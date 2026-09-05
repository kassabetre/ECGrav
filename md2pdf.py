#!/usr/bin/env python3
"""Render the repo's Markdown specs to PDF with real LaTeX maths.

    python3 md2pdf.py Theory.md              # -> Theory.pdf
    python3 md2pdf.py                        # -> every *.md at the repo root
    python3 md2pdf.py --check FacetLabeledCount.md   # build + verify, no install

There is no pandoc on this machine (and no Homebrew to get one), but there is a
full TeX distribution at /Library/TeX/texbin, so this converts Markdown to LaTeX
directly and runs pdflatex.  Ghostscript (/usr/local/bin/gs) is used to extract
text back out of the PDF for the verification pass.

Handles exactly what these documents use: ATX headings, $$..$$ display maths with
\\tag, $..$ inline maths, pipe tables, blockquotes, bullet/numbered/checkbox
lists, fenced code, thematic breaks, links, **bold**, *em*, `code`.

Two traps this is built around, both of which produce a PDF that compiles with
exit 0 and is wrong:

  1. ESCAPING ORDER.  Escape the ASCII LaTeX specials FIRST, then map non-ASCII
     to LaTeX commands.  Doing it the other way round re-escapes the backslashes
     of the commands just inserted, and every em dash renders as the literal text
     "\\textemdash{}".  Nothing in the compile log mentions it.
  2. PIPES INSIDE CODE SPANS.  A table cell like `1 / |Aut|` contains a pipe, so
     the row must have its code spans stashed away BEFORE it is split on "|",
     or the table silently gains phantom columns.

So: never trust the exit code.  --check greps the extracted text for escaping
leakage and confirms every heading and equation tag survived.
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile

TEXBIN = "/Library/TeX/texbin"
GS = "/usr/local/bin/gs"

# LaTeX specials, escaped before any non-ASCII mapping runs (see trap 1).
SPECIAL = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
           "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
           "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}

# Non-ASCII -> LaTeX, for prose.  Every character occurring in the repo's specs.
UNI = {
    "—": r"\textemdash{}", "–": r"\textendash{}", "§": r"\S{}", "·": r"\textperiodcentered{}",
    "…": r"\dots{}", "′": r"$'$", "½": r"$\tfrac12$", "²": r"$^2$", "³": r"$^3$", "⁷": r"$^7$",
    "₀": r"$_0$", "₁": r"$_1$", "₄": r"$_4$",
    "ö": r'\"o', "ó": r"\'o", "ü": r'\"u', "é": r"\'e",
    "Γ": r"$\Gamma$", "Σ": r"$\Sigma$", "Δ": r"$\Delta$", "σ": r"$\sigma$", "λ": r"$\lambda$",
    "χ": r"$\chi$", "ρ": r"$\rho$", "β": r"$\beta$", "ν": r"$\nu$", "δ": r"$\delta$",
    "π": r"$\pi$", "ℓ": r"$\ell$", "𝒜": r"$\mathcal{A}$",
    "×": r"$\times$", "→": r"$\to$", "↔": r"$\leftrightarrow$", "−": r"$-$", "±": r"$\pm$",
    "≤": r"$\le$", "≥": r"$\ge$", "≠": r"$\ne$", "≈": r"$\approx$", "≡": r"$\equiv$",
    "∝": r"$\propto$", "∈": r"$\in$", "∩": r"$\cap$", "∪": r"$\cup$", "⋂": r"$\bigcap$",
    "⊆": r"$\subseteq$", "√": r"$\sqrt{\,}$", "∎": r"$\blacksquare$",
    "⟨": r"$\langle$", "⟩": r"$\rangle$",
    "✅": r"$\checkmark$", "►": r"$\blacktriangleright$",
    # box drawing, only ever inside fenced diagrams
    "─": "-", "┐": "+", "┘": "+", "┼": "+", "├": "+", "│": "|", "└": "+", "┌": "+",
}
# Inside verbatim nothing may be a command, so fall back to ASCII there.
MONO = dict(UNI)
MONO.update({k: v for k, v in {
    "—": "--", "–": "-", "§": "S", "·": ".", "…": "...", "′": "'", "½": "1/2",
    "²": "^2", "³": "^3", "⁷": "^7", "₀": "_0", "₁": "_1", "₄": "_4",
    "ö": "oe", "ó": "o", "ü": "ue", "é": "e",
    "Γ": "Gamma", "Σ": "Sigma", "Δ": "Delta", "σ": "sigma", "λ": "lambda", "χ": "chi",
    "ρ": "rho", "β": "beta", "ν": "nu", "δ": "delta", "π": "pi", "ℓ": "l", "𝒜": "A",
    "×": "x", "→": "->", "↔": "<->", "−": "-", "±": "+/-", "≤": "<=", "≥": ">=",
    "≠": "!=", "≈": "~=", "≡": "==", "∝": "prop", "∈": "in", "∩": "^", "∪": "v",
    "⋂": "^", "⊆": "subset", "√": "sqrt", "∎": "QED", "⟨": "<", "⟩": ">",
    "✅": "[x]", "►": ">",
}.items()})

PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb}
\usepackage{booktabs,tabularx,longtable}
\usepackage[margin=1in]{geometry}
\usepackage{fancyvrb}
\usepackage[hidelinks]{hyperref}
\usepackage{microtype}
\sloppy
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt}
\setlength{\emergencystretch}{3em}
\renewcommand{\arraystretch}{1.15}
\author{}\date{}
\title{%s}
\begin{document}
\maketitle
"""


class Vault:
    """Stash verbatim spans behind opaque keys so escaping cannot reach them."""

    def __init__(self):
        self.d = {}

    def put(self, s):
        k = "ZQV" + hashlib.md5(s.encode()).hexdigest()[:14] + "ZQV"
        self.d[k] = s
        return k

    def restore(self, s):
        for _ in range(6):                      # nested placeholders
            if not any(k in s for k in self.d):
                break
            for k, v in self.d.items():
                s = s.replace(k, v)
        return s


def esc(s, mono=False):
    """ASCII specials first, then non-ASCII -> commands.  Order matters (trap 1)."""
    if mono:
        for u, r in MONO.items():
            s = s.replace(u, r)
        return s
    s = "".join(SPECIAL.get(c, c) for c in s)
    for u, r in UNI.items():
        s = s.replace(u, r)
    return s


def breakable(s):
    """Let long paths in \\texttt wrap."""
    return s.replace("/", r"/\allowbreak{}").replace(r"\_", r"\_\allowbreak{}")


def inline(s, v):
    """Inline markup -> LaTeX.  Maths and code are vaulted before escaping."""
    s = re.sub(r"\$([^$\n]+)\$", lambda m: v.put("$" + m.group(1) + "$"), s)
    s = re.sub(r"`([^`]+)`",
               lambda m: v.put(r"\texttt{" + breakable(esc(m.group(1))) + "}"), s)
    def link(m):
        text, url = m.group(1), m.group(2)
        if url.startswith(("http://", "https://")):
            return v.put(r"\href{" + url.replace("%", r"\%").replace("#", r"\#")
                         + "}{") + text + v.put("}")
        return text + v.put(r"\,\texttt{\footnotesize(" + breakable(esc(url)) + ")}")
    s = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", link, s)
    s = re.sub(r"\*\*(.+?)\*\*", lambda m: v.put(r"\textbf{") + m.group(1) + v.put("}"), s)
    s = re.sub(r"(?<![\w*])\*([^*\n]+?)\*(?![\w*])",
               lambda m: v.put(r"\emph{") + m.group(1) + v.put("}"), s)
    return esc(s)


def split_cells(row, v):
    """Split a table row on '|' AFTER vaulting code spans, so a pipe inside a
    code span cannot invent a column (trap 2).

    The span must be vaulted as its FINAL LaTeX, not as the raw backticked
    text: a vaulted token is invisible to inline(), so stashing "`x`" here
    means nothing ever converts it and the cell renders with literal
    backticks around unescaped source."""
    return [c.strip() for c in
            re.sub(r"`([^`]+)`",
                   lambda m: v.put(r"\texttt{" + breakable(esc(m.group(1))) + "}"),
                   row).strip().strip("|").split("|")]


def table(rows, v):
    hdr = split_cells(rows[0], v)
    body = [split_cells(r, v) for r in rows[2:]]
    n = len(hdr)
    body = [(r + [""] * n)[:n] for r in body if any(c.strip() for c in r)]
    widest = [max([len(hdr[i])] + [len(r[i]) for r in body] or [0]) for i in range(n)]
    wide = [i for i, w in enumerate(widest) if w > 30]
    spec = "".join("X" if i in wide else "l" for i in range(n))
    if not wide:
        spec = "l" * n
    env = "tabularx"
    out = [r"\begin{center}\small",
           r"\begin{%s}{\linewidth}{%s}" % (env, spec) if wide
           else r"\begin{tabular}{%s}" % spec,
           r"\toprule",
           " & ".join(r"\textbf{" + inline(h, v) + "}" for h in hdr) + r" \\",
           r"\midrule"]
    for r in body:
        out.append(" & ".join(inline(c, v) for c in r) + r" \\")
    out += [r"\bottomrule",
            r"\end{%s}" % env if wide else r"\end{tabular}",
            r"\end{center}"]
    return out


def convert(md, title):
    v = Vault()

    # Fenced code -> Verbatim, ASCII-folded so pdflatex never sees a stray byte.
    def fence(m):
        body = "\n".join(esc(l, mono=True) for l in m.group(2).rstrip("\n").split("\n"))
        return v.put("\n\\begin{Verbatim}[fontsize=\\small,frame=leftline,"
                     "framesep=6pt,xleftmargin=4pt]\n" + body + "\n\\end{Verbatim}\n")
    md = re.sub(r"```[ \t]*(\w*)\n(.*?)```", fence, md, flags=re.S)

    # Display maths.  \tag only works in a numbered environment.
    def disp(m):
        body = m.group(1).strip()
        env = "equation" if r"\tag{" in body else "equation*"
        return v.put("\n\\begin{%s}\n%s\n\\end{%s}\n" % (env, body, env))
    md = re.sub(r"\$\$(.*?)\$\$", disp, md, flags=re.S)

    out, lines, i = [], md.split("\n"), 0
    while i < len(lines):
        ln = lines[i]
        if not ln.strip():
            out.append("")
            i += 1
            continue
        if re.match(r"^\s*(---+|\*\*\*+|___+)\s*$", ln):
            out.append(r"\vspace{2mm}\hrule\vspace{2mm}")
            i += 1
            continue
        m = re.match(r"^(#{1,6})\s+(.*)$", ln)
        if m:
            lvl, txt = len(m.group(1)), inline(m.group(2), v)
            if lvl == 1:
                pass                                    # H1 is the title
            else:
                cmd = {2: r"\section*{%s}", 3: r"\subsection*{%s}"}.get(
                    lvl, r"\subsubsection*{%s}")
                out.append(cmd % txt)
            i += 1
            continue
        if ln.lstrip().startswith("|"):
            j = i
            while j < len(lines) and lines[j].lstrip().startswith("|"):
                j += 1
            if j - i >= 2:
                out += table(lines[i:j], v)
            i = j
            continue
        if ln.startswith(">"):
            j = i
            while j < len(lines) and lines[j].startswith(">"):
                j += 1
            body = " ".join(re.sub(r"^>\s?", "", l).strip() for l in lines[i:j])
            out += [r"\begin{quote}\itshape", inline(body, v), r"\end{quote}"]
            i = j
            continue
        m = re.match(r"^(\s*)([-*+]|\d+[.)])\s+(.*)$", ln)
        if m:
            env = "itemize" if m.group(2) in "-*+" else "enumerate"
            out.append(r"\begin{%s}" % env)
            j = i
            while j < len(lines):
                mm = re.match(r"^(\s*)([-*+]|\d+[.)])\s+(.*)$", lines[j])
                if mm:
                    txt = mm.group(3)
                    box = re.match(r"^\[([ xX])\]\s*(.*)$", txt)
                    prefix = ""
                    if box:
                        prefix = (r"$\checkmark$~" if box.group(1).lower() == "x"
                                  else r"$\square$~")
                        txt = box.group(2)
                    j += 1
                    while (j < len(lines) and lines[j].strip()
                           and lines[j].startswith((" ", "\t"))
                           and not re.match(r"^\s*([-*+]|\d+[.)])\s+", lines[j])
                           and not lines[j].lstrip().startswith(("|", ">", "#"))):
                        txt += " " + lines[j].strip()
                        j += 1
                    out.append(r"\item " + v.put(prefix) + inline(txt, v))
                elif (not lines[j].strip() and j + 1 < len(lines)
                      and re.match(r"^\s*([-*+]|\d+[.)])\s+", lines[j + 1])):
                    j += 1
                else:
                    break
            out.append(r"\end{%s}" % env)
            i = j
            continue
        para = [ln]
        i += 1
        while (i < len(lines) and lines[i].strip()
               and not re.match(r"^(#|>|\s*(---+|\*\*\*+)\s*$|\s*([-*+]|\d+[.)])\s)", lines[i])
               and not lines[i].lstrip().startswith("|")):
            para.append(lines[i])
            i += 1
        out.append(inline(" ".join(para), v))

    return PREAMBLE % title + v.restore("\n".join(out)) + "\n\\end{document}\n"


def build(src, install=True, check=True):
    """Convert, compile, verify.  Returns (ok, message)."""
    md = open(src, encoding="utf-8").read()
    m = re.search(r"^#\s+(.*)$", md, re.M)
    title = esc(re.sub(r"[`*]", "", m.group(1))) if m else esc(
        os.path.splitext(os.path.basename(src))[0])
    stem = os.path.splitext(os.path.basename(src))[0]
    env = dict(os.environ, PATH=TEXBIN + ":" + os.environ.get("PATH", ""))

    with tempfile.TemporaryDirectory() as d:
        tex = os.path.join(d, stem + ".tex")
        open(tex, "w", encoding="utf-8").write(convert(md, title))
        for _ in range(2):                              # twice, for hyperref
            subprocess.run(["pdflatex", "-interaction=nonstopmode", stem + ".tex"],
                           cwd=d, env=env, capture_output=True)
        pdf = os.path.join(d, stem + ".pdf")
        if not os.path.exists(pdf):
            log = open(os.path.join(d, stem + ".log"), encoding="utf-8",
                       errors="replace").read()
            errs = [l for l in log.split("\n") if l.startswith("!")][:3]
            return False, "pdflatex produced no PDF: " + "; ".join(errs)

        log = open(os.path.join(d, stem + ".log"), encoding="utf-8",
                   errors="replace").read()
        pages = re.search(r"Output written .*?\((\d+) page", log)
        pages = pages.group(1) if pages else "?"
        overfull = log.count("Overfull \\hbox")
        missing = len(re.findall(r"Missing character", log))
        notes = []

        if check and os.path.exists(GS):
            txt = os.path.join(d, "out.txt")
            subprocess.run([GS, "-dNOPAUSE", "-dBATCH", "-sDEVICE=txtwrite",
                            "-sOutputFile=" + txt, pdf], cwd=d, capture_output=True)
            body = open(txt, encoding="utf-8", errors="replace").read()
            flat = re.sub(r"\s+", " ", body)
            leaks = [t for t in ("ZQV", "textbackslash", "textemdash", r"\tag{",
                                 "begin{equation", "$$") if t in flat]
            if leaks:
                notes.append("LEAK:" + ",".join(leaks))
            heads = re.findall(r"^#{2,3}\s+(.*)$", md, re.M)
            miss = 0
            # Normalise BOTH sides to bare lowercase alphanumerics before
            # comparing. Anything less fails spuriously: punctuation sits
            # between the words on one side but not the other, pdflatex
            # hyphenates across line breaks, and maths normalises differently
            # in the source than in the extracted text.
            # Squash the WHOLE heading on both sides. Selecting words fails:
            # drop the digits and "Phase 0 - Foundation" stops matching
            # "phase0foundation"; keep the punctuation and "Two conditions four"
            # stops matching "Two conditions, four". Only removing every
            # non-alphanumeric from both sides is symmetric.
            #
            # Two artifacts of gs text extraction have to be absorbed, or this
            # reports failures against PDFs that are perfectly correct:
            #   - LaTeX ligatures come back as a single unmapped glyph, so
            #     "coefficient" extracts as "coefIcient";
            #   - subscripts are emitted after their base in the content stream,
            #     so "S_n, not S_n x S_M" extracts as "S , not S xS n n M".
            # The first is repaired below; the second is why the match is an
            # ordered subsequence inside a local window rather than contiguous.
            for lig, rep in (("Ï", "ffi"), ("ﬀ", "ff"), ("ﬁ", "fi"),
                             ("ﬂ", "fl"), ("ﬃ", "ffi"), ("ﬄ", "ffl")):
                flat = flat.replace(lig, rep)
            squash = lambda s: re.sub(r"[^a-z0-9]+", "", s.lower())
            flatn = squash(flat)

            def present(key):
                if key in flatn:
                    return True
                head, span = key[:8], len(key) * 3 + 16
                start = flatn.find(head)
                while start != -1:
                    win, k = flatn[start:start + span], 0
                    for ch in win:
                        if k < len(key) and ch == key[k]:
                            k += 1
                    if k == len(key):
                        return True
                    start = flatn.find(head, start + 1)
                return False

            for h in heads:
                key = squash(re.sub(r"\\[a-zA-Z]+", " ", h))
                if key and not present(key):
                    miss += 1
            if miss:
                # ADVISORY, never a failure. Maths-bearing headings cannot be
                # checked this way: gs emits every base glyph before its
                # subscripts, so "$S_n$, not $S_n x S_M$" extracts in a
                # genuinely permuted order ("snotssnnm" against "snnotsnsm").
                # No ordered match bridges that, and a character-multiset
                # comparison would be too weak to mean anything. Treat a hit
                # here as a prompt to eyeball the page, not as a broken build.
                notes.append("headings unmatched (advisory): %d/%d" % (miss, len(heads)))
            tags = re.findall(r"\\tag\{([^}]*)\}", md)
            lost = [t for t in tags if "(" + t + ")" not in flat]
            if lost:
                notes.append("tags missing: " + ",".join(lost))

        if missing:
            notes.append("missing glyphs: %d" % missing)
        if install:
            shutil.copy(pdf, os.path.join(os.path.dirname(os.path.abspath(src)) or ".",
                                          stem + ".pdf"))
    # Hard failures only. The heading scan is advisory (see above); leakage,
    # dropped equation tags and missing glyphs are reliable and are not.
    ok = not any(n.startswith(("LEAK", "tags", "missing glyphs")) for n in notes)
    msg = "%s pages, %d overfull" % (pages, overfull)
    if notes:
        msg += " | " + "; ".join(notes)
    return ok, msg


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    install = "--check" not in sys.argv
    srcs = args or sorted(f for f in os.listdir(".") if f.endswith(".md"))
    width = max(len(s) for s in srcs) + 1
    bad = 0
    for s in srcs:
        ok, msg = build(s, install=install)
        print("%-*s %s  %s" % (width, s, "ok  " if ok else "FAIL", msg))
        bad += not ok
    print("\n%d/%d built cleanly" % (len(srcs) - bad, len(srcs)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
