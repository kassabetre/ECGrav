#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*WLSToNotebook -- turn a .wls script into an evaluatable .nb*)

(*
	Usage:

		wolframscript -file WLSToNotebook.wls  in.wls  [out.nb]

	or from a notebook:

		Get["WLSToNotebook.wls"];
		WLSToNotebook["in.wls"]          (* writes in.nb next to it *)

	The front end already opens a .wls directly (double-click, or File > Open), renders the
	(* ::Section:: *) markers as real section cells, and will File > Save As it to .nb.  This is
	the same conversion done in the kernel, for when you want it scripted or the file is on a
	machine with no front end.

	The kernel has no built-in for it: Import[file, "Notebook"] on a .wls just returns the text,
	because the marker-to-cell-style mapping lives in the front end, not in Import.

	What it does: a (* ::Style:: *) marker followed by a (*text*) line becomes a cell of that
	style; the leading block comment becomes a Text cell; everything else is grouped into Input
	cells, split at blank lines between top-level expressions.
*)

BeginPackage["WLSToNotebook`"];

WLSToNotebook::usage =
"WLSToNotebook[\"in.wls\"] writes in.nb beside it and returns the path.  \
WLSToNotebook[\"in.wls\", \"out.nb\"] writes to out.nb.";

Begin["`Private`"];

(*A marker line is (* ::Section:: *) and the (*title*) line under it carries the text.*)
markerStyle[line_String] :=
	Replace[StringCases[line, "(* ::" ~~ s : Except[":"] .. ~~ ":: *)" :> s],
		{{s_} :> s, _ -> None}];

commentText[line_String] :=
	Replace[StringCases[line, "(*" ~~ s___ ~~ "*)" :> s], {{s_} :> StringTrim[s], _ -> None}];

(*Input cells hold plain strings.  The front end parses them on evaluation, which keeps the code
	byte-identical to the script -- no MakeBoxes round trip to lose a comment or reformat a line.*)
inputCell[code_String] :=
	If[StringTrim[code] === "", Nothing, Cell[StringTrim[code, "\n" | " "], "Input"]];

WLSToNotebook[in_String] :=
	WLSToNotebook[in, StringReplace[in, ".wls" ~~ EndOfString -> ".nb"]];

WLSToNotebook[in_String, out_String] :=
Module[{lines, cells = {}, buf = {}, i, n, style, txt, depth, j, block},
	lines = StringSplit[Import[in, "Text"], "\n"];
	lines = DeleteCases[lines, _?(StringStartsQ[#, "#!"] &), 1, 1];

	flush[] := (cells = Join[cells, {inputCell[StringRiffle[buf, "\n"]]}]; buf = {});

	n = Length[lines];
	i = 1;
	While[i <= n,
		style = markerStyle[lines[[i]]];
		Which[
			(*a style marker: the next non-blank line is its (*text*)*)
			style =!= None && style =!= "Package",
				flush[];
				j = i + 1;
				While[j <= n && StringTrim[lines[[j]]] === "", j++];
				txt = If[j <= n, commentText[lines[[j]]], None];
				If[txt =!= None,
					AppendTo[cells, Cell[txt, style]]; i = j + 1,
					i = i + 1],

			style === "Package", i++,

			(*a standalone block comment becomes prose, but only OUTSIDE code -- one already
				buffered means it is a comment on the code, so leave it in the Input cell*)
			buf === {} && StringStartsQ[StringTrim[lines[[i]]], "(*"],
				depth = 0; j = i;
				While[j <= n,
					depth += StringCount[lines[[j]], "(*"] - StringCount[lines[[j]], "*)"];
					If[depth <= 0, Break[]];
					j++];
				block = StringRiffle[lines[[i ;; Min[j, n]]], "\n"];
				AppendTo[cells,
					Cell[StringTrim@StringReplace[block,
						{StartOfString ~~ "(*" -> "", "*)" ~~ EndOfString -> ""}], "Text"]];
				i = Min[j, n] + 1,

			(*a blank line ends an Input cell only when brackets are balanced, so a blank line
				inside a Module or a Table does not split the expression in half*)
			StringTrim[lines[[i]]] === "",
				If[buf =!= {} && balancedQ[StringRiffle[buf, "\n"]], flush[]];
				i++,

			True, AppendTo[buf, lines[[i]]]; i++]];
	flush[];

	Export[out, Notebook[cells, WindowSize -> {900, 800}], "Notebook"];
	out
];

(*Brackets outside strings and comments.  Cheap and good enough to decide where a cell ends.*)
balancedQ[code_String] :=
Module[{s = StringReplace[code, {"(*" ~~ Shortest[___] ~~ "*)" -> "",
			"\"" ~~ Shortest[___] ~~ "\"" -> "\"\""}]},
	And @@ Table[StringCount[s, p[[1]]] === StringCount[s, p[[2]]],
		{p, {{"[", "]"}, {"{", "}"}, {"(", ")"}}}]
];

End[];
EndPackage[];

If[!$Notebooks && Length[$ScriptCommandLine] >= 2,
	Print[WLSToNotebook @@ Rest[$ScriptCommandLine]]];
