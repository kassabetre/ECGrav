#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Continuous-integration runner for the ECGrav .wlt regression suite.

   Runs every test file under Tests/ through MUnit`TestReport and sets the process exit
   code: 0 if every test succeeded, 1 otherwise -- so CI (e.g. GitHub Actions) turns the
   job red on any regression. Invoke with:

       wolframscript -file Tests/ci-run.wls

   The repo root is taken from the ECGRAV_ROOT environment variable if it is already set,
   otherwise from this script's own location (Tests/ci-run.wls -> repo root). It is then
   exported back into ECGRAV_ROOT so TestPrelude.wl can locate the package from source
   regardless of where the repository is checked out. *)

Needs["MUnit`"];

(* Pick the first candidate that actually contains the package: an explicit ECGRAV_ROOT,
   this script's grandparent (Tests/ci-run.wls -> repo root), or the current directory
   (which GitHub Actions sets to the checkout root). *)
hasPackage[dir_] := StringQ[dir] && dir =!= "" &&
    FileExistsQ[FileNameJoin[{dir, "Kernel", "ECGrav.wl"}]];
candidates = {
    Environment["ECGRAV_ROOT"],
    If[$InputFileName =!= "", ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]], Nothing],
    Directory[]};
root = SelectFirst[candidates, hasPackage, Directory[]];
SetEnvironment["ECGRAV_ROOT" -> root];

Print["ECGrav repo root: ", root];

(* Environment preflight. The log has to be able to answer "was it the code or the box?"
   on its own, because a red run here is usually the box: the container image is
   :latest and floats between runs, and CountHoles wraps ResourceFunction["BettiNumbers"],
   which a fresh container must fetch from the Wolfram Cloud every time -- there is no
   cached copy in the image. LaunchKernels is checked because TestPrelude calls it and
   LowEnergyStates divides by $KernelCount. *)
(* GitHub Actions workflow commands. Job logs need push access to read through the API,
   but ANNOTATIONS are readable unauthenticated, so the failure summary is emitted as
   ::error:: lines too. That is the difference between diagnosing a red build and asking
   someone to paste the log. Newlines must be percent-encoded, and the "%" itself first.
   GitHub renders at most ten annotations per level per step, hence the cap below. *)
inActions = Environment["GITHUB_ACTIONS"] === "true";
ghEnc[t_String] := StringReplace[t, {"%" -> "%25", "\r" -> "%0D", "\n" -> "%0A"}];
ghCmd[kind_String, msg_String] := If[inActions, Print["::" <> kind <> "::" <> ghEnc[msg]]];

pre = {};
addPre[label_String, value_] := (
    Print[StringPadRight[label, 17], ": ", value];
    AppendTo[pre, label <> "=" <> ToString[value]]);
addPre["Wolfram version", $Version];
addPre["ProcessorCount", $ProcessorCount];
Module[{k = Quiet[LaunchKernels[]]},
    addPre["LaunchKernels", ToString[Length[If[ListQ[k], k, {}]]] <>
        " kernels, $KernelCount=" <> ToString[$KernelCount]]];
(* No network probe any more. Until 1.14.1 CountHoles wrapped
   ResourceFunction["BettiNumbers"], the package's only network dependency, and this
   preflight probed it -- which is how the 15 red CountHoles tests of 2026-09-05 were
   diagnosed from the annotations alone. Betti numbers are now computed locally, so
   there is nothing left here that needs the Wolfram Cloud. If that ever changes, probe
   it here again: a red build whose cause is the box should say so in its own log. *)
addPre["Network needed", "no (Betti numbers computed locally since 1.14.1)"];
ghCmd["notice", "preflight: " <> StringRiffle[pre, " | "]];

testFiles = Sort@FileNames["*.wlt", FileNameJoin[{root, "Tests"}]];
Print["Test files: ", FileNameTake /@ testFiles];

reports = TestReport /@ testFiles;
outcomes = Catenate[Through[Values[#["TestResults"]]["Outcome"]] & /@ reports];

(* Per file, so a red run names the file before anything else. *)
Print[""];
Do[Print["  ", FileNameTake[testFiles[[i]]], ": ",
        Tally[Through[Values[reports[[i]]["TestResults"]]["Outcome"]]]],
    {i, Length[testFiles]}];
Print["Test outcomes: ", Tally[outcomes]];

(* Name every failing test. Tally alone says "8 failures" and never which, which makes a
   red log almost as uninformative as no log. Details are capped so that a wholesale
   breakage cannot bury the summary. *)
(* Rendering, for a log rather than a notebook. ToString[Short[x, n]] prints the literal
   Short[...] wrapper instead of shortening, so cap the rendered string. Plain ToString
   uses OutputForm, which typesets 1/0 across three lines -- unreadable interleaved with
   other CI output -- so force InputForm at infinite page width for one line per value.
   The HoldForm wrappers that survive are TestResultObject's own and are left alone: they
   are noise, but resolving them would mean evaluating the very expression that failed. *)
clip[x_] := With[{s = ToString[x, InputForm, PageWidth -> Infinity]},
    If[StringLength[s] > 220, StringTake[s, 200] <> " ...(" <> ToString[StringLength[s]] <> " chars)", s]];

Module[{failed, shown, cap = 25},
    failed = Catenate[Function[i,
            Function[t, <|"file" -> FileNameTake[testFiles[[i]]], "test" -> t|>] /@
                Select[Values[reports[[i]]["TestResults"]], #["Outcome"] =!= "Success" &]] /@
        Range[Length[testFiles]]];
    If[failed =!= {},
        Print[""];
        Print["FAILING TESTS (", Length[failed], "):"];
        shown = Take[failed, UpTo[cap]];
        Do[Module[{t = f["test"]},
            Print["  [", f["file"], "] ", t["TestID"], "  (", t["Outcome"], ")"];
            Print["      expected : ", clip[t["ExpectedOutput"]]];
            Print["      actual   : ", clip[t["ActualOutput"]]];
            If[t["ActualMessages"] =!= {},
                Print["      messages : ", clip[t["ActualMessages"]]]]],
            {f, shown}];
        If[Length[failed] > cap,
            Print["  ... and ", Length[failed] - cap, " more, not shown"]];

        (* Same information as annotations, which are publicly readable. One summary
           naming every failing test, then a few with detail, staying under the
           ten-per-step render limit. *)
        ghCmd["error", "FAILING (" <> ToString[Length[failed]] <> "): " <>
            StringRiffle[ToString[#["test"]["TestID"]] & /@ failed, ", "]];
        Do[Module[{t = f["test"]},
            ghCmd["error", f["file"] <> "  " <> ToString[t["TestID"]] <> "  (" <> ToString[t["Outcome"]] <> ")" <>
                "\nexpected: " <> clip[t["ExpectedOutput"]] <>
                "\nactual  : " <> clip[t["ActualOutput"]] <>
                If[t["ActualMessages"] =!= {}, "\nmessages: " <> clip[t["ActualMessages"]], ""]]],
            {f, Take[failed, UpTo[8]]}]]];

Print[""];
If[outcomes =!= {} && AllTrue[outcomes, # === "Success" &],
    Print["All ", Length[outcomes], " tests passed."];
    Exit[0],
    Print["Test failures detected (or no tests ran)."];
    Exit[1]];
