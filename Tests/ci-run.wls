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

testFiles = Sort@FileNames["*.wlt", FileNameJoin[{root, "Tests"}]];
Print["Test files: ", FileNameTake /@ testFiles];

reports = TestReport /@ testFiles;
outcomes = Catenate[Through[Values[#["TestResults"]]["Outcome"]] & /@ reports];
Print["Test outcomes: ", Tally[outcomes]];

If[outcomes =!= {} && AllTrue[outcomes, # === "Success" &],
    Print["All ", Length[outcomes], " tests passed."];
    Exit[0],
    Print["Test failures detected (or no tests ran)."];
    Exit[1]];
