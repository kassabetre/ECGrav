#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Reproducible build script for the ECGrav paclet.

   Replaces the manual "Building a paclet" steps in BuildingAndInstallingPaclets.nb.
   It cleans macOS .DS_Store junk from the content directories, builds the paclet with
   PacletTools`PacletBuild (which also regenerates the documentation search / spell
   indexes under build/ECGrav/), and then VERIFIES the artifact:
     - the documentation build reported no failures,
     - the archive is named for the version declared in PacletInfo.wl,
     - the bundled PacletInfo.wl + Kernel sources are byte-identical to the working tree,
     - no .DS_Store leaked into the archive.
   A summary is printed and the process exit code is set (0 = success, 1 = any failure),
   so this can gate a release or run in CI.

   Run:     wolframscript -file build.wls
   Output:  build/ECGrav-<version>.paclet   (<version> comes from PacletInfo.wl)

   The repo root is ECGRAV_ROOT if set, else this script's own directory (or its parent),
   else the current directory.

   This script deliberately does NOT install the paclet. PacletUninstall removes a paclet
   by NAME (all versions), so install-testing in a shared Wolfram environment can delete an
   already-installed copy. To install the freshly built paclet:
       PacletInstall["<repo>/build/ECGrav-<version>.paclet"]
*)

Needs["PacletTools`"];

hasPackage[dir_] := StringQ[dir] && dir =!= "" &&
    FileExistsQ[FileNameJoin[{dir, "PacletInfo.wl"}]] &&
    FileExistsQ[FileNameJoin[{dir, "Kernel", "ECGrav.wl"}]];

scriptDir = If[$InputFileName =!= "", DirectoryName[ExpandFileName[$InputFileName]], ""];
root = SelectFirst[
    {Environment["ECGRAV_ROOT"], scriptDir, ParentDirectory[scriptDir], Directory[]},
    hasPackage, $Failed];
If[root === $Failed,
    Print["ERROR: could not locate the ECGrav paclet root (need PacletInfo.wl and Kernel/ECGrav.wl)."];
    Exit[1]];

buildDir = FileNameJoin[{root, "build"}];
version = Get[FileNameJoin[{root, "PacletInfo.wl"}]]["Version"];
Print["ECGrav paclet root : ", root];
Print["Declared version   : ", version];

(* 1. Remove .DS_Store from the content dirs so the build manifest stays clean. *)
junk = FileNames[".DS_Store",
    Select[{FileNameJoin[{root, "Kernel"}], FileNameJoin[{root, "Documentation"}]}, DirectoryQ], Infinity];
If[junk =!= {}, DeleteFile[junk]];
Print["Removed .DS_Store  : ", Length[junk]];

(* 2. Build: regenerates build/ECGrav/ and writes build/ECGrav-<version>.paclet. *)
Print["Building ..."];
res = PacletTools`PacletBuild[root, buildDir];
If[Head[res] =!= Success,
    Print["ERROR: PacletBuild did not succeed: ", res]; Exit[1]];
archive = res["PacletArchive"];
docFailures = Total[(#["FailedFilesCount"] &) /@ res["DocumentationBuildResults"]];
Print["Archive            : ", archive];
Print["Doc build failures : ", docFailures];

(* 3. Verify the artifact against the working tree (reads entries straight from the ZIP). *)
expectedName = "ECGrav-" <> version;
identityFiles = {"PacletInfo.wl", "Kernel/ECGrav.wl",
    "Kernel/Subpackages/MCSims.wl", "Kernel/Subpackages/PureComplexes.wl"};
mismatched = Select[identityFiles,
    Import[archive, {"ZIP", expectedName <> "/" <> #, "Byte"}] =!= Import[FileNameJoin[{root, #}], "Byte"] &];
dsInArchive = Select[Import[archive, {"ZIP", "FileNames"}], StringContainsQ[#, "DS_Store"] &];
nameOK = FileBaseName[archive] === expectedName;

Print["Archive name OK    : ", nameOK];
Print["Sources match tree : ", mismatched === {},
    If[mismatched =!= {}, " (MISMATCH: " <> ToString[mismatched] <> ")", ""]];
Print[".DS_Store in paclet: ", Length[dsInArchive]];

Print[""];
If[FileExistsQ[archive] && docFailures === 0 && nameOK && mismatched === {} && dsInArchive === {},
    Print["BUILD OK -> ", archive];
    Print["Install with: PacletInstall[\"", archive, "\"]"];
    Exit[0],
    Print["BUILD FAILED (see the checks above)."];
    Exit[1]];
