#!/usr/bin/env wolframscript

(* 
Objective: plot the non convex tilings on
https://en.wikipedia.org/wiki/Uniform_tiling#Expanded_lists_of_uniform_tilings
 *)


frameCount = 5;
imageSize = {4, 3} * 720 / 3;
shape = "labs_septic";
imageFolder = "output";
exportToPov = True;
plotPoints = 60;
maxRecursion = 2;
amax = 0.15;
amin = -0.45;

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

pFunction[x_, y_, z_] := 
  x (x^6 - 3*7*x^4 y^2 + 5*7*x^2 y^4 - 7*y^6) + 7 z ((x^2 + y^2)^3 - 2^3 z^2 (x^2 + y^2)^2 + 2^4 z^4 (x^2 + y^2)) - 2^6 z^7;
a1[a_] := -12/7 a^2 - 384/49 a - 8/7;
a2[a_] := -32/7 a^2 + 24/49 a - 4;
a3[a_] := -4 a^2 + 24/49 a - 4;
a4[a_] := -8/7 a^2 + 8/49 a - 8/7;
a5[a_] := 49 a^2 - 7 a + 50;
uFunction[x_, y_, z_, a_, w_] := (z + a5[a] w) ((z + w) (x^2 + y^2) + a1[a] z^3 + a2[a] z^2 w + a3[a] z w^2 + a4[a] w^3)^2;
sFunction[x_, y_, z_, a_, w_] := pFunction[x, y, z] - uFunction[x, y, z, a, w];

shapeWithParams = shape <> "_" <> IntegerString[plotPoints] <> "_" <> IntegerString[maxRecursion];

exportFrame[aValue_, frameIndex_, outputFolderName_] := Module[
	{},

	range = 2;
	renderRange = 3 * range;

(* 	graph = ContourPlot3D[
	 x^2+y^2+100 z^2 - range^2 == 0, {x, -range, range}, {y, -range, 
	  range}, {z, -range, range}, 
	 RegionFunction -> Function[{x, y, z}, x^2 + y^2 + z^2 <= range^2], 
	 RegionBoundaryStyle -> None, PlotPoints -> plotPoints, 
	 Boxed -> False, Axes -> False, MaxRecursion -> maxRecursion, 
	 MeshFunctions -> {0 &}, MeshShading -> { { Pink } },
	 ViewPoint -> {50, -50, 50}
	]; *)

	graph = ContourPlot3D[
	 sFunction[x, y, z, aValue, 1] == 0, {x, -renderRange, renderRange}, {y, -renderRange, 
	  renderRange}, {z, -renderRange, renderRange}, 
	 RegionBoundaryStyle -> None, PlotPoints -> plotPoints, 
	 Boxed -> False, Axes -> False, MaxRecursion -> maxRecursion, 
	 MeshFunctions -> {0 &}, MeshShading -> { { Pink } },
	 ViewPoint -> {70, -40, 90}
	];

	graph = Show[graph, PlotRange->{{-range, range},{-range, range},{-range, range}}];

	outputFolder = FileNameJoin[{imageFolder, outputFolderName}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	outputFileName = "spinning_frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
	outputFileName = FileNameJoin[{outputFolder, outputFileName}];

	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	Print["Saving an image to " <> outputFileName];
	If[exportToPov,
	  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
	  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
	  Run[cropCommand];
	  ,
	  Export[outputFileName, graph, ImageSize -> imageSize[[2]]];
	];
	Print["Saved an image to " <> outputFileName];
];

exportAnimation[] := Module[
	{
		outputFolder,
		gifFileName,
		gifFolder,
		maxExplosionFactor,
		explodeStep,
		explosionFactors,
		frameIndex,
		explosionFactor,
		gifCommand,
		edges
	},

	outputFolderName = shapeWithParams <> "_spinning_frames";
	outputFolder = FileNameJoin[{imageFolder, outputFolderName}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	gifFileName = shapeWithParams <> "_spinning";
	gifFolder = FileNameJoin[{imageFolder, shape}];
	If[!DirectoryQ[gifFolder], CreateDirectory[gifFolder]];

	aValues = Table[a, {a, amax, amin, -(amax - amin)/(frameCount-1) }];
	Print[aValues];

	For[aid = 1, aid <= Length[aValues], aid++,
		aValue = aValues[[aid]];
		Print[{aid, aValue}];
		exportFrame[aValue, aid, outputFolderName];
	];
	
	delay = If[frameCount >= 100, 3, Floor[300/frameCount]];

	gifCommand = "magick convert -quiet +repage -alpha remove -loop 0 -delay " <> IntegerString[delay] <> " -dispose 2 " <> outputFolder <> "\\spinning_frame*.png " <> gifFolder <> "\\" <> gifFileName <> ".gif";

	Print[gifCommand];
	Run[gifCommand];
	Print["Saved the animation to " <> gifFolder <> "\\" <> gifFileName <> ".gif"];
	
  mp4Command = "ffmpeg -i " <> gifFolder <> "\\" <> gifFileName <> ".gif" 
      <> " -movflags faststart -pix_fmt yuv420p -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" "
      <> gifFolder <> "\\" <> gifFileName <> ".mp4";

  mp4FileName = gifFolder <> "\\" <> gifFileName <> ".mp4";
  If[FileExistsQ[mp4FileName], DeleteFile[mp4FileName]];
  Print[mp4Command];
  Run[mp4Command];
  Print["Saved the video to " <> gifFolder <> "\\" <> gifFileName <> ".mp4"];
];

exportAnimation[];

(* exportFrame[-0.14, 0, "test"]; *)



(* 

for ($num = 69 ; $num -le 99 ; $num++) { cp spinning_frame_0068.png "spinning_frame_00$num.png" }

 *)