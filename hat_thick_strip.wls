#!/usr/bin/env wolframscript

(* 
Objective: plot the non convex tilings on
https://en.wikipedia.org/wiki/Uniform_tiling#Expanded_lists_of_uniform_tilings
 *)


(* nGenerations = 13; *)
(* nGenerations = 20; *)
nGenerations = 30;
(* nGenerations = 1; *)
(* nGenerations = 7; *)

frameCount = 4;

exportToPov = True;
highlightCenter = False;

(* projectToBall = True; *)
projectToBall = False;

topView = True;

trackVertices = False;

verticesOnly = True;
(* verticesOnly = False; *)

imageSize = {4, 3} * 720 / 3;

stepsPerTurn = If[trackVertices, 35, 35];
(* stepsPerTurn = 3; *)

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 2, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]],
	p = 3; q = 6; (* default *)
];

(* radius = 15; *)
radius = 5;

periodStart = -3; 
periodEnd = 1;

periodStart = 0; 
periodEnd = 0;

maxStellationLevel = If[Length[argv] >= 3, ToExpression[argv[[3]]], 1];
maxRatio = Tan[maxStellationLevel Pi/p]/Tan[Pi/p];

dihedralSpherical = 0;
dihedralHyperbolic = 0;

shape = "hat_thick_strip";
opacity = 1;
lightingOrigin = {50, -50, 50};
lighting = {{"Point", White, lightingOrigin}};
colors = {Blue, Green, Cyan, Yellow, Red, Purple, Brown, Orange, Lighter[Red, 0.5] };

splitEdgeParts = 8;

imageFolder = "output";

viewPointBeforeRotation = If[topView, {0, 0, 1}, {-3, -2.5, 3.5}];
(* rotationLighting = Transpose[RotationMatrix[{{0.8, -1.0, 1.2}, viewPointBeforeRotation}]]; *)
(* rotationLighting = Transpose[RotationMatrix[{{0.8, -1.0, 1.2}, viewPointBeforeRotation}]]; *)
(* rotationLighting = Transpose[RotationMatrix[{ {0, 0, 1}, lightingOrigin}]]; *)
rotationLighting = IdentityMatrix[3];
(* rotationLighting = RotationMatrix[Pi/4, {0, 1, 1}]; *)

(* viewVerticalBeforeRotation = viewPointBeforeRotation.RotationMatrix[Pi/2, {0, 1, 1}]; *)
viewVerticalBeforeRotation = {0, 1, 0};

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

diagonal1[v_] := Round[v[[2]] Sqrt[3] + v[[3]], epsilon];
diagonal2[v_] := Round[v[[2]] Sqrt[3] - v[[3]], epsilon];
diagonal3[v_] := Round[v[[3]], epsilon];
HInner[v_, u_, sig_] := If[sig > 0, 2 * v[[1]] * u[[1]] - Dot[v, u], Dot[v, u]];
HNormSquare[v_, sig_] := HInner[v, v, sig];
HNorm[v_, sig_] := Sqrt[HNormSquare[v, sig]];
rotated[vertices_, center_, axis_, angle_] := Map[# + center &, Map[# - center &, vertices].RotationMatrix[angle, axis]];
Rotation[t_] :=
{
  {1,  0, 0},
  {0, Cos[t], -Sin[t]},
  {0, Sin[t], Cos[t]}
};
epsilon = 0.00000001;
HReflect[point_, mirror_, sig_] := If[sig!= 0, 
	point - 2 * HInner[point, mirror, sig]/HInner[mirror, mirror, sig]*mirror,
	point - 2 * (mirror[[2]]*point[[2]] + mirror[[3]]*point[[3]] + mirror[[1]])/(mirror[[2]]^2 + mirror[[3]]^2) * {0, mirror[[2]], mirror[[3]]}
	];
HDoubleReflect[point_, mirror1_, mirror2_, sig_] := HReflect[HReflect[point, mirror1, sig], mirror2, sig];
(* ApproxSamePoint[point1_, point2_] := Norm[N[point1 - point2]] < 0.00001; *)
ApproxSamePoint[point1_, point2_] := Round[point1, epsilon] == Round[point2, epsilon];

getEdgesFromFace[face_] := Table[{face[[i+1]], face[[Mod[i+1, Length[face]] + 1]]}, {i, 0, Length[face] - 1}];
sameCenter[set1_, set2_] := ApproxSamePoint[Total[set1], Total[set2]];
(* sameEdge[e1_, e2_] := sameCenter[e1, e2] && ((ApproxSamePoint[e1[[1]], e2[[1]]] && ApproxSamePoint[e1[[2]], e2[[2]]] ) || (ApproxSamePoint[e1[[1]], e2[[2]]] && ApproxSamePoint[e1[[2]], e2[[1]]] )); *)
sameEdge[e1_, e2_] := sameCenter[e1, e2];

getKlein[v_] := {0, v[[2]], v[[3]]}/v[[1]];
getPoincare[v_] := {0, v[[2]], v[[3]]}/(1+v[[1]]);
splitEdge[edge_, n_] := Table[{k edge[[1]] + (1-k) edge[[2]], (k+1/n) edge[[1]] + (1-k-1/n) edge[[2]] }, {k, 0, 1 - 1/n, 1/n}];
HNormalize[v_, norm_] := v/HNorm[v] * norm;
getHyperboloid[v_] := v;

scale[edge_, ratio_] := Module[
	{center},
	center = Mean[edge];
	Map[(# - center) * ratio + center&, edge]
];

trimFaces[faces_] := Module[
	{},
	allEdges = Flatten[Map[getEdgesFromFace, faces], 1];
	allEdgesTally = Tally[allEdges, ApproxSamePoint[Total[#1], Total[#2]]&];
	commonEdges = Map[#[[1]]&, Select[allEdgesTally, #[[2]] > 1 &]];
	trimmedFaces = Select[faces, Length[Intersection[getEdgesFromFace[#], commonEdges, SameTest -> sameCenter]] > 1 & ];
	trimmedFaces
];

stereo3[p_] := {
	2 p[[1]]/(1 + p[[1]]^2 + p[[2]]^2), 
	2 p[[2]]/(1 + p[[1]]^2 + p[[2]]^2),
	(1 - p[[1]]^2 - p[[2]]^2)/(1 +	p[[1]]^2 + p[[2]]^2)
};

scaleFace[face_, ratio_, stellationLevel_] := Module[
	{edges, edgeOrders},
	edges = Map[scale[#, ratio]&, getEdgesFromFace[face]];
	edgeOrders = Table[Mod[stellationLevel * k, p]+1, {k, 0, p-1}];
	(* edges[[edgeOrders]] *)
	(* Map[#[[1]]&, edges[[edgeOrders]]] *)
	Flatten[edges[[edgeOrders]], 1]
];

boost[t_] :=
{
  {Cosh[t], 0, Sinh[t]},
  {0, 1, 0},
  {Sinh[t], 0, Cosh[t]}
};

(* only works if axis[[1]] == 0 *)
boostMatrix[angle_, axis_] := RotationMatrix[Arg[axis[[2]] + I axis[[3]]], {1, 0, 0}].boost[angle].RotationMatrix[-Arg[axis[[2]] + I axis[[3]]], {1, 0, 0}];

boosted[vertices_, center_, axis_, angle_] := Map[# + center &, Map[# - center &, vertices].boostMatrix[angle, axis]];

shapeWithParams = shape <> "_" <> IntegerString[frameCount] <> If[projectToBall, "", "_no_projection"];

elements3D[faces_] := Module[
	{allEdges, edges, allVertices, vertices, selectedFaces, displayFaces},

	selectedFaces = Map[(#.rotationLighting) &, faces, {2}];
	displayFaces = Map[getHyperboloid, selectedFaces, {2}];

	allEdges = Flatten[Map[getEdgesFromFace, selectedFaces], 1];
	edges = Union[allEdges , SameTest -> sameEdge];

	displayEdges = Map[getHyperboloid, edges, {2}];

	allVertices = Round[Flatten[displayEdges, 1], epsilon];
	vertices = Union[allVertices , SameTest -> ApproxSamePoint];
	verticesTally = Tally[allVertices];

	elements = {
			Table[
				{
					Opacity[opacity], FaceForm[ allColors[[k]] ], Polygon[displayFaces[[k]]]
				}, 
				{k, Length[displayFaces]}
			],
			{
				Opacity[opacity], FaceForm[White], Map[Cylinder[#, thickness*opacity]&, displayEdges]
			},
			{
				Opacity[opacity], FaceForm[White], Map[Sphere[#[[1]], thickness*opacity]&, verticesTally]
			}
	};

	elements
];

hatPrototype =
	{
		{0,0},
		{0,Sqrt[3]},
		{1,Sqrt[3]},
		{3/2,(3 Sqrt[3])/2},
		{3,Sqrt[3]},
		{3,0},
		{4,0},
		{9/2,-(Sqrt[3]/2)},
		{3,-Sqrt[3]},
		{3/2,-(Sqrt[3]/2)},
		{1,-Sqrt[3]},
		{-1,-Sqrt[3]},
		{-(3/2),-(Sqrt[3]/2)}
	};

getTriangle[hat_] := hat[[ {1, 5, 9} ]];
getCenter[hat_] := Mean[getTriangle[hat]];

hatPrototype = Map[# - getCenter[hatPrototype] &, hatPrototype];
hatTriangle = getTriangle[hatPrototype];

hatTriangle = getTriangle[hatPrototype];
triangleCenterVectors = Reverse[Sort[2 * Map[Mean, getEdgesFromFace[hatTriangle]]]];

(* 
	Three transforms in sequence:
	1. flip. If True, flip around the y axis. x -> -x. y unchanged
	2. rotation, in degrees. Positive: CW
	3. offset, 3 dimensions. {right, up left, down left}
*)

getHat[flip_, deg_, offset_] := Map[(# + offset.triangleCenterVectors) &, hatPrototype.{{If[flip, -1, 1], 0}, {0, 1}}.RotationMatrix[deg/180 * Pi]]//N;

getOffsetShapes[shapes_, offset_] := Map[ (# + offset.triangleCenterVectors)&, shapes, {2} ];
hatSignature[hat_] := Round[Simplify[ {(hat[[1]] - getCenter[hat]), (hat[[2]] - getCenter[hat])} ], epsilon];

period = {7, 0, 2};
periodOffset = period.triangleCenterVectors;
rotationCenter = Mean[getHat[True, 60, {0, 0, 0}][[ {8, 13} ]]];

getShapes3D[periodStart_, periodEnd_, periodicOffsetRatio_] := Module[
	{centerRow, upperHalf, lowerHalf, fundamentalDomain, signatures, colorMap, shapes2D, rotatedShapes2D},

	centerRow = {
		getHat[True, 60, {0, 0, 0}]
		, getHat[False, 120, {2, 0, 1}]
		, getHat[False, 120, {4, 0, 2}]
		, getHat[False, 120, {6, 0, 3}]
	};

	upperHalf = {
		(* blue row *)
		getHat[False, -120, {1, 2, 0}]
		, getHat[False, 60, {3, 1, 0}]
		, getHat[False, 180, {4, 0, 0}]

		(* green row *)
		, getHat[False, 120, {3, 3, 0}]
		, getHat[False, 60, {5, 2, 0}]
		, getHat[False, 0, {7, 2, 0}]

		(* purple row *)
		, getHat[False, 180, {5, 5, 0}]
		, getHat[True, 0, {6, 4, 0}]
		, getHat[False, 60, {8, 5, 0}]
		, getHat[False, -60, {9, 4, 0}]

		(* yellow row *)
		, getHat[False, 180, {6, 7, 0}]
		, getHat[False, 120, {8, 7, 0}]
		, getHat[False, 60, {10, 6, 0}]

		(* cyan row *)
		, getHat[False, 180, {7, 9, 0}]
	};

	lowerHalf = getOffsetShapes[ Map[( (2 rotationCenter - #) )&, upperHalf, {2}], period];

	fundamentalDomain = Join[ centerRow, upperHalf, lowerHalf ];
	signatures = Sort[Map[#[[1]]&, Tally[Map[hatSignature, fundamentalDomain ] ]]];
	colorMap = <||>;

	For[sid = 1, sid<=Length[signatures], sid++,
		colorMap[signatures[[sid]] ] = colors[[sid]]
	];

	(* shapes2D = fundamentalDomain; *)
	shapes2D = Apply[Join, Table[getOffsetShapes[fundamentalDomain, k * period], {k, periodStart, periodEnd}]];

	allColors = Map[( colorMap[hatSignature[#]] )&, shapes2D];

	shiftedShapes2D = getOffsetShapes[shapes2D, periodicOffsetRatio * period];
	rotatedShapes2D = Map[#.RotationMatrix[{{1, 0}, periodOffset}]&, shiftedShapes2D];

	Map[{#[[1]], #[[2]], 0}&, rotatedShapes2D, {2}]
];

exportFrame[faces_, frameIndex_, outputFolderName_] := Module[
	{},
	(* If[frameIndex < 100, Return[]]; *)
	(* If[Mod[frameIndex, 2] == 0, Return[]]; *)
	(* Return[]; *)

	elements = elements3D[faces];

	viewCenter = {0, 0, 0};
	range = radius * 0.8;
	thickness = range * 0.003;

	graph = Graphics3D[
		elements,
		ViewPoint -> viewPointBeforeRotation.Transpose[rotationLighting] * 100,
		ViewVertical -> viewVerticalBeforeRotation.Transpose[rotationLighting],
		PlotRange -> {{viewCenter[[1]]-range, viewCenter[[1]]+range}, {viewCenter[[2]]-range, viewCenter[[2]]+range}, {viewCenter[[3]]-range, viewCenter[[3]]+range}},
		Lighting -> lighting
	];

	outputFolder = FileNameJoin[{imageFolder, outputFolderName}];
	If[!DirectoryQ[outputFolder], CreateDirectory[outputFolder]];

	outputFileName = "spinning_frame_" <> IntegerString[frameIndex, 10, 4] <> ".png";
	outputFileName = FileNameJoin[{outputFolder, outputFileName}];

	cropCommand = "mogrify -crop " <> IntegerString[imageSize[[2]]] <> "x" <> IntegerString[imageSize[[2]]] <> "+" <> IntegerString[(imageSize[[1]]-imageSize[[2]])/2] <> "+0 " <> outputFileName;

	If[exportToPov,
	  If[FileExistsQ[outputFileName], DeleteFile[outputFileName]];
	  POVRayRender[graph, Method -> "Triangle", OutputPath -> outputFileName, ImageSize -> imageSize, RenderOptions -> "+A0.001 -J +UA"];
	  Run[cropCommand];
	  ,
	  Export[outputFileName, graph];
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


	For[frameIndex = 1, frameIndex <= frameCount, frameIndex++,
		Print["frameIndex/frameCount: " <> IntegerString[frameIndex] <> " / " <> IntegerString[frameCount] ];
		shapes3D = getShapes3D[periodStart, periodEnd, frameIndex/frameCount];
		
		If[projectToBall,
			shapes3D = Map[radius stereo3[# / radius]&, shapes3D, {2}]
		];

		exportFrame[shapes3D, frameIndex, outputFolderName];
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

(* exportAnimation[]; *)

(* For[fid=0, fid<=10, fid++,
	exportFrame[getShapes3D[-2, 1, fid/10.0], 100+fid, "test"];
];
 *)

exportFrame[getShapes3D[0, 0, 0], 100, "test"];

(* 

for ($num = 69 ; $num -le 99 ; $num++) { cp spinning_frame_0068.png "spinning_frame_00$num.png" }

 *)