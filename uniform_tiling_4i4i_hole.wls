#!/usr/bin/env wolframscript

(* 
Objective: plot the non convex tilings on
https://en.wikipedia.org/wiki/Uniform_tiling#Expanded_lists_of_uniform_tilings
 *)


(* nGenerations = 13; *)
(* nGenerations = 10; *)
nGenerations = 20;
(* nGenerations = 1; *)

framesPerStage = 10;

stillStages = 0;

holeCenter2D = {1, -1, 0};

exportToPov = True;
highlightCenter = False;

projectToBall = True;
(* projectToBall = False; *)

topView = True;
(* topView = False; *)
If[projectToBall == False, topView = True];

trackVertices = False;

verticesOnly = True;
(* verticesOnly = False; *)

imageSize = {4, 3} * 720 / 3;

stepsPerTurn = If[trackVertices, 35, 35];
(* stepsPerTurn = 3; *)

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 2, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]],
	p = 4; q = 4; (* default *)
];

(* radius = 2.5; *)
radius = 5;

maxStellationLevel = If[Length[argv] >= 3, ToExpression[argv[[3]]], 1];
maxRatio = Tan[maxStellationLevel Pi/p]/Tan[Pi/p];

dihedralSpherical = 0;
dihedralHyperbolic = 0;

shape = "uniform_4i4i_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm];

colors = {Yellow, Blue, Green, Cyan,  Red, Purple  };

splitEdgeParts = 8;

imageFolder = "output";

viewPointBeforeRotation = If[topView, {1, 0, 0}, {-3, -2.5, 3.5}];
(* rotationLighting = Transpose[RotationMatrix[{{0.8, -1.0, 1.2}, viewPointBeforeRotation}]]; *)
(* rotationLighting = Transpose[RotationMatrix[{{0.8, -1.0, 1.2}, viewPointBeforeRotation}]]; *)
rotationLighting = RotationMatrix[- Pi/4, {0, 1, 1}];
viewVerticalBeforeRotation = viewPointBeforeRotation.RotationMatrix[Pi/2, {0, 1, 1}];

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

diagonal1[v_] := Round[v[[2]] + v[[3]], epsilon];
diagonal2[v_] := Round[v[[2]] - v[[3]], epsilon];

centerMod2[face_] := Mod[ Round[Mean[face]][[3]], 2];
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
splitFaceToQuads[face_] := Module[
	{pieces},
	pieces = {};
	For[id = 2, id <= Length[face] - 1, id = id + 2,
		AppendTo[pieces, {face[[1]], face[[id]], face[[id + 1]], face[[Min[id + 2, Length[face]]]]}]
	];
	pieces
];
getKlein[v_] := {0, v[[2]], v[[3]]}/v[[1]];
getPoincare[v_] := {0, v[[2]], v[[3]]}/(1+v[[1]]);
splitEdge[edge_, n_] := Table[{k edge[[1]] + (1-k) edge[[2]], (k+1/n) edge[[1]] + (1-k-1/n) edge[[2]] }, {k, 0, 1 - 1/n, 1/n}];
HNormalize[v_, norm_] := v/HNorm[v] * norm;
getHyperboloid[v_] := {-v[[1]], v[[2]], v[[3]]};

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

stereo[p_] := {
	(-1 + p[[2]]^2 + p[[3]]^2)/(1 +	p[[2]]^2 + p[[3]]^2),
	2 p[[2]]/(1 + p[[2]]^2 + p[[3]]^2), 
	2 p[[3]]/(1 + p[[2]]^2 + p[[3]]^2)	
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

shapeWithParams = shape <> "_hole_" <> IntegerString[nGenerations];

getGenerator[p_, q_] := Module[
	{signature, p0, p1, p2, otherp2, ch2psi},

	mirror1 = {0, -Sin[Pi/p], Cos[Pi/p]};
	mirror2 = {0, 0, 1};
	mirror3 = {-Cos[Pi/p], Cos[Pi/p], Sin[Pi/p]};

	vertex1 = {1, 1, 0};

	mirrors = {mirror1, mirror2, mirror3};
	Print["mirrors"];
	Print[mirrors//N];

	{vertex1, mirrors, 0}
];

generateFaces[p_, q_, option_] := Module[
	{sol, x, y, generator, vertex1, mirrors, signature, image1, image2, image3},

	Print["GenerateFaces"];
	Print[{p,q,option}];
	generator = getGenerator[p, q];
	(* vertex1 = generator[[1]]; *)
	mirrors = generator[[2]];
	signature = generator[[3]];

	vertex1 = {1, x, y};
	image1 = HReflect[vertex1, mirrors[[1]], signature];
	image2 = HReflect[vertex1, mirrors[[2]], signature];
	image3 = HReflect[vertex1, mirrors[[3]], signature];

	sol = {};
	(* 4,8,8 https://en.wikipedia.org/wiki/Truncated_square_tiling *)
	If[option == 1 || option == 3,
		sol = Solve[Norm[image1-vertex1]==Norm[image2 - vertex1] && Norm[vertex1 - image3] == 0 && x < 1 && x >= 0 && y >= 0, {x, y}];
	];

	(* large octagon. Related to 3,4,6,4 but not sure what exactly this is *)
	If[option == 2,
		sol = {{x -> 1 + Sqrt[2]/2, y -> Sqrt[2]/2}};
	];

	Print[sol];
	
	vertex1 = If[option == 0, generator[[1]], vertex1/.sol[[1]]];

	Print["vertex1"];
	Print[vertex1//N];

	(* Exit[]; *)

	faces = {};
	
	For[mirror1Index = 1, mirror1Index <= Length[mirrors], mirror1Index++,
		mirror2Index = mirror1Index + 1;
		If[mirror2Index > Length[mirrors], mirror2Index = 1];

		face1 = {vertex1}//N;
		finished = False;
		For[i = 0, i < 20 && finished == False, i++,
		  beforeCount = Length[face1];
		  newVertices1 = Map[HReflect[#, mirrors[[mirror1Index]], signature] &, face1];
		  newVertices2 = Map[HReflect[#, mirrors[[mirror2Index]], signature] &, face1];
		  face1 = Union[face1, newVertices1, newVertices2, SameTest -> ApproxSamePoint];
		  Print[{i, Length[face1]}];
		  If[beforeCount == Length[face1], finished = True]
		];

		faceCenter = Mean[face1];
		proj = Orthogonalize[{{1, 2, 3}, {1.23, 4.34, 0}}];
		p2d = Map[(# - faceCenter) . Transpose[proj] &, face1];
		args = Map[Arg[#[[1]] + #[[2]] I] &, p2d];
		face1 = face1[[Ordering[args]]];

		If[option == 3 && mirror1Index == 1,
			face1 = face1[[ Table[Mod[3 k + 1, Length[face1]] + 1, {k, 1, Length[face1]}] ]];
		];

		Print["Each face contains " <> IntegerString[Length[face1]] <> " vertices. mirror1Index = " <> IntegerString[mirror1Index]];

		If[Length[face1] > 2,
			AppendTo[faces, face1]
		];	
	];


	(* faces = {face1}; *)
	countStart = -1;
	For[i = 0, i <= nGenerations && countStart != Length[faces], i++,
		(* Print["i = " <> IntegerString[i]]; *)
		countStart = Length[faces];
		For[mirrorIndex = 3, mirrorIndex >= 1, mirrorIndex--,
			newFaces = Map[HReflect[#, mirrors[[mirrorIndex]] , signature ]&, faces, {2}];
			faces = Union[faces, newFaces, SameTest -> sameCenter];
		];
	];

	(* Print["Making it more symmetric"]; *)

	(* Making it symmetric *)
	centeredMirrors = Select[mirrors, #[[1]] == 0 &];

	countStart = -1;
	For[i = 0, i < 20 && countStart != Length[faces], i++,
		(* Print["i = " <> IntegerString[i]]; *)
		countStart = Length[faces];
		For[mirrorIndex = Length[centeredMirrors], mirrorIndex >= 1, mirrorIndex--,
			newFaces = Map[HReflect[#, centeredMirrors[[mirrorIndex]], signature ]&, faces, {2}];
			faces = Union[faces, newFaces, SameTest -> sameCenter];
		];
	];

	faces
];

elements3D[faces_, spinAngle_,  stage_, sr_, color_] := Module[
	{allEdges, edges, allVertices, vertices, selectedFaces, displayFaces},

	rotation = RotationMatrix[spinAngle, {1, 0, 0}];	
	(* selectedFaces = Select[faces, (Total[#]/Length[#])[[1]] < 10 &]; *)
	selectedFaces = Map[(#.rotation.rotationLighting) &, faces, {2}];
	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	holeCenter = holeCenter3D.rotation.rotationLighting;
	displayHoleCenter = N[getHyperboloid[holeCenter]];

	selectedFacesNoHole = Select[selectedFaces,
		Length[Intersection[#, {holeCenter}, SameTest -> ApproxSamePoint]] == 0 &
	];
	allFacesAfterTruncation = selectedFacesNoHole;

	selectedFacesWithHole = Select[selectedFaces,
		Length[Intersection[#, {holeCenter}, SameTest -> ApproxSamePoint]] > 0 &
	];

	truncatedEdges = {};

	For[fid = 1, fid <= Length[selectedFacesWithHole], fid++,
		faceWithHoles = selectedFacesWithHole[[fid]];
		vid = Select[Range[Length[faceWithHoles]], ApproxSamePoint[faceWithHoles[[#]], holeCenter]&][[1]];
		faceWithHoles = RotateLeft[faceWithHoles, vid - 1];
		truncatedFaceWithHoles = Join[
			{(faceWithHoles[[1]] + faceWithHoles[[ Length[faceWithHoles] ]])/2},
			{(faceWithHoles[[1]] + faceWithHoles[[2]])/2},
			faceWithHoles[[Range[2, Length[faceWithHoles] ]]]
		];
		AppendTo[allFacesAfterTruncation, truncatedFaceWithHoles];
		AppendTo[truncatedEdges, 
			{
				(faceWithHoles[[1]] + faceWithHoles[[ Length[faceWithHoles] ]])/2,
				(faceWithHoles[[1]] + faceWithHoles[[2]])/2
			}
		];
	];

	displayFaces = Map[getHyperboloid, allFacesAfterTruncation, {2}];
	allEdges = Flatten[Map[getEdgesFromFace, allFacesAfterTruncation], 1];
	edges = Union[allEdges , SameTest -> sameEdge];

	displayEdges = Map[getHyperboloid, edges, {2}];
	displayTruncatedEdges = Map[getHyperboloid, truncatedEdges, {2}];

	allVertices = Round[Flatten[displayEdges, 1], epsilon];
	truncatedVertices = Round[Flatten[displayTruncatedEdges, 1], epsilon];

	vertices = Union[allVertices , SameTest -> ApproxSamePoint];
	verticesTally = Tally[allVertices];
	verticesTallyCounts = Map[#[[2]]&, verticesTally];
	tallyTally = Tally[verticesTallyCounts];
	(* Print["There are " <> IntegerString[Length[vertices]] <> " vertices in total."]; *)
	Print[Map[Length, {faces, allEdges, edges, allVertices, vertices}]];
	Print[tallyTally];

	centeredFaces = Select[faces, (Abs[Total[#][[2]]]<0.001)&&(Abs[Total[#][[3]]]<0.001)&];
	displayCenteredFaces = Map[getHyperboloid, centeredFaces, {2}];
	displayCenteredFaces = Map[(# + {0.001, 0, 0})&, displayCenteredFaces, {2}];

	centeredEdges = Flatten[Map[getEdgesFromFace, centeredFaces], 1];
	
	(* splitCenteredEdges = Map[HNormalize[#, norm]&, splitCenteredEdges, {2}]; *)
	(* Print["There are " <> IntegerString[Length[centeredEdges]] <> " centered edges after split."]; *)
	displayCenteredEdges = Map[getHyperboloid, centeredEdges, {2}];

	(* centeredVertices = displayCenteredFaces[[1]]; *)

	opacity = 1;
	lighting = {{"Point", White, {50, -50, 50}}};

	range = radius * 0.8;
	thickness = range * 0.008;

	elements = {
			{
				Opacity[opacity], FaceForm[color], Map[Polygon, displayFaces]
			},
			{
				Opacity[1], FaceForm[White], Map[Cylinder[#, thickness]&, displayEdges]
			},
			{
				Opacity[1], FaceForm[White], Map[Sphere[#[[1]], thickness]&, verticesTally]
			},
			{
				Opacity[1], FaceForm[Red], Map[Cylinder[#, thickness * 1.1]&, displayTruncatedEdges]
			},
			{
				Opacity[1], FaceForm[Red], Map[Sphere[#[[1]], thickness * 1.1]&, truncatedVertices]
			}
	};

	elements
];

exportFrame[components_, spinAngle_, stage_, sr_, frameIndex_, outputFolderName_] := Module[
	{},
	(* Return[]; *)
	facesWithHoles = {};
	elements = Union[ Table[elements3D[components[[k]], spinAngle, stage, sr, colors[[k]] ], {k, 1, Length[components]}] ];

	(* rangeZoom = Sqrt[maxRatio] * 1.2 * If[p < 6, 1, 0.5]; *)
	
	(* range = 6 * rangeZoom; *)
	(* thickness = range / 80.0 / rangeZoom * If[p < 6, 1, 0.7]; *)
	viewCenter = {0, 0, 0};

	graph = Graphics3D[
		elements,
		ViewPoint -> viewPointBeforeRotation.Transpose[rotationLighting] * 100,
		ViewVertical -> viewVerticalBeforeRotation.Transpose[rotationLighting],
		PlotRange -> {{viewCenter[[1]]-range, viewCenter[[1]]+range}, {viewCenter[[2]]-range, viewCenter[[2]]+range}, {viewCenter[[3]]-range, viewCenter[[3]]+range}},
		Lighting -> lighting
	];

	If[Length[facesWithHoles] > 1,
		(* graph = Apply[Show, Join[{graph}, facesWithHoles]]; *)
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

	faces1 = generateFaces[4, 4, 0];
	groupedFaces = GroupBy[faces1, Length];

	faces4half = {Select[groupedFaces[4], centerMod2[#] == 0 &]};
	faces4halfCenters = Map[Mean, faces4half[[1]]];
	faces4halfCenters = Select[faces4halfCenters, Abs[ #[[3]] ]<epsilon &];
	(* Print[faces4halfCenters]; *)
	(* Exit[]; *)

 	vertices = Union[Flatten[faces4half[[1]], 1], SameTest -> ApproxSamePoint];
 	Print[Length[vertices]];

 	degenFaces1 = Values[GroupBy[vertices, diagonal1]];
 	degenFaces1 = Select[degenFaces1, Length[#]>2 &];
 	degenFaces2 = Values[GroupBy[vertices, diagonal2]];
 	degenFaces2 = Select[degenFaces2, Length[#]>2 &];

	components = Join[faces4half, {degenFaces1}, {degenFaces2}];

	stageCount = 1;

	For[stage = 1, stage <= stageCount * (1 + stillStages), stage++,
		For[srId = If[stage==1, 0, 1], srId <= framesPerStage, srId ++,
		(* For[srId = 3, srId <= 3, srId ++, *)
			sr = srId/framesPerStage;
			Print["sr = " <> ToString[N[sr]]];
			yOffset = 2 Mod[sr, 1];
			movedComponents = Map[(# + {0, yOffset, 0})&, components, {3}];

			truncateStage = Ceiling[stage/(1 + stillStages)];
			truncateSr = If[truncateStage * (1 + stillStages) == stage, sr, 0];

			componentsStereo = If[projectToBall, Map[radius stereo[# / radius]&, movedComponents, {3}], movedComponents];

			movedHoleCenter2D = holeCenter2D + {0, yOffset, 0};
			holeCenter3D = radius stereo[movedHoleCenter2D / radius];

			frameIndex = (stage - 1) * framesPerStage + srId;
			exportFrame[componentsStereo, -Pi/4, truncateStage, truncateSr, frameIndex, outputFolderName];
		];
	];

	frameCount = frameIndex;
	delay = If[frameCount >= 100, 3, Floor[300/framesPerStage]];

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

(* exportFrame[componentsStereo, 0, 1, "test"]; *)



(* 

for ($num = 69 ; $num -le 99 ; $num++) { cp spinning_frame_0068.png "spinning_frame_00$num.png" }

 *)