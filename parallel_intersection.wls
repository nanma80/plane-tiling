#!/usr/bin/env wolframscript

(* nGenerations = 10; *)
nGenerations = 7;
(* nGenerations = 12; *)

framesPerStage = 100;

stillStages = 0;

exportToPov = True;
highlightCenter = False;
topView = True;
topView = False;

projectToBall = True;
(* projectToBall = False; *)

trackVertices = False;

verticesOnly = True;
(* verticesOnly = False; *)

imageSize = {4, 3} * 720 / 3;

stepsPerTurn = If[trackVertices, 35, 35];
(* stepsPerTurn = 3; *)

argv = Rest @ $ScriptCommandLine;
If[Length[argv] >= 2, 
	p = ToExpression[argv[[1]]]; q = ToExpression[argv[[2]]],
	p = 6; q = 3; (* default *)
];

radius = 2.5;

maxStellationLevel = If[Length[argv] >= 3, ToExpression[argv[[3]]], 1];
maxRatio = Tan[maxStellationLevel Pi/p]/Tan[Pi/p];

dihedralSpherical = 0;
dihedralHyperbolic = 0;

shape = "para_int_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm];
colors = {Green, Blue, Red};

splitEdgeParts = 8;

imageFolder = "output";
viewPointBeforeRotation = If[topView, {1, 0, 0}, {-1, -2.5, 3.5}];
rotationLighting = Transpose[RotationMatrix[{{0.5, -1.0, 1.5}, viewPointBeforeRotation}]];

Needs["POVRayRender`"];
ConfigurePOVRayRender["POVRayPath" -> "C:\\Program Files\\POV-Ray\\v3.7\\bin\\pvengine64.exe"];

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

shapeWithParams = shape <> If[topView, "_top", ""] <> "_" <> IntegerString[nGenerations] <> If[projectToBall, "", "_no_projection"];

getGenerator[p_, q_] := Module[
	{inner01, inner02, inner12, signature, p0, p1, p2, otherp2, ch2psi},

	(* chi: p0 and p2 *)
	inner02 = Cot[Pi/p]*Cot[Pi/q];
	(* psi: p1 and p2 *)
	inner12 = Cos[Pi/q]/Sin[Pi/p];
	(* phi: p0 and p1 *)
	inner01 = Cos[Pi/p]/Sin[Pi/q];

	signature = Sign[1/2-1/p-1/q];
	(* signature = 1 for hyperbolic, -1 for spherical, 0 for planar *)

	p0 = {inner02, 
		Sqrt[signature*(inner02^2-1)], 
		0}//Simplify;
	p1 = {inner12, 
		Sqrt[signature*(inner12^2-1)] Cos[Pi/p], 
		Sqrt[signature*(inner12^2-1)] Sin[Pi/p]}//Simplify;
	(* other p2 such that the midpoint of p2 and otherp2 is p1 *)
	p2 = {1, 0, 0};

	ch2psi = 2 * inner12 ^ 2 - 1;
	Print["ch2psi"];
	Print[ch2psi//N];
	(* is 1 for {3,6} *)
	otherp2 = {ch2psi, 
		Sqrt[signature*(ch2psi^2-1)] Cos[Pi/p], 
		Sqrt[signature*(ch2psi^2-1)] Sin[Pi/p]};
	

	mirror1 = {0, 0, 1};
	mirror2 = {0, -Sin[Pi/p], Cos[Pi/p]};
	mirror3 = otherp2 - p2//Simplify;
	p0 = {1, 1, 0};
	p1 = {1, Cos[Pi/p], Sin[Pi/p]}*Cos[Pi/p];
	p2 = {1, 0, 0};
	otherp2 = 2 p1 - p2;
	mirror3 = {-Cos[Pi/p], Cos[Pi/p], Sin[Pi/p]};

	vertex1 = p0;

	Print[p0//N];
	Print[p1//N];
	Print[otherp2//N];

	mirrors = {mirror1, mirror2, mirror3}//N;
	Print["mirrors"];
	Print[mirrors//N];

	{vertex1, {mirror1, mirror2, mirror3}//N, signature}
];

generateFaces[p_, q_] := Module[
	{},

	generator = getGenerator[p, q];
	vertex1 = generator[[1]];
	mirrors = generator[[2]];
	signature = generator[[3]];

	face1 = {vertex1}//N;
	finished = False;
	For[i = 0, i < 20 && finished == False, i++,
		newVertex = FullSimplify[HDoubleReflect[face1[[-1]], mirrors[[1]], mirrors[[2]], signature]];
		If[!(ApproxSamePoint[newVertex, face1[[1]]]),
			AppendTo[face1, newVertex],
			finished = True;
		]
	];

	Print["Each face contains " <> IntegerString[Length[face1]] <> " vertices."];

	faces = {face1};
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

	scaledFaces = Map[(# / vertex1[[2]] + {1- vertex1[[1]]/vertex1[[2]], 0, 0} )&, faces, {2}];
	scaledFaces
]

findIntersection[sp11_, sp12_, sp21_, sp22_] := Module[
	{sol},
	sol = Minimize[Norm[(l1 sp11 + (1 - l1) sp12) - (l2 sp21 + (1 - l2) sp22)]^2, {l1, l2}];
	If[sol[[1]] < epsilon,
		l1 sp11 + (1 - l1) sp12 /. sol[[2]],
		{}
	]
];

getNewEdges[sp11_, sp12_, sp21_, sp22_] := Module[
	{},
	int = findIntersection[sp11, sp12, sp21, sp22];
	If[Length[int] == 0,
		{},
		{
			{int, sp11},
			{int, sp12},
			{int, sp21},
			{int, sp22}
		}
	]
];

elements3D[faces_, spinAngle_, color_] := Module[
	{allEdges, edges, allVertices, vertices, selectedFaces, displayFaces},

	rotation = RotationMatrix[spinAngle, {1, 0, 0}];	
	(* selectedFaces = Select[faces, (Total[#]/Length[#])[[1]] < 10 &]; *)
	selectedFaces = Map[(#.rotation.rotationLighting) &, faces, {2}];
	(* kleinFaces = Map[getKlein, faces, {2}]; *)
	displayFaces = Map[getHyperboloid, selectedFaces, {2}];

	(* Print["There are " <> IntegerString[Length[selectedFaces]] <> " faces."]; *)
	allEdges = Flatten[Map[getEdgesFromFace, selectedFaces], 1];
	edges = Union[allEdges , SameTest -> sameEdge];

	If[Length[selectedFaces[[1]]] == 4, 
		newEdges1 = Map[getNewEdges[#[[1]], #[[2]], #[[4]], #[[3]]]&, selectedFaces];
		newEdges2 = Map[getNewEdges[#[[2]], #[[3]], #[[1]], #[[4]]]&, selectedFaces];
		edges = Join[edges, Flatten[newEdges1, 1], Flatten[newEdges2, 1]]
	];

	If[Length[selectedFaces[[1]]] == 6, 
		newEdges1 = Map[getNewEdges[#[[1]], #[[2]], #[[5]], #[[4]]]&, selectedFaces];
		newEdges2 = Map[getNewEdges[#[[2]], #[[3]], #[[6]], #[[5]]]&, selectedFaces];
		newEdges3 = Map[getNewEdges[#[[3]], #[[4]], #[[1]], #[[6]]]&, selectedFaces];
		edges = Join[edges, Flatten[newEdges1, 1], Flatten[newEdges2, 1], Flatten[newEdges3, 1]]
	];

	(* Print["There are " <> IntegerString[Length[edges]] <> " edges in total."]; *)
	(* norm = HNorm[edges[[1]][[1]]]; *)

	(* splitEdges = Flatten[ Map[splitEdge[#, splitEdgeParts]&, edges], 1]; *)
	(* splitEdges = Map[HNormalize[#, norm]&, splitEdges, {2}]; *)
	(* Print["There are " <> IntegerString[Length[splitEdges]] <> " edges after split."]; *)
	displayEdges = Map[getHyperboloid, edges, {2}];

	allVertices = Round[Flatten[displayEdges, 1], epsilon];
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

	(* planar *)
	opacity = If[topView, 1, 1];
	lighting = {{"Point", White, {50, -50, 50}}};

	elements = {
			{
				Opacity[opacity], FaceForm[color], Map[Polygon, displayFaces]
			},
			{
				Opacity[1], FaceForm[White], Map[Cylinder[#, thickness]&, displayEdges]
			},
			{
				Opacity[1], FaceForm[Purple], Map[Sphere[#[[1]], thickness * 2]&, verticesTally]
			}
	};

	If[highlightCenter,
		highlightColor = White;
		elements = Union[
			elements,
			{
				{
					Opacity[1], FaceForm[highlightColor], Map[Polygon, displayCenteredFaces]
				},
				{
					Opacity[1], FaceForm[highlightColor], Map[Cylinder[#, thickness * 2.2]&, displayCenteredEdges]
				}
			}
		];
	];
	elements
];

exportFrame[components_, spinAngle_, frameIndex_, outputFolderName_] := Module[
	{},

	elements = Union[ Table[elements3D[components[[k]], spinAngle, colors[[k]] ], {k, 1, Length[components]}] ];

	(* rangeZoom = Sqrt[maxRatio] * 1.2 * If[p < 6, 1, 0.5]; *)
	
	(* range = 6 * rangeZoom; *)
	(* thickness = range / 80.0 / rangeZoom * If[p < 6, 1, 0.7]; *)
	viewCenter = {0, 0, 0};
	range = radius * 3;
	thickness = range * 0.008;

	graph = Graphics3D[
		elements,
		ViewPoint -> viewPointBeforeRotation.Transpose[rotationLighting] * 100,
		ViewVertical -> If[topView, {0, 1, 0}, {-1, 0, 0}].Transpose[rotationLighting],
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

	faces = generateFaces[p, q];

	(* ensure all faces are rotating to the same direction *)
	centeredFaces = Table[Map[# - Mean[faces[[k]]]&, faces[[k]]], {k, 1, Length[faces]}];
	argDiff =
		Table[
			Mod[(Arg[centeredFaces[[k]][[1]][[2]] + I centeredFaces[[k]][[1]][[3]]] - Arg[centeredFaces[[k]][[2]][[2]] + I centeredFaces[[k]][[2]][[3]]])/Pi, 2]
			, {k, 1, Length[faces]} ];

	For[k = 1, k <= Length[faces], k++,
		If[argDiff[[k]] < 1,
			faces[[k]] = Reverse[faces[[k]]];
		];
	];

	(* move a bit to avoid singularity *)
	faces = Map[(# + {0, Sqrt[7]/100., Sqrt[11]/100. })&, faces, {2}];

	stageCount = 1;
	For[stage = 1, stage <= stageCount * (1 + stillStages), stage++,
		For[srId = 1, srId <= framesPerStage, srId ++,
			sr = srId/framesPerStage;
			Print["sr = " <> ToString[N[sr]]];
			zOffset = If[{p, q}=={4,4}, 2, Sqrt[3]];
			movedFaces = Map[(# + {0, 0, -zOffset Mod[sr, 1] })&, faces, {2}];

			components = {movedFaces};

			componentsStereo = If[projectToBall, Map[radius stereo[# / radius]&, components, {3}], components];

			frameIndex = (stage - 1) * framesPerStage + srId;
			exportFrame[componentsStereo, 0, frameIndex, outputFolderName];
		];
	];

	frameCount = frameIndex;
	delay = If[frameCount >= 50, 3, Floor[500/frameCount]];

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

(* centers = Map[Mean, faces]; *)

(* xCenters = Select[centers, Abs[#[[2]]] < epsilon || Abs[#[[3]]] < epsilon & ];
Print[xCenters];
 *)
(* Exit[]; *)

(* component1 = Map[(# + {0, 0, 0.05})&, component1, {2}]; *)


exportAnimation[];





(* 

for ($num = 69 ; $num -le 99 ; $num++) { cp spinning_frame_0068.png "spinning_frame_00$num.png" }

 *)