#!/usr/bin/env wolframscript

(* 
Objective: plot the non convex tilings on
https://en.wikipedia.org/wiki/Uniform_tiling#Expanded_lists_of_uniform_tilings
 *)
(* nGenerations = 10; *)
nGenerations = 7;
(* nGenerations = 13; *)
(* nGenerations = 0; *)

framesPerStage = 1;

stillStages = 1;

exportToPov = True;
highlightCenter = False;
topView = True;
(* topView = False; *)

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

radius = 4;

maxStellationLevel = If[Length[argv] >= 3, ToExpression[argv[[3]]], 1];
maxRatio = Tan[maxStellationLevel Pi/p]/Tan[Pi/p];

dihedralSpherical = 0;
dihedralHyperbolic = 0;

shape = "uniform_" <> ToString[p, InputForm] <> "_" <> ToString[q, InputForm];
colors = {Yellow, Blue, Red};

splitEdgeParts = 8;

imageFolder = "output";
viewPointBeforeRotation = If[topView, {1, 0, 0}, {-3, -2.5, 3.5}];
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

generateFaces[p_, q_] := Module[
	{},

	generator = getGenerator[p, q];
	vertex1 = generator[[1]];
	mirrors = generator[[2]];
	signature = generator[[3]];

	vertex1 = {1, x, y};
	image1 = HReflect[vertex1, mirrors[[1]], signature];
	image2 = HReflect[vertex1, mirrors[[2]], signature];
	image3 = HReflect[vertex1, mirrors[[3]], signature];
	
	(* 3,4,6,4 https://en.wikipedia.org/wiki/Rhombitrihexagonal_tiling *)
	(* sol = Solve[Norm[image1-vertex1]==Norm[image3 - vertex1] && Norm[vertex1 - image2] == 0 && x < 1 && x >= 0 && y >= 0, {x, y}]; *)

	(* 4,6,12 https://en.wikipedia.org/wiki/Truncated_trihexagonal_tiling *)
	(* sol = Solve[Norm[image1-vertex1]==Norm[image3 - vertex1] && Norm[image1-vertex1]==Norm[image2 - vertex1] && x < 1 && x >= 0 && y >= 0, {x, y}]; *)

	(* Gruenbaum 18, 19, 17? https://en.wikipedia.org/wiki/Uniform_tiling#Expanded_lists_of_uniform_tilings *)
	(* sol = Solve[Norm[image1-vertex1]==Norm[image3 - vertex1] && Norm[vertex1 - image2] == 0 && x > 1 && x >= 0 && y >= 0, {x, y}]; *)

	(* 12, 3, 4. Related to 3,4,6,4 but not sure what exactly this is *)
	sol = Solve[Norm[image1-vertex1]==Norm[image3 - vertex1] && Norm[image1-vertex1]==Norm[image2 - vertex1] && x > 1 && x >= 0 && y >= 0, {x, y}];

	Print[sol];
	vertex1 = vertex1/.sol[[1]];
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
		Print["Each face contains " <> IntegerString[Length[face1]] <> " vertices"];
		AppendTo[faces, face1];
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

	scaledFaces = Map[(# / vertex1[[2]] + {1- vertex1[[1]]/vertex1[[2]], 0, 0} )&, faces, {2}];
	scaledFaces
]

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
				Opacity[1], FaceForm[White], Map[Sphere[#[[1]], thickness]&, verticesTally]
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
	range = radius * 0.8;
	thickness = range * 0.008;

	graph = Graphics3D[
		elements,
		ViewPoint -> viewPointBeforeRotation.Transpose[rotationLighting] * 100,
		ViewVertical -> If[topView, {0, 1, 0}, {1, 0, 0}].Transpose[rotationLighting],
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

	faces = If[{p, q}=={3,6}, generateFaces[6, 3], generateFaces[p, q]];

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
		If[{p, q} == {3,6},
			mod = Mod[faces[[k]][[1]][[2]] + 1, 3/2];
			If[mod > epsilon && mod < 3/2 - epsilon, faces[[k]] = RotateLeft[faces[[k]]]];
		];
	];

	(* Map[Print[Round[#, 0.001]]&, faces]; *)



	stageCount = If[{p,q}=={6,3}, 3, 2];

	For[stage = 1, stage <= stageCount * (1 + stillStages), stage++,
		For[srId = 1, srId <= framesPerStage, srId ++,
			sr = srId/framesPerStage;
			Print["sr = " <> ToString[N[sr]]];
			zOffset = If[{p, q}=={4,4}, 2, Sqrt[3]];
			movedFaces = Map[(# + {0, 0, -zOffset Mod[sr, 1] })&, faces, {2}];

			truncateStage = Ceiling[stage/(1 + stillStages)];
			truncateSr = If[truncateStage * (1 + stillStages) == stage, sr, 0];
			component1 = movedFaces;

			If[{p, q}=={6,3}, 
				component1 = Flatten[ Map[truncateFace63[#, truncateStage, truncateSr]&, movedFaces], 1];
			];

			If[{p,q}=={4,4},
				component1 = Flatten[ Map[truncateFace44[#, truncateStage, truncateSr]&, movedFaces], 1]; 
			];

			If[{p,q}=={3,6},
				component1 = Flatten[ Map[truncateFace36[#, truncateStage, truncateSr]&, movedFaces], 1];
			];

			components = {component1};

			If[{p, q}=={3,6} || {p, q}=={6,3},
				component2 = Map[(# + {0, 1/2, Sqrt[3]/2})&, component1, {2}];
				component3 = Map[(# + {0, -1/2, Sqrt[3]/2})&, component1, {2}];

				If[{p, q}=={3,6}, component3 = Map[rotateFace, component3]];
				components = If[{p, q}=={6,3}, {component1, component2, component3}, {component1, component3}];

				If[{p,q}=={4,4},
					component2 = Map[(# + {0, 1, 0})&, component1, {2}];
					components = {component1, component2}
				];
			];

			Print[components//Length];
			Print[components[[1]]//Length];

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

truncateFace63[f_, stage_, sr_] := 
	If[stage == 1,
		{
			{
				f[[1]] (sr/3 + 2/3) + (1/3 - sr/3) f[[3]],
				f[[3]] (sr/3 + 2/3) + (1/3 - sr/3) f[[1]],
				f[[3]] (sr/3 + 2/3) + (1/3 - sr/3) f[[5]],
				f[[5]] (sr/3 + 2/3) + (1/3 - sr/3) f[[3]],
				f[[5]] (sr/3 + 2/3) + (1/3 - sr/3) f[[1]],
				f[[1]] (sr/3 + 2/3) + (1/3 - sr/3) f[[5]]
			},
			{
				f[[2]] (sr/3 + 2/3) + (1/3 - sr/3) f[[4]],
				f[[4]] (sr/3 + 2/3) + (1/3 - sr/3) f[[2]],
				f[[4]] (sr/3 + 2/3) + (1/3 - sr/3) f[[6]],
				f[[6]] (sr/3 + 2/3) + (1/3 - sr/3) f[[4]],
				f[[6]] (sr/3 + 2/3) + (1/3 - sr/3) f[[2]],
				f[[2]] (sr/3 + 2/3) + (1/3 - sr/3) f[[6]]
			}	
		},
		If[stage == 2, {
			{
				Mean[f] + (f[[1]] - Mean[f]).Rotation[sr Pi/3],
				f[[1]],
				Mean[f] + (f[[3]] - Mean[f]).Rotation[sr Pi/3],
				f[[3]],
				Mean[f] + (f[[5]] - Mean[f]).Rotation[sr Pi/3],
				f[[5]]
			},
			{
				Mean[f] + (f[[2]] - Mean[f]).Rotation[sr Pi/3],
				f[[2]],
				Mean[f] + (f[[4]] - Mean[f]).Rotation[sr Pi/3],
				f[[4]],
				Mean[f] + (f[[6]] - Mean[f]).Rotation[sr Pi/3],
				f[[6]]
			}
		},
		{
			{
				f[[1]] (1-sr/3) + (sr/3) f[[3]],
				f[[2]] (1-sr/3) + (sr/3) f[[4]],
				f[[3]] (1-sr/3) + (sr/3) f[[5]],
				f[[4]] (1-sr/3) + (sr/3) f[[6]],
				f[[5]] (1-sr/3) + (sr/3) f[[1]],
				f[[6]] (1-sr/3) + (sr/3) f[[2]]
			}
		}]
	];

truncateFace44[f_, stage_, sr_] := 
	If[stage == 1,
		{
			{
				f[[1]] (1+sr)/2 + (1-sr)/2 f[[2]],
				f[[2]] (1+sr)/2 + (1-sr)/2 f[[3]],
				f[[3]] (1+sr)/2 + (1-sr)/2 f[[4]],
				f[[4]] (1+sr)/2 + (1-sr)/2 f[[1]]
			}
		},
		{
			{
				f[[2]] (2-sr)/2 + (sr)/2 f[[1]],
				f[[3]] (2-sr)/2 + (sr)/2 f[[2]],
				f[[4]] (2-sr)/2 + (sr)/2 f[[3]],
				f[[1]] (2-sr)/2 + (sr)/2 f[[4]]
			}
		}
	];

rotateFace[f_] := Map[Mean[f] + (# - Mean[f]).Rotation[Pi/3]&, f];

truncateFace36[f_, stage_, sr_] := 
	If[stage == 1,
		{
			{
				Mean[f] + (f[[1]] - Mean[f]).Rotation[sr Pi/3],
				f[[1]],
				Mean[f] + (f[[3]] - Mean[f]).Rotation[sr Pi/3],
				f[[3]],
				Mean[f] + (f[[5]] - Mean[f]).Rotation[sr Pi/3],
				f[[5]]
			}
		},
		{
			{
				Mean[f] + (f[[1]] - Mean[f]).Rotation[(sr + 1)Pi/3],
				f[[1]],
				Mean[f] + (f[[3]] - Mean[f]).Rotation[(sr + 1) Pi/3],
				f[[3]],
				Mean[f] + (f[[5]] - Mean[f]).Rotation[(sr + 1) Pi/3],
				f[[5]]
			}
		}
	];


(* exportAnimation[]; *)

faces = generateFaces[6, 3];
(* components = {faces}; *)
groupedFaces = GroupBy[faces, Length];
sortedKey = Sort[Keys[groupedFaces]];
components = Map[groupedFaces, sortedKey];

componentsStereo = If[projectToBall, Map[radius stereo[# / radius]&, components, {3}], components];
exportFrame[componentsStereo, 0, 1, "test"];



(* 

for ($num = 69 ; $num -le 99 ; $num++) { cp spinning_frame_0068.png "spinning_frame_00$num.png" }

 *)