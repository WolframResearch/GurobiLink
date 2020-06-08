(* Mathematica Package *)

BeginPackage["GUROBILink`"]

$GUROBILinkLibrary::usage  = "$GUROBILinkLibrary is the full path to the GUROBILink library loaded."
$GUROBILinkDirectory::usage = "$GUROBILinkDirectory gives the location of the GUROBILink library."
LoadGUROBILink::usage  = "LoadGUROBILink[] loads the GUROBILink library."
GUROBISetVariableTypesAndObjectiveVector::usage = "GUROBISetVariableTypesAndObjectiveVector[data, vartypes, objvector], vartypes is a string with \"C\" at positions of continuous variables and \"I\" for integer variables."
GUROBISetVariableTypesAndBoundsAndObjectiveVector::usage = "GUROBISetVariableTypesAndBoundsAndObjectiveVector[data, vartypes, lowerbounds, upperbounds, objvector], vartypes is a string with \"C\" at positions of continuous variables and \"I\" for integer variables."
GUROBIAddQuadraticObjectiveMatrix::usage = "GUROBIAddQuadraticObjectiveMatrix[data, saQmat]"
GUROBIAddLinearConstraintIndices::usage = "GUROBIAddLinearConstraintIndices[data, indices, values, sense, rhs], where sense is \">\", \"<\" or \"=\"."
GUROBIAddLinearConstraint::usage = "GUROBIAddLinearConstraint[data, a, sense, b] -> a.x sense b, where sense is \">\", \"<\" or \"=\"."
GUROBIAddLinearConstraints::usage = "GUROBIAddLinearConstraints[data, mat, sense, rhs], where sense is \">\", \"<\" or \"=\"."
GUROBIAddQuadraticConstraintIndices::usage = "GUROBIAddQuadraticConstraintIndices[data, linind, linvals, quadrow, quadcol, quadvals, sense, rhs], where sense is \">\", \"<\" or \"=\"."
GUROBIAddQuadraticConstraint::usage = "GUROBIAddQuadraticConstraint[data, Qmat, qvec, sense, rhs] --> 1/2 x^T.Q.x + q.x sense b, where sense is \">\", \"<\" or \"=\"."
GUROBIAddSOCMembershipConstraint::usage = "GUROBIAddSOCMembershipConstraint[data, ind]"
GUROBIAddSOCAffineConstraint::usage = "GUROBIAddSOCAffineConstraint[data, A, b] adds the constaraint VectorGreaterEqual[{Ax + b, 0}, {\"NormCone\", n}]"

GUROBIOptimize::usage = "GUROBIOptimize[data]"

GUROBISolutionStatus::usage ="GUROBISolutionStatus[data]"
GUROBIx::usage ="GUROBIx[data]"
GUROBIObjectiveValue::usage = "GUROBIObjectiveValue[data]"
GUROBISlack::usage = "GUROBISlack[data]"	

GUROBILink::usage = "GUROBILink is used as a symbol for message heads from GUROBI error codes."




Begin["`Private`"]
(* Implementation of the package *)

$GUROBILinkDirectory = DirectoryName[$InputFileName];
$targetDir = FileNameJoin[{$GUROBILinkDirectory, "LibraryResources", $SystemID}];
(* $targetDir ="C:\\Users\\ninad\\Desktop\\dt\\GUROBI";*)

$GUROBILinkLibrary = Block[{$LibraryPath = $targetDir}, FindLibrary["GUROBILink"]];
(* $GUROBILinkLibrary = "C:\\Users\\ninad\\Desktop\\dt\\GUROBI\\GUROBILink.dll"; *)

$GUROBILibrariesToPreload = Switch[$SystemID,
	"Windows-x86-64",
		Map[Function[FileNameJoin[{$targetDir, #}]], {"GUROBI90.dll"}],
	_,
		{}
]

(*
 Load all the functions from the GUROBILink library
*)

$GUROBIPrintLevel = 0;
dPrint = Optimization`Debug`OptimizationDebugPrint;
(* Set print level 5 with Optimization`Debug`SetPrintLevel[5] *)
GUROBILink`Private`dPrint = Print;
(*GUROBILink`Private`pPrint = Print;*)

needInitialization = True;

$GUROBIInfinity = 1.*10^30;

(*
GUROBISetVariableTypesAndObjectiveVector[data, vartypes, objvector]
GUROBISetVariableTypesAndBoundsAndObjectiveVector[data, vartypes, lowerbounds, upperbounds, objvector]
GUROBIAddQuadraticObjectiveMatrix[data, Qmat]
GUROBIAddLinearConstraint[data, vec, sense, rhs]
GUROBIAddLinearConstraintIndices[data, indices, values, sense, rhs]
GUROBIAddLinearConstraints[data, mat, sense, rhs]
GUROBIAddQuadraticConstraint[data, Qmat, qvec, sense, rhs]
GUROBIAddQuadraticConstraintIndices[data, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]
GUROBIAddSOCMembershipConstraint[data, ind]
GUROBIAddSOCAffineConstraint[data, mat, vec]
GUROBIOptimize[data]
GUROBISolutionStatus[data]
GUROBIObjectiveValue[data]
GUROBIx[data]
GUROBISlack[data]
*)


LoadGUROBILink[] :=
Block[{$LibraryPath = $targetDir}, 
	Map[LibraryLoad, $GUROBILibrariesToPreload];
	GUROBISetVariableTypesAndObjectiveVector0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetVariableTypesAndObjectiveVector", {Integer, UTF8String, {Real, 1}}, Integer];
	GUROBISetVariableTypesAndBoundsAndObjectiveVector0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetVariableTypesAndBoundsAndObjectiveVector", {Integer, UTF8String, {Real, 1}, {Real, 1}, {Real, 1}}, Integer];
	GUROBIAddQuadraticObjectiveMatrix0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddQuadraticObjectiveMatrix", {Integer, LibraryDataType[SparseArray, Real, 2]}, Integer];
	GUROBIAddLinearConstraint0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddLinearConstraint", {Integer, {Integer, 1}, {Real, 1}, UTF8String, Real}, Integer];
	GUROBIAddLinearConstraint1 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddLinearConstraint1", {Integer, LibraryDataType[SparseArray, Real, 1], UTF8String, Real}, Integer];
	GUROBIAddLinearConstraints0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddLinearConstraints", {Integer, LibraryDataType[SparseArray, Real, 2], UTF8String, {Real, 1}}, Integer];
	GUROBIAddQuadraticConstraint0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddQuadraticConstraint", {Integer, {Integer, 1}, {Real, 1}, {Integer, 1}, {Integer, 1}, {Real, 1}, UTF8String, Real}, Integer];
	GUROBIAddQuadraticConstraint1 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddQuadraticConstraint1", {Integer, LibraryDataType[SparseArray, Real, 2], LibraryDataType[SparseArray, Real, 1], UTF8String, Real}, Integer];

	GUROBIOptimize0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_OptimizeModel", {Integer}, Integer];

	(*GUROBIGetSolution0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_GetSolution", {Integer, Integer}, Integer];*)
	GUROBISolutionStatus0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_GetSolutionStatus", {Integer}, Integer];
	GUROBIObjectiveValue0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_GetObjectiveValue", {Integer}, Real];
	GUROBIx0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_Getx", {Integer}, {Real, 1}];
	GUROBISlack0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_GetSlack", {Integer}, {Real, 1}];

	GUROBIDataDelete0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIDataMap_delete", {Integer}, Integer];
	GUROBIDataIDList = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIDataMap_retIDList", {}, {Integer, 1}];
	needInitialization = False;
]

LoadGUROBILink[]


(* GUROBIData expression (GUROBISolMap) related: *)

GUROBIDataID[e_GUROBIData] := ManagedLibraryExpressionID[e, "GUROBI_data_instance_manager"];

GUROBIDataQ[e_GUROBIData] := ManagedLibraryExpressionQ[e, "GUROBI_data_instance_manager"];
GUROBIDataQ[_] := False;

testGUROBIData[][e_] := testGUROBIData[GUROBIData][e];
testGUROBIData[mhead_Symbol][e_] :=
If[TrueQ[GUROBIDataQ[e]],
	True,
	Message[MessageName[mhead, "GUROBIinst"], e]; False
];
testGUROBIData[_][e_] := TrueQ[GUROBIDataQ[e]];

General::GUROBIinst = "`1` does not represent an active GUROBIData object.";

GUROBIDataCreate[] :=
Module[{},
	If[needInitialization, LoadGUROBILink[]];
	CreateManagedLibraryExpression["GUROBI_data_instance_manager", GUROBIData]
];

GUROBIDataDelete[GUROBIData[id_]?(testGUROBIData[GUROBIDataDelete])] := GUROBIDataDelete0[id];

GUROBIDataDelete[l:{_GUROBIData..}] := GUROBIDataDelete /@ l;

GUROBIDataExpressions[] :=
Module[{list},
	If[needInitialization, LoadGUROBI[]];
	list = GUROBIDataIDList[];
	If[!ListQ[list],
	   $Failed,
	   Map[GUROBIData, list]]
]


(* Functions to set up the problem *)


GUROBISetVariableTypesAndObjectiveVector[GUROBIData[id_]?(testGUROBIData[GUROBISetVariableTypesAndObjectiveVector]), vartypes_, objvector_]:=
Module[{},
	GUROBISetVariableTypesAndObjectiveVector0[id, vartypes, objvector]
];

GUROBISetVariableTypesAndBoundsAndObjectiveVector[GUROBIData[id_]?(testGUROBIData[GUROBISetVariableTypesAndBoundsAndObjectiveVector]), vartypes_, lb_, ub_, objvector_]:=
Module[{},
	GUROBISetVariableTypesAndBoundsAndObjectiveVectort0[id, vartypes, lb, ub, objvector]
];

GUROBIAddQuadraticObjectiveMatrix[GUROBIData[id_]?(testGUROBIData[GUROBIAddQuadraticObjectiveMatrix]), Qmat_SparseArray]:=
Module[{},
	GUROBIAddQuadraticObjectiveMatrix0[id, Qmat]
];

GUROBIAddLinearConstraintIndices[GUROBIData[id_]?(testGUROBIData[GUROBIDataDelete]), indices_, values_, sense_, rhs_]:=
Module[{},
	GUROBIAddLinearConstraint0[id, indices, values, sense, rhs]
];

GUROBIAddLinearConstraint[GUROBIData[id_]?(testGUROBIData[GUROBIDataDelete]), vector_SparseArray, sense_, rhs_]:=
Module[{},
	GUROBIAddLinearConstraint1[id, vector, sense, rhs]
];

GUROBIAddLinearConstraints[GUROBIData[id_]?(testGUROBIData[GUROBIDataDelete]), mat_SparseArray, sense_, rhs_]:=
Module[{},
Print["in GUROBIAddLinearConstraints"];
Print[xGUROBIAddLinearConstraints0[id, mat, sense, rhs]];
	GUROBIAddLinearConstraints0[id, mat, sense, rhs]
];


GUROBIAddQuadraticConstraintIndices[GUROBIData[id_]?(testGUROBIData[GUROBIAddQuadraticConstraint]), linind_, linvals_, quadrow_, quadcol_, quadvals_, sense_, rhs_] :=
Module[{},
	GUROBIAddQuadraticConstraint0[data, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]
]

GUROBIAddQuadraticConstraint[GUROBIData[id_]?(testGUROBIData[GUROBIAddQuadraticOptimizationConstraint]), Q_SparseArray, q_SparseArray, sense_, b_] :=
Module[{},
	(* 1/2 x^T.Q.x + q.x <= b *)
	(* accept any matrix Q and vector q inputs and make them sparse if they are not *)
	GUROBIAddQuadraticConstraint1[id, Q, q, sense, b]
]

GUROBIAddSOCMembershipConstraint[GUROBIData[id_]?(testGUROBIData[GUROBIAddSOCMembershipConstraint]), ind_]:=
Module[{n, linind, linvals, quadrow, quadcol, quadvals, sense, rhs},
	(* x1^2+...+x(n-1)^2-xn^2 <= 0, xn>=0 *)
Print["In GUROBIAddSOCMembershipConstraint"];
	n = Length[ind];
Print[xGUROBIAddLinearConstraint0[id, {ind[[-1]]}, {1}, ">", 0]];
	GUROBIAddLinearConstraint0[id, {ind[[-1]]}, {1}, ">", 0];
	linind = {};
	linvals = {};
	quadrow = ind;
	quadcol = ind;
	quadvals = Append[ConstantArray[1, n-1], -1];
	sense = "<";
	rhs = 0;
	(*{Integer, {Integer, 1}, {Real, 1}, {Integer, 1}, {Integer, 1}, {Real, 1}, UTF8String, Real}, Integer];*)
Print[xGUROBIAddQuadraticConstraint0[id, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]];
  	GUROBIAddQuadraticConstraint0[id, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]
]

GUROBIAddSOCAffineConstraint[GUROBIData[id_]?(testGUROBIData[GUROBIAddSOCConstraintAsQuadratic]), A_, b_]:=
Module[{n, Q, q, b1},
	(* it is not clear if this will be accepted by the GUROBI solver for convex quadratic constraint,
	if not one may try to set the "NonConvex" parameter to 2 but this should not be very efficient *)
	Print["In GUROBIAddSOCAffineConstraint"];
	
	(*VectorGreaterEqual[{Ax + b, 0}, {"NormCone", n}] --> 1/2 x^T.Q.x + q.x <= b *)
	(*
	an.x + bn >= 0,
	||{a1.x + b1, ..., a (n - 1).x + b (n - 1)}|| <= an.x + bn
	Q = a1 a1^T + ... + a (n - 1).a (n - 1)^T - an^2
	b = -(b1^2 + ... + b (n - 1)^2 - bn^2)
	q = 2 (b1*a1 + ... b (n - 1)*a (n - 1) - bn*an) 
	*)
	n = Length[b];
Print["b: ", b];
Print["n: ", n];
Print[ A[[n]] ];
Print[ b[[n]] ];
Print[xGUROBIAddLinearConstraint1[id, SparseArray[A[[n]]], ">", -b[[n]]]];
	GUROBIAddLinearConstraint1[id, SparseArray[A[[n]]], ">", -b[[n]]];
	Q = Sum[Transpose[{A[[i]]}].{A[[i]]}, {i, 1, n-1}] - Transpose[{A[[n]]}].{A[[n]]};
	b1 = -Sum[b[[i]]^2, {i, 1, n-1}] + b[[n]]^2;
Print["b1: ", b1];
	q = 2*(Sum[b[[i]]*A[[i]], {i, 1, n-1}] - b[[n]]*A[[n]]);
Print[xGUROBIAddQuadraticConstraint1[id, SparseArray[Q], SparseArray[q], "<", b1]];
	GUROBIAddQuadraticConstraint1[id, SparseArray[Q], SparseArray[q], "<", b1]
]


(* Function to solve the problem *)

Options[GUROBIOptimize] = {Method->Automatic, MaxIterations->Automatic, Tolerance->Automatic, PerformanceGoal:>$PerformanceGoal(*, "Verbose"->False*)}

GUROBIOptimize[GUROBIData[id_]?(testGUROBIData[GUROBIOptimize]), OptionsPattern[GUROBIOptimize]] :=
Module[{tol = OptionValue[Tolerance], maxiter = OptionValue[MaxIterations], status, error},
	If[SameQ[tol, Automatic], tol = 10.^-8, tol = N[tol]];
	If[SameQ[maxiter, Automatic], maxiter = 400];
	dPrint[1, "In GUROBISolve"];
	dPrint[3, "tol: ", tol];
	dPrint[3, "maxiter: ", maxiter];

	dPrint[1, "Solving with GUROBI..."];
	(*error = GUROBIOptimize0[id, tol, maxiter];*)
	error = GUROBIOptimize0[id];
Print["error: ", error];
	status = GUROBISolutionStatus0[id];
Print["status: ", status]; status
	(* TODO: make GUROBIStringStatus string status 
	and fill them in in Optimization`SolutionData`*)
];

(* Selectors *)

GUROBISolutionStatus[GUROBIData[id_]?(testGUROBIData[GUROBISolutionStatus])] := GUROBISolutionStatus0[id];

GUROBIx[GUROBIData[id_]?(testGUROBIData[GUROBIx])] := GUROBIx0[id];

GUROBIObjectiveValue[GUROBIData[id_]?(testGUROBIData[GUROBIObjectiveValue])] := GUROBIObjectiveValue0[id];

GUROBISlack[GUROBIData[id_]?(testGUROBIData[GUROBISlack])] := GUROBISlack0[id];


(* Register convex method *)

Optimization`ConvexSolvers`RegisterConvexMethod["GUROBI1", 
	Association[
		"SolveFunction" -> GUROBISolve1,
		"ObjectiveSupport" -> "Quadratic",
  		"ConicConstraintSupport" -> Thread[{"EqualityConstraint", "NonNegativeCone", "NormCone"}->"Affine"],
		"MixedIntegerSupport"->True
  	]
]

Optimization`ConvexSolvers`RegisterConvexMethod["GUROBI2", 
	Association[
		"SolveFunction" -> GUROBISolve2,
		"ObjectiveSupport" -> "Quadratic",
  		"ConicConstraintSupport" -> {"EqualityConstraint"->"Affine", "NonNegativeCone"->"Affine", "NormCone"->"Membership"},
 		"MixedIntegerSupport"->True
 	]
]

Optimization`ConvexSolvers`RegisterConvexMethod["GUROBI", 
	Association[
		"Submethods"->{"GUROBI1", "GUROBI2"}
  	]
]

GUROBISolve1[problemData_, pmopts___] :=
Module[{t0, data, objvec, objmat, nvars, affine, coneSpecifications, lpos, intvars, vartypes},
	(* only for fully suppoted by GUROBI condstraints *)
	t0 = AbsoluteTime[];

	dPrint[1, "In GUROBISolve"];
	dPrint[3, "pmopts: ", pmopts];
Print[InputForm[problemData]];

	data = GUROBIDataCreate[];
	
	objvec = problemData["ObjectiveVector"];
	nvars = Length[objvec];
	objmat = problemData["ObjectiveMatrix"];
	intvars = problemData["IntegerVariableColumns"];
	vartypes = StringJoin[Table[If[MemberQ[intvars, i], "I", "C"], {i, nvars}]];
Print["nvars: ", nvars];
Print["intvars: ", intvars];
Print["vartypes: ", vartypes];
	GUROBISetVariableTypesAndObjectiveVector[data, vartypes, objvec];
	GUROBIAddQuadraticObjectiveMatrix[data, SparseArray[objmat]];
	
	affine = problemData["ConicConstraintAffineLists"];
	coneSpecifications = problemData["ConicConstraintConeSpecifications"];
Print[coneSpecifications];

	(* "EqualityConstraint" will always come first and then "NonNegativeCone" *)
	lpos = 1;
	If[Length[coneSpecifications]>=1 && MatchQ[coneSpecifications[[1]],{"EqualityConstraint", _}],
		{a, b} = affine[[1]];
		GUROBIAddLinearConstraints[data, SparseArray[a], "=", -b];
		lpos = 2;
		Print["eq {a,b} -> ", {a,b}];
	];

	If[Length[coneSpecifications]>=lpos && MatchQ[coneSpecifications[[lpos]],{"NonNegativeCone", _}],
		{a, b} = affine[[lpos]];
		GUROBIAddLinearConstraints[data, SparseArray[a], ">", -b];
		lpos += 1;
		Print["ineq {a,b}-> ", {a,b}];
	];
	If[Length[coneSpecifications]>=lpos && MatchQ[coneSpecifications[[lpos]], {"NormCone", _}],
		{a, b} = affine[[lpos]];
		GUROBIAddSOCAffineConstraint[data, a, b]
		Print["norm {a,b}-> ", {a,b}];
	];


	(*
	If GUROBIAddSOCAffineConstraint won't work, use (method GURPBI2 with)
		GUROBIAddSOCMembershipConstraint[data, ind];
	*)

	status = GUROBIOptimize[data, pmopts];
	(*If[!StringQ[status], Return[$Failed]];*)
	status = "Solved";
	status = Optimization`SolutionData`WrapStatus[status];
Print["obj val: ", GUROBIObjectiveValue[data]];
Print["x: ",GUROBIx[data]];

	dPrint[1, "status: ", status];
	{status, GUROBI1Data[data, {}]}
]

GUROBI1Data[data_, _]["PrimalMinimumValue"] := GUROBIObjectiveValue[data];
GUROBI1Data[data_, _]["PrimalMinimizerVector"] := GUROBIx[data];
GUROBI1Data[data_, _]["Slack"] := GUROBISlack[data];
GUROBI1Data[data_, _]["DualMaximumValue"] := Missing["NotAvailable"];
GUROBI1Data[data_, _]["DualityGap"] := Missing["NotAvailable"];
GUROBI1Data[data_, _]["DualMaximizer"] := Missing["NotAvailable"];


GUROBISolve2[problemData_, pmopts___] :=
Module[{a, b, objvec, data, nvars, status, lpos, nextra, norig, integerColumns, t0, 
	affine, coneSpecifications, coneVariableIndexes, vartypes},
	(* corresponds to MOSEKPrimal *)
	t0 = AbsoluteTime[];

	dPrint[1, "In GUROBISolve2"];
	dPrint[3, "pmopts: ", pmopts];
	data = GUROBIDataCreate[];
	
	objvec = problemData["ObjectiveVector"];
	nvars = Length[objvec];
	objmat = problemData["ObjectiveMatrix"];
	nextra = Lookup[problemData, "ExtraColumns", 0];
	norig = nvars - nextra;

	affine = problemData["ConicConstraintAffineLists"];
	coneSpecifications = problemData["ConicConstraintConeSpecifications"];
 	coneVariableIndexes = problemData["ConeVariableColumns"];
	coneLengths = Map[Length, affine];
Print["objvec: ", objvec];
Print["affine: ", affine];
Print["coneSpecifications: ", coneSpecifications];
	If[MatchQ[coneVariableIndexes, _Missing], coneVariableIndexes = ConstantArray[None, Length[coneSpecifications]]];

Print["coneVariableIndexes: ", coneVariableIndexes];

	integerColumns = problemData["IntegerVariableColumns"];

	vartypes = StringJoin[Table[If[MemberQ[integerColumns, i], "I", "C"], {i, nvars}]];
Print["vartypes: ", vartypes];

	
	GUROBISetVariableTypesAndObjectiveVector[data, vartypes, objvec];
	GUROBIAddQuadraticObjectiveMatrix[data, SparseArray[objmat]];

	(* add variable bounds *)
	(*pleq = Flatten[Position[coneSpecifications, {"EqualityConstraint", _}, {1}]];
	plin = Flatten[Position[coneSpecifications, {"NonNegativeCone", _}, {1}]];
	nleq = Total[coneLengths[[pleq]]];
	nlin = Total[coneLengths[[plin]]];
	Print["pleq: ", pleq];
	Print["plin: ", plin];
	Print["nleq: ", nleq];
	Print["nlin: ", nlin];
	GUROBISetVariableTypesAndZeroLowerBoundsForLinIneqDualVarsAndObjectiveVector[data, vartypes, nleq, nlin, objvector];*)
	
	(* "EqualityConstraint" will always come first and then "NonNegativeCone" *)
	lpos = 1;
	If[Length[coneSpecifications]>=1 && MatchQ[coneSpecifications[[1]],{"EqualityConstraint", _}],
		{a, b} = affine[[1]];
		GUROBIAddLinearConstraints[data, SparseArray[a], "=", -b];
		lpos = 2;
	];
	If[Length[coneSpecifications]>=lpos && MatchQ[coneSpecifications[[lpos]],{"NonNegativeCone", _}],
		{a, b} = affine[[lpos]];
		GUROBIAddLinearConstraints[data, SparseArray[a], ">", -b];
		lpos += 1;
	];
	If[Length[coneSpecifications]>=lpos && MatchQ[coneSpecifications[[lpos]], {"NormCone", _}],
		GUROBIAddSOCMembershipConstraint[data, coneVariableIndexes[[lpos]]]];

	(* GUROBISolve: *)
	status = GUROBIOptimize[data, pmopts];
	(*If[!StringQ[status], Return[$Failed]];*)
	status = "Solved";
	status = Optimization`SolutionData`WrapStatus[status];
	dPrint[1, "status: ", status];
Print["obj val: ", GUROBIObjectiveValue[data]];
Print["x: ",GUROBIx[data]];
	{status, GUROBI2Data[data, {integerColumns}]}
];

GUROBI2Data[data_, _]["PrimalMinimumValue"] := GUROBIObjectiveValue[data];
(*GUROBI2Data[data_, {integerColumns_}]["PrimalMinimizerVector"] :=  GUROBIx[data];*)
GUROBI2Data[data_, {integerColumns_}]["PrimalMinimizerVector"] := 
Block[{x = GUROBIx[data], res},
Print["In PrimalMInimizerVector"];
	If[Length[integerColumns] > 0,
		x[[integerColumns]] = Round[x[[integerColumns]]]
	];
	x
];
GUROBI2Data[data_, _]["DualMaximumValue"] := Missing["NotAvailable"];
GUROBI2Data[data_, _]["DualityGap"] := Missing["NotAvailable"];
GUROBI2Data[data_, _]["DualMaximizer"] := Missing["NotAvailable"];
GUROBI2Data[data_, _]["Slack"] := Missing["NotAvailable"];


End[]
EndPackage[]