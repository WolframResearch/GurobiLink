(* Mathematica Package *)

BeginPackage["GUROBILink`"]

$GUROBILinkLibrary::usage  = "$GUROBILinkLibrary is the full path to the GUROBILink library loaded."
$GUROBILinkDirectory::usage = "$GUROBILinkDirectory gives the location of the GUROBILink library."
LoadGUROBILink::usage  = "LoadGUROBILink[] loads the GUROBILink library."
GUROBILink::usage = "GUROBILink is used as a symbol for message heads from GUROBI error codes."
GUROBICheckLicense::usage = "GUROBICheckLicense[data]"
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

GUROBIStatusValues::usage ="GUROBIStatusValues[data]"
GUROBIx::usage ="GUROBIx[data]"
GUROBIObjectiveValue::usage = "GUROBIObjectiveValue[data]"
GUROBISlack::usage = "GUROBISlack[data]"	

GUROBIData::usage = "GUROBIData[id] represents an instance of an GUROBIData expression created by GUROBIDataCreate."
GUROBIDataID::usage = "GUROBIDataID[data] gives the instance id of an GUROBIData expression data."
GUROBIDataQ::usage = "GUROBIDataQ[expr] gives True if expr represents an active instance of an GUROBIData object."
GUROBIDataCreate::usage = "GUROBIDataCreate[] creates an instance of an GUROBIData expression."
GUROBIDataExpressions::ussage = "GUROBIDataExpressions[] shows all active GUROBIData expression instances."
GUROBIDataDelete::usage = "GUROBIDataDelete[expr] removes an instance of an GUROBIData expression, freeing up memory."


(*TODO:
add initial point
add more messages
add statuses to kernel
*) 

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
GUROBISetPArameters[data, maxit, tol, nonconvex]
GUROBIOptimize[data]
GUROBIStatusValue[data]
GUROBIObjectiveValue[data]
GUROBIx[data]
GUROBISlack[data]
*)


LoadGUROBILink[] :=
Block[{$LibraryPath = $targetDir}, 
	Map[LibraryLoad, $GUROBILibrariesToPreload];
	GUROBICheckLicense0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_CheckLicense", {Integer}, Integer];
	GUROBISetVariableTypesAndObjectiveVector0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetVariableTypesAndObjectiveVector", {Integer, UTF8String, {Real, 1}}, Integer];
	GUROBISetVariableTypesAndBoundsAndObjectiveVector0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetVariableTypesAndBoundsAndObjectiveVector", {Integer, UTF8String, {Real, 1}, {Real, 1}, {Real, 1}}, Integer];
	GUROBIAddQuadraticObjectiveMatrix0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddQuadraticObjectiveMatrix", {Integer, LibraryDataType[SparseArray, Real, 2]}, Integer];
	GUROBIAddLinearConstraint0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddLinearConstraint", {Integer, {Integer, 1}, {Real, 1}, UTF8String, Real}, Integer];
	GUROBIAddLinearConstraint1 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddLinearConstraint1", {Integer, LibraryDataType[SparseArray, Real, 1], UTF8String, Real}, Integer];
	GUROBIAddLinearConstraints0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddLinearConstraints", {Integer, LibraryDataType[SparseArray, Real, 2], UTF8String, {Real, 1}}, Integer];
	GUROBIAddQuadraticConstraint0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddQuadraticConstraint", {Integer, {Integer, 1}, {Real, 1}, {Integer, 1}, {Integer, 1}, {Real, 1}, UTF8String, Real}, Integer];
	GUROBIAddQuadraticConstraint1 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddQuadraticConstraint1", {Integer, LibraryDataType[SparseArray, Real, 2], LibraryDataType[SparseArray, Real, 1], UTF8String, Real}, Integer];
	(*GUROBISetNumberOfConstraints0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetNumberOfConstraints", {Integer, Integer}, Integer];*)
	GUROBISetParameters0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetParameters", {Integer, Integer, Real, Integer}, Integer];
	GUROBISetStartingPoint0 =  LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetStartingPoint", {Integer, {Real, 1}}, Integer];

	GUROBIOptimize0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_OptimizeModel", {Integer}, Integer];

	(*GUROBIGetSolution0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_GetSolution", {Integer, Integer}, Integer];*)
	GUROBIStatusValue0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_GetStatusValue", {Integer}, Integer];
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

General::gurobiinst = "`1` does not represent an active GUROBIData object.";

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
(* Check License *)

GUROBICheckLicense[GUROBIData[id_]?(testGUROBIData[GUROBICheckLicense])]:=
Module[{},
	dPrint[1, "Checking for GUROBI license..."];
	res = GUROBICheckLicense0[id];
	If[res === 0,
		dPrint[1, "... found license."];
		True
		, (* else *)
		Message[GUROBILink::license];
		False
	]
]

GUROBILink::license = "Cannot find a valid GUROBI license for version 9.0 or greater." 

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
	dPrint[5, Normal[xGUROBIAddLinearConstraints0[id, mat, sense, rhs]]];
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
	dPrint[5, Normal[xGUROBIAddLinearConstraint0[id, {ind[[-1]]}, {1}, ">", 0]]];
	GUROBIAddLinearConstraint0[id, {ind[[-1]]}, {1}, ">", 0];
	linind = {};
	linvals = {};
	quadrow = ind;
	quadcol = ind;
	quadvals = Append[ConstantArray[1, n-1], -1];
	sense = "<";
	rhs = 0;
	(*{Integer, {Integer, 1}, {Real, 1}, {Integer, 1}, {Integer, 1}, {Real, 1}, UTF8String, Real}, Integer];*)
	dPrint[5, Normal[xGUROBIAddQuadraticConstraint0[id, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]]];
  	GUROBIAddQuadraticConstraint0[id, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]
]

GUROBIAddSOCAffineConstraint[GUROBIData[id_]?(testGUROBIData[GUROBIAddSOCConstraintAsQuadratic]), A_, b_]:=
Module[{n, Q, q, b1},
	(* This is only accepted for SOC constraint for diagonal A, with non-diagonal A
	one may try to set the "NonConvex" parameter to 2 but this should not be very efficient *)
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
	dPrint[5, Normal[xGUROBIAddLinearConstraint1[id, SparseArray[A[[n]]], ">", -b[[n]]]]];
	GUROBIAddLinearConstraint1[id, SparseArray[A[[n]]], ">", -b[[n]]];
	Q = Sum[Transpose[{A[[i]]}].{A[[i]]}, {i, 1, n-1}] - Transpose[{A[[n]]}].{A[[n]]};
	b1 = -Sum[b[[i]]^2, {i, 1, n-1}] + b[[n]]^2;
	q = 2*(Sum[b[[i]]*A[[i]], {i, 1, n-1}] - b[[n]]*A[[n]]);
	dPrint[5, Normal[xGUROBIAddQuadraticConstraint1[id, SparseArray[Q], SparseArray[q], "<", b1]]];
	GUROBIAddQuadraticConstraint1[id, SparseArray[Q], SparseArray[q], "<", b1]
]

GUROBISetNumberOfConstraints[GUROBIData[id_]?(testGUROBIData[GUROBISetNumberOfConstraints]), ncons_]:=
Module[{},
	GUROBISetNumberOfConstraints0[id, ncons];
]

GUROBISetParameters[GUROBIData[id_]?(testGUROBIData[GUROBISetNumberOfConstraints]), maxit_, tol_, nonconvex_]:=
Module[{}
	error = GUROBISetParameters0[id, maxit, tol, nonconvex]
]

GUROBISetStartingPoint[GUROBIData[id_]?(testGUROBIData[GUROBISetNumberOfConstraints]), startpt_]:=
Module[{}
	error = GUROBISetStartingPoint0[id, N[startpt]];
]

(* Function to solve the problem *)

INTMAX = 2^31-1;

Options[GUROBIOptimize] = {Method->Automatic, MaxIterations->Automatic, Tolerance->Automatic, 
"NonConvex"->Automatic, "StartingPoint"->Automatic, "Caller"-> Automatic,
	PerformanceGoal:>$PerformanceGoal, WorkingPrecision->MachinePrecision}

GUROBIOptimize[GUROBIData[id_]?(testGUROBIData[GUROBIOptimize]), OptionsPattern[GUROBIOptimize]] :=
Module[{tol, maxiter, nonconvex, mhead, verbose, status, error, data=GUROBIData[id]},
	dPrint[1, "In GUROBIOptimize"];
	tol = OptionValue[Tolerance]; 
	maxiter = OptionValue[MaxIterations];
	nonconvex = OptionValue["NonConvex"];
	mhead = OptionValue["Caller"];
	startpt = OptionValue["StartingPoint"];
	If[SameQ[tol, Automatic], tol = 10.^-6, tol = N[tol]];
	If[SameQ[maxiter, Automatic], maxiter = 1000, If[SameQ[maxiter, Infinity], maxiter = MAXINT]];
	If[TrueQ[nonconvex] || SameQ[nonconvex, 2], nonconvex = 2, nonconvex = 1];
	If[SameQ[mhead, Automatic], mhead = ConicOptimization];
	If[!SameQ[startpt, Automatic],
		dPrint[3, "Setting starting point ", startpt];
		error = GUROBISetStartingPoint0[id, N[startpt]];
	];
	dPrint[3, "Setting parameters {maxiter, tol, nonconvex} -> ", {maxiter, tol, nonconvex}];
	error = GUROBISetParameters0[id, maxiter, tol, nonconvex];

	dPrint[1, "Solving with GUROBI..."];
	error = GUROBIOptimize0[id];
	dPrint[1, "error: ", error];
	status = GUROBIStatusValue0[id];
	dPrint["status: ", status];

	status = GUROBIStringStatus[data];
	dPrint[1, "string status: ", status];
	If[!StringQ[status], Return[$Failed]];
	ReportError[status, mhead, maxiter];
	status
];

(* https://www.gurobi.com/documentation/9.0/refman/optimization_status_codes.html *)
GUROBIStringStatus[data_] :=
	Switch[GUROBIStatusValue[data],
		1, "Loaded",
		2, "Solved",
		3, "PrimalInfeasible",
		4, "DualInfeasible",
		5, "Unbounded",
		6, "Cutoff",
		7, "MaxIterationsReached",
		8, "NodeLimitReached",
		9, "MaxCPUTimeReached",
		10, "SolutionLimit",
		11, "Interrupted",
		12, "Numeric",
		13, "Solved/Inaccurate",
		14, "InProgress",
		15, "ObjectiveLimitReached",
		_, "Unknown"
	]

ReportError[status_, mhead_, maxiter_] :=
	Switch[status,
		"Solved/Inaccurate", Message[mhead::inac, maxiter],
		"MaxIterationsReached", Message[mhead::maxit, maxiter],
		"PrimalInfeasible", Message[mhead::nsolc],
		"DualInfeasible", Message[mhead::dinfeas],
		"Interrupted", Message[mhead::userstop],
		"MaxCPUTimeReached", Message[mhead::maxcput]
	];


(* Selectors *)

GUROBIStatusValue[GUROBIData[id_]?(testGUROBIData[GUROBIStatusValues])] := GUROBIStatusValue0[id];

GUROBIx[GUROBIData[id_]?(testGUROBIData[GUROBIx])] := GUROBIx0[id];

GUROBIObjectiveValue[GUROBIData[id_]?(testGUROBIData[GUROBIObjectiveValue])] := GUROBIObjectiveValue0[id];

GUROBISlack[GUROBIData[id_]?(testGUROBIData[GUROBISlack])] := GUROBISlack0[id];


(* Register convex method *)

(* Use GUROBI1 method with membership SOC constraint 'x in K', or 'Ax+b in K' where A is diagonal *)
Optimization`ConvexSolvers`RegisterConvexMethod["GUROBI1", 
	Association[
		"SolveFunction" -> GUROBISolve1,
		"ObjectiveSupport" -> "Quadratic",
  		"ConicConstraintSupport" -> Thread[{"EqualityConstraint", "NonNegativeCone", "NormCone"}->"Affine"],
		"MixedIntegerSupport"->True
  	]
]

(* GUROBI2 is for affine SOC constraint 'Ax+b in K'. It works by adding extra variables y = Ax+b *)
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
		"Submethods"->{"GUROBI2", "GUROBI1"}
  	]
]

GUROBISolve1[problemData_, pmopts___] :=
Module[{t0, data, objvec, objmat, nvars, ncons, affine, coneSpecifications, lpos, intvars, vartypes},
	(* only for condstraints suppoted by GUROBI -- linear, quadratic and soc membership 
	for others can try the option "NonConvex"->True *)
	t0 = AbsoluteTime[];

	dPrint[1, "In GUROBISolve1"];
	dPrint[3, "pmopts: ", pmopts];

	data = GUROBIDataCreate[];
	If[!GUROBIDataQ[data] || !GUROBICheckLicense[data], Return[$Failed]];

	objvec = problemData["ObjectiveVector"];
	nvars = Length[objvec];
	objmat = problemData["ObjectiveMatrix"];
	intvars = problemData["IntegerVariableColumns"];
	vartypes = StringJoin[Table[If[MemberQ[intvars, i], "I", "C"], {i, nvars}]];
	dPrint[3, "nvars: ", nvars];
	dPrint[3, "intvars: ", intvars];
	dPrint[3, "vartypes: ", vartypes];
	GUROBISetVariableTypesAndObjectiveVector[data, vartypes, objvec];
	GUROBIAddQuadraticObjectiveMatrix[data, SparseArray[objmat]];
	
	affine = problemData["ConicConstraintAffineLists"];
	coneSpecifications = problemData["ConicConstraintConeSpecifications"];
	dPrint[3, "coneSpecifications: ", coneSpecifications];

	(* "EqualityConstraint" will always come first and then "NonNegativeCone" *)
	lpos = 1;
	If[Length[coneSpecifications]>=1 && MatchQ[coneSpecifications[[1]],{"EqualityConstraint", _}],
		{a, b} = affine[[1]];
		GUROBIAddLinearConstraints[data, SparseArray[a], "=", -b];
		lpos = 2;
		dPrint[5, "eq {a,b} -> ", {a,b}];
	];

	If[Length[coneSpecifications]>=lpos && MatchQ[coneSpecifications[[lpos]],{"NonNegativeCone", _}],
		{a, b} = affine[[lpos]];
		GUROBIAddLinearConstraints[data, SparseArray[a], ">", -b];
		lpos += 1;
		dPrint[5, "ineq {a,b}-> ", {a,b}];
	];
	If[Length[coneSpecifications]>=lpos && MatchQ[coneSpecifications[[lpos]], {"NormCone", _}],
		{a, b} = affine[[lpos]];
		GUROBIAddSOCAffineConstraint[data, a, b];
		dPrint[5, "norm {a,b}-> ", {a, b}];
		(*
		GUROBIAddSOCAffineConstraint seems to work only with diagonal matrix a, 
		otherwise use GUROBI2 which uses GUROBIAddSOCMembershipConstraint[data, ind];
		*)
		If[!DiagonalMatrixQ[a], Return[$Failed]];
	];

	status = GUROBIOptimize[data, pmopts];
	(*status = "Solved";*)
	status = GUROBIStringStatus[data];
	If[!StringQ[status], Return[$Failed]];
	status = Optimization`SolutionData`WrapStatus[status];
	(*Print["obj val: ", GUROBIObjectiveValue[data]];
	Print["x: ",GUROBIx[data]];
	Print["slack: ", GUROBISlack[data]];*)

	dPrint[1, "status: ", status];
	{status, GUROBI1Data[data, {}]}
]

GUROBI1Data[data_, _]["PrimalMinimumValue"] := GUROBIObjectiveValue[data];
GUROBI1Data[data_, _]["PrimalMinimizerVector"] := GUROBIx[data];
(*GUROBI1Data[data_, _]["Slack"] := GUROBISlack[data];*)
GUROBI1Data[data_, _]["Slack"] := Missing["NotAvailable"];
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
	If[!GUROBIDataQ[data] || !GUROBICheckLicense[data], Return[$Failed]];
	
	objvec = problemData["ObjectiveVector"];
	nvars = Length[objvec];
	objmat = problemData["ObjectiveMatrix"];
	nextra = Lookup[problemData, "ExtraColumns", 0];
	norig = nvars - nextra;

	affine = problemData["ConicConstraintAffineLists"];
	coneSpecifications = problemData["ConicConstraintConeSpecifications"];
 	coneVariableIndexes = problemData["ConeVariableColumns"];
	dPrint[5, "objvec: ", objvec];
	dPrint[5, "affine: ", affine];
	dPrint[3, "coneSpecifications: ", coneSpecifications];
	If[MatchQ[coneVariableIndexes, _Missing], coneVariableIndexes = ConstantArray[None, Length[coneSpecifications]]];

	dPrint[5, "coneVariableIndexes: ", coneVariableIndexes];

	integerColumns = problemData["IntegerVariableColumns"];

	vartypes = StringJoin[Table[If[MemberQ[integerColumns, i], "I", "C"], {i, nvars}]];
	Print[3, "vartypes: ", vartypes];

	
	GUROBISetVariableTypesAndObjectiveVector[data, vartypes, objvec];
	GUROBIAddQuadraticObjectiveMatrix[data, SparseArray[objmat]];
	
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
	(*status = "Solved";*)
	status = GUROBIStringStatus[data];
	If[!StringQ[status], Return[$Failed]];
	status = Optimization`SolutionData`WrapStatus[status];
	(*Print["obj val: ", GUROBIObjectiveValue[data]];
	Print["x: ",GUROBIx[data]];
	Print["slack: ", GUROBISlack[data]];*)

	{status, GUROBI2Data[data, {integerColumns}]}
];

GUROBI2Data[data_, _]["PrimalMinimumValue"] := GUROBIObjectiveValue[data];
(*GUROBI2Data[data_, {integerColumns_}]["PrimalMinimizerVector"] :=  GUROBIx[data];*)
GUROBI2Data[data_, {integerColumns_}]["PrimalMinimizerVector"] := 
Block[{x = GUROBIx[data], res},
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