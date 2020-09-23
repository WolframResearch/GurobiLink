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
GUROBIDataQ::usage = "GUROBIDataQ[data] gives True if expr represents an active instance of an GUROBIData object."
GUROBIDataCreate::usage = "data = GUROBIDataCreate[] creates an instance of an GUROBIData expression."
GUROBIDataExpressions::ussage = "GUROBIDataExpressions[] shows all active GUROBIData expression instances."
GUROBIDataDelete::usage = "GUROBIDataDelete[data] removes an instance of an GUROBIData expression, freeing up memory."

GUROBIEnvironmentCreate::usage = "env = GUROBIEnvironmentCreate[] creates GUROBI environment."
GUROBIEnvironmentDelete::usage = "GUROBIEnvironmentDelete[env] deletes GUROBI environment."

(*TODO:
add statuses to kernel
*) 

Begin["`Private`"]
(* Implementation of the package *)

$GUROBILinkDirectory = DirectoryName[$InputFileName];
$targetDir = FileNameJoin[{$GUROBILinkDirectory, "LibraryResources", $SystemID}];

$GUROBILinkLibrary = Block[{$LibraryPath = $targetDir}, FindLibrary["GUROBILink"]];

$GUROBILibrariesToPreload = Switch[$SystemID,
	"Windows-x86-64",
		Map[Function[FileNameJoin[{$targetDir, #}]], {"GUROBI90.dll"}],
	_,
		{}
]

(*
 Load all the functions from the GUROBILink library
*)

dPrint = Optimization`Debug`OptimizationDebugPrint;
pReset = Optimization`Debug`OptimizationProfileReset;
pPrint = Optimization`Debug`OptimizationProfilePrint;
(* Set print level 5 with Optimization`Debug`SetPrintLevel[5]
   For GUROBILink only prints set:
   GUROBILink`Private`dPrint = Print;*)
(* Set profile print level 5 with Optimization`Debug`SetProfilePrintLevel[5]
   and use Optimization`Debug`OptimizationProfile[...] *)

needInitialization = True;

$GUROBIInfinity = 1.*10^30;


LoadGUROBILink[] :=
Block[{$LibraryPath = $targetDir}, 
	Map[LibraryLoad, $GUROBILibrariesToPreload];
	GUROBICheckLicense0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_CheckLicense", {Integer}, Integer];
	GUROBICheckModel0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_CheckModel", {Integer}, Integer];
	GUROBISetVariableTypesAndObjectiveVector0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetVariableTypesAndObjectiveVector", {Integer, {Integer, 1}, {Real, 1}}, Integer];
	GUROBISetVariableTypesAndBoundsAndObjectiveVector0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetVariableTypesAndBoundsAndObjectiveVector", {Integer, UTF8String, {Real, 1}, {Real, 1}, {Real, 1}}, Integer];
	GUROBIAddQuadraticObjectiveMatrix0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddQuadraticObjectiveMatrix", {Integer, LibraryDataType[SparseArray, Real, 2]}, Integer];
	GUROBIAddLinearConstraint0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddLinearConstraint", {Integer, {Integer, 1}, {Real, 1}, UTF8String, Real}, Integer];
	GUROBIAddLinearConstraint1 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddLinearConstraint1", {Integer, LibraryDataType[SparseArray, Real, 1], UTF8String, Real}, Integer];
	GUROBIAddLinearConstraints0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddLinearConstraints", {Integer, LibraryDataType[SparseArray, Real, 2], UTF8String, {Real, 1}}, Integer];
	GUROBIAddQuadraticConstraint0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddQuadraticConstraint", {Integer, {Integer, 1}, {Real, 1}, {Integer, 1}, {Integer, 1}, {Real, 1}, UTF8String, Real}, Integer];
	GUROBIAddQuadraticConstraint1 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_AddQuadraticConstraint1", {Integer, LibraryDataType[SparseArray, Real, 2], LibraryDataType[SparseArray, Real, 1], UTF8String, Real}, Integer];
	GUROBISetParameters0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetParameters", {Integer, Integer, Real, Integer}, Integer];
	GUROBISetStartingPoint0 =  LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_SetStartingPoint", {Integer, {Real, 1}}, Integer];

	GUROBIOptimize0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_OptimizeModel", {Integer}, Integer];

	GUROBIStatusValue0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_GetStatusValue", {Integer}, Integer];
	GUROBIObjectiveValue0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_GetObjectiveValue", {Integer}, Real];
	GUROBIx0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_Getx", {Integer}, {Real, 1}];
	GUROBISlack0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIData_GetSlack", {Integer}, {Real, 1}];

	GUROBIDataDelete0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIDataMap_delete", {Integer}, Integer];
	GUROBIDataIDList = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIDataMap_retIDList", {}, {Integer, 1}];
	GUROBIEnvironmentDelete0 = LibraryFunctionLoad[$GUROBILinkLibrary, "GUROBIEnvironmentMap_delete", {Integer}, Integer];
	needInitialization = False;
]

LoadGUROBILink[]

(* GUROBIEnvironment related *)
GUROBIEnvironmentID[e_GUROBIEnvironment] := ManagedLibraryExpressionID[e, "GUROBI_environment_instance_manager"];

GUROBIEnvironmentQ[e_GUROBIEnvironment] := ManagedLibraryExpressionQ[e, "GUROBI_environment_instance_manager"];
GUROBIEnvironmentQ[_] := False;

testGUROBIEnvironment[][e_] := testGUROBIEnvironment[GUROBIEnvironment][e];
testGUROBIEnvironment[mhead_Symbol][e_] :=
If[TrueQ[GUROBIEnvironmentQ[e]],
	True,
	Message[MessageName[mhead, "gurobienvinst"], e]; False
];
testGUROBIEnvironment[_][e_] := TrueQ[GUROBIEnvironmentQ[e]];

General::gurobienvinst = "`1` does not represent an active GUROBIEnvironment object.";

GUROBIEnvironmentCreate[] :=
Module[{},
	If[needInitialization, LoadGUROBILink[]];
	CreateManagedLibraryExpression["GUROBI_environment_instance_manager", GUROBIEnvironment]
];

GUROBIEnvironmentDelete[GUROBIEnvironment[id_]?(testGUROBIEnvironment[GUROBIEnvironmentDelete])] := GUROBIEnvironmentDelete0[id];

(* Create one environment for the entire session, if not already created *)
If[!GUROBIEnvironmentQ[env],
	env = GUROBIEnvironmentCreate[];
	If[GUROBIEnvironmentQ[env],
		dPrint[5, env, " was created"],
		dPrint[5, "Failed to create GUROBI environment"]
	];
	,
	dPrint[5, "GUROBI environment is still ", env];
];

(* Check License *)

GUROBICheckLicense[GUROBIEnvironment[id_]?(testGUROBIEnvironment[GUROBICheckLicense])]:=
Module[{error},
	dPrint[1, "Checking for GUROBI license..."];
	error = GUROBICheckLicense0[id];
	If[error === 0,
		dPrint[1, "...................license found."];
		True
		,
		dPrint[1, "Creating GUROBI environment failed with error ", error];
		If[error === 10009,
			Message[GUROBILink::license];
			False,
			$Failed
		]
	]
]

GUROBILink::license = "Cannot find a valid GUROBI license for version 9.0 or greater." 

(* Check for license, just once *)
hasGUROBILicense = TrueQ[GUROBICheckLicense[env]];

(* Register convex methods, only if license was found *)
If[hasGUROBILicense,
	(* Method "GUROBI1" is only for condstraints suppoted by GUROBI -- linear, quadratic and soc membership.
		For general affine SOC constraints 'Ax+b in K', unless A is diagonal, use method "GUROBI",
		or try Method -> {"GUROBI1", "NonConvex" -> 2} *)
	Optimization`MethodFramework`RegisterOptimizationMethod["GUROBI1",
		Association[
			"SolveFunction" -> GUROBISolve1,
			"ObjectiveSupport" -> "Quadratic",
			"ConstraintSupport" -> Association[{"EqualityConstraint" -> "Affine", "NonNegativeCone" -> "Affine",
			"NormCone" -> "Affine", "QuadraticConstraint" -> "Affine"}],
			"MixedIntegerSupport"->True
		]
	];
	(* "GUROBI" method solves problems with linear or quadratic objective and
		linear, quadratic and second order cone affine constraints.
		In order to handle affine SOC constraints it adds new variables y = A.x+b *)
	Optimization`MethodFramework`RegisterOptimizationMethod["GUROBI",
		Association[
			"SolveFunction" -> GUROBISolve2,
			"ObjectiveSupport" -> "Quadratic",
			"ConstraintSupport" -> Association[{"EqualityConstraint"->"Affine", "NonNegativeCone"->"Affine",
			"NormCone"->"Membership", "QuadraticConstraint" -> "Affine"}],
			"MixedIntegerSupport"->True
		]
	];

	, (* else *)
	Clear[env];
];

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

GUROBICheckModel[GUROBIData[id_]?(testGUROBIData[GUROBICheckModel])]:=
Module[{error},
	error = GUROBICheckModel0[id];
	If[!SameQ[error, 0], dPrint[1, "Model creation failed with error ", error];
		If[SameQ[error, 10001], Message[GUROBILink::nomem]];];
	error
];

(* Functions to set up the problem *)

GUROBISetVariableTypesAndObjectiveVector[GUROBIData[id_]?(testGUROBIData[GUROBISetVariableTypesAndObjectiveVector]), intvars_, objvector_]:=
Module[{},
	dPrint[5, "Setting variable types and objective vector..."];
	GUROBISetVariableTypesAndObjectiveVector0[id, intvars, Normal[objvector]]
];

GUROBISetVariableTypesAndBoundsAndObjectiveVector[GUROBIData[id_]?(testGUROBIData[GUROBISetVariableTypesAndBoundsAndObjectiveVector]), vartypes_, lb_, ub_, objvector_]:=
Module[{},
	GUROBISetVariableTypesAndBoundsAndObjectiveVectort0[id, vartypes, lb, ub, objvector]
];

GUROBIAddQuadraticObjectiveMatrix[GUROBIData[id_]?(testGUROBIData[GUROBIAddQuadraticObjectiveMatrix]), Qmat_SparseArray]:=
Module[{QGmat = Qmat/2},
	dPrint[5, " In GUROBIAddQuadraticObjectiveMatrix"];
	(*the GUROBI Q matrix absorbs the 1/2 coeefficient*)
	dPrint[5, "Adding quadratic objective matrix ", QGmat];
	GUROBIAddQuadraticObjectiveMatrix0[id, QGmat]
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
	dPrint[5, xGUROBIAddLinearConstraints0[id, mat, sense, rhs]];
	GUROBIAddLinearConstraints0[id, SparseArray[mat], sense, Normal[rhs]]
];

GUROBIAddQuadraticConstraintIndices[GUROBIData[id_]?(testGUROBIData[GUROBIAddQuadraticConstraint]), linind_, linvals_, quadrow_, quadcol_, quadvals_, sense_, rhs_] :=
Module[{},
	GUROBIAddQuadraticConstraint0[data, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]
]

GUROBIAddQuadraticConstraint[GUROBIData[id_]?(testGUROBIData[GUROBIAddQuadraticOptimizationConstraint]), Q_SparseArray, q_SparseArray, sense_, b_] :=
Module[{},
	(* x^T.Q.x + q.x <= b *)
	GUROBIAddQuadraticConstraint1[id, Q, q, sense, b]
]

GUROBIAddSOCMembershipConstraint[GUROBIData[id_]?(testGUROBIData[GUROBIAddSOCMembershipConstraint]), ind_]:=
Module[{n, linind, linvals, quadrow, quadcol, quadvals, sense, rhs},
	(* x1^2+...+x(n-1)^2-xn^2 <= 0, xn>=0 *)
	dPrint[3, "In GUROBIAddSOCMembershipConstraint"];
	n = Length[ind];
	dPrint[5, xGUROBIAddLinearConstraint0[id, {ind[[-1]]}, {1}, ">", 0]];
	GUROBIAddLinearConstraint0[id, {ind[[-1]]}, {1}, ">", 0];
	linind = {};
	linvals = {};
	quadrow = ind;
	quadcol = ind;
	quadvals = Append[ConstantArray[1, n-1], -1];
	sense = "<";
	rhs = 0;
	dPrint[5, xGUROBIAddQuadraticConstraint0[id, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]];
  	GUROBIAddQuadraticConstraint0[id, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]
]

GUROBIAddSOCAffineConstraint[GUROBIData[id_]?(testGUROBIData[GUROBIAddSOCConstraintAsQuadratic]), A_, b_]:=
Module[{n, Q, q, b1, psd},
	(* This is only accepted for SOC constraint for diagonal A, with non-diagonal A
	one may try to set the "NonConvex" parameter to 2 but this should not be very efficient *)
	dPrint[3, "In GUROBIAddSOCAffineConstraint"];
	(*VectorGreaterEqual[{Ax + b, 0}, {"NormCone", n}] --> x^T.Q.x + q.x <= b *)
	(*
	an.x + bn >= 0,
	||{a1.x + b1, ..., a (n - 1).x + b (n - 1)}|| <= an.x + bn
	Q = a1 a1^T + ... + a (n - 1).a (n - 1)^T - an^2
	b = -(b1^2 + ... + b (n - 1)^2 - bn^2)
	q = 2 (b1*a1 + ... b (n - 1)*a (n - 1) - bn*an) 
	*)
	n = Length[b];
	dPrint[5, xGUROBIAddLinearConstraint1[id, SparseArray[A[[n]]], ">", -b[[n]]]];
	GUROBIAddLinearConstraint1[id, SparseArray[A[[n]]], ">", -b[[n]]];
	Q = Sum[Transpose[{A[[i]]}].{A[[i]]}, {i, 1, n-1}] - Transpose[{A[[n]]}].{A[[n]]};
	b1 = -Sum[b[[i]]^2, {i, 1, n-1}] + b[[n]]^2;
	q = 2*(Sum[b[[i]]*A[[i]], {i, 1, n-1}] - b[[n]]*A[[n]]);
	If[!DiagonalMatrixQ[A] && !PositiveSemidefiniteMatrixQ[Q],
		dPrint[3, "Matrix Q is not positive semidefinite and matrix A is not diagonal."];
		Message[GUROBILink::badmethod];
	];
	dPrint[5, xGUROBIAddQuadraticConstraint1[id, SparseArray[Q], SparseArray[q], "<", b1]];
	GUROBIAddQuadraticConstraint1[id, SparseArray[Q], SparseArray[q], "<", b1]
]

GUROBILink::badmethod = "Method \"GUROBI1\" is not suitable for general affine SOC constraints. \
It is better to use method \"GUROBI\" instead. \
If looking for adventure, try Method -> {\"GUROBI1\", \"NonConvex\" -> 2}."

GUROBISetNumberOfConstraints[GUROBIData[id_]?(testGUROBIData[GUROBISetNumberOfConstraints]), ncons_]:=
Module[{},
	GUROBISetNumberOfConstraints0[id, ncons];
]

GUROBISetParameters[GUROBIData[id_]?(testGUROBIData[GUROBISetNumberOfConstraints]), maxit_, tol_, nonconvex_]:=
Module[{},
	error = GUROBISetParameters0[id, maxit, tol, nonconvex]
]

GUROBISetStartingPoint[GUROBIData[id_]?(testGUROBIData[GUROBISetNumberOfConstraints]), startpt_]:=
Module[{},
	error = GUROBISetStartingPoint0[id, N[startpt]];
]

(* Function to solve the problem *)

INTMAX = 2^31-1;

Options[GUROBIOptimize] = {Method->Automatic, MaxIterations->Automatic, Tolerance->Automatic, 
"NonConvex"->Automatic, "StartingPoint"->Automatic, "Caller"-> Automatic,
	PerformanceGoal:>$PerformanceGoal, WorkingPrecision->MachinePrecision}

GUROBIOptimize[GUROBIData[id_]?(testGUROBIData[GUROBIOptimize]), OptionsPattern[GUROBIOptimize]] :=
Module[{tol, maxiter, nonconvex, mhead, verbose, status, error, data=GUROBIData[id]},
	dPrint[3, "In GUROBIOptimize"];
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
	pReset[5];
	error = GUROBIOptimize0[id];
	pPrint[5, "GUROBI solver"];
	dPrint[1, "GUROBI solver error: ", error];
	status = GUROBIStatusValue0[id];
	dPrint["status: ", status];
	status = GUROBIStringStatus[data];
	dPrint[1, "string status: ", status];
	ReportError[error, status, mhead, maxiter];
	If[error=!=0 || !StringQ[status], Return[$Failed]];
	status
];

(* https://www.gurobi.com/documentation/9.0/refman/optimization_status_codes.html *)
GUROBIStringStatus[data_] :=
	Switch[GUROBIStatusValue[data],
		1, "NoSolutionAvailable",
		2, "Solved",
		3, "PrimalInfeasible",
		4, "DualInfeasible",
		5, "Unbounded",
		6, "ObjectiveCutoffReached",
		7, "MaxIterationsReached",
		8, "NodeLimitReached",
		9, "MaxCPUTimeReached",
		10, "SolutionLimit",
		11, "Interrupted",
		12, "NumericalError",
		13, "Solved/Inaccurate",
		14, "AsynchronousOptimizationInProgress",
		15, "ObjectiveLimitReached",
		_, "Unknown"
	]

ReportError[error_, status_, mhead_, maxiter_] := (
	Switch[error,
		10020, Message[mhead::npsd]
	];
	Switch[status,
		"Solved/Inaccurate", Message[mhead::inac, maxiter],
		"MaxIterationsReached", Message[mhead::maxit, maxiter],
		"PrimalInfeasible", Message[mhead::nsolc],
		"DualInfeasible", Message[mhead::dinfeas],
		"Interrupted", Message[mhead::userstop],
		"MaxCPUTimeReached", Message[mhead::maxcput]
	];)

General::npsd = "The objective matrix is not positive semi-definite."

(* Selectors *)

GUROBIStatusValue[GUROBIData[id_]?(testGUROBIData[GUROBIStatusValues])] := GUROBIStatusValue0[id];

GUROBIx[GUROBIData[id_]?(testGUROBIData[GUROBIx])] := GUROBIx0[id];

GUROBIObjectiveValue[GUROBIData[id_]?(testGUROBIData[GUROBIObjectiveValue])] := GUROBIObjectiveValue0[id];

GUROBISlack[GUROBIData[id_]?(testGUROBIData[GUROBISlack])] := GUROBISlack0[id];

(* Method functions *)

GUROBISolve1[problemData_, pmopts___] :=
Module[{data, objvec, objmat, ncons, coefficients, coneSpecifications, lpos, intvars, error},

	pReset[5];
	dPrint[1, "In GUROBISolve1"];
	dPrint[3, "pmopts: ", pmopts];

	data = GUROBIDataCreate[];
	If[!GUROBIDataQ[data], Return[$Failed]];
	If[!GUROBICheckModel[data], Return[$Failed]];

	objvec = problemData["ObjectiveVector"];

	objmat = problemData["ObjectiveMatrix"];
	intvars = problemData["IntegerVariableColumns"];
	dPrint[3, "intvars: ", intvars];

	error = GUROBISetVariableTypesAndObjectiveVector[data, intvars, objvec];
	If[error=!=0, Return[$Failed]];

	error = GUROBIAddQuadraticObjectiveMatrix[data, SparseArray[objmat]];
	If[error=!=0, Return[$Failed]];
	
	coefficients = problemData["ConstraintCoefficientArrays"];
	coneSpecifications = problemData["ConstraintSpecifications"];
	ncons = Length[coneSpecifications];
	dPrint[3, "coneSpecifications: ", coneSpecifications];

	(* "EqualityConstraint" will always come first, then "NonNegativeCone",
		any number of "NormCone"s, any number of quadratic constraints *)
	lpos = 1;
	If[ncons>=1 && MatchQ[coneSpecifications[[1]], {"EqualityConstraint", _}],
		{b, a} = coefficients[[1]];
		dPrint[5, "eq {a, b} -> ", {a, b}];
		error = GUROBIAddLinearConstraints[data, a, "=", -b];
		If[error=!=0, Return[$Failed]];
		lpos = 2;
	];
	If[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], {"NonNegativeCone", _}],
		{b, a} = coefficients[[lpos]];
		dPrint[5, "ineq {a, b}-> ", {a, b}];
		error = GUROBIAddLinearConstraints[data, a, ">", -b];
		If[error=!=0, Return[$Failed]];
		lpos += 1;
	];
	While[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], {"NormCone", _}],
			{b, a} = coefficients[[lpos]];
			dPrint[5, "norm {a, b}-> ", {a, b}];
			error = GUROBIAddSOCAffineConstraint[data, a, b];
			If[error=!=0, Return[$Failed]];
			lpos += 1;
	];
	While[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], "QuadraticConstraint"],
			{d, c, q} = coefficients[[lpos]];
			dPrint[5, "quad {q, c, d}-> ", {q, c, d}];
			error = GUROBIAddQuadraticConstraint[data, SparseArray[q], SparseArray[c], "<", -d];
			If[error=!=0, Return[$Failed]];
			lpos += 1;
	];

	status = GUROBIOptimize[data, pmopts];
	status = GUROBIStringStatus[data];
	If[!StringQ[status], Return[$Failed]];
	status = Optimization`SolutionData`WrapStatus[status];
	(*Print["obj val: ", GUROBIObjectiveValue[data]];
	Print["x: ",GUROBIx[data]];*)

	dPrint[1, "status: ", status];
	{status, GUROBI1Data[data, {}]}
]

GUROBI1Data[data_, _]["PrimalMinimumValue"] := GUROBIObjectiveValue[data];
GUROBI1Data[data_, _]["PrimalMinimizerVector"] := GUROBIx[data];
GUROBI1Data[data_, _]["Slack"] := Missing["NotAvailable"];
GUROBI1Data[data_, _]["DualMaximumValue"] := Missing["NotAvailable"];
GUROBI1Data[data_, _]["DualityGap"] := Missing["NotAvailable"];
GUROBI1Data[data_, _]["DualMaximizer"] := Missing["NotAvailable"];

GUROBISolve2[problemData_, pmopts___] :=
Module[{a, b, objvec, data, nvars, status, lpos, nextra, norig, integerColumns,
	coefficients, coneSpecifications, coneVariableIndexes},
	(* For method GUROBI *)

	dPrint[1, "In GUROBISolve2"];
	dPrint[3, "pmopts: ", pmopts];
	pReset[5];

	data = GUROBIDataCreate[];
	If[!GUROBIDataQ[data], Return[$Failed]];
	If[!GUROBICheckModel[data], Return[$Failed]];

	pPrint[5, "Setting up GUROBI: create and check data"];
	objvec = problemData["ObjectiveVector"];
	nvars = Length[objvec];
	objmat = problemData["ObjectiveMatrix"];
	nextra = Lookup[problemData, "ExtraColumns", 0];
	norig = nvars - nextra;
	coefficients = problemData["ConstraintCoefficientArrays"];
	coneSpecifications = problemData["ConstraintSpecifications"];
	ncons = Length[coneSpecifications];
 	coneVariableIndexes = problemData["ConeVariableColumns"];

	dPrint[5, "objvec: ", objvec];
	dPrint[5, "coefficients: ", coefficients];
	dPrint[3, "coneSpecifications: ", coneSpecifications];
	dPrint[3, "ncons: ", ncons];
	If[MatchQ[coneVariableIndexes, _Missing], coneVariableIndexes = ConstantArray[None, ncons]];
	dPrint[5, "coneVariableIndexes: ", coneVariableIndexes];

	integerColumns = problemData["IntegerVariableColumns"];
	dPrint[5, "integerColumns: ", integerColumns];
	pPrint[5, "Setting up GUROBI: get problem data"];

	error = GUROBISetVariableTypesAndObjectiveVector[data, integerColumns, objvec];
	If[error=!=0, Return[$Failed]];
	pPrint[5, "Setting up GUROBI: set vars and linear objective"];

	error = GUROBIAddQuadraticObjectiveMatrix[data, SparseArray[objmat]];
	If[error=!=0, Return[$Failed]];
	pPrint[5, "Setting up GUROBI: set quad objective"];

	(* "EqualityConstraint" will always come first, then "NonNegativeCone",
		any number of "NormCone"s, any number of quadratic constraints *)
	lpos = 1;
	If[ncons >= 1 && MatchQ[coneSpecifications[[1]],{"EqualityConstraint", _}],
		{b, a} = coefficients[[1]];
		error = GUROBIAddLinearConstraints[data, SparseArray[a], "=", -b];
		If[error=!=0, Return[$Failed]];
		lpos = 2;
	];

	If[ncons >= lpos && MatchQ[coneSpecifications[[lpos]],{"NonNegativeCone", _}],
		{b, a} = coefficients[[lpos]];
		error = GUROBIAddLinearConstraints[data, SparseArray[a], ">", -b];
		If[error=!=0, Return[$Failed]];
		lpos += 1;
	];
	pPrint[5, "Setting up GUROBI: set linear constraints"];

	While[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], {"NormCone", _}],
			dPrint[5, "NormCone indices-> ", coneVariableIndexes[[lpos]]];
			error = GUROBIAddSOCMembershipConstraint[data, coneVariableIndexes[[lpos]]];
			If[error=!=0, Return[$Failed]];
			lpos += 1;
	];
	pPrint[5, "Setting up GUROBI: set SOC constraints"];

	While[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], "QuadraticConstraint"],
			{d, c, q} = coefficients[[lpos]];
			dPrint[5, "quad {q, c, d}-> ", {q, c, d}];
			error = GUROBIAddQuadraticConstraint[data, SparseArray[q], SparseArray[c], "<", -d];
			If[error=!=0, Return[$Failed]];
			lpos += 1;
	];
	pPrint[5, "Setting up GUROBI: set quadratic constraints"];

	(* GUROBISolve: *)
	status = GUROBIOptimize[data, pmopts];
	status = GUROBIStringStatus[data];
	If[!StringQ[status], Return[$Failed]];
	status = Optimization`SolutionData`WrapStatus[status];
	(*Print["obj val: ", GUROBIObjectiveValue[data]];
	Print["x: ",GUROBIx[data]];*)

	{status, GUROBI2Data[data, {integerColumns}]}
];

GUROBI2Data[data_, _]["PrimalMinimumValue"] := GUROBIObjectiveValue[data];
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