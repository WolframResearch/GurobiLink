
BeginPackage["GurobiLink`"]

$GurobiLinkLibrary::usage  = "$GurobiLinkLibrary is the full path to the loaded GurobiLink library."
$GurobiLinkDirectory::usage = "$GurobiLinkDirectory gives the location of the GurobiLink library."
LoadGurobiLink::usage  = "LoadGurobiLink[] loads the GurobiLink library."
GurobiLink::usage = "GurobiLink is a symbol used for message heads, e.g. in messages triggered by Gurobi error codes."
GurobiTestLicense::usage = "GurobiTestLicense[] returns True if a valid Gurobi license was found."
GurobiSetVariableTypesAndObjectiveVector::usage = "GurobiSetVariableTypesAndObjectiveVector[data, vartypes, objvector] sets objective vector and variable types as a string with \"C\" at positions of continuous variables and \"I\" for integer variables."
GurobiSetVariableTypesAndBoundsAndObjectiveVector::usage = "GurobiSetVariableTypesAndBoundsAndObjectiveVector[data, vartypes, lowerbounds, upperbounds, objvector] sets objective vector, vectors of lower and upper bounds for the variables and variable types as a string with \"C\" at positions of continuous variables and \"I\" for integer variables."
GurobiAddQuadraticObjectiveMatrix::usage = "GurobiAddQuadraticObjectiveMatrix[data, Q] adds a quadratic objective matrix Q that must be sparse, symmetric and positive semi-definite."
GurobiAddLinearConstraint::usage = "GurobiAddLinearConstraint[data, avec, sense, bnum] adds a constraint avec.x sense bnum, where sense is \">\", \"<\" or \"=\"."
GurobiAddLinearConstraints::usage = "GurobiAddLinearConstraints[data, amat, sense, bvec] adds a constraint amat.x sense bvec, where sense is \">\", \"<\" or \"=\"."
GurobiAddLinearConstraintIndices::usage = "GurobiAddLinearConstraintIndices[data, nzindices, nzvalues, sense, bnum] add a constraint avec.x sense bnum, where sense is \">\", \"<\" or \"=\" and avec is described by its nonzero values and their positions as {nzindices, nzvalues}."
GurobiAddQuadraticConstraint::usage = "GurobiAddQuadraticConstraint[data, Qmat, qvec, sense, rhs] adds a constraint 1/2 x^T.Q.x + q.x sense b, where sense is \">\", \"<\" or \"=\"."
GurobiAddQuadraticConstraintIndices::usage = "GurobiAddQuadraticConstraintIndices[data, linind, linvals, quadrow, quadcol, quadvals, sense, rhs] adds a constraint 1/2 x^T.Q.x + q.x sense b, where sense is \">\", \"<\" or \"=\", Q is represented by its nonzero values and their positions {quadrow, quadcol, quadvals} and q is represented by its non-zero values and their positioons {linind, linvals}."
GurobiAddSOCMembershipConstraint::usage = "GurobiAddSOCMembershipConstraint[data, ind]"
GurobiAddSOCAffineConstraint::usage = "GurobiAddSOCAffineConstraint[data, A, b] adds the constaraint VectorGreaterEqual[{Ax + b, 0}, {\"NormCone\", n}]"

GurobiOptimize::usage = "GurobiOptimize[data]"

GurobiStatusValues::usage ="GurobiStatusValues[data]"
Gurobix::usage ="Gurobix[data]"
GurobiObjectiveValue::usage = "GurobiObjectiveValue[data]"
GurobiSlack::usage = "GurobiSlack[data]"

GurobiData::usage = "GurobiData[id] represents an instance of an GurobiData expression created by GurobiDataCreate."
GurobiDataID::usage = "GurobiDataID[data] gives the instance id of an GurobiData expression data."
GurobiDataQ::usage = "GurobiDataQ[data] gives True if expr represents an active instance of an GurobiData object."
GurobiDataCreate::usage = "data = GurobiDataCreate[] creates an instance of an GurobiData expression."
GurobiDataExpressions::ussage = "GurobiDataExpressions[] shows all active GurobiData expression instances."
GurobiDataDelete::usage = "GurobiDataDelete[data] removes an instance of an GurobiData expression, freeing up memory."

GurobiEnvironmentCreate::usage = "env = GurobiEnvironmentCreate[] creates Gurobi environment."
GurobiEnvironmentDelete::usage = "GurobiEnvironmentDelete[env] deletes Gurobi environment."

Begin["`Private`"]
(* Implementation of the package *)

$GurobiLinkDirectory = DirectoryName[$InputFileName];
$targetDir = FileNameJoin[{$GurobiLinkDirectory, "LibraryResources", $SystemID}];

$GurobiLinkLibrary = Block[{$LibraryPath = $targetDir}, FindLibrary["GurobiLink"]];

$GurobiLibrariesToPreload = Switch[$SystemID,
	"Windows-x86-64",
		FileNames["gurobi*.dll", $targetDir],
	_,
		{}
]

(*
 Load all the functions from the GurobiLink library
*)

dPrint = Optimization`Debug`OptimizationDebugPrint;
pReset = Optimization`Debug`OptimizationProfileReset;
pPrint = Optimization`Debug`OptimizationProfilePrint;

needInitialization = True;

$GurobiInfinity = 1.*10^30;


LoadGurobiLink[] :=
Block[{$LibraryPath = $targetDir}, 
	Map[LibraryLoad, $GurobiLibrariesToPreload];
	GurobiCheckLicense0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_CheckLicense", {Integer}, Integer];
	GurobiCheckModel0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_CheckModel", {Integer}, Integer];
	GurobiSetVariableTypesAndObjectiveVector0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_SetVariableTypesAndObjectiveVector", {Integer, {Integer, 1}, {Real, 1, "Constant"}}, Integer];
	GurobiSetVariableTypesAndBoundsAndObjectiveVector0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_SetVariableTypesAndBoundsAndObjectiveVector", {Integer, UTF8String, {Real, 1, "Constant"}, {Real, 1, "Constant"}, {Real, 1, "Constant"}}, Integer];
	GurobiAddQuadraticObjectiveMatrix0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_AddQuadraticObjectiveMatrix", {Integer, {LibraryDataType[SparseArray, Real, 2], "Constant"}}, Integer];
	GurobiAddLinearConstraint0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_AddLinearConstraint", {Integer, {Integer, 1}, {Real, 1, "Constant"}, UTF8String, Real}, Integer];
	GurobiAddLinearConstraint1 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_AddLinearConstraint1", {Integer, {LibraryDataType[SparseArray, Real, 1], "Constant"}, UTF8String, Real}, Integer];
	GurobiAddLinearConstraints0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_AddLinearConstraints", {Integer, {LibraryDataType[SparseArray, Real, 2], "Constant"}, UTF8String, {Real, 1, "Constant"}}, Integer];
	GurobiAddQuadraticConstraint0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_AddQuadraticConstraint", {Integer, {Integer, 1}, {Real, 1, "Constant"}, {Integer, 1}, {Integer, 1}, {Real, 1, "Constant"}, UTF8String, Real}, Integer];
	GurobiAddQuadraticConstraint1 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_AddQuadraticConstraint1", {Integer, {LibraryDataType[SparseArray, Real, 2], "Constant"}, {LibraryDataType[SparseArray, Real, 1], "Constant"}, UTF8String, Real}, Integer];
	GurobiSetParameters0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_SetParameters", {Integer, Integer, Real, Integer}, Integer];
	GurobiSetStartingPoint0 =  LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_SetStartingPoint", {Integer, {Real, 1, "Constant"}}, Integer];

	GurobiOptimize0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_OptimizeModel", {Integer}, Integer];

	GurobiStatusValue0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_GetStatusValue", {Integer}, Integer];
	GurobiObjectiveValue0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_GetObjectiveValue", {Integer}, Real];
	Gurobix0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_Getx", {Integer}, {Real, 1}];
	GurobiSlack0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiData_GetSlack", {Integer}, {Real, 1}];

	GurobiDataDelete0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiDataMap_delete", {Integer}, Integer];
	GurobiDataIDList = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiDataMap_retIDList", {}, {Integer, 1}];
	GurobiEnvironmentDelete0 = LibraryFunctionLoad[$GurobiLinkLibrary, "GurobiEnvironmentMap_delete", {Integer}, Integer];
	needInitialization = False;
]

LoadGurobiLink[]

(* GurobiEnvironment related *)
GurobiEnvironmentID[e_GurobiEnvironment] := ManagedLibraryExpressionID[e, "Gurobi_environment_instance_manager"];

GurobiEnvironmentQ[e_GurobiEnvironment] := ManagedLibraryExpressionQ[e, "Gurobi_environment_instance_manager"];
GurobiEnvironmentQ[_] := False;

testGurobiEnvironment[][e_] := testGurobiEnvironment[GurobiEnvironment][e];
testGurobiEnvironment[mhead_Symbol][e_] :=
If[TrueQ[GurobiEnvironmentQ[e]],
	True,
	Message[MessageName[mhead, "gurobienvinst"], e]; False
];
testGurobiEnvironment[_][e_] := TrueQ[GurobiEnvironmentQ[e]];

General::gurobienvinst = "`1` does not represent an active GurobiEnvironment object.";

GurobiEnvironmentCreate[] :=
Module[{},
	If[needInitialization, LoadGurobiLink[]];
	CreateManagedLibraryExpression["Gurobi_environment_instance_manager", GurobiEnvironment]
];

GurobiEnvironmentDelete[GurobiEnvironment[id_]?(testGurobiEnvironment[GurobiEnvironmentDelete])] := GurobiEnvironmentDelete0[id];

(* Create one environment for the entire session, if not already created *)
If[!GurobiEnvironmentQ[env],
	env = GurobiEnvironmentCreate[];
	If[GurobiEnvironmentQ[env],
		dPrint[5, env, " was created"],
		dPrint[5, "Failed to create Gurobi environment"]
	];
	,
	dPrint[5, "Gurobi environment is still ", env];
];

(* Check License *)

GurobiCheckLicense[GurobiEnvironment[id_]?(testGurobiEnvironment[GurobiCheckLicense])]:=
Module[{error},
	dPrint[1, "Checking for Gurobi license..."];
	error = GurobiCheckLicense0[id];
	If[error === 0,
		dPrint[1, "...................license found."];
		True
		,
		dPrint[1, "Creating Gurobi environment failed with error ", error];
		If[error === 10009,
			Message[GurobiLink::license];
			False,
			$Failed
		]
	]
]

GurobiLink::license = "Cannot find a valid Gurobi license for version 9.0 or greater."

(* Check for license, just once *)
GurobiTestLicense[] :=
	(GurobiTestLicense[] = TrueQ[GurobiCheckLicense[env]])

licenseData = <|"TestFunction"->GurobiTestLicense, "WorkflowName"->"Gurobi"|>;

(* Register optimization methods *)

(* Method "Gurobi1" is only for condstraints suppoted by Gurobi -- linear, quadratic and second order cone(SOC) membership.
   For general affine SOC constraints 'Ax+b in K', unless A is diagonal, use method "Gurobi",
   or try Method -> {"Gurobi1", "NonConvex" -> 2} *)
Optimization`MethodFramework`RegisterOptimizationMethod["Gurobi1",
	Association[
		"SolveFunction" -> GurobiSolve1,
		"ObjectiveSupport" -> "Quadratic",
		"ConstraintSupport" -> Association[{"EqualityConstraint" -> "Affine", "NonNegativeCone" -> "Affine",
		"NormCone" -> "Affine", "QuadraticConstraint" -> "Affine"}],
		"MixedIntegerSupport"->True,
		"License"->licenseData
	]
];
(* Method "Gurobi" solves problems with linear or quadratic objective and
    linear, quadratic and second order cone(SOC) affine constraints.
    In order to handle affine SOC constraints it adds new variables y = A.x+b *)
Optimization`MethodFramework`RegisterOptimizationMethod["Gurobi",
	Association[
	"SolveFunction" -> GurobiSolve2,
		"ObjectiveSupport" -> "Quadratic",
		"ConstraintSupport" -> Association[{"EqualityConstraint"->"Affine", "NonNegativeCone"->"Affine",
		"NormCone"->"Membership", "QuadraticConstraint" -> "Affine"}],
		"MixedIntegerSupport"->True,
		"License"->licenseData
	]
];

Optimization`MethodFramework`RegisterOptimizationMethodAliases["Gurobi", {"Gurobi"}];

(* GurobiData expression (GurobiSolMap) related: *)

GurobiDataID[e_GurobiData] := ManagedLibraryExpressionID[e, "Gurobi_data_instance_manager"];

GurobiDataQ[e_GurobiData] := ManagedLibraryExpressionQ[e, "Gurobi_data_instance_manager"];
GurobiDataQ[_] := False;

testGurobiData[][e_] := testGurobiData[GurobiData][e];
testGurobiData[mhead_Symbol][e_] :=
If[TrueQ[GurobiDataQ[e]],
	True,
	Message[MessageName[mhead, "Gurobiinst"], e]; False
];
testGurobiData[_][e_] := TrueQ[GurobiDataQ[e]];

General::gurobiinst = "`1` does not represent an active GurobiData object.";

GurobiDataCreate[] :=
Module[{},
	If[needInitialization, LoadGurobiLink[]];
	CreateManagedLibraryExpression["Gurobi_data_instance_manager", GurobiData]
];

GurobiDataDelete[GurobiData[id_]?(testGurobiData[GurobiDataDelete])] := GurobiDataDelete0[id];

GurobiDataDelete[l:{_GurobiData..}] := GurobiDataDelete /@ l;

GurobiDataExpressions[] :=
Module[{list},
	If[needInitialization, LoadGurobi[]];
	list = GurobiDataIDList[];
	If[!ListQ[list],
	   $Failed,
	   Map[GurobiData, list]]
]

GurobiCheckModel[GurobiData[id_]?(testGurobiData[GurobiCheckModel])]:=
Module[{error},
	error = GurobiCheckModel0[id];
	If[!SameQ[error, 0], dPrint[1, "Model creation failed with error ", error];
		If[SameQ[error, 10001], Message[GurobiLink::nomem]];];
	error
];

(* Functions to set up the problem *)

GurobiSetVariableTypesAndObjectiveVector[GurobiData[id_]?(testGurobiData[GurobiSetVariableTypesAndObjectiveVector]), intvars_, objvector_]:=
Module[{},
	dPrint[5, "Setting variable types and objective vector..."];
	GurobiSetVariableTypesAndObjectiveVector0[id, intvars, Normal[objvector]]
];

GurobiSetVariableTypesAndBoundsAndObjectiveVector[GurobiData[id_]?(testGurobiData[GurobiSetVariableTypesAndBoundsAndObjectiveVector]), vartypes_, lb_, ub_, objvector_]:=
Module[{},
	GurobiSetVariableTypesAndBoundsAndObjectiveVectort0[id, vartypes, lb, ub, objvector]
];

GurobiAddQuadraticObjectiveMatrix[GurobiData[id_]?(testGurobiData[GurobiAddQuadraticObjectiveMatrix]), Qmat_SparseArray]:=
Module[{QGmat = Qmat/2},
	dPrint[5, " In GurobiAddQuadraticObjectiveMatrix"];
	(*the Gurobi Q matrix absorbs the 1/2 coeefficient*)
	dPrint[5, "Adding quadratic objective matrix ", QGmat];
	GurobiAddQuadraticObjectiveMatrix0[id, QGmat]
];

GurobiAddLinearConstraintIndices[GurobiData[id_]?(testGurobiData[GurobiAddLinearConstraintIndices]), indices_, values_, sense_, rhs_]:=
Module[{},
	GurobiAddLinearConstraint0[id, indices, values, sense, rhs]
];

GurobiAddLinearConstraint[GurobiData[id_]?(testGurobiData[GurobiAddLinearConstraint]), vector_SparseArray, sense_, rhs_]:=
Module[{},
	GurobiAddLinearConstraint1[id, vector, sense, rhs]
];

GurobiAddLinearConstraints[GurobiData[id_]?(testGurobiData[GurobiAddLinearConstraints]), mat_SparseArray, sense_, rhs_]:=
Module[{},
	dPrint[5, xGurobiAddLinearConstraints0[id, mat, sense, rhs]];
	GurobiAddLinearConstraints0[id, SparseArray[mat], sense, Normal[rhs]]
];

GurobiAddQuadraticConstraintIndices[GurobiData[id_]?(testGurobiData[GurobiAddQuadraticConstraintIndices]), linind_, linvals_, quadrow_, quadcol_, quadvals_, sense_, rhs_] :=
Module[{},
	GurobiAddQuadraticConstraint0[id, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]
]

GurobiAddQuadraticConstraint[GurobiData[id_]?(testGurobiData[GurobiAddQuadraticOptimizationConstraint]), Q_SparseArray, q_SparseArray, sense_, b_] :=
Module[{},
	(* x^T.Q.x + q.x <= b *)
	GurobiAddQuadraticConstraint1[id, Q, q, sense, b]
]

GurobiAddSOCMembershipConstraint[GurobiData[id_]?(testGurobiData[GurobiAddSOCMembershipConstraint]), ind_]:=
Module[{n, linind, linvals, quadrow, quadcol, quadvals, sense, rhs},
	(* x1^2+...+x(n-1)^2-xn^2 <= 0, xn>=0 *)
	dPrint[3, "In GurobiAddSOCMembershipConstraint"];
	n = Length[ind];
	dPrint[5, xGurobiAddLinearConstraint0[id, {ind[[-1]]}, {1}, ">", 0]];
	GurobiAddLinearConstraint0[id, {ind[[-1]]}, {1}, ">", 0];
	linind = {};
	linvals = {};
	quadrow = ind;
	quadcol = ind;
	quadvals = Append[ConstantArray[1, n-1], -1];
	sense = "<";
	rhs = 0;
	dPrint[5, xGurobiAddQuadraticConstraint0[id, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]];
  	GurobiAddQuadraticConstraint0[id, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]
]

GurobiAddSOCAffineConstraint[GurobiData[id_]?(testGurobiData[GurobiAddSOCAffineConstraint]), A_, b_]:=
Module[{n, Q, q, b1, psd},
	(* This is only accepted for SOC constraint for diagonal A, with non-diagonal A
	one may try to set the "NonConvex" parameter to 2 but this should not be very efficient *)
	dPrint[3, "In GurobiAddSOCAffineConstraint"];
	(*VectorGreaterEqual[{Ax + b, 0}, {"NormCone", n}] --> x^T.Q.x + q.x <= b *)
	(*
	an.x + bn >= 0,
	||{a1.x + b1, ..., a (n - 1).x + b (n - 1)}|| <= an.x + bn
	Q = a1 a1^T + ... + a (n - 1).a (n - 1)^T - an^2
	b = -(b1^2 + ... + b (n - 1)^2 - bn^2)
	q = 2 (b1*a1 + ... b (n - 1)*a (n - 1) - bn*an) 
	*)
	n = Length[b];
	dPrint[5, xGurobiAddLinearConstraint1[id, SparseArray[A[[n]]], ">", -b[[n]]]];
	GurobiAddLinearConstraint1[id, SparseArray[A[[n]]], ">", -b[[n]]];
	Q = Sum[Transpose[{A[[i]]}].{A[[i]]}, {i, 1, n-1}] - Transpose[{A[[n]]}].{A[[n]]};
	b1 = -Sum[b[[i]]^2, {i, 1, n-1}] + b[[n]]^2;
	q = 2*(Sum[b[[i]]*A[[i]], {i, 1, n-1}] - b[[n]]*A[[n]]);
	If[!DiagonalMatrixQ[A] && !PositiveSemidefiniteMatrixQ[Q],
		dPrint[3, "Matrix Q is not positive semidefinite and matrix A is not diagonal."];
		Message[GurobiLink::badmethod];
	];
	dPrint[5, xGurobiAddQuadraticConstraint1[id, SparseArray[Q], SparseArray[q], "<", b1]];
	GurobiAddQuadraticConstraint1[id, SparseArray[Q], SparseArray[q], "<", b1]
]

GurobiLink::badmethod = "Method \"Gurobi1\" is not suitable for general affine SOC constraints. \
Use method \"Gurobi\"."

GurobiSetParameters[GurobiData[id_]?(testGurobiData[GurobiSetParameters]), maxit_, tol_, nonconvex_]:=
Module[{error},
	error = GurobiSetParameters0[id, maxit, tol, nonconvex]
]

GurobiSetStartingPoint[GurobiData[id_]?(testGurobiData[GurobiSetStartingPoint]), startpt_]:=
Module[{error},
	error = GurobiSetStartingPoint0[id, N[startpt]];
]

(* Function to solve the problem *)

INTMAX = 2^31-1;

Options[GurobiOptimize] = {Method->Automatic, MaxIterations->Automatic, Tolerance->Automatic,
"NonConvex"->Automatic, "StartingPoint"->Automatic, "Caller"-> Automatic,
	PerformanceGoal:>$PerformanceGoal, WorkingPrecision->MachinePrecision}

GurobiOptimize[GurobiData[id_]?(testGurobiData[GurobiOptimize]), OptionsPattern[GurobiOptimize]] :=
Module[{tol, maxiter, nonconvex, mhead, verbose, status, error, startpt, data=GurobiData[id]},
	dPrint[3, "In GurobiOptimize"];
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
		error = GurobiSetStartingPoint0[id, N[startpt]];
	];

	dPrint[3, "Setting parameters {maxiter, tol, nonconvex} -> ", {maxiter, tol, nonconvex}];
	error = GurobiSetParameters0[id, maxiter, tol, nonconvex];

	dPrint[1, "Solving with Gurobi..."];
	pReset[5];
	error = GurobiOptimize0[id];
	pPrint[5, "Gurobi solver"];
	dPrint[1, "Gurobi solver error: ", error];
	status = GurobiStatusValue0[id];
	dPrint["status: ", status];
	status = GurobiStringStatus[data];
	dPrint[1, "string status: ", status];
	ReportError[error, status, mhead, maxiter];
	If[error=!=0 || !StringQ[status], Return[$Failed]];
	status
];

(* https://www.gurobi.com/documentation/9.0/refman/optimization_status_codes.html *)
GurobiStringStatus[data_] :=
	Switch[GurobiStatusValue[data],
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

GurobiStatusValue[GurobiData[id_]?(testGurobiData[GurobiStatusValues])] := GurobiStatusValue0[id];

Gurobix[GurobiData[id_]?(testGurobiData[Gurobix])] := Gurobix0[id];

GurobiObjectiveValue[GurobiData[id_]?(testGurobiData[GurobiObjectiveValue])] := GurobiObjectiveValue0[id];

GurobiSlack[GurobiData[id_]?(testGurobiData[GurobiSlack])] := GurobiSlack0[id];

(* Determine the NonConvex setting for the Gurobi optimizer:
ttps://www.gurobi.com/documentation/9.0/refman/nonconvex.html*)

getNonConvexSetting[problemData_] :=
Module[{convexity = problemData["Convexity"]},
	dPrint[5, "convexity: ", convexity];
	If[MatchQ[convexity, "Convex"], 1, 2]
]

(* Method functions *)

(* Solve function for method Gurobi1 *)
GurobiSolve1[problemData_, pmopts___] :=
Module[{a, b, c, d, q, data, objvec, objmat, ncons, coefficients, coneSpecifications, lpos, intvars, nonconvex, error, status},
	pReset[5];
	dPrint[1, "In GurobiSolve1"];
	dPrint[3, "pmopts: ", pmopts];

	(* Retrieve and set up the problem data *)
	data = GurobiDataCreate[];
	If[!GurobiDataQ[data], Return[$Failed]];
	If[!GurobiCheckModel[data], Return[$Failed]];

	objvec = problemData["ObjectiveVector"];

	objmat = problemData["ObjectiveMatrix"];
	intvars = problemData["IntegerVariableColumns"];
	dPrint[3, "intvars: ", intvars];

	error = GurobiSetVariableTypesAndObjectiveVector[data, intvars, objvec];
	If[error=!=0, Return[$Failed]];

	error = GurobiAddQuadraticObjectiveMatrix[data, SparseArray[objmat]];
	If[error=!=0, Return[$Failed]];
	
	coefficients = problemData["ConstraintCoefficientArrays"];
	coneSpecifications = problemData["ConstraintSpecifications"];
	ncons = Length[coneSpecifications];
	dPrint[3, "coneSpecifications: ", coneSpecifications];

	nonconvex = getNonConvexSetting[problemData];
	dPrint[5, "nonconvex: ", nonconvex];

	(* "EqualityConstraint" will always come first, then "NonNegativeCone",
		any number of "NormCone"s, any number of quadratic constraints *)
	lpos = 1;
	If[ncons>=1 && MatchQ[coneSpecifications[[1]], {"EqualityConstraint", _}],
		{b, a} = coefficients[[1]];
		dPrint[5, "eq {a, b} -> ", {a, b}];
		error = GurobiAddLinearConstraints[data, a, "=", -b];
		If[error=!=0, Return[$Failed]];
		lpos = 2;
	];
	If[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], {"NonNegativeCone", _}],
		{b, a} = coefficients[[lpos]];
		dPrint[5, "ineq {a, b}-> ", {a, b}];
		error = GurobiAddLinearConstraints[data, a, ">", -b];
		If[error=!=0, Return[$Failed]];
		lpos += 1;
	];
	While[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], {"NormCone", _}],
			{b, a} = coefficients[[lpos]];
			dPrint[5, "norm {a, b}-> ", {a, b}];
			error = GurobiAddSOCAffineConstraint[data, a, b];
			If[error=!=0, Return[$Failed]];
			lpos += 1;
	];
	While[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], "QuadraticConstraint"],
			{d, c, q} = coefficients[[lpos]];
			dPrint[5, "quad {q, c, d}-> ", {q, c, d}];
			error = GurobiAddQuadraticConstraint[data, SparseArray[q], SparseArray[c], "<", -d];
			If[error=!=0, Return[$Failed]];
			lpos += 1;
	];

	(* Solve *)
	status = GurobiOptimize[data, {pmopts, "NonConvex"->nonconvex}];

	status = GurobiStringStatus[data];
	If[!StringQ[status], Return[$Failed]];
	status = Optimization`SolutionData`WrapStatus[status];
	(*Print["obj val: ", GurobiObjectiveValue[data]];
	Print["x: ", Gurobix[data]];*)

	dPrint[1, "status: ", status];
	{status, Gurobi1Data[data, {}]}
]

Gurobi1Data[data_, _]["PrimalMinimumValue"] := GurobiObjectiveValue[data];
Gurobi1Data[data_, _]["PrimalMinimizerVector"] := Gurobix[data];
Gurobi1Data[data_, _]["Slack"] := Missing["NotAvailable"];
Gurobi1Data[data_, _]["DualMaximumValue"] := Missing["NotAvailable"];
Gurobi1Data[data_, _]["DualityGap"] := Missing["NotAvailable"];
Gurobi1Data[data_, _]["DualMaximizer"] := Missing["NotAvailable"];

(* Solve function for method Gurobi *)
GurobiSolve2[problemData_, pmopts___] :=
Module[{a, b, c, d, q, objvec, objmat, data, nvars, status, lpos, nextra, norig, integerColumns,
	coefficients, coeff, coneSpecifications, ncons, coneVariableIndexes, nonconvex, error},

	dPrint[1, "In GurobiSolve2"];
	dPrint[3, "pmopts: ", pmopts];
	pReset[5];

	(* Retrieve and set up the problem data: *)
	data = GurobiDataCreate[];
	If[!GurobiDataQ[data], Return[$Failed]];
	If[!GurobiCheckModel[data], Return[$Failed]];

	pPrint[5, "Setting up Gurobi: create and check data"];
	objvec = problemData["ObjectiveVector"];
	dPrint[5, "objvec: ", objvec];
	nvars = Length[objvec];
	objmat = problemData["ObjectiveMatrix"];
	nextra = Lookup[problemData, "ExtraColumns", 0];
	norig = nvars - nextra;
	coneSpecifications = problemData["ConstraintSpecifications"];
	dPrint[3, "coneSpecifications: ", coneSpecifications];
	coefficients = problemData["ConstraintCoefficientArrays"];
	dPrint[5, "coefficients: ", coefficients];
	ncons = Length[coneSpecifications];
	dPrint[3, "ncons: ", ncons];
	coneVariableIndexes = problemData["ConeVariableColumns"];
	If[MatchQ[coneVariableIndexes, _Missing], coneVariableIndexes = ConstantArray[None, ncons]];
	dPrint[5, "coneVariableIndexes: ", coneVariableIndexes];
	nonconvex = getNonConvexSetting[problemData];
	dPrint[5, "nonconvex: ", nonconvex];
	integerColumns = problemData["IntegerVariableColumns"];
	dPrint[5, "integerColumns: ", integerColumns];
	pPrint[5, "Setting up Gurobi: get problem data"];

	error = GurobiSetVariableTypesAndObjectiveVector[data, integerColumns, objvec];
	If[error=!=0, Return[$Failed]];
	pPrint[5, "Setting up Gurobi: set vars and linear objective"];

	error = GurobiAddQuadraticObjectiveMatrix[data, SparseArray[objmat]];
	If[error=!=0, Return[$Failed]];
	pPrint[5, "Setting up Gurobi: set quad objective"];

	(* "EqualityConstraint" will always come first, then "NonNegativeCone",
		any number of "NormCone"s, any number of quadratic constraints *)
	lpos = 1;
	If[ncons >= 1 && MatchQ[coneSpecifications[[1]], {"EqualityConstraint", _}],
		coeff = coefficients[[lpos]];
		If[ListQ[coeff] && Length[coeff] === 2,
			{b, a} = coefficients[[1]],
			Return[$Failed]
		];
		error = GurobiAddLinearConstraints[data, SparseArray[a], "=", -b];
		If[error=!=0, dPrint[3, "add lin = error: ", error]; Return[$Failed]];
		lpos = 2;
	];
	If[ncons >= lpos && MatchQ[coneSpecifications[[lpos]], {"NonNegativeCone", _}],
		coeff = coefficients[[lpos]];
		If[ListQ[coeff] && Length[coeff] === 2,
			{b, a} = coefficients[[lpos]],
			Return[$Failed]
		];
		error = GurobiAddLinearConstraints[data, SparseArray[a], ">", -b];
		If[error =!= 0, dPrint[3, "add lin > error: ", error]; Return[$Failed]];
		lpos += 1;
	];
	pPrint[5, "Setting up Gurobi: set linear constraints"];

	While[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], {"NormCone", _}],
			dPrint[5, "NormCone indices-> ", coneVariableIndexes[[lpos]]];
			error = GurobiAddSOCMembershipConstraint[data, coneVariableIndexes[[lpos]]];
			If[error=!=0, dPrint[3, "add soc const error: ", error]; Return[$Failed]];
			lpos += 1;
	];
	pPrint[5, "Setting up Gurobi: set SOC constraints"];

	While[ncons>=lpos && MatchQ[coneSpecifications[[lpos]], "QuadraticConstraint"],
			{d, c, q} = coefficients[[lpos]];
			dPrint[5, "quad {q, c, d}-> ", {q, c, d}];
			error = GurobiAddQuadraticConstraint[data, SparseArray[q], SparseArray[c], "<", -d];
			If[error=!=0, dPrint[3, "add quad const error: ", error]; Return[$Failed]];
			lpos += 1;
	];
	pPrint[5, "Setting up Gurobi: set quadratic constraints"];

	(* Solve: *)
	status = GurobiOptimize[data, {pmopts, "NonConvex"->nonconvex}];
	status = GurobiStringStatus[data];
	If[!StringQ[status], Return[$Failed]];
	status = Optimization`SolutionData`WrapStatus[status];
	(*Print["obj val: ", GurobiObjectiveValue[data]];
	Print["x: ", Gurobix[data]];*)

	{status, Gurobi2Data[data, {integerColumns}]}
];

Gurobi2Data[data_, _]["PrimalMinimumValue"] := GurobiObjectiveValue[data];
Gurobi2Data[data_, {integerColumns_}]["PrimalMinimizerVector"] :=
Block[{x = Gurobix[data]},
	If[Length[integerColumns] > 0,
		x[[integerColumns]] = Round[x[[integerColumns]]]
	];
	x
];
Gurobi2Data[data_, _]["DualMaximumValue"] := Missing["NotAvailable"];
Gurobi2Data[data_, _]["DualityGap"] := Missing["NotAvailable"];
Gurobi2Data[data_, _]["DualMaximizer"] := Missing["NotAvailable"];
Gurobi2Data[data_, _]["Slack"] := Missing["NotAvailable"];

End[]
EndPackage[]
