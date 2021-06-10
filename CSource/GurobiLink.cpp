#include <algorithm>
#include <cstring>
#include "GurobiSolution.h"

EXTERN_C DLLEXPORT mint WolframLibrary_getVersion()
{
	return WolframLibraryVersion;
}

EXTERN_C DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData)
{

	int err;
	err = (*libData->registerLibraryExpressionManager)("Gurobi_data_instance_manager", GurobiDataMap_manage);
	if (err)
		return err;
	err = (*libData->registerLibraryExpressionManager)("Gurobi_environment_instance_manager", GurobiEnvironmentMap_manage);

	return err;
}

EXTERN_C DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData)
{
	// remove environment
	GurobiEnvironmentMap_manage(libData, 1, 1);
	(*libData->unregisterLibraryCallbackManager)("Gurobi_data_instance_manager");
	(*libData->unregisterLibraryCallbackManager)("Gurobi_environment_instance_manager");
}

/************************************************************************/
/*                GurobiData_CheckLicense                               */
/*                                                                      */
/*             GurobiCheckLicense[environment]                          */
/************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_CheckLicense(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	mint envID;
	GurobiEnvironment Gurobienvironment;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}

	envID = MArgument_getInteger(Args[0]);
	Gurobienvironment = GurobiEnvironmentMap_get(envID);

	// load return values
	MArgument_setInteger(Res, Gurobienvironment->error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GurobiData_CheckModel                                 */
/*                                                                      */
/*                GurobiCheckModel[data]                                */
/************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_CheckModel(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	mint dataID;
	GurobiData Gurobidata;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}

	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	// return the error from GRBnewmodel in GurobiEData_new
	MArgument_setInteger(Res, (mint)(Gurobidata->error));

	return LIBRARY_NO_ERROR;
}

/******************************************************************************/
/*             GurobiData_SetVariableTypesAndObjectiveVector                  */
/*                                                                            */
/*     GurobiSetVariableTypesAndObjectiveVector[data, intvars, objvector]     */
/******************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_SetVariableTypesAndObjectiveVector(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint i, nvars, nintvars, dataID, *intvars = nullptr;
	char* vartypes;
	double *objvec, *lbounds = nullptr;
	MTensor pT;
	GurobiData Gurobidata;

	if (Argc != 3)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	pT = MArgument_getMTensor(Args[1]);
	nintvars = libData->MTensor_getFlattenedLength(pT);
	intvars = libData->MTensor_getIntegerData(pT);

	pT = MArgument_getMTensor(Args[2]);
	nvars = libData->MTensor_getFlattenedLength(pT);
	objvec = libData->MTensor_getRealData(pT);
	Gurobidata->nvars = nvars;

	lbounds = (double*)malloc(nvars * sizeof(double));
	// set lb to be a vector of -GRB_INFINITY; if we leave it nullptr the default is a vector of 0.0
	std::fill_n(lbounds, nvars, -GRB_INFINITY);

	vartypes = (char*)malloc(nvars * sizeof(char));
	std::memset(vartypes, 'C', nvars);
	for (i = 0; i < nintvars; i++)
	{
		vartypes[intvars[i] - 1] = 'I';
	}

	// Add objective vector and variable types
	error = GRBaddvars(Gurobidata->model, (int)nvars, 0, nullptr, nullptr, nullptr, objvec, lbounds, nullptr, vartypes, nullptr);

	// free stuff
	if (lbounds)
		free(lbounds);
	if (vartypes)
		free(vartypes);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/*****************************************************************************************************/
/*             GurobiData_SetVariableTypesAndBoundsAndObjectiveVector                                */
/*                                                                                                   */
/*  GurobiSetVariableTypesAndBoundsAndObjectiveVector[data, vartypes, lbounds, ubounds, objvector]   */
/*****************************************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_SetVariableTypesAndBoundsAndObjectiveVector(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, nvars;
	char* vartypes = nullptr;
	double *objvec, *lbounds, *ubounds;
	mint dataID;
	MTensor pT;
	GurobiData Gurobidata;

	if (Argc != 5)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	vartypes = MArgument_getUTF8String(Args[1]);
	nvars = (int)strlen(vartypes);

	Gurobidata->nvars = nvars;

	pT = MArgument_getMTensor(Args[2]);
	lbounds = libData->MTensor_getRealData(pT);

	pT = MArgument_getMTensor(Args[3]);
	ubounds = libData->MTensor_getRealData(pT);

	pT = MArgument_getMTensor(Args[4]);
	objvec = libData->MTensor_getRealData(pT);

	// Add objective vector and variable types
	error = GRBaddvars(Gurobidata->model, nvars, 0, nullptr, nullptr, nullptr, objvec, lbounds, ubounds, vartypes, nullptr);

	// free stuff
	libData->UTF8String_disown(vartypes);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/**********************************************************************/
/*                GurobiData_AddQuadraticObjectiveMatrix              */
/*                                                                    */
/*             GurobiAddQuadraticObjectiveMatrix[data, Qmat]          */
/* adds the term 1/2 x^T.Q.x to already specified lin objective  c.x  */
/**********************************************************************/

EXTERN_C DLLEXPORT int GurobiData_AddQuadraticObjectiveMatrix(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error = 0, numquadnz, *quadrow = nullptr, *quadcol = nullptr;
	double *quadval = nullptr, *explicitValues = nullptr;
	mint i, err, dataID, *explicitPositions;
	MTensor *pT = nullptr, TQ = nullptr;
	GurobiData Gurobidata;
	MSparseArray Qmat;

	if (Argc != 2)
		return LIBRARY_FUNCTION_ERROR;
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	Qmat = MArgument_getMSparseArray(Args[1]);

	pT = (*(libData->sparseLibraryFunctions->MSparseArray_getExplicitValues))(Qmat);
	quadval = libData->MTensor_getRealData(*pT);
	if (*pT != nullptr)
		numquadnz = (int)(libData->MTensor_getFlattenedLength(*pT));
	else
		numquadnz = 0;

	if (numquadnz > 0)
	{
		err = (*(libData->sparseLibraryFunctions->MSparseArray_getExplicitPositions))(Qmat, &TQ);
		if (err)
			return LIBRARY_FUNCTION_ERROR;
		explicitPositions = libData->MTensor_getIntegerData(TQ);

		quadrow = (int*)malloc(numquadnz * sizeof(int));
		quadcol = (int*)malloc(numquadnz * sizeof(int));
		for (i = 0; i < numquadnz; i++)
		{
			quadrow[i] = (int)explicitPositions[2 * i] - 1;
			quadcol[i] = (int)explicitPositions[2 * i + 1] - 1;
		}

		error = GRBaddqpterms(Gurobidata->model, numquadnz, quadrow, quadcol, quadval);

		if (quadrow)
			free(quadrow);
		if (quadcol)
			free(quadcol);
	}

	// free more stuff
	if (TQ)
		libData->MTensor_free(TQ);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/*****************************************************************************/
/*                GurobiData_AddLinearConstraint                             */
/*                                                                           */
/*    GurobiAddLinearConstraint[data, indices, values, sense, rhs]           */
/* sense: GRB_EQUAL -> '=', GRB_LESS_EQUAL -> '<', GRB_GREATER_EQUAL -> '>'  */
/*****************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_AddLinearConstraint(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int j, error, nz, *indices = nullptr;
	double *values = nullptr, rhs;
	mint dataID, *idata = nullptr;
	MTensor pT;
	GurobiData Gurobidata;
	char sense, *senseString;

	if (Argc != 5)
	{
		return LIBRARY_FUNCTION_ERROR;
	}

	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	pT = MArgument_getMTensor(Args[1]);
	nz = (int)(libData->MTensor_getFlattenedLength(pT));
	idata = libData->MTensor_getIntegerData(pT);
	indices = (int*)malloc(nz * sizeof(int));
	for (j = 0; j < nz; j++)
	{
		indices[j] = (int)idata[j] - 1;
	}
	pT = MArgument_getMTensor(Args[2]);
	values = libData->MTensor_getRealData(pT);

	senseString = MArgument_getUTF8String(Args[3]);
	sense = senseString[0];

	rhs = MArgument_getReal(Args[4]);

	error = GRBaddconstr(Gurobidata->model, nz, indices, values, sense, rhs, nullptr);

	// free stuff
	if (indices)
		free(indices);
	libData->UTF8String_disown(senseString);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*	       SpArrayData: convenient data wrapper for MSparseArray        */
/************************************************************************/
typedef struct
{
	mint nentries;
	mint rank;
	mint const* dims;
	mint* explicitPositions;
	mreal* implicitValue;
	mreal* explicitValues;
} SpArrayData;

SpArrayData SpArrayData_fromMSparseArray(WolframLibraryData libData, MSparseArray S, MTensor* pT1)
{
	WolframSparseLibrary_Functions spFuns = libData->sparseLibraryFunctions;
	MTensor* pT = nullptr;
	SpArrayData sad;
	mint error;

	sad.rank = (*(spFuns->MSparseArray_getRank))(S);
	sad.dims = (*(spFuns->MSparseArray_getDimensions))(S);

	pT = (*(spFuns->MSparseArray_getImplicitValue))(S);
	sad.implicitValue = libData->MTensor_getRealData(*pT);

	pT = (*(spFuns->MSparseArray_getExplicitValues))(S);
	sad.explicitValues = libData->MTensor_getRealData(*pT);
	if (*pT != nullptr)
	{
		sad.nentries = libData->MTensor_getFlattenedLength(*pT);
	}
	else
	{
		sad.nentries = 0;
	}

	error = (*(spFuns->MSparseArray_getExplicitPositions))(S, pT1);
	sad.explicitPositions = libData->MTensor_getIntegerData(*pT1);

	return sad;
}

/************************************************************************/
/*                GurobiData_AddLinearConstraint1                       */
/*                                                                      */
/*         GurobiAddLinearConstraint1[data, vec, sense, rhs]            */
/*     sense: GRB_EQUAL, GRB_LESS_EQUAL, GRB_GREATER_EQUAL              */
/************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_AddLinearConstraint1(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, nz, *indices = nullptr;
	char sense, *senseString = nullptr;
	double *values = nullptr, rhs;
	mint i, dataID;
	MTensor Tvec = nullptr;
	GurobiData Gurobidata;
	SpArrayData vector;
	MSparseArray vectorSA;

	if (Argc != 4)
	{
		return LIBRARY_FUNCTION_ERROR;
	}

	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	/* find nz, indices and values */
	vectorSA = MArgument_getMSparseArray(Args[1]);
	vector = SpArrayData_fromMSparseArray(libData, vectorSA, &Tvec);
	nz = (int)vector.nentries;
	indices = (int*)malloc(nz * sizeof(int));
	for (i = 0; i < nz; i++)
	{
		indices[i] = (int)vector.explicitPositions[i] - 1;
	}
	values = vector.explicitValues;

	senseString = MArgument_getUTF8String(Args[2]);
	sense = senseString[0];

	rhs = MArgument_getReal(Args[3]);

	error = GRBaddconstr(Gurobidata->model, nz, indices, values, sense, rhs, nullptr);

	// free stuff
	if (Tvec)
		libData->MTensor_free(Tvec);
	if (indices)
		free(indices);
	libData->UTF8String_disown(senseString);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GurobiData_AddLinearConstraints                       */
/*                                                                      */
/*          GurobiAddLinearConstraints[data, mat, sense, rhs]           */
/*         sense: GRB_EQUAL, GRB_LESS_EQUAL, GRB_GREATER_EQUAL          */
/************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_AddLinearConstraints(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	char sense, *senseString = nullptr, *sensevec = nullptr;
	int error, ncons, Anz, *Ap = nullptr, *Ai = nullptr;
	mint j, dataID, *idata = nullptr;
	double *rhs = nullptr, *Ax = nullptr;
	MTensor pT, *Tp = nullptr;
	GurobiData Gurobidata;
	MSparseArray AMat = nullptr;

	if (Argc != 4)
		return LIBRARY_FUNCTION_ERROR;

	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	pT = MArgument_getMTensor(Args[3]);
	rhs = libData->MTensor_getRealData(pT);
	ncons = (int)(libData->MTensor_getFlattenedLength(pT));

	AMat = MArgument_getMSparseArray(Args[1]);

	Tp = (*libData->sparseLibraryFunctions->MSparseArray_getRowPointers)(AMat);
	if (libData->MTensor_getFlattenedLength(*Tp) != ncons + 1)
		return LIBRARY_FUNCTION_ERROR;
	idata = libData->MTensor_getIntegerData(*Tp);
	Ap = (int*)malloc(ncons * sizeof(int));

	for (j = 0; j < ncons; j++)
		Ap[j] = (int)idata[j];
	Anz = (int)idata[ncons];
	if (Anz == 0)
	{
		free(Ap);
		Ap = nullptr;
		Ai = nullptr;
		Ax = nullptr;
	}
	else
	{
		Ai = (int*)malloc(sizeof(int) * Anz);
		Tp = (*libData->sparseLibraryFunctions->MSparseArray_getColumnIndices)(AMat);
		if (libData->MTensor_getFlattenedLength(*Tp) != Anz)
			return LIBRARY_FUNCTION_ERROR;
		idata = libData->MTensor_getIntegerData(*Tp);
		for (j = 0; j < Anz; j++)
			Ai[j] = (int)(idata[j] - 1);

		Tp = (*libData->sparseLibraryFunctions->MSparseArray_getExplicitValues)(AMat);
		if (libData->MTensor_getFlattenedLength(*Tp) != Anz)
			return LIBRARY_FUNCTION_ERROR;
		Ax = libData->MTensor_getRealData(*Tp);
	}

	senseString = MArgument_getUTF8String(Args[2]);
	sense = senseString[0];

	sensevec = (char*)malloc(sizeof(char) * ncons);
	for (j = 0; j < ncons; j++)
		sensevec[j] = (char)sense;

	pT = MArgument_getMTensor(Args[3]);
	rhs = libData->MTensor_getRealData(pT);

	error = GRBaddconstrs(Gurobidata->model, ncons, Anz, Ap, Ai, Ax, sensevec, rhs, nullptr);

	// free stuff
	if (Ap)
		free(Ap);
	if (Ai)
		free(Ai);
	if (sensevec)
		free(sensevec);
	libData->UTF8String_disown(senseString);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/****************************************************************************************************/
/*                GurobiData_AddQuadraticConstraint                                                 */
/*                                                                                                  */
/*    GurobiAddQuadraticConstraint[data, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]   */
/*     sense: GRB_EQUAL, GRB_LESS_EQUAL, GRB_GREATER_EQUAL                                          */
/****************************************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_AddQuadraticConstraint(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, numlinnz, numquadnz, *linind = nullptr, *quadrow = nullptr, *quadcol = nullptr;
	char sense, *senseString = nullptr;
	double *linval, *quadval, rhs;
	mint j, dataID, *idata;
	MTensor pT;
	GurobiData Gurobidata;

	if (Argc != 8)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	pT = MArgument_getMTensor(Args[1]);
	numlinnz = (int)(libData->MTensor_getFlattenedLength(pT));

	if (numlinnz == 0)
	{
		linind = nullptr;
		linval = nullptr;
	}
	else
	{
		idata = libData->MTensor_getIntegerData(pT);
		linind = (int*)malloc(numlinnz * sizeof(int));
		for (j = 0; j < numlinnz; j++)
		{
			linind[j] = (char)idata[j] - 1;
		}
		pT = MArgument_getMTensor(Args[2]);
		linval = libData->MTensor_getRealData(pT);
	}

	pT = MArgument_getMTensor(Args[3]);
	numquadnz = (int)(libData->MTensor_getFlattenedLength(pT));
	idata = libData->MTensor_getIntegerData(pT);
	quadrow = (int*)malloc(numquadnz * sizeof(int));
	quadcol = (int*)malloc(numquadnz * sizeof(int));
	for (j = 0; j < numquadnz; j++)
	{
		quadrow[j] = (int)idata[j] - 1;
	}
	pT = MArgument_getMTensor(Args[4]);
	idata = libData->MTensor_getIntegerData(pT);
	for (j = 0; j < numquadnz; j++)
	{
		quadcol[j] = (int)idata[j] - 1;
	}
	pT = MArgument_getMTensor(Args[5]);
	quadval = libData->MTensor_getRealData(pT);

	senseString = MArgument_getUTF8String(Args[6]);
	sense = senseString[0];

	rhs = MArgument_getReal(Args[7]);

	error = GRBaddqconstr(Gurobidata->model, numlinnz, linind, linval, numquadnz, quadrow, quadcol, quadval, sense, rhs, "qc");

	// free stuff
	if (linind)
		free(linind);
	if (quadrow)
		free(quadrow);
	if (quadcol)
		free(quadcol);
	libData->UTF8String_disown(senseString);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/**********************************************************************/
/*                GurobiData_AddQuadraticConstraint1                  */
/*                                                                    */
/*    GurobiAddQuadraticConstraint1[data, mat, vec, sense, rhs]       */
/*                   1/2 x^T.Q.x + q.x sense b                        */
/**********************************************************************/

EXTERN_C DLLEXPORT int GurobiData_AddQuadraticConstraint1(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, numlinnz, numquadnz, *linind, *quadrow, *quadcol;
	char sense, *senseString = nullptr;
	double *linval, *quadval, rhs;
	mint i, dataID;
	MTensor TQ = nullptr, Tq = nullptr;
	GurobiData Gurobidata;
	SpArrayData Qmat, qvec;
	MSparseArray QmatSA, qvecSA;

	if (Argc != 5)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	QmatSA = MArgument_getMSparseArray(Args[1]);
	Qmat = SpArrayData_fromMSparseArray(libData, QmatSA, &TQ);
	numquadnz = (int)Qmat.nentries;
	quadrow = (int*)malloc(numquadnz * sizeof(int));
	quadcol = (int*)malloc(numquadnz * sizeof(int));
	for (i = 0; i < numquadnz; i++)
	{
		quadrow[i] = (int)Qmat.explicitPositions[2 * i] - 1;
		quadcol[i] = (int)Qmat.explicitPositions[2 * i + 1] - 1;
	}
	quadval = Qmat.explicitValues;
	qvecSA = MArgument_getMSparseArray(Args[2]);
	qvec = SpArrayData_fromMSparseArray(libData, qvecSA, &Tq);
	numlinnz = (int)qvec.nentries;
	if (numlinnz == 0)
	{
		linind = nullptr;
		linval = nullptr;
	}
	else
	{
		linind = (int*)malloc(numlinnz * sizeof(int));
		for (i = 0; i < numlinnz; i++)
		{
			linind[i] = (int)qvec.explicitPositions[i] - 1;
		}
		linval = qvec.explicitValues;
	}

	senseString = MArgument_getUTF8String(Args[3]);
	sense = senseString[0];

	rhs = MArgument_getReal(Args[4]);

	error = GRBaddqconstr(Gurobidata->model, numlinnz, linind, linval, numquadnz, quadrow, quadcol, quadval, sense, rhs, "qc1");

	// free stuff
	if (Tq)
		libData->MTensor_free(Tq);
	if (TQ)
		libData->MTensor_free(TQ);
	if (linind)
		free(linind);
	if (quadrow)
		free(quadrow);
	if (quadcol)
		free(quadcol);
	libData->UTF8String_disown(senseString);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GurobiData_SetParameters                              */
/*                                                                      */
/*        GurobiSetParameters[data, maxit, tol, nonconvex]              */
/************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_SetParameters(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, maxit, nonconvex;
	mint dataID;
	double tol;
	GurobiData Gurobidata;

	if (Argc != 4)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	maxit = (int)MArgument_getInteger(Args[1]);
	error = GRBsetintparam(GRBgetenv(Gurobidata->model), "BarIterLimit", maxit);
	if (error)
		return LIBRARY_FUNCTION_ERROR;

	tol = MArgument_getReal(Args[2]);
	if (tol > 1)
	{
		tol = 1.;
	}
	error = GRBsetdblparam(GRBgetenv(Gurobidata->model), "BarConvTol", tol / 100);
	if (error)
		return LIBRARY_FUNCTION_ERROR;
	error = GRBsetdblparam(GRBgetenv(Gurobidata->model), "BarQCPConvTol", tol);
	if (error)
		return LIBRARY_FUNCTION_ERROR;

	nonconvex = (int)MArgument_getInteger(Args[3]);
	error = GRBsetintparam(GRBgetenv(Gurobidata->model), "NonConvex", nonconvex);
	if (error)
		return LIBRARY_FUNCTION_ERROR;

	/* load return values */
	MArgument_setInteger(Res, 0);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GurobiData_SetStartingPoint                           */
/*                                                                      */
/*             GurobiSetStartingPoint[data, initpt]                     */
/************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_SetStartingPoint(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	double* start;
	MTensor pT;
	GurobiData Gurobidata;

	if (Argc != 2)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	pT = MArgument_getMTensor(Args[1]);
	start = libData->MTensor_getRealData(pT);

	error = GRBsetdblattrarray(Gurobidata->model, "Start", 0, (int)(Gurobidata->nvars), start);

	// load return values
	MArgument_setInteger(Res, error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GurobiData_OptimizeModel                              */
/*                                                                      */
/*                  GurobiOptimize[data]                                */
/************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_OptimizeModel(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	GurobiData Gurobidata;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	/* optimize model */
	error = GRBoptimize(Gurobidata->model);

	/* load return values */
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/*****************************************************************/
/*                GurobiData_GetStatusValue                      */
/*                                                               */
/*                 GurobiStatusValue[data]                       */
/*****************************************************************/

EXTERN_C DLLEXPORT int GurobiData_GetStatusValue(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	GurobiData Gurobidata;
	int optimstatus;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	/* get solution status */
	error = GRBgetintattr(Gurobidata->model, GRB_INT_ATTR_STATUS, &optimstatus);

	if (error)
		return LIBRARY_FUNCTION_ERROR;

	/* load return values */
	MArgument_setInteger(Res, (mint)optimstatus);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GurobiData_GetObjectiveValue                          */
/*                                                                      */
/*                 GurobiObjectiveValue[data]                           */
/************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_GetObjectiveValue(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	GurobiData Gurobidata;
	double objval = -1;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	/* get optimal objective value */
	error = GRBgetdblattr(Gurobidata->model, GRB_DBL_ATTR_OBJVAL, &objval);

	if (error)
		return LIBRARY_FUNCTION_ERROR;

	/* load return values */
	MArgument_setReal(Res, objval);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                    GurobiData_Getx                                   */
/*                                                                      */
/*                      Gurobix[data]                                   */
/************************************************************************/

inline mint setTensor(WolframLibraryData libData, const mreal* t, mint tdim, MTensor& T)
{
	if (T)
	{
		// free the tensor if it is already set
		libData->MTensor_free(T);
		T = nullptr;
	}
	// allocate the new tensor
	mint dims[1];
	dims[0] = tdim;
	libData->MTensor_new(MType_Real, 1, dims, &T);
	// fill the tensor
	mreal* tdata = libData->MTensor_getRealData(T);
	for (int i = 0; i < tdim; i++)
		tdata[i] = t[i];
	return 0;
}

EXTERN_C DLLEXPORT int GurobiData_Getx(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	GurobiData Gurobidata;
	double* sol = nullptr;
	MTensor solT = nullptr;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	/* get the minimizer vector */
	sol = (double*)malloc(sizeof(double) * (Gurobidata->nvars));
	error = GRBgetdblattrarray(Gurobidata->model, GRB_DBL_ATTR_X, 0, (int)(Gurobidata->nvars), sol);

	if (error)
		return LIBRARY_FUNCTION_ERROR;

	setTensor(libData, sol, Gurobidata->nvars, solT);

	// load return values
	MArgument_setMTensor(Res, solT);

	// free stuff
	if (sol)
		free(sol);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                    GurobiData_GetSlack                               */
/*                                                                      */
/*                     GurobiSlack[data]                                */
/************************************************************************/

EXTERN_C DLLEXPORT int GurobiData_GetSlack(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, ncons, nqcons, allcons;
	mint dataID;
	GurobiData Gurobidata;
	double* slack = nullptr;
	MTensor slackT = nullptr;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	Gurobidata = GurobiDataMap_get(dataID);

	error = GRBgetintattr(Gurobidata->model, "NumConstrs", &ncons);
	if (error)
		return LIBRARY_FUNCTION_ERROR;
	error = GRBgetintattr(Gurobidata->model, "NumQConstrs", &nqcons);
	if (error)
		return LIBRARY_FUNCTION_ERROR;
	allcons = ncons + nqcons;

	/* get the minimizer vector */
	slack = (double*)malloc(sizeof(double) * (allcons));
	error = GRBgetdblattrarray(Gurobidata->model, "Slack", 0, allcons, slack);

	if (error)
		return LIBRARY_FUNCTION_ERROR;

	setTensor(libData, slack, allcons, slackT);

	/* load return values */
	MArgument_setMTensor(Res, slackT);

	// free stuff
	if (slack)
		free(slack);

	return LIBRARY_NO_ERROR;
}
