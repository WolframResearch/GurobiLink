#include <cstdlib>
#include <cstring>
#include <algorithm>
#include "WolframLibrary.h"
#include "WolframSparseLibrary.h"
#include "GUROBISolution.h"
#include "gurobi_c.h"

EXTERN_C DLLEXPORT mint WolframLibrary_getVersion()
{
	return WolframLibraryVersion;
}

EXTERN_C DLLEXPORT int WolframLibrary_initialize(WolframLibraryData libData)
{

	int err;
	err = (*libData->registerLibraryExpressionManager)("GUROBI_data_instance_manager", GUROBIDataMap_manage);
	if (err)
		return err;
	err = (*libData->registerLibraryExpressionManager)("GUROBI_environment_instance_manager", GUROBIEnvironmentMap_manage);

	return err;
}

EXTERN_C DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData libData)
{

	(*libData->unregisterLibraryCallbackManager)("GUROBI_data_instance_manager");
	(*libData->unregisterLibraryCallbackManager)("GUROBI_environment_instance_manager");
}

/************************************************************************/
/*                GUROBIData_CheckLicense                               */
/*                                                                      */
/*             GUROBICheckLicense[environment]                          */
/************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_CheckLicense(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	mint envID;
	GUROBIEnvironment GUROBIenvironment;

	if (Argc != 1) {
		return LIBRARY_FUNCTION_ERROR;
	}

	envID = MArgument_getInteger(Args[0]);
	GUROBIenvironment = GUROBIEnvironmentMap_get(envID);

	// load return values
	MArgument_setInteger(Res, GUROBIenvironment->error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GUROBIData_CheckModel                                 */
/*                                                                      */
/*                GUROBICheckModel[data]                                */
/************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_CheckModel(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	mint dataID;
	GUROBIData GUROBIdata;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}

	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	// return the error from GRBnewmodel in GUROBIEData_new
	MArgument_setInteger(Res, (mint)(GUROBIdata->error));

	return LIBRARY_NO_ERROR;
}

/******************************************************************************/
/*             GUROBIData_SetVariableTypesAndObjectiveVector                  */
/*                                                                            */
/*     GUROBISetVariableTypesAndObjectiveVector[data, intvars, objvector]     */
/******************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_SetVariableTypesAndObjectiveVector(WolframLibraryData libData, mint Argc,
																	 MArgument* Args, MArgument Res)
{
	int error;
	mint i, nvars, nintvars, dataID, *intvars = nullptr;
	char* vartypes;
	double *objvec, *lbounds = nullptr;
	MTensor pT;
	GUROBIData GUROBIdata;

	if (Argc != 3)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	pT = MArgument_getMTensor(Args[1]);
	nintvars = libData->MTensor_getFlattenedLength(pT);
	intvars = libData->MTensor_getIntegerData(pT);

	pT = MArgument_getMTensor(Args[2]);
	nvars = libData->MTensor_getFlattenedLength(pT);
	objvec = libData->MTensor_getRealData(pT);
	GUROBIdata->nvars = nvars;

	lbounds = (double*)malloc(nvars * sizeof(double));
	// set lb to be a vector of -GRB_INFINITY, if we leave it NULL the default is 0.0
	std::fill_n(lbounds, nvars, -GRB_INFINITY);
	/*for (i = 0; i < nvars; i++) {
		lbounds[i] = -GRB_INFINITY;
	}*/

	vartypes = (char*)malloc(nvars * sizeof(char));
	memset(vartypes, 'C', nvars);
	/*for (i = 0; i < nvars; i++) {
		vartypes[i] = 'C';
	}*/
	for (i = 0; i < nintvars; i++)
	{
		vartypes[intvars[i] - 1] = 'I';
	}

	// Add objective vector and variable types
	error = GRBaddvars(GUROBIdata->model, (int)nvars, 0, nullptr, nullptr, nullptr, objvec, lbounds, nullptr, vartypes, nullptr);

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
/*             GUROBIData_SetVariableTypesAndBoundsAndObjectiveVector                                */
/*                                                                                                   */
/*  GUROBISetVariableTypesAndBoundsAndObjectiveVector[data, vartypes, lbounds, ubounds, objvector]   */
/*****************************************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_SetVariableTypesAndBoundsAndObjectiveVector(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, nvars;
	char* vartypes = nullptr;
	double *objvec, *lbounds, *ubounds;
	mint dataID;
	MTensor pT;
	GUROBIData GUROBIdata;

	if (Argc != 5)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	vartypes = MArgument_getUTF8String(Args[1]);
	nvars = (int)strlen(vartypes);

	GUROBIdata->nvars = nvars;

	pT = MArgument_getMTensor(Args[2]);
	lbounds = libData->MTensor_getRealData(pT);

	pT = MArgument_getMTensor(Args[3]);
	ubounds = libData->MTensor_getRealData(pT);

	pT = MArgument_getMTensor(Args[4]);
	objvec = libData->MTensor_getRealData(pT);

	// Add objective vector and variable types
	error = GRBaddvars(GUROBIdata->model, nvars, 0, nullptr, nullptr, nullptr, objvec, lbounds, ubounds, vartypes, nullptr);

	// free stuff

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

/**********************************************************************/
/*                GUROBIData_AddQuadraticObjectiveMatrix              */
/*                                                                    */
/*             GUROBIAddQuadraticObjectiveMatrix[data, Qmat]          */
/* adds the term 1/2 x^T.Q.x to already specified lin objective  c.x  */
/**********************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_AddQuadraticObjectiveMatrix(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error = 0, numquadnz, *quadrow = nullptr, *quadcol = nullptr;
	double* quadval = nullptr;
	mint i, dataID;
	MTensor TQ = nullptr;
	GUROBIData GUROBIdata;
	SpArrayData Qmat;
	MSparseArray QmatSA;

	if (Argc != 2)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	QmatSA = MArgument_getMSparseArray(Args[1]);
	Qmat = SpArrayData_fromMSparseArray(libData, QmatSA, &TQ);
	numquadnz = (int)Qmat.nentries;

	if (numquadnz > 0)
	{
		quadrow = (int*)malloc(numquadnz * sizeof(int));
		quadcol = (int*)malloc(numquadnz * sizeof(int));
		quadval = (double*)malloc(numquadnz * sizeof(double));
		for (i = 0; i < numquadnz; i++)
		{
			quadrow[i] = (int)Qmat.explicitPositions[2 * i] - 1;
			quadcol[i] = (int)Qmat.explicitPositions[2 * i + 1] - 1;
			quadval[i] = Qmat.explicitValues[i];
		}

		error = GRBaddqpterms(GUROBIdata->model, numquadnz, quadrow, quadcol, quadval);

		if (quadrow)
			free(quadrow);
		if (quadcol)
			free(quadcol);
		if (quadval)
			free(quadval);
	}

	// free stuff
	if (TQ)
		libData->MTensor_free(TQ);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/*****************************************************************************/
/*                GUROBIData_AddLinearConstraint                             */
/*                                                                           */
/*    GUROBIAddLinearConstraint[data, indices, values, sense, rhs]           */
/* sense: GRB_EQUAL -> '=', GRB_LESS_EQUAL -> '<', GRB_GREATER_EQUAL -> '>'  */
/*****************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_AddLinearConstraint(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int j, error, nz, *indices = nullptr;
	double *values = nullptr, rhs;
	mint dataID, *idata = nullptr;
	MTensor pT;
	GUROBIData GUROBIdata;
	char sense, *senseString;

	if (Argc != 5)
	{
		return LIBRARY_FUNCTION_ERROR;
	}

	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);
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

	// Add constraint coefficients indices and values, equal/ineq sign and right hand side
	// example: x + 2 y + 3 z <= 4
	// nz = 3;
	// indices[0] = 0; indices[1] = 1; indices[2] = 2;
	// indices[0] = 1; indices[1] = 2; indices[2] = 3;
	// sense = GRB_LESS_EQUAL
	// rhs = 4

	error = GRBaddconstr(GUROBIdata->model, nz, indices, values, sense, rhs, nullptr);

	// free stuff
	if (indices)
		free(indices);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GUROBIData_AddLinearConstraint1                       */
/*                                                                      */
/*         GUROBIAddLinearConstraint1[data, vec, sense, rhs]            */
/*     sense: GRB_EQUAL, GRB_LESS_EQUAL, GRB_GREATER_EQUAL              */
/************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_AddLinearConstraint1(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, nz, *indices = nullptr;
	char sense, *senseString = nullptr;
	double *values = nullptr, rhs;
	mint i, dataID;
	MTensor Tvec = nullptr;
	GUROBIData GUROBIdata;
	SpArrayData vector;
	MSparseArray vectorSA;

	if (Argc != 4)
	{
		return LIBRARY_FUNCTION_ERROR;
	}

	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	/* find nz, indices and values */
	vectorSA = MArgument_getMSparseArray(Args[1]);
	vector = SpArrayData_fromMSparseArray(libData, vectorSA, &Tvec);
	nz = (int)vector.nentries;
	indices = (int*)malloc(nz * sizeof(int));
	values = (double*)malloc(nz * sizeof(double));
	for (i = 0; i < nz; i++)
	{
		indices[i] = (int)vector.explicitPositions[i] - 1;
		values[i] = vector.explicitValues[i];
	}

	senseString = MArgument_getUTF8String(Args[2]);
	sense = senseString[0];

	rhs = MArgument_getReal(Args[3]);

	// Add constraint coefficients indices and values, equal/ineq sign and right hand side
	// example: x + 2 y + 3 z <= 4
	// nz = 3;
	// indices[0] = 0; indices[1] = 1; indices[2] = 2;
	// indices[0] = 1; indices[1] = 2; indices[2] = 3;
	// sense = GRB_LESS_EQUAL
	// rhs = 4

	error = GRBaddconstr(GUROBIdata->model, nz, indices, values, sense, rhs, nullptr);

	// free stuff
	if (Tvec)
		libData->MTensor_free(Tvec);
	if (indices)
		free(indices);
	if (values)
		free(values);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GUROBIData_AddLinearConstraints                       */
/*                                                                      */
/*          GUROBIAddLinearConstraints[data, mat, sense, rhs]           */
/*         sense: GRB_EQUAL, GRB_LESS_EQUAL, GRB_GREATER_EQUAL          */
/************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_AddLinearConstraints(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, ncons, Anz;
	char sense, *senseString = nullptr, *sensevec = nullptr;
	double *rhs = nullptr, *Ax = nullptr;
	mint j, dataID, *idata = nullptr;
	MTensor pT;
	GUROBIData GUROBIdata;
	MSparseArray AMat = nullptr;
	MTensor* Tp = nullptr;
	int *Ap = nullptr, *Ai = nullptr;

	if (Argc != 4)
	{
		return LIBRARY_FUNCTION_ERROR;
	}

	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

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
		return LIBRARY_FUNCTION_ERROR;

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

	senseString = MArgument_getUTF8String(Args[2]);
	sense = senseString[0];

	sensevec = (char*)malloc(sizeof(char) * ncons);
	for (j = 0; j < ncons; j++)
		sensevec[j] = (char)sense;

	pT = MArgument_getMTensor(Args[3]);
	rhs = libData->MTensor_getRealData(pT);

	error = GRBaddconstrs(GUROBIdata->model, ncons, Anz, Ap, Ai, Ax, sensevec, rhs, nullptr);

	// free stuff
	if (Ap)
		free(Ap);
	if (Ai)
		free(Ai);
	if (sensevec)
		free(sensevec);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/****************************************************************************************************/
/*                GUROBIData_AddQuadraticConstraint                                                 */
/*                                                                                                  */
/*    GUROBIAddQuadraticConstraint[data, linind, linvals, quadrow, quadcol, quadvals, sense, rhs]   */
/*     sense: GRB_EQUAL, GRB_LESS_EQUAL, GRB_GREATER_EQUAL                                          */
/****************************************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_AddQuadraticConstraint(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, numlinnz, numquadnz, *linind = nullptr, *quadrow = nullptr, *quadcol = nullptr;
	char sense, *senseString = nullptr;
	double *linval, *quadval, rhs;
	mint j, dataID, *idata;
	MTensor pT;
	GUROBIData GUROBIdata;

	if (Argc != 8)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

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

	// Add constraint coefficients indices and values, equal/ineq sign and right hand side
	// example: 2 x0^2 + x0 x1 + x1^2 + 2 x1 + x2 <= 1
	// int    numlinnz = 2;
	// int    linind[] = {1, 2};
	// double linval[] = {2.0, 1.0};
	// int    numquadnz = 3;
	// int    quadrow[] = {0, 0, 1};
	// int    quadcol[] = {0, 1, 1};
	// double quadval[] = {2.0, 1.0, 1.0};

	error = GRBaddqconstr(GUROBIdata->model, numlinnz, linind, linval, numquadnz, quadrow, quadcol, quadval, sense, rhs, "qc");

	// free stuff
	if (linind)
		free (linind);
	if (quadrow)
		free (quadrow);
	if (quadcol)
		free (quadcol);

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/**********************************************************************/
/*                GUROBIData_AddQuadraticConstraint1                  */
/*                                                                    */
/*    GUROBIAddQuadraticConstraint1[data, mat, vec, sense, rhs]       */
/*                   1/2 x^T.Q.x + q.x sense b                        */
/**********************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_AddQuadraticConstraint1(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, numlinnz, numquadnz, *linind, *quadrow, *quadcol;
	char sense, *senseString = nullptr;
	double *linval, *quadval, rhs;
	mint i, dataID;
	MTensor TQ = nullptr, Tq = nullptr;
	GUROBIData GUROBIdata;
	SpArrayData Qmat, qvec;
	MSparseArray QmatSA, qvecSA;

	if (Argc != 5)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

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

	// Add constraint coefficients indices and values, equal/ineq sign and right hand side
	// example: 2 x0^2 + x0 x1 + x1^2 + 2 x1 + x2 <= 1
	// int    numlinnz = 2;
	// int    linind[] = {1, 2};
	// double linval[] = {2.0, 1.0};
	// int    numquadnz = 3;
	// int    quadrow[] = {0, 0, 1};
	// int    quadcol[] = {0, 1, 1};
	// double quadval[] = {2.0, 1.0, 1.0};

	error = GRBaddqconstr(GUROBIdata->model, numlinnz, linind, linval, numquadnz, quadrow, quadcol, quadval, sense, rhs, "qc1");

	// error = GRBaddqconstr(GUROBIdata->model, 0, NULL, NULL, 0, NULL, NULL, NULL, sense, rhs, NULL);

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

	// load return values
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GUROBIData_SetParameters                              */
/*                                                                      */
/*        GUROBISetParameters[data, maxit, tol, nonconvex]              */
/************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_SetParameters(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, maxit, nonconvex;
	mint dataID;
	double tol;
	GUROBIData GUROBIdata;

	if (Argc != 4)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	maxit = (int)MArgument_getInteger(Args[1]);
	error = GRBsetintparam(GRBgetenv(GUROBIdata->model), "BarIterLimit", maxit);
	if (error)
		return LIBRARY_FUNCTION_ERROR;

	tol = MArgument_getReal(Args[2]);
	if (tol > 1)
	{
		tol = 1.;
	}
	error = GRBsetdblparam(GRBgetenv(GUROBIdata->model), "BarConvTol", tol / 100);
	if (error)
		return LIBRARY_FUNCTION_ERROR;
	error = GRBsetdblparam(GRBgetenv(GUROBIdata->model), "BarQCPConvTol", tol);
	if (error)
		return LIBRARY_FUNCTION_ERROR;

	nonconvex = (int)MArgument_getInteger(Args[3]);
	error = GRBsetintparam(GRBgetenv(GUROBIdata->model), "NonConvex", nonconvex);
	if (error)
		return LIBRARY_FUNCTION_ERROR;

	/* load return values */
	MArgument_setInteger(Res, 0);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GUROBIData_SetStartingPoint                           */
/*                                                                      */
/*             GUROBISetStartingPoint[data, initpt]                     */
/************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_SetStartingPoint(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	double* start;
	MTensor pT;
	GUROBIData GUROBIdata;

	if (Argc != 2)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	pT = MArgument_getMTensor(Args[1]);
	start = libData->MTensor_getRealData(pT);

	error = GRBsetdblattrarray(GUROBIdata->model, "Start", 0, (int)(GUROBIdata->nvars), start);

	// free stuff

	// load return values
	MArgument_setInteger(Res, error);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GUROBIData_OptimizeModel                              */
/*                                                                      */
/*                  GUROBIOptimize[data]                                */
/************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_OptimizeModel(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	GUROBIData GUROBIdata;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	/* optimize model */
	error = GRBoptimize(GUROBIdata->model);

	/* load return values */
	MArgument_setInteger(Res, (mint)error);

	return LIBRARY_NO_ERROR;
}

/*****************************************************************/
/*                GUROBIData_GetStatusValue                      */
/*                                                               */
/*                 GUROBIStatusValue[data]                       */
/*****************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_GetStatusValue(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	GUROBIData GUROBIdata;
	int optimstatus;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	/* get solution status */
	error = GRBgetintattr(GUROBIdata->model, GRB_INT_ATTR_STATUS, &optimstatus);

	if (error)
		return LIBRARY_FUNCTION_ERROR;

	/* load return values */
	MArgument_setInteger(Res, (mint)optimstatus);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                GUROBIData_GetObjectiveValue                          */
/*                                                                      */
/*                 GUROBIObjectiveValue[data]                           */
/************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_GetObjectiveValue(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	GUROBIData GUROBIdata;
	double objval = -1;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	/* get optimal objective value */
	error = GRBgetdblattr(GUROBIdata->model, GRB_DBL_ATTR_OBJVAL, &objval);

	if (error)
		return LIBRARY_FUNCTION_ERROR;

	/* load return values */
	MArgument_setReal(Res, objval);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                    GUROBIData_Getx                                   */
/*                                                                      */
/*                      GUROBIx[data]                                   */
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

EXTERN_C DLLEXPORT int GUROBIData_Getx(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error;
	mint dataID;
	GUROBIData GUROBIdata;
	double* sol = nullptr;
	MTensor solT = nullptr;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	/* get the minimizer vector */
	sol = (double*)malloc(sizeof(double) * (GUROBIdata->nvars));
	error = GRBgetdblattrarray(GUROBIdata->model, GRB_DBL_ATTR_X, 0, (int)(GUROBIdata->nvars), sol);

	if (error)
		return LIBRARY_FUNCTION_ERROR;

	setTensor(libData, sol, GUROBIdata->nvars, solT);

	// load return values
	MArgument_setMTensor(Res, solT);

	// free stuff
	if (sol)
		free(sol);

	return LIBRARY_NO_ERROR;
}

/************************************************************************/
/*                    GUROBIData_GetSlack                               */
/*                                                                      */
/*                     GUROBISlack[data]                                */
/************************************************************************/

EXTERN_C DLLEXPORT int GUROBIData_GetSlack(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument Res)
{
	int error, ncons, nqcons, allcons;
	mint dataID;
	GUROBIData GUROBIdata;
	double* slack = nullptr;
	MTensor slackT = nullptr;

	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	dataID = MArgument_getInteger(Args[0]);
	GUROBIdata = GUROBIDataMap_get(dataID);

	error = GRBgetintattr(GUROBIdata->model, "NumConstrs", &ncons);
	if (error)
		return LIBRARY_FUNCTION_ERROR;
	error = GRBgetintattr(GUROBIdata->model, "NumQConstrs", &nqcons);
	if (error)
		return LIBRARY_FUNCTION_ERROR;
	allcons = ncons + nqcons;

	// It is not clear what length to expect for the slack and in any case it doesn't seem to work for quadratic constraints

	/* get the minimizer vector */
	slack = (double*)malloc(sizeof(double) * (allcons));
	error = GRBgetdblattrarray(GUROBIdata->model, "Slack", 0, allcons, slack);

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
