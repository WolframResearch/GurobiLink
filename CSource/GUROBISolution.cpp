#include <stdlib.h>
#include <stdio.h>
#include "WolframLibrary.h"
#include "WolframSparseLibrary.h"
#include "GUROBISolution.h"
#include <unordered_map>
#include "gurobi_c.h"

int GUROBIData_initialize(GUROBIData GUROBIdata) {
	GUROBIdata->env  = NULL;
	GUROBIdata->model = NULL;
	GUROBIdata->nvars = 0;
	GUROBIdata->ncons = 0;
	GUROBIdata->error = 0;
	return 0;
}

GUROBIData GUROBIData_new(void)
{
	GUROBIData GUROBIdata;
	GUROBIdata = (GUROBIData)malloc(sizeof(*GUROBIdata));
	GUROBIData_initialize(GUROBIdata);
	GUROBIdata->error = GRBemptyenv(&(GUROBIdata->env));
	if (!(GUROBIdata->error)) {
		/* 0 variables, no problem info yet */
		// GUROBIdata->error = GRBsetstrparam(GUROBIdata->env, "LogFile", "Mo.log");
		GUROBIdata->error = GRBstartenv(GUROBIdata->env);
		if (!(GUROBIdata->error)) {
			GUROBIdata->error = GRBnewmodel(GUROBIdata->env, &(GUROBIdata->model), "Mo", 0, NULL, NULL, NULL, NULL, NULL);
		}
	}
	return GUROBIdata;
}

int GUROBIData_delete(GUROBIData gurobidata) {
	/* delete env, model */
	GRBfreemodel(gurobidata->model);
	GRBfreeenv(gurobidata->env);
	free(gurobidata);
	gurobidata = NULL;
	return 0;
}

static std::unordered_map<mint, GUROBIData> GUROBIDataMap;

EXTERN_C DLLEXPORT void GUROBIDataMap_manage(WolframLibraryData libData, mbool mode, mint id)
{
	if (mode == 0) {
		GUROBIDataMap[id] = GUROBIData_new();
	}
	else if (GUROBIDataMap[id] != NULL) {
		GUROBIData_delete(GUROBIDataMap[id]);
		GUROBIDataMap.erase(id);
	}
}

EXTERN_C DLLEXPORT int GUROBIDataMap_delete(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument res)
{
	mint id;
	if (Argc != 1) {
		return LIBRARY_FUNCTION_ERROR;
	}
	id = MArgument_getInteger(Args[0]);
	if (GUROBIDataMap[id] == NULL) {
		return LIBRARY_FUNCTION_ERROR;
	}
	GUROBIData_delete(GUROBIDataMap[id]);
	return (*libData->releaseManagedLibraryExpression)("GUROBI_solution_instance_manager", id);
}

GUROBIData GUROBIDataMap_get(mint id)
{
	return GUROBIDataMap[id];
}

EXTERN_C DLLEXPORT int GUROBIDataMap_retIDList(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument res)
{
	mint i, num = GUROBIDataMap.size();
	mint dims[1];
	MTensor resTen;

	dims[0] = num;
	int err = libData->MTensor_new(MType_Integer, 1, dims, &resTen);
	if (err)
		return err;
	mint* elems = libData->MTensor_getIntegerData(resTen);
	std::unordered_map<mint, GUROBIData>::const_iterator iter = GUROBIDataMap.begin();
	std::unordered_map<mint, GUROBIData>::const_iterator end = GUROBIDataMap.end();
	for (i = 0; i < num; i++) {
		elems[i] = iter->first;
		if (iter != end) {
			iter++;
		}
	}
	MArgument_setMTensor(res, resTen);
	return err;
}

GUROBIData GUROBIData_set(WolframLibraryData libData, mint solID)
{
	GUROBIData GUROBIdata = GUROBIDataMap[solID];

	return GUROBIdata;
}
