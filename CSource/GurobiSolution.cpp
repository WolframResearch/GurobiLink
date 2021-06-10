#include "GurobiSolution.h"
#include <unordered_map>

int GurobiData_initialize(GurobiData Gurobidata)
{
	Gurobidata->model = nullptr;
	Gurobidata->nvars = 0;
	Gurobidata->optstatus = 0;
	Gurobidata->error = 0;
	return 0;
}

int GurobiEnvironment_initialize(GurobiEnvironment Gurobienvironment)
{
	Gurobienvironment->env = nullptr;
	Gurobienvironment->error = 0;
	return 0;
}

GurobiData GurobiData_new()
{
	GurobiData Gurobidata;
	Gurobidata = (GurobiData)malloc(sizeof(*Gurobidata));
	GurobiData_initialize(Gurobidata);
	GurobiEnvironment Gurobienvironment = GurobiEnvironmentMap_get(1);
	/* Create new model with 0 variables, no problem info yet */
	Gurobidata->error =
		GRBnewmodel(Gurobienvironment->env, &(Gurobidata->model), "Mo", 0, nullptr, nullptr, nullptr, nullptr, nullptr);
	return Gurobidata;
}

GurobiEnvironment GurobiEnvironment_new()
{
	GurobiEnvironment Gurobienvironment;
	Gurobienvironment = (GurobiEnvironment)malloc(sizeof(*Gurobienvironment));
	GurobiEnvironment_initialize(Gurobienvironment);
	Gurobienvironment->error = GRBemptyenv(&(Gurobienvironment->env));
	if (!Gurobienvironment->error)
	{
		// Gurobienvironment->error = GRBsetstrparam(Gurobienvironment->env, "LogFile", "Mo.log");
		Gurobienvironment->error = GRBstartenv(Gurobienvironment->env);
	}
	return Gurobienvironment;
}

int GurobiData_delete(GurobiData gurobidata)
{
	GRBfreemodel(gurobidata->model);
	free(gurobidata);
	gurobidata = nullptr;
	return 0;
}

int GurobiEnvironment_delete(GurobiEnvironment Gurobienvironment)
{
	GRBfreeenv(Gurobienvironment->env);
	free(Gurobienvironment);
	Gurobienvironment = nullptr;
	return 0;
}

static std::unordered_map<mint, GurobiData> GurobiDataMap;
static std::unordered_map<mint, GurobiEnvironment> GurobiEnvironmentMap;

EXTERN_C DLLEXPORT void GurobiDataMap_manage(WolframLibraryData libData, mbool mode, mint id)
{
	if (mode == 0)
	{
		GurobiDataMap[id] = GurobiData_new();
	}
	else if (GurobiDataMap[id] != nullptr)
	{
		GurobiData_delete(GurobiDataMap[id]);
		GurobiDataMap.erase(id);
	}
}

EXTERN_C DLLEXPORT void GurobiEnvironmentMap_manage(WolframLibraryData libData, mbool mode, mint id)
{
	if (mode == 0)
	{
		GurobiEnvironmentMap[id] = GurobiEnvironment_new();
	}
	else if (GurobiEnvironmentMap[id] != nullptr)
	{
		GurobiEnvironment_delete(GurobiEnvironmentMap[id]);
		GurobiEnvironmentMap.erase(id);
	}
}

EXTERN_C DLLEXPORT int GurobiDataMap_delete(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res)
{
	mint id;
	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	id = MArgument_getInteger(Args[0]);
	if (GurobiDataMap[id] == nullptr)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	GurobiData_delete(GurobiDataMap[id]);
	return (*libData->releaseManagedLibraryExpression)("Gurobi_data_instance_manager", id);
}

EXTERN_C DLLEXPORT int GurobiEnvironmentMap_delete(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res)
{
	mint id;
	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	id = MArgument_getInteger(Args[0]);
	if (GurobiEnvironmentMap[id] == nullptr)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	GurobiEnvironment_delete(GurobiEnvironmentMap[id]);
	return (*libData->releaseManagedLibraryExpression)("Gurobi_environment_instance_manager", id);
}

GurobiData GurobiDataMap_get(mint id)
{
	return GurobiDataMap[id];
}

GurobiEnvironment GurobiEnvironmentMap_get(mint id)
{
	return GurobiEnvironmentMap[id];
}

EXTERN_C DLLEXPORT int GurobiDataMap_retIDList(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res)
{
	mint i, num = GurobiDataMap.size();
	mint dims[1];
	MTensor resTen;

	dims[0] = num;
	int err = libData->MTensor_new(MType_Integer, 1, dims, &resTen);
	if (err)
		return err;
	mint* elems = libData->MTensor_getIntegerData(resTen);
	std::unordered_map<mint, GurobiData>::const_iterator iter = GurobiDataMap.begin();
	std::unordered_map<mint, GurobiData>::const_iterator end = GurobiDataMap.end();
	for (i = 0; i < num; i++)
	{
		elems[i] = iter->first;
		if (iter != end)
		{
			iter++;
		}
	}
	MArgument_setMTensor(res, resTen);
	return err;
}
