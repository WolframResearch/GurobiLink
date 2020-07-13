#include <cstdlib>
#include "WolframLibrary.h"
#include "GUROBISolution.h"
#include <unordered_map>
#include "gurobi_c.h"

int GUROBIData_initialize(GUROBIData GUROBIdata)
{
	GUROBIdata->model = nullptr;
	GUROBIdata->nvars = 0;
	GUROBIdata->optstatus = 0;
	GUROBIdata->error = 0;
	return 0;
}

int GUROBIEnvironment_initialize(GUROBIEnvironment GUROBIenvironment)
{
	GUROBIenvironment->env = nullptr;
	GUROBIenvironment->error = 0;
	return 0;
}

GUROBIData GUROBIData_new()
{
	GUROBIData GUROBIdata;
	GUROBIdata = (GUROBIData)malloc(sizeof(*GUROBIdata));
	GUROBIData_initialize(GUROBIdata);
	GUROBIEnvironment GUROBIenvironment = GUROBIEnvironmentMap_get(1);
	GUROBIdata->error =
		GRBnewmodel(GUROBIenvironment->env, &(GUROBIdata->model), "Mo", 0, nullptr, nullptr, nullptr, nullptr, nullptr);
	return GUROBIdata;
}

GUROBIEnvironment GUROBIEnvironment_new()
{
	GUROBIEnvironment GUROBIenvironment;
	GUROBIenvironment = (GUROBIEnvironment)malloc(sizeof(*GUROBIenvironment));
	GUROBIEnvironment_initialize(GUROBIenvironment);
	GUROBIenvironment->error = GRBemptyenv(&(GUROBIenvironment->env));
	if (!GUROBIenvironment->error)
	{
		/* 0 variables, no problem info yet */
		GUROBIenvironment->error = GRBsetstrparam(GUROBIenvironment->env, "LogFile", "Mo.log");
		GUROBIenvironment->error = GRBstartenv(GUROBIenvironment->env);
	}
	return GUROBIenvironment;
}

int GUROBIData_delete(GUROBIData gurobidata)
{
	GRBfreemodel(gurobidata->model);
	free(gurobidata);
	gurobidata = nullptr;
	return 0;
}

int GUROBIEnvironment_delete(GUROBIEnvironment GUROBIenvironment)
{
	GRBfreeenv(GUROBIenvironment->env);
	free(GUROBIenvironment);
	GUROBIenvironment = nullptr;
	return 0;
}

static std::unordered_map<mint, GUROBIData> GUROBIDataMap;
static std::unordered_map<mint, GUROBIEnvironment> GUROBIEnvironmentMap;

EXTERN_C DLLEXPORT void GUROBIDataMap_manage(WolframLibraryData libData, mbool mode, mint id)
{
	if (mode == 0)
	{
		GUROBIDataMap[id] = GUROBIData_new();
	}
	else if (GUROBIDataMap[id] != nullptr)
	{
		GUROBIData_delete(GUROBIDataMap[id]);
		GUROBIDataMap.erase(id);
	}
}

EXTERN_C DLLEXPORT void GUROBIEnvironmentMap_manage(WolframLibraryData libData, mbool mode, mint id)
{
	if (mode == 0)
	{
		GUROBIEnvironmentMap[id] = GUROBIEnvironment_new();
	}
	else if (GUROBIEnvironmentMap[id] != nullptr)
	{
		GUROBIEnvironment_delete(GUROBIEnvironmentMap[id]);
		GUROBIEnvironmentMap.erase(id);
	}
}

EXTERN_C DLLEXPORT int GUROBIDataMap_delete(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res)
{
	mint id;
	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	id = MArgument_getInteger(Args[0]);
	if (GUROBIDataMap[id] == nullptr)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	GUROBIData_delete(GUROBIDataMap[id]);
	return (*libData->releaseManagedLibraryExpression)("GUROBI_data_instance_manager", id);
}

EXTERN_C DLLEXPORT int GUROBIEnvironmentMap_delete(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res)
{
	mint id;
	if (Argc != 1)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	id = MArgument_getInteger(Args[0]);
	if (GUROBIEnvironmentMap[id] == nullptr)
	{
		return LIBRARY_FUNCTION_ERROR;
	}
	GUROBIEnvironment_delete(GUROBIEnvironmentMap[id]);
	return (*libData->releaseManagedLibraryExpression)("GUROBI_environment_instance_manager", id);
}

GUROBIData GUROBIDataMap_get(mint id)
{
	return GUROBIDataMap[id];
}

GUROBIEnvironment GUROBIEnvironmentMap_get(mint id)
{
	return GUROBIEnvironmentMap[id];
}

EXTERN_C DLLEXPORT int GUROBIDataMap_retIDList(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res)
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
