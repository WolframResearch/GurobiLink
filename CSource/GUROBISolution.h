
#ifndef GUROBI_SOLUTION_H
#define GUROBI_SOLUTION_H

#include <cstdlib>
#include <cstdio>

#include "WolframLibrary.h"
#include "WolframSparseLibrary.h"
#include "gurobi_c.h"

typedef struct GUROBIData_struct
{
	// GRBenv *env;
	GRBmodel* model;
	mint nvars;
	int optstatus;
	int error;

} * GUROBIData;

typedef struct GUROBIEnvironment_struct
{
	GRBenv* env;
	int error;
} * GUROBIEnvironment;

// Adds a data instance to the solution map if mode is 0 and deletes a solution with given id if mode is 1
EXTERN_C DLLEXPORT void GUROBIDataMap_manage(WolframLibraryData libData, mbool mode, mint id);

// Deletes a data instance with given id;
EXTERN_C DLLEXPORT int GUROBIDataMap_delete(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res);

// Retrieves a data instance from the solution map with a specified solution id
EXTERN_C DLLEXPORT GUROBIData GUROBIDataMap_get(mint id);

// Adds an environment instance to the solution map if mode is 0 and deletes a solution with given id if mode is 1
EXTERN_C DLLEXPORT void GUROBIEnvironmentMap_manage(WolframLibraryData libData, mbool mode, mint id);

// Deletes an environment instance with given id;
EXTERN_C DLLEXPORT int GUROBIEnvironmentMap_delete(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res);

// Retrieves a data instance from the solution map with a specified solution id
EXTERN_C DLLEXPORT GUROBIEnvironment GUROBIEnvironmentMap_get(mint id);

#endif
