
#ifndef GUROBI_SOLUTION_H
#define GUROBI_SOLUTION_H

#include <cstdlib>
#include <cstdio>

#include "WolframLibrary.h"
#include "WolframSparseLibrary.h"
#include "gurobi_c.h"

typedef struct GurobiData_struct
{
	// GRBenv *env;
	GRBmodel* model;
	mint nvars;
	int optstatus;
	int error;

} * GurobiData;

typedef struct GurobiEnvironment_struct
{
	GRBenv* env;
	int error;
} * GurobiEnvironment;

// Adds a data instance to the solution map if mode is 0 and deletes a solution with given id if mode is 1
EXTERN_C DLLEXPORT void GurobiDataMap_manage(WolframLibraryData libData, mbool mode, mint id);

// Deletes a data instance with given id;
EXTERN_C DLLEXPORT int GurobiDataMap_delete(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res);

// Retrieves a data instance from the solution map with a specified solution id
EXTERN_C DLLEXPORT GurobiData GurobiDataMap_get(mint id);

// Adds an environment instance to the solution map if mode is 0 and deletes a solution with given id if mode is 1
EXTERN_C DLLEXPORT void GurobiEnvironmentMap_manage(WolframLibraryData libData, mbool mode, mint id);

// Deletes an environment instance with given id;
EXTERN_C DLLEXPORT int GurobiEnvironmentMap_delete(WolframLibraryData libData, mint Argc, MArgument* Args, MArgument res);

// Retrieves a data instance from the solution map with a specified solution id
EXTERN_C DLLEXPORT GurobiEnvironment GurobiEnvironmentMap_get(mint id);

#endif
