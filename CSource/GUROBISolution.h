
#ifndef GUROBI_SOLUTION_H
#define GUROBI_SOLUTION_H

#include <stdlib.h>
#include <stdio.h>

#include "WolframLibrary.h"
#include "WolframSparseLibrary.h"
#include "gurobi_c.h"

typedef struct GUROBIData_struct {
	GRBenv *env;
	GRBmodel *model;
	mint nvars;
	mint ncons;
	int optstatus;
	int error;

} *GUROBIData;


extern GUROBIData GUROBIData_set(WolframLibraryData libData, mint solID);

// Adds a solution instance to the solution map if mode is 0 and deletes a solution with given id if mode is 1
EXTERN_C DLLEXPORT void GUROBIDataMap_manage(WolframLibraryData libData, mbool mode, mint id);

// Deletes a solution instance with given id;
EXTERN_C DLLEXPORT int GUROBIDataMap_delete(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument res);

// Retrieves a solution instance from the solution map with a specified solution id
EXTERN_C DLLEXPORT GUROBIData GUROBIDataMap_get(mint id);

#endif
