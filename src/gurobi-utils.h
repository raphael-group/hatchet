#ifndef _UTILS_GUROBI_H_
#define _UTILS_GUROBI_H_

#include "utils.h"
#include <gurobi_c++.h>

/// Gurobi variable array
typedef std::vector<GRBVar> VarArray;
/// Gurobi variable matrix
typedef std::vector<VarArray> VarMatrix;
/// Gurobi variable 3D matrix
typedef std::vector<VarMatrix> Var3Matrix;
/// Gurobi variable 4D matrix
typedef std::vector<Var3Matrix> Var4Matrix;
/// Gurobi variable 5D matrix
typedef std::vector<Var4Matrix> Var5Matrix;
/// Gurobi variable 6D matrix
typedef std::vector<Var5Matrix> Var6Matrix;

#define NON_ZERO 0.01

#endif // _UTILS_GUROBI_H_
