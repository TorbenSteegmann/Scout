#pragma once

#include <set>
#include <vector>

#include "Common.hpp"
#include "Relation.hpp"

namespace scout
{
namespace MatrixOperations
{

// these functions refer to the algorithms shown in the thesis
matrix ParametricFloydWarshallAlgorithm(matrix m, bool tighten);

cell MinTerms(cell const& c1, cell const& c2, cell const& c3);

cell MinCell(std::set<std::pair<int, int>> const& omniSet);

matrix CalculateTightClosure(matrix m);

cell HalfTerms(cell const& c);

int HalfInt(int val);

int IDash(int val);

void PrintMatrix(matrix const& m);

matrix IntegerMatrixSubtraction(matrix m1, matrix const& m2);

matrix MatrixAddition(matrix m1, matrix const& m2);

matrix MatrixComposition(matrix const& m1, matrix const& m2, bool tighten);
matrix CalcExtremalPaths(matrix m);

} // namespace MatrixOperations
} // namespace scout
