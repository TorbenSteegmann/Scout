#pragma once

#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <vector>

#include "Common.hpp"
#include "MatrixOperations.hpp"

namespace scout
{
class Relation
{
public:
    // algorithm 1 of the thesis
    void CalculateTransitiveClosure();

    // maxConsistent
    std::optional<int> MaxConsistent(int b, matrix LambdaB);
    // minGamma in the thesis
    static std::optional<int> ParametricConsistencyCheck(int alpha, int beta);

    // maxPeriodic
    std::optional<int> MaxPeriodic(matrix const& LambdaB, int c);
    // min Kappa in the thesis
    static std::optional<int> CheckPeriod(matrix const& LambdaB, matrix const& m2, int const l);

    // prints the calculatedTransitiveClosure; if you set MAX_SMT in common.h to true it will also make an SMTLIB file (requires you to manually change the path)
    void PrintTransitiveClosure();

    void PrintCell(cell c, int valI, int valJ, char signOne, char signTwo, std::string& fileString, bool& containsK);

    std::string SearchVariable(int value);

    void SetVariableMap(std::map<int, std::string> const& variableMap);

    std::map<int, std::string> GetVariableMap();

    void SetIsOctagonal(bool isOctagonal);

    [[nodiscard]] bool GetIsOctagonal() const;

    void AddPowerOfRelation(int power, matrix const& m);

    void CalcAddPowerOfRelation(int power);

    matrix CalcNextPowerOfRelation(matrix m);

    std::pair<int, matrix> SearchPowerOfRelation(int power);

    static bool ConsistencyCheck(matrix m);


private:
    std::map<int, std::string> variableMap;
    std::map<int, matrix> powersOfRelation;
    std::vector<matrix> transitiveClosure;
    int prefix = 0;
    bool isOctagonal;
    matrix test;
};
} // namespace scout
