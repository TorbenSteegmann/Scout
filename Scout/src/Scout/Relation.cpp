#include "Relation.hpp"

#include <fstream>
#include <utility>

void scout::Relation::CalculateTransitiveClosure()
{
    int b = 1;
    int b_jump = 1;
    while (true)
    {
        std::optional<int> K;
        for (int c = 1; c <= b; ++c)
        {
            for (int l = 0; l <= 2; ++l)
            {
                CalcAddPowerOfRelation(b + l * c);
                if (!ConsistencyCheck(powersOfRelation[b + l * c]))
                {
                    if (b == 1)
                    {
                        this->transitiveClosure.emplace_back(this->powersOfRelation[1]);
                        ++prefix;
                    }
                    for (int i = b + 1; i < b + l * c; ++i)
                    {
                        CalcAddPowerOfRelation(i);
                        this->transitiveClosure.emplace_back(this->powersOfRelation[i]);
                        ++this->prefix;
                    }
                    return;
                }
            }
            auto Lambda = MatrixOperations::IntegerMatrixSubtraction(powersOfRelation[b + c], powersOfRelation[b]);
            if (Lambda == MatrixOperations::IntegerMatrixSubtraction(powersOfRelation[b + 2 * c], powersOfRelation[b + c]))
            {
                auto LambdaB = MatrixOperations::MatrixAddition(powersOfRelation[b], Lambda);
                LambdaB = MatrixOperations::ParametricFloydWarshallAlgorithm(LambdaB, false);

                K = MaxConsistent(b, LambdaB);

                std::optional<int> L = MaxPeriodic(LambdaB, c);
                if (K)
                {
                    if (L)
                    {
                        L = std::min(*K, *L);
                        L = *L == 0 ? 2 : *L;
                    }
                    else
                    {
                        L = *K;
                    }
                }

                if (!L)
                {
                    this->transitiveClosure.emplace_back(LambdaB);
                    for (int j = 1; j < c; ++j)
                    {
                        CalcAddPowerOfRelation(j);
                        matrix LambdaBJ = MatrixOperations::MatrixComposition(LambdaB, this->powersOfRelation[j], true);
                        LambdaBJ = MatrixOperations::CalcExtremalPaths(LambdaBJ);
                        this->transitiveClosure.emplace_back(LambdaBJ);
                    }
                    return;
                }
                b_jump = std::max(b_jump, b + c * (*L + 1));
            }
        }
        int b_next = std::max(b + 1, b_jump);
        for (int i = b; i < b_next; ++i)
        {
            CalcAddPowerOfRelation(i);
            this->transitiveClosure.emplace_back(this->powersOfRelation[i]);
            ++this->prefix;
        }
        b = std::max(b + 1, b_jump);
    }
}


std::optional<int> scout::Relation::MaxConsistent(int b, matrix LambdaB)
{
    int const l = 2;
    std::optional<int> minGammaDB;
    std::optional<int> gamma;
    auto size = LambdaB.size();
    for (int i = 0; i < size; ++i)
    {
        for (auto const& term : LambdaB[i][i])
        {
            gamma = ParametricConsistencyCheck(term.first, term.second);
            if (!gamma)
            {
                continue;
            }
            if (*gamma == 3)
            {
                return l;
            }
            minGammaDB = minGammaDB ? std::min(*minGammaDB, *gamma) : *gamma;
        }
    }
    if (!this->isOctagonal)
    {
        minGammaDB = !minGammaDB ? minGammaDB : *minGammaDB - 1;
        return minGammaDB;
    }

    std::set<std::pair<int, int>> L;
    std::set<std::pair<int, int>> U;
    for (int i = 0; i < size; ++i)
    {
        auto cell1 = LambdaB[i][MatrixOperations::IDash(i)];
        auto cell2 = LambdaB[MatrixOperations::IDash(i)][i];
        if (cell1.empty() || cell2.empty())
        {
            continue;
        }
        for (auto const& term_i : cell1)
        {
            for (auto const& term_j : cell2)
            {
                L.emplace(term_i.first + term_j.first, MatrixOperations::HalfInt(term_i.second) + MatrixOperations::HalfInt(term_j.second));
                U.emplace(term_i.first + term_j.first,
                          MatrixOperations::HalfInt(term_i.first + term_i.second) + MatrixOperations::HalfInt(term_j.first + term_j.second));
            }
        }
    }
    std::optional<int> minGammaL;
    std::optional<int> gammaL;
    std::optional<int> minGammaU;
    std::optional<int> gammaU;
    for (auto const& term : L)
    {
        gammaL = ParametricConsistencyCheck(term.first, term.second);
        if (!gammaL)
        {
            continue;
        }
        if (*gammaL == l)
        {
            return gammaL;
        }
        minGammaL = minGammaL ? std::min(*minGammaL, *gammaL) : *gammaL;
    }
    for (auto const& term : U)
    {
        gammaU = ParametricConsistencyCheck(term.first, term.second);
        if (!gammaU)
        {
            continue;
        }
        if (*gammaU == l)
        {
            return gammaU;
        }
        minGammaU = minGammaU ? std::min(*minGammaU, *gammaU) : *gammaU;
    }

    std::vector<int> v;
    std::optional<int> minGamma;
    if (minGammaDB)
    {
        v.emplace_back(*minGammaDB);
    }
    if (minGammaL)
    {
        v.emplace_back(*minGammaL * 2);
    }
    if (minGammaU)
    {
        v.emplace_back(*minGammaU * 2 - 1);
    }
    for (int& i : v)
    {
        minGamma = !minGamma ? i : std::min(*minGamma, i);
    }
    minGamma = !minGamma ? minGamma : *minGamma - 1;
    return minGamma;
}

std::optional<int> scout::Relation::MaxPeriodic(matrix const& LambdaB, int c)
{
    auto const l = 0;
    CalcAddPowerOfRelation(c);
    auto LambdaBC = MatrixOperations::MatrixComposition(LambdaB, this->powersOfRelation[c], false);
    std::optional<int> kappa;

    if (!this->isOctagonal)
    {
        LambdaBC = MatrixOperations::CalcExtremalPaths(LambdaBC);
        auto ret = CheckPeriod(LambdaB, LambdaBC, l);
        if (ret)
        {
            *ret -= *ret;
        }
        return ret;
    }

    auto size = LambdaBC.size();
    auto size2 = LambdaB.size();
    matrix M_1_L = LambdaBC;
    matrix M_1_U = LambdaBC;
    matrix M_2_L = LambdaB;
    matrix M_2_U = LambdaB;

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            auto cell0 = LambdaBC[i][j];
            auto cell1 = LambdaBC[i][MatrixOperations::IDash(i)];
            auto cell2 = LambdaBC[MatrixOperations::IDash(j)][j];
            std::set<std::pair<int, int>> minTermsL;
            std::set<std::pair<int, int>> minTermsU;
            for (auto const& term : cell0)
            {
                // Split univariate cell_ij into L and U using Lemma 4.23
                minTermsL.emplace(2 * term.first, term.second);
                minTermsU.emplace(2 * term.first, term.first + term.second);
            }
            for (auto const& term_i : cell1)
            {
                for (auto const& term_j : cell2)
                {
                    // Add Halfterms in cell_ii' + cell_j'j using Lemma 4.23
                    minTermsL.emplace(term_i.first + term_j.first, MatrixOperations::HalfInt(term_i.second) + MatrixOperations::HalfInt(term_j.second));
                    minTermsU.emplace(term_i.first + term_j.first,
                                      MatrixOperations::HalfInt(term_i.first + term_i.second) + MatrixOperations::HalfInt(term_j.first + term_j.second));
                }
            }
            // ensure only min terms remain in both cells
            M_1_L[i][j] = MatrixOperations::MinCell(minTermsL);
            M_1_U[i][j] = MatrixOperations::MinCell(minTermsU);
        }
    }

    M_1_L = MatrixOperations::CalcExtremalPaths(M_1_L);
    M_1_U = MatrixOperations::CalcExtremalPaths(M_1_U);

    for (int i = 0; i < size2; ++i)
    {
        for (int j = 0; j < size2; ++j)
        {
            // M2 is already tightly closed, we just split it into L and U using Lemma 4.23
            if (LambdaB[i][j].empty())
            {
                continue;
            }
            M_2_L[i][j][0] = std::make_pair(2 * LambdaB[i][j][0].first, LambdaB[i][j][0].first + LambdaB[i][j][0].second); //+alpha da l+1
            M_2_U[i][j][0] = std::make_pair(2 * LambdaB[i][j][0].first, 2 * LambdaB[i][j][0].first + LambdaB[i][j][0].second);
        }
    }

    auto P_L = CheckPeriod(M_2_L, M_1_L, l);
    auto P_U = CheckPeriod(M_2_U, M_1_U, l);

    std::vector<int> v;
    std::optional<int> maxPeriod;
    if (P_L)
    {
        v.emplace_back(*P_L * 2 + 1);
    }
    if (P_U)
    {
        v.emplace_back(*P_U * 2);
    }
    for (int& i : v)
    {
        maxPeriod = !maxPeriod ? i : std::min(*maxPeriod, i);
    }

    return maxPeriod;
}

std::optional<int> scout::Relation::CheckPeriod(matrix const& m1, matrix const& m2, int const l)
{
    std::optional<int> kappa;
    auto size = m1.size();
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            if (m1[i][j].empty())
            {
                continue;
            }
            auto t_0 = m1[i][j][0];
            for (auto t_i : m2[i][j])
            {
                if (t_i.first == t_0.first)
                {
                    continue;
                }
                int candidate = (t_i.second - t_0.second) / (t_0.first - t_i.first);
                if (candidate <= l + 1) // needed in case of negative cylce
                    continue;
                kappa = !kappa ? candidate : std::min(candidate, *kappa);
            }
        }
    }
    return kappa;
}

std::pair<int, scout::matrix> scout::Relation::SearchPowerOfRelation(int power)
{
    if (power <= 0)
    {
        throw std::invalid_argument("Searched for negative power of relation.");
    }

    auto lowerBound = this->powersOfRelation.lower_bound(power);

    auto closestPower = lowerBound->first == power ? lowerBound->second : (--lowerBound)->second;

    return std::make_pair(lowerBound->first, closestPower);
}

void scout::Relation::SetVariableMap(std::map<int, std::string> const& variableMap) { this->variableMap = variableMap; }

std::map<int, std::string> scout::Relation::GetVariableMap() { return this->variableMap; }

void scout::Relation::SetIsOctagonal(bool isOctagonal) { this->isOctagonal = isOctagonal; }

bool scout::Relation::GetIsOctagonal() const { return this->isOctagonal; }

void scout::Relation::AddPowerOfRelation(int power, matrix const& m) { this->powersOfRelation.emplace(power, m); }

bool scout::Relation::ConsistencyCheck(matrix m)
{
    for (int i = 0; i < m.size(); ++i)
    {
        if (m[i][i][0].second < 0)
        { // only works because center diagonal is always < INF
            return false;
        }
    }
    return true;
}

std::optional<int> scout::Relation::ParametricConsistencyCheck(int alpha, int beta)
{
    std::optional<int> gamma;
    if (alpha < 0)
    {
        gamma = std::max(2, ((beta / (-1 * alpha)) + 1));
    }
    else if (alpha * 0 + beta < 0)
    {
        return 0;
    }
    return gamma;
}

void scout::Relation::CalcAddPowerOfRelation(int power)
{
    auto closestPower = SearchPowerOfRelation(power);
    while (closestPower.first < power)
    {
        ++closestPower.first;
        closestPower.second = CalcNextPowerOfRelation(closestPower.second);
        AddPowerOfRelation(closestPower.first, closestPower.second);
    }
}

scout::matrix scout::Relation::CalcNextPowerOfRelation(matrix m)
{
    m = MatrixOperations::MatrixComposition(m, this->powersOfRelation[1], this->isOctagonal);

    m = MatrixOperations::CalcExtremalPaths(m);

    return m;
}

void scout::Relation::PrintTransitiveClosure()
{
    auto numberOfRelations = transitiveClosure.size();
    std::string fileString;
    auto containsK = false;


    fileString.append(";Variable declarations\n");
    for (auto const& var : variableMap)
    {
        fileString.append("(declare-fun |" + var.second + "| () Int)\n");
        fileString.append("(declare-fun |" + var.second + "'| () Int)\n");
    }
    fileString.append("(declare-fun |$k| () Int)\n");
    fileString.append("(declare-fun k () Int)\n");
    fileString.append("(declare-fun t0 () bool)\n");
    fileString.append("(declare-fun t1 () bool)\n");
    for (auto i = 0; i < numberOfRelations; ++i)
    {
        fileString.append("(declare-fun s" + std::to_string(i) + " () bool)\n");
    }
    fileString.append("\n;Constraints\n");

    int k = 0;
    for (auto matrix : this->transitiveClosure)
    {
        ++k;
        std::cout << "(";
        fileString.append("(assert (= s" + std::to_string(k - 1) + " (and ");
        if (this->isOctagonal)
        {
            int relationHalfSize = (int)matrix.size() / 2;
            for (int i = 0; i < relationHalfSize; ++i)
            {
                int valI = i;
                int inputI = i;
                for (int j = i; j < relationHalfSize; ++j)
                {
                    int valJ = j;
                    int inputJ = j;

                    if (i == j)
                    {
                        PrintCell(matrix[2 * valI][2 * valI + 1], inputI, inputJ, '\0', '+', fileString, containsK);
                        PrintCell(matrix[2 * valI + 1][2 * valI], inputI, inputJ, '-', '-', fileString, containsK);
                        continue;
                    }

                    if (matrix[2 * valI][2 * valJ] == matrix[2 * valJ + 1][2 * valI + 1])
                    { // def 2.18
                        PrintCell(matrix[2 * valI][2 * valJ], inputI, inputJ, '\0', '-', fileString, containsK);
                    }
                    if (matrix[2 * valJ][2 * valI] == matrix[2 * valI + 1][2 * valJ + 1])
                    {
                        PrintCell(matrix[2 * valJ][2 * valI], inputI, inputJ, '-', '+', fileString, containsK);
                    }
                    if (matrix[2 * valI + 1][2 * valJ] == matrix[2 * valJ + 1][2 * valI])
                    {
                        PrintCell(matrix[2 * valI + 1][2 * valJ], inputI, inputJ, '-', '-', fileString, containsK);
                    }
                    if (matrix[2 * valI][2 * valJ + 1] == matrix[2 * valJ][2 * valI + 1])
                    {
                        PrintCell(matrix[2 * valI][2 * valJ + 1], inputI, inputJ, '\0', '+', fileString, containsK);
                    }
                }
            }
        }
        else
        {
            for (int i = 0; i < matrix.size(); ++i)
            {
                for (int j = 0; j < matrix.size(); ++j)
                {
                    int valI = i;
                    int valJ = j;
                    if (i == j)
                    {
                        continue;
                    }
                    PrintCell(matrix[i][j], valI, valJ, '\0', '-', fileString, containsK);
                }
            }
        }
        if (containsK)
        {
            fileString.append(" (>= |$k| 0))))\n");
            std::cout << " k >= 0)";
            containsK = false;
        }
        else
        {
            fileString.append(" )))\n");
            std::cout << " )";
        }

        if (k < numberOfRelations)
        {
            std::cout << " || \n";
        }
    }
    fileString.append("\n(assert (= t0 (or");
    for (int i = 0; i < numberOfRelations; ++i)
    {
        fileString.append(" s" + std::to_string(i));
    }
    fileString.append(")))");

    if (MAKE_SMT)
    {
        std::ofstream file("/path/to/file");
        file << fileString;
        file.close();
    }
}

void scout::Relation::PrintCell(cell c, int valI, int valJ, char signOne, char signTwo, std::string& fileString, bool& containsK)
{
    if (c.empty())
    {
        return;
    }
    fileString.append("(<= (");
    fileString += signTwo;
    fileString += ' ';

    if (signOne == '-')
    {
        std::cout << signOne << SearchVariable(valI);
        fileString += '(';
        fileString += signOne;
        fileString.append(" |" + SearchVariable(valI) + "|)");
    }
    else
    {
        std::cout << SearchVariable(valI);
        fileString.append("|" + SearchVariable(valI) + "|");
    }
    std::cout << signTwo << SearchVariable(valJ) << "<=";
    fileString.append(" |" + SearchVariable(valJ) + "|) ");

    auto alpha = c[0].first;
    auto beta = c[0].second;

    if (beta != 0 && alpha != 0)
    {
        containsK = true;
        std::cout << alpha << "k";
        if (beta > 0)
        {
            std::cout << "+";
        }
        std::cout << beta;
        fileString.append("(+ (* " + std::to_string(alpha) + " |$k|) " + std::to_string(beta) + ")");
    }
    else if (alpha != 0)
    {
        containsK = true;
        std::cout << alpha << "k";
        fileString.append("(* " + std::to_string(alpha) + " |$k|)");
    }
    else
    {
        std::cout << beta;
        fileString.append(std::to_string(beta));
    }
    std::cout << ",";
    fileString.append(") ");
}

std::string scout::Relation::SearchVariable(int value)
{
    ++value;
    for (auto const& pair : this->variableMap)
    {
        if (pair.first == value)
        {
            return pair.second;
        }
        else if (pair.first == value - this->variableMap.size())
        {
            return pair.second + "'";
        }
    }
    return "null";
}
