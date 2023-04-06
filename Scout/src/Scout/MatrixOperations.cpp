#include "MatrixOperations.hpp"

#include <chrono>
#include <iostream>

static bool isTerm1DominatedBy2(int a1, int b1, int a2, int b2)
{
    if (b1 < b2)
        return false; // for x = 0, term 1 is smaller than 2

    if (a1 < a2)
        return false; // for some large enough x, term 2 becomes larger than term 1

    return true;
}

static void updateCellWithTerm(scout::cell& c, int a_new, int b_new)
{
    // 1. check if (a_new, b_new) is minterm
    for (auto& [a, b] : c)
    {
        if (isTerm1DominatedBy2(a_new, b_new, a, b))
        {
            return;
        }
    }
    // 2. remove all previous terms that are not min now
    auto is_non_min = [a_new, b_new](std::pair<int, int> v) -> bool
    {
        auto [a, b] = v;
        return isTerm1DominatedBy2(a, b, a_new, b_new);
    };
    c.erase(std::remove_if(c.begin(), c.end(), is_non_min), c.end());

    // 3. actuall add term
    c.emplace_back(a_new, b_new);
}
scout::matrix scout::MatrixOperations::ParametricFloydWarshallAlgorithm(matrix m, bool tighten)
{
    auto size = int(m.size());

    cell tmp;

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < size; ++k)
    {
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                if (i == k || j == k)
                    continue;
                auto& c1 = m[i][j];
                auto const& c2 = m[i][k];
                auto const& c3 = m[k][j];

                tmp.clear();
                for (auto [a2, b2] : c2)
                {
                    for (auto [a3, b3] : c3)
                    {
                        tmp.emplace_back(a2 + a3, b2 + b3);
                    }
                }
                if (tmp.empty())
                {
                    continue;
                }

                std::sort(tmp.begin(), tmp.end());

                auto prev_a = tmp[0].first - 1;
                for (auto [a, b] : tmp)
                {
                    if (a != prev_a)
                    {
                        prev_a = a;

                        updateCellWithTerm(c1, a, b);
                    }
                }

                // verification if every term is minimal for some n
                cell tmp2;
                auto numberOfLoops = c1.size();
                for (int l = 0; l < numberOfLoops; ++l)
                {
                    auto term1 = c1[l];
                    bool minTerm = true;

                    for (int n = l + 1; n < numberOfLoops; ++n)
                    {
                        if (term1.first >= 0)
                            break;
                        auto term2 = c1[n];
                        if (term1.first == term2.first || term2.first >= 0)
                            continue;

                        auto minAlpha = term1.first < term2.first ? term1 : term2;
                        auto midAlpha = term1;
                        auto maxAlpha = term1.first > term2.first ? term1 : term2;

                        for (int o = n + 1; o < numberOfLoops; ++o)
                        {
                            auto term3 = c1[o];
                            if (term1.first == term3.first || term2.first == term3.first || term3.first >= 0)
                            {
                                continue;
                            }
                            if (term3.first < minAlpha.first)
                            {
                                midAlpha = minAlpha;
                                minAlpha = term3;
                            }
                            else if (term3.first > maxAlpha.first)
                            {
                                midAlpha = maxAlpha;
                                maxAlpha = term3;
                            }
                            else
                            {
                                midAlpha = term3;
                            }

                            auto s1 = (maxAlpha.second - minAlpha.second) / (minAlpha.first - maxAlpha.first);
                            auto s2 = (maxAlpha.second - midAlpha.second) / (midAlpha.first - maxAlpha.first);

                            if (s1 <= s2)
                            {
                                if (term1 == midAlpha)
                                {
                                    minTerm = false;
                                    break;
                                }
                                else if (term2 == midAlpha)
                                {
                                    c1[n] = term1;
                                    c1[l] = term2;
                                    minTerm = false;
                                    break;
                                }
                                else
                                {
                                    c1[o] = term1;
                                    c1[l] = term3;
                                    minTerm = false;
                                    break;
                                }
                            }
                        }
                        if (!minTerm)
                        {
                            break;
                        }
                    }
                    if (minTerm)
                    {
                        tmp2.emplace_back(term1);
                    }
                }
                m[i][j] = tmp2;
            }
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    // tighten
    if (tighten)
    {
        m = MatrixOperations::CalculateTightClosure(m);
    }

    return m;
}

scout::matrix scout::MatrixOperations::CalculateTightClosure(matrix m)
{
    auto size = m.size();
    matrix tightClosure = m;

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            auto cell1 = m[i][j];
            auto cell2 = MatrixOperations::HalfTerms(m[i][MatrixOperations::IDash(i)]);
            auto cell3 = MatrixOperations::HalfTerms(m[MatrixOperations::IDash(j)][j]);

            tightClosure[i][j] = MatrixOperations::MinTerms(cell1, cell2, cell3);
        }
    }

    return tightClosure;
}

scout::cell scout::MatrixOperations::MinTerms(cell const& c1, cell const& c2, cell const& c3)
{
    std::set<std::pair<int, int>> omniSet(c1.begin(), c1.end());
    std::set<std::pair<int, int>> minTerms;
    for (auto termInCell2 : c2)
    {
        for (auto termInCell3 : c3)
        {
            omniSet.emplace(termInCell2.first + termInCell3.first, termInCell2.second + termInCell3.second);
        }
    }

    return MinCell(omniSet);
}

scout::cell scout::MatrixOperations::MinCell(std::set<std::pair<int, int>> const& omniSet)
{
    cell minTerms;
    for (auto it1 = omniSet.begin(); it1 != omniSet.end(); ++it1)
    {
        bool isMinTerm = true;
        for (auto it2 = omniSet.begin(); it2 != omniSet.end(); ++it2)
        {
            if (it1 == it2)
            {
                continue;
            }
            if ((it2->first < it1->first || it2->second < it1->second) && (it2->first <= it1->first && it2->second <= it1->second))
            {
                isMinTerm = false;
                break;
            }
        }
        if (isMinTerm)
        {
            minTerms.emplace_back(*it1);
        }
    }
    return minTerms;
}

scout::cell scout::MatrixOperations::HalfTerms(cell const& c)
{
    cell halfTerms;

    for (auto term : c)
    {
        halfTerms.emplace_back(MatrixOperations::HalfInt(term.first), MatrixOperations::HalfInt(term.second));
    }

    return halfTerms;
}

int scout::MatrixOperations::HalfInt(int val)
{
    auto value = ((val < 0) && (val % 2 != 0)) ? (val / 2) - 1 : val / 2;
    return value;
}

int scout::MatrixOperations::IDash(int val)
{
    auto i = (val % 2 == 0) ? val + 1 : val - 1;
    return i;
}

[[maybe_unused]] void scout::MatrixOperations::PrintMatrix(matrix const& m)
{
    for (auto const& row : m)
    {
        for (auto const& column : row)
        {
            if (column.empty())
            {
                std::cout << "INFTY\t";
                continue;
            }
            std::cout << column[0].first << "k+" << column[0].second << '\t';
        }
        std::cout << '\n';
    }
    std::cout << std::endl;
}

scout::matrix scout::MatrixOperations::IntegerMatrixSubtraction(matrix m1, matrix const& m2)
{
    auto size = m1.size();
    if (size != m2.size())
    {
        throw std::invalid_argument("Illegal Matrix Subtraction");
    }
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < m1[i].size(); ++j)
        {
            if ((m1[i][j].empty() && !m2[i][j].empty()))
            {
                m1[i][j].emplace_back(1, 0);
                return m1;
            }
            if ((!m1[i][j].empty() && m2[i][j].empty()))
            {
                m1[i][j][0].first = -1;
                return m1;
            }
            if (m1[i][j].empty() || m2[i][j].empty())
            {
                continue;
            }
            m1[i][j][0].second -= m2[i][j][0].second;
        }
    }

    return m1;
}

scout::matrix scout::MatrixOperations::MatrixAddition(matrix m1, matrix const& m2)
{
    auto size = m1.size();
    if (size != m2.size())
    {
        throw std::invalid_argument("Illegal Matrix Addition");
    }
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            if (m1[i][j].empty() || m2[i][j].empty())
            {
                continue;
            }
            m1[i][j][0].first = m2[i][j][0].second;
        }
    }
    return m1;
}

scout::matrix scout::MatrixOperations::MatrixComposition(matrix const& m1, matrix const& m2, bool tighten)
{
    auto baseSize = m1.size();
    auto halfSize = baseSize / 2;
    auto size = baseSize + halfSize;
    matrix res(size);

    for (int i = 0; i < size; ++i)
    {
        std::vector<cell> row(size);
        for (int j = 0; j < size; ++j)
        {
            int i2 = i - (int)halfSize;
            int j2 = j - (int)halfSize;

            if (i < baseSize && j < baseSize)
            {
                for (auto const& term : m1[i][j])
                {
                    row[j].emplace_back(term);
                }
            }
            if (i2 >= 0 && j2 >= 0)
            {
                for (auto const& term : m2[i2][j2])
                {
                    row[j].emplace_back(term);
                }
            }
        }
        res[i] = row;
    }
    res = ParametricFloydWarshallAlgorithm(res, tighten);
    return res;
}

scout::matrix scout::MatrixOperations::CalcExtremalPaths(matrix m)
{
    auto baseSize = m.size() / 3 * 2;
    auto halfSize = baseSize / 2;
    matrix reduced(baseSize);

    for (int i = 0; i < baseSize; ++i)
    {
        for (int j = 0; j < baseSize; ++j)
        {
            int posI = i >= halfSize ? i + (int)halfSize : i;
            int posJ = j >= halfSize ? j + (int)halfSize : j;
            reduced[i].emplace_back(m[posI][posJ]);
        }
    }
    return reduced;
}
