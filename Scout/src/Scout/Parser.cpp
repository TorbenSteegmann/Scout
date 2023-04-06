#include "Parser.hpp"
#include "MatrixOperations.hpp"

#include <string_view>

// Allowed Symbols
static constexpr std::string_view VARIABLE_NAME_SYMBOLS = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_$";
static constexpr std::string_view DIGITS = "1234567890";
static constexpr std::string_view LOGIC_SYMBOLS = "+-*/<>=&";

scout::Relation scout::Parser::RetrieveRelation(std::string const& filePath)
{
    Relation r;
    auto relationAsString = Parser::ReadRelation(filePath);
    auto tokens = Parser::TokenizeRelation(relationAsString, r);
    auto tokenizedFormula = Parser::MakeFormula(tokens);
    tokenizedFormula = Parser::AddTokens(tokenizedFormula);
    tokenizedFormula = Parser::NormalizeTokens(tokenizedFormula);
    r.SetIsOctagonal(Parser::VerifyValidity(tokenizedFormula));
    Parser::MakeRelation(tokenizedFormula, r);
    return r;
}

std::string scout::Parser::ReadRelation(std::string const& filePath)
{
    std::ifstream file;
    file.open(filePath);
    if (!file)
    {
        throw std::invalid_argument("No such File exists");
    }
    std::string line, input, relation;
    while (getline(file, line))
    {
        input.append(line);
    }
    auto pos1 = input.find(':');
    if (pos1 != std::string::npos)
    {
        for (auto pos2 = pos1 + 1; pos2 + 1 != std::string::npos; ++pos2)
        {
            if (input[pos2] == ';')
            {
                break;
            }
            if (input[pos2] != ' ')
            {
                relation += input[pos2];
            }
        }
    }
    return relation;
}

std::vector<scout::token> scout::Parser::TokenizeRelation(std::string const& relation, Relation& r)
{
    std::pair<std::vector<token>, std::string> tokenStream = std::make_pair(std::vector<token>{}, relation);
    return Parser::ConsumeToken(tokenStream, r).first;
}

std::pair<std::vector<scout::token>, std::string> scout::Parser::ConsumeToken(std::pair<std::vector<token>, std::string> tokenStream, Relation& r)
{
    if (tokenStream.second.length() == 0)
    {
        return tokenStream;
    }

    if (LOGIC_SYMBOLS.find(tokenStream.second[0]) != std::string::npos)
    {
        tokenStream = Parser::ConsumeSymbol(tokenStream);
    }
    else
    {
        tokenStream = Parser::ConsumeVarOrConst(tokenStream, r);
    }

    return Parser::ConsumeToken(tokenStream, r);
}


std::pair<std::vector<scout::token>, std::string> scout::Parser::ConsumeVarOrConst(std::pair<std::vector<token>, std::string> tokenStream, Relation& r)
{
    std::string numberOrName;
    auto isVariable = false;
    auto primed = false;
    auto succeeded = false;

    for (auto character : tokenStream.second)
    {
        if (DIGITS.find(character) != std::string::npos)
        {
            numberOrName += character;
        }
        else if (VARIABLE_NAME_SYMBOLS.find(character) != std::string::npos)
        {
            numberOrName += character;
            isVariable = true;
        }
        else if (character == '\'')
        {
            if (primed || !isVariable)
            {
                throw std::invalid_argument("More than one \' in Variable or pure numbered Name.");
            }
            primed = true;
        }
        else if (LOGIC_SYMBOLS.find(character) != std::string::npos)
        {
            succeeded = true;
        }
        else
        {
            throw std::invalid_argument("Character not allowed: ");
        }
        if (succeeded)
        {
            break;
        }
    }
    if (numberOrName.length() == 0)
    {
        throw std::invalid_argument("Parsing Error: Empty Variable Name");
    }
    if (!isVariable)
    {
        auto number = std::stoi(numberOrName);

        tokenStream.first.emplace_back(token{token::CONST, variable{.factor = number}});
        tokenStream.second.erase(0, numberOrName.length());
        return tokenStream;
    }

    // has to be a isVariable here
    int number = 0;
    for (auto const& pair : r.GetVariableMap())
    {
        if (pair.second == numberOrName)
        {
            number = pair.first;
        }
    }
    if (number == 0)
    {
        auto newMap = r.GetVariableMap();
        number = (int)newMap.size() + 1;
        newMap.emplace(number, numberOrName);
        r.SetVariableMap(newMap);
    }

    tokenStream.first.emplace_back(token{token::VAR, variable{numberOrName, number, primed, 1}});
    tokenStream.second.erase(0, numberOrName.length() + primed);

    return tokenStream;
}

std::pair<std::vector<scout::token>, std::string> scout::Parser::ConsumeSymbol(std::pair<std::vector<token>, std::string> tokenStream)
{
    int eraseLength = 0;
    auto failed = false;
    for (auto i = 0; i < tokenStream.second.length() && !failed; ++i)
    {
        switch (tokenStream.second[i])
        {
        case '+':
            tokenStream.first.emplace_back(token{token::PLUS, variable{.factor = 1}});
            break;
        case '-':
            tokenStream.first.emplace_back(token{token::MINUS});
            break;
        case '*':
            tokenStream.first.emplace_back(token{token::MUL});
            break;
        case '/':
            tokenStream.first.emplace_back(token{token::DIV});
            break;
        case '<':
            if (tokenStream.second[i + 1] == '=')
            {
                tokenStream.first.emplace_back(token{token::SMALLER_EQ});
                ++i;
                ++eraseLength;
            }
            else
            {
                tokenStream.first.emplace_back(token{token::SMALLER});
            }
            break;
        case '>':
            if (tokenStream.second[i + 1] == '=')
            {
                tokenStream.first.emplace_back(token{token::GREATER_EQ});
                ++i;
                ++eraseLength;
            }
            else
            {
                tokenStream.first.emplace_back(token{token::GREATER});
            }
            break;
        case '=':
            tokenStream.first.emplace_back(token{token::EQUALS});
            break;
        case '&':
            if (tokenStream.second[i + 1] == '&')
            {
                tokenStream.first.emplace_back(token{token::LAND});
                ++i;
                ++eraseLength;
            }
            else
            {
                throw std::invalid_argument("No singular \'&\' allowed");
            }
            break;
        default:
            failed = true;
            break;
        }
        if (!failed)
        {
            ++eraseLength;
        }
    }
    if (eraseLength != 0)
    {
        tokenStream.second.erase(0, eraseLength);
    }
    return tokenStream;
}


std::vector<scout::conjunct> scout::Parser::MakeFormula(std::vector<token> const& tokens)
{
    std::vector<conjunct> rel;
    conjunct f;
    bool newFormula = false;
    int factor = 1;
    bool rightSide = false;
    bool addNormal = false;
    bool addInverted = false;
    bool multiplication = false;
    int rightSideFactor = 1;
    for (auto const& t : tokens)
    {
        if (newFormula)
        {
            if (addNormal)
            {
                rel.emplace_back(f);
            }
            if (addInverted)
            {
                rel.emplace_back(Parser::NegateFormula(f));
            }
            else if (!addNormal)
            {
                throw std::invalid_argument("There is no comparator in formula");
            }
            f = {};
            newFormula = false;
            factor = 1;
            rightSide = false;
            rightSideFactor = 1;
            addNormal = false;
            addInverted = false;
            multiplication = false;
        }
        if (t.type > 5 && t.type != token::LAND)
        {
            if (rightSide)
            {
                throw std::invalid_argument("Multiple Comparators");
            }
            rightSide = true;
        }
        if (rightSide)
        {
            rightSideFactor = -1;
        }
        switch (t.type)
        {
        case token::VAR:
            if (multiplication)
            {
                if (f.back().name)
                {
                    throw std::invalid_argument("No variable multiplication allowed");
                }
                f.back().name = t.var->name;
                f.back().number = t.var->number;
                f.back().primed = t.var->primed;
                f.back().factor *= -1 * factor * t.var->factor;
                factor = 1;
                multiplication = false;
                break;
            }
            f.emplace_back(variable{t.var->name, t.var->number, t.var->primed, rightSideFactor * factor * (t.var->factor)});
            factor = 1;
            break;
        case token::CONST:
            if (multiplication)
            {
                f.back().factor *= factor * t.var->factor;
                factor = 1;
                multiplication = false;
                break;
            }
            f.emplace_back(variable{.factor = -1 * rightSideFactor * factor * (t.var->factor)});
            factor = 1;
            break;
        case token::PLUS:
            continue;
        case token::MINUS:
            factor = -1 * factor;
            break;
        case token::MUL:
            multiplication = true;
            break;
        case token::DIV:
            throw std::invalid_argument("No division allowed");
        case token::SMALLER:
            f.emplace_back(variable{.factor = -1});
            rightSide = true;
            addNormal = true;
            break;
        case token::GREATER:
            f.emplace_back(variable{.factor = 1});
            rightSide = true;
            addInverted = true;
            break;
        case token::EQUALS:
            addNormal = true;
            addInverted = true;
            break;
        case token::SMALLER_EQ:
            addNormal = true;
            break;
        case token::GREATER_EQ:
            addInverted = true;
            break;
        case token::LAND:
            newFormula = true;
            break;
        }
    }
    if (addNormal)
    {
        rel.emplace_back(f);
    }
    if (addInverted)
    {
        rel.emplace_back(Parser::NegateFormula(f));
    }
    else if (!addNormal)
    {
        throw std::invalid_argument("There is no comparator in formula");
    }

    return rel;
}

scout::conjunct scout::Parser::NegateFormula(conjunct f)
{
    for (auto& var : f)
    {
        var.factor = -1 * var.factor;
    }
    return f;
}

void scout::Parser::PrintTokens(std::vector<conjunct> const& tokenizedRelation)
{
    for (auto const& conjunct : tokenizedRelation)
    {
        std::cout << "\nConjunct:";
        for (auto var : conjunct)
        {
            std::cout << " ";
            if (var.name)
            {
                std::cout << "(" << *var.name << " " << *var.number << " " << *var.primed << " " << var.factor << ")";
            }
            else
            {
                std::cout << "(" << var.factor << ")";
            }
        }
    }
    std::cout << std::endl;
}

std::vector<scout::conjunct> scout::Parser::AddTokens(std::vector<conjunct> const& tokenizedFormula)
{
    std::vector<conjunct> resultingFormula;
    conjunct resultingConjunct;
    std::set<std::pair<int, bool>> AddedValues;
    bool foundConstant;

    for (auto const& conjunct : tokenizedFormula)
    {
        AddedValues = {};
        resultingConjunct = {};
        foundConstant = false;
        for (int i = 0; i < conjunct.size(); ++i)
        {
            variable v = conjunct[i];
            if (v.number)
            {
                if (AddedValues.count(std::make_pair(*v.number, *v.primed)) != 0)
                { // doesnt respect prime
                    continue;
                }
                else
                {
                    AddedValues.insert(std::make_pair(*v.number, *v.primed));
                }
            }
            else
            {
                if (AddedValues.count(std::make_pair(-1, 0)) != 0)
                {
                    continue;
                }
                else
                {
                    AddedValues.insert(std::make_pair(-1, 0));
                    foundConstant = true;
                }
            }
            for (int j = i + 1; j < conjunct.size(); ++j)
            {
                if (v.name == conjunct[j].name && v.primed == conjunct[j].primed)
                {
                    v.factor += conjunct[j].factor;
                }
            }
            if (v.factor != 0 || !v.name)
            {
                resultingConjunct.emplace_back(v);
            }
        }
        if (!foundConstant)
        {
            resultingConjunct.emplace_back(variable{.factor = 0});
        }
        resultingFormula.emplace_back(resultingConjunct);
    }


    return resultingFormula;
}

bool scout::Parser::VerifyValidity(std::vector<conjunct> const& resolvedTokens)
{
    bool octagon = true;
    bool dbr = true;
    for (auto const& conjunct : resolvedTokens)
    {
        int numberOfVariables = 0;
        int numberOfPositiveVariables = 0;
        for (auto const& var : conjunct)
        {
            if (var.name)
            {
                if (var.factor > 0)
                {
                    numberOfPositiveVariables++;
                }
                numberOfVariables++;
            }
        }
        switch (numberOfVariables)
        {
        case 1:
            dbr = false;
            break;
        case 2:
            if (numberOfPositiveVariables != 1)
            {
                dbr = false;
            }
            break;
        default:
            octagon = false;
            dbr = false;
            break;
        }
    }
    if (!octagon)
    {
        throw std::invalid_argument("Invalid Formula");
    }


    return !dbr;
}

std::vector<scout::conjunct> scout::Parser::NormalizeTokens(std::vector<conjunct> const& tokenizedFormula)
{
    std::vector<conjunct> unifiedFormulas;
    conjunct unifiedFormula;
    variable var1;
    variable var2;
    variable constant;
    int highestFactor;

    for (auto const& conjunct : tokenizedFormula)
    {
        unifiedFormula = {};
        var1 = {};
        var2 = {};
        highestFactor = 1;
        for (auto const& token : conjunct)
        {
            if (token.name)
            {
                var2 = var1.name ? token : var2;
                var1 = var1.name ? var1 : token;
                highestFactor = std::max(std::abs(highestFactor), std::abs(token.factor));
            }
            else
            {
                constant = token;
            }
        }
        if (!var2.name)
        {
            var2 = var1;
            constant.factor *= 2;
        }
        var1.factor = var1.factor % highestFactor == 0 ? var1.factor / highestFactor : throw std::invalid_argument("Can't Normalize Variable Factors.");
        var2.factor = var2.factor % highestFactor == 0 ? var2.factor / highestFactor : throw std::invalid_argument("Can't Normalize Variable Factors.");
        constant.factor = constant.factor / highestFactor;

        unifiedFormula.emplace_back(var1);
        unifiedFormula.emplace_back(var2);
        unifiedFormula.emplace_back(constant);

        unifiedFormulas.emplace_back(unifiedFormula);
    }

    return unifiedFormulas;
}

void scout::Parser::MakeRelation(std::vector<conjunct> const& tokenizedFormula, Relation& r)
{
    int size = r.GetIsOctagonal() ? 4 * (int)r.GetVariableMap().size() : 2 * (int)r.GetVariableMap().size();
    std::vector<cell> row(size);
    matrix m;
    for (int i = 0; i < size; ++i)
    {
        m.emplace_back(row);
    }
    for (auto const& conjunct : tokenizedFormula)
    {
        if (!r.GetIsOctagonal())
        {
            int posI, posJ;
            if (conjunct[0].factor > 0)
            {
                posI = *conjunct[0].number + *conjunct[0].primed * size / 2 - 1;
                posJ = *conjunct[1].number + *conjunct[1].primed * size / 2 - 1;
            }
            else
            {
                posJ = *conjunct[0].number + *conjunct[0].primed * size / 2 - 1;
                posI = *conjunct[1].number + *conjunct[1].primed * size / 2 - 1;
            }
            auto weight = {std::make_pair(0, conjunct[2].factor)};
            m[posI][posJ] = weight;
            continue;
        }

        int posIOne = conjunct[0].factor > 0 ? 2 * *conjunct[0].number - 2 + *conjunct[0].primed * size / 2
                                             : 2 * *conjunct[0].number - 1 + *conjunct[0].primed * size / 2;
        int posITwo = conjunct[1].factor > 0 ? 2 * *conjunct[1].number - 2 + *conjunct[1].primed * size / 2
                                             : 2 * *conjunct[1].number - 1 + *conjunct[1].primed * size / 2;
        int posJOne = conjunct[1].factor > 0 ? 2 * *conjunct[1].number - 1 + *conjunct[1].primed * size / 2
                                             : 2 * *conjunct[1].number - 2 + *conjunct[1].primed * size / 2;
        int posJTwo = conjunct[0].factor > 0 ? 2 * *conjunct[0].number - 1 + *conjunct[0].primed * size / 2
                                             : 2 * *conjunct[0].number - 2 + *conjunct[0].primed * size / 2;

        auto weight = {std::make_pair(0, conjunct[2].factor)};

        m[posIOne][posJOne] = weight;
        m[posITwo][posJTwo] = weight;
    }
    for (int i = 0; i < size; ++i)
    {
        m[i][i] = {std::make_pair(0, 0)};
    }
    //    MatrixOperations::PrintMatrix(m);
    //    r.transitiveClosure = {m};
    //    r.PrintTransitiveClosure();
    //    std::cout << "here" << std::endl;
    m = MatrixOperations::ParametricFloydWarshallAlgorithm(m, r.GetIsOctagonal());
    //    MatrixOperations::PrintMatrix(m);
    // r.settest(m);
    // std::cout << std::endl;

    r.AddPowerOfRelation(1, m);
}
