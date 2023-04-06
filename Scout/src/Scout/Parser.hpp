#pragma once

#include <fstream>
#include <iostream>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "Common.hpp"
#include "Relation.hpp"

namespace scout
{

struct variable
{
    std::optional<std::string> name;
    std::optional<int> number;
    std::optional<bool> primed;
    int factor;
};

struct token
{
    enum tokenType
    {
        VAR,
        CONST,
        PLUS,
        MINUS,
        MUL,
        DIV,
        SMALLER,
        GREATER,
        EQUALS,
        SMALLER_EQ,
        GREATER_EQ,
        LAND
    } type{};
    std::optional<variable> var;
};

typedef std::vector<variable> conjunct;

namespace Parser
{
// Wrapper function of the parser
Relation RetrieveRelation(std::string const& filePath);

// Filters out the relation in between first : and ; of a file
std::string ReadRelation(std::string const& filePath);

// Converts string into token format
std::vector<token> TokenizeRelation(std::string const& relation, Relation& r);

// helper functions for TokenizeRelation
std::pair<std::vector<token>, std::string> ConsumeToken(std::pair<std::vector<token>, std::string> tokenStream, Relation& r);

std::pair<std::vector<token>, std::string> ConsumeVarOrConst(std::pair<std::vector<token>, std::string> tokenStream, Relation& r);

std::pair<std::vector<token>, std::string> ConsumeSymbol(std::pair<std::vector<token>, std::string> tokenStream);

// first step of turning tokens into proper form
std::vector<conjunct> MakeFormula(std::vector<token> const& tokens);

// helper function for MakeFormula
conjunct NegateFormula(conjunct f);

std::vector<conjunct> AddTokens(std::vector<conjunct> const& tokenizedFormula);

// second step of turning tokens into proper form
std::vector<conjunct> NormalizeTokens(std::vector<conjunct> const& tokenizedFormula);

// check if the resulting form is octagonal, dbr or invalid
bool VerifyValidity(std::vector<conjunct> const& resolvedTokens);

// converts tokens into matrix
void MakeRelation(std::vector<conjunct> const& tokenizedFormula, Relation& r);

// prints tokens for debugging
[[maybe_unused]] static void PrintTokens(std::vector<conjunct> const& tokenizedRelation);

} // namespace Parser
} // namespace scout
