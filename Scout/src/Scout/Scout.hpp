#pragma once

#include "Parser.hpp"
#include "Relation.hpp"

/* Usage Example:
 *
 *   scout::Relation r = scout::Parser::RetrieveRelation(filePath);
 *   r.CalculateTransitiveClosure();
 *   r.PrintTransitiveClosure();
 */
