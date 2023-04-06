#include "Scout/Scout.hpp"

int main()
{
    scout::Relation r = scout::Parser::RetrieveRelation("dataset/1567523039294457.koat.rel");
    r.CalculateTransitiveClosure();
    r.PrintTransitiveClosure();
}
