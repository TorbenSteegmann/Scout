#pragma once

#include <string>
#include <vector>

namespace scout
{
constexpr bool MAKE_SMT = false;

// typedefs
typedef std::vector<std::pair<int, int>> cell;
typedef std::vector<std::vector<cell>> matrix;
} // namespace scout
