#pragma once
#include <vector>
#include <array>
#include "Parameters.h"

namespace postProcessor
{
    void writeSolution(const std::vector<std::array<double, 3>> &U,
                       const std::vector<double> &x,
                       const caseParameters &params);
    void printProgress(const caseParameters &params);
}
