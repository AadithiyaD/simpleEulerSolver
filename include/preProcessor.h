#pragma once
#include <vector>
#include <array>
#include "parameters.h"

/*  the stuff where we've writte caseParameters&, std::vecotr<double> x , are like type hinting
    When I hover over these functions, I can see the variables i need to give, and their expected data type
    Pretty neat*/

namespace preProcessor
{
    void initializeParameters(caseParameters &params);
    void initializeMesh(std::vector<double> &x, const caseParameters &params);
    void initiazlieSolution(std::vector<std::array<double, 3>> &U,
                            const std::vector<double> &x,
                            const caseParameters &params);
} // namespace preProcessor
