#pragma once
#include <vector>
#include <array>
#include "Parameters.h"

namespace solver
{
    /* Start the solution. We first specify the conserved variable vector U, then U at faces
    , conserved fluxes at faces, number of gridpoints, case params, reconstruction scheme and
    flux limiter*/
    void advance(std::vector<std::array<double, 3>> &U,
                 std::vector<std::array<std::array<double, 3>, 2>> &Ufaces,
                 std::vector<std::array<std::array<double, 3>, 2>> &Ffaces,
                 const std::vector<double> &x,
                 caseParameters &params,
                 SCHEME scheme,
                 LIMITER limiter);
}
