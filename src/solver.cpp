#include "solver.h"
#include <cmath>
#include <algorithm>

namespace solver
{
    void advance(std::vector<std::array<double, 3>> &U,
                 std::vector<std::array<std::array<double, 3>, 2>> &Ufaces,
                 std::vector<std::array<std::array<double, 3>, 2>> &Ffaces,
                 const std::vector<double> &x,
                 caseParameters &params,
                 SCHEME scheme,
                 LIMITER limiter)
    {
        // Preparing solution update (store old soln and calc stable timestep)
        auto UOld = U;

        // calc stable time step
        double speedMax = 0.0;
        for (int i = 0; i < params.numberOfPoints; i++)
        {
            // Getting primitive var in order to calc wave speed based on speed of sound and local velo
            auto rho = U[i][0];
            auto u = U[i][1] / rho;
            auto p = (params.gamma - 1.0) * (U[i][2] - 0.5 * rho * std::pow(u, 2));

            // calculate wave speed for each cell
            double speedOfSound = std::sqrt(params.gamma * p / rho);
            if (speedOfSound + std::fabs(u) > speedMax)
            {
                speedMax = speedOfSound + std::fabs(u);
            }
        }
        double dt = (params.CFL * params.dx) / speedMax;

        // solve equations
        if (scheme == SCHEME::CONSTANT)
        {
            // piecewise constant reconstruction
            for (int i = 0; i < params.numberOfPoints; i++)
            {
                // use this temp "variable" to iterate through each of the conserved variables
                /* the value of the conserved variable at the centroid is copied to the right and left faces */
                for (int variable = 0; variable < 3; ++variable)
                {
                    Ufaces[i][FACE::WEST][variable] = U[i][variable];
                    Ufaces[i][FACE::EAST][variable] = U[i][variable];
                }
            }
        }
        else if (scheme == SCHEME::MUSCL)
        {
            // MUSCL scheme
            /* use lower-order scheme near boundaries. The 0 is at the left, 1 is at right boundary*/
            for (int variable = 0; variable < 3; ++variable)
            {
                Ufaces[0][FACE::WEST][variable] = U[0][variable];
                Ufaces[0][FACE::EAST][variable] = U[0][variable];

                Ufaces[params.numberOfPoints - 1][FACE::WEST][variable] = U[params.numberOfPoints - 1][variable];
                Ufaces[params.numberOfPoints - 1][FACE::EAST][variable] = U[params.numberOfPoints - 1][variable];
            }

            // use high-reso MUSCL Scheme on interior nodes / cells
            for (int i = 1; i < params.numberOfPoints - 1; i++)
            {
                for (int variable = 0; variable < 3; ++variable)
                {
                    auto du_i_plus_half = U[i + 1][variable] - U[i][variable];
                    auto du_i_minus_half = U[i - 1][variable] - U[i][variable];

                    double rL = du_i_minus_half / (du_i_plus_half + 1e-8);
                    double rR = du_i_plus_half / (du_i_minus_half + 1e-8);

                    double psiL = 1.0;
                    double psiR = 1.0;

                    // apply limiter to make scheme TCD (total variation diminishing)
                    if (limiter == LIMITER::MINMOD)
                    {
                        psiL = std::max(0.0, std::min(1.0, rL));
                        psiR = std::max(0.0, std::min(1.0, rR));
                    }
                    else if (limiter == LIMITER::OSHER)
                    {
                        double beta = 1; // Beta value arbitrarily selected??
                        psiL = std::max(0.0, std::min(rL, beta));
                        psiR = std::max(0.0, std::min(rR, beta));
                    }

                    else if (limiter == LIMITER::VANLEER)
                    {
                        psiL = (rL + std::fabs(rL)) / (1.0 + std::fabs(rL));
                        psiR = (rR + std::fabs(rR)) / (1.0 + std::fabs(rR));
                    }

                    Ufaces[i][FACE::WEST][variable] = U[i][variable] - 0.5 * psiL * du_i_plus_half;
                    Ufaces[i][FACE::EAST][variable] = U[i][variable] + 0.5 * psiR * du_i_minus_half;
                }
            }
        }

        // compute fluxes at faces
        std::array<double, 3> fluxL;
        std::array<double, 3> fluxR;

        for (int i = 1; i < params.numberOfPoints - 1; i++)
        {
            for (int face = FACE::WEST; face <= FACE::EAST; ++face)
            {
                int indexOffset = 0;
                if (face == FACE::WEST)
                    indexOffset = 0;
                else if (face == FACE::EAST)
                    indexOffset = 1;

                // get primitive vars at east and west face for each cell
                auto rhoL = Ufaces[i - 1 + indexOffset][FACE::EAST][0];
                auto uL = Ufaces[i - 1 + indexOffset][FACE::EAST][1] / rhoL;
                auto EL = Ufaces[i - 1 + indexOffset][FACE::EAST][2];
                auto pL = (params.gamma - 1.0) * (EL - 0.5 * rhoL * std::pow(uL, 2));
                auto aL = std::sqrt(params.gamma * pL / rhoL);

                auto rhoR = Ufaces[i - 1 + indexOffset][FACE::WEST][0];
                auto uR = Ufaces[i - 1 + indexOffset][FACE::WEST][1] / rhoR;
                auto ER = Ufaces[i - 1 + indexOffset][FACE::WEST][2];
                auto pR = (params.gamma - 1.0) * (ER - 0.5 * rhoR * std::pow(uR, 2));
                auto aR = std::sqrt(params.gamma * pR / rhoR);

                // compute flux vector
                fluxL[0] = rhoL * uL;
                fluxL[1] = pL + rhoL * std::pow(uL, 2);
                fluxL[2] = uL * (EL + pL);

                fluxR[0] = rhoR * uR;
                fluxR[1] = pR + rhoR * std::pow(uR, 2);
                fluxR[2] = uR * (ER + pR);

                // Rusanov Riemannn solver
                auto speedMax = std::max(std::fabs(uL) + aL, std::fabs(uR) + aR); // calc wavespeed max
                for (int variable = 0; variable < 3; ++variable)
                {
                    const auto &qL = Ufaces[i - 1 + indexOffset][FACE::EAST][variable];
                    const auto &qR = Ufaces[i + indexOffset][FACE::WEST][variable];
                    const auto &fL = fluxL[variable];
                    const auto &fR = fluxR[variable];
                    Ffaces[i][face][variable] = 0.5 * (fL + fR) - speedMax * (qR - qL);
                }
            }
        }

        // integrate eqn in time (IN LISTING THIS WASNT OPENED {})
        // this doesnt integrate at boundary points
        for (int i = 1; i < params.numberOfPoints - 1; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                const auto &dF = Ffaces[i][FACE::EAST][j] - Ffaces[i][FACE::WEST][j];
                U[i][j] = UOld[i][j] - (dt / params.dx) * dF;
            }
        }

        // update boundary conditions
        // dirichlet for velo, neumann for everything else
        // neumann is basically copying whatever value is in adjacent cell onto bounadry
        // i.e partial (phi)/ partial (n) = 0
        auto rhoL = 1.0;
        auto uL = 0.0;

        auto rhoR = 0.125;
        auto uR = 0.0;

        U[0][0] = U[1][0];
        U[0][1] = rhoL * uL;
        U[0][2] = U[1][2];

        U[params.numberOfPoints - 1][0] = U[params.numberOfPoints - 2][0];
        U[params.numberOfPoints - 1][1] = rhoR * uR;
        U[params.numberOfPoints - 1][2] = U[params.numberOfPoints - 2][2];

        params.time += dt;
        params.timeStep++;
    }
}
