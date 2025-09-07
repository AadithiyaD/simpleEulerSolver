#include "parameters.h"
#include "preProcessor.h"
#include "solver.h"
#include "postProcessor.h"
#include <iostream>

int main()
{
    // create the param struct and initialize them with default values
    caseParameters params;
    preProcessor::initializeParameters(params);

    // the below variables are not covered in initializeParameters, so we do them here
    // Centroids of our cells
    std::vector<double> x(params.numberOfPoints);

    /*  The conserved vector U. This is a vector of arrays, each containing 3 elements:
       rho, rho * U, E */
    std::vector<std::array<double, 3>> U(params.numberOfPoints);

    /* The conserved variables at each of our cell faces
    It is a vector of arrays of 2 elements, representing the left and right faces of the current mesh cell in this 1D case
    Each sub element represents one of the 3 conserved variables in U */
    std::vector<std::array<std::array<double, 3>, 2>> Ufaces(params.numberOfPoints);

    /* The conserved fluxes at each cell face, similar explanation of strcuture to Ufaces */
    std::vector<std::array<std::array<double, 3>, 2>> Ffaces(params.numberOfPoints);

    // Ini mesh and soln
    preProcessor::initializeMesh(x, params);
    preProcessor::initiazlieSolution(U, x, params);

    // specify reconstruction scheme and limiter we want to use
    auto scheme = SCHEME::MUSCL;
    auto limiter = LIMITER::OSHER;

    while (params.time < params.endTime)
    {
        solver::advance(U, Ufaces, Ffaces, x, params, scheme, limiter);
        postProcessor::writeSolution(U, x, params);
        postProcessor::printProgress(params);
    }

    std::cout << "\n Simulation finished \n"
              << std::endl;
    return 0;
}