#include "preProcessor.h"
#include <cmath>

void preProcessor::initializeParameters(caseParameters &params)
{
    params.numberOfPoints = 101;
    params.gamma = 1.4;
    params.domainLength = 1.0;
    params.endTime = 0.2;
    params.CFL = 0.1;

    params.dx = params.domainLength / (params.numberOfPoints - 1);
    params.time = 0.0;
    params.timeStep = 0;
}

void preProcessor::initializeMesh(std::vector<double> &x, const caseParameters &params)
{
    for (int i = 0; i < params.numberOfPoints; i++)
    {
        x[i] = i * params.dx;
    }
}

void preProcessor::initiazlieSolution(std::vector<std::array<double, 3>> &U,
                                      const std::vector<double> &x,
                                      const caseParameters &params)
{
    double rho = 0.0;
    double u = 0.0;
    double p = 0.0;

    for (int i = 0; i < params.numberOfPoints; i++)
    {
        if (x[i] < 0.5)
        {
            rho = 1.0;
            u = 0.0;
            p = 1.0;
        }
        else
        {
            rho = 0.125;
            u = 0.0;
            p = 0.1;
        }

        U[i][0] = rho;
        U[i][1] = rho * u;
        U[i][2] = p / (params.gamma - 1.0) + 0.5 * rho * std::pow(u, 2);
    }
}