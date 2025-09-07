#include "postProcessor.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <cmath>

void postProcessor::writeSolution(
    const std::vector<std::array<double, 3>> &U, const std::vector<double> &x,
    const caseParameters &params)
{
    std::filesystem::create_directory("../solnData");

    // convert times step and points into 6 digit strings with leading zeros
    std::ostringstream timeStepTemp, pointsTemp;

    timeStepTemp << std::setfill('0') << std::setw(6);
    timeStepTemp << params.timeStep;
    auto timeStep = timeStepTemp.str();

    pointsTemp << std::setfill('0') << std::setw(6);
    pointsTemp << params.numberOfPoints;
    auto points = pointsTemp.str();

    // building file name string
    std::string filename = "./solnData/solution_" + points + "_" + timeStep + ".csv";

    // open file
    std::ofstream outputFile(filename);

    // check to see if output file is open
    if (!outputFile.is_open())
    {
        throw std::runtime_error("Congratulations! Your file: " + filename + " failed to open!");
    }

    outputFile << "x,rho,u,p" << std::endl;

    for (int i = 0; i < params.numberOfPoints; i++)
    {
        auto rho = U[i][0];
        auto u = U[i][1] / rho;
        auto p = (params.gamma - 1.0) * (U[i][2] - 0.5 * rho * std::pow(u, 2));
        outputFile << x[i] << "," << rho << "," << u << "," << p << std::endl;
    }
}

void postProcessor::printProgress(const caseParameters &params)
{
    std::cout << "Current time: " << std::scientific << std::setw(10) << std::setprecision(3) << params.time;
    std::cout << ", End time: " << std::scientific << std::setw(10) << std::setprecision(3) << params.endTime;
    std::cout << ", Current time step: " << std::fixed << std::setw(7) << params.timeStep;
    std::cout << "\r";
}