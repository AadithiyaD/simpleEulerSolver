#pragma once

// global enum for easy variable access

/* Specfiy which reconstruction scheme to use.
 Can be CONSTANT, MUSCL. Const => piecewise reconstruction */
enum SCHEME
{
    CONSTANT = 0,
    MUSCL
};

/* Flux limiter specification. Can be NONE, MINMOD, VANLEER, OSHER */
enum LIMITER
{
    NONE = 0,
    MINMOD,
    VANLEER,
    OSHER
};

// Specifies the face currently under consideration
enum FACE
{
    WEST = 0,
    EAST
};

// definition for case parameter structure to hold case-specific strings
struct caseParameters
{
    int numberOfPoints;
    double gamma;
    double domainLength;
    double endTime;
    double CFL;
    double dx;
    double time;
    int timeStep;
};