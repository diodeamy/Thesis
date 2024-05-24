#ifndef ADDITIONAL_COMPUTATIONS_HEADER
#define ADDITIONAL_COMPUTATIONS_HEADER


#include <array.h>
#include <k_space.h>
#include <FFTW.h>
#include <MMF_header.h>

typedef Array<Real,NO_DIM>             ArrayReal3D;
typedef Array<FFTW_COMPLEX,NO_DIM>     ArrayComplex3D;



void densityLogarithm(ArrayReal3D *result,
                      Real const minValue,
                      bool VERBOSE = true);
void densityLogarithm(ArrayReal3D &input,
                      Real const minValue,
                      ArrayReal3D *result,
                      bool VERBOSE = true);
void densityExponential(ArrayReal3D *result,
                        Real const factor,
                        bool VERBOSE = true);

void computeFT_gravitationalField(ArrayComplex3D *result,
                                  size_t *grid,
                                  Real const boxLength[],
                                  Real const factor);

void solvePoisson(ArrayReal3D *result,
                  Real boxLength[],
                  Real const factor,
                  string variableName);

void velocityDivergenceLogarithm(ArrayReal3D *result,
                                 Real const minValue);
void velocityDivergenceLogarithm(ArrayReal3D *result,
                                 Real const minValue,
                                 Real const offset);

void maskInputField(ArrayReal3D *result,
                    MMF_header &mmfHeader,
                    Real maskValue,
                    string maskFilename,
                    string maskName);


#endif
