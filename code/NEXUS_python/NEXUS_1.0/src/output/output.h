#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <list>
#include <utility>

#include <gadget_header.h>
#include <density_header.h>
#include <miscellaneous.h>
#include <hopGroup.h>
#include <array.h>
#include <defines.h>



//! Template functions
#include "raw.h"



//! Functions in "visivo.cc"
void outputVisivo(float *delta,
                  Int const gridSize[3],
                  std::string rootName);

void logarithm(Real *delta,
               Int const totalSize,
               Real const delta0,
               Real const minimumDelta);


//! Functions in "ascii.cc"
void outputASCII_matrix(Real *densityMatrix,
                        Int const nGrid[],
                        std::string outputFileName,
                        double const limitDown,
                        double const limitUp,
                        std::string const inputFileName);

void outputASCII_Haloes(std::list<HopGroup> &blobHaloes,
                        std::string outputFileName,
                        std::string descriptionOfData,
                        std::string programOptions);
void outputASCII_Haloes(HopGroup *blobHaloes,
                        Int const noGroups,
                        std::string outputFileName,
                        std::string descriptionOfData,
                        std::string programOptions);
void outputASCII_Haloes(std::list<HopGroup> &blobHaloes,
                        std::string outputFileName);
void outputASCII_Haloes(HopGroup *blobHaloes,
                        Int const noGroups,
                        std::string outputFileName);

void outputASCII_HaloesPositions(std::list<HopGroup> &blobHaloes,
                                 std::string outputFileName);
void outputASCII_HaloesPositions(HopGroup *blobHaloes,
                                 Int const noGroups,
                                 std::string outputFileName);


//! Functions in "densitySlices.cc"
void densitySliceComputation(Real *delta3D,
                             Real *d2D,
                             Int const nGrid[3],
                             int const axis,
                             Real posMin,
                             Real posMax,
                             Real const boxLength);

Int compute_2D_index(Int const n1,
                     Int const n2,
                     Int n3,
                     int const choice,
                     Int const nGrid[]);
int perpendicularDirection(int const direction,
                           int const choice);




//! Functions in "compareDensityOutput.cc"
void outputDensityHistogram(Array<Real,3> &density,
                            std::vector<Real> range,
                            std::string outputFilename,
                            std::string const programOptions);

void outputPointRatioDistribution(Array<Real,3> &density,
                                  Array<Real,3> &pointRatios,
                                  Real const invalid,
                                  std::vector<Real> xRange,
                                  std::vector<Real> yRange,
                                  std::string outputFilename,
                                  std::string const programOptions);

void outputDensityRatio(Array<Real,3> &density,
                        Array<Real,3> &mainDensity,
                        std::vector<Real> range,
                        std::vector<Real> ratioValues,
                        std::string outputFilename,
                        std::string const programOptions);

void outputPointRatio(Array<Real,3> &density,
                      Array<Real,3> &pointRatios,
                      Real const invalid,
                      std::string outputFilename,
                      std::string const programOptions);

void outputPointStatistics(Array< std::pair< std::pair<Real,size_t> ,std::pair<Real,Real> >, 1 >  &statistics,
                           std::string outputFilename,
                           std::string const programOptions);

void outputcellCount(Array<Real,3> &density,
                     Array<int,1> &cellCountPoints,
                     std::vector<Real> range,
                     std::string outputFilename,
                     std::string const programOptions);









#endif


