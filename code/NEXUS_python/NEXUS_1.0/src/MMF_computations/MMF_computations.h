#ifndef MMF_COMPUTATIONS_HEADER
#define MMF_COMPUTATIONS_HEADER
#include <list>

#include <defines.h>
#include <miscellaneous.h>
#include <array.h>
#include <vector.h>
#include <FFTW.h>
using namespace std;


typedef Array<Real,NO_DIM>             ArrayReal3D;
typedef Vector<Real,NO_DIM>            VectorReal3D;    // real vector with 3 coordinates
typedef Vector<Real,NO_DIM*(NO_DIM-1)> VectorReal9D;    // real vector with 6 coordinates (in 3D)
typedef Array<VectorReal3D,NO_DIM>     ArrayVector3D;
typedef Array<VectorReal9D,NO_DIM>     ArrayVector9D;
typedef Array<FFTW_COMPLEX,NO_DIM>     ArrayComplex3D;


#include "threshold.h"
#include "hessian.h"
#include "objects.h"
#include "response.h"
#include "additionalComputations.h"



//! Functions in "eigenvectors.cc"
void rescaleEigenvalues(ArrayVector3D *eigenvalues,
                        Real const radius,
                        Real const bias);

void eigenMaximumResponse(Array<int,3> const &scale,
                          ArrayVector3D *eigenvalues,
                          ArrayVector9D *eigenvectors,
                          ArrayVector3D &eigenvalues2,
                          ArrayVector9D &eigenvectors2,
                          int const scale2);

void featureEigenvector(ArrayVector9D &eigenvectors,
                        ArrayVector3D *eigenvalues,
                        int const feature);

VectorReal3D lastEigenvector(VectorReal9D &eigenvectors);



//! Functions in 'filters.c'
void computeResponse(ArrayVector3D &eigenvalue,
                     ArrayReal3D &response,
                     int const filterNumber);

string filterDescription(int const filterNumber);








#endif

