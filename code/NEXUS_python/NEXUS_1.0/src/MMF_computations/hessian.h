#ifndef HESSIAN_HEADER
#define HESSIAN_HEADER


#include <defines.h>
#include <array.h>
#include <vector.h>
#include <FFTW.h>


typedef Array<Real,NO_DIM>             ArrayReal3D;
typedef Vector<Real,NO_DIM>            VectorReal3D;    // real vector with 3 coordinates
typedef Vector<Real,NO_DIM*(NO_DIM-1)> VectorReal9D;    // real vector with 6 coordinates (in 3D)
typedef Array<VectorReal3D,NO_DIM>     ArrayVector3D;
typedef Array<VectorReal9D,NO_DIM>     ArrayVector9D;
typedef Array<FFTW_COMPLEX,NO_DIM>     ArrayComplex3D;



// returns the a 1D gaussian in Fourier space
void gaussianFilterFT(Real result[],
                      size_t const gridSize,
                      Real const radius,
                      Real const boxLength);


// computes the Hessian eigenvalues and eigenvectors for the FT of a density field. It also smooths the input field at radius 'r'.
void hessianMatrix(ArrayComplex3D &density,
                   ArrayVector3D &eigenvalues,
                   ArrayVector9D &eigenvectors,
                   Real const radius,
                   Real const sideLength[],
                   Real const bias,
                   bool const computeEigenvectors);

// computes the Hessian eigenvalues and eigenvectors for the SMOOTHED density field. NOTE: It does not smooth the field.
void hessianMatrix(ArrayReal3D &density,
                   ArrayVector3D &eigenvalues,
                   ArrayVector9D &eigenvectors,
                   Real const radius,
                   Real const sideLength[],
                   Real const bias,
                   bool const computeEigenvectors);


// computes the actual entries of the Hessian matrix. NOTE: this function smooths the input density field with a Gaussian filter.
void computeHessian(ArrayComplex3D &density,
                    ArrayVector9D &eigenvectors,
                    Real const radius,
                    Real const sideLength[],
                    Real const bias);

// computes each of the 6-independent entries of the Hessian matrix
void hessianEntry(FFTW_COMPLEX density[],
                  Real gauss[],
                  Int const gridSize[],
                  Real const radius,
                  Real const sideLength[],
                  Real const bias,
                  int const noDerivativeX, int const noDerivativeY, int const noDerivativeZ,
                  FFTW_COMPLEX result[]);


// return the eigenvalues and eigenvectors of the Hessian matrix
void eigenvaluesEigenvectors(ArrayVector3D &eigenvalues,
                             ArrayVector9D &eigenvectors,
                             bool const computeEigenvectors);


// applies a Gaussian filter to a density field
void gaussianFilter(ArrayReal3D &density,
                    Real const radius,
                    Real const sideLength[],
                    ArrayReal3D *result);



// computes the eigenvalues and eigenvectors for the 'reduced' tidal field (=\phi_{ij} - 1/3\delta). The input field 'density' is the gravitational potential.
void reducedTidalFieldEigenvectors(ArrayReal3D &density,
                                   ArrayVector3D &eigenvalues,
                                   ArrayVector9D &eigenvectors,
                                   Real const radius,
                                   Real const sideLength[],
                                   Real const bias,
                                   bool const computeEigenvectors);


#endif
