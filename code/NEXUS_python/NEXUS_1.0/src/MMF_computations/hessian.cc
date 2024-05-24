#include <omp.h>
#include <cmath>
#include <ctime>


#include "hessian.h"
#include <miscellaneous.h>
#include <matrix.h>
#include <k_space.h>
using namespace std;



/* Computes the 1D Gaussian in FT-space. */
void gaussianFilterFT(Real result[],
                      size_t const gridSize,
                      Real const radius,
                      Real const boxLength)
{
    K_space<Real,1> momentum( gridSize, boxLength );
    Real factor = -2.*PI*PI / (boxLength*boxLength) * radius*radius; // -1/2 * momentum^2 * radius^2
    
    for (size_t i=0; i<gridSize; ++i)
    {
        Real temp = momentum.nx(i);
        result[i] = exp( factor * temp*temp );
    }
}




/* Computes:
   the Hessian matrix
   the eigenvalues and eigenvectors of the Hessian
*/
void hessianMatrix(ArrayComplex3D &density,
                   ArrayVector3D &eigenvalues,
                   ArrayVector9D &eigenvectors,
                   Real const radius,
                   Real const sideLength[],
                   Real const bias,
                   bool const computeEigenvectors)
{
    Timer t;
    
    // compute the Hessian and save its values in the array 'eigenvectors'
    t.start();
    computeHessian( density, eigenvectors, radius, sideLength, bias );  //! works only for 3 spatial dimensions
    t.printTime("hessian");

    
    // compute the eigenvalues and eigenvectors
    t.start();
    eigenvaluesEigenvectors( eigenvalues, eigenvectors, computeEigenvectors ); //! works for any dimension of the space
    t.printTime("eigenvalues");
}


/* This function computes the eigenvalues and eigenvectors of an array of Hessian matrices.
*/
void eigenvaluesEigenvectors(ArrayVector3D &eigenvalues,
                             ArrayVector9D &eigenvectors,
                             bool const computeEigenvectors)
{
    cout << "Computing the eigenvalues and eigenvectors of the Hessian matrix ... " << flush;
    
    
    Int const totalSize = eigenvalues.totalSize();
    MATRIX::Symmetric<Real,NO_DIM> mat;
    if ( computeEigenvectors )
    {
        #pragma omp parallel for private(mat)
        for (Int i=0; i<totalSize; ++i)
        {
            mat.data( &(eigenvectors[i][0]) );
            mat.eigenvalues( &(eigenvalues[i][0]), &(eigenvectors[i][0]), false ); // returns only the NO_DIM-1 eigenvectors, the remaining one can be computed from these values
        }
    }
    else // do not return eigenvectors - faster computation
    {
        #pragma omp parallel for private(mat)
        for (Int i=0; i<totalSize; ++i)
        {
            mat.data( &(eigenvectors[i][0]) );
            mat.eigenvalues( &(eigenvalues[i][0]), NULL, false ); // returns only the eigenvalues
        }
    }

    cout << "Done.\n";
}


/* Computes the values of the Hessian matrix by multiplying in Fourier space the FT of the Gaussian filter and its derivatives with the FT of the density/gravitational filed.
*/
void computeHessian(ArrayComplex3D &density,
                    ArrayVector9D &eigenvectors,
                    Real const radius,
                    Real const sideLength[],
                    Real const bias)
{
    intervalCheck( NO_DIM, 3, 3, "'NO_DIM' in function 'computeHessian'" );  // at the moment this method is implemented only for 3 dimensions
    cout << "Computing the entries of the Hessian matrix using 'bias'=" << bias << " ... " << flush;
    
    // define some constants
    Array<Int,1> gridSize = density.getSize<Int>();
    Array<Int,1> tempGrid = eigenvectors.getSize<Int>();
    Int const nx = eigenvectors.getSize(0);
    Int const ny = eigenvectors.getSize(1);
    Int const nz = eigenvectors.getSize(2);
    if ( nx!=gridSize(0) or ny!=gridSize(1) or (nz/2+1)!=gridSize(2) )
        throwError( "In function 'computeHessian' the variables 'density' and 'eigenvectors' do not have the same dimensions." );
    
    
    // define temporary variables to store the Gaussian filter along each direction
    Array<Real,1> gauss( nx+ny+nz );
    gaussianFilterFT( &(gauss[0]), nx, radius, sideLength[0] );          // get the values of the FT of the Gaussian filter
    gaussianFilterFT( &(gauss[nx]), ny, radius, sideLength[1] );
    gaussianFilterFT( &(gauss[nx+ny]), nz, radius, sideLength[2] );
    
    // create variables for the FFTW
    FFTW_class fftw( (eigenvectors.getSize<int>()).ptrData(), NO_DIM, FFTW_class::C2R, omp_get_max_threads() );
    fftw.plan();
    FFTW_COMPLEX *in = fftw.ptrComplexInput();  //points to the input complex data for the FFTW
    
    
    // Now compute each of the 6 independent entries of the Hessian matrix - the entries are saved in the order {h_xx, h_xy, h_xz, h_yy, h_yz, h_zz} in the memory location 0 to 5 starting at eigenvectors[i][0]
    hessianEntry( density.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 2, 0, 0, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][0]), NO_DIM*(NO_DIM-1) );  // writes the FFTW output in memory starting with '&(eigenvectors[0][0])' at every 'NO_DIM2' steps further
    
    hessianEntry( density.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 1, 1, 0, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][1]), NO_DIM*(NO_DIM-1) );
    
    hessianEntry( density.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 1, 0, 1, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][2]), NO_DIM*(NO_DIM-1) );
    
    hessianEntry( density.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 0, 2, 0, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][3]), NO_DIM*(NO_DIM-1) );
    
    hessianEntry( density.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 0, 1, 1, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][4]), NO_DIM*(NO_DIM-1) );
    
    hessianEntry( density.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 0, 0, 2, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][5]), NO_DIM*(NO_DIM-1) );
    
    
    cout << "Done.\n";
}


/* Computes FT of entry 00 of the Hessian. */
void hessianEntry(FFTW_COMPLEX density[],
                  Real gauss[],
                  Int const gridSize[],
                  Real const radius,
                  Real const sideLength[],
                  Real const bias,
                  int const noDerivativeX, int const noDerivativeY, int const noDerivativeZ,
                  FFTW_COMPLEX result[])
{
    if ( noDerivativeX<0 or noDerivativeY<0 or noDerivativeZ<0 ) throwError( "At least one of the 'noDerivativeX/Y/Z' argument of function 'hessianEntry' has a negative value. This is not allowed since those arguments give the number of derivatives applied along each axis." );
    if ( noDerivativeX+noDerivativeY+noDerivativeZ!=2 ) throwError( "The sum of the number of derivatives along all 3 directions must be 2 in function 'hessianEntry'. This wasn't the case now." );
    
    // define some constants
    Int const nx = gridSize[0];
    Int const ny = gridSize[1];
    Int const nz = gridSize[2]/2 + 1;
    K_space<Real,1> mX( gridSize[0], sideLength[0] ),
                    mY( gridSize[1], sideLength[1] ),
                    mZ( gridSize[2], sideLength[2] );
    Int totalGrid = gridSize[0] * gridSize[1] * gridSize[2];
    Real const factor = Real(-1.) / Real(totalGrid) * pow(double(radius),double(2.*bias));    // 'totalGrid' is the 1/N from the DFT, while the last term is the hessian normalization for different scales\
    
    Array<Real,1> tempRes1(nx);
    Array<Real,1> tempRes2(ny);
    Array<Real,1> tempRes3(nz);
    for (Int i1=0; i1<nx; ++i1)
    {
        Real temp = 1.;
        if (noDerivativeX==1) temp = mX.kx(i1);
        else if (noDerivativeX==2) temp = mX.kSquare(i1);
        tempRes1(i1) = temp * factor * gauss[i1];
    }
    for (Int i1=0; i1<ny; ++i1)
    {
        Real temp = 1.;
        if (noDerivativeY==1) temp = mY.kx(i1);
        else if (noDerivativeY==2) temp = mY.kSquare(i1);
        tempRes2(i1) = temp * gauss[nx+i1];
    }
    for (Int i1=0; i1<nz; ++i1)
    {
        Real temp = 1.;
        if (noDerivativeZ==1) temp = mZ.kx(i1);
        else if (noDerivativeZ==2) temp = mZ.kSquare(i1);
        tempRes3(i1) = temp * gauss[nx+ny+i1];
    }
    
    
    #pragma omp parallel for
    for (Int i1=0; i1<nx; ++i1)
        for (Int i2=0; i2<ny; ++i2)
            for (Int i3=0; i3<nz; ++i3)
            {
                Int const index = i1*ny*nz + i2*nz + i3;
                result[index] = tempRes1[i1] * tempRes2[i2] * tempRes3[i3] * density[index];
            }
}




/*  Applies Gaussian filter to the input field. */
void gaussianFilter(ArrayReal3D &density,
                    Real const radius,
                    Real const sideLength[],
                    ArrayReal3D *result)
{
    size_t nx = density.getSize(0);
    size_t ny = density.getSize(1);
    size_t nz = density.getSize(2);
    size_t nz2 = density.getSize(2)/2 + 1;
    Real const factor = Real(1.) / Real(nx*ny*nz);    // 1/'totalGrid' is the 1/N from the DFT
    
    
    // compute the FT of the input signal
    ArrayComplex3D denComplex( nx, ny, nz2 );
    computeFFTW( density, denComplex, FFTW_class::R2C, omp_get_max_threads() );
    
    
    // get the values of the FT of the Gaussian function for the given filter radius
    Array<Real,1> gauss( nx+ny+nz );
    gaussianFilterFT( &(gauss[0]), nx, radius, sideLength[0] );          // get the values of the FT of the Gaussian filter
    gaussianFilterFT( &(gauss[nx]), ny, radius, sideLength[1] );
    gaussianFilterFT( &(gauss[nx+ny]), nz, radius, sideLength[2] );
    
    
    // multiply the FT of the input signal and Gaussian filter
    for (Int i1=0; i1<nx; ++i1)
        for (Int i2=0; i2<ny; ++i2)
            for (Int i3=0; i3<nz2; ++i3)
            {
                Int const index = i1*ny*nz2 + i2*nz2 + i3;
                denComplex[index] = factor * gauss[i1] * gauss[nx+i2] * gauss[nx+ny+i3] * denComplex[index];
            }
    
    
    // compute the inverse FT
    computeFFTW( denComplex, *result, FFTW_class::C2R, omp_get_max_threads() );
}



/* Computes the Hessian matrix for a given filtered field. */
void computeHessian(ArrayReal3D &density,
                    ArrayVector9D &eigenvectors,
                    Real const radius,
                    Real const sideLength[],
                    Real const bias)
{
    size_t nx = density.getSize(0);
    size_t ny = density.getSize(1);
    size_t nz = density.getSize(2);
    size_t nz2 = density.getSize(2)/2 + 1;
    Array<Int,1> tempGrid = eigenvectors.getSize<Int>();
    Real const factor = Real(1.) / Real(nx*ny*nz);    // 1/'totalGrid' is the 1/N from the DFT
    cout << "Computing the entries of the Hessian matrix ... " << flush;
    
    
    // compute the FT of the input signal
    ArrayComplex3D denComplex( nx, ny, nz2 );
    computeFFTW( density, denComplex, FFTW_class::R2C, omp_get_max_threads(), false );
    
    Array<Real,1> gauss( nx+ny+nz );
    gauss.assign( Real(1.) );
    
    
    // create variables for the FFTW
    FFTW_class fftw( (eigenvectors.getSize<int>()).ptrData(), NO_DIM, FFTW_class::C2R, omp_get_max_threads() );
    fftw.plan();
    FFTW_COMPLEX *in = fftw.ptrComplexInput();  //points to the input complex data for the FFTW
    
    
    // Now compute each of the 6 independent entries of the Hessian matrix - the entries are saved in the order {h_xx, h_xy, h_xz, h_yy, h_yz, h_zz} in the memory location 0 to 5 starting at eigenvectors[i][0]
    hessianEntry( denComplex.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 2, 0, 0, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][0]), NO_DIM*(NO_DIM-1) );  // writes the FFTW output in memory starting with '&(eigenvectors[0][0])' at every 'NO_DIM2' steps further
    
    hessianEntry( denComplex.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 1, 1, 0, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][1]), NO_DIM*(NO_DIM-1) );
    
    hessianEntry( denComplex.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 1, 0, 1, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][2]), NO_DIM*(NO_DIM-1) );
    
    hessianEntry( denComplex.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 0, 2, 0, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][3]), NO_DIM*(NO_DIM-1) );
    
    hessianEntry( denComplex.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 0, 1, 1, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][4]), NO_DIM*(NO_DIM-1) );
    
    hessianEntry( denComplex.ptrData(), gauss.ptrData(), tempGrid.ptrData(), radius, sideLength, bias, 0, 0, 2, in );
    fftw.FFTW();
    fftw.copyOutput( &(eigenvectors[0][5]), NO_DIM*(NO_DIM-1) );
    
    cout << "Done.\n";
}



/* Computes the hessian matrix, eignevalues and eigenvectors for an already filtered field. */
void hessianMatrix(ArrayReal3D &density,
                   ArrayVector3D &eigenvalues,
                   ArrayVector9D &eigenvectors,
                   Real const radius,
                   Real const sideLength[],
                   Real const bias,
                   bool const computeEigenvectors)
{
    Timer t;
    
    // compute the Hessian and save its values in the array 'eigenvectors'
    t.start();
    computeHessian( density, eigenvectors, radius, sideLength, bias );  //! works only for 3 spatial dimensions
    t.printTime("hessian");

    
    // compute the eigenvalues and eigenvectors
    t.start();
    eigenvaluesEigenvectors( eigenvalues, eigenvectors, computeEigenvectors ); //! works for any dimension of the space
    t.printTime("eigenvalues");
}





/* Computes the eigenvalues and eigenvectors of a 'reduced' tidal field using the filtered density field. */
void reducedTidalFieldEigenvectors(ArrayReal3D &density,
                                   ArrayVector3D &eigenvalues,
                                   ArrayVector9D &eigenvectors,
                                   Real const radius,
                                   Real const sideLength[],
                                   Real const bias,
                                   bool const computeEigenvectors)
{
    Timer t;
    
    // compute the Hessian and save its values in the array 'eigenvectors'
    t.start();
    computeHessian( density, eigenvectors, radius, sideLength, bias );  //! works only for 3 spatial dimensions
    t.printTime("hessian");
    
    
    // take out from the diagonal terms the smoothed density
    for (Int i=0; i<density.size(); ++i)
    {
        eigenvectors[i][0] -= density[i] / Real(3.);
        eigenvectors[i][3] -= density[i] / Real(3.);
        eigenvectors[i][5] -= density[i] / Real(3.);
    }
    
    
    // compute the eigenvalues and eigenvectors
    t.start();
    eigenvaluesEigenvectors( eigenvalues, eigenvectors, computeEigenvectors ); //! works for any dimension of the space
    t.printTime("eigenvalues");
}



