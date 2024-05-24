#include <iostream>
#include <cmath>
#include <vector>

#include "MMF_computations.h"
#include <miscellaneous.h>
#include <array.h>
using namespace std;



/* Rescales the eigenvalue according to the bias value. */
void rescaleEigenvalues(ArrayVector3D *eigenvalues,
                        Real const radius,
                        Real const bias)
{
    Real factor = pow( radius, Real(2.*bias) );
    Int const gridSize = eigenvalues->totalSize();
    cout << "Rescaling Hessian matrix eigenvalues using 'bias'=" << bias << " (multiplication by " << factor << ") ... " << flush;
    
    for (Int i=0; i<gridSize; ++i)
        (*eigenvalues)[i] *= factor;
    
    cout << "Done.\n";
}


/* Selects the maximum response eigenvalues and eigenvectors for each grid cell. */
void eigenMaximumResponse(Array<int,3> const &scale,
                          ArrayVector3D *eigenvalues,
                          ArrayVector9D *eigenvectors,
                          ArrayVector3D &eigenvalues2,
                          ArrayVector9D &eigenvectors2,
                          int const scale2)
{
    cout << "Extracting maximum response eigenvalues and eigenvectors from scale '" << scale2 << "' ... " << flush;
    Int const gridSize = scale.totalSize();
    
    for (Int i=0; i<gridSize; ++i)
        if ( scale.at(i)==scale2 )
        {
            (*eigenvalues)[i] = eigenvalues2[i];
            (*eigenvectors)[i] = eigenvectors2[i];
        }
    cout << "Done.\n";
}



/* This function selects the eigenvector corresponding to the given feature:
        feature = 3 - filament - : selects the 3rd eigenvector (in 3D) that gives the filament direction
        feature = 2 - wall - : selects the 1st eigenvector that gives the wall normal
*/
void featureEigenvector(ArrayVector9D &eigenvectors,
                        ArrayVector3D *eigenvalues,
                        int const feature)
{
    cout << "Computing the feature eigenvector for " << (feature==3? "filaments" : "walls") << " ... " << flush;
    if ( feature==3 )
    {
        for (size_t i=0; i<eigenvectors.size(); ++i)
        {
            (*eigenvalues)[i] = lastEigenvector( eigenvectors[i] );
            if ( (*eigenvalues)[i][0]<Real(0.) ) (*eigenvalues)[i] *= Real(-1.);
        }
    }
    else if ( feature==2 )
    {
        for (size_t i=0; i<eigenvectors.size(); ++i)
        {
            for (size_t j=0; j<NO_DIM; ++j)
                (*eigenvalues)[i][j] = eigenvectors[i][j];
            if ( (*eigenvalues)[i][0]<Real(0.) ) (*eigenvalues)[i] *= Real(-1.);
        }
    }
    else throwError( "Unknown value for the 'feature' argument in function 'featureEigenvector'." );
    cout << "Done.\n";
//     for (int i=0; i<10; ++i)
//     {
//         cout << "( ";
//         for (int j=0; j<NO_DIM; ++j)
//             cout << (*eigenvalues)[i][j] << ",";
//         cout << ")\n";
//         Real temp1 = eigenvectors[i][0]*(*eigenvalues)[i][0] + eigenvectors[i][1]*(*eigenvalues)[i][1] + eigenvectors[i][2]*(*eigenvalues)[i][2];
//         Real temp2 = eigenvectors[i][3]*(*eigenvalues)[i][0] + eigenvectors[i][4]*(*eigenvalues)[i][1] + eigenvectors[i][5]*(*eigenvalues)[i][2];
//         cout << temp1 << "  " << temp2 << "\n";
//     }
}



/* This function computes the N-th eigenvector given N-1 eigenvectors. Works for 2D and 3D cases. */
VectorReal3D lastEigenvector(VectorReal9D &e)
{
    VectorReal3D temp;
#if NO_DIM==2
    temp[0] = -e[1];
    temp[1] = e[0];
#elif NO_DIM==3
    temp[0] = e[2]*e[4] - e[1]*e[5];
    temp[1] = e[0]*e[5] - e[2]*e[3];
    temp[2] = e[1]*e[3] - e[0]*e[4];
    
    temp.normalize();   // normalize the vector to magnitude 1
#endif
    return temp;
}


