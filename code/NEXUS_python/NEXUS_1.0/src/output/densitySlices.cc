#include <iostream>

#include "output.h"
#include <miscellaneous.h>
using namespace std;




/* This function computes a slice trough the 3D density map.
    axis = 0,1 or 2 gives the direction along which to take the slice
    posMin, posMax = the min and max position of the slice in Mpc
*/
void densitySliceComputation(Real *delta3D,
                             Real *d2D,
                             Int const nGrid[3],
                             int const axis,
                             Real posMin,
                             Real posMax,
                             Real const boxLength)
{
    if ( posMin<0. ) posMin += boxLength;
    if ( posMax<0. ) posMax += boxLength;
    if ( posMin>=boxLength ) posMin -= boxLength;
    if ( posMax>=boxLength ) posMax -= boxLength;
    cout << "Computing a density slice along direction " << axis << " between " << posMin << " to " << posMax << " Mpc ...  " << flush;
    
    // the next part computes the position of the bins that have to be summed for the density slice along the "axis" direction
    Real const densityCell = boxLength / nGrid[axis];//thickness of 3D density cell in Mpc
    Int const gridMin = int( posMin/densityCell );    // lowest 3D map bin along the "axis" direction
    Int const gridMax = int( posMax/densityCell );    // highest 3D map bin along the "axis" direction
    intervalCheck( gridMin, Int(0), nGrid[axis]-1, "'gridMin' in function 'outputDensitySlice'"); //checks that indeed 0<= gridMin <nGrid
    intervalCheck( gridMax, Int(0), nGrid[axis]-1, "'gridMax' in function 'outputDensitySlice'"); //checks that indeed 0<= gridMin <nGrid
    // the outermost bins may not be included fully, hence define a fraction of how much of it to include
    Real const fractionMin = gridMin+1. - posMin/densityCell;    // fraction for lowest bin
    Real const fractionMax = posMax/densityCell - gridMax;    // fraction for highest bin
    Real normalization = gridMax - gridMin -1. + fractionMin + fractionMax;   // normalization constant such that the density is expressed in units of average density
    if ( gridMax==gridMin ) normalization = fractionMin + fractionMax;  //if the slice is very thin and is all inside a density map layer
    
    
    // do the calculations; first add the outermost bins after which the bins in between
    int const j1 = perpendicularDirection( axis, 1 );
    int const j2 = perpendicularDirection( axis, 2 );
    for (Int i1=0; i1<nGrid[j1]; ++i1)
    {
        for (Int i2=0; i2<nGrid[j2]; ++i2)
        {
            Int index1 = i1*nGrid[j2] + i2; //index for the 2D density slice
            
            //add the lowest most bin
            Int index2 = compute_2D_index( i1, i2, gridMin, axis, nGrid );  //index for the 3D density map
            d2D[index1] += delta3D[index2] * fractionMin;
            //add the bins in between the lowest and highest
            for (Int i3=gridMin+1; i3<gridMax; ++i3)
            {
                index2 = compute_2D_index( i1, i2, i3, axis, nGrid );
                d2D[index1] += delta3D[index2];
            }
            //add the highest most bin
            index2 = compute_2D_index( i1, i2, gridMax, axis, nGrid );
            d2D[index1] += delta3D[index2] * fractionMax;
            
            //now we still need to divide by the number of bins summed above such that we express the result with respect to the average density of the universe
            d2D[index1] /= normalization;
        }
    }
    cout << "Done.\n";
}


/* Computes the index of the 3D density map grid cells that have to be merged to obtain the 2D density map. 'n1' and 'n2' give the grid cell of the 2D density map while while n3 gives the position along the 3rd axis that has to be stacked to get the 2D map. 'choice' tells along which direction we make the slices (which direction we stack up).
*/
Int compute_2D_index(Int const n1,
                     Int const n2,
                     Int const n3,
                     int const choice,
                     Int const nGrid[])
{
    Int n3Temp = n3;
    if ( n3Temp<0 ) n3Temp += nGrid[choice];
    else if ( n3Temp>=nGrid[choice] ) n3Temp -= nGrid[choice];
    
    switch (choice)
    {
        case 0:
            return n3Temp*nGrid[1]*nGrid[2] + n1*nGrid[2] + n2;
            break;
        case 1:
            return n1*nGrid[1]*nGrid[2] + n3Temp*nGrid[2] + n2;
            break;
        case 2:
            return n1*nGrid[1]*nGrid[2] + n2*nGrid[2] + n3Temp;
            break;
        default:
            throwError( "Invalid input for the 4th argument in function 'compute_2D_index'." );
    }
    return -1;
}



/* Returns what are the directions that span the slices. 'direction' gives the direction perpendicular to the slice and 'choice' tells the function to return the 1st or 2nd direction that gives the plane of the slice. For example if the slice is along the z (direction=2) direction, than it returns the x (==0) and y (==1) directions.
*/
int perpendicularDirection(int const direction,
                           int const choice)
{
    if ( choice!=1 and choice!=2 )
        throwError( "Invalid input for the 2nd argument in function 'perpendicularDirection'. This argument can take only the values 1 or 2." );
    
    switch( direction )
    {
        case 0:
            if (choice==1) return 1;
            else if (choice==2) return 2;
            break;
        case 1:
            if (choice==1) return 0;
            else if (choice==2) return 2;
            break;
        case 2:
            if (choice==1) return 0;
            else if (choice==2) return 1;
            break;
        default:
            throwError( "Invalid input for the 1st argument in function 'perpendicularDirection'.");
    }
    return -1;
}


