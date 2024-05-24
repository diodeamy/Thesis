#include <string>
#include <cmath>


#include "MMF_computations.h"
#include <miscellaneous.h>
using namespace std;

// functions for node detection
void nodeFilter40(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);
void nodeFilter41(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);
void nodeFilter42(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);
void nodeFilter43(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);
void nodeFilter44(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);

// functions for filament detection
void filamentFilter30(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response);
void filamentFilter31(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response);
void filamentFilter32(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response);
void filamentFilter33(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response);
void filamentFilter34(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response);

// functions for wall detection
void wallFilter20(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);
void wallFilter21(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);
void wallFilter22(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);
void wallFilter23(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);
void wallFilter24(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response);





/* Computes the response function using the selected filter. */
void computeResponse(ArrayVector3D &eigenvalue,
                     ArrayReal3D &response,
                     int const filterNumber)
{
    switch(filterNumber)
    {
        case 40:
            nodeFilter40( eigenvalue, response );
            return;
        case 41:
            nodeFilter41( eigenvalue, response );
            return;
        case 42:
            nodeFilter42( eigenvalue, response );
            return;
        case 43:
            nodeFilter43( eigenvalue, response );
            return;
        case 44:
            nodeFilter44( eigenvalue, response );
            return;
        
        case 30:
            filamentFilter30( eigenvalue, response );
            return;
        case 31:
            filamentFilter31( eigenvalue, response );
            return;
        case 32:
            filamentFilter32( eigenvalue, response );
            return;
        case 33:
            filamentFilter33( eigenvalue, response );
            return;
        case 34:
            filamentFilter34( eigenvalue, response );
            return;
        
        case 20:
            wallFilter20( eigenvalue, response );
            return;
        case 21:
            wallFilter21( eigenvalue, response );
            return;
        case 22:
            wallFilter22( eigenvalue, response );
            return;
        case 23:
            wallFilter23( eigenvalue, response );
            return;
        case 24:
            wallFilter24( eigenvalue, response );
            return;
        default:
            throwError( "Unknown filter type in the function 'computeResponse'." );
    }
}


/* Outputs the name of the coresponding filter and a short description of it. */
string filterDescription(int const filterNumber)
{
    switch(filterNumber)
    {
        case 40:
             return "Filter 40 -- node filter: described in the algorithm paper.";
        case 41:
             return "Filter 41 -- none.";
        case 42:
             return "Filter 42 -- none.";
        case 43:
             return "Filter 43 -- none.";
        case 44:
             return "Filter 44 -- node filter: no strength feature, using only the eigenvalue lambda_3 as characteristic of the environment.";
        
        case 30:
             return "Filter 30 -- filament filter: described in the algorithm paper.";
        case 31:
             return "Filter 31 -- filament filter: described in the algorithm paper - but NOT the requirement that |lamda_2|>|lambda_3|.";
        case 32:
             return "Filter 32 -- filament filter: described in the algorithm paper - but NOT the requirement that |lamda_2|>|lambda_3| and also uses (1-lambda_3/lambda_2) instead of (1-abs(lambda_3/lambda_2)).";
        case 33:
             return "Filter 33 -- filament filter: no strength feature, using logarithm(eigenvalue lambda_2) as characteristic of the environment.";
        case 34:
             return "Filter 34 -- filament filter: no strength feature, using only the eigenvalue lambda_2 as characteristic of the environment.";
        
        case 20:
             return "Filter 20 -- wall filter: described in the algorithm paper.";
        case 21:
             return "Filter 21 -- wall filter: described in the algorithm paper - but NOT the requirements that |lamda_1|>|lambda_3| and |lamda_1|>|lambda_2|.";
        case 22:
             return "Filter 22 -- wall filter: described in the algorithm paper - but NOT the requirements that |lamda_2|>|lambda_3|and |lamda_1|>|lambda_2|. Also uses (1-lambda_(2,3)/lambda_1) instead of (1-abs(lambda_(2,3)/lambda_1))";
        case 23:
             return "Filter 23 -- wall filter: no strength feature, using logarithm(eigenvalue lambda_1) as characteristic of the environment.";
        case 24:
             return "Filter 24 -- wall filter: no strength feature, using only the eigenvalue lambda_1 as characteristic of the environment.";
        
        default:
            return "Unknown filter type.";
    }
}





/********************************************************************
The filters are implemented for the eigenvalues of the Hessian
matrix of the density. In this case the eigenvalues are negative and
they are ordered according to:  l_1 < l_2 < l_3
(with l_1, l_2 and l_3 the eigenvelues of the Hessian matrix)
*********************************************************************/

//! Filter functions for nodes.
/* The node filter function 40 - the one described in the algorithm paper. */
void nodeFilter40(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    Real const beta = Real( 0.5 );
    Real const beta2 = 2. * beta * beta;
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'nodeFilter40' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the node response function using filter 40 ... " << flush;
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][2]>=0 ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        Real temp = eigenvalue[i][2] / eigenvalue[i][0];
        temp = temp*temp / beta2;
        response[i] = (1. - exp(-temp)) * fabs( eigenvalue[i][2] );
    }
    cout << "Done.\n";
}
/* The node filter 41 - None implemented. */
void nodeFilter41(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'nodeFilter41' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    
    throwError( "No node filter is implemented in function 'nodeFilter41'. There is no filter implemented for filter number '41'." );
    
    cout << "Computing the node response function using filter 41 ... " << flush;
    cout << "Done.\n";
}
/* The node filter 42 - none implemented. */
void nodeFilter42(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'nodeFilter42' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    
    throwError( "No node filter is implemented in function 'nodeFilter42'. There is no filter implemented for filter number '42'." );
    
    cout << "Computing the node response function using filter 42 ... " << flush;
    cout << "Done.\n";
}
/* The node filter 43 -none implemented . */
void nodeFilter43(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'nodeFilter43' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    
    throwError( "No node filter is implemented in function 'nodeFilter43'. There is no filter implemented for filter number '43'." );
    
    cout << "Computing the node response function using filter 43 ... " << flush;
    cout << "Done.\n";
}
/* The node filter 44 - uses only |lambda_3| as response for the node environment. */
void nodeFilter44(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'nodeFilter44' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the node response function using filter 44 ... " << flush;
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][2]>=0 ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        response[i] = fabs( eigenvalue[i][2] );
    }
    cout << "Done.\n";
}



//! Filter functions for filaments.
/* The filament filter function 30 - the one described in the algorithm paper. */
void filamentFilter30(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'filamentFilter30' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the filament response function using filter 30 ... " << flush;
    
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][1]>=0 or fabs(eigenvalue[i][1])<fabs(eigenvalue[i][2]) ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        Real temp = (eigenvalue[i][1]/eigenvalue[i][0]) * (1. - fabs(eigenvalue[i][2]/eigenvalue[i][1]) );
        response[i] =  temp * fabs(eigenvalue[i][1]);
    }
    cout << "Done.\n";
}
/* The filament filter function 31 - the one described in the algorithm paper - but NOT the requirement that |lamda_2|>|lambda_3|. */
void filamentFilter31(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    Real const beta = Real( 0.5 );
    Real const beta2 = 2. * beta * beta;
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'filamentFilter31' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the filament response function using filter 31 ... " << flush;
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][1]>=0 or fabs(eigenvalue[i][1])<fabs(eigenvalue[i][2]) ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        Real temp = (eigenvalue[i][1]/eigenvalue[i][0]) * (1. - fabs(eigenvalue[i][2]/eigenvalue[i][1]) );
        temp = temp*temp / beta2;
        response[i] = (1. - exp(-temp)) * fabs( eigenvalue[i][1] );
    }
    cout << "Done.\n";
}
/* The filament filter function 32 - the one described in the algorithm paper - but NOT the requirement that |lamda_2|>|lambda_3| and also uses (1-lambda_3/lambda_2) instead of (1-abs(lambda_3/lambda_2)). */
void filamentFilter32(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    Real const beta = Real( 0.5 );
    Real const beta2 = 2. * beta * beta;
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'filamentFilter32' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the filament response function using filter 32 ... " << flush;
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][1]>=0 ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        Real temp = (eigenvalue[i][1]/eigenvalue[i][0]) * (1. - fabs(eigenvalue[i][2]/eigenvalue[i][1]) );
        temp = temp*temp / beta2;
        response[i] = (1. - exp(-temp)) * fabs( eigenvalue[i][1] );
    }
    cout << "Done.\n";
}
/* The filament filter function 33 - no strength feature, using logarithm(eigenvalue lambda_2) as characteristic of the environment. */
void filamentFilter33(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'filamentFilter33' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the filament response function using filter 33 ... " << flush;
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][1]>=0 ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        response[i] = log10( 1.e10 * fabs( eigenvalue[i][1] ) );
    }
    cout << "Done.\n";
}
/* The filament filter function 34 - no strength feature, uses eigenvalue lambda_2 as characteristic of the environment. */
void filamentFilter34(ArrayVector3D &eigenvalue,
                      ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'filamentFilter34' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the filament response function using filter 34 ... " << flush;
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][1]>=0 ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        response[i] = fabs( eigenvalue[i][1] );
    }
    cout << "Done.\n";
}




//! Filter functions for walls.
/* The wall filter function 20 - described in the algorithm paper. */
void wallFilter20(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'wallFilter20' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the wall response function using filter 20 ... " << flush;
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][0]>=0 or fabs(eigenvalue[i][0])<fabs(eigenvalue[i][2]) or fabs(eigenvalue[i][0])<fabs(eigenvalue[i][1]) ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        Real temp = ( fabs(eigenvalue[i][0]) - fabs(eigenvalue[i][1]) ) * (fabs(eigenvalue[i][0])- fabs(eigenvalue[i][2]) );
        if ( fabs( eigenvalue[i][0] )>1.e-5 )
            response[i] = temp / fabs( eigenvalue[i][0] );
        else
            response[i] = temp;
    }
    cout << "Done.\n";
}
/* The wall filter function 21 - described in the algorithm paper - but NOT the requirements that |lamda_1|>|lambda_3| and |lamda_1|>|lambda_2|. */
void wallFilter21(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    Real const beta = Real( 0.5 );
    Real const beta2 = 2. * beta * beta;
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'wallFilter21' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the wall response function using filter 21 ... " << flush;
    
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][0]>=0 or fabs(eigenvalue[i][0])<fabs(eigenvalue[i][2]) or fabs(eigenvalue[i][0])<fabs(eigenvalue[i][1]) ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        Real temp = (1.- fabs(eigenvalue[i][1]/eigenvalue[i][0]) ) * (1.- fabs(eigenvalue[i][2]/eigenvalue[i][0]) );
        temp = temp*temp / beta2;
        response[i] = (1. - exp(-temp)) * fabs( eigenvalue[i][0] );
    }
    cout << "Done.\n";
}
/* The wall filter function 22 - described in the algorithm paper - but NOT the requirements that |lamda_1|>|lambda_3| and |lamda_1|>|lambda_2|. Also uses (1-lambda_(2,3)/lambda_1) instead of (1-abs(lambda_(2,3)/lambda_1)).*/
void wallFilter22(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    Real const beta = Real( 0.5 );
    Real const beta2 = 2. * beta * beta;
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'wallFilter22' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the wall response function using filter 22 ... " << flush;
    
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][0]>=0 ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        Real temp = (1.- fabs(eigenvalue[i][1]/eigenvalue[i][0]) ) * (1.- fabs(eigenvalue[i][2]/eigenvalue[i][0]) );
        temp = temp*temp / beta2;
        response[i] = (1. - exp(-temp)) * fabs( eigenvalue[i][0] );
    }
    cout << "Done.\n";
}
/* The wall filter function 23 - no strength feature, using logarithm(eigenvalue lambda_1) as characteristic of the environment. */
void wallFilter23(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'wallFilter23' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the wall response function using filter 23 ... " << flush;
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][0]>=0. ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        response[i] = log10( 1.e10 * fabs( eigenvalue[i][0] ) );
    }
    cout << "Done.\n";
}

/* The wall filter function 23 - no strength feature, using eigenvalue lambda_1 as characteristic of the environment. */
void wallFilter24(ArrayVector3D &eigenvalue,
                  ArrayReal3D &response)
{
    Int const totalSize = response.totalSize();
    if ( totalSize!=eigenvalue.totalSize() ) throwError( "In function 'wallFilter24' the 'eigenvalue' and 'response' arrays have different sizes. Both must have the same dimensions." );
    cout << "Computing the wall response function using filter 24 ... " << flush;
    
    for (Int i=0; i<totalSize; ++i)
    {
        if ( eigenvalue[i][0]>=0 ) // apply morphology mask
        {
            response[i] = 0.;
            continue;
        }
        
        // morphology response
        response[i] = fabs( eigenvalue[i][0] );
    }
    cout << "Done.\n";
}
