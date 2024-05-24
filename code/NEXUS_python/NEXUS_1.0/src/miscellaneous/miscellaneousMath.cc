#include <iostream>
#include <string>
#include <cmath>
#include <vector>


#include "miscellaneous.h"
using namespace std;



/* This function takes root "rootPower" from an integer and returns also an integer. If the input integer is input = result^rootPower, than it returns result, otherwise it terminates the program.
*/
int rootN(int const input,
          int const rootPower)
{
    lowerBoundCheck( input, 1, "'input' in function 'rootN'" );
    lowerBoundCheck( rootPower, 0, "'rootPower' in function 'rootN'" );
    
    int result;
    double temp = pow( input-1., 1./rootPower ); //since pow returns a double, we take the root from 'input-1' to be sure that we get a double smaller than the integer we hope to find
    result = int( temp ) + 1;	//returns integer part of temp + 1
    
    
    // check that indeed input = result^rootPower
    int temp2 = 1;
    for (int i=0; i<rootPower; ++i)
        temp2 *= result;
    if (temp2==input)
        return result;
    else
    {
        cout << "~~~ ERROR ~~~ The first argument of the function 'rootN' is " << input << " which cannot be written as " << result << "^" << rootPower << "=" << temp <<". The program ended unsuccesfully!\n";
        exit( EXIT_FAILURE );
    }
}

/* This function checks if the root "rootPower" from an integer is an integer (returns true), otherwise returns false.
*/
bool isRootN(int const input,
             int const rootPower)
{
    lowerBoundCheck( input, 1, "'input' in function 'rootN'" );
    lowerBoundCheck( rootPower, 0, "'rootPower' in function 'rootN'" );
    
    int result;
    double temp = pow( input-1., 1./rootPower ); //since pow returns a double, we take the root from 'input-1' to be sure that we get a double smaller than the integer we hope to find
    result = int( temp ) + 1;	//returns integer part of temp + 1
    
    
    // check that indeed input = result^rootPower
    int temp2 = 1;
    for (int i=0; i<rootPower; ++i)
        temp2 *= result;
    if (temp2==input)
        return true;
    return false;
}



/* Computes the magnitude of a 3 dimensional vector.
*/
double vectorMagnitude(int x1,
                       int x2,
                       int x3)
{
    int temp = x1*x1 + x2*x2 + x3*x3;
    return sqrt( double(temp) );
}



/* Computes the magnitude of a vector.
*/
double vectorMagnitude(double x[], int const n)
{
    double temp = 0.;
    for (int i=0; i<n; ++i)
        temp += x[i]*x[i];
    return sqrt( temp );
}



/* Computes the magnitude of a 3 dimensional vector in momentum space using the conventions of the FFTW library.
*/
double momentumMagnitude(int x1,
                       int x2,
                       int x3,
                       int const nGrid)
{
    if (x1>nGrid/2) x1 -= nGrid;
    if (x2>nGrid/2) x2 -= nGrid;
    if (x3>nGrid/2) x3 -= nGrid;
    
    int temp = x1*x1 + x2*x2 + x3*x3;
    return sqrt( double(temp) );
}


/* Computes the magnitude of a 3 dimensional vector in momentum space using the conventions of the FFTW library.
*/
double momentumMagnitude(int x1, int const n1,
                         int x2, int const n2,
                         int x3, int const n3)
{
    if (x1>n1/2)
        x1 -= n1;
    if (x2>n2/2)
        x2 -= n2;
    if (x3>n3/2)
        x3 -= n3;
    
    int temp = x1*x1 + x2*x2 + x3*x3;
    return sqrt( double(temp) );
}



/* This function outputs an array whose entries are at equal logarithmic steps between "minimun" and "maximum".
*/
void log10Steps(double array[],
                double const minimum,
                double const maximum,
                int const size)
{
    double const delta = log10( maximum/minimum ) / (size-1.);
    
    for (int i=0; i<size; ++i)
        array[i] = minimum * pow ( 10., delta*i );
}



/* This function computes the sinc function:
    sin(x)/x
where:
    x = \pi k/n
*/
double sinc(int const k,
            int const n)
{
   if ( k==0 )
       return 1.;
   else
       return sin( PI *double(k)/double(n) ) / ( PI *double(k)/double(n) );
}



/* This function computes the sinc function:
    sin(x)/x
*/
double sinc(double const x,
            double const lowerValue = 1.e-10)
{
    if ( x<lowerValue )
        return 1.;
    else
        return sin(x) / x;
}


/* This function computes the sinc function on a period lattice where k can take values from -n/2+1 to n/2:
    sin(x)/x
where:
    x = \pi k/n
and with k normalized to values from -n/2+1 to n/2.
*/
double sinc_k(int const k,
              int const n)
{
    if ( k==0 )
        return 1.;
    else if (k>n/2)
        return sin( PI *double(k-n)/double(n) ) / ( PI *double(k-n)/double(n) );
    else
        return sin( PI *double(k)/double(n) ) / ( PI *double(k)/double(n) );
}



/* Computes the numerical integral of the function 'f' between the integration boundaries 'xMin' to 'xMax' using the Romberg method for numerical integration. 'maxError' dictates the accuracy of the result, such that the relative error of the numerical integration is smaller than 'maxError'.

This method is suitable for integration from '0.' to 'some finite value, not very large'. If the function goes to zero at infinity, than one can integrate the 2nd part of the function using the function "RombergIntInfinity" given below, which takes the step 'h' proportional to 1/x. In this case 'xMax' will be a very large number.

For details see "Numerical Recipies" or:
http://en.wikipedia.org/wiki/Romberg%27s_method
http://math.fullerton.edu/mathews/n2003/RombergMod.html
*/
double RombergInt(double (*f) (double),
                  double const xMin,
                  double const xMax,
                  double const maxError)
{
    size_t const maxSteps = 20;	//maximum iterations of the method
    double h = xMax - xMin;	//step size
	
    //allocate memory for the triangular matrix that stores the intermediate results of the Romberg method
    vector< vector<double> > res;	//vector for the triangular matrix
    res.resize(maxSteps);
    for (size_t i=0; i<res.size(); ++i)
        res.at(i).resize(i+1);
	
	//now the actual calculation
    size_t c1 = 1;	//a variable that stores 2^(i-1) in the next loop, needed to compute res[i][0]
    size_t c2 = 1;	//a variable that stores 4^j, needed when computing res[i][j] with j>0
	// the 1st element is just:
    res.at(0).at(0) = .5 * h * ( (*f)(xMin) + (*f)(xMax) );
	// the next elements:
    for (size_t i=1; i<res.size(); ++i)
    {
        h *= .5;
        double x = xMin - h;
		
		// to compute the first column must use the trapeizodal rule
        res.at(i).at(0) = .5 * res.at(i-1).at(0);
        for (size_t j=1; j<=c1; ++j)
        {
            x += 2.*h;	//now we have "x = xMin + (2*j-1) * h"
            res.at(i).at(0) += h * (*f)(x);
        }
		
		// compute res[i][j] with j>0
        for (size_t j=1; j<res.at(i).size(); ++j)
        {
            c2 *= 4;	//now we have "c2 = 4^j"
            res.at(i).at(j) = res.at(i).at(j-1) + 1./(c2-1.) * (res.at(i).at(j-1) - res.at(i-1).at(j-1));
        }
		
		// check if convergence is fulfilled
        if ( abs( res.at(i).at(i)-res.at(i-1).at(i-1) )<abs( maxError*res.at(i).at(i) ) )
            return res.at(i).at(i);
		
        c1 *= 2;	//now "c1 = 2^(i-1)"
        c2 = 1;		//initialize c2=1 to prepare for next loop
    }
	
	// if the method wasn't convergent fast enough, print error message and stop
    double relativeError = abs( ( res.at(maxSteps-1).at(maxSteps-1) - res.at(maxSteps-2).at(maxSteps-2) ) / res.at(maxSteps-1).at(maxSteps-1) );
    cout.precision(4);
    cout << scientific;
	
    cout << "~~~ ERROR ~~~ The Romberg integration method failled to converge in " << maxSteps << " steps. Using " << maxSteps << " steps, the Romberg method yields 'integral'=" << res.at(maxSteps-1).at(maxSteps-1) << " with the relative change commpared to the result at the previos step being 'relative change'=" << relativeError << ". The condition for the result to be deemed acceptable is 'relative error'<" << maxError << ". The program ended unsuccesfully! \n";
    exit( EXIT_FAILURE );
}



/* Computes the numerical integral of the function 'f' between the integration boundaries 'xMin' to 'xMax' using the Romberg method for numerical integration. 'maxError' dictates the accuracy of the result, such that the relative error of the numerical integration is smaller than 'maxError'.

This is the same as the function "RombergInt" given above, just that the x-steps taken are proportional to 1/x. This method is suitable for integration from 'some nozero value' to 'infinity'. In this case 'xMax' will be a very large number.

*/
double RombergIntInfinity(double (*f) (double),
                          double const xMin,
                          double const xMax,
                          double const maxError)
{
    size_t const maxSteps = 20;	//maximum iterations of the method
    double const xmin = 1./xMax;    //since we integrate using the 1/x variable, this are the new integration limits
    double const xmax = 1./xMin;
    double h = xmax - xmin;	//step size
	
    //allocate memory for the triangular matrix that stores the intermediate results of the Romberg method
    vector< vector<double> > res;	//vector for the triangular matrix
    res.resize(maxSteps);
    for (size_t i=0; i<res.size(); ++i)
        res.at(i).resize(i+1);
	
	//now the actual calculation
    size_t c1 = 1;	//a variable that stores 2^(i-1) in the next loop, needed to compute res[i][0]
    size_t c2 = 1;	//a variable that stores 4^j, needed when computing res[i][j] with j>0
	// the 1st element is just:
    res.at(0).at(0) = .5 * h * ( (*f)(1./xmin)/(xmin*xmin) + (*f)(1./xmax)/(xmax*xmax) );
	// the next elements:
    for (size_t i=1; i<res.size(); ++i)
    {
        h *= .5;
        double x = xmin - h;
	
		// to compute the first column must use the trapeizodal rule
        res.at(i).at(0) = .5 * res.at(i-1).at(0);
        for (size_t j=1; j<=c1; ++j)
        {
            x += 2.*h;	//now we have "x = xmin + (2*j-1) * h"
            res.at(i).at(0) += h * (*f)(1./x)/(x*x);
        }
	
		// compute res[i][j] with j>0
        for (size_t j=1; j<res.at(i).size(); ++j)
        {
            c2 *= 4;	//now we have "c2 = 4^j"
            res.at(i).at(j) = res.at(i).at(j-1) + 1./(c2-1.) * (res.at(i).at(j-1) - res.at(i-1).at(j-1));
        }
	
		// check if convergence is fulfilled
        if ( abs( res.at(i).at(i)-res.at(i-1).at(i-1) )<abs( maxError*res.at(i).at(i) ) )
            return res.at(i).at(i);
		
        c1 *= 2;	//now "c1 = 2^(i-1)"
        c2 = 1;		//initialize c2=1 to prepare for next loop
    }
	
	// if the method wasn't convergent fast enough, print error message and stop
    double relativeError = abs( ( res.at(maxSteps-1).at(maxSteps-1) - res.at(maxSteps-2).at(maxSteps-2) ) / res.at(maxSteps-1).at(maxSteps-1) );
    cout.precision(4);
    cout << scientific;
	
    cout << "~~~ ERROR ~~~ The Romberg integration 1/x method failled to converge in " << maxSteps << " steps. Using " << maxSteps << " steps, the Romberg method yields 'integral'=" << res.at(maxSteps-1).at(maxSteps-1) << " with the relative change commpared to the result at the previos step being 'relative change'=" << relativeError << ". The condition for the result to be deemed acceptable is 'relative error'<" << maxError << ". The program ended unsuccesfully! \n";
    exit( EXIT_FAILURE );
}








