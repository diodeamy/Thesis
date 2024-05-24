
#include <string>
#include <fstream>
#include <cmath>
#include <limits>

#include <defines.h>
#include <array.h>
#include <threshold.h>
#include <miscellaneous.h>
#include <linear_fit.h>
#include <arrayProperties.h>
using namespace std;

namespace LSM = LEAST_SQUARE_METHOD;
typedef LSM::FitPowerLawFunction<Real>                       powerLaw;
typedef LSM::FitLinearFunction<Real>                         linear;
typedef LSM::FitParabolaFunction<Real>                       parabola;
typedef LSM::FitLogNormalFunction<Real>                      logNormal;
typedef LSM::FitMultipleFunctions<Real,powerLaw,powerLaw>    twoPowerLaw;



/* Finds the threshold x-value where the y value == 'optimalValue'. */
bool optimalThreshold(vector<Real> &x,
                      vector<Real> &y,
                      Real optimalValue,
                      bool const highToLow,
                      Real *optimalX)
{
    Real minI = -1, maxI = -1;
    if ( highToLow )    //start the search with the larger x-values and descrease x
    {
        for (size_t i=0; i<x.size()-1; ++i)
            if ( optimalValue>=y[i+1] and optimalValue<y[i] )
            {
                minI = i+1;
                maxI = i;
                break;
            }
    }
    else
    {
        for (size_t i=x.size()-1; i>0; --i)
            if ( optimalValue>=y[i-1] and optimalValue<y[i] )
            {
                minI = i-1;
                maxI = i;
                break;
            }
    }
    if (minI==-1) return false; //optimal value could not be found
    
    
    // interpolate the threshold value between the two values found
    Real slope = (y[maxI]-y[minI]) / (x[maxI]-x[minI]);
    Real constant = y[minI] - slope*x[minI];
    *optimalX = (optimalValue - constant) / slope;
    return true;
}

/* Fits a parabola (2nd degree polynomial) for y=f(x) and returns the maximum of the function. */
bool parabolaFitMaximum(vector<Real> &x,
                        vector<Real> &y,
                        Real fitIntervalFraction,
                        Real *optimalX,
                        Real fitParams[])
{
    // find all the points that have y-values between yMax*fitIntervalFraction < y < yMax
    intervalCheck<Real>( fitIntervalFraction, 0., 1., "'fitIntervalFraction' argument in function 'parabolaFitMaximum'" );
    Real yMax = Max<Real>( &(y[0]), y.size() );
    int minI = -1, maxI = -1;
    for (size_t i=1; i<x.size(); ++i)   //first find 'minI'
        if ( y[i]>yMax*fitIntervalFraction )
        {
            minI = i-1;
            break;
        }
    for (size_t i=x.size()-2; i>0; --i)   //next find 'maxI'
        if ( y[i]>yMax*fitIntervalFraction )
        {
            maxI = i+1;
            break;
        }
    // check that there are enough data points to have a decent fit
    int minNoPoints = 7;
    minNoPoints = y.size()<size_t(minNoPoints) ? y.size() : minNoPoints;
    while ( (maxI-minI)<minNoPoints )   //if not enough points, add more at the ends of the interval
    {
        if ( minI>0 ) --minI;
        if ( maxI<int(y.size()-1) ) ++maxI;
    }
    
    
    // fit a 2nd order polynomial to the data
    // make the arrays used for the fit
    int const noFitElem = maxI-minI;
    Real xData[noFitElem], yData[noFitElem];
    for (int i=minI; i<=maxI; ++i)
    {
        xData[i-minI] = x[i];
        yData[i-minI] = y[i];
    }
    parabola fit;
    fit.dataToFit( xData, yData, noFitElem );
    size_t const noParameters = fit.noFitParamaters();
    for (size_t i=0; i<noParameters; ++i)
        fitParams[i] = fit.param_n( i );
    
    // get the maximum
    *optimalX = -fitParams[1] / (2.*fitParams[2]);    // the maximum of the fit function
    if ( *optimalX<Min(x[minI],x[maxI]) or *optimalX>Max(x[minI],x[maxI]) )
        return false;   //the maximum is outside the fit interval - something went wrong
    return true;
}


/* Computes the derivatives of the mass and volume increase as a function of threshold, derivative with respect to log10(threshold). */
void thresholdDerivatives(Array<Real,2> &data,
                          int const i[][2],
                          int const noVariables)
{
    // derivative for first point
    int N = data.getSize(0);    // the number of data rows in 'data'
    Real dx = log10( data(1,0) / data(0,0) ); //the size of the threshold variation - user for derivative normalization
    for (int j=0; j<noVariables; ++j)
        data(0,i[j][0]) = ( data(0,i[j][1]) - data(1,i[j][1]) ) / dx;  //derivative
    
    // derivatives for middle points
    for (int k=1; k<N-1; ++k)
    {
        dx = log10( data(k+1,0) / data(k-1,0) ); //the size of the threshold variation - user for derivative normalization
        for (int j=0; j<noVariables; ++j)
            data(k,i[j][0]) = ( data(k-1,i[j][1]) - data(k+1,i[j][1]) ) / dx;  //derivative
    }
    
    //derivative for last point
    dx = log10( data(N-1,0) / data(N-2,0) ); //the size of the threshold variation - user for derivative normalization
    for (int j=0; j<noVariables; ++j)
        data(N-1,i[j][0]) = ( data(N-2,i[j][1]) - data(N-1,i[j][1]) ) / dx;  //derivative
}


void writeData(Array<Real,2> &data,
               fstream &outputFile)
{
    for (size_t i=0; i<data.getSize(0); ++i)
    {
        outputFile << data(i,0);
        for (size_t j=1; j<data.getSize(1); ++j)
            outputFile << "\t" << data(i,j);
        outputFile << "\n";
    }
}


/* Runs over the data and uses a simple filter to damp out any large changes in the input data. */
void smoothData(vector<Real> &y,
                Real maxRatio)
{
    size_t noElem = y.size();
    
    if (y[0]>(1.+maxRatio)*y[1])
        y[0] = y[1];
    
    for (size_t i=1; i<noElem-1; ++i)
    {
        Real temp1 = fabs( y[i-1] - y[i] );
        Real temp2 = fabs( y[i+1] - y[i] );
        temp1 = temp1>temp2 ? temp2 : temp1;
        Real avg = (y[i-1] + y[i+1]) / 2.;
        if ( temp1>maxRatio*avg )
            y[i] = avg;
    }
    
    if (y[noElem-1]>(1.+maxRatio)*y[noElem-2])
        y[noElem-1] = y[noElem-2];
}



/* Computes the optimal node threshold. We define significant nodes the ones with mass higher than 'minimumSize'.
The optimal node treshold is defined as the response value where the fraction of significant nodes with average density larger than the virial density is 0.5.
It prints the data and the fit details to an ASCII file.
*/
Real outputNodeThreshold(Array<Real,3> &density,
                         Array<Real,3> &response,
                         int const neighborType,
                         double const virialDensity,
                         double const minimumSize,
                         Array<Real,1> &threshold,
                         string outputFilename,
                         string const programOptions)
{
    // Find the fraction of virial objects for a given threshold value
    size_t const noThresholds = threshold.size();
    cout << "\n\nComputing the optimal threshold for nodes using " << noThresholds << " thresholds values, with 'virial density' = " << virialDensity << " and 'minSize' = " << minimumSize << " in units of average mass in voxel:\n" << flush;
    cout << "  threshold index    threshold value     valid nodes(\%)\n";
    
    
    // assign memory for the quantities needed by the computation
    Array<Real,2> thresholdResponse( noThresholds, 11 );    // keeps track of the final results
    Array<int,3>  mask( density.getSize(0), density.getSize(1), density.getSize(2) );         //mask that keeps track of the objects
    mask.assign( -1 );  //initialize the mask to empty cells
    int trackStatistics[] = {0,0};  // array needed by the 'fractionVirialNodes' function
    
    
    // loop over the thresholds and get the data
    Real norm = Real( density.size() ); //normalization factor for mass and volume data - this gives the total mass and volume in the box
    Real minThreshold = threshold[0], maxThreshold = numeric_limits<Real>::max();
    for (size_t i=0; i<noThresholds; ++i)
    {
        if (i%10==0 or i==noThresholds-1)
            cout << "   " << setw(4) << i+1 << " / " << noThresholds << "\t\t" << setprecision(2) << setw(4) << minThreshold << "\t\t   " << flush;
        Real properties[6];
        fractionVirialNodes( response, minThreshold, maxThreshold, &mask, neighborType, density, minimumSize, virialDensity, properties, trackStatistics );
        
        thresholdResponse(i,0) = minThreshold;       //threshold value
        thresholdResponse(i,1) = properties[2];      //fraction valid nodes (mass>minimumSize and averageDensity > virialDensity)
        thresholdResponse(i,3) = properties[0]/norm; //total mass in significant nodes (mass>minimumSize)
        thresholdResponse(i,5) = properties[1]/norm; //total volume in significant nodes
        thresholdResponse(i,6) = properties[5];      //number of valid objects
        thresholdResponse(i,8) = properties[3]/norm; //total mass in all nodes
        thresholdResponse(i,10) = properties[4]/norm;//total volume in all nodes
        
        maxThreshold = threshold[i];
        if (i<noThresholds-1) minThreshold = threshold[i+1];
        if (i%10==0 or i==noThresholds-1)
            cout << setprecision(3) << setw(4) << thresholdResponse(i,1)*100. << "\n";
    }
    
    
    // get the mass and volume change between two different thresholds
    int indices[][2] = { {2,3}, {4,5}, {7,8}, {9,10} };  //each pair of indices gives the values (index quantity variation, index total quantity) - values used in the function 'thresholdDerivatives'
    // {2,3} = mass in significant objects, {4,5} = volume in significant objects, {7,8} = mass in all objects, {9,10} = volume in all objects
    thresholdDerivatives( thresholdResponse, indices, 4 );
    
    
    // compute the optimal threshold using the fraction of valid objects
    vector<Real> x, y;
    for (size_t i=0; i<noThresholds; ++i)
    {
        x.push_back( log10( thresholdResponse(i,0) ) );
        y.push_back( thresholdResponse(i,1) );
    }
    Real optimalThreshold_validFraction, optimalFraction = 0.5;
    bool succesfulThreshold = optimalThreshold( x, y, optimalFraction, true, &optimalThreshold_validFraction );
    optimalThreshold_validFraction = pow( double(10.), double(optimalThreshold_validFraction) ); 
    
    
    // compute the optimal threshold using the peak of the mass change
    y.clear();
    for (size_t i=0; i<noThresholds; ++i)
        y.push_back( thresholdResponse(i,2)*thresholdResponse(i,3) );
    Real optimalThreshold_massChange, fitParams[3];
    smoothData( y, .1 );
    parabolaFitMaximum( x, y, .5, &optimalThreshold_massChange, fitParams );
    optimalThreshold_massChange = pow( double(10.), double(optimalThreshold_massChange) ); 
    
    
    if ( succesfulThreshold ) cout << "\nOptimal threshold value: " << setprecision(4) << optimalThreshold_validFraction << "\n";
    else cout << "# Optimal threshold value: NOT FOUND\n";
    
    
    
    // open the output file and write the data
    cout << "Writting the results to the ASCII file '" << outputFilename << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    
    if ( succesfulThreshold )
        outputFile << "# Optimal threshold value: " << optimalThreshold_validFraction << "\n";
    else outputFile << "# Optimal threshold value: NOT_FOUND\n";
    outputFile << "# The optimal node threshold value was computed by requiring that " << optimalFraction << " of the significant nodes (i.e. objects of mass larger than the threshold) have a density larger than the virial one.\n";
    outputFile << "# Alternative value for optimal threshold value: " << optimalThreshold_massChange << "\n"
               << "# The above value is the peak of 3rd*4th columns as a function of threshold (only for significant objects). The results were obtain by fitting a 2nd degree polynomial to the data.\n"
               << "# Best fit function: " << fitParams[0] << " + " << fitParams[1] << "*x + " << fitParams[2] << "*x**2\n";
    outputFile << "# The file contains the following data (the mass and volume units are in fractions of total box mass and volume):\n"
               << "#\t 1st column = the threshold value\n"
               << "#\t 2nd column = the fraction of valid nodes ( significant objects with average density > viraila density / number of significant objects )\n"
               << "#\t 3rd column = the mass change of significant objects with threshold \n"
               << "#\t 4th column = the mass in significant objects\n"
               << "#\t 5th column = the volume change of significant objects with threshold\n"
               << "#\t 6th column = the volume in significant objects\n"
               << "#\t 7th column = the number of significant nodes (mass higher than the threshold one) \n"
               << "#\t 8th column = the mass change of all objects with threshold\n"
               << "#\t 9th column = the mass in all objects\n"
               << "#\t 10th column= the volume change of all objects with threshold\n"
               << "#\t 11th column= the volume in all objects\n";
    outputFile << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
    writeData( thresholdResponse, outputFile ); //write the data to file
    outputFile.close();
    cout << "Done.\n";
    
    if ( not succesfulThreshold )
        throwError( "Could not find the optimal threshold for nodes in function 'outputNodeThreshold'." );
    return optimalThreshold_validFraction;
}



/* Computes the optimal filament and wall threshold. We define significant fialments/walls the ones with volume higher than 'minimumSize'.
The optimal filament/wall treshold is defined as the value when we have percolation - defined as when the most massive object takes 50% of the total volume in significant objects.
It prints the data and the fit details to an ASCII file.
*/
Real outputFilamentThreshold(Array<Real,3> &density,
                             Array<Real,3> &response,
                             int const neighborType,
                             Int const minimumSize,
                             Array<Real,1> &threshold,
                             string outputFilename,
                             string const programOptions)
{
    size_t const noThresholds = threshold.size();
    cout << "\n\nComputing the optimal threshold for filaments/walls using " << noThresholds << " thresholds values and 'minimum object volume' = " << minimumSize << ":\n" << flush;
    cout << "  threshold index   threshold value   largest object(\%)\n";
    
    
    // assign memory for the quantities needed by the computation
    Array<Real,2> thresholdResponse( noThresholds, 11 );    // keeps track of the final results
    Array<int,3>  mask( density.getSize(0), density.getSize(1), density.getSize(2) );         //mask that keeps track of the objects
    mask.assign( -1 );  //initialize the mask to empty cells
    int trackStatistics[] = {0,0};  // array needed by the 'fractionVirialNodes' function
    
    
    // loop over the thresholds and get the data
    Real norm = Real( density.size() ); //normalization factor for mass and volume data - this gives the total mass and volume in the box
    Real minThreshold = threshold[0], maxThreshold = numeric_limits<Real>::max();
    for (size_t i=0; i<noThresholds; ++i)
    {
        if (i%10==0 or i==noThresholds-1)
            cout << "   " << setw(4) << i+1 << " / " << noThresholds << "\t\t" << setprecision(2) << setw(4) << minThreshold << "\t\t" << flush;
        
        Real properties[6];
        thresholdProperties( response, minThreshold, maxThreshold, &mask, neighborType, density, minimumSize, properties, trackStatistics );
        
        thresholdResponse(i,0) = threshold[i];       //threshold value
        thresholdResponse(i,1) = properties[2];      //volume fraction of largest object compared to the significant objects volume
        thresholdResponse(i,3) = properties[0]/norm; //total mass in significant objects (volume>minimumSize)
        thresholdResponse(i,5) = properties[1]/norm; //total volume in significant objects
        thresholdResponse(i,6) = properties[5];      //volume fraction of largest object compared to the all objects volume
        thresholdResponse(i,8) = properties[3]/norm; //total mass in all nodes
        thresholdResponse(i,10)= properties[4]/norm; //total volume in all nodes
        
        maxThreshold = threshold[i];
        if (i<noThresholds-1) minThreshold = threshold[i+1];
        if (i%10==0 or i==noThresholds-1)
            cout << setprecision(2) << setw(4) << thresholdResponse(i,1)*100. << "\n";
    }
    
    
    // get the mass and volume change between two different thresholds
    int indices[][2] = { {2,3}, {4,5}, {7,8}, {9,10} };  //each pair of indices gives the values (index quantity variation, index total quantity) - values used in the function 'thresholdDerivatives'
    // {1,2} = mass in significant objects, {3,4} = volume in significant objects, {6,7} = mass in all objects, {8,9} = volume in all objects
    thresholdDerivatives( thresholdResponse, indices, 4 );
    
    
    // compute the optimal threshold using the peak of the mass change
    vector<Real> x, y;
    for (size_t i=0; i<noThresholds; ++i)
    {
        x.push_back( log10( thresholdResponse(i,0) ) );
        y.push_back( thresholdResponse(i,1) );
    }
    Real optimalThreshold_volumeFraction, optimalFraction = 0.5;
    smoothData( y, .1 );
    bool succesfulThreshold = optimalThreshold( x, y, optimalFraction, false, &optimalThreshold_volumeFraction );
    optimalThreshold_volumeFraction = pow( double(10.), double(optimalThreshold_volumeFraction) );
    
    
    // compute the optimal threshold using the fraction of volume occupied by the largest object
    y.clear();
    for (size_t i=0; i<noThresholds; ++i)
        y.push_back( thresholdResponse(i,2)*thresholdResponse(i,3) );
    Real optimalThreshold_massChange, fitParams[3];
    parabolaFitMaximum( x, y, .75, &optimalThreshold_massChange, fitParams );
    optimalThreshold_massChange = pow( double(10.), double(optimalThreshold_massChange) ); 
    
    
    if ( succesfulThreshold ) cout << "\nOptimal threshold value: " << setprecision(4) << optimalThreshold_volumeFraction << "\n";
    else cout << "# Optimal threshold value: NOT FOUND\n";
    
    
    
    // open the output file and write the data
    cout << "Writting the results to the ASCII file '" << outputFilename << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    
    if ( succesfulThreshold )
        outputFile << "# Optimal threshold value: " << optimalThreshold_volumeFraction << "\n";
    else outputFile << "# Optimal threshold value: NOT_FOUND\n";
    outputFile << "# The above value is the percolation threshold, i.e. the threshold when the largest object takes half of the volume in significant objects.\n";
    outputFile << "# Alternative value for optimal threshold value: " << optimalThreshold_massChange << "\n"
               << "# The optimal filament/wall threshold value was computed by finding the peak of the 3rd*4th columns as a function of threshold (only for significant objects). The results were obtain by fitting a 2nd degree polynomial to the data.\n"
               << "# Best fit function: " << fitParams[0] << " + " << fitParams[1] << "*x + " << fitParams[2] << "*x**2\n";
    outputFile << "# The file contains the following data (the mass and volume units are in fractions of total box mass and volume):\n"
               << "#\t 1st column = the threshold value\n"
               << "#\t 2nd column = the volume fraction of the largest object (fraction with respect to volume in significant objects)\n"
               << "#\t 3rd column = the mass change of significant objects with threshold \n"
               << "#\t 4th column = the mass in significant objects\n"
               << "#\t 5th column = the volume change of significant objects with threshold\n"
               << "#\t 6th column = the volume in significant objects\n"
               << "#\t 7th column = the volume fraction of the largest object (fraction with respect to volume in all objects)\n"
               << "#\t 8th column = the mass change of all objects with threshold\n"
               << "#\t 9th column = the mass in all objects\n"
               << "#\t 10th column= the volume change of all objects with threshold\n"
               << "#\t 11th column= the volume in all objects\n";
    outputFile << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
    writeData( thresholdResponse, outputFile ); //write the data to file
    outputFile.close();
    cout << "Done.\n";
    
    if ( not succesfulThreshold )
        throwError( "Could not find the optimal threshold for filaments/walls in function 'outputFilamentThreshold'." );
    return optimalThreshold_volumeFraction;
}




/* Computes the optimal filament and wall threshold. We define significant fialments/walls the ones with volume higher than 'minimumSize'.
The optimal filament/wall treshold is defined as the response value where the is the peak of the mass derivative as a function of threshold. We find the peak by fitting a 2nd order polynomial to the data.
It prints the data and the fit details to an ASCII file.
*/
Real outputWallThreshold(Array<Real,3> &density,
                         Array<Real,3> &response,
                         int const neighborType,
                         Int const minimumSize,
                         Array<Real,1> &threshold,
                         string outputFilename,
                         string const programOptions)
{
    size_t const noThresholds = threshold.size();
    cout << "\n\nComputing the optimal threshold for filaments/walls using " << noThresholds << " thresholds values and 'minimum object volume' = " << minimumSize << ":\n" << flush;
    cout << "  threshold index   threshold value   largest object(\%)    mass_derivative\n";
    
    
    // assign memory for the quantities needed by the computation
    Array<Real,2> thresholdResponse( noThresholds, 11 );    // keeps track of the final results
    Array<int,3>  mask( density.getSize(0), density.getSize(1), density.getSize(2) );         //mask that keeps track of the objects
    mask.assign( -1 );  //initialize the mask to empty cells
    int trackStatistics[] = {0,0};  // array needed by the 'fractionVirialNodes' function
    
    
    // loop over the thresholds and get the data
    Real norm = Real( density.size() ); //normalization factor for mass and volume data - this gives the total mass and volume in the box
    Real minThreshold = threshold[0], maxThreshold = numeric_limits<Real>::max();
    for (size_t i=0; i<noThresholds; ++i)
    {
        if (i%10==0 or i==noThresholds-1)
            cout << "   " << setw(4) << i+1 << " / " << noThresholds << "\t\t" << setprecision(2) << setw(4) << minThreshold << "\t\t" << flush;
        
        Real properties[6];
        thresholdProperties( response, minThreshold, maxThreshold, &mask, neighborType, density, minimumSize, properties, trackStatistics );
        
        thresholdResponse(i,0) = threshold[i];       //threshold value
        thresholdResponse(i,2) = properties[0]/norm; //total mass in significant objects (volume>minimumSize)
        thresholdResponse(i,4) = properties[1]/norm; //total volume in significant objects
        thresholdResponse(i,5) = properties[2];      //volume fraction of largest object compared to the significant objects volume
        thresholdResponse(i,7) = properties[3]/norm; //total mass in all nodes
        thresholdResponse(i,9) = properties[4]/norm; //total volume in all nodes
        thresholdResponse(i,10)= properties[5];      //volume fraction of largest object compared to the all objects volume
        
        maxThreshold = threshold[i];
        if (i<noThresholds-1) minThreshold = threshold[i+1];
        if (i%10==0 or i==noThresholds-1)
            cout << setprecision(2) << setw(4) << thresholdResponse(i,5)*100. << "\t\t" << (i==0 ? 0. : thresholdResponse(i,2)-thresholdResponse(i-1,2)) << "\n";
    }
    
    
    // get the mass and volume change between two different thresholds
    int indices[][2] = { {1,2}, {3,4}, {6,7}, {8,9} };  //each pair of indices gives the values (index quantity variation, index total quantity) - values used in the function 'thresholdDerivatives'
    // {1,2} = mass in significant objects, {3,4} = volume in significant objects, {6,7} = mass in all objects, {8,9} = volume in all objects
    thresholdDerivatives( thresholdResponse, indices, 4 );
    
    
    // compute the optimal threshold using the peak of the mass change
    vector<Real> x, y;
    for (size_t i=0; i<noThresholds; ++i)
    {
        x.push_back( log10( thresholdResponse(i,0) ) );
        y.push_back( thresholdResponse(i,1)*thresholdResponse(i,2) );
    }
    Real optimalThreshold_massChange, fitParams[3];
    smoothData( y, .1 );
    bool succesfulThreshold = parabolaFitMaximum( x, y, .75, &optimalThreshold_massChange, fitParams );
    optimalThreshold_massChange = pow( double(10.), double(optimalThreshold_massChange) ); 
    
    
    // compute the optimal threshold using the fraction of volume occupied by the largest object
    y.clear();
    for (size_t i=0; i<noThresholds; ++i)
        y.push_back( thresholdResponse(i,5) );
    Real optimalThreshold_volumeFraction, optimalFraction = 0.5;
    optimalThreshold( x, y, optimalFraction, false, &optimalThreshold_volumeFraction );
    optimalThreshold_volumeFraction = pow( double(10.), double(optimalThreshold_volumeFraction) );
    
    
    if ( succesfulThreshold ) cout << "\nOptimal threshold value: " << setprecision(4) << optimalThreshold_massChange << "\n";
    else cout << "# Optimal threshold value: NOT FOUND\n";
    
    
    
    // open the output file and write the data
    cout << "Writting the results to the ASCII file '" << outputFilename << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    
    if ( succesfulThreshold )
        outputFile << "# Optimal threshold value: " << optimalThreshold_massChange << "\n";
    else outputFile << "# Optimal threshold value: NOT_FOUND\n";
    outputFile << "# The optimal filament/wall threshold value was computed by finding the peak of the 3rd*4th columns (the mass derivative * mass) as a function of threshold (only for significant objects). The results were obtain by fitting a 2nd degree polynomial to the data.\n"
               << "# Best fit function: " << fitParams[0] << " + " << fitParams[1] << "*x + " << fitParams[2] << "*x**2\n";
    outputFile << "# Alternative value for optimal threshold value: " << optimalThreshold_volumeFraction << "\n"
               << "# The above value is the percolation threshold, i.e. the threshold when the largest object takes half of the volume in significant objects.\n";
    outputFile << "# The file contains the following data (the mass and volume units are in fractions of total box mass and volume):\n"
               << "#\t 1st column = the threshold value\n"
               << "#\t 2nd column = the mass change of significant objects with threshold \n"
               << "#\t 3rd column = the mass in significant objects\n"
               << "#\t 4th column = the volume change of significant objects with threshold\n"
               << "#\t 5th column = the volume in significant objects\n"
               << "#\t 6th column = the volume fraction of the largest object (fraction with respect to volume in significant objects)\n"
               << "#\t 7th column = the mass change of all objects with threshold\n"
               << "#\t 8th column = the mass in all objects\n"
               << "#\t 9th column = the volume change of all objects with threshold\n"
               << "#\t 10th column= the volume in all objects\n"
               << "#\t 11th column= the volume fraction of the largest object (fraction with respect to volume in all objects)\n";
    outputFile << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
    writeData( thresholdResponse, outputFile ); //write the data to file
    outputFile.close();
    cout << "Done.\n";
    
    if ( not succesfulThreshold )
        throwError( "Could not find the optimal threshold for filaments/walls in function 'outputWallThreshold'." );
    return optimalThreshold_massChange;
}







// /* Fits the blob response function with two power laws and finds the intersection point of the two functions. It prints the data and the fit details to an ASCII file. */
//Real outputNodeThreshold(Array<Real,3> &density,
//                         Array<Real,3> &response,
//                         int const neighborType,
//                         Real const virialDensity,
//                         Real const minimumSize,
//                         Array<Real,1> &threshold,
//                         string outputFilename,
//                         string const programOptions)
//{
//    // Find the fraction of virial objects for a given threshold value
//    size_t const noThresholds = threshold.size();
//    cout << "\n\nComputing the optimal threshold for nodes using " << noThresholds << " thresholds values, with 'virial density' = " << virialDensity << " and 'minSize' = " << minimumSize << " in units of average mass in voxel:\n" << flush;
//    cout << "  threshold index    threshold value     valid nodes(\%)\n";
//    Array<Real,1> thresholdResponse( noThresholds );
//    for (size_t i=0; i<noThresholds; ++i)
//    {
//        cout << "   " << setw(4) << i+1 << " / " << noThresholds << "\t\t" << setprecision(2) << setw(4) << threshold[i] << "\t\t" << flush;
//        thresholdResponse[i] = fractionVirialBlobs( density, response, threshold[i], neighborType, virialDensity, minimumSize );
//        cout << setprecision(2) << setw(4) << thresholdResponse[i]*100. << "\n";
//    }
    
    
//     // compute the optimal threshold
//    Array<Real,1> tempFraction(noThresholds);
//    for (size_t i=0; i<noThresholds; ++i) tempFraction[i] = thresholdResponse[i];
//    int position = -1;
//    Real optimalFraction = Real(0.5);
//    for (size_t i=0; i<noThresholds-1; ++i)
//        if ( tempFraction[i]<=optimalFraction and tempFraction[i+1]>optimalFraction )
//        {
//            position = i;
//            break;
//        }
//    Real optimalThreshold = Real(0.);
//    if ( position!=-1 ) // approximate as linear the tempFraction function between the two thresholds below and above 'optimalFraction'
//    {
//        Real slope = (tempFraction[position] - tempFraction[position-1]) / (threshold[position] - threshold[position-1]);
//        Real constant = tempFraction[position] - slope * threshold[position];
//        optimalThreshold = (optimalFraction - constant) / slope;
//        cout << "\nOptimal threshold value: " << setprecision(4) << optimalThreshold << "\n";
//    }
//    else cout << "# Optimal threshold value: None\n";
    
    
    
//    // open the output file and write the data
//    cout << "Writting the results to the ASCII file '" << outputFilename << "' ... " << flush;
//    fstream outputFile;
//    openOutputTextFile( outputFile, outputFilename );   //open the output file
//    BOOST_ASSERT( noParameters==4 );
    
//    if ( position!=-1 )
//        outputFile << "# Optimal threshold value: " << optimalThreshold << "\n";
//    else outputFile << "# Optimal threshold value: None\n";
//    outputFile << "# The file contains the results of finding the optimal threshold for node detection with MMF. The 1st column gives the threshold values while the second gives the fraction of nodes which have an average density above the virial density.\n"
//            << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
//    for (size_t i=0; i<noThresholds; ++i)
//        outputFile << threshold[i] << "\t" << thresholdResponse[i] << "\n";
//    outputFile.close();
//    cout << "Done.\n";
//}



// /* Fits the filament response function with a quadratic function around the maximum value of 'thresholdResponseFunction'. It prints the data and the fit details to an ASCII file. */
//void outputFilamentThreshold(Array<Real,3> &response,
//                             int const neighborType,
//                             Real const minimumSize,
//                             Array<Real,1> &threshold,
//                             string outputFilename,
//                             string const programOptions)
//{
//    // Find the optimal filament threshold
//    size_t const noThresholds = threshold.size();
//    size_t gridSize = response.size();
//    cout << "\n\nComputing the optimal threshold for filaments using " << noThresholds << " thresholds values and 'minimum object size' = " << minimumSize << ":\n" << flush;
//    cout << "  threshold index   threshold value   largest filament(\%)\n";
//    Array<Real,1> largestObjectVolumeFraction(noThresholds), // largest object volume fraction compared to large objects
//            volumeFill(noThresholds);                        // volume fill of all objects
//    Array<size_t,1> largestObjectVolume(noThresholds),       // largest object volume
//            noLargeObjects(noThresholds);                    // number large objects above 'minimumSize'
    
//    for (size_t i=0; i<noThresholds; ++i)
//    {
//        cout << "   " << setw(4) << i+1 << " / " << noThresholds << "\t\t  " << setprecision(2) << setw(4) << threshold[i] << "\t\t" << flush;
//        largestObjectVolumeFraction[i] = countObjects( response, threshold[i], neighborType, minimumSize, &(largestObjectVolume[i]), &(noLargeObjects[i]), &(volumeFill[i]) );
//        cout << setprecision(2) << setw(4) << largestObjectVolume[i]/(volumeFill[i]*gridSize)*100. << "\n";
//    }
    
    
//    // compute the optimal threshold
//    Array<Real,1> tempFraction(noThresholds);
//    for (size_t i=0; i<noThresholds; ++i) tempFraction[i] = largestObjectVolume[i]/(volumeFill[i]*gridSize);
//    int position = -1;
//    Real percolationFraction = Real(0.5);
//    for (size_t i=noThresholds-1; i>0; --i)
//        if ( tempFraction[i]<=percolationFraction and tempFraction[i-1]>percolationFraction )
//        {
//            position = i;
//            break;
//        }
//    Real optimalThreshold = Real(0.);
//    if ( position!=-1 ) // approximate as linear the tempFraction function between the two thresholds below and above 'percolationFraction'
//    {
//        Real slope = (tempFraction[position] - tempFraction[position-1]) / (threshold[position] - threshold[position-1]);
//        Real constant = tempFraction[position] - slope * threshold[position];
//        optimalThreshold = (percolationFraction - constant) / slope;
//        cout << "\nOptimal threshold value: " << setprecision(4) << optimalThreshold << "\n";
//    }
//    else cout << "# Optimal threshold value: None\n";
    
    
    
//    ArrayProperties<size_t,size_t> prop( noLargeObjects.ptrData(), noLargeObjects.size() );
//    Real maxNoLargeObjects = Real(prop.valueMaximum());
    
    
//    // open the output file and write the data
//    cout << "Writting the results to the ASCII file '" << outputFilename << "' ... " << flush;
//    fstream outputFile;
//    openOutputTextFile( outputFile, outputFilename );   //open the output file
//    BOOST_ASSERT( noParameters==3 );
    
//    if ( position!=-1 )
//        outputFile << "# Optimal threshold value: " << optimalThreshold << "\n";
//    else outputFile << "# Optimal threshold value: None\n";
//    outputFile << "# The file contains the results of finding the optimal threshold for filament detection with MMF.The columns give the following information:\n"
//            << "#\t 1st = threshold values\n"
//            << "#\t 2nd = ratio of the volume of the largest object to the volume in all objects\n"
//            << "#\t 3rd = ratio of the volume of the largest object to the volume in large objects\n"
//            << "#\t 4th = the normalized number of large objects for the given threshold\n"
//            << "#\t 5th = the volume fraction filled by all objects for the given threshold\n"
//            << "#\t 6th = the absolute size of the largest object\n"
//            << "#\t 7th = the absolute numbers for column 4\n"
//            << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
//    for (size_t i=0; i<noThresholds; ++i)
//        outputFile << threshold[i] << "\t" << largestObjectVolume[i]/(volumeFill[i]*gridSize) << "\t" << largestObjectVolumeFraction[i] << "\t" << noLargeObjects[i]/maxNoLargeObjects << "\t" << volumeFill[i] << "\t" << largestObjectVolume[i] << "\t" << noLargeObjects[i] << "\n";
//    outputFile.close();
//    cout << "Done.\n";
//}



// /* Fits the wall response function with a quadratic function around the maximum value of 'thresholdResponseFunction'. It prints the data and the fit details to an ASCII file. */
//void outputWallThreshold(Array<Real,3> &response,
//                         int const neighborType,
//                         Real const minimumSize,
//                         Array<Real,1> &threshold,
//                         string outputFilename,
//                         string const programOptions)
//{
//    // Find the optimal wall threshold
//    size_t const noThresholds = threshold.size();
//    size_t gridSize = response.size();
//    cout << "\n\nComputing the optimal threshold for walls using " << noThresholds << " thresholds values and 'minimum object size' = " << minimumSize << ":\n" << flush;
//    cout << "  threshold index   threshold value   largest wall(\%)\n";
//    Array<Real,1> largestObjectVolumeFraction(noThresholds), // largest object volume fraction compared to large objects
//            volumeFill(noThresholds);                        // volume fill of all objects
//    Array<size_t,1> largestObjectVolume(noThresholds),       // largest object volume
//            noLargeObjects(noThresholds);                    // number large objects above 'minimumSize'
    
//    for (size_t i=0; i<noThresholds; ++i)
//    {
//        cout << "   " << setw(4) << i+1 << " / " << noThresholds << "\t\t  " << setprecision(2) << setw(4) << threshold[i] << "\t\t" << flush;
//        largestObjectVolumeFraction[i] = countObjects( response, threshold[i], neighborType, minimumSize, &(largestObjectVolume[i]), &(noLargeObjects[i]), &(volumeFill[i]) );
//        cout << setprecision(2) << setw(4) << largestObjectVolume[i]/(volumeFill[i]*gridSize)*100. << "\n";
//    }
    
    
//    // compute the optimal threshold
//    Array<Real,1> tempFraction(noThresholds);
//    for (size_t i=0; i<noThresholds; ++i) tempFraction[i] = largestObjectVolume[i]/(volumeFill[i]*gridSize);
//    int position = -1;
//    Real percolationFraction = Real(0.5);
//    for (size_t i=noThresholds-1; i>0; --i)
//        if ( tempFraction[i]<=percolationFraction and tempFraction[i-1]>percolationFraction )
//    {
//        position = i;
//        break;
//    }
//    Real optimalThreshold = Real(0.);
//    if ( position!=-1 ) // approximate as linear the tempFraction function between the two thresholds below and above 'percolationFraction'
//    {
//        Real slope = (tempFraction[position] - tempFraction[position-1]) / (threshold[position] - threshold[position-1]);
//        Real constant = tempFraction[position] - slope * threshold[position];
//        optimalThreshold = (percolationFraction - constant) / slope;
//        cout << "\nOptimal threshold value: " << optimalThreshold << "\n";
//    }
//    else cout << "# Optimal threshold value: None\n";
    
    
    
//    ArrayProperties<size_t,size_t> prop( noLargeObjects.ptrData(), noLargeObjects.size() );
//    Real maxNoLargeObjects = Real(prop.valueMaximum());
    
    
//    // open the output file and write the data
//    cout << "Done.\nWritting the results to the ASCII file '" << outputFilename << "' ... " << flush;
//    fstream outputFile;
//    openOutputTextFile( outputFile, outputFilename );   //open the output file
//    BOOST_ASSERT( noParameters==3 );
    
//    if ( position!=-1 )
//        outputFile << "# Optimal threshold value: " << optimalThreshold << "\n";
//    else outputFile << "# Optimal threshold value: None\n";
//    outputFile << "# The file contains the results of finding the optimal threshold for wall detection with MMF.The columns give the following information:\n"
//            << "#\t 1st = threshold values\n"
//            << "#\t 2nd = ratio of the volume of the largest object to the volume in all objects\n"
//            << "#\t 3rd = ratio of the volume of the largest object to the volume in large objects\n"
//            << "#\t 4th = the normalized number of large objects for the given threshold\n"
//            << "#\t 5th = the volume fraction filled by all objects for the given threshold\n"
//            << "#\t 6th = the absolute size of the largest object\n"
//            << "#\t 7th = the absolute numbers for column 4\n"
//            << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
//    for (size_t i=0; i<noThresholds; ++i)
//        outputFile << threshold[i] << "\t" << largestObjectVolume[i]/(volumeFill[i]*gridSize) << "\t" << largestObjectVolumeFraction[i] << "\t" << noLargeObjects[i]/maxNoLargeObjects << "\t" << volumeFill[i] << "\t" << largestObjectVolume[i] << "\t" << noLargeObjects[i] << "\n";
//    outputFile.close();
//    cout << "Done.\n";
//}


/* This function computes the histogram for a array and outputs the results to an ASCII file. */
void outputResponseHessian(Real *ptr,
                           Int const size,
                           Real const *histogramProperties,
                           string outputFilename,
                           string const programOptions)
{
    // first compute the histogram
    cout << "Computing the histogram of the response data using " << int(histogramProperties[2]) << " logarithmic bins in the interval " << histogramProperties[0] << " to " << histogramProperties[1] << " of the maximum response value:\n";
    ArrayProperties<Real,Int> prop( ptr, size );
    Real maximum = prop.valueMaximum();       // get the maximum response
    Int const histogramBins = Int( histogramProperties[2] );
    cout << "\t maximum response: " << maximum << "\n" << flush;
    prop.histogram( histogramProperties[0]*maximum, histogramProperties[1]*maximum, histogramBins, LOGARITHMIC_BIN );
    vector<Int> histValues = prop.histogramValues();
    vector<Real> histBins = prop.histogramBins();
    
    // compute the cumulative sum of the histogram values
    ArrayProperties<Int,Int> prop2( &(histValues[0]), histValues.size() );
    vector<double> histCumulativeSum = prop2.cumulativeSum( ORDER_DESCENDING );
    cout << "\t filled volume: " << histCumulativeSum[0] / Real(size)*100. << "\%\n";
    
    
    cout << "Fitting the data with a log-normal distribution and writing the results to the ASCII file '" << outputFilename << "' ..." << flush;
    Real histogramMaxValue = Real( prop2.valueMaximum() );
    Int indexMax = prop2.indexMaximum();
    // fit around the maximum of histogram - take 4 elements before maximum to reduce the effects of the noise for smaller values
    size_t fitExtent = 5;
    size_t indexFitStart = indexMax-fitExtent, indexFitEnd = indexMax+3*fitExtent ;
    if ( indexFitStart<0 )
    {
        throwWarning( "You should increase the histogram minimum range since it is close to the maximum value of the histogram and will give wrong results for the lognormal fit." );
        indexFitStart = 0;
    }
    if ( indexFitEnd>histogramBins )
    {
        throwWarning( "You should increase the histogram maximum range since it is close to the maximum value of the histogram and will give wrong results for the lognormal fit." );
        indexFitEnd = histogramBins;
    }
    size_t const fitSize = indexFitEnd - indexFitStart;
    
    // do the fit
    Real tempY[fitSize];
    for (Int i=indexFitStart; i<indexFitEnd; ++i)
        tempY[i-indexFitStart] = Real( histValues[i] );
    logNormal fit;
    fit.dataToFit( &(histBins[indexFitStart]), tempY, fitSize );
    size_t const noParameters = fit.noFitParamaters();
    Real parameters[noParameters];
    for (size_t i=0; i<noParameters; ++i)
        parameters[i] = fit.param_n( i );
    
    
    // open the output file and write the data
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    BOOST_ASSERT( noParameters==3 );
    
    outputFile << "# The file contains the histogram of a MMF response. The columns give the following information:\n"
            << "#\t 1st = the value at the center of the bin\n"
            << "#\t 2nd = the count number in each bin\n"
            << "#\t 3rd = the normalized histogram to the largest histogram count value\n"
            << "#\t 4th = the descending cumulative sum of the histogram\n"
            << "#\t 5th = same as 4th column, but normalized to the total number of grid cells\n"
            << "# The results were obtained using the command:  " << programOptions << "\n\n"
            << "# The histogram was fitted with a lognormal distribution - only " << fitExtent << " points to the left  and " << 2*fitExtent << " points to the right of the maximum where considered for the fit to reduce noise influence; the fit is:\n"
            << "# f1(x) = " << parameters[0] << " * exp( " << parameters[1] <<  " * (log10(x)-" << parameters[2] << ")**2 )\n"
            << "# The maximum of the fit function is at x = " << pow( Real(10.), parameters[2] ) << " .\n\n";
    
    for (size_t i=0; i<histogramBins; ++i)
        outputFile << histBins[i] << "\t" << histValues[i] << "\t" << histValues[i]/histogramMaxValue << "\t" << histCumulativeSum[i] << "\t" << histCumulativeSum[i] / Real(size) << "\n";
    outputFile.close();
    cout << "Done.\n";
}

