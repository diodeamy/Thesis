#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>

#include <defines.h>
#include <array.h>
#include <miscellaneous.h>
#include <arrayProperties.h>
using namespace std;



/* Outputs the density histogram. */
void outputDensityHistogram(Array<Real,3> &density,
                            vector<Real> range,
                            string outputFilename,
                            string const programOptions)
{
    ArrayProperties<Real,Int> prop( density.ptrData(), density.size() );
    //check if the user supplied a range
    if ( range[3]==0. ) //default range
    {
        range[0] = prop.valueMinimum(); // gets the minimum value of density
        range[1] = prop.valueMaximum(); // gets the maximum value of density
        range[2] = 100.;                // sets 100 bins
    }
    size_t const noBins = size_t(range[2]);
    Real const averageDensity = prop.mean();
    cout << "\nAverage density: " << averageDensity << "\n";
    
    
    // compute the histogram
    cout << "Computing the density histogram:\n"
            << "\t range   : " << range[0] << " to " << range[1] << "\n"
            << "\t bins    : " << noBins << "\n" << flush;
    prop.histogram( range[0], range[1]+1., noBins, LOGARITHMIC_BIN );
    vector<Int> histValues = prop.histogramValues();
    vector<Real> histBins = prop.histogramBins();
    
    
    // compute the cumulative sum of the density histogram values
    ArrayProperties<Int,Int> prop2( &(histValues[0]), histValues.size() );
    Real maxHistogramValue = Real( prop2.valueMaximum() );
    vector<double> histCumulativeSum = prop2.cumulativeSum( ORDER_ASCENDING );
    
    
    // compute the mass vs. density graph
    prop.histogramWeighted( range[0], range[1]+1., noBins, density.ptrData(), LOGARITHMIC_BIN );
    vector<double> massHistogram = prop.histogramWeightedValues();
    ArrayProperties<double,Int> prop3( &(massHistogram[0]), massHistogram.size() );
    double maxMassHistogramValue = prop3.valueMaximum();
    vector<double> histMassCumulativeSum = prop3.cumulativeSum( ORDER_ASCENDING );
    Real totalSize = Real( density.size() );
    
    
    // output the results
//     Real medianApprox = prop.medianApproximative( 100, LOGARITHMIC_BIN );
    cout << "Writting the density histogram to the ASCII file '" << outputFilename << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    
    outputFile << "# The file contains the histogram of a density field:\n"
            << "#\t 1st column: the bin density values.\n"
            << "#\t 2nd column: the histogram value.\n"
            << "#\t 3rd column: histogram values normalized to largest value.\n"
            << "#\t 4th column: the cumulative sum of histogram values (divided by number of density cells).\n"
            << "#\t 5th column: delta(mass) / delta(density) - i.e. mass histogram vs. density.\n"
            << "#\t 6th column: delta(mass) / delta(density) - normalized by the highest value.\n"
            << "#\t 7th column: cummulative mass histogram vs. density (divided by number of density cells).\n"
            << "# The results were obtained using the command:  " << programOptions << "\n"
            << "# Number of density cells: " << density.size() << "\n"
            << "# Average density: " << averageDensity << "\n" << "\n";
//             << "# Median density : " << prop.median( LOGARITHMIC_BIN ) << "\n"
//             << "# Median density : " << medianApprox << "  with error: " << prop.medianApproximativeError() << "\n\n";
    
    for (size_t i=0; i<histBins.size(); ++i)
        outputFile << histBins[i] << "\t" << histValues[i] << "\t" << histValues[i]/maxHistogramValue << "\t" << histCumulativeSum[i]/totalSize << "\t" << massHistogram[i] << "\t" << massHistogram[i]/maxMassHistogramValue << "\t" << histMassCumulativeSum[i]/totalSize << "\n";
    
    outputFile.close();
    cout << "Done.\n";
}



/* Outputs the value of the density and the distribution of the ratios for the given density cell if the density cell != 'invalid'. */
void outputPointRatioDistribution(Array<Real,3> &density,
                                  Array<Real,3> &pointRatios,
                                  Real const invalid,
                                  vector<Real> xRange,
                                  vector<Real> yRange,
                                  string outputFilename,
                                  string const programOptions)
{
    int const xBins = int( xRange[2] );
    int const yBins = int( yRange[2] );
    cout << "Computing the distribution of ratio values:\n"
        << "\t x range   : " << xRange[0] << " to " << xRange[1] << "   using " << xBins << "bins\n"
        << "\t y range   : " << yRange[0] << " to " << yRange[1] << "   using " << yBins << "bins\n" << flush;
    
    //compute the distribution
    Array<Real,2> distribution( xBins, yBins ); // will keep track of the distribution
    distribution.assign( 0. );
    Real dx = log10(xRange[1]/xRange[0]) / xBins;
    Real dy = (yRange[1]-yRange[0]) / yBins;
    for (Int i=0; i<density.size(); ++i)
    {
        if ( pointRatios[i]==invalid )
            continue;
        int posX = int( log10(density[i]/xRange[0]) / dx );
        int posY = int( (pointRatios[i]-yRange[0]) / dy );
        if ( posX>=0 and posX<xBins and posY>=0 and posY<yBins )
            distribution(posX,posY) += 1.;
    }
    
    // normalize the distribution to the largest value in each x-bin
    Array<Real,2> normalizedDistribution( xBins, yBins );
    for (int i=0; i<xBins; ++i)
    {
        ArrayProperties<Real,size_t> prop( &(distribution(i,0)), yBins );
        Real maxValue = prop.valueMaximum();
        for (int j=0; j<yBins; ++j)
            normalizedDistribution(i,j) = distribution(i,j) / maxValue;
    }
    
    // compute the x and y bins values
    Array<Real,1> xBinValues(xBins), yBinValues(yBins);
    for (int i=0; i<xBins; ++i)
        xBinValues[i] = xRange[0] * pow( 10., (i+0.5)*dx );
    for (int i=0; i<yBins; ++i)
        yBinValues[i] = yRange[0] + (i+0.5) * dy;
    
    
    
    cout << "Writting the distribution of ratio values to the ASCII file '" << outputFilename << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    
    outputFile << "# The file contains the distribution of ratio of two values for each grid cell of a density file:\n"
            << "#\t 1st column: density bin value.\n"
            << "#\t 2nd column: ratio bin value.\n"
            << "#\t 3rd column: number of grid cells with density and ratio in the given bins.\n"
            << "#\t 4th column: normalized column 3 values to the largest value for each density bin.\n"
            << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
    for (int i=0; i<xBins; ++i)
    {
        for (int j=0; j<yBins; ++j)
            outputFile << xBinValues[i] << "\t" << yBinValues[j] << "\t" << distribution(i,j) << "\t" << normalizedDistribution(i,j) << "\n";
        outputFile << "\n";
    }
    
    
    outputFile.close();
    cout << "Done.\n";
}



/* Outputs the statistics of comparing two different density distributions. */
void outputDensityRatio(Array<Real,3> &density,
                        Array<Real,3> &mainDensity,
                        vector<Real> range,
                        vector<Real> ratioValues,
                        string outputFilename,
                        string const programOptions)
{
    if ( density.size()!=mainDensity.size() ) throwError( "The two density distributions in function 'outputDensityRatio' have different sizes. The comparison cannot be made." );
    sort( ratioValues.begin(), ratioValues.end() );
    cout << "\nComputing the ratio of two density distributions. The output file will contain statistics for the full distribution as well as the number of cells with error values larger than: ";
    for (size_t i=0; i<ratioValues.size(); ++i)
        cout << ratioValues[i] << "%, ";
    cout << "when compared with the expected value of 1.\n" << flush;
    
    
    ArrayProperties<Real,Int> prop( mainDensity.ptrData(), mainDensity.size() );
    if ( range[3]==0. ) //default range
    {
        range[0] = prop.valueMinimum(); // gets the minimum value of density
        range[1] = prop.valueMaximum(); // gets the maximum value of density
        range[2] = 100.;                // sets 100 bins
        if ( range[0]==0. )
            range[0] = range[1] * 1.e-6;
    }
    size_t const noBins = size_t(range[2]);
    cout << "Computing the density ratio in the density range:\n"
            << "\t range   : " << range[0] << " to " << range[1] << "\n"
            << "\t bins    : " << noBins << "\n" << flush;
    
    
    // now compute how many grid cells have an error larger than a required value (given in 'ratioValues')
    size_t yBins = ratioValues.size() + 1 ;
    Array<Real,2> distribution( noBins, yBins );
    distribution.assign( Real(0.) );
    Real dx = log10(range[1]/range[0]) / noBins;
    
#define ABSOLUTE(x) (((x)<0)?(-(x)):(x))
    for (size_t i=0; i<density.size(); ++i)
    {
        if ( mainDensity[i]<=Real(0.) )
            continue;
        
        Real temp = ABSOLUTE(density[i] / mainDensity[i] -1.) * Real(100.);
        int posX = int( log10(mainDensity[i]/range[0]) / dx );
        if ( not (posX>=0 and posX<int(noBins)) )
            continue;
        distribution(posX,0) += Real(1.);
        for (size_t j=0; j<ratioValues.size()-1; ++j)
            if ( temp>=ratioValues[j] and temp<ratioValues[j+1] )
            {
                distribution(posX,j+1) += Real(1.);
                break;
            }
        if ( temp>=ratioValues[yBins-2] )
            distribution(posX,yBins-1) += Real(1.);
    }
    
    
    // find the descending cumulative sum
    for (size_t i=0; i<noBins; ++i)
    {
        for (size_t j=ratioValues.size()-1; j>0; --j)
            distribution(i,j) += distribution(i,j+1);
        for (size_t j=ratioValues.size(); j>0; --j)
            if ( distribution(i,0)!=Real(0.) )
                distribution(i,j) /= (distribution(i,0) / 100.);  // normalize the values to percentange of total bin grids
    }
    
    
    // compute the x bins values
    Array<Real,1> xBinValues(noBins);
    for (size_t i=0; i<noBins; ++i)
        xBinValues[i] = range[0] * pow( 10., (i+0.5)*dx );
    
    
    cout << "Writting the distribution of density ratio errors to the ASCII file '" << outputFilename << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    
    outputFile << "# The file contains the distribution of density ratio errors for each grid cell of a density file:\n"
            << "#\t 1st column: density bin value.\n"
            << "#\t 2nd column: total number of density cells in the given density bin.\n"
            << "#\t 3rd column to last column: relative number of grid cells (to total number of grid cells in the given bin) with error of the density ratio larger than the values:  ";
    for (size_t i=0; i<ratioValues.size(); ++i)
        outputFile << ratioValues[i] << "%, ";
    outputFile << "\n# The results were obtained using the command:  " << programOptions << "\n\n";
    
    for (size_t i=0; i<noBins; ++i)
    {
        outputFile << xBinValues[i] << "\t" << distribution(i,0);
        for (size_t j=1; j<yBins; ++j)
            outputFile << "\t" << distribution(i,j);
        outputFile << "\n";
    }
    
    
    outputFile.close();
    cout << "Done.\n";
}




/* Outputs the value of the density and the ratio for the given density cell if the density cell != 'invalid'. */
void outputPointRatio(Array<Real,3> &density,
                      Array<Real,3> &pointRatios,
                      Real const invalid,
                      string outputFilename,
                      string const programOptions)
{
    cout << "Writting the ratio values for each grid point to the ASCII file '" << outputFilename << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    
    outputFile << "# The file contains the ratio of two values for each grid cell of a density file. The 1st column gives the density value while the second gives the value of the ratio.\n"
            << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
    for (Int i=0; i<density.size(); ++i)
        if ( pointRatios[i]!=invalid )
            outputFile << density[i] << "\t" << pointRatios[i] << "\n";
    
    outputFile.close();
    cout << "Done.\n";
}


/* Outputs the statistics of a density comparison. */
void outputPointStatistics(Array< pair< pair<Real,size_t> ,pair<Real,Real> >, 1 >  &statistics,
                           string outputFilename,
                           string const programOptions)
{
    cout << "Writting the statistics for the density comparison to the ASCII file '" << outputFilename << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    
    outputFile << "# The file contains the statistics of a density comparison. The 1st column gives the density value for the corresponding bin, 2nd column gives the number of cells in that bin, 3rd column gives the average value of the result in the given bin while the 4th column gives the standard deviation of the result.\n"
            << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
    for (Int i=0; i<statistics.size(); ++i)
        outputFile << statistics[i].first.first << "\t" << statistics[i].first.second << "\t" << statistics[i].second.first << "\t" << statistics[i].second.second << "\n";
    
    outputFile.close();
    cout << "Done.\n";
}


/* Output statistics of how many random points are in each density grid cell; for each density bin. */
void outputcellCount(Array<Real,3> &density,
                     Array<int,1> &cellCountPoints,
                     vector<Real> range,
                     string outputFilename,
                     string const programOptions)
{
    size_t noBins = size_t( range[2] );
    vector<size_t> yRange;
    size_t const middle = 100, yBins = middle;
    yRange.push_back( middle*0 ); yRange.push_back( middle*3 );
    Array<size_t,1> counts( noBins );
    counts.assign( 0 );
    
    
    // get the density bins
    Real dx = log10(range[1]/range[0]) / noBins;
    Array<Real,1> xBinValues(noBins);
    Array<Real,1> xBinBoundaries(noBins+1);
    for (size_t i=0; i<noBins; ++i)
    {
        xBinValues[i] = range[0] * pow( Real(10.), Real((i+0.5)*dx) );
        xBinBoundaries[i] = range[0] * pow( Real(10.), Real(i*dx) );
    }
    xBinBoundaries[noBins] = range[1];
    // find how many grid cells are in each density bin
    for (size_t i=0; i<density.size(); ++i)
    {
        if ( density[i]>=range[0] and density[i]<range[1] )
        {
            int temp = int(log10(density[i]/range[0]) / dx);
            ++counts[temp];
        }
    }
    
    
    // loop over all density bins and compute the distribution of random points in each density bin
    Array<Int,2> distribution(noBins,yBins);
    vector<double> average;
    vector<int> yBinValues;
    for (size_t i=0; i<noBins; ++i)
    {
        Array<int,1> randomPoints( counts[i] );
        size_t k = 0;
        for (size_t j=0; j<density.size(); ++j)
            if ( density[j]>=xBinBoundaries[i] and density[j]<xBinBoundaries[i+1] )
                randomPoints[k++] = cellCountPoints[j];
        
        ArrayProperties<int,Int> prop( randomPoints.ptrData(), randomPoints.size() );
        prop.histogram( yRange[0], yRange[1], yBins, LINEAR_BIN );
        prop.histogramValues( &(distribution(i,0)), yBins );
        
        average.push_back( prop.mean() );
        if ( i==0 ) yBinValues = prop.histogramBins();
    }
    
    
    cout << "Writting the statistics for the number of random points in each density grid cell to the ASCII file '" << outputFilename << "' ... " << flush;
    fstream outputFile;
    openOutputTextFile( outputFile, outputFilename );   //open the output file
    
    outputFile << "# The file contains the statistics for the number of random points in each density grid cell. The 1st column gives the number of random points, while the rest of the columns give the point distribution for the density bins: \n # density bin   average count\n";
    for (size_t i=0; i<xBinValues.size(); ++i)
        outputFile << "# " << xBinValues[i] << "\t " << average[i] << "\n";
    outputFile << "# The results were obtained using the command:  " << programOptions << "\n\n";
    
    for (size_t i=0; i<yBins; ++i)
    {
        outputFile << yBinValues[i];
        for (size_t j=0; j<noBins; ++j)
            outputFile << "\t" << distribution( j, i );
        outputFile << "\n";
    }
    
    outputFile.close();
    cout << "Done.\n";
}


