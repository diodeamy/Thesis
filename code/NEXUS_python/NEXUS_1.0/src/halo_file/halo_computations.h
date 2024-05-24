#ifndef HALO_COMPUTATIONS_HEADER
#define HALO_COMPUTATIONS_HEADER

#include <iostream>
#include <miscellaneous.h>
#include <defines.h>
#include <array.h>
#include <binarySearch.h>
#include <halo_header.h>


/* Select only a certain number of columns for the output data. 
The function takes the following arguments:
        inputData = a 2D array with the input data
        columns = an array with the indices of the columns from 'inputData' to be copied
        noColumns = number of columns in the 'columns' array
        outputData = a 2D array with the output data
        shiftColumn = should be !=0 only if to write the data in the 'outputData' starting at a different column than 0. Then this gives the index where the first column will be written.
*/
template <typename T, typename T2>
void selectColumns(Array<T,2> &inputData,
                   T2 *columns,
                   int const noColumns,
                   Array<T,2> *outputData,
                   int const shiftColumn = 0)
{
    Int noRows = inputData.axisSize(0);
    if ( noRows>outputData->axisSize(0) )
        throwError( "In function 'selectColumns'. The number of rows in the 'outputData' array is smaller than the number of rows in the input data 'inputData'. Cannot perform the required operation." );
    if ( noColumns>outputData->axisSize(1)-shiftColumn )
        throwError( "In function 'selectColumns'. The number of columns in the 'outputData' array is smaller than the number of columns requested to copy from the input data 'inputData'. Cannot perform the required operation." );
    
    for (Int i=0; i<noRows; ++i)
        for (int j=0; j<noColumns; ++j)
            (*outputData)( i,shiftColumn+j ) = inputData( i,columns[j] );
}

template <typename T, typename T2>
void selectColumns(Array<T,2> &inputData,
                   T2 *columns,
                   int const noColumns,
                   Array<T,1> *outputData)
{
    Int noRows = inputData.axisSize(0);
    if ( noRows>outputData->axisSize(0) )
        throwError( "In function 'selectColumns'. The number of rows in the 'outputData' array is smaller than the number of rows in the input data 'inputData'. Cannot perform the required operation." );
    if ( noColumns>1 )
        throwError( "In function 'selectColumns'. The array 'outputData' has a single column which is smaller than the number of columns requested to copy from the input data 'inputData'. Cannot perform the required operation." );
    
    for (Int i=0; i<noRows; ++i)
        for (int j=0; j<noColumns; ++j)
            (*outputData)[ i ] = inputData( i,columns[j] );
}


/* Select only a certain number of rows for the output data. 
The function takes the following arguments:
        inputData = a 2D array with the input data
        rows = an array with the indices of the rows from 'inputData' to be copied
        noRows = number of rows in the 'rows' array
        outputData = a 2D array with the output data
        shiftRow = should be !=0 only if to write the data in the 'outputData' starting at a different row than 0. Then this gives the index where the first row will be written.
*/
template <typename T, typename T2>
void selectRows(Array<T,2> &inputData,
               T2 *rows,
               Int const noRows,
               Array<T,2> *outputData,
               Int const shiftRow = 0)
{
    int noColumns = inputData.axisSize(1);
    if ( noRows>outputData->axisSize(0)-shiftRow )
        throwError( "In function 'selectColumns'. The number of rows in the 'outputData' array is smaller than the number of rows requested to copy from the input data 'inputData'. Cannot perform the required operation." );
    if ( noColumns>outputData->axisSize(1) )
        throwError( "In function 'selectColumns'. The number of columns in the 'outputData' array is smaller than the number of columns in the input data 'inputData'. Cannot perform the required operation." );
    
    for (Int i=0; i<noRows; ++i)
        for (int j=0; j<noColumns; ++j)
            (*outputData)( i+shiftRow,j ) = inputData( rows[i],j );
}

template <typename T, typename T2>
void selectRows(Array<T,1> &inputData,
               T2 *rows,
               Int const noRows,
               Array<T,1> *outputData,
               Int const shiftRow = 0)
{
    if ( noRows>outputData->axisSize(0)-shiftRow )
        throwError( "In function 'selectColumns'. The number of rows in the 'outputData' array is smaller than the number of rows requested to copy from the input data 'inputData'. Cannot perform the required operation." );
    
    for (Int i=0; i<noRows; ++i)
        (*outputData)[ i+shiftRow ] = inputData[ rows[i] ];
}


/* Select only a certain number of rows and columns for the output data. 
The function takes the following arguments:
        inputData = a 2D array with the input data
        rows = an array with the indices of the rows from 'inputData' to be copied
        noRows = number of rows in the 'rows' array
        columns = an array with the indices of the columns from 'inputData' to be copied
        noColumns = number of columns in the 'columns' array
        outputData = a 2D array with the output data
        shiftRow = should be !=0 only if to write the data in the 'outputData' starting at a different row than 0. Then this gives the index where the first row will be written
        shiftColumn = should be !=0 only if to write the data in the 'outputData' starting at a different column than 0. Then this gives the index where the first column will be written.
*/
template <typename T, typename T2, typename T3>
void selectRowsColumns(Array<T,2> &inputData,
                       T2 *rows,
                       Int const noRows,
                       T3 *columns,
                       int const noColumns,
                       Array<T,2> *outputData,
                       Int const shiftRow = 0,
                       int const shiftColumn = 0)
{
    if ( noRows>outputData->axisSize(0)-shiftRow )
        throwError( "In function 'selectColumns'. The number of rows in the 'outputData' array is smaller than the number of rows requested to copy from the input data 'inputData'. Cannot perform the required operation." );
    if ( noColumns>outputData->axisSize(1)-shiftColumn )
        throwError( "In function 'selectColumns'. The number of columns in the 'outputData' array is smaller than the number of columns requested to copy from the input data 'inputData'. Cannot perform the required operation." );
    
    for (Int i=0; i<noRows; ++i)
        for (int j=0; j<noColumns; ++j)
            (*outputData)( i+shiftRow,j+shiftColumn ) = inputData( rows[i], columns[j] );
}


/* Selects only the indices of the elements in the given range. */
template <typename T>
Int indicesInRange(T *inputData,
                   Int const size,
                   T const minValue, T const maxValue,
                   Int *indices)
{
    Int count = 0;
    for (Int i=0; i<size; ++i)
    {
        if ( inputData[i]>=minValue and inputData[i]<maxValue )
            indices[count++] = i;
    }
    return count;
}




/* Computes the grid indices associated to each halo. */
template <typename T1, typename T2>
void pointGridIndices(Array<T1,2> &position,
                      Box<T1,NO_DIM> dataBox,
                      T2 *gridSize,
                      Array<int,1> *indices)
{
    Real boxLength[NO_DIM], offset[NO_DIM], dx[NO_DIM];
    dataBox.length( boxLength );
    for (int i=0; i<NO_DIM; ++i)
    {
        offset[i] = dataBox[2*i];
        dx[i]     = boxLength[i] / gridSize[i];
    }
    
    // compute the grid index of each point
    for (Int i=0; i<position.axisSize(0); ++i)
    {
        int temp[NO_DIM];
        for (int j=0; j<NO_DIM; ++j)
            temp[j] = std::floor( (position(i,j)-offset[j]) / dx[j] );
        (*indices)[i] = temp[0]*gridSize[1]*gridSize[2] + temp[1]*gridSize[2] + temp[2];
    }
}




/* Return only the unique halos. Duplicate halos are the ones that have the same position and mass within a given tolerance. */
template <typename T1, typename T2>
void uniqueHalos(Halo_header &header,
                 Array<T1,2> &dataIntegers,
                 Array<T2,2> &dataFloats,
                 bool const relabel = true,
                 T2 const tolerance = T2(1.e-4) )
{
    std::cout << "Searching for duplicate halos ...\n" << std::flush;
    T2 tol2 = tolerance * tolerance;
    
    // get the halo position and mass
    Int const noHalos = header.noHalos;
    int columnIndices[] = { header.positionColumns[0], header.positionColumns[1], header.positionColumns[1], header.massColumn };
    Array<Real,2> haloPos( noHalos, NO_DIM );
    Array<Real,1> haloMass( noHalos );
    Array<bool,1> valid( noHalos );
    selectColumns( dataFloats,   &(columnIndices[0]), 3, &haloPos );
    selectColumns( dataFloats,   &(columnIndices[3]), 1, &haloMass );
    
    
    // some variables needed for the binary tree
    Box<Real,NO_DIM> dataBox;
    dataBox.assign( header.box, 2*NO_DIM );
    Real boxLength = dataBox[1]-dataBox[0];     // the boxLength along x-coordinate
    Real radius = boxLength/256;                // the grid spacing for the binary tree
    
    
    // construct a binary tree for fast search
    BinarySearch<T2,NO_DIM> tree( dataBox, radius, false );
    tree.buildTree( haloPos.ptrData(), noHalos );
    
    
    // loop over all halos within the binary tree grid cell and find if any is a duplicate of the halo in question
    Int noValidHalos = 0;
    #pragma omp parallel for reduction(+ :noValidHalos)
    for (Int i=0; i<noHalos; ++i)   //implement the search right here to speed things up
    {
        // get the grid index of the point
        valid[i] = true;
        T2 *pos = &( haloPos(i,0) );
        size_t indexPointCell = tree.getPointGridIndex( pos );
        if ( indexPointCell>=tree.gridSize ) continue;
        
        // loop over all the points in the grid cell
        for (size_t j=tree.offsets[indexPointCell]; j<tree.offsets[indexPointCell]+tree.noCount[indexPointCell]; ++j )
        {
            Int index = tree.indices[j];
            if ( index==i ) continue;                  // skip if neighbor is the halo itself
            
            T2 *pos1 = &( tree.data[ index*NO_DIM ] );
            T2 dis = tree.distance( pos, pos1 );
            T2 deltaMass = std::fabs( haloMass[i] - haloMass[index] ) / haloMass[i];
            if ( dis<tol2 and deltaMass<tolerance and i>index )     //if two halos are duplicate, than get rid of the one with a larger index
                valid[i] = false;
        }
        if ( valid[i] ) ++noValidHalos;
    }
    cout << "\t found " << noValidHalos << " unique halos out of " << noHalos << " input halos -- " << std::setprecision(4) << double(100.*noValidHalos/noHalos) << " % are unique.\n" << std::flush;
    
    if ( noValidHalos==noHalos ) return;    //no duplicates were found
    cout << "\t keeping only valid halos and relabeling the remaining ones ...\n" << std::flush;
    
    
    // get the indices of the valid halos
    Array<Int,1> indices( noValidHalos );
    Int count = 0;
    for (Int i=0; i<noHalos; ++i)
        if ( valid[i] )
            indices[count++] = i;
    
    
    // copy the data to two new arrays
    header.noHalos = noValidHalos;
    Array<T1,2> data1( noValidHalos, header.noColumnsIntegers );
    Array<T2,2> data2( noValidHalos, header.noColumnsFloats );
    selectRows<T1,Int>( dataIntegers, indices.ptrData(),  noValidHalos, &data1 );
    dataIntegers = data1;
    if ( relabel )
        for (Int i=0; i<noValidHalos; ++i)
            dataIntegers(i,0) = T1(i);       //relabel the halos
    selectRows<T2,Int>( dataFloats, indices.ptrData(),  noValidHalos, &data2 );
    dataFloats = data2;
}




#endif
