

#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>
#include <boost/math/special_functions/fpclassify.hpp>


#include "MMF_contraction.h"
#include <vector.h>
#include <matrix.h>
#include <binarySearch.h>
#include <periodic.h>

#ifdef ENABLE_OPENMP
#include <omp.h>
#endif




using namespace std;





// returns the total number of valid 'object' cells
size_t noObjectCells(Array<int,1> &objectSize)
{
    int noCells = 0;
    for (int i=0; i<(int)objectSize.axisSize(0); ++i)
        noCells += objectSize(i);
    return noCells;
}
size_t noObjectCells(Array<int,3> &mask)
{
    size_t noCells = 0;
    Int NX = mask.axisSize(0), NY = mask.axisSize(1), NZ = mask.axisSize(2);
    for (Int i1=0; i1<NX; ++i1)
        for (Int i2=0; i2<NY; ++i2)
            for (Int i3=0; i3<NZ; ++i3)
                if ( mask(i1,i2,i3)>=0 )
                    ++noCells;
                return noCells;
}
size_t noObjectCells(Array<shortInt,3> &mask)
{
    size_t noCells = 0;
    Int NX = mask.axisSize(0), NY = mask.axisSize(1), NZ = mask.axisSize(2);
    for (Int i1=0; i1<NX; ++i1)
        for (Int i2=0; i2<NY; ++i2)
            for (Int i3=0; i3<NZ; ++i3)
                if ( mask(i1,i2,i3)>0 )
                    ++noCells;
                return noCells;
}






// returns the cell indices and positions corresponding to valid object grid cells
void validObjectCell(Array<shortInt,3> &mask,
                     Array<int,2> &cellIndices,
                     Array<Real,2> &spinePosition,
                     Box<Real,NO_DIM> &box)
{
    int object = 0;
    int Nx = mask.axisSize(0), Ny = mask.axisSize(1), Nz = mask.axisSize(2);
    Real dx[NO_DIM] = { (box[1]-box[0])/Nx, (box[3]-box[2])/Ny, (box[5]-box[4])/Nz };
    for (int i1=0; i1<Nx; ++i1)
        for (int i2=0; i2<Ny; ++i2)
            for (int i3=0; i3<Nz; ++i3)
                if ( mask(i1,i2,i3)>0 )
                {
                    cellIndices(object,0) = i1;
                    cellIndices(object,1) = i2;
                    cellIndices(object,2) = i3;
                    spinePosition(object,0) = box[0] + i1*dx[0];
                    spinePosition(object,1) = box[2] + i2*dx[1];
                    spinePosition(object,2) = box[4] + i3*dx[2];
                    ++object;
                }
}
void validObjectCell(Array<shortInt,3> &mask,
                     Array<Real,2> &spinePosition,
                     Box<Real,NO_DIM> &box)
{
    int object = 0;
    int Nx = mask.axisSize(0), Ny = mask.axisSize(1), Nz = mask.axisSize(2);
    Real dx[NO_DIM] = { (box[1]-box[0])/Nx, (box[3]-box[2])/Ny, (box[5]-box[4])/Nz };
    for (int i1=0; i1<Nx; ++i1)
        for (int i2=0; i2<Ny; ++i2)
            for (int i3=0; i3<Nz; ++i3)
                if ( mask(i1,i2,i3)>0 )
                {
                    spinePosition(object,0) = box[0] + i1*dx[0];
                    spinePosition(object,1) = box[2] + i2*dx[1];
                    spinePosition(object,2) = box[4] + i3*dx[2];
                    ++object;
                }
}





/* Takes a distribution of points and computes the shape of that point distribution. After that it uses the shape to compute the eigenvalues and eigenvectors of the shape.
It takes as input:
        point - the point for which we compute the shape of the cloud of points around it
        positions - a vector with the positions
        radius - the length of the shape
        shapeSize - number of cells for the shape along one direction for the 'radius' distance
        maxExtent - the maximum size associated to a point in the input array
*/
void shapeDirections(Real point[],
                     vector<Real> &positions,
                     Real const radius,
                     int const shapeSize,
                     Real const maxExtent,
                     Real *eigenvalues,
                     Real *eigenvectors)
{
    if ( positions.empty() or positions.size()==NO_DIM ) return;
    
    // define an array used to compute the shape of the point cloud
    int const bufferSize = int(maxExtent / radius * shapeSize) + 1;
    int const sSize2 = 2* (shapeSize+bufferSize) + 1;
    bool shapeArray[sSize2][sSize2][sSize2];
    Real const shapeDX = radius / shapeSize;    //width of a shape cell
    Real const shapeDisplacement = shapeDX *  sSize2 / 2.;  //left displacement of the shape array with respect to point position
    
    int const noPoints = positions.size() / NO_DIM;
    Real pointSize = pow( double(radius*radius*radius/noPoints), 1./3. ); //compute the size of a point if they would fill the volume radius^3 uniformly
    if ( pointSize>maxExtent )  // in the case there are only very few points in the neighborhood
        pointSize = maxExtent;
    
    // find the shape of the cloud of points
    for (int i0=0; i0<sSize2; ++i0)
        for (int i1=0; i1<sSize2; ++i1)
            for (int i2=0; i2<sSize2; ++i2)
                shapeArray[i0][i1][i2] = false;
    for (int i=0; i<noPoints; ++i)
    {
        int iStart[NO_DIM], iEnd[NO_DIM];
        for (int j=0; j<NO_DIM; ++j)
        {
            Real pos = positions[i*NO_DIM+j] - point[j] + shapeDisplacement;
            iStart[j] = int( (pos-pointSize) / shapeDX );
            iEnd[j] = int( (pos+pointSize) / shapeDX );
            if ( iStart[j]<0 ) iStart[j] = 0;
            if ( iEnd[j]>=sSize2 ) iEnd[j] = sSize2-1;
        }
        for (int i0=iStart[0]; i0<=iEnd[0]; ++i0)
            for (int i1=iStart[1]; i1<=iEnd[1]; ++i1)
                for (int i2=iStart[2]; i2<=iEnd[2]; ++i2)
                    shapeArray[i0][i1][i2] = true;
    }
    
    
    // compute the reduced momentum of inertia
    Real massCenter[] = {0.,0.,0.};  //the center of the shape
    Real inertiaTensor[] = {0.,0.,0.,0.,0.,0.};//the inertia tensor
    int noShapeCells = 0;
    for (int i0=0; i0<sSize2; ++i0)
        for (int i1=0; i1<sSize2; ++i1)
            for (int i2=0; i2<sSize2; ++i2)
                if ( shapeArray[i0][i1][i2] )
                {
                    massCenter[0] += i0;
                    massCenter[1] += i1;
                    massCenter[2] += i2;
                    ++noShapeCells;
                }
    for (int j=0; j<NO_DIM; ++j)
        massCenter[j] /= noShapeCells;
    for (int i0=0; i0<sSize2; ++i0)
        for (int i1=0; i1<sSize2; ++i1)
            for (int i2=0; i2<sSize2; ++i2)
                if ( shapeArray[i0][i1][i2] )
                {
                    Real tempP[] = { i0-massCenter[0], i1-massCenter[1], i2-massCenter[2] };
                    Real tempD = tempP[0]*tempP[0] + tempP[1]*tempP[1] + tempP[2]*tempP[2];
                    inertiaTensor[0] += tempP[0]*tempP[0] / tempD;
                    inertiaTensor[1] += tempP[0]*tempP[1] / tempD;
                    inertiaTensor[2] += tempP[0]*tempP[2] / tempD;
                    inertiaTensor[3] += tempP[1]*tempP[1] / tempD;
                    inertiaTensor[4] += tempP[1]*tempP[2] / tempD;
                    inertiaTensor[5] += tempP[2]*tempP[2] / tempD;
                }
    
    
    // get the eigenvalues and eigenvectors
    MATRIX::Symmetric<Real,NO_DIM> mat;
    mat.data( inertiaTensor );
    mat.eigenvalues( eigenvalues, eigenvectors, true );
}

// Computes the momentum of inertia for a cloud of points
bool cloudDirections(vector<Real> &positions,
                     Real *massCenter,
                     Real *eigenvalues,
                     Real *eigenvectors)
{
    if ( positions.empty() or positions.size()==NO_DIM ) return false;
    
    // compute the reduced momentum of inertia
    Real inertiaTensor[] = {0.,0.,0.,0.,0.,0.};//the inertia tensor
    int const noPoints = positions.size() / NO_DIM;
    for (int i=0; i<noPoints; ++i)
    {
        Real tempP[] = { positions[i*NO_DIM]-massCenter[0], positions[i*NO_DIM+1]-massCenter[1], positions[i*NO_DIM+2]-massCenter[2] };
        //Real tempD = tempP[0]*tempP[0] + tempP[1]*tempP[1] + tempP[2]*tempP[2];
        Real tempD = 1.;
        inertiaTensor[0] += tempP[0]*tempP[0] / tempD;
        inertiaTensor[1] += tempP[0]*tempP[1] / tempD;
        inertiaTensor[2] += tempP[0]*tempP[2] / tempD;
        inertiaTensor[3] += tempP[1]*tempP[1] / tempD;
        inertiaTensor[4] += tempP[1]*tempP[2] / tempD;
        inertiaTensor[5] += tempP[2]*tempP[2] / tempD;
    }
    
    // get the eigenvalues and eigenvectors
    MATRIX::Symmetric<Real,NO_DIM> mat;
    mat.data( inertiaTensor );
    mat.eigenvalues( eigenvalues, eigenvectors, true );
    return true;
}
bool cloudDirections(vector<Real> &positions,
                     Real *eigenvalues,
                     Real *eigenvectors)
{
    if ( positions.empty() or positions.size()==NO_DIM ) return false;
    
    // compute the reduced momentum of inertia
    Real massCenter[] = {0.,0.,0.};  //the center of the shape
    int const noPoints = positions.size() / NO_DIM;
    for (int i=0; i<noPoints; ++i)
        for (int j=0; j<NO_DIM; ++j)
            massCenter[j] += positions[i*NO_DIM+j];
    for (int j=0; j<NO_DIM; ++j)
        massCenter[j] /= noPoints;
    return cloudDirections( positions, massCenter, eigenvalues, eigenvectors );
}


// Checks to see if the point cloud represents the intersection of several 'well-defined' filaments - up to 4 filaments
bool validIntersection(Datas &data,
                       int const Feature,
                       Real const eigenvalueRatioThreshold,
                       Real const maxAngle = Real(20.))
{
    int const maxDirections2 = 4;
    int const noPoints = data.neighbors.size()/NO_DIM;
    
    // find the directions with the most number of particles
    vector<Real> tempDirs;
    vector< pair<int,int> > trackCount;
    Real cosThreshold = cos( 2.*maxAngle*PI/180. );
    
    int count = 0;
    for (int i=0; i<noPoints; ++i)
    {
        bool newDirection = true;
        if ( !tempDirs.empty() )
        {
            for (int k=0; k<tempDirs.size()/NO_DIM; ++k)
                if ( fabs(cosAngle( &(tempDirs[k*NO_DIM]), &(data.directions[i*NO_DIM]) )) >= cosThreshold )   //angle within the given bin
                {
                    for (int j=0; j<NO_DIM; ++j)
                        tempDirs[k*NO_DIM+j] += data.directions[i*NO_DIM+j];
                    newDirection = false;
                    ++(trackCount[k].first);
                }
        }
        if ( newDirection )
        {
            for (int j=0; j<NO_DIM; ++j)
                tempDirs.push_back( data.directions[i*NO_DIM+j] );
            trackCount.push_back( make_pair(1,count) );
            ++count;
        }
    }
    
    
    // find the 'maxDirections' with the largest number of particles
    sort( trackCount.begin(), trackCount.end() );
    reverse( trackCount.begin(), trackCount.end() );
    int const maxDirections = count>maxDirections2 ? maxDirections2 : count;
    Real avgDirs[NO_DIM*maxDirections];
    for (int i=0; i<maxDirections; ++i)
        for (int j=0; j<NO_DIM; ++j)
            avgDirs[i*NO_DIM+j] = tempDirs[(trackCount[i].second)*NO_DIM+j];
    
    
    // get the points for each of the directions
    vector< vector<Real> > tempPos;
    tempPos.resize(maxDirections+1);
    cosThreshold = cos( maxAngle*PI/180. );
    int objectIndex = -1;
    for (int i=0; i<noPoints; ++i)
    {
        bool newDirection = true;
        for (int k=0; k<maxDirections; ++k)
            if ( fabs(cosAngle( &(avgDirs[k*NO_DIM]), &(data.directions[i*NO_DIM]) )) >= cosThreshold )   //angle within the given bin
            {
                for (int j=0; j<NO_DIM; ++j)
                    tempPos[k].push_back( data.neighbors[i*NO_DIM+j] );
                newDirection = false;
                if (data.objectIndex==i) objectIndex = k;
                break;
            }
        if ( newDirection )
        {
            for (int j=0; j<NO_DIM; ++j)
                tempPos[maxDirections].push_back( data.neighbors[i*NO_DIM+j] );
            if (data.objectIndex==i) objectIndex = maxDirections;
        }
    }
    
    
    // compute the shape of the points for each direction and see if they correspond to well defined filaments/walls
    for (int k=0; k<maxDirections+1; ++k)
    {
        if ( tempPos[k].size()<=1 ) continue;
        Real eigenvalues[NO_DIM], eigenvectors[NO_DIM*NO_DIM];
        cloudDirections( tempPos[k], eigenvalues, eigenvectors);
        Real eigenvalueRatio = Feature==3 ? eigenvalues[1]/eigenvalues[2] : eigenvalues[0]/eigenvalues[1];
        int offset = Feature==3 ? 6 : 0;
        if ( eigenvalueRatio>eigenvalueRatioThreshold )
            return false;
        else if ( k==objectIndex and data.updateDirection )
            for (int j=0; j<NO_DIM; ++j)
                data.spineDirection[j] = eigenvectors[offset+j];
    }
    return true;
}

// Returns the direction at a the intersection of 2 or more objects
void intersectionDirection(Datas &data,
                           int const Feature,
                           Real const radius,
                           Real const eigenvalueRatioThreshold)
{
    int const noPoints = data.neighbors.size()/NO_DIM;
    int const noIter = 8;
    Real fraction[noIter] = { 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2 };
    Real *pos = &( data.neighbors[data.objectIndex*NO_DIM] );
    int offset = Feature==3 ? 6 : 0;
    
    for (int i=0; i<noIter; ++i)
    {
        Real tempDis = radius * fraction[i];
        vector<Real> tempPos;
        
        for (int j=0; j<noPoints; ++j)
        {
            if ( distance(pos,&(data.neighbors[j*NO_DIM]))<=tempDis )
                for (int k=0; k<NO_DIM; ++k)
                    tempPos.push_back( data.neighbors[j*NO_DIM+k] );
        }
        if (tempPos.size()/NO_DIM<5)
            break;
        
        // get the eigenvalues
        Real eigenvalues[NO_DIM], eigenvectors[NO_DIM*NO_DIM];
        cloudDirections( tempPos, eigenvalues, eigenvectors );
        Real eigenvalueRatio = Feature==3 ? eigenvalues[1]/eigenvalues[2] : eigenvalues[0]/eigenvalues[1];
        if ( eigenvalueRatio<=eigenvalueRatioThreshold )    //valid direction found
        {
            for (int j=0; j<NO_DIM; ++j)
                data.spineDirection[j] = eigenvectors[offset+j];
            break;
        }
    }
}


// this function moves a point according to a given direction and also checks if the point reached convergence
inline bool movePoint(Datas &data,
                      int const Feature,
                      Real const stepFraction,
                      Real const distanceThreshold,
                      Real const eigenvalueRatioThreshold)
{
    int const noPoints = data.neighbors.size() / NO_DIM;
    Real *pos = data.spinePosition;
    Real *directionVector = data.spineDirection;
    
    //compute the center of where the particles are in the given radius
    Real massCenter[] = {0.,0.,0.}; //will keep track of the mass center
    for (int i=0; i<noPoints; ++i)
        for (int j=0; j<NO_DIM; ++j)
            massCenter[j] += data.neighbors[i*NO_DIM+j];
    for (int j=0; j<NO_DIM; ++j)
        massCenter[j] /= noPoints;
    
    Real moveDirection[NO_DIM], distance = 0.;
    if (Feature==3)
        distance = distanceAlongPerpendicularDirection( pos, massCenter, directionVector, moveDirection );
    else if (Feature==2)
        distance = distanceAlongDirection( pos, massCenter, directionVector, moveDirection );
    for (int j=0; j<NO_DIM; ++j)
        if ( boost::math::isfinite(moveDirection[j]) )
            pos[j] += stepFraction * moveDirection[j];
    
    if (distance<distanceThreshold) // check if the point converged to its final position
    {
        Real eigenvalues[NO_DIM], eigenvectors[NO_DIM*NO_DIM];
        if (not cloudDirections( data.neighbors, massCenter, eigenvalues, eigenvectors) )
            return false;
        Real eigenvalueRatio = Feature==3 ? eigenvalues[1]/eigenvalues[2] : eigenvalues[0]/eigenvalues[1];
        if (eigenvalueRatio<eigenvalueRatioThreshold)
            return true;
    }
    return false;
}
inline bool movePoint2(Datas &data,
                       int const Feature,
                       Real const stepFraction,
                       Real const distanceThreshold,
                       Real const eigenvalueRatioThreshold)
{
    int const noPoints = data.neighbors.size() / NO_DIM;
    Real *pos = data.spinePosition;
    Real *directionVector = data.spineDirection;
    
    // compute the center of where the particles are in the given radius
    Real massCenter[] = {0.,0.,0.}; //will keep track of the mass center
    for (int i=0; i<noPoints; ++i)
        for (int j=0; j<NO_DIM; ++j)
            massCenter[j] += data.neighbors[i*NO_DIM+j];
    for (int j=0; j<NO_DIM; ++j)
        massCenter[j] /= noPoints;
    
    // get the direction associated to the point cloud
    Real eigenvalues[NO_DIM], eigenvectors[NO_DIM*NO_DIM];
    if (not cloudDirections( data.neighbors, massCenter, eigenvalues, eigenvectors) )
        return false;
    int offset = Feature==3 ? 6 : 0;
    for (int j=0; j<NO_DIM; ++j)
            directionVector[j] = eigenvectors[offset+j];
    
    // move the point perpendicular to the main direction of the point cloud
    Real moveDirection[NO_DIM], distance = Real(0.);
    if (Feature==3)
        distance = distanceAlongPerpendicularDirection( pos, massCenter, directionVector, moveDirection );
    else if (Feature==2)
        distance = distanceAlongDirection( pos, massCenter, directionVector, moveDirection );
    for (int j=0; j<NO_DIM; ++j)
        if ( boost::math::isfinite(moveDirection[j]) )
            pos[j] += stepFraction * moveDirection[j];
    
    // check if point reached convergence criterion
    Real eigenvalueRatio = Feature==3 ? eigenvalues[1]/eigenvalues[2] : eigenvalues[0]/eigenvalues[1];
    if ( distance<distanceThreshold and eigenvalueRatio<eigenvalueRatioThreshold )
        return true;
    return false;
}





/* This function computes the inertia tensor of the cloud of points around a given point. Returns the numbers of cells that changed direction more than a given angle. It takes the following arguments:
 *              spinePosition   = the input point positions
 *              oldDirections   = the input onld directions associated to the points
 *              directions      = array to keep track of the output directions
 *              options         = keep track of the box dimensions, radius, feature and convergence options
 */
size_t directionDetection(Array<Real,2> &spinePosition,
                          Array<Real,2> &oldDirections,
                          Array<Real,2> &directions,
                          ContractionOptions &options,
                          bool const VERBOSE)
{
    ProgressMessage message;
    if (VERBOSE)
    {
        cout << "Finding the feature direction for each valid object voxel.\n\tDone: ";
        //message.updateProgress( 0 );
    }
    // set some variables
    int const Feature = options.feature;
    Real const radius = options.radius;
    Real const eigenvalueRatioThreshold = options.eigenvalueThreshold;
    Real const cosTolerance = options.cosDirectionThreshold;
    size_t const minimumNoNeigbors = 5;
    
    
    // take into account the periodic boundaries => copy the particles on the edge
    vector<Real> tempSpinePosition;
    vector<Int>  tempIndices;       //vector that will keep track of the indices of the extra points added due to periodic boundary conditions
    Real *data = &( spinePosition(0,0) );
    size_t const noBaseCells = spinePosition.axisSize(0);
    size_t noTotalCells = noBaseCells;
    if ( options.periodic ) // if periodic boundaries - add padding particles
    {
        noTotalCells = periodicBoundary( spinePosition, radius, options.box, tempSpinePosition, &tempIndices, false );
        data = &( tempSpinePosition[0] );
    }
    
    
    // build the tree used to search the data
    std::vector< BSR<Real> > result;
    BinarySearch<Real,NO_DIM> tree( options.box, radius, options.periodic );
    tree.buildTree( data, noTotalCells );
    
    
    // loop over all the cells
    size_t noDifferentDirections = 0, noNoDirection = 0;
#ifdef ENABLE_OPENMP
    int const numThreads = omp_get_max_threads();
    #pragma omp parallel for private(result) reduction(+ :noDifferentDirections) reduction(+ :noNoDirection)
#endif
    for (unsigned long long object=0; object<noBaseCells; ++object)
    {
        tree.pointNeighbors_radius( &(data[object*NO_DIM]), &result );
        if ( result.empty() or result.size()<minimumNoNeigbors )
        {
            for (int i=0; i<NO_DIM; ++i)
                directions(object,i) = oldDirections(object,i);
            if ( magnitude( &(directions(object,0)) )<Real(0.5) )
                directions(object,0) = 1.;
            ++noNoDirection;
            continue;
        }
        
        // copy the neighbors data for further computations
        Datas dataProxy;                            // structure to keep track of the data that need to be transfered between the different functions
        dataProxy.neighbors.clear();
        dataProxy.directions.clear();
        dataProxy.spinePosition = &(spinePosition(object,0));
        dataProxy.spineDirection = &(directions(object,0));
        dataProxy.updateDirection = true;
        dataProxy.objectIndex = -1;
        
        for (size_t i=0; i<result.size(); ++i)
        {
            // get the neighbors of the point
            int index = result[i].id;
            for (int j=0; j<NO_DIM; ++j)
                dataProxy.neighbors.push_back( data[index*NO_DIM+j] );
            // keep track of which particle is the particle in question
            if (index==object)
                dataProxy.objectIndex = i;
        }
        
        // compute the eigenvalues and eigenvectors for the given point cloud
        Real eigenvalues[NO_DIM], eigenvectors[NO_DIM*NO_DIM];
        cloudDirections( dataProxy.neighbors, eigenvalues, eigenvectors );
        
        Real eigenvalueRatio = Feature==3 ? eigenvalues[1]/eigenvalues[2] : eigenvalues[0]/eigenvalues[1];
        int offset = Feature==3 ? 6 : 0;
        for (int j=0; j<NO_DIM; ++j)
            dataProxy.spineDirection[j] = eigenvectors[offset+j];
        if (eigenvalueRatio>eigenvalueRatioThreshold)   //if intersection of filaments or walls
            intersectionDirection( dataProxy, Feature, radius, eigenvalueRatioThreshold );
        Real tempM = magnitude( dataProxy.spineDirection ); //magnitude
        Real factor = dataProxy.spineDirection[0]>0. ? 1. : -1.;   //choose the same direction for all eigenvectors
        for (int j=0; j<NO_DIM; ++j)
            dataProxy.spineDirection[j] *= factor/tempM;    //normalize the direction to 1 and give the same pointing along positive x-direction
        if ( tempM<1.e-5 or not boost::math::isfinite(tempM) )  //keep track if no eigenvector is set
            for (int j=0; j<NO_DIM; ++j)
                dataProxy.spineDirection[j] = 1./sqrt(3.); //set some random direction
        
        // compare the new direction with the old one
        if ( fabs( scalarProduct(&(oldDirections(object,0)),dataProxy.spineDirection) )<cosTolerance )
            ++noDifferentDirections;
        for (int j=0; j<NO_DIM; ++j)
            oldDirections(object,j) = directions(object,j);
    }
    cout << "   Not converged: " << noDifferentDirections << " - " << std::setprecision(3) << 100.*noDifferentDirections/noBaseCells << "\%.  For " << noNoDirection << " points could not find a direction.\n";
    return noDifferentDirections;
}






/* Contracts the filaments/walls to their spine/principal wall. It takes the following arguments:
 *          spinePosition   = the initial position of each voxel that is part of the filament/wall - on exit will return the contracted spine/wall position for that point
 *          directionVector = the direction of the filament / normal to the wall for each of the above points
 *          options         = keep track of the box dimensions, radius, feature and convergence options
 * NOTE: For a point to be considered converged it needs to satisfy all 3 threshold criteria: distanceThreshold, eigenvalueThreshold and cosDirectionThreshold.
 *  */
void MMFContraction(Array<Real,2> &spinePosition,
                    Array<Real,2> &directionVector,
                    ContractionOptions options,
                    Real const convergeFraction,
                    Real const distanceThreshold,
                    Real const eigenvalueThreshold,
                    int const maxLoop)
{
    cout << "Contracting the Cosmic Web objects to a spine for filaments and plane for walls:\n" << flush;
    int const noBaseCells = spinePosition.axisSize(0);
    
    
    // define variable to keep track of the convergence procedure
    Array<bool,1> convergedPoints(noBaseCells); // true if point has converged to final position and will not be displaced anymore
    convergedPoints = false;
    Array<bool,1> hasCriteria(noBaseCells);     //true if point has satisfied all convergence criteria except that all its neighbors should have converged
    hasCriteria = false;
    
    
    // the loop that contracts the structures
    Real const stepFraction = .99;          //what fraction of the distance to move every step
    bool continueLoop = true;
    int loopNo = 0, noNotConverged = 0;
    Real notConvergedFraction = 1.;
    while (continueLoop and loopNo<maxLoop)
    {
        // take into account the periodic boundaries => copy the particles on the edge
        vector<Real> tempSpinePosition;
        vector<Int>  tempIndices;       //vector that will keep track of the indices of the extra points added due to periodic boundary conditions
        Real *data = &( spinePosition(0,0) );
        size_t noTotalCells = noBaseCells;
        if ( options.periodic ) // if periodic boundaries - add padding particles
        {
            noTotalCells = periodicBoundary( spinePosition, options.radius, options.box, tempSpinePosition, &tempIndices, false );
            data = &( tempSpinePosition[0] );
        }
        
        // construct the tree for fast searching
        std::vector< BSR<Real> > result;
        BinarySearch<Real,NO_DIM> tree( options.box, options.radius, options.periodic );
        tree.buildTree( data, noTotalCells );
        
        
        // search all the neighbors within 'radius' and move the point accordingly
        noNotConverged = 0;
#ifdef ENABLE_OPENMP
        #pragma omp parallel for private(result) reduction(+ :noNotConverged)
#endif
        for (int object=0; object<noBaseCells; ++object)
        {
            if ( convergedPoints(object) ) continue;
            
            tree.pointNeighbors_radius( &(data[object*NO_DIM]), &result );
            if ( result.empty() ) continue;
            
            // copy the neighbors data for further computations
            Datas dataProxy;                            // structure to keep track of the data that need to be transfered between the different functions
            dataProxy.neighbors.clear();
            dataProxy.directions.clear();
            dataProxy.spinePosition = &(spinePosition(object,0));
            dataProxy.spineDirection = &(directionVector(object,0));
            dataProxy.updateDirection = false;
            dataProxy.objectIndex = -1;
            
            bool neighborsConverged = true;
            for (size_t i=0; i<result.size(); ++i)
            {
                // get the neighbors of the point
                int index = result[i].id;
                for (int j=0; j<NO_DIM; ++j)
                    dataProxy.neighbors.push_back( data[index*NO_DIM+j] );
                // find if all the neighbors have converged
                if ( index<noBaseCells and not hasCriteria(index) )
                    neighborsConverged = false;
            }
            
            // move the points acording to the point gradient
            hasCriteria(object) = movePoint( dataProxy, options.feature, stepFraction, distanceThreshold, eigenvalueThreshold );
            if ( neighborsConverged and hasCriteria(object) )
                convergedPoints(object) = true;
            else
                ++noNotConverged;
        }
        Real tempNotConverged = 1.*noNotConverged/noBaseCells;
        if ( notConvergedFraction-tempNotConverged < convergeFraction and loopNo>=5 )
            continueLoop = false;
        notConvergedFraction = tempNotConverged;
        if ( options.periodic ) // fold particles outside the box boundary
            periodicFolding(spinePosition, options.box, false);
        ++loopNo;
    }
    cout << "\tDid " << loopNo << " iterations out of a maximum of " << maxLoop << ".   Fraction not converged: "  << std::setprecision(3) << noNotConverged*100./noBaseCells << "\%\n";
}
void MMFContraction(Array<Real,2> &spinePosition,
                    ContractionOptions options,
                    Real const convergeFraction,
                    Real const distanceThreshold,
                    Real const eigenvalueThreshold,
                    int const maxLoop)
{
    cout << "Contracting the Cosmic Web objects:\n" << flush;
    int const noBaseCells = spinePosition.axisSize(0);
    ProgressMessage message;
    
    
    // define variable to keep track of the convergence procedure
    Array<bool,1> convergedPoints(noBaseCells); // true if point has converged to final position and will not be displaced anymore
    convergedPoints = false;
    Array<bool,1> hasCriteria(noBaseCells);     //true if point has satisfied all convergence criteria except that all its neighbors should have converged
    hasCriteria = false;
    
    
    // the loop that contracts the structures
    Real const stepFraction = .99;          //what fraction of the distance to move every step
    bool continueLoop = true;
    int loopNo = 0, noNotConverged = 0;
    Real notConvergedFraction = 1.;
    while (continueLoop and loopNo<maxLoop)
    {
        // take into account the periodic boundaries => copy the particles on the edge
        vector<Real> tempSpinePosition;
        vector<Int>  tempIndices;       //vector that will keep track of the indices of the extra points added due to periodic boundary conditions
        Real *data = &( spinePosition(0,0) );
        size_t noTotalCells = noBaseCells;
        if ( options.periodic ) // if periodic boundaries - add padding particles
        {
            noTotalCells = periodicBoundary( spinePosition, options.radius, options.box, tempSpinePosition, &tempIndices, false );
            data = &( tempSpinePosition[0] );
        }
        
        // construct the tree for fast searching
        std::vector< BSR<Real> > result;
        BinarySearch<Real,NO_DIM> tree( options.box, options.radius, options.periodic );
        tree.buildTree( data, noTotalCells );
        
        
        // search all the neighbors within 'radius' and move the point accordingly
        noNotConverged = 0;
#ifdef ENABLE_OPENMP
        #pragma omp parallel for private(result) reduction(+ :noNotConverged)
#endif
        for (int object=0; object<noBaseCells; ++object)
        {
            if ( convergedPoints(object) ) continue;
            
            tree.pointNeighbors_radius( &(data[object*NO_DIM]), &result );
            if ( result.empty() ) continue;
            
            // copy the neighbors data for further computations
            Real directionVector[NO_DIM];
            Datas dataProxy;                            // structure to keep track of the data that need to be transfered between the different functions
            dataProxy.neighbors.clear();
            dataProxy.directions.clear();
            dataProxy.spinePosition = &(spinePosition(object,0));
            dataProxy.spineDirection = directionVector;
            dataProxy.updateDirection = true;
            dataProxy.objectIndex = -1;
            
            bool neighborsConverged = true;
            for (size_t i=0; i<result.size(); ++i)
            {
                // get the neighbors of the point
                int index = result[i].id;
                for (int j=0; j<NO_DIM; ++j)
                    dataProxy.neighbors.push_back( data[index*NO_DIM+j] );
                // find if all the neighbors have converged
                if ( index<noBaseCells and not hasCriteria(index) )
                    neighborsConverged = false;
                // keep track of which particle is the particle in question
                if (index==object)
                    dataProxy.objectIndex = i;
            }
            
            // move the points acording to the point gradient and update the direction of the filament/wall
            hasCriteria(object) = movePoint2( dataProxy, options.feature, stepFraction, distanceThreshold, eigenvalueThreshold );
            if ( neighborsConverged and hasCriteria(object) )
                convergedPoints(object) = true;
            else
                ++noNotConverged;
        }
        Real tempNotConverged = 1.*noNotConverged/noBaseCells;
        if ( notConvergedFraction-tempNotConverged < convergeFraction and loopNo>=5 )
            continueLoop = false;
        notConvergedFraction = tempNotConverged;
        if ( options.periodic ) // fold particles outside the box boundary
            periodicFolding(spinePosition, options.box, false);
        ++loopNo;
        cout << loopNo << "   Fraction not converged: "  << std::setprecision(3) << noNotConverged*100./noBaseCells << "\%\n";
    }
    cout << "\tDid " << loopNo << " iterations out of a maximum of " << maxLoop << ".     Fraction not converged: "  << std::setprecision(3) << noNotConverged*100./noBaseCells << "\%\n";
}


/* Uses an iterative method to get the direction of the filament/wall.
It finds the neighbors of the point within a given radius 'radius' and than it moves it to that center of mass following the point density gradient. It repeats this procedure until all points lie in pronunced filaments/wall. Once it does so for all the points, computes the characteristic direction of the resulting point set.
* In the next step it repeats the above iterative procedure but starting with a gues for the filament/wall direction. The gues for the filament/wall direction is refined after each contraction.*/
void MMFObjectDirection(Array<Real,2> &spinePosition,
                        Array<Real,2> &directionVector,
                        ContractionOptions &options,
                        bool hasInitialDirection)
{
    cout << "Contracting the MMF objects to a spine for filaments and plane for walls ...\n" << flush;
    int const noBaseCells = spinePosition.axisSize(0);
    
    
    // first contract the points along the volume center direction
    cout << "\tContracting the objects using the volume center in the imediate neighborhood:\n" << flush;
    Array<Real,2> initialSpinePosition(noBaseCells,NO_DIM);
    initialSpinePosition = spinePosition;
    if ( not hasInitialDirection )
        MMFContraction( spinePosition, options, options.convergeFraction, options.radius, options.eigenvalueThreshold, options.maxLoop );
    
    
    // compute the eigenvector of the cloud points and the direction associated to them
    Array<Real,2> tempDirections(noBaseCells,NO_DIM);
    tempDirections.assign( Real(0.) );
    if ( not hasInitialDirection )
        directionDetection( spinePosition, tempDirections, directionVector, options, true );
    
    
    // now iterate until reaching convergence in the filament/wall direction 
    cout << "\tIteration computation until achieving convergence :\n" << flush;
    Real notConvergedFraction = 1.;
    for (int i=0; i<options.maxIterations; ++i)
    {
        cout << "ITERATION:  <<<  " << i+1 << "  >>>\n";
        spinePosition = initialSpinePosition;
        tempDirections = directionVector;
        
        MMFContraction( spinePosition, directionVector, options, options.convergeFraction, options.distanceThreshold, options.eigenvalueThreshold, options.maxLoop );
        
        size_t const notConverged = directionDetection( spinePosition, tempDirections, directionVector, options, true );
        Real tempNotConverged = 1.*notConverged/noBaseCells;
        if ( notConvergedFraction-tempNotConverged < options.convergeFraction and i>=3 )
            break;
        notConvergedFraction = tempNotConverged;
    }
}






// searches to which bin corresponds the current values
template <typename T>
inline bool binIndex(Array<T,2> &bins, T position, int *index)
{
    int last = bins.axisSize(0) - 1;
    if ( position<bins(0,0) )
        return false;
    else if ( position>=bins(last,0) )
        return false;
    else
    {
        for (int i=0; i<last; ++i)
            if ( position<bins(i+1,0) and position>=bins(i,0) )
            {
                *index = i;
                return true;
            }
    }
    return false;
}

// Takes as input a value that corresponds to a given interval. It finds the interval (binMin, binMax) in the given search array bins and than distributes the input value in equal weight in the (binMin, binMax) interval.
template <typename T>
inline void binIndex_distribute(Array<T,2> &bins, T binMin, T binMax, T value, int jIndex)
{
    int last = bins.axisSize(0) - 1;
    if ( binMin>=binMax ) throwError( "In function 'binIndex_distribute'. The parameter 'binMin' must be smaller than 'binMax'." );
    if ( binMax<bins(0,0) or binMin>bins(last,0) )
        return;
    
    //search the bin index where binMin and binMax are located
    int indexMin = 0, indexMax = last;
    for (int i=0; i<last; ++i)
    {
        if ( binMin<bins(i+1,0) and binMin>=bins(i,0) )
            indexMin = i;
        if ( binMax<bins(i+1,0) and binMax>=bins(i,0) )
        {
            indexMax = i;
            break;
        }
    }
    
    // if both the (binMin,binMax) interval is in a single bin
    if ( indexMin==indexMax )
    {
        bins(indexMin,jIndex) += value;
        return;
    }
    
    // distribute the value within the given bins when the (binMin,binMax) interval is in multiple bins
    T temp = value / (binMax - binMin);
    T weight = binMin>bins(indexMin,0) ? (bins(indexMin+1,0)-binMin) : (bins(indexMin+1,0)-bins(indexMin,0)) ;
    bins(indexMin,jIndex) += temp * weight;
    
    for (int i=indexMin+1; i<indexMax; ++i)
        bins(i,jIndex) += temp * ( bins(i+1,0)-bins(i,0) );
    
    weight = binMax<bins(indexMax+1,0) ? (binMax-bins(indexMax,0)) : (bins(indexMax+1,0)-bins(indexMax,0)) ;
    bins(indexMax,jIndex) += temp * weight;
}





// function that computes different quantities along the filament and plane of the wall
void environmentProperties(Array<Real,1> &density,
                           Array<Real,2> &spinePosition,
                           ContractionOptions &options,
                           Array<Real,2> &properties,
                           Array<double,2> &sizeData,
                           Array<double,2> &massData,
                           Array<Real,1> &averageData)
{
    // define some variables
    int const feature = options.feature;
    Real const radius = options.radius;
    int const noBaseCells = spinePosition.axisSize(0);
    Real dx[NO_DIM];
    for (int i=0; i<NO_DIM; ++i) dx[i] = options.length[i] / options.grid[i];
    
    cout << "Computing properties (e.g. thickness, mass density) of the NEXUS filaments and walls smoothed on a " << radius << " Mpc/h scale. There are " << noBaseCells << " valid cells (" << std::setprecision(3) << (noBaseCells*100.)/(options.grid[0]*options.grid[1]*options.grid[2]) << "\% of the total volume):\n" << flush;
    
    
    // take into account the periodic boundaries => copy the particles on the edge
    vector<Real> tempSpinePosition;
    vector<Int>  tempIndices;       //vector that will keep track of the indices of the extra points added due to periodic boundary conditions
    Real *data = &( spinePosition(0,0) );
    size_t noTotalCells = noBaseCells;
    if ( options.periodic ) // if periodic boundaries - add padding particles
    {
        noTotalCells = periodicBoundary( spinePosition, options.radius, options.box, tempSpinePosition, &tempIndices, false );
        data = &( tempSpinePosition[0] );
    }
    
    // construct the tree for fast searching
    std::vector< BSR<Real> > result;
    BinarySearch<Real,NO_DIM> tree( options.box, options.radius, options.periodic );
    tree.buildTree( data, noTotalCells );
    
    
    // loop over the voxels which are part of the object
    Real R2 = radius*radius;
    Real factorF = 4./3.14 * dx[0]*dx[1]*dx[2] / radius;    // multiplication factor for filament diameter computation
    Real factorW = 1./3.14 * dx[0]*dx[1]*dx[2] / R2;        // multiplication factor for wall thickness computation
    int minNeighbors = feature==3 ? int(radius/dx[0]+1.) :  int(R2/dx[0]/dx[0]+1.);    //minimum number of neighbors for valid cells
    size_t badCells = 0, badCellsVolume = 0;
    size_t averageVolume = 0;
    double totalLength = 0., averageDiameter = 0.;    // for filaments: totalLength=total length of filaments, averageDiameter=average diameter
    //for walls: totalLength = total surface area, averageDiameter = average thickness for walls
    cout << "Discarding regions with less than " << minNeighbors << " neighbors within the search radius\n" << flush;
    for (int object=0; object<noBaseCells; ++object)
    {
        tree.pointNeighbors_radius( &(data[object*NO_DIM]), &result );
        size_t volume = result.size();
        if ( volume<minNeighbors )
        {
            ++badCells;
            badCellsVolume += volume;
            continue;
        }
        
        averageVolume += volume;
        Real mass = 0.;
        //computes the mass in the given smoothing length
        for (size_t j=0; j<volume; ++j)
            if ( result[j].id<noBaseCells )
                mass += density[ result[j].id ];
            else
                mass += density[ tempIndices[result[j].id] ];
        
        
        Real len = 0.;              // filament length / wall area assigned to this point
        Real diameter[2] = {0.,0.}; // stores the filament diameter / wall thickness
        Real linDensity = 0.;       // stores the linear mass density for filaments / surface mass density for walls
        
        if (feature==3) //for filaments
        {
            len = radius / volume;             // length of the filament associated to the given point
            totalLength += len;
            diameter[0] = sqrt( (volume-.5)*factorF ); // filament diameter - minimum value
            diameter[1] = sqrt( (volume+.5)*factorF ); // filament diameter - maximum value
            linDensity = mass / radius;        // linear mass density of the filament
        }
        else if (feature==2) //for walls
        {
            len = 3.14*R2 / volume;         // area of the wall associated to the given point
            totalLength += len;
            diameter[0] = (volume-.5) * factorW;    // wall thickness - minimum value
            diameter[1] = (volume+.5) * factorW;    // wall thickness - maximum value
            linDensity = mass / (3.14*R2);  // surface mass density of the wall
        }
        averageDiameter += (diameter[0] + diameter[1]) * len/2.;
        
        int index = 0;
        binIndex_distribute( sizeData, double(diameter[0]), double(diameter[1]), double(len), 1 );
        binIndex_distribute( sizeData, double(diameter[0]), double(diameter[1]), double( density[object] ), 3 );
        binIndex_distribute( sizeData, double(diameter[0]), double(diameter[1]), double(1.), 4 );
        if ( binIndex(massData,double(linDensity),&index) ) massData(index,1) += len;
        
        properties(object,0) = (diameter[0] + diameter[1])/2.;    // filament diameter / wall thickness
        properties(object,1) = linDensity;                        // filament linear mass density / wall surface mass density
    }
    
    for (int i=0; i<sizeData.axisSize(0); ++i)
        sizeData(i,2) = sizeData(i,1) / totalLength;
    for (int i=0; i<massData.axisSize(0); ++i)
        massData(i,2) = massData(i,1) / totalLength;
    
    // compute the cumulative sums for the mass and volume distribution
    double totalMass = 0., totalVolume = 0.;
    for (Int i=0; i<sizeData.axisSize(0); ++i)
    {
        totalMass += sizeData(i,3);
        totalVolume += sizeData(i,4);
        sizeData(i,3) = totalMass;
        sizeData(i,4) = totalVolume;
    }
    for (Int i=0; i<sizeData.axisSize(0); ++i)
    {
        sizeData(i,3) /= totalMass;
        sizeData(i,4) /= totalVolume;
    }
    
    
    averageData(0) = averageDiameter / totalLength;      // average diameter
    averageData(1) = totalLength/(options.length[0]*options.length[1]*options.length[2]); //average length per unit volume
    
    
    cout << "\t Average volume per cell : " << std::setprecision(4) << 1.*averageVolume/noBaseCells << " cells \n";
    cout << "\t Average " << (feature==3 ? "filament diameter" : "wall thickness") << " :  " << std::setprecision(4) << averageData(0) << "  Mpch\n";
    cout << "\t Total " << (feature==3 ? "filament length" : "wall area") << " :  "<< std::setprecision(4) << averageData(1)*1.e6 << (feature==3 ? "  Mpch" : "  Mpch^2") << " / (100 Mpch)^3\n";
    if (badCells!=0) cout << "\t Number of cells below the limit : " << badCells << "(" << std::setprecision(4) << float(badCells)/noBaseCells*100. << "\% of valid volume) with average volume: " << std::setprecision(4) << float(badCellsVolume)/badCells << "\n" << flush;
}


