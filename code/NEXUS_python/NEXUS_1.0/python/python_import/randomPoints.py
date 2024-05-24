
import numpy as np
import math
from scipy.weave import inline, converters
from miscellaneous import throwError, throwWarning


def randomPointsOnSphereSurface(N):
    """ Generates a numpy array of random points on a sphere's surface - sphere of radius 1. """
    sph = np.empty( (N,3), np.float )
    
    sph[:,2] = np.random.uniform(-1.0,1.0,N) # z-coordinates
    z2 = np.sqrt(1.0 - sph[:,2]**2)
    phi = (2.0 * math.pi) * np.random.random( N )
    sph[:,0] = z2 * np.cos(phi) # x 
    sph[:,1] = z2 * np.sin(phi) # y
    
    return sph


def randomPointsInSphere(N):
    """ Generates a numpy array of random points inside a sphere of radius 1.  """
    sph = np.empty( (N,3), np.float )
    
    cube = np.random.uniform( -1.0, 1.0, (10*N,3) )
    count = 0
    for i in range(3*N):
        if not (cube[i,0]**2 + cube[i,1]**2 + cube[i,2]**2)>1.:
            sph[count,:] = cube[i,:]
            ++count
        if count >= N: break
    
    return sph




def uniformPointsOnSphereSurface(N):
    """ Generates a numpy array of uniformly distributed points on a sphere's surface - sphere of radius 1. """
    uniformPointsSuportCode = """
#line 1000 "uniformPointsSuportCode"
#include <iostream>  // sync_with_stdio
#include <cstdio>    // printf
#include <cstdlib>   // I believe sqrt is in here
#include <ctime>     // Used to salt the random number generator
#include <cmath>     // Various mathematical needs
#include <vector>    // Used to contain many points
#include <valarray>  // Used to implement a true XYZ vector
using namespace std; // Necessary to gain access to many C++ names

//******************************************************************************
typedef valarray<double>coordinates; // To simplify declarations
//******************************************************************************
class XYZ { // This class contains operations for placing and moving points
  public:
    double change_magnitude;
    coordinates xyz;   // This holds the coordinates of the point.
    coordinates dxyz;  // This holds the summed force vectors from other points
    inline double random() { return(double((rand()%1000)-500)); } // ordinates
    inline double square(const double& n) { return(n*n); }
    inline coordinates square(const coordinates& n) { return(n*n); }
    inline double inverse(const double& n) { return(1.0/n); }
    XYZ& inverse_square() { xyz*=inverse(square(magnitude())); return *this; }
    inline double magnitude() { return(sqrt((xyz*xyz).sum())); }
    void normalize() { xyz/=magnitude(); } // unit vector
    XYZ(): xyz(3), dxyz(3) {
      xyz[0]=random(); xyz[1]=random(); xyz[2]=random(); normalize();
    }
    XYZ(const double& x,const double& y,const double& z) : xyz(3), dxyz(3) {
      xyz[0]=x; xyz[1]=y; xyz[2]=z;
    }
    XYZ(const coordinates& p) : xyz(3), dxyz(3) {
      xyz=p;
    }
    ~XYZ() { }
    coordinates& array() { return xyz; }
    void zero_force() { dxyz=0.0; }
    double change() { return(change_magnitude); }
    double magnitude(XYZ& b) { // Return length of vector.  (not const)
      return(sqrt( square(b.array()-xyz).sum() ));
    }
    void sum_force(XYZ& b) { // Cause force from each point to sum.  (not const)
      dxyz+=(XYZ(b.array()-xyz).inverse_square().array()); // Calculate and add
    }
    void move_over_sphere() { // Cause point to move due to force
      coordinates before=xyz;                       // Save previous position
      xyz-=dxyz;                                    // Follow movement vector
      normalize();                                  // Project back to sphere
      before-=xyz;                                  // Calculate traversal
      change_magnitude=sqrt((before*before).sum()); // Record largest
    }
};

//******************************************************************************
class points { // This class organizes expression of relations between points
  public:
    const size_t N;   // Number of point charges on surface of sphere
    const size_t R;   // Number of rounds after which to stop
    const double L;   // Threshold of movement below which to stop
    char        *S;   // Name of this vertex set
    size_t rounds;    // Index of rounds processed
    vector<XYZ>V;     // List of point charges
    vector<double>H;  // List of minimum distances
    double maximum_change; // The distance traversed by the most moved point
    double minimum_radius; // The radius of the smallest circle
    time_t T0;        // Timing values

    void relax() { // Cause all points to sum forces from all other points
      size_t i, j;
      rounds=0;
      do {
        maximum_change=0.0;
        for(i=1;i<N;i++) {   // for all points other than the fixed point
	  V[i].zero_force();                       // Initialize force vector
	  for(j=  0;j<i;j++) V[i].sum_force(V[j]); // Get contributing forces
	  // Skip i==j
	  for(j=i+1;j<N;j++) V[i].sum_force(V[j]); // Get contributing forces
	}
        for(i=1;i<N;i++) {  // React to summed forces except for the fixed point
	  V[i].move_over_sphere();
	  if(V[i].change()>maximum_change) maximum_change=V[i].change();
	}
      } while(maximum_change>L&&++rounds<R); // Until small or too much movement
    }
    points(char *s,const size_t& n,const double& l,const size_t& r) :
      N(n), L(l), R(r)
    {
      S=s;
      T0=time(0L);                   // Get the current time
      srand(T0);                     // Salt the random number generator.
      V.push_back(XYZ(1.0,0.0,0.0)); // Create Anchored first point V[0] (1,0,0)
      H.push_back(2.0);
      while(V.size()<N) {   // For all other points, until we have enough
	V.push_back(XYZ()); // Create randomized position
        H.push_back(2.0);
	coordinates& last=V.back().array(); // Remember this position
	for(size_t i=V.size()-1;i--;) { // And check to see if it is occupied
	  coordinates& temp=V[i].array();
	  if(temp[0]==last[0]&&temp[1]==last[1]&&temp[2]==last[2]) {
	    V.pop_back(); // Remove the position if it is already occupied
	    break;
	  }
	}
      }
      relax();  // After vector construction, start the relaxation process
      size_t i, j;
      minimum_radius=1.0; // On a unit sphere, the maximum circle radius is 1.0
      for(i=0;i<V.size();i++) { // Discover the minimum distance between points.
	for(j=0;j<V.size();j++) {
	  if(j==i) continue;
	  double rtemp=V[i].magnitude(V[j])/2.0;
	  if(rtemp<minimum_radius) minimum_radius=rtemp; // Record when smaller.
	  if(rtemp<H[i]) H[i]=rtemp;
	}
      }
    }
    ~points() {}
    coordinates& operator[](const size_t& i) { // Caller access to positions
      return(V[i].array());
    }
};
    """
    
    uniformPointsCode = """
    #line 1000 "uniformPointsCode"
size_t Rounds=10000;        // Stop relaxation after this number of rounds.
double Jiggle = 0.03/sqrt(double(NumberPoints));    // Stop relaxation when movement drops below.

points P("default",NumberPoints,Jiggle,Rounds); // compute the unfirom points on the sphere

// write the points in the output array
for (int i=0; i<NumberPoints; ++i)
    for (int j=0; j<3; ++j)
        output(i,j) = P.V[i].xyz[j];
    """
    
    output = np.zeros( (N,3), np.float )
    NumberPoints = N
    inline( uniformPointsCode, ['NumberPoints', 'output'], type_converters=converters.blitz, support_code=uniformPointsSuportCode )
    
    return output