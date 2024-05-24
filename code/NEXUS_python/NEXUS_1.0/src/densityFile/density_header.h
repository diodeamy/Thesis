#ifndef DENSITY_HEADER
#define DENSITY_HEADER

#include <string>

#ifndef EXTERNAL_COMPILATION
#include <gadget_header.h>
#include <defines.h>

static const size_t noVelComp = NO_DIM;                     // number of independent velocity components
static const size_t noShearComp = (NO_DIM*(NO_DIM+1))/2;    // number of independent velocity shear components (3 for 2D and 6 for 3D)
static const size_t noVortComp = ((NO_DIM-1)*NO_DIM)/2;
#endif



// constants to keep track of how the density was computed
static size_t const DTFE_METHOD = 1;
static size_t const TSC_METHOD = 2;
static size_t const SPH_METHOD = 3;
static size_t const UNKNOW_METHOD = -1;
// constants to keep track what fields the file contains
static int const DENSITY_FILE = 1;
static int const VELOCITY_FILE = 11;
static int const VELOCITY_GRADIENT_FILE = 12;
static int const VELOCITY_DIVERGENCE_FILE = 13;
static int const VELOCITY_SHEAR_FILE = 14;
static int const VELOCITY_VORTICITY_FILE = 15;
static int const VELOCITY_STD_FILE = 16;
static int const SCALAR_FIELD_FILE = 20;
static int const SCALAR_FIELD_GRADIENT_FILE = 21;
static int const UNKNOW_FILE = -1;


static int const fillSize = 1024 - 13*8 - 8*18 - 8*2;

// data structure that stores information in the header of a density file output
struct Density_header
{
    // Information about the density computations
    size_t  gridSize[3];    // the size of the density grid along all 3D directions
    size_t  totalGrid;      // the total size of the grid = gridSize[0]*gridSize[1]*gridSize[2]
    int     fileType;       // keeps track of the field type stored in the file: DENSITY_FILE, VELOCITY_FILE, etc...
    int     noDensityFiles; // the number of density files corresponding to this run
    int     densityFileGrid[3];  // if density is saved within multiple files, each file stores only a patch of the full density (the patches are given by a regular grid with this variable giving the dimensions of that grid)
    int     indexDensityFile;    // keeps track of this density file index compared to the rest for multiple files (this tells the program what region of the density is saved in this file)
    double  box[6];         // keep track of the box coordinates in which the density was computed (xMin, xMax, yMin, yMax, zMin, zMax)
    
    
    // the following are the same as for the gadget file header and are meant to store information about the input file used to compute the density
    size_t   npartTotal[6];// gives the total number of particles in the gadget simulation for which we compute the density profile
    double   mass[6];      // the mass of the particles in the N-body code
    double   time;         // the expansion parameter 'a' of the snapshot file from which we computed the density
    double   redshift;     // the corresponding redshift
    double   BoxSize;      // the box size in kpc
    double   Omega0;       // Omega_matter
    double   OmegaLambda;  // Omega_Lambda
    double   HubbleParam;  // Hubble parameter h (where H=100 h km/s /Mpc )
    
    
    // additional information about files
    size_t  method;        // the method used to compute the density DTFE_METHOD, TSC_METHOD or SPH_METHOD
    char    fill[fillSize];// fill to 1024 bytes - left 760 - used to keep track of information on how the file was obtained
    size_t  FILE_ID;       // keep a unique id for this type of file
    
    
    Density_header();       // constructor - initializes to 0 or to non-assigned value (=-1) 
    void print();           // print the content of the density header to the standard output
    void copyGadgetHeader(Gadget_header const &gadgetHeader);  // copy the shared information from the Gadget header
    void updateGridSize(size_t size[3]);                       // update the 'gridSize' and 'totalGrid' fields
    void updateGridSize(size_t const sizeX,
                        size_t const sizeY,
                        size_t const sizeZ);                   // update the 'gridSize' and 'totalGrid' fields
    size_t dataSize();                                         // returns the number of data points pointed by the header
    double boxVolume();                                        // returns the volume of the box for the given data
    bool emptyBox();                                           // returns true if there are no values set for the box coordinates for the given density data
    double massInCell();                                       // returns the average mass in a density=1 cell
    
    void checkSimilar(Density_header &other);                  // check if 'other' has the same grid dimension as this header
    bool isSimilar(Density_header &other);                     // check if 'other' has the same grid dimension as this header
    
    void updateDensityFilename(std::string densityFileName,
                               size_t const fileGrid[]);
    void noneFileGrid();                                       // don't use multiple files
    void updateObservations(char **obs, int const size);       // update the observation field
    void updateObservations(std::string newObservations);      // update the observation field
    void overwriteObservations(std::string newObservations);   // overwrite in the observation field
    
    // the following are functions used with multiple density files
    std::string filename(std::string &rootName, int const i);  // return the filename of file 'i' when the density is saved in multiple files
    void secondarySize(int const fileNo, size_t *ptr);         // returns the size along each direction of a density file from a density written in multiple files
    size_t secondaryTotalSize(int const fileNo);               // returns the total length of the data to be saved in file 'fileNo' when density is saved in multiple files
    void boundariesIndices(int const fileNo, size_t *ptr);     // returns the indices delimiting the boundaries of file 'fileNo' when density is saved in multiple files
    void boundariesUnits(int const fileNo, Real *ptr);         // returns the boundaries (in units) of file 'fileNo' when density is saved in multiple files
    void boundariesBoxLength(int const fileNo, Real *ptr);     // returns the boundaries (as a fraction of boxLength) of file 'fileNo' when density is saved in multiple files
    void copyDensityToMain(Real *mainArray,
                           Real const *secondaryArray,
                           int const fileNo);                  // copies to the mainArray (which stores the density in the full box) the density corresponding to file 'fileNo' for a density saved in multiple files
    void copyDensityToSecondary(Real *secondaryArray,
                                int const fileNo,
                                Real const *mainArray);        // copies the corresponding density to file 'fileNo' from the main density array
};


#endif

