#ifndef MMF_HEADER
#define MMF_HEADER

#include <string>
#include <boost/assert.hpp>

#include <defines.h>
#include <density_header.h>


static int const fillSize2 = 1024- 16*8- 18*8- 8;

// types of binary MMF files (each one contains a different result)
static int const MMF_RESPONSE = 1;                          // contains the response value at a given scale
static int const MMF_EIGEN = 5;                             // contains the eigenvalues and eigenvectors for a given scale
static int const MMF_EIGENVECTOR = 6;                       // contains the significant eigenvectors for a given features (eigenvector 3 for filaments -the direction along the filament- and eigenvector 1 for wall -the direction perpendicular to the wall-)
static int const MMF_MAX_RESPONSE = 10;                     // contains the values of the maximum MMF response
static int const MMF_MAX_RESPONSE_SCALE = 11;             // contains the values of the maximum MMF response and also the scale corresponding to the maximum response
static int const MMF_MAX_EIGEN = 15;                        // the eigenvalues and eigenvectors corresponding to the maximum response scale
static int const MMF_MAX_EIGENVECTOR = 16;                  // same as 6, but corresponding to the maximum response scale
static int const MMF_CLEAN_RESPONSE = 20;                   // contains a short int data with values of 0 or 1, depending if the pixel is a valid feature or not (e.g. for filaments: 1=pixel is part of a filament, 0=pixel is not part of the filament)
static int const MMF_CLEAN_RESPONSE_COMBINED = 21;          // contains the above information, but for all the environments: 0 = void, 2 = wall, 3 = filament and 4 = node
static int const MMF_OBJECTS = 30;                          // all the pixels corresponding to the same object have the same value, the object tag/id
static int const MMF_DIRECTIONS = 40;                       // gives the geometrical direction of the filament/wall ONLY at valid pixels positions
static int const MMF_PROPERTIES = 50;                       // gives the properties of filaments / environments using teh contraction procedure

// methods used to compute the MMF response
static int const MMF_METHOD_DENSITY = 1;                    // MMF environments computed from the density field
static int const MMF_METHOD_LOGARITHM_FILTERING = 100;      // MMF environments computed by gaussian smoothing of the density field in logarithmic space
static int const MMF_METHOD_DENSITY_LOGARITHM = 5;          // MMF environments computed from the logarithm of the density field
static int const MMF_METHOD_GRAVITY = 10;                   // MMF environments computed from the gravitational potential (potential obtain from the density map using Poisson's equation)
static int const MMF_METHOD_VELOCITY = 20;                  // MMF environments computed from the velocity divergence field
static int const MMF_METHOD_VELOCITY_LOGARITHM = 25;        // MMF environments computed from the logarithm of the velocity divergence field
static int const MMF_METHOD_VELOCITY_POTENTIAL = 30;        // MMF environments computed from the potential of the velocity divergence field


// data structure that stores information in the header of a MMF file output
struct MMF_header
{
    // Information about the MMF process
    size_t   gridSize[3];   // grid size along each dimension
    size_t   totalGrid;     // total number of grid cells
    int      feature;       // to which feature corresponds this result (4=blob, 3=filament, 2=walls)
    int      scale;         // the index of the current scale (if '-1' -combined result of several scales)
    float    radius;        // Gaussian smoothing radius in cell units (if '-1.' -combined result of several scales)
    float    bias;          // the value of the bias factor '\gamma' from Miguel's thesis
    int      filter;        // which filter was used to obtain this result
    int      fileType;      // which data is written to this file: - see the variables defined above under "types of MMF binary files"
    int      noMMF_files;   // the number of MMF files corresponding to this run
    int      MMF_fileGrid[3]; // if the MMF is saved within multiple files, each file stores only a patch of the full grid (the patches are given by a regular grid with this variable giving the dimensions of that grid)
    int      indexMMF_file; // keeps track of this MMF file index compared to the rest for multiple files (this tells the program what region of the density is saved in this file)
    int      method;        // method used to compute the MMF response (1=density, 2=gravity, )
    double   box[6];         // keep track of the box coordinates in which the MMF response was computed (xMin, xMax, yMin, yMax, zMin, zMax)
    
    
    // Information about the initial data file
    size_t   npartTotal[6]; // the total number of particles
    double   mass[6];   // the mass of the particles in the N-body code
    double   time;      // the expansion parameter 'a' of the snapshot file from which we computed the density
    double   redshift;  // the corresponding redshift
    double   BoxSize;   // the box size in kpc
    double   Omega0;    // Omega_matter
    double   OmegaLambda;   // Omega_Lambda
    double   HubbleParam;   // Hubble parameter h (where H=100 h km/s /Mpc )
    
    
    // The following are non-zero only if this file contains the combined result of several scales
    char    fill[fillSize2];// fill to 1024 bytes - left 744 - used to keep track of information on how the file was obtained
    size_t  FILE_ID;       // keep a unique id for this type of file
    
    
    
    MMF_header();
    void print();                                 // print the information contained in the MMF header
    void copyDensityHeader(Density_header const &densityHeader );     // copy the shared information from the density header
    
    void compatible(MMF_header &header);          // check if the two MMF files have the same grid dimensions
    void compatible(size_t *grid);                // check if the input grid has the same dimensions as the current grid size
    void checkScale(int const scaleNo);           // check if the scale of the MMF header is the same as the argument
    Int totalSize();                              // returns the total number of grid cells of the response
    Int dataSize();                               // returns the total number of grid cells of the response
    
    static float filterRadius(float const radius0,
                              float const base,
                              int const scaleNo); // computes the filter radius corresponding to scale 'scale' via 'baseRadius^scale'
    void updateScale(int const scaleNo,
                     float const radius0,
                     float const base);           // updates both 'scale' and 'radius'
    void updateFeature(int const filterType);     // updates both 'feature' and 'filterType'
    void updateFileType(int const file_type);     // updates 'fileType'
    void updateBias(float const biasValue);       // updates 'bias'
    void updateMMF_fileGrid(int const grid[], int const size);    // update 'MMF_fileGrid' and 'noMMF_files'
    void updateMMF_fileGrid(int const nX,
                            int const nY,
                            int const nZ);        // update 'MMF_fileGrid' and 'noMMF_files'
    void updateMMF_fileIndex(int const fileIndex);// update 'indexMMF_file'
    void updateMaximumResponse();                 // update 'scale' = -1 and 'radius'=-1.
    double boxVolume();                           // returns the volume of the box
    double gridCellVolume();                      // returns the volume of the grid cell (voxel)
    bool emptyBox();                              // returns true if there are no box coordinates
    
    
    // Functions that deal with filename for the different outputs of the MMF computation
    std::string responseFilename(std::string const &rootName,
                                 int const scale,
                                 int const fileIndex=1);    //get the name of the MMF response file for a given scale
    std::string maxResponseFilename(std::string const &rootName,
                                    int const fileIndex=1); //get the name of the maximum MMF response file
    std::string eigenFilename(std::string const &rootName,
                              int const scale,
                              int const fileIndex=1);       //get the name of the file storing the eigenvalues and eigenvectors
    std::string maxEigenFilename(std::string const &rootName,
                                 int const fileIndex=1);    //get the name of the file storing the eigenvalues and eigenvectors for the maximum MMF response
    void checkResponseFiles(std::string const &mmfRootName,
                            int scale[],
                            int const noScales);            //check that the MMF files given by the array of scales exist
    void checkEigenFiles(std::string const &mmfRootName,
                         int scale[],
                         int const noScales);               //check that the Hessian eigevalues and eigenvector files given by the array of scales exist
    void checkFiles(std::string const &mmfRootName);        //check if all the files exist for a result saved in multiple files
    
    // Miscellanous functions
    void updateObservations(std::string newObservations);    // update the observation field
    void updateObservations(char **obs, int const size);
    void overwriteObservations(std::string newObservations); // overwrite in the observation field
    void overwriteObservations(char **obs, int const size);
};




#endif

