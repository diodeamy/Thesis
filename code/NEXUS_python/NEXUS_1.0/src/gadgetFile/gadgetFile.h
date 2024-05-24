#include <string>
#include "gadget_header.h"

using namespace std;




//! Functions in "readGadgetFile.cc"
void readGadgetHeader(string const &inputFileName,
                      Gadget_header *gadgetHeader);

void readGadgetData(Particle p[],
                    int noParticles,
                    Gadget_header & gadgetHeader,
                    string const & inputFileName);

void readGadgetData(Particle_pvm p[],
                    int noParticles,
                    Gadget_header & gadgetHeader,
                    string const & inputFileName);

void readGadgetData(Particle_pm p[],
                    int noParticles,
                    Gadget_header & gadgetHeader,
                    string const & inputFileName);

void readGadgetData(Particle_p p[],
                    int noParticles,
                    Gadget_header & gadgetHeader,
                    string const & inputFileName);


//!  Functions in "writeGadgetFile.cc"
void writeGadgetFile(string &outputFileName,
		     Gadget_header &gadgetHeader,
                     Particle p[]);
void writeGadgetFile_DM(string &outputFileName,
		        Gadget_header gadgetHeader,
                        Particle p[]);

void writeGadgetFile(string &outputFileName,
                     Gadget_header &gadgetHeader,
                     Particle_pm p[]);
void writeGadgetFile_DM(string &outputFileName,
                        Gadget_header gadgetHeader,
                        Particle_pm p[]);





