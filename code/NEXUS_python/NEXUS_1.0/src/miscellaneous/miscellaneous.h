#ifndef MISCELLANEOUS_HEADER
#define MISCELLANEOUS_HEADER


#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>


using std::cout;
using std::string;
using std::vector;


#include "misc.h"





//! Functions in "miscellaneousUserInput.cc"
int getIntegerFromUser();

int getIntegerFromUser(int minValue);

int getIntegerFromUser(int minValue,
               int maxValue);

double getDoubleFromUser();

double getDoubleFromUser(double minValue);

double getDoubleFromUser(double minValue, double maxValue);

bool isStringInteger(string const str,
                     int *number);

bool isStringDouble(string const str,
                    double *result);


//! Functions in "miscellaneousInput.cc"
void getIntegerVector(char input[],
                      size_t const length,
                      vector<int> *scales,
                      char interval,
                      char split);
void getIntegerVector(string input,
                      vector<int> *scales,
                      char interval,
                      char split);

void copyStringFront(char const *oldString,
                     size_t const oldSize,
                     char *newString,
                     size_t const newSize,
                     string message);
void copyStringEnd(char const *oldString,
                   size_t const oldSize,
                   char *newString,
                   size_t const newSize,
                   string message);

void copyStringFront(string const old,
                     char *newString,
                     size_t const newSize,
                     string message);
void copyStringEnd(string const old,
                   char *newString,
                   size_t const newSize,
                   string message);

bool duplicateEntries(int *array,
                      int const size);

string fullFilePath(string filename);
void initializeCharArray(char *array, size_t const size);


//! Functions in "miscellaneousMath.cc"
int rootN( int const input, int const rootPower);
bool isRootN( int const input, int const rootPower);

double vectorMagnitude(int x1,
                       int x2,
                       int x3);

double vectorMagnitude(double x[], int const n);

double momentumMagnitude(int x1,
                         int x2,
                         int x3,
                         int const nGrid);

double momentumMagnitude(int x1, int const n1,
                         int x2, int const n2,
                         int x3, int const n3);

void log10Steps(double array[],
                double const minimun,
                double const maximum,
                int const size);

double sinc(int const k,
            int const n);

double sinc(double const x,
            double const lowerValue);

double sinc_k(int const k,
              int const n);

double RombergInt(double (*f) (double),
                  double const xMin,
                  double const xMax,
                  double const maxError=1.0e-10 );

double RombergIntInfinity(double (*f) (double),
                          double const xMin,
                          double const xMax,
                          double const maxError=1.0e-10 );



//! Functions in "miscellaneousFile.cc"
void openInputBinaryFile(std::fstream & inputFile,
                         string & fileName);

void openOutputBinaryFile(std::fstream & outputFile,
                          string & fileName);

void openInputTextFile(std::fstream & inputFile,
                       string & fileName);

void openOutputTextFile(std::fstream & outputFile,
                        string & fileName);

void existentFile(string const & fileName);

void getFullPath(string & fileName,
                 char output[],
                 int const maxSize);

std::string getCurrentPath();
std::string getFilename(std::string fileName);






#endif
