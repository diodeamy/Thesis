#ifndef THRESHOLD_HEADER
#define THRESHOLD_HEADER

#include <iostream>
#include <iomanip>
#include <cmath>
#include <list>
#include <utility>
#include <algorithm>

#include "../defines/defines.h"
#include "../classes/array.h"
#include "objects.h"


// returns filament and wall properties needed for computing the optimal threshold
void thresholdProperties(Array<Real,3> &response,
                         Real const minThreshold,
                         Real const maxThreshold,
                         Array<int,3> *mask,
                         int const neighborFindingMethod,
                         Array<Real,3> &mass,
                         Int const minimumSize,
                         Real *properties,
                         int *trackStatistics);

// returns node properties needed for computing the optimal threshold
Real fractionVirialNodes(Array<Real,3> &response,
                         Real const minThreshold,
                         Real const maxThreshold,
                         Array<int,3> *mask,
                         int const neighborFindingMethod,
                         Array<Real,3> &mass,
                         double const minSize,
                         double const virialDensity,
                         Real *properties,
                         int *trackStatistics);



#endif


