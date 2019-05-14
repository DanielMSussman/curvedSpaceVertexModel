#ifndef convexHullCGAL_H
#define convexHullCGAL_H

#include <cmath>
#include <algorithm>
#include <memory>
#include <ctype.h>
#include <random>
#include <stdio.h>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/time.h>
#include <string.h>
#include <stdexcept>
#include <cassert>

#include "vector_types.h"
#include <vector>

#include "dDimensionalVectorTypes.h"
class convexHullCGALInterface
    {
    public:
        std::vector< std::vector<int> > allNeighs; //!<The list of neighbors of every point in the convex hull
        std::vector<int> numNeighs;
        void sphericalConvexHull(dVec *points, int n);
    };
#endif
