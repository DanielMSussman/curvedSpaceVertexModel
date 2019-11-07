#ifndef vectorValue_h
#define vectorValue_h

#include "baseNetCDFDatabase.h"
#include "simpleDatabase.h"

class vectorValueNetCDF : public baseNetCDFDatabase, public simpleDatabase
    {
    public:
        vectorValueNetCDF(string fn="temp.nc", int N=2, NcFile::FileMode mode = NcFile::ReadOnly);

        void writeState(vector<scalar> &vec, scalar val);
        void readState(vector<scalar> &vec, scalar &val, int rec);
    };

#endif
