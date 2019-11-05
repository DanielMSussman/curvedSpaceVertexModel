#include "vectorValueNetCDF.h"

vectorValueNetCDF::vectorValueNetCDF(string fn, int N, NcFile::FileMode mode) : baseNetCDFDatabase(fn,mode), simpleDatabase(fn)
    {
    addDimension("unit",1);
    addDimension("N",N);
    addVariable("value",1);
    addVariable("vector",2);

    switch(mode)
        {
        case NcFile::ReadOnly:
            getDimVar();
            break;
        case NcFile::Replace:
            setDimVar();
            break;
        case NcFile::Write:
            setDimVar();
            break;
        case NcFile::New:
            setDimVar();
            break;
        default:
            ;
        };

    }

void vectorValueNetCDF::writeState(vector<scalar> &vec,scalar val)
    {
    int rec = dimPointers[0]->size();
    varPointers[0]->put_rec(&val,rec);
    varPointers[1]->put_rec(&vec[0],rec);
    file.sync();
    };

void vectorValueNetCDF::readState(vector<scalar> &vec,scalar &val, int rec)
    {
    /*
    int totalRecords = GetNumRecs();
    if (rec >= totalRecords)
        {
        printf("Trying to read a database entry that does not exist\n");
        throw std::exception();
        };
        vecVar->set_cur(rec);
        vecVar->get(&vec[0],1,dofDim->size());
        valVar->set_cur(rec);
        valVar->get(&val,1,1);
        */
    };

