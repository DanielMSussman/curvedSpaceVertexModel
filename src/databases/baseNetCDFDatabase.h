#ifndef baseNetCDF_H
#define baseNetCDF_H

#include <netcdfcpp.h>
#include "std_include.h"
/*! \file baseNetCDFDatabase.h */
//! A base class that implements a details-free  netCDF4-based data storage system
/*!
BaseDatabase just provides an interface to a file and a mode of operation.
*/
class baseNetCDFDatabase
    {
    public:
        //!The NcFile itself
        NcFile file;

        //!The default constructor starts a bland filename in readonly mode
        baseNetCDFDatabase(string fn="temp.nc", NcFile::FileMode mode=NcFile::ReadOnly);

        vector<NcDim*> dimPointers;
        vector<NcVar*> varPointers;

        vector<pair<string, int> > dimSpecifiers;
        vector<pair<string, pair<int,int> > >varSpecifiers;

        int getNumberOfRecords()
            {
            NcDim *rd = file.get_dim("record");
            return rd->size();
            }
        void addDimension(string _name, int _size);
        void addVariable(string _name, int dimIdx1, int dimIdx2= -1);
        void setDimVar();
        void getDimVar();
    };

#endif
