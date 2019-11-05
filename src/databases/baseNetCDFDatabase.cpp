#include "baseNetCDFDatabase.h"

baseNetCDFDatabase::baseNetCDFDatabase(string fn, NcFile::FileMode mode)
     : file(fn.c_str(), mode)
{
    NcError err(NcError::silent_nonfatal);

    addDimension("record",0);
}

void baseNetCDFDatabase::addVariable(string _name, int _dimIdx1, int _dimIdx2)
    {
    pair<string,pair<int,int>> varPair;
    varPair.first = _name;
    varPair.second.first = _dimIdx1;
    varPair.second.second = _dimIdx2;
    varSpecifiers.push_back(varPair);
    }
void baseNetCDFDatabase::addDimension(string _name, int _size)
    {
    pair<string,int> dimPair;
    dimPair.first = _name;
    dimPair.second = _size;
    dimSpecifiers.push_back(dimPair);
    }

void baseNetCDFDatabase::getDimVar()
    {
    NcDim *fDim;
    for (int dd = 0; dd < dimSpecifiers.size(); ++dd)
        {
        fDim = file.get_dim(dimSpecifiers[dd].first.c_str());
        dimPointers.push_back(fDim);
        }
    NcVar *fVar;
    for (int dd = 0; dd < varSpecifiers.size(); ++dd)
        {
        fVar = file.get_var(varSpecifiers[dd].first.c_str());
        varPointers.push_back(fVar);
        }
    }

void baseNetCDFDatabase::setDimVar()
    {
    NcDim *fDim;
    for (int dd = 0; dd < dimSpecifiers.size(); ++dd)
        {
        if(dimSpecifiers[dd].second ==0)
            {
            fDim = file.add_dim(dimSpecifiers[dd].first.c_str());
            dimPointers.push_back(fDim);
            }
        else
            {
            fDim = file.add_dim(dimSpecifiers[dd].first.c_str(),dimSpecifiers[dd].second);
            dimPointers.push_back(fDim);
            }
        }
    NcVar *fVar;
    for (int dd = 0; dd < varSpecifiers.size(); ++dd)
        {
        int i1 = varSpecifiers[dd].second.first;
        int i2 = varSpecifiers[dd].second.second;
        string vName = varSpecifiers[dd].first;
        NcDim *d1 = dimPointers[i1];
        if(i2 > -1)
            {
            NcDim *d2 = dimPointers[i2];
            fVar = file.add_var(vName.c_str(),ncDouble,d1,d2);
            }
        else
            {
            NcDim *d2 = dimPointers[0];
            fVar = file.add_var(vName.c_str(),ncDouble,d2,d1);
            }
        varPointers.push_back(fVar); 
        }
    file.sync();
    }
