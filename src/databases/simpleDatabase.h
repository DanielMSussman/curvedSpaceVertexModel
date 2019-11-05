#ifndef simpleDatabase_h
#define simpleDatabase_h

#include "std_include.h"
#include "simpleModel.h"

/*! \file simpleDatabase.h */
//! A base class defining the operations of a database save and read scheme

class simpleDatabase
    {
    protected:
        typedef shared_ptr<simpleModel> STATE;

    public:
        //! Base constructure takes a bland filename in readonly mode
        simpleDatabase(string fn="temp.txt", int _mode=-1):filename(fn), mode(_mode),records(0){};
        //!The name of the file
        string filename;
        //!The desired mode (integer representation of replace, new, write, readonly, etc)
        const int mode;
        //!The number of saved records in the database
        int records;

        //!Write the current state; if the default value of rec=-1 is used, add a new entry
        virtual void writeState(STATE c, scalar time = -1.0, int rec = -1){};
        //Read the rec state of the database. If geometry = true, call computeGeomerty routines (instead of just reading in the d.o.f.s)
        virtual void readState(STATE c, int rec, bool geometry = true){};
    
    };
#endif

