#include "std_include.h" // std library includes, definition of scalar, etc.. has a "using namespace std" in it, because I'm lazy

//we'll use TCLAP as our command line parser
#include <tclap/CmdLine.h>
#include "cuda_profiler_api.h"

#include "functions.h"
#include "gpuarray.h"
#include "periodicBoundaryConditions.h"
#include "simulation.h"
#include "simpleModel.h"
#include "sphericalVectorialVicsek.h"
#include "baseUpdater.h"
#include "energyMinimizerFIRE.h"
#include "velocityVerlet.h"
#include "noseHooverNVT.h"
#include "noiseSource.h"
#include "harmonicRepulsion.h"
#include "lennardJones6_12.h"
#include "indexer.h"
#include "hyperrectangularCellList.h"
#include "neighborList.h"
#include "poissonDiskSampling.h"
#include "sphericalVoronoi.h"

using namespace std;
using namespace TCLAP;

//!What, after all, *is* the volume of a d-dimensional sphere?
scalar sphereVolume(scalar radius, int dimension)
    {
    if(dimension == 1)
        return 2*radius;
    else
        if(dimension == 2)
            return PI*radius*radius;
        else
            return (2.*PI*radius*radius)/((scalar) dimension)*sphereVolume(radius,dimension-2);
    };

void getMeanVar(vector<int> &vec, scalar &mean, scalar &var)
    {
    mean = 0.0;
    var = 0.0;
    for (int ii = 0; ii < vec.size(); ++ii)
        mean +=vec[ii];
    if (vec.size() > 0)
        mean /= vec.size();
    for (int ii = 0; ii < vec.size(); ++ii)
        var += (vec[ii]-mean)*(vec[ii]-mean);
    if(vec.size() >0)
        var /=(vec.size()-1);
    }

/*!
This file runs some dynamics on particles interacting according to some
potential... when this repository is ever meant to be used this all should be
updated.
*/
int main(int argc, char*argv[])
{
    // wrap tclap in a try block
    try
    {
    //First, we set up a basic command line parser...
    //      cmd("command description message", delimiter, version string)
    CmdLine cmd("basic testing of dDimSim", ' ', "V0.1");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programSwitchArg("z","programSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> gpuSwitchArg("g","USEGPU","an integer controlling which gpu to use... g < 0 uses the cpu",false,-1,"int",cmd);
    ValueArg<int> nSwitchArg("n","Number","number of particles in the simulation",false,100,"int",cmd);
    ValueArg<int> maxIterationsSwitchArg("i","iterations","number of timestep iterations",false,100,"int",cmd);
    ValueArg<scalar> lengthSwitchArg("l","sideLength","size of simulation domain",false,10.0,"double",cmd);
    ValueArg<scalar> temperatureSwitchArg("t","temperature","temperature of simulation",false,.001,"double",cmd);

    //allow setting of system size by either volume fraction or density (assuming N has been set)
    scalar phiDest = 1.90225*exp(-(scalar)DIMENSION / 2.51907);
    ValueArg<scalar> phiSwitchArg("p","phi","volume fraction",false,phiDest,"double",cmd);
    ValueArg<scalar> rhoSwitchArg("r","rho","density",false,-1.0,"double",cmd);
    ValueArg<scalar> dtSwitchArg("e","dt","timestep",false,0.1,"double",cmd);
    ValueArg<scalar> v0SwitchArg("v","v0","v0",false,0.5,"double",cmd);
    //parse the arguments
    cmd.parse( argc, argv );

    int programSwitch = programSwitchArg.getValue();
    int N = nSwitchArg.getValue();
    int maximumIterations = maxIterationsSwitchArg.getValue();
    scalar L = lengthSwitchArg.getValue();
    scalar Temperature = temperatureSwitchArg.getValue();
    scalar phi = phiSwitchArg.getValue();
    scalar rho = rhoSwitchArg.getValue();
    scalar dt = dtSwitchArg.getValue();
    scalar v0 = v0SwitchArg.getValue();

    int gpuSwitch = gpuSwitchArg.getValue();
    bool GPU = false;
    if(gpuSwitch >=0)
        GPU = chooseGPU(gpuSwitch);

    if(phi >0)
        {
        L = pow(N*sphereVolume(.5,DIMENSION) / phi,(1.0/(scalar) DIMENSION));
        rho = N/pow(L,(scalar)DIMENSION);
        }
    else
        phi = N*sphereVolume(.5,DIMENSION) / pow(L,(scalar)DIMENSION);

    if(rho >0)
        {
        L = pow(((scalar)N/rho),(1.0/(scalar) DIMENSION));
        phi = rho * sphereVolume(.5,DIMENSION);
        }
    else
        rho = N/pow(L,(scalar)DIMENSION);


    int dim =DIMENSION;
    noiseSource noise(true);
    shared_ptr<simpleModel> Configuration = make_shared<simpleModel>(N,noise,GPU,!GPU);
    scalar boxL = pow(2.,1./3.)*pow(N,1./3.);
    Configuration->setRadius(boxL*0.5);
    printf("new box side length = %f\n",boxL);
    Configuration->getNeighbors();
            //Configuration->setSoftRepulsion();
    shared_ptr<vectorialVicsek> vicsek = make_shared<vectorialVicsek>();
    vicsek->setEta(0.1);
    vicsek->setV0(v0);
    vicsek->setDeltaT(dt);

    shared_ptr<Simulation> sim = make_shared<Simulation>();
    sim->setConfiguration(Configuration);
    sim->addUpdater(vicsek,Configuration);

    if(gpuSwitch >=0)
        {
        sim->setCPUOperation(false);
        };

    char filename[256];
    sprintf(filename,"../data/GNFTesting_N%i_rho%.3f.txt",N,N/(1.0*boxL*boxL*boxL));
    ofstream myfile;
    myfile.open(filename);myfile.setf(ios_base::scientific);myfile << setprecision(10);

    vector<hyperrectangularCellList> cls;
    vector<vector<int> > vs;
    scalar lastSize = 0.0;
    for (scalar cellLengthSize = 1.0; cellLengthSize < boxL/5; cellLengthSize += 0.25)
        {
        hyperrectangularCellList CL1(cellLengthSize,boxL/2);
        CL1.setGPU(gpuSwitch >=0);
        CL1.computeAdjacentCells(1);
        scalar currentSize = CL1.getCellSize().x[0];
        if(currentSize > lastSize)
            {
            printf("%f\n",currentSize);
            vector<int> v1(CL1.totalCells);
            cls.push_back(CL1);
            vs.push_back(v1);
            lastSize = currentSize;
            }
        }

    vector<scalar3> smallCellPositions(cls[0].totalCells);
    vector<vector< int> > smallBins;

    for (int ii = 0; ii < maximumIterations; ++ii)
        sim->performTimestep();

    int rate = floor(boxL/(v0*dt));

    dVec meanDir;
    for (int ii = 0; ii < maximumIterations; ++ii)
        {
        sim->performTimestep();
        if(ii%rate == 0 )
            {
            Configuration->getMeanDirection(meanDir);
            myfile << ii;
            for (int dd = 0; dd < DIMENSION; ++dd)
                myfile << "\t" << meanDir.x[dd];
            for (int cc = 0; cc < cls.size(); ++cc)
                {
                cls[cc].computeCellList(Configuration->returnPositions());
                ArrayHandle<unsigned int> particlesPerCell(cls[cc].elementsPerCell);
                int totalCells = cls[cc].totalCells;
                for (int currentCell = 0; currentCell < totalCells; ++currentCell)
                    {
                    int particlesInBin =  particlesPerCell.data[currentCell];
                    vs[cc][currentCell] = particlesInBin;
                    };
                if (cc == 0)
                    smallBins.push_back(vs[cc]);
                scalar mean, var;
                getMeanVar(vs[cc],mean,var);
                printf("timestep %i\t m,v = %f , %f\n",ii, mean, var);
                myfile << "\t" << mean << "\t" << var;
                }
            myfile << "\n";
            }
        }

    myfile.close();

    //read out the fluctuation statistics of the smallest bins
    sprintf(filename,"../data/smallBinTesting_N%i_rho%.3f.txt",N,N/(1.0*boxL*boxL*boxL));
    ofstream myfile2;
    myfile2.open(filename);myfile2.setf(ios_base::scientific);myfile2 << setprecision(10);


    dVec cellSizes = cls[0].getCellSize();
    for (int bb = 0; bb < smallBins[0].size();++bb)
        {
        iVec cellPos = cls[0].indexToiVec(bb);
        vector<int> counts(smallBins.size());
        for (int tt =0; tt < counts.size();++tt)
            counts[tt] = smallBins[tt][bb];
        scalar mean, var;
        getMeanVar(counts,mean,var);
        for (int dd = 0; dd < DIMENSION; ++dd)
            myfile2 << cellSizes.x[dd]*cellPos.x[dd] << "\t";
        myfile2 << mean << "\t" << var << "\n";
        }

    myfile2.close();
//
//The end of the tclap try
//
    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
    return 0;
};
