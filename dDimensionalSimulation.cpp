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
    //parse the arguments
    cmd.parse( argc, argv );

    int programSwitch = programSwitchArg.getValue();
    int N = nSwitchArg.getValue();
    int maximumIterations = maxIterationsSwitchArg.getValue();
    scalar L = lengthSwitchArg.getValue();
    scalar Temperature = temperatureSwitchArg.getValue();
    scalar phi = phiSwitchArg.getValue();
    scalar rho = rhoSwitchArg.getValue();

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
    shared_ptr<simpleModel> Configuration = make_shared<simpleModel>(N,noise);
    scalar boxL = pow(0.5,1./3.)*pow(N,1./3.);
    Configuration->setRadius(boxL*0.5);
    Configuration->getNeighbors();
            //Configuration->setSoftRepulsion();
    shared_ptr<vectorialVicsek> vicsek = make_shared<vectorialVicsek>();
    vicsek->setEta(0.1);
    vicsek->setV0(0.5);
    vicsek->setDeltaT(.1);

    shared_ptr<Simulation> sim = make_shared<Simulation>();
    sim->setConfiguration(Configuration);
    sim->addUpdater(vicsek,Configuration);

    if(gpuSwitch >=0)
        {
        sim->setCPUOperation(false);
        };

    int rate = 100;

vector<hyperrectangularCellList> cls;
vector<vector<int> > vs;
for (scalar cellLengthSize = 1.0; cellLengthSize < 3.0; cellLengthSize += 0.5)
    {
    hyperrectangularCellList CL1(cellLengthSize,boxL/2);
    CL1.setGPU(gpuSwitch >=0);
    CL1.computeAdjacentCells(1);
    vector<int> v1(CL1.totalCells);
    cls.push_back(CL1);
    vs.push_back(v1);
    }

    for (int ii = 0; ii < maximumIterations; ++ii)
        {
        sim->performTimestep();
        if(ii%rate == 0 )
            {
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
                scalar mean, var;
                getMeanVar(vs[cc],mean,var);
                printf("timestep %i\t m,v = %f , %f\n",ii, mean, var);
                }
            }
        }


//
//The end of the tclap try
//
    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
    return 0;
};
