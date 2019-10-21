#include "std_include.h" // std library includes, definition of scalar, etc.. has a "using namespace std" in it, because I'm lazy

//we'll use TCLAP as our command line parser
#include <tclap/CmdLine.h>
#include "cuda_profiler_api.h"

#include "functions.h"
#include "gpuarray.h"
#include "periodicBoundaryConditions.h"
#include "simulation.h"
#include "simpleModel.h"
#include "sphericalVertexModel.h"
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


/*!
core of spherical vertexmodel
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
    ValueArg<int> maxIterationsSwitchArg("i","iterations","number of timestep iterations",false,0,"int",cmd);
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

    /*
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
    */

    int dim =DIMENSION;
    noiseSource noise(true);
    shared_ptr<sphericalVertexModel> Configuration = make_shared<sphericalVertexModel>(N,noise,GPU,!GPU);

    printf("sphere size  = %f\n",Configuration->sphere.radius);
    shared_ptr<Simulation> sim = make_shared<Simulation>();
    sim->setConfiguration(Configuration);
    scalar temperature  =0.0;
    shared_ptr<noseHooverNVT> NVT = make_shared<noseHooverNVT>(Configuration,temperature);
    sim->setIntegrationTimestep(0.001);
    sim->addUpdater(NVT,Configuration);

    if(gpuSwitch >=0)
        {
        sim->setCPUOperation(false);
        };

    Configuration->computeForces(true);

    for (int ii = 0; ii <maximumIterations; ++ii)
        sim->performTimestep();
//
//The end of the tclap try
//
    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
    return 0;
};
