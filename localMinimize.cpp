#include "std_include.h" // std library includes, definition of scalar, etc.. has a "using namespace std" in it, because I'm lazy

//we'll use TCLAP as our command line parser
#include <tclap/CmdLine.h>
#include "cuda_profiler_api.h"

#include "functions.h"
#include "profiler.h"
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
#include "brownianDynamics.h"
#include "noiseSource.h"
#include "harmonicRepulsion.h"
#include "lennardJones6_12.h"
#include "indexer.h"
#include "hyperrectangularCellList.h"
#include "neighborList.h"
#include "poissonDiskSampling.h"
#include "sphericalVoronoi.h"
#include "vectorValueNetCDF.h"
#include "simpleUtilities.h"
#include "analysisPackage.h"

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
    ValueArg<int> maxIterationsSwitchArg("i","iterations","number of timestep iterations",false,100000,"int",cmd);
    ValueArg<int> fileIdxSwitch("f","file","file Index",false,-1,"int",cmd);
    ValueArg<scalar> forceToleranceSwitchArg("t","fTarget","target minimization threshold for norm of residual forces",false,0.000000000000001,"scalar",cmd);
    ValueArg<scalar> lengthSwitchArg("l","sideLength","size of simulation domain",false,10.0,"double",cmd);
    ValueArg<scalar> krSwitchArg("k","springRatio","kA divided by kP",false,1.0,"double",cmd);

    //allow setting of system size by either volume fraction or density (assuming N has been set)
    ValueArg<scalar> a0SwitchArg("a","a0","preferred area per cell",false,1.0,"double",cmd);
    ValueArg<scalar> rhoSwitchArg("r","rho","density",false,-1.0,"double",cmd);
    ValueArg<scalar> dtSwitchArg("e","dt","timestep",false,0.001,"double",cmd);
    ValueArg<scalar> v0SwitchArg("v","v0","v0",false,0.5,"double",cmd);
    ValueArg<scalar> p0SwitchArg("p","p0","preferred perimeter",false,3.78,"double",cmd);
    ValueArg<scalar> pMaxSwitchArg("q","pmax","preferred perimeter",false,3.79,"double",cmd);
    ValueArg<scalar> dpSwitchArg("w","dp","preferred perimeter",false,.01,"double",cmd);
    //parse the arguments
    cmd.parse( argc, argv );

    int programSwitch = programSwitchArg.getValue();
    int fIdx = fileIdxSwitch.getValue();
    int N = nSwitchArg.getValue();
    int maximumIterations = maxIterationsSwitchArg.getValue();
    scalar L = lengthSwitchArg.getValue();
    scalar rho = rhoSwitchArg.getValue();
    scalar dt = dtSwitchArg.getValue();
    scalar v0 = v0SwitchArg.getValue();
    scalar p0 = p0SwitchArg.getValue();
    scalar pMax = pMaxSwitchArg.getValue();
    scalar dp = dpSwitchArg.getValue();
    scalar a0 = a0SwitchArg.getValue();
    scalar kr = krSwitchArg.getValue();
    scalar forceCutoff = forceToleranceSwitchArg.getValue();

    int gpuSwitch = gpuSwitchArg.getValue();
    bool GPU = false;
    if(gpuSwitch >=0)
        GPU = chooseGPU(gpuSwitch);

    int dim =DIMENSION;
    bool reproducible = fIdx <=0 ? true : false;
    noiseSource noise(reproducible);
    vector<scalar2> results;
for (scalar pVal = p0; pVal < pMax; pVal += dp)
    {
    vector<scalar> outputVec(3);
    char fname[256];
    sprintf(fname,"data/minimizationResults_N%i_p%.4f_kr%.3f.nc",N,pVal,kr);
    string outFile(fname);
    vectorValueNetCDF vvdat(outFile,3,NcFile::Replace);

    int numAtZero = 0;
    int minimizedStates = 0;
    for (int ff = 0; ff < fIdx; ++ff)
        {
        noise.setReproducibleSeed(ff+1137);

        shared_ptr<sphericalVertexModel> Configuration = make_shared<sphericalVertexModel>(N,noise,a0,pVal,GPU,!GPU);
        Configuration->setScalarModelParameter(kr);
        printf("sphere size  = %f\n",Configuration->sphere->radius);

        shared_ptr<energyMinimizerFIRE> fire =  make_shared<energyMinimizerFIRE>(Configuration);
        shared_ptr<Simulation> sim = make_shared<Simulation>();
        sim->setConfiguration(Configuration);
        fire->setMaximumIterations(maximumIterations);

        scalar alphaStart=.99; scalar deltaTInc=1.1; scalar deltaTDec=0.95;
        scalar alphaDec=0.9; int nMin=4;scalar alphaMin = .0;
        scalar deltaTMax = 50*dt;
        fire->setFIREParameters(dt,alphaStart,deltaTMax,deltaTInc,deltaTDec,alphaDec,nMin,forceCutoff,alphaMin);

        sim->addUpdater(fire,Configuration);
        sim->setCPUOperation(!GPU);//have cpu and gpu initialized the same...for debugging

        sim->performTimestep();
        scalar energy = sim->computeEnergy();
        scalar meanForce = fire->getMaxForce();
        int iters = fire->getCurrentIterations();
        printf("iterations %i\t energy %g\t force %g\n",iters,energy,meanForce);
        outputVec[0] = energy;
        outputVec[1] = meanForce;
        outputVec[2] = iters;
        vvdat.writeState(outputVec,ff);
        if(meanForce < forceCutoff)
            {
            minimizedStates +=1;
            if(energy < 1e-10)
                numAtZero +=1;
            }
        }
        printf("fraction of mechanically unstable states: %f\n",(scalar)numAtZero/(scalar)minimizedStates);
        scalar2 ans; ans.x = pVal; ans.y = (scalar)numAtZero/(scalar)minimizedStates;
        results.push_back(ans);
    }
    printf("{");
    for(int rr = 0; rr < results.size()-1; ++rr)
        printf("{%f,%f},",results[rr].x,results[rr].y);
    printf("{%f,%f}}\n",results[results.size()-1].x,results[results.size()-1].y);
//
//The end of the tclap try
//
    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
    return 0;
};
