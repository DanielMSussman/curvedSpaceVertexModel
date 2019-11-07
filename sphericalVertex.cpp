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
    ValueArg<int> maxIterationsSwitchArg("i","iterations","number of timestep iterations",false,0,"int",cmd);
    ValueArg<int> fileIdxSwitch("f","file","file Index",false,0,"int",cmd);
    ValueArg<scalar> lengthSwitchArg("l","sideLength","size of simulation domain",false,10.0,"double",cmd);
    ValueArg<scalar> temperatureSwitchArg("t","temperature","temperature of simulation",false,.001,"double",cmd);

    //allow setting of system size by either volume fraction or density (assuming N has been set)
    ValueArg<scalar> p0SwitchArg("p","p0","preferred perimeter",false,3.78,"double",cmd);
    ValueArg<scalar> a0SwitchArg("a","a0","preferred area per cell",false,1.0,"double",cmd);
    ValueArg<scalar> rhoSwitchArg("r","rho","density",false,-1.0,"double",cmd);
    ValueArg<scalar> dtSwitchArg("e","dt","timestep",false,0.001,"double",cmd);
    ValueArg<scalar> v0SwitchArg("v","v0","v0",false,0.5,"double",cmd);
    //parse the arguments
    cmd.parse( argc, argv );

    int programSwitch = programSwitchArg.getValue();
    int fIdx = fileIdxSwitch.getValue();
    int N = nSwitchArg.getValue();
    int maximumIterations = maxIterationsSwitchArg.getValue();
    scalar L = lengthSwitchArg.getValue();
    scalar Temperature = temperatureSwitchArg.getValue();
    scalar rho = rhoSwitchArg.getValue();
    scalar dt = dtSwitchArg.getValue();
    scalar v0 = v0SwitchArg.getValue();
    scalar p0 = p0SwitchArg.getValue();
    scalar a0 = a0SwitchArg.getValue();

    int gpuSwitch = gpuSwitchArg.getValue();
    bool GPU = false;
    if(gpuSwitch >=0)
        GPU = chooseGPU(gpuSwitch);

    int dim =DIMENSION;
    bool reproducible = fIdx ==0 ? true : false;
    noiseSource noise(reproducible);
    shared_ptr<sphericalVertexModel> Configuration = make_shared<sphericalVertexModel>(N,noise,a0,p0,GPU,!GPU);
    printf("sphere size  = %f\n",Configuration->sphere->radius);
    
    shared_ptr<brownianDynamics> BD = make_shared<brownianDynamics>(reproducible);
    BD->setT(Temperature);
    shared_ptr<Simulation> sim = make_shared<Simulation>();
    sim->setConfiguration(Configuration);
    sim->addUpdater(BD,Configuration);
    sim->setIntegrationTimestep(dt);

    if(gpuSwitch >=0)
        {
        sim->setCPUOperation(false);
        };
    sim->setReproducible(reproducible);

    int stepsPerTau = floor(1./dt);
    //initialize
    cout << "initialization..." << endl;
    for (int ii = 0; ii < min(maximumIterations,1000*stepsPerTau); ++ii)
        {
        sim->performTimestep();
        }
DEBUGCODEHELPER;
    dynamicalFeatures dynFeat(Configuration->cellPositions,Configuration->sphere);
DEBUGCODEHELPER;
    logSpacedIntegers lsi(0,0.05);
    lsi.update();

    char fname[256];
    sprintf(fname,"cellPositions_N%i_p%.4f_T%.5f_fidx%i.nc",N,p0,Temperature,fIdx);
    string outFile(fname);
    sprintf(fname,"msd_N%i_p%.4f_T%.5f_fidx%i.nc",N,p0,Temperature,fIdx);
    string outFile2(fname);
    sprintf(fname,"overlap_N%i_p%.4f_T%.5f_fidx%i.nc",N,p0,Temperature,fIdx);
    string outFile3(fname);

    vectorValueNetCDF vvdat(outFile,N*3,NcFile::Replace);
    vectorValueNetCDF msddat(outFile2,2,NcFile::Replace);
    vectorValueNetCDF overlapdat(outFile3,2,NcFile::Replace);
    vector<scalar> outputVec(2);
    vector<scalar> cellPositions(3*N);
    for (int ii = 0; ii <maximumIterations; ++ii)
        {
        if(ii == lsi.nextSave)
            {
            lsi.update();
            scalar e = sim->computeEnergy();
            ArrayHandle<dVec> cp(Configuration->cellPositions);
            for(int cc = 0; cc < N; ++cc)
                for (int dd=0;dd<3;++dd)
                    cellPositions[3*cc+dd] = cp.data[cc][dd];
            vvdat.writeState(cellPositions,e);
            scalar msd = dynFeat.computeMSD(Configuration->cellPositions);
            scalar overlap = dynFeat.computeOverlapFunction(Configuration->cellPositions);
            outputVec[0] = ii*dt;
            outputVec[1] = msd;
            msddat.writeState(outputVec,ii);
            outputVec[1] = overlap;
            overlapdat.writeState(outputVec,ii);
            cout << "wrote state at timestep " << ii << endl;
            }
        sim->performTimestep();
        }


//
//The end of the tclap try
//
    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
    return 0;
};
