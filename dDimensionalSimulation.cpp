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

class runningStats
    {
    public:
        runningStats() : m_n(0) {};
        
        void clear(){m_n = 0;};
        void push(scalar x)
            {
            m_n++;
            if(m_n ==1)
                {
                m_oldM = m_newM = x;
                m_oldS = 0.0;
                }
            else
                {
                m_newM = m_oldM+(x-m_oldM)/m_n;
                m_newS = m_oldS+(x-m_oldM)*(x-m_newM);

                m_oldM = m_newM;
                m_oldS = m_newS;
                };
            }

        int numEntries() const {return m_n;};
        scalar mean() const {return (m_n>0) ? m_newM : 0.0;};
        scalar variance() const {return (m_n>1) ? m_newS/(m_n-1) : 0.0;};

    protected:
        int m_n;
        scalar m_oldM, m_oldS,m_newM,m_newS;
    };

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

void getMeanVar(vector<unsigned int> &vec, scalar &mean, scalar &var,int group)
    {
    mean = 0.0;
    var = 0.0;
    for (int ii = 0; ii < floor(vec.size()/group); ++ii)
        {
        scalar partialSum = 0.0;
        for (int gg = 0; gg <group; ++gg)
            partialSum += vec[group*ii+gg];
        mean +=partialSum;
        };
    if (vec.size() > 0)
        mean /= floor(vec.size()/group);

    for (int ii = 0; ii < floor(vec.size()/group); ++ii)
        {
        scalar partialSum = 0.0;
        for (int gg = 0; gg <group; ++gg)
            partialSum += vec[group*ii+gg];
        var += (partialSum-mean)*(partialSum-mean);
        };
    if(vec.size() >0)
        var /=(floor(vec.size()/group)-1);
    }

void computeRelativeFluctuationScaling(shared_ptr<hyperrectangularCellList> cl, vector<vector< runningStats> > &stats, scalar drStep, scalar dtStep, dVec meanDir)
    {
    ArrayHandle<unsigned int> particlesPerCell(cl->elementsPerCell);
    dVec cellSizes = cl->getCellSize();
    scalar aveN = 1;
    for (int dd = 0; dd < DIMENSION; ++dd)
        aveN *= cellSizes[dd];
    aveN*=0.5;//set density correctly

    printf("total cells %i \n" ,cl->totalCells-1);
    for (int b1 = 0; b1 < cl->totalCells-1;++b1)
        {
        iVec cellPos1 = cl->indexToiVec(b1);
        dVec pos1;
        for (int dd = 0; dd < DIMENSION;++dd)
            pos1[dd] = cellSizes[dd]*cellPos1[dd];
        for (int b2 = b1 +1; b2 < cl->totalCells; ++b2)
            {
            dVec distanceBetweenBoxes;
            iVec cellPos2 = cl->indexToiVec(b2);
            dVec pos2;
            for (int dd = 0; dd < DIMENSION;++dd)
                pos2[dd] = cellSizes[dd]*cellPos2[dd];

            cl->Box->minDist(pos1,pos2,distanceBetweenBoxes);
            scalar dispNorm = norm(distanceBetweenBoxes);
            if(dispNorm < drStep*stats.size())
                {
                int rBin = floor(dispNorm / drStep);
                scalar theta = acos(dot(meanDir,distanceBetweenBoxes)/(dispNorm*norm(meanDir)));
                int tBin = floor(theta/dtStep);
                if(tBin < stats[0].size() && rBin < stats.size())
                    {
                    int bc1 = particlesPerCell.data[b1]-aveN;
                    int bc2 = particlesPerCell.data[b2]-aveN;
                    stats[rBin][tBin].push(bc1*bc2);
                    };
                }
            }
        }
    };

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
//    Configuration->setSoftRepulsion(0.2,1.0);
    shared_ptr<vectorialVicsek> vicsek = make_shared<vectorialVicsek>();
    vicsek->setEta(0.2);
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

/*
    vector<hyperrectangularCellList> cls;
    vector<vector<int> > vs;
    scalar lastSize = 0.0;
    int targetCellList = 0;
    for (scalar cellLengthSize = 1.0; cellLengthSize < 2.0; cellLengthSize += 0.25)
        {
        hyperrectangularCellList CL1(cellLengthSize,boxL/2);
        CL1.setGPU(gpuSwitch >=0);
        CL1.computeAdjacentCells(1);
        scalar currentSize = CL1.getCellSize().x[0];
        if(currentSize > lastSize)
            {
            printf("%f grid size with %i total cells \n",currentSize,CL1.totalCells);
            vector<int> v1(CL1.totalCells);
            cls.push_back(CL1);
            vs.push_back(v1);
            lastSize = currentSize;
            if(CL1.totalCells > 1001)
                targetCellList += 1;
            }
        }
    targetCellList = 2;
*/
    shared_ptr<hyperrectangularCellList> CL = make_shared<hyperrectangularCellList>(1.25,boxL/2);
    CL->setGPU(gpuSwitch >=0);
    CL->computeAdjacentCells(1);
    scalar currentSize = CL->getCellSize().x[0];

    scalar drStep = Configuration->metricNeighbors.cellList->getCellSize().x[0];
    drStep = currentSize;

    ///let's say we're thinking of tonerTest[distance][angle]
    int rBins = floor(boxL/2./drStep);
    int nAngles = 40;
    scalar dthetaStep = PI/(nAngles);
    vector<runningStats> angles(nAngles);
    vector<vector< runningStats> > tonerTest(rBins,angles);
    for (int rr = 0; rr < rBins; ++rr)
        for (int tt = 0; tt < nAngles; ++tt)
            tonerTest[rr][tt].clear();
    printf("setting up a grid with %i angular steps and %i radial steps of size %f\n",(int)angles.size(),rBins,drStep);

    //initialize
    int rate = floor(boxL/(v0*dt));
    for (int ii = 0; ii < maximumIterations; ++ii)
        {
        if(ii%rate == 0 )
            printf("%i out of %i initialization steps\n",ii,maximumIterations);

        sim->performTimestep();
        }

    vector<int> cellListGroupings;
    int last = 1;
    while(last *0.5 < 0.15*boxL*boxL*boxL)
        {
        cellListGroupings.push_back(last);
        last = last*2;
        }

    vector<unsigned int> clData;
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
            Configuration->getNeighbors();
            copyGPUArrayData(Configuration->metricNeighbors.cellList->elementsPerCell,clData);
            for (int gg = 0 ; gg < cellListGroupings.size(); ++gg)
                {
                scalar mean, var;
                getMeanVar(clData,mean,var,cellListGroupings[gg]);
                printf("timestep %i\t m,v = %f , %f\t\t %f\n",ii, mean, var,norm(meanDir));
                myfile << "\t" << mean << "\t" << var;
                }

            /*
            for (int cc = 0; cc < cls.size(); ++cc)
                {
                cls[cc].computeCellList(Configuration->returnPositions());
                copyGPUArrayData(cls[cc].elementsPerCell,clData);
                ArrayHandle<unsigned int> particlesPerCell(cls[cc].elementsPerCell);
                int totalCells = cls[cc].totalCells;
                for (int currentCell = 0; currentCell < totalCells; ++currentCell)
                    {
                    int particlesInBin =  particlesPerCell.data[currentCell];
                    vs[cc][currentCell] = particlesInBin;
                    };
                scalar mean, var;
                getMeanVar(vs[cc],mean,var);
                printf("timestep %i\t m,v = %f , %f\t\t %f\n",ii, mean, var,norm(meanDir));
                myfile << "\t" << mean << "\t" << var;
                }
            */
            myfile << "\n";
            cout << "starting fluctation processing...";cout.flush();
            CL->computeCellList(Configuration->returnPositions());
            computeRelativeFluctuationScaling(CL,tonerTest,drStep,dthetaStep,meanDir);
            cout <<" finished." << endl;cout.flush();
            }
        }

    myfile.close();

    //read out the fluctuation statistics of the smallest bins
    sprintf(filename,"../data/smallBinTesting_N%i_rho%.3f.txt",N,N/(1.0*boxL*boxL*boxL));
    ofstream myfile2;
    myfile2.open(filename);myfile2.setf(ios_base::scientific);myfile2 << setprecision(10);
    for (int rr = 0; rr < tonerTest.size(); ++rr)
        {
        for (int tt = 0; tt < tonerTest[0].size(); ++tt)
            {
            scalar m = tonerTest[rr][tt].mean();
            scalar v = tonerTest[rr][tt].variance();
            scalar ans = 0.0;
//            if(m!= 0)
                ans = m;
            myfile2 << ans << "\t";
            }
        myfile2 << "\n";
        };



    myfile2.close();

//
//The end of the tclap try
//
    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
    return 0;
};
