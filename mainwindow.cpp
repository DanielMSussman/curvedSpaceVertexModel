
#include <QMainWindow>
#include <QGuiApplication>

//#include <Qt3DCore/QEntity>
//#include <Qt3DRender/QCamera>
//#include <Qt3DRender/QCameraLens>
//#include <Qt3DCore/QTransform>
//#include <Qt3DCore/QAspectEngine>

//#include <Qt3DInput/QInputAspect>

//#include <Qt3DRender/QRenderAspect>
//#include <Qt3DExtras/Qt3DWindow>
//#include <Qt3DExtras/QForwardRenderer>
//#include <Qt3DExtras/QPhongMaterial>
//#include <Qt3DExtras/QCylinderMesh>
//#include <Qt3DExtras/QSphereMesh>
//#include <Qt3DExtras/QTorusMesh>

#include <QPropertyAnimation>
#include <chrono>

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::~MainWindow()
{
    delete ui;
}


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

    ui->setupUi(this);

    ui->fileSaveWidget->hide();
    /*
    ui->setPhaseConstants->hide();
    ui->setDistortionConstants1->hide();
    ui->setDistortionConstants2->hide();
    ui->setDistortionConstants3->hide();
    ui->fireParametersWidget->hide();
    ui->addObjectsWidget->hide();
    ui->fileImportWidget->hide();

    ui->multithreadingWidget->hide();
    ui->nesterovWidget->hide();
    ui->applyFieldWidget->hide();
    ui->moveObjectWidget->hide();
    ui->colloidalEvolutionWidget->hide();
    ui->colloidalTrajectoryWidget->hide();
    ui->LOLBFGSWidget->hide();
    */
    connect(ui->displayZone,SIGNAL(xRotationChanged(int)),ui->xRotSlider,SLOT(setValue(int)));
    connect(ui->displayZone,SIGNAL(zRotationChanged(int)),ui->zRotSlider,SLOT(setValue(int)));

    hideControls();
    QString printable = QStringLiteral("Welcome!");
    ui->testingBox->setText(printable);
}

void MainWindow::hideControls()
{
    /*
    ui->label_41->hide();ui->label_43->hide();ui->label_44->hide();ui->label_45->hide(); ui->label_40->hide(); ui->label_42->hide();
    ui->label_12->hide();ui->label_13->hide();ui->label_56->hide();ui->label_57->hide(); ui->label_7->hide(); ui->label_39->hide();
    ui->resetQTensorsButton->hide();
    ui->minimizeButton->hide();
    ui->addObjectButton->hide();
    ui->minimizationParametersButton->hide();
    ui->addIterationsButton->hide();
    ui->addIterationsBox->hide();
    ui->displayZone->hide();
    ui->drawStuffButton->hide();
    ui->latticeSkipBox->hide();
    ui->directorScaleBox->hide();
    ui->xRotSlider->hide();
    ui->zRotSlider->hide();
    ui->zoomSlider->hide();
    ui->visualProgressCheckBox->hide();
    ui->defectThresholdBox->hide();
    ui->defectDrawCheckBox->hide();
    ui->progressBar->hide();
    ui->reprodicbleRNGBox->hide();
    ui->globalAlignmentCheckBox->hide();
    ui->builtinBoundaryVisualizationBox->hide();
    ui->boundaryFromFileButton->hide();
    ui->nesterovMinimizationButton->hide();
    ui->computeEnergyButton->hide();
    ui->lolbfgsMinimizationButton->hide();
    */
}
void MainWindow::showControls()
{
    /*
    ui->label_41->show();ui->label_43->show();ui->label_44->show();ui->label_45->show();
    ui->label_39->show();ui->label_42->show();ui->label_40->show();ui->label_7->show();
    ui->label_12->show();ui->label_13->show();ui->label_56->show();ui->label_57->show();
    ui->defectDrawCheckBox->show();
    ui->resetQTensorsButton->show();
    ui->minimizeButton->show();
    ui->addObjectButton->show();
    ui->minimizationParametersButton->show();
    ui->addIterationsButton->show();
    ui->addIterationsBox->show();
    ui->displayZone->show();
    ui->drawStuffButton->show();
    ui->latticeSkipBox->show();
    ui->directorScaleBox->show();
    ui->xRotSlider->show();
    ui->zRotSlider->show();
    ui->zoomSlider->show();
    ui->visualProgressCheckBox->show();
    ui->defectThresholdBox->show();
    ui->progressBar->show();
    ui->reprodicbleRNGBox->show();
    ui->globalAlignmentCheckBox->show();
    ui->builtinBoundaryVisualizationBox->show();
    ui->boundaryFromFileButton->show();
    ui->nesterovMinimizationButton->show();
    ui->computeEnergyButton->show();
    //ui->lolbfgsMinimizationButton->show();
    */
}

void MainWindow::on_initializeButton_released()
    {
    N = ui->boxNTotalSize->text().toInt();
    v0 =ui->initialSpeed->text().toDouble();
    eta =ui->initialEta->text().toDouble();
    dt =ui->initialDt->text().toDouble();
    noise.Reproducible= ui->reproducibleButton->isChecked();
    ui->initializationFrame->hide();
    if(noise.Reproducible)
        {
        ui->reprodicbleRNGBox->setChecked(true);
        }
    else
        {
        ui->reprodicbleRNGBox->setChecked(false);
        }

    simulationInitialize();

    sim->setCPUOperation(true);
    ui->progressBar->setValue(50);


    QString printable = QStringLiteral("N %1").arg(Configuration->getNumberOfParticles());
    ui->testingBox->setText(printable);
    ui->progressBar->setValue(100);
    on_drawStuffButton_released();
}

void MainWindow::simulationInitialize()
{
    Configuration = make_shared<sphericalVoronoi>(N,noise);
    sim = make_shared<Simulation>();
    sim->setConfiguration(Configuration);
    /*

     landauLCForce = make_shared<landauDeGennesLC>();

     landauLCForce->setPhaseConstants(A,B,C);
     landauLCForce->setModel(Configuration);
     sim->addForce(landauLCForce);
     on_fireParamButton_released();
     ui->reproducibleButton->setEnabled(true);
     */
}
   /*
void MainWindow::on_minimizeButton_released()
{

    bool graphicalProgress = ui->visualProgressCheckBox->isChecked();
    auto upd = sim->updaters[0].lock();

    ui->progressBar->setValue(0);
    QString printable1 = QStringLiteral("minimizing");
    ui->testingBox->setText(printable1);
    auto t1 = chrono::system_clock::now();
    int initialIterations = upd->getCurrentIterations();
    if(!graphicalProgress)
        sim->performTimestep();
    else
    {
        int stepsToTake = upd->getMaxIterations();
        for (int ii = 1; ii <= 10; ++ii)
        {
            upd->setMaximumIterations(upd->getCurrentIterations()+stepsToTake/10);
            sim->performTimestep();
            on_drawStuffButton_released();
            ui->progressBar->setValue(10*ii);
            QString printable2 = QStringLiteral("minimizing");
            ui->testingBox->setText(printable2);
        };
    };
    int iterationsTaken = upd->getCurrentIterations() - initialIterations;
    ui->progressBar->setValue(50);
    auto t2 = chrono::system_clock::now();
    chrono::duration<scalar> diff = t2-t1;
    ui->progressBar->setValue(75);

    ui->progressBar->setValue(80);
    scalar maxForce = sim->getMaxForce();
    QString printable = QStringLiteral("minimization iterations took %2 total time for %3 steps...<f> = %4 ")
                .arg(diff.count()).arg(iterationsTaken).arg(maxForce);
    ui->testingBox->setText(printable);
    ui->progressBar->setValue(100);

}
  */
void MainWindow::on_resetQTensorsButton_released()
{
    /*
    ui->progressBar->setValue(0);
    QString printable1 = QStringLiteral("resetting q tensors ");
    ui->progressBar->setValue(20);
    ui->testingBox->setText(printable1);
    ui->progressBar->setValue(40);
    scalar S0 = (-B+sqrt(B*B-24*A*C))/(6*C);
    ui->progressBar->setValue(60);
    if(noise.Reproducible)
        noise.setReproducibleSeed(13377);
    bool globalAlignment = ui->globalAlignmentCheckBox->isChecked();
    Configuration->setNematicQTensorRandomly(noise,S0,globalAlignment);

    ui->progressBar->setValue(80);
    QString printable = QStringLiteral("Qtensor values reset at S0=%1...").arg(S0);
    ui->testingBox->setText(printable);
    ui->progressBar->setValue(100);
    if(ui->visualProgressCheckBox->isChecked())
        on_drawStuffButton_released();
    */
}

void MainWindow::on_addIterationsButton_released()
{
    /*
    bool graphicalProgress = ui->visualProgressCheckBox->isChecked();
    ui->progressBar->setValue(0);

    int additionalIterations = ui->addIterationsBox->text().toInt();
    maximumIterations = additionalIterations;
    int subdivisions =  10;

    if (iterationsPerColloidalEvolution >0)
        subdivisions = additionalIterations / iterationsPerColloidalEvolution;

    int stepsPerSubdivision = additionalIterations / subdivisions;
    vector<int3> moveChain;
    for (int ii = 0; ii < subdivisions; ++ii)
        {
        {
        auto upd = sim->updaters[0].lock();
        int curIterations = upd->getCurrentIterations();
        upd->setMaximumIterations(curIterations+stepsPerSubdivision);
        }
        sim->performTimestep();
        if(iterationsPerColloidalEvolution > 0)
            {
            for (int bb = 0; bb < Configuration->boundaryState.size();++bb)
                {
                if(Configuration->boundaryState[bb]==1)
                    landauLCForce->computeObjectForces(bb);
                }
            for (int bb = 0; bb < Configuration->boundaryState.size();++bb)
                {
                if(Configuration->boundaryState[bb]==1)
                    {
                    scalar3 currentForce = Configuration->boundaryForce[bb];
                    int dirx = -10;int diry = -10;int dirz = -10;
                    bool moved = false;
                    if(colloidalEvolutionPrefactor*currentForce.x > 1 || colloidalEvolutionPrefactor*currentForce.x < -1)
                        {
                        dirx = (currentForce.x > 0) ? 1 : 0;
                        currentForce.x = 1;moved = true;
                        Configuration->boundaryForce[bb].x = 0;
                        };
                    if(colloidalEvolutionPrefactor*currentForce.y > 1|| colloidalEvolutionPrefactor*currentForce.y < -1)
                        {
                        diry = (currentForce.y > 0) ? 3 : 2;
                        currentForce.y = 1;moved = true;
                        Configuration->boundaryForce[bb].y = 0;
                        };
                    if(colloidalEvolutionPrefactor*currentForce.z > 1|| colloidalEvolutionPrefactor*currentForce.z < -1)
                        {
                        dirz = (currentForce.z > 0) ? 5 : 4;
                        currentForce.z = 1;moved = true;
                        Configuration->boundaryForce[bb].z = 0;
                        };
                    if(moved)
                        {
                        int xm=0; int ym = 0; int zm = 0;
                        if(dirx >= 0)
                            {
                            Configuration->displaceBoundaryObject(bb, dirx,1);
                            xm = (dirx ==0 ) ? -1 : 1;
                            }
                        if(diry >= 0)
                            {
                            Configuration->displaceBoundaryObject(bb, diry,1);
                            ym = (diry ==2) ? -1 : 1;
                            }
                        if(dirz >= 0)
                            {
                            Configuration->displaceBoundaryObject(bb, dirz,1);
                            zm = (dirz ==4) ? -1 : 1;
                            }
                        int3 thisMove; thisMove.x = xm; thisMove.y=ym;thisMove.z=zm;
                        moveChain.push_back(thisMove);
                        }
                    }
                }
            }//end check of colloidal moves

        if(graphicalProgress) on_drawStuffButton_released();
        int progress = ((1.0*ii/(1.0*subdivisions))*100);
        QString printable2 = QStringLiteral("evolving... %1 percent done").arg(progress);
        ui->testingBox->setText(printable2);
        ui->progressBar->setValue(progress);
        }
    scalar maxForce = sim->getMaxForce();
    QString printable3 = QStringLiteral("system evolved...mean force is %1").arg(maxForce);
    ui->testingBox->setText(printable3);
    ui->progressBar->setValue(100);
    printf("move chain:\n");
    for(int ii = 0; ii < moveChain.size(); ++ii)
        printf("{%i,%i,%i},",moveChain[ii].x,moveChain[ii].y,moveChain[ii].z);
    printf("\nmove chain end :\n");
    */
}

void MainWindow::on_drawStuffButton_released()
{
    ArrayHandle<dVec> p(Configuration->returnPositions(),access_location::host,access_mode::read);
    ArrayHandle<dVec> v(Configuration->returnVelocities(),access_location::host,access_mode::read);
    scalar scale = ui->directorScaleBox->text().toDouble();
    vector<scalar3> lineSegments(2*N);
    vector<scalar3> defects;
    scalar3 director;
    QString printable1 = QStringLiteral("finding vectors ");
    ui->testingBox->setText(printable1);
    for (int ii = 0; ii < N; ++ii)
        {
        scalar3 pos;
        pos.x = p.data[ii].x[0];
        pos.y = p.data[ii].x[1];
        pos.z = p.data[ii].x[2];

        director.x=v.data[ii].x[0];
        director.y=v.data[ii].x[1];
        director.z=v.data[ii].x[2];

        scalar3 lineSegment1;
        scalar3 lineSegment2;

        lineSegment1.x = pos.x;
        lineSegment2.x = pos.x+0.5*scale*director.x;
        lineSegment1.y = pos.y;
        lineSegment2.y = pos.y+0.5*scale*director.y;
        lineSegment1.z = pos.z;
        lineSegment2.z = pos.z+0.5*scale*director.z;

        lineSegments.push_back(lineSegment1);
        lineSegments.push_back(lineSegment2);
    }
    int3 zero; zero.x=0;zero.y=0;zero.z=0;
    ui->displayZone->setLines(lineSegments,zero);
    ui->displayZone->update();
    /*


    bool goodVisualization = ui->builtinBoundaryVisualizationBox->isChecked();
    if(goodVisualization)
        {
        ui->displayZone->setSpheres(Configuration->latticeIndex.sizes);
        ui->displayZone->drawBoundaries = true;
        }
    else
        {
        on_builtinBoundaryVisualizationBox_released();
        };
    QString printable3 = QStringLiteral("drawing stuff ");
    ui->testingBox->setText(printable3);
    ui->displayZone->update();
    */
}

void MainWindow::on_xRotSlider_valueChanged(int value)
{
    ui->displayZone->setXRotation(value);
}

void MainWindow::on_zRotSlider_valueChanged(int value)
{
    ui->displayZone->setZRotation(value);
}

void MainWindow::on_zoomSlider_valueChanged(int value)
{
    zoom = value;
    ui->displayZone->zoom = zoom;
    //ui->displayZone->setLines(ui->displayZone->lines,Configuration->latticeIndex.sizes);
    //on_builtinBoundaryVisualizationBox_released();
    on_drawStuffButton_released();
}

/*
void MainWindow::on_addSphereButton_released()
{

    scalar3 spherePos;
    spherePos.x = ui->xSpherePosBox->text().toDouble()*BoxX;
    spherePos.y = ui->ySpherePosBox->text().toDouble()*BoxX;
    spherePos.z = ui->zSpherePosBox->text().toDouble()*BoxX;
    scalar rad = ui->sphereRadiusBox->text().toDouble()*BoxX;


    scalar W0 = ui->boundaryEnergyBox->text().toDouble();
    scalar s0b = ui->boundaryS0Box->text().toDouble();

    QString homeotropic ="homeotropic anchoring";
    QString planarDegenerate="planar degenerate anchoring";
    if(ui->anchoringComboBox->currentText() ==homeotropic)
        {
        boundaryObject homeotropicBoundary(boundaryType::homeotropic,W0,s0b);
        sim->createSphericalColloid(spherePos,rad,homeotropicBoundary);
        //Configuration->createSimpleSpherialColloid(spherePos,rad, homeotropicBoundary);
        }
    else if(ui->anchoringComboBox->currentText() ==planarDegenerate)
        {
        boundaryObject planarDegenerateBoundary(boundaryType::degeneratePlanar,W0,s0b);
        sim->createSphericalColloid(spherePos,rad,planarDegenerateBoundary);
        //Configuration->createSimpleSpherialColloid(spherePos,rad, planarDegenerateBoundary);
        }
    spherePositions.push_back(spherePos);
    sphereRadii.push_back(rad);
    QString printable1 = QStringLiteral("sphere added ");
    ui->testingBox->setText(printable1);
    ui->displayZone->addSphere(spherePos,rad);

}
  */

void MainWindow::on_actionReset_the_system_triggered()
{
    hideControls();
    ui->displayZone->clearObjects();
    ui->initializationFrame->show();
    QString printable1 = QStringLiteral("system reset");
    ui->testingBox->setText(printable1);
}

void MainWindow::on_reprodicbleRNGBox_stateChanged(int arg1)
{
    bool repro = ui->reprodicbleRNGBox->isChecked();
    noise.Reproducible= repro;
    if(repro)
        noise.setReproducibleSeed(13377);
    sim->setReproducible(repro);
}


void MainWindow::on_saveFileNowButton_released()
{
    /*
    QString fname = ui->saveFileNameBox->text();
    string fileName = fname.toStdString();

    ArrayHandle<dVec> pp(Configuration->returnPositions());
    ArrayHandle<int> tt(Configuration->returnTypes());
    ofstream myfile;
    myfile.open (fileName.c_str());
    for (int ii = 0; ii < Configuration->getNumberOfParticles();++ii)
        {
        int3 pos = Configuration->latticeIndex.inverseIndex(ii);
        myfile << pos.x <<"\t"<<pos.y<<"\t"<<pos.z;
        for (int dd = 0; dd <DIMENSION; ++dd)
            myfile <<"\t"<<pp.data[ii][dd];
        myfile << "\t"<<tt.data[ii]<<"\n";
        };

    sim->saveState("../data/saveTesting.txt");

    myfile.close();
    QString printable1 = QStringLiteral("File saved");
    ui->testingBox->setText(printable1);
    ui->fileSaveWidget->hide();
    */
}

void MainWindow::on_computeEnergyButton_released()
{
    /*
     ui->progressBar->setValue(0);
    landauLCForce->computeEnergy();
     ui->progressBar->setValue(90);
    scalar totalEnergy = 0.0;
    for(int ii = 0; ii < landauLCForce->energyComponents.size();++ii)
        totalEnergy+=landauLCForce->energyComponents[ii];
    QString energyString = QStringLiteral("Total energy: %1, components (phase, distortion, anchoring, E, H):  ").arg(totalEnergy);
    for(int ii = 0; ii < landauLCForce->energyComponents.size();++ii)
        energyString += QStringLiteral(" %1,  ").arg(landauLCForce->energyComponents[ii]);
    ui->testingBox->setText(energyString);
    ui->progressBar->setValue(100);
    */
}
