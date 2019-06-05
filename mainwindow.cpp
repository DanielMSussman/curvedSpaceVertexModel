#include <QMainWindow>
#include <QGuiApplication>
#include <QPropertyAnimation>
#include <chrono>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "profiler.h"
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
    ui->evolutionParametersWidget->hide();
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

    ui->initialSpeedSet->setText(ui->initialSpeed->text());
    ui->initialEtaSet->setText(ui->initialEta->text());
    ui->initialDtSet->setText(ui->initialDt->text());

    radius = ui->boxRadius->text().toDouble();
    density = ui->boxDensity->text().toDouble();
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

void MainWindow::on_boxNTotalSize_textChanged(const QString &arg1)
{
    N = ui->boxNTotalSize->text().toInt();
    density = ui->boxDensity->text().toDouble();
    radius = 0.5*sqrt((double)N/(density*PI));
    QString valueAsString = QString::number(radius);
    ui->boxRadius->setText(valueAsString);
};


void MainWindow::on_boxRadius_textEdited(const QString &arg1)
{
    N = ui->boxNTotalSize->text().toInt();
    radius = ui->boxRadius->text().toDouble();
    density = ((double) N) / (4.0*PI*radius*radius) ;
    QString valueAsString = QString::number(density);
    ui->boxDensity->setText(valueAsString);
    QString valueAsString2 = QString::number(0.25*PI*N/(4.0*PI*radius*radius));
    ui->boxPackingFraction->setText(valueAsString2);
}

void MainWindow::on_boxDensity_textEdited(const QString &arg1)
{
    N = ui->boxNTotalSize->text().toInt();
    density = ui->boxDensity->text().toDouble();
    radius = 0.5*sqrt((double)N/(density*PI));
    QString valueAsString = QString::number(radius);
    ui->boxRadius->setText(valueAsString);
    QString valueAsString2 = QString::number(0.25*PI*N/(4.0*PI*radius*radius));
    ui->boxPackingFraction->setText(valueAsString2);
}

void MainWindow::on_boxPackingFraction_textEdited(const QString &arg1)
{
     N = ui->boxNTotalSize->text().toInt();
     scalar phi =  ui->boxPackingFraction->text().toDouble();
     radius = sqrt(0.0625*N/phi);
     QString valueAsString = QString::number(radius);
     ui->boxRadius->setText(valueAsString);
     QString valueAsString2 = QString::number(N/(4.0*PI*radius*radius));
     ui->boxDensity->setText(valueAsString2);
}


void MainWindow::on_setParametersButton_released()
{
    v0 =ui->initialSpeedSet->text().toDouble();
    eta =ui->initialEtaSet->text().toDouble();
    dt =ui->initialDtSet->text().toDouble();
    vicsek->setEta(eta);
    vicsek->setV0(v0);
    vicsek->setDeltaT(dt);
    ui->evolutionParametersWidget->hide();
    ui->initialSpeed->setText(ui->initialSpeedSet->text());
    ui->initialEta->setText(ui->initialEtaSet->text());
    ui->initialDt->setText(ui->initialDtSet->text());

}

void MainWindow::simulationInitialize()
{
    if(ui->sphericalModel->isChecked())
        {
        if(ui->topologicalModel->isChecked())
            Configuration = make_shared<sphericalVoronoi>(N,noise);
        else
            Configuration = make_shared<sphericalModel>(N,noise);
        Configuration->setRadius(radius);
        Configuration->getNeighbors();
        if(ui->softRepulsion->isChecked())
            Configuration->setSoftRepulsion();
        vicsek = make_shared<sphericalVectorialVicsek>();
        }
    else
        {
        Configuration = make_shared<simpleModel>(N,noise);
        Configuration->setRadius(radius);
        Configuration->getNeighbors();
        if(ui->softRepulsion->isChecked())
            Configuration->setSoftRepulsion();
        vicsek = make_shared<vectorialVicsek>();
        }
    printf("%f %f %f\n",eta,v0,dt);
    vicsek->setEta(eta);
    vicsek->setV0(v0);
    vicsek->setDeltaT(dt);
    sim = make_shared<Simulation>();
    sim->setConfiguration(Configuration);
    sim->addUpdater(vicsek,Configuration);
    /*

     landauLCForce = make_shared<landauDeGennesLC>();

     landauLCForce->setPhaseConstants(A,B,C);
     landauLCForce->setModel(Configuration);
     sim->addForce(landauLCForce);
     on_fireParamButton_released();
     ui->reproducibleButton->setEnabled(true);
     */
}

void MainWindow::on_resetSystemButton_released()
{
    Configuration->setParticlePositionsRandomly(noise);
    Configuration->getNeighbors();
    on_drawStuffButton_released();
}

void MainWindow::on_resetSystemBandButton_released()
{
    int angularExtent = ui->angularExtent->value();
    scalar ae = angularExtent*0.01*0.5;
    Configuration->setParticlePositionsBandedRandomly(noise,ae);
    Configuration->getNeighbors();
    on_drawStuffButton_released();
}

void MainWindow::on_addIterationsButton_released()
{
    bool graphicalProgress = ui->visualProgressCheckBox->isChecked();
    ui->progressBar->setValue(0);

    int additionalIterations = ui->addIterationsBox->text().toInt();

    int stepsPerSubdivision = 1 / dt;
    int subdivisions = additionalIterations/stepsPerSubdivision;
    profiler prof1("drawing");
    profiler prof2("evolving");
    for (int ii = 0; ii < subdivisions; ++ii)
        {
        for (int jj = 0; jj < stepsPerSubdivision; ++jj)
            {
            prof2.start();
            sim->performTimestep();
            prof2.end();
            }
        if(graphicalProgress)
            {
            prof1.start();
            on_drawStuffButton_released();
            prof1.end();
            }
        int progress = ((1.0*ii/(1.0*subdivisions))*100);
        QString printable2 = QStringLiteral("evolving... %1 percent done").arg(progress);
        ui->testingBox->setText(printable2);
        ui->progressBar->setValue(progress);
        }
    prof2.print();
    prof1.print();

    QString printable3 = QStringLiteral("system evolved...");
    ui->testingBox->setText(printable3);
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
    ArrayHandle<dVec> v(Configuration->returnDirectors(),access_location::host,access_mode::read);
    scalar scale = ui->directorScaleBox->text().toDouble();
    vector<scalar3> lineSegments;
    vector<scalar3> defects;
    scalar3 director;
    QString printable1 = QStringLiteral("finding vectors ");
    ui->testingBox->setText(printable1);
    double factor = sqrt(2./N)/pow(radius,.25);
    for (int ii = 0; ii < N; ++ii)
        {
        scalar3 pos;
        pos.x = p.data[ii].x[0]/radius;
        pos.y = p.data[ii].x[1]/radius;
        pos.z = p.data[ii].x[2]/radius;

        director.x=v.data[ii].x[0];
        director.y=v.data[ii].x[1];
        director.z=v.data[ii].x[2];

        scalar3 lineSegment1;
        scalar3 lineSegment2;

        lineSegment1.x = pos.x;
        lineSegment2.x = pos.x+factor*scale*director.x;
        lineSegment1.y = pos.y;
        lineSegment2.y = pos.y+factor*scale*director.y;
        lineSegment1.z = pos.z;
        lineSegment2.z = pos.z+factor*scale*director.z;

        lineSegments.push_back(lineSegment1);
        lineSegments.push_back(lineSegment2);
    }
    int3 zero; zero.x=1;zero.y=1;zero.z=1;
    ui->displayZone->setLines(lineSegments,zero);

    ArrayHandle<unsigned int> nNeighs(Configuration->numberOfNeighbors);
    ArrayHandle<int> neighs(Configuration->neighbors);
    vector<scalar3> connections;
    if(ui->showNeighborsCheckbox->isChecked())
        {
        for (int ii = 0; ii < N; ++ii)
            {
            for (int jj = 0; jj < nNeighs.data[ii]; ++jj)
                {
                int neighbor = neighs.data[Configuration->neighborIndex(jj,ii)]; 
                if(neighbor > ii)
                    {
                    scalar3 start,end;
                    start.x = p.data[ii][0]/radius;
                    start.y = p.data[ii][1]/radius;
                    start.z = p.data[ii][2]/radius;
                    end.x = p.data[neighbor][0]/radius;
                    end.y = p.data[neighbor][1]/radius;
                    end.z = p.data[neighbor][2]/radius;
                    if(!ui->sphericalModel->isChecked())
                        {
                        scalar len2 = radius*((start.x-end.x)*(start.x-end.x) + (start.y-end.y)*(start.y-end.y)+(start.z-end.z)*(start.z-end.z));
                        if(len2 < radius)
                            {
                            connections.push_back(start);
                            connections.push_back(end);
                            };
                        }
                    else
                        {
                        connections.push_back(start);
                        connections.push_back(end);
                        }
                    }
                }
            }

        }
    ui->displayZone->setConnections(connections,zero);

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
    QString energyString = QStringLiteral("This button doesn't do anything at the moment. But thanks for checking");
    ui->testingBox->setText(energyString);
    ui->progressBar->setValue(100);
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




void MainWindow::on_computeEnergyButton_2_released()
{
    ui->evolutionParametersWidget->show();
}

void MainWindow::on_cancelEvolutionParametersButton_pressed()
{
    ui->evolutionParametersWidget->hide();
}




