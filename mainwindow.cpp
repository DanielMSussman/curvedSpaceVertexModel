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
    connect(ui->displayZone,SIGNAL(xRotationChanged(int)),ui->xRotSlider,SLOT(setValue(int)));
    connect(ui->displayZone,SIGNAL(zRotationChanged(int)),ui->zRotSlider,SLOT(setValue(int)));

    hideControls();
    QString printable = QStringLiteral("Welcome!");
    ui->testingBox->setText(printable);
}

void MainWindow::hideControls()
{
    ui->evolutionParametersWidget->hide();
}
void MainWindow::showControls()
{
}

void MainWindow::on_initializeButton_released()
    {
    N = ui->boxNTotalSize->text().toInt();

    //eta =ui->initialEta->text().toDouble();
    dt =ui->initialDt->text().toDouble();

    //ui->initialEtaSet->setText(ui->initialEta->text());
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

void MainWindow::on_boxRadius_textEdited(const QString &arg1)
{
    N = ui->boxNTotalSize->text().toInt();
    radius = ui->boxRadius->text().toDouble();
        density = ((double) N) / (4.0*PI*radius*radius) ;
        QString valueAsString = QString::number(density);
        ui->boxDensity->setText(valueAsString);
}

void MainWindow::on_boxDensity_textEdited(const QString &arg1)
{
    N = ui->boxNTotalSize->text().toInt();
    density = ui->boxDensity->text().toDouble();
        radius = 0.5*sqrt((double)N/(density*PI));
        QString valueAsString = QString::number(radius);
        ui->boxRadius->setText(valueAsString);
}



void MainWindow::on_boxNTotalSize_textChanged(const QString &arg1)
{
    N = ui->boxNTotalSize->text().toInt();
    radius = sqrt(N/(4.0*PI));
    QString valueAsString2 = QString::number(N/(4.0*PI*radius*radius));
    ui->boxDensity->setText(valueAsString2);
    QString valueAsString = QString::number(radius);
    ui->boxRadius->setText(valueAsString);
}

void MainWindow::on_setParametersButton_released()
{
    //v0 =ui->initialSpeedSet->text().toDouble();
    //eta =ui->initialEtaSet->text().toDouble();
    dt =ui->initialDtSet->text().toDouble();
    //vicsek->setEta(eta);
    //vicsek->setV0(v0);
    //vicsek->setDeltaT(dt);
    ui->evolutionParametersWidget->hide();
    //ui->initialSpeed->setText(ui->initialSpeedSet->text());
    //ui->initialEta->setText(ui->initialEtaSet->text());
    ui->initialDt->setText(ui->initialDtSet->text());

}

void MainWindow::simulationInitialize()
{
    scalar preferredP = ui->p0Box->text().toDouble();
    scalar preferredA = ui->a0Box->text().toDouble();
    scalar initialT = ui->initialT->text().toDouble();
    Configuration = make_shared<sphericalVertexModel>(N,noise,preferredA,preferredP,false,true);
    scalar temperature = initialT;
    BD = make_shared<brownianDynamics>();
    BD->setT(temperature);
    N = Configuration->getNumberOfParticles();

    sim = make_shared<Simulation>();
    sim->setConfiguration(Configuration);
    sim->addUpdater(BD,Configuration);
    sim->setIntegrationTimestep(dt);
    sim->setCPUOperation(true);
    
    scalar3 zero; zero.x = zero.y= zero.z=0;
    int3 one; one.x = one.y=one.z=1;
    scalar rad = 1.0;
    ui->displayZone->addSphere(zero,rad);
    ui->displayZone->setSpheres(one);
}

void MainWindow::on_resetSystemButton_released()
{
    Configuration->setParticlePositionsRandomly(noise);
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
    if(subdivisions ==0)
        {
        prof2.start();
        for (int ii = 0; ii < additionalIterations;++ii)
           sim->performTimestep();
        prof2.end();
        on_drawStuffButton_released();
        }
    prof2.print();
    prof1.print();

    QString printable3 = QStringLiteral("system evolved...");
    ui->testingBox->setText(printable3);
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
#if DIMENSION == 3
        pos.z = p.data[ii].x[2]/radius;
        director.z=v.data[ii].x[2];
#else
        pos.z = 0.;
        director.z=0.;
#endif
        director.x=v.data[ii].x[0];
        director.y=v.data[ii].x[1];

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
                    end.x = p.data[neighbor][0]/radius;
                    end.y = p.data[neighbor][1]/radius;
                    start.z = p.data[ii][2]/radius;
                    end.z = p.data[neighbor][2]/radius;
                    scalar len = norm(p.data[neighbor]-p.data[ii]);
                    connections.push_back(start);
                    connections.push_back(end);
                    };

                }
            }

        }
    ui->displayZone->setConnections(connections,zero);

    ui->displayZone->update();
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
}


void MainWindow::on_computeEnergyButton_2_released()
{
    ui->evolutionParametersWidget->show();
}

void MainWindow::on_cancelEvolutionParametersButton_pressed()
{
    ui->evolutionParametersWidget->hide();
}

