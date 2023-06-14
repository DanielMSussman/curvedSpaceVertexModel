#include <QApplication>
#include <QMainWindow>
#include <QSplashScreen>
#include <QDesktopWidget>
#include <QTimer>
#include <QGuiApplication>

#include <QPropertyAnimation>
#include "mainwindow.h"
#include <tclap/CmdLine.h>


using namespace TCLAP;
int main(int argc, char*argv[])
{

    QApplication a(argc, argv);

    QSplashScreen *splash = new QSplashScreen;

    QRect screenGeometry = QApplication::desktop()->screenGeometry();


    splash->setPixmap(QPixmap("../examples/splashWithText.jpeg"));
    splash->move(0,0*screenGeometry.height() / 2);
    splash->show();

    MainWindow w;
    int x = 0*(screenGeometry.width()-w.width())/2;
    int y = (screenGeometry.height()-w.height())/2;

    w.move(x,y);

    QTimer::singleShot(750,splash,SLOT(close()));

    QTimer::singleShot(750,&w,SLOT(show()));

    return a.exec();
};
