#include "oglwidget.h"
#include <QOpenGLFunctions>


OGLWidget::OGLWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    xRot = 45;
    zRot = 135;
    scalar3 zero; zero.x = zero.y= zero.z=0;
//    baseSpherePositions.push_back(zero);
//    baseSphereRadii.push_back(1.);
    int3 one; one.x = one.y=one.z=1;
//    setSpheres(one);
}

OGLWidget::~OGLWidget()
{

}

void OGLWidget::clearObjects()
{
    baseSpherePositions.clear();
    baseSphereRadii.clear();
    spherePositions.clear();
    sphereRadii.clear();
    walls.clear();
}
void OGLWidget::initializeGL()
{
    glClearColor(0,0,0,1);
    glEnable(GL_DEPTH_TEST);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
}

void OGLWidget::setLines(vector<scalar3> &lineSegments, int3 sizes)
{
    float zm = 0.1*zoom;
    Sizes=sizes;
    lines =lineSegments;
    for (int ii = 0; ii < lines.size(); ++ii)
    {
        lines[ii].x = zm*((lines[ii].x-0.*sizes.x)/sizes.z);
        lines[ii].y = zm*((lines[ii].y-0.*sizes.y)/sizes.z);
        lines[ii].z = zm*((lines[ii].z-0.*sizes.z)/sizes.z);
    }
}

void OGLWidget::setConnections(vector<scalar3> &lineSegments, int3 sizes)
{
    float zm = 0.1*zoom;
    Sizes=sizes;
    connections =lineSegments;
    for (int ii = 0; ii < connections.size(); ++ii)
    {
        connections[ii].x = zm*((connections[ii].x-0.*sizes.x)/sizes.z);
        connections[ii].y = zm*((connections[ii].y-0.*sizes.y)/sizes.z);
        connections[ii].z = zm*((connections[ii].z-0.*sizes.z)/sizes.z);
    }
}

void OGLWidget::setDefects(vector<scalar3> &def, int3 sizes)
{
    Sizes=sizes;
    defects =def;
    float zm = 0.1*zoom;
    for (int ii = 0; ii < defects.size(); ++ii)
    {
        defects[ii].x = zm*((defects[ii].x-0.5*sizes.x)/sizes.z);
        defects[ii].y = zm*((defects[ii].y-0.5*sizes.y)/sizes.z);
        defects[ii].z = zm*((defects[ii].z-0.5*sizes.z)/sizes.z);
    }
}

void OGLWidget::setSpheres(int3 sizes)
{
    spherePositions.resize(baseSpherePositions.size());
    sphereRadii.resize(baseSpherePositions.size());
    float zm = 0.1*zoom;

    for (int ii = 0; ii < baseSpherePositions.size(); ++ii)
    {
        spherePositions[ii].x = zm*((baseSpherePositions[ii].x-0.*sizes.x)/sizes.z);
        spherePositions[ii].y = zm*((baseSpherePositions[ii].y-0.*sizes.y)/sizes.z);
        spherePositions[ii].z = zm*((baseSpherePositions[ii].z-0.*sizes.z)/sizes.z);
        sphereRadii[ii] = zm*baseSphereRadii[ii]/sizes.z;
    }
}

void OGLWidget::addSphere(scalar3 &pos,scalar &radii)
{
    baseSpherePositions.push_back(pos);
    baseSphereRadii.push_back(radii);
}

void OGLWidget::addWall(int3 planeAndNormalAndType)
{
    walls.push_back(planeAndNormalAndType);
}
void OGLWidget::draw()
{
    glBegin(GL_LINES);
    for (int ii = 0; ii < lines.size(); ++ii)
    {
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(lines[ii].x,lines[ii].y,lines[ii].z);
    }
    for (int ii = 0; ii < connections.size(); ++ii)
    {
        glColor4f(.0, 0.0, .9,.2);
        glVertex3f(connections[ii].x,connections[ii].y,connections[ii].z);
    }

    glEnd();
    float zm = 0.1*zoom;
    float empiricallyNiceRadius = 1./(pow(zm,0.25));
    for (int ii = 0; ii < defects.size(); ++ii)
    {
        glColor3f(0.0, 0.0, 1.0);
        GLUquadric *quad;
        quad = gluNewQuadric();
        glTranslatef(defects[ii].x,defects[ii].y,defects[ii].z);
        gluSphere(quad,empiricallyNiceRadius,20,20);
        glTranslatef(-defects[ii].x,-defects[ii].y,-defects[ii].z);
    }

    glEnd();
}

void OGLWidget::drawBoundarySites()
{
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    float empiricallyNiceRadius = 0.2/(pow(zoom,0.15));
    for (int ii = 0; ii < boundarySites.size(); ++ii)
        {
            glColor4f(0.,0.,0.8,0.2);
            GLUquadric *quad;
            quad = gluNewQuadric();
            glTranslatef(boundarySites[ii].x,boundarySites[ii].y,boundarySites[ii].z);
            gluSphere(quad,empiricallyNiceRadius,10,10);
            glTranslatef(-boundarySites[ii].x,-boundarySites[ii].y,-boundarySites[ii].z);
        }
        glDisable (GL_BLEND);
}

void OGLWidget::drawSpheres()
{
    glEnable (GL_BLEND);

    glBlendFunc (GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    for (int ii = 0; ii < spherePositions.size(); ++ii)
    {
        glColor4f(0.,0.,0.7,0.3);
        GLUquadric *quad;
        quad = gluNewQuadric();
        glTranslatef(spherePositions[ii].x,spherePositions[ii].y,spherePositions[ii].z);
        gluSphere(quad,sphereRadii[ii],100,20);
        glTranslatef(-spherePositions[ii].x,-spherePositions[ii].y,-spherePositions[ii].z);
    }
    glDisable (GL_BLEND);
}

void OGLWidget::drawWalls()
{

}

void OGLWidget::resizeGL(int w, int h)
{
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, (float)w/h, 0.01, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(5,5,0,0,0,0,0,0,1);
}

void OGLWidget::paintGL()
{
    scalar dirx,diry,dirz;

    scalar theta = zRot*1.0*PI/360.;
    scalar phi = xRot*2.0*PI/360.;
    dirx = sin(theta)*cos(phi);
    diry = sin(theta)*sin(phi);
    dirz = cos(theta);

    glLoadIdentity();
    gluLookAt(8*dirx,8*diry,8*dirz,0,0,0,0,0,1);
    if(drawBoundaries)
        {
        draw();
        drawSpheres();
        }
    else
        {
        draw();
        drawBoundarySites();
        }
}

void OGLWidget::setAllBoundarySites(vector<int3> &sites)
{
    boundarySites.resize(sites.size());
    float zm = ((float)zoom)*0.1;
    for(int ii = 0; ii < sites.size();++ii)
        {
        boundarySites[ii].x = zm*((sites[ii].x-0.5*Sizes.x)/Sizes.z);
        boundarySites[ii].y = zm*((sites[ii].y-0.5*Sizes.y)/Sizes.z);
        boundarySites[ii].z = zm*((sites[ii].z-0.5*Sizes.z)/Sizes.z);
        };
}

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360;
    while (angle > 360)
        angle -= 360;
}

void OGLWidget::setXRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(xRot);
        update();
    }
}

void OGLWidget::setZRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(zRot);
        update();
    }
}

void OGLWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}

void OGLWidget::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();

    if(event->buttons() & Qt::LeftButton)
    {
        setXRotation(xRot + 1.*dx);
        setZRotation(zRot + 1.*dy);
    }

    lastPos =event->pos();
}
