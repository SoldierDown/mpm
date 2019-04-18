#include <iostream>
#include <sys/stat.h> 
#include <string.h>
#include <GL/glut.h>
#include "MPM_Solver.h"


using namespace std;

void timer(int i_);
int time_ = 0;
float x = 0;
float y = 0;

string dir_name("boundary_");

MPM_Solver mpm;

void Screendump(char *tga_file, short W, short H) 
{
    FILE   *out = fopen(tga_file, "w");
    char   pixel_data[3*W*H];
    short  TGAhead[] = {0, 2, 0, 0, 0, 0, W, H, 24};
    
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, W, H, GL_BGR, GL_UNSIGNED_BYTE, pixel_data);
    fwrite(TGAhead, sizeof(TGAhead), 1, out);
    fwrite(pixel_data, 3*W*H, 1, out);
    fclose(out); 
}
void init (void)
{
    glClear (GL_COLOR_BUFFER_BIT);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity();
    glOrtho(OFF_SET, WIDTH, OFF_SET, WIDTH, -1.0, 1.0);
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();
}


 void drawSquare()
{
    double density = 1.0;
    double radius = 2.0;
    glClear(GL_COLOR_BUFFER_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    double x, y;
    for(auto p : mpm.particles)
    {
        x = WIDTH * p.pos(0);
        y = WIDTH * p.pos(1);
        glColor3f(density,density,density);
        glBegin(GL_POLYGON);
            glVertex3f(x + radius, y + radius, 0.0);
            glVertex3f(x + radius, y - radius, 0.0);
            glVertex3f(x - radius, y - radius, 0.0);
            glVertex3f(x - radius, y + radius, 0.0);
        glEnd();
    }

    glFlush();
  }

void display()
{   
    mpm.Step();
    if(time_ % int(mpm.frame_dt/mpm.dt) == 0) 
    {
        char path[256];
        sprintf(path, "./%s/Frame%05d.tga", dir_name.c_str(), time_);
        Screendump(path, WIDTH, WIDTH);
        drawSquare();
    }
    glFlush();
    glutSwapBuffers();
}

int main (int argc, char **argv)
{
    mkdir(dir_name.c_str(), 0777);
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize (WIDTH, WIDTH);
    glutCreateWindow ("MPM");

    init();
    glutDisplayFunc (display);
    glutTimerFunc(1000.0/FRAME, timer, 0);
    glDisable (GL_TEXTURE_2D);
    glutMainLoop();
    return 0;
}

void timer(int i_)
{
    time_++;
    glutPostRedisplay();
    glutTimerFunc(1000.0 / FRAME, timer, 0);
}