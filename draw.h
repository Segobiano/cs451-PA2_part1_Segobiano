//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma once

#include <GL/gli.h>
#include <time.h>

#include <map>
using namespace std;

#include "Basic.h"
#include "model.h"
using namespace mathtool;

//-----------------------------------------------------------------------------
//variables used in rendering

bool showWire=false; //on/off wireframe
bool showSharpEdgeOnly=false;
bool randomColor=false;
bool background=true; //on/off background
bool light=true; //on/off background

//Store display IDs and model colors
map<model*,int> model_solid_gids;
map<model*,int> model_wire_gids;
map<model*,Vector3d> model_colors;
vector<Point4d> SelectedNode;
static int Selectedx=-1;
static int Selectedy = -1;
static int Selectedz = -1;
double storedz = 0;
double oldx;
double oldy;
double oldz;
double deltay;
double deltax;

double deltaz;

//-----------------------------------------------------------------------------
//this defines the size of the lattice
extern unsigned int lattice_nx, lattice_ny, lattice_nz;
extern vector<Point3d> FFD_lattice; //This stores all lattice nodes, FFD_lattice has size = (lattice_nx X lattice_ny X lattice_nz)
extern vector<Point3d> FFD_parameterization; //This stores all parameterized coordinates of all vertices from the model
                                             //FFD_parameterization has size = models.front().v_size

//-----------------------------------------------------------------------------
static bool toggleTurn() { 

}
inline void DisplayLattice()
{
	glEnable(GL_LIGHTING);
	for (int i = 0; i < BBmatrix.size(); i++){
		for (int j = 0; j < BBmatrix[0].size(); j++){
			for (int h = 0; h < BBmatrix[0][0].size(); h++){
				//if the point is selected
				if (Selectedx == h && Selectedy == j && Selectedz == i){
					glPointSize(20.0);
					glBegin(GL_POINTS);
					glColor3f(0.0f, 1.0f, 0.0f);
					glVertex3f(BBmatrix[i][j][h][0], BBmatrix[i][j][h][1], BBmatrix[i][j][h][2]);
					glEnd();
					glColor3f(0.0f, 0.0f, 1.0f);
					glPointSize(8.0);

				}
				else{
					glPointSize(8.0);
					glBegin(GL_POINTS);
					glColor3f(0.95f, 0.207, 0.031f);
					glVertex3f(BBmatrix[i][j][h][0], BBmatrix[i][j][h][1], BBmatrix[i][j][h][2]);
					glEnd();
					glColor3f(0.0f, 0.0f, 1.0f);
				}
				//draws lines between points
				glBegin(GL_LINES);
				if (h != lattice_nx - 1){
					glVertex3f(BBmatrix[i][j][h][0], BBmatrix[i][j][h][1], BBmatrix[i][j][h][2]);
					glVertex3f(BBmatrix[i][j][h + 1][0], BBmatrix[i][j][h + 1][1], BBmatrix[i][j][h + 1][2]);
				}
				if (j != lattice_ny - 1){
					glVertex3f(BBmatrix[i][j][h][0], BBmatrix[i][j][h][1], BBmatrix[i][j][h][2]);
					glVertex3f(BBmatrix[i][j + 1][h][0], BBmatrix[i][j + 1][h][1], BBmatrix[i][j + 1][h][2]);
				}
				if (i != lattice_nz - 1){
					glVertex3f(BBmatrix[i][j][h][0], BBmatrix[i][j][h][1], BBmatrix[i][j][h][2]);
					glVertex3f(BBmatrix[i + 1][j][h][0], BBmatrix[i + 1][j][h][1], BBmatrix[i + 1][j][h][2]);
				}
				glEnd();
			}
		}
	}
	
	//glEnd();
	//TODO: draw lattice nodes using FFD_lattice


	//TODO: draw lattice edges using FFD_lattice

}

inline void DisplayModel(model& M, bool randcolor=false)
{
	//draw
	if (randcolor){
		if (model_colors.find(&M) == model_colors.end())
			model_colors[&M] = Vector3d(drand48() + 0.5, drand48() + 0.5, drand48(), +0.5).normalize() + Vector3d(0.25, 0.25, 0.25);
		glColor3dv(model_colors[&M].get());
	}
	
	//Draw facets
	glEnable( GL_POLYGON_OFFSET_FILL );
	glPolygonOffset( 0.5f, 0.5f );
	//for(list<polygon>::iterator i=M.polys.begin();i!=M.polys.end();i++)
	glBegin(GL_TRIANGLES);
	for(unsigned int i=0;i<M.t_size;i++)
	{
        const triangle & t=M.tris[i];
        glNormal3dv(M.tris[i].n.get());
        for(int k=0;k<3;k++)
		{
            const Point3d& pt=M.vertices[t.v[k]].p;

			//send pt to OpenGL
            glVertex3d(pt[0],pt[1],pt[2]);
        }
	}
	glEnd();
	glDisable( GL_POLYGON_OFFSET_FILL );
}

inline void DisplayModelWireFrame(model& M, bool randcolor=false)
{
    //Draw Edges
    if(showWire)
	{
		glBegin(GL_LINES);
        for(uint i=0;i<M.e_size;i++){
            glColor3f(0,0,0);
            const edge & e=M.edges[i];
            if(e.fid.size()==2){//normal case, check if e is sharp
                triangle& f1=M.tris[e.fid.front()];
                triangle& f2=M.tris[e.fid.back()];
                if(fabs(1-f1.n*f2.n)<1e-2){
                    if(showSharpEdgeOnly) continue; //not sharp
                    else
                        glColor3f(0.7f,0.7f,0.7f);
                }
            }

            Point3d& p1=M.vertices[e.vid[0]].p;
            Point3d& p2=M.vertices[e.vid[1]].p;
            glVertex3d(p1[0],p1[1],p1[2]);
            glVertex3d(p2[0],p2[1],p2[2]);
        }
        glEnd();
    }
}


//copied from meshlab
void DisplayBackground(void)
{
	float topcolor[]={1,1,1};
	float bottomcolor[]={1,1,0.5};
	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-1,1,-1,1,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glDisable(GL_TEXTURE_2D);
	glBegin(GL_TRIANGLE_STRIP);
		glColor3fv(topcolor);  	glVertex2f(-1, 1);
		glColor3fv(bottomcolor);	glVertex2f(-1,-1);
		glColor3fv(topcolor);	glVertex2f( 1, 1);
		glColor3fv(bottomcolor);	glVertex2f( 1,-1);
	glEnd();
	
	glPopAttrib();
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

void drawAll()
{
  
    //show the inputs
    glColor3f(1,1,1);
    if(light) glEnable(GL_LIGHTING);
    else glDisable(GL_LIGHTING);

    for(list<model>::iterator i=models.begin();i!=models.end();i++){
        DisplayModel(*i,randomColor);
    }
    for(list<model>::iterator i=models.begin();i!=models.end();i++){
        DisplayModelWireFrame(*i);
    }

	//draw lattice
	DisplayLattice();
}

//-----------------------------------------------------------------------------
void Display( void )
{
    //Init Draw
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    if(background) DisplayBackground();
    
    
    glPushMatrix();
    glLoadIdentity();
    static GLfloat light_position1[] = {  100, 100, 100.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
    static GLfloat light_position2[] = { -100, -100, 50.0f, 1.0f };
    glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
    glPopMatrix();

    drawAll();

    glDisable(GL_LIGHTING);
}


//-----------------------------------------------------------------------------
// regular openGL callback functions
bool InitGL()
{
    // transparent
    glShadeModel(GL_SMOOTH);
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

	glEnable( GL_LINE_SMOOTH );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );

    // others
    glEnable( GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
    glClearColor(1.0, 1.0, 1.0, 0.0);

    //Let's have light!
    GLfloat Diffuse[] =  { 0.9f, 0.9f, 0.9f, 1.0f };
    glMaterialfv(GL_FRONT, GL_DIFFUSE, Diffuse);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    GLfloat WhiteLight[] =  { 0.75f, 0.75f, 0.75f, 1.0f };
    glLightfv(GL_LIGHT0,GL_DIFFUSE,WhiteLight);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,WhiteLight);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    return true;
}

void Reshape( int w, int h)
{
    glViewport( 0, 0, (GLsizei)w, (GLsizei)h );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();

	//Perspective view
    gluPerspective( 60, 1.0f*w/h, R/100, R*100 );

	//Othogonal view
	//glOrtho(-R * 1.5, R * 1.5, -R * 1.5, R * 1.5, -R * 100, R * 100);


    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
}

void Mouse(int button, int state, int x, int y){
	//if the mouse press is down otherwise nothing is clicked
	if (!state){
		
		GLdouble* modelmax = new GLdouble[16];
		GLdouble* project = new GLdouble[16];
		GLint* view = new GLint[4];
		glGetDoublev(GL_MODELVIEW_MATRIX, modelmax);
		glGetDoublev(GL_PROJECTION_MATRIX, project);
		glGetIntegerv(GL_VIEWPORT, view);

		GLdouble* scx = new GLdouble(0);
		GLdouble* scy = new GLdouble(0);
		GLdouble* scz = new GLdouble(0);
		double maxdistance = 4000.0;
		double pointdistance = 4000.0;
		double numberholder = 0;
		//data structure that is 3d array of points
		for (int i = 0; i < BBmatrix.size(); i++){
			for (int j = 0; j < BBmatrix[0].size(); j++){
				for (int h = 0; h < BBmatrix[0][0].size(); h++){
				
					gluProject(BBmatrix[i][j][h][0], BBmatrix[i][j][h][1], BBmatrix[i][j][h][2], modelmax, project, view, scx, scy, scz);
					*scx = *scx;
					*scy = (view[3] - *scy) - 10;
					//find the closet node to the click
					pointdistance = sqrt((*scx - x)*(*scx - x) + (*scy - y)*(*scy - y));
					if ((pointdistance < maxdistance) && (pointdistance < 6)){
						gli::onSelect = true;
						maxdistance = pointdistance;
						Selectedx = h;
						Selectedy = j;
						Selectedz = i;
						storedz = *scz;
					}
				}
			}
		}
	}
	else {
		//unselect when not clicked
		Selectedx = -1;
		Selectedy = -1;
		Selectedz = -1;
		//tells gli that a node is not selected
		gli::onSelect = false;
	}

}


void Motion(int x, int y)
{
	//only if you have a selected node
	if (Selectedx != -1 && Selectedy != -1 && Selectedz != -1){
		GLdouble* modelmax = new GLdouble[16];
		GLdouble* project = new GLdouble[16];
		GLint* view = new GLint[4];
		glGetDoublev(GL_MODELVIEW_MATRIX, modelmax);
		glGetDoublev(GL_PROJECTION_MATRIX, project);
		glGetIntegerv(GL_VIEWPORT, view);
		GLdouble* scx = new GLdouble(0);
		GLdouble* scy = new GLdouble(0);
		GLdouble* scz = new GLdouble(0);
		//convert back to world coords
		// invert the y direction
		y = (view[3] - y) - 10;
		gluUnProject(x, y, storedz, modelmax, project, view, scx, scy, scz);
		Point3d tempPoint = Point3d(*scx, *scy, *scz);
		//set the nodes new location
		BBmatrix[Selectedz][Selectedy][Selectedx] = tempPoint;
		
	}

	model& m = models.front();
	glutPostRedisplay();
}

void PassiveMotion(int x, int y)
{
	// This handles mouse motion when mouse button is NOT pressed
	// does nothing now...

	glutPostRedisplay();
}

//Used for simulation/anitmation. 
void TimerCallback(int value)
{
    //in simuation state
    glutPostRedisplay();
    glutTimerFunc(30, TimerCallback,value);
}

//Handle keyboard events
void Keyboard( unsigned char key, int x, int y )
{
    // find closest colorPt3D if ctrl is pressed...
    switch( key ){
        case 27: exit(0);
        case 'w' : showWire=!showWire; break;
        case 'r' : randomColor=!randomColor; break;
		case 'R' : model_colors.clear(); break;
		case 'L' : light=!light; break;
		case 'b' : background=!background; break;
		case 'S' : showSharpEdgeOnly=!showSharpEdgeOnly;
		           for(map<model*,int>::iterator i=model_wire_gids.begin();i!=model_wire_gids.end();i++) glDeleteLists(i->second,1);
		           model_wire_gids.clear();
		           break;
    }
    glutPostRedisplay();
}



