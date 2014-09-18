
//------------------------------------------------------------------------------
//  Copyright 2007-2014 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------


#pragma once

#include "objReader.h"
#include "model.h"
#include <list>
#include <float.h>
#include "gli.h"
using namespace std;

//-----------------------------------------------------------------------------
// INPUTS
list<string> input_filenames;

//this defines the size of the lattice
unsigned int lattice_nx = 4, lattice_ny = 3, lattice_nz = 3; 

//-----------------------------------------------------------------------------
// Intermediate data
list<model> models; //NOTE: only the first model of the list is used in this code
float R=0;          //radius
Point3d COM;        //center of mass

//TODO: fill FFD_lattice in computeCOM_R() below
vector<Point3d> FFD_lattice; //This stores all lattice nodes, FFD_lattice has size = (lattice_nx X lattice_ny X lattice_nz)
vector <vector <vector <Point3d > > > BBmatrix;//3d Array of Nodes or clicking points


//TODO: fill FFD_parameterization in parameterizeModel() below
vector<Point3d> FFD_parameterization; //This stores all parameterized coordinates of all vertices from the model
                                      //FFD_parameterization has size = models.front().v_size

//-----------------------------------------------------------------------------
bool readfromfile();
void computeCOM_R();
void parameterizeModel();

//-----------------------------------------------------------------------------
bool parseArg(int argc, char ** argv)
{
    for(int i=1;i<argc;i++){
        if(argv[i][0]=='-')
        {
			if (string(argv[i]) == "-nx")      lattice_nx = atoi(argv[++i]);
			else if (string(argv[i]) == "-ny") lattice_ny = atoi(argv[++i]);
			else if (string(argv[i]) == "-nz") lattice_nz = atoi(argv[++i]);
			else
				return false; //unknown
        }
        else{
            input_filenames.push_back(argv[i]);
        }
    }

    return true;
}

void printUsage(char * name)
{
    //int offset=20;
    cerr<<"Usage: "<<name<<" [options] -nx integer -ny integer -nz integer *.obj \n"
        <<"options:\n\n";
    cerr<<"\n-- Report bugs to: Jyh-Ming Lien jmlien@gmu.edu"<<endl;
}

//-----------------------------------------------------------------------------

bool readfromfiles()
{
	if (input_filenames.empty())
	{
		cerr << "! Error: No input model" << endl;
		return false;
	}

    long vsize=0;
    long fsize=0;

    uint id=0;
    for(list<string>::iterator i=input_filenames.begin();i!=input_filenames.end();i++,id++){
        cout<<"- ["<<id<<"/"<<input_filenames.size()<<"] Start reading "<<*i<<endl;
        model m;
        if(!m.build(*i)) continue;
        cout<<"- ["<<id<<"/"<<input_filenames.size()<<"] Done reading "<<m.v_size<<" vertices and "<<m.t_size<<" facets"<<endl;
        vsize+=m.v_size;
        fsize+=m.t_size;
        models.push_back(m);
    }
    cout<<"- Total: "<<vsize<<" vertices, "<<fsize<<" triangles, and "<<input_filenames.size()<<" models"<<endl;
    computeCOM_R();
	parameterizeModel();

    return true;
}

void computeCOM_R()
{
    //compute a bbox
    double box[6]={FLT_MAX,-FLT_MAX,FLT_MAX,-FLT_MAX,FLT_MAX,-FLT_MAX};
    //-------------------------------------------------------------------------
    for(list<model>::iterator i=models.begin();i!=models.end();i++){
        for(unsigned int j=0;j<i->v_size;j++){
            Point3d& p=i->vertices[j].p;
            if(p[0]<box[0]) box[0]=p[0];
            if(p[0]>box[1]) box[1]=p[0];
            if(p[1]<box[2]) box[2]=p[1];
            if(p[1]>box[3]) box[3]=p[1];
            if(p[2]<box[4]) box[4]=p[2];
            if(p[2]>box[5]) box[5]=p[2];
        }//j
    }//i
	//scales bounding box
	box[0] = box[0] * 1.2;
	box[1] = box[1] * 1.2;
	box[2] = box[2] * 1.2;
	box[3] = box[3] * 1.2;
	box[4] = box[4] * 1.2;
	box[5] = box[5] * 1.2;
	//calculates the armount of space between nodes
	double xspace = (abs(box[1] - box[0])) / (lattice_nx - 1);
	double yspace = (abs(box[3] - box[2])) / (lattice_ny - 1);
	double zspace = (abs(box[5] - box[4])) / (lattice_nz - 1);
	double tempx = box[0];
	double tempy = box[2];
	double tempz = box[4];
	

	for (int i = 0; i < lattice_nz; i++){
		//holds planes
		vector < vector <Point3d >> tempplane;
		for (int j = 0; j < lattice_ny; j++){
			//holds rows
			vector <Point3d > temprow;
			for (int h = 0; h < lattice_nx; h++){
				//fills row
				temprow.push_back(Point3d(tempx, tempy, tempz));
				FFD_lattice.push_back(Vector3d(tempx, tempy, tempz));
				tempx = tempx + xspace;
				
			}
			//fills planes
			tempplane.push_back(temprow);
			temprow.clear();
			tempx = box[0];
			tempy = tempy + yspace;
		}
		//fill matrix layers
		BBmatrix.push_back(tempplane);
		tempplane.clear();
		tempy = box[2];
		tempz = tempz + zspace;
	}


	//TODO: scale the bounding box "box[6]" 1.2 times 
	//TODO: build FFD_lattice here using the scaled bounding box

    //-------------------------------------------------------------------------
    // compute center of mass and R...
    COM.set( (box[1]+box[0])/2,(box[3]+box[2])/2,(box[5]+box[4])/2);

    //-------------------------------------------------------------------------
	R=0;
    for(list<model>::iterator i=models.begin();i!=models.end();i++){
        for(unsigned int j=0;j<i->v_size;j++){
            Point3d& p=i->vertices[j].p;
            float d=(float)(p-COM).normsqr();
            if(d>R) R=d;
        }//j
    }//i

    R=sqrt(R);
}

//convert each vertex of the model into parameterized space
//and store the parameterization in FFD_parameterization
void parameterizeModel()
{
	model& m = models.front();

	for (unsigned int i = 0; i < m.v_size; i++)
	{
		vertex & v = m.vertices[i];
		
		//TODO: convert v.p (the position of p) into the parameterized space using FFD_lattice
		Point3d pc; //parameterized coordinate

		FFD_parameterization.push_back(pc);
	}
	//end i

	//done
}

//-----------------------------------------------------------------------------
//
//
//
//  Open GL stuff below
//
//
//-----------------------------------------------------------------------------

#include <draw.h>



