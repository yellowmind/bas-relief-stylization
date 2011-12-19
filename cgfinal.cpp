#include <math.h>
#include <stdlib.h>

#include "GLee.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <cv.h>
//#include <cvaux.h>
#include <highgui.h>

#include "swgl.h"
#include "math3d.h"
#include "glm.h"

#include <iostream>
#include <vector>

#ifndef bool
#define bool int
#define false 0
#define true 1
#endif

float maxf(float a, float b)
{
    return (a < b) ? b : a;
} 
float minf(float a, float b)
{
    return (a < b) ? a : b;
}
  
#ifndef M_PI
#define M_PI 3.14159
#endif
using namespace std;
//using namespace cv;

int 	winWidth, winHeight;

float 	angle = 0.0, axis[3], trans[3];
bool 	trackingMouse = false;
bool 	redrawContinue = false;
bool    trackballMove = false;
GLdouble TRACKM[16]={1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

GLdouble DEBUG_M[16];

GLdouble Angle1=0, Angle2=0;
GLint TICK=0;

GLMmodel *MODEL;
GLfloat lightPos0[] = { -15.f, 15.0f, 15.0, 1.0f };
GLfloat lightPos1[] = { -25.f, 15.0f, 7.5, 1.0f };
GLfloat lightPos2[] = { -27.f, 15.0f, 7.5, 1.0f };

GLdouble MODELSCALE = 1.0;
GLdouble LIGHTP = 15;

bool relief=true,relief2=true;
GLdouble angleX = 0,angleY = 0,angleZ = 0,scale =2.5;
GLdouble reliefAngleX = 0, reliefAngleY = 0, reliefAngleZ = 0;
GLdouble outputHeight = 0.1	;//0.075;
GLfloat threshold = 0.02;
vector< vector<GLfloat> > heightList, laplaceList;
const int pyrLevel = 4;
vector< vector<GLfloat > > heightPyr(pyrLevel);
CvMat **imgPyr;
//vector<GLfloat> height;
vector<GLfloat> equalizeHeight;
int boundary =20;
int disp=0;

GLdouble *pThreadRelief = NULL;
GLdouble *pThreadNormal = NULL;
GLdouble *pThreadEqualizeRelief = NULL;
GLdouble *pThreadEqualizeNormal = NULL;
GLint vertCount = 1;

int DRAWTYPE = 1;// 0:hw1, 1:hw2, 2:Gouraud shading, 3: Phong Shading

/*----------------------------------------------------------------------*/
/*
** Draw the wireflame cube.
*/
GLfloat vertices[][3] = {
    {-1.0,-1.0,-1.0},{1.0,-1.0,-1.0}, {1.0,1.0,-1.0}, {-1.0,1.0,-1.0}, 
    {-1.0,-1.0,1.0}, {1.0,-1.0,1.0}, {1.0,1.0,1.0}, {-1.0,1.0,1.0}
};

GLfloat colors[][3] = {
    {0.0,0.0,0.0},{1.0,0.0,0.0}, {1.0,1.0,0.0}, {0.0,1.0,0.0}, 
    {0.0,0.0,1.0}, {1.0,0.0,1.0}, {1.0,1.0,1.0}, {0.0,1.0,1.0}
};


inline void SwglTri(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble nx1=1, GLdouble ny1=0, GLdouble nz1=0, 
			 GLdouble nx2=1, GLdouble ny2=0, GLdouble nz2=0, 
			 GLdouble nx3=1, GLdouble ny3=0, GLdouble nz3=0,
			 GLdouble r1=1, GLdouble g1=1, GLdouble b1=1,
			 GLdouble r2=1, GLdouble g2=1, GLdouble b2=1,
			 GLdouble r3=1, GLdouble g3=1, GLdouble b3=1)
{
	//copy to homogenous coordinate
	GLdouble h1[4]={x1, y1, z1, 1.0};
	GLdouble h2[4]={x2, y2, z2, 1.0};
	GLdouble h3[4]={x3, y3, z3, 1.0};

	//window coordinate
	GLdouble w1[4]={x1, y1, 0, 1.0}; 
	GLdouble w2[4]={x2, y2, 0, 1.0};
	GLdouble w3[4]={x3, y3, 0, 1.0};

	//implement the opengl pipeline here
	swTransformation(h1, w1);
	swTransformation(h2, w2);
	swTransformation(h3, w3);

	switch(DRAWTYPE) {
		case 0:
		{
			//copy to homogenous coordinate
			GLdouble h1[4]={x1, y1, z1, 1.0};
			GLdouble h2[4]={x2, y2, z2, 1.0};
			GLdouble h3[4]={x3, y3, z3, 1.0};

			//window coordinate
			GLdouble w1[4]={x1, y1, 0, 1.0}; 
			GLdouble w2[4]={x2, y2, 0, 1.0};
			GLdouble w3[4]={x3, y3, 0, 1.0};

			//implement the opengl pipeline here
			swTransformation(h1, w1);
			swTransformation(h2, w2);
			swTransformation(h3, w3);

			writepixel(w1[0], w1[1], r1, g1, b1);
			writepixel(w2[0], w2[1], r2, g2, b2);
			writepixel(w3[0], w3[1], r3, g3, b3);
		}
		break;

		case 1:
		{
			//copy to homogenous coordinate
			GLdouble h1[4]={x1, y1, z1, 1.0};
			GLdouble h2[4]={x2, y2, z2, 1.0};
			GLdouble h3[4]={x3, y3, z3, 1.0};

			//window coordinate
			GLdouble w1[4]={x1, y1, 0, 1.0}; 
			GLdouble w2[4]={x2, y2, 0, 1.0};
			GLdouble w3[4]={x3, y3, 0, 1.0};

			//implement the opengl pipeline here
			swTransformation(h1, w1);
			swTransformation(h2, w2);
			swTransformation(h3, w3);

			swTriangle(w1[0], w1[1], w1[2],
					   w2[0], w2[1], w2[2],
					   w3[0], w3[1], w3[2],
					   r1, g1, b1);
		}
		break;

		case 2:
		{
			swTriangleG(x1, y1, z1,
						x2, y2, z2,
						x3, y3, z3,
						nx1, ny1, nz1,
						nx2, ny2, nz2,
						nx3, ny3, nz3,
						r1, g1, b1,
						r2, g2, b2,
						r3, g3, b3);
		}
		break;

		case 3:
		{
			swTriangleP(x1, y1, z1,
						x2, y2, z2,
						x3, y3, z3,
						nx1, ny1, nz1,
						nx2, ny2, nz2,
						nx3, ny3, nz3,
						r1, g1, b1,
						r2, g2, b2,
						r3, g3, b3);
		}
		break;

	}

}

void SwglTri(int index1, int index2, int index3)
{
	SwglTri( vertices[index1][0], vertices[index1][1], vertices[index1][2],
		     vertices[index2][0], vertices[index2][1], vertices[index2][2],
			 vertices[index3][0], vertices[index3][1], vertices[index3][2]);
}


void SwglLine(GLdouble x1, GLdouble y1, GLdouble z1, GLdouble x2, GLdouble y2, GLdouble z2)
{
	//copy to homogenous coordinate
	GLdouble h1[4]={x1, y1, z1, 1.0};
	GLdouble h2[4]={x2, y2, z2, 1.0};
	
	GLdouble w1[4]={x1, y1, 0, 1.0}; //window coordinate
	GLdouble w2[4]={x2, y2, 0, 1.0};

	//implement the opengl pipeline here
	swTransformation(h1, w1);
	swTransformation(h2, w2);
	
	////draw the 2D line
	//glBegin(GL_LINES);
	//	//glColor3fv(colors[index1]);
	//	glVertex2f(w1[0], w1[1]);
	//	//glColor3fv(colors[index2]);
	//	glVertex2f(w2[0], w2[1]);
	//glEnd();

	//implement 
	switch(DRAWTYPE) {
		case 0:
		{
			writepixel(w1[0], w1[1], 1, 0, 0);
			writepixel(w2[0], w2[1], 0, 1, 0);
		}
		break;

		case 1: case 2: case 3:
		{
			GLdouble col[4];
			glGetDoublev(GL_CURRENT_COLOR, col);
			BresenhamLine(w1[0], w1[1], w2[0], w2[1], col[0], col[1], col[2]);
		}
		break;
	}
}

void SwglLine(int index1, int index2)
{
	SwglLine(vertices[index1][0], vertices[index1][1], vertices[index1][2],
		     vertices[index2][0], vertices[index2][1], vertices[index2][2]);
}


void SolidQuad(int a, int b, int c, int d, bool USINGOPENGL)
{
	if(USINGOPENGL) {
		glBegin(GL_TRIANGLES);
			glVertex3fv(vertices[a]);
			glVertex3fv(vertices[b]);
			glVertex3fv(vertices[c]);

			glVertex3fv(vertices[c]);
			glVertex3fv(vertices[d]);
			glVertex3fv(vertices[a]);
		glEnd();
	} else {
		SwglTri(a, b, c);
		SwglTri(c, d, a);
	}
}

void swSolidCube(void)
{
    // map vertices to faces */
    SolidQuad(1,0,3,2, false);
    SolidQuad(3,7,6,2, false);
    SolidQuad(7,3,0,4, false);
    SolidQuad(2,6,5,1, false);
    SolidQuad(4,5,6,7, false);
    SolidQuad(5,4,0,1, false);
}

void glSolidCube(void)
{
    // map vertices to faces */
    SolidQuad(1,0,3,2, true);
    SolidQuad(3,7,6,2, true);
    SolidQuad(7,3,0,4, true);
    SolidQuad(2,6,5,1, true);
    SolidQuad(4,5,6,7, true);
    SolidQuad(5,4,0,1, true);
}



void OpenglLine(GLdouble x1, GLdouble y1, GLdouble z1, GLdouble x2, GLdouble y2, GLdouble z2)
{
	glBegin(GL_LINES);
		glVertex3f(x1, y1, z1);
		glVertex3f(x2, y2, z2);
	glEnd();
}

void OpenglLine(int index1, int index2)
{
	OpenglLine(vertices[index1][0], vertices[index1][1], vertices[index1][2],
		       vertices[index2][0], vertices[index2][1], vertices[index2][2]);
}

void WireQuad(int a, int b, int c , int d, bool USINGOPENGL)
{
	if(USINGOPENGL) {
		OpenglLine(a, b);
		OpenglLine(b, c);
		OpenglLine(c, d);
		OpenglLine(d, a);
	} else {
		SwglLine(a, b);
		SwglLine(b, c);
		SwglLine(c, d);
		SwglLine(d, a);
	}
}

void swWireCube(void)
{
    // map vertices to faces */
    WireQuad(1,0,3,2, false);
    WireQuad(3,7,6,2, false);
    WireQuad(7,3,0,4, false);
    WireQuad(2,6,5,1, false);
    WireQuad(4,5,6,7, false);
    WireQuad(5,4,0,1, false);
}

void glWireCube(void)
{
    // map vertices to faces */
    WireQuad(1,0,3,2, true);
    WireQuad(3,7,6,2, true);
    WireQuad(7,3,0,4, true);
    WireQuad(2,6,5,1, true);
    WireQuad(4,5,6,7, true);
    WireQuad(5,4,0,1, true);
}

void polygon(int a, int b, int c , int d, int face)
{
    /* draw a polygon via list of vertices */
    glBegin(GL_POLYGON);
  	glColor3fv(colors[a]);
  	glVertex3fv(vertices[a]);
  	glColor3fv(colors[b]);
  	glVertex3fv(vertices[b]);
  	glColor3fv(colors[c]);
  	glVertex3fv(vertices[c]);
  	glColor3fv(colors[d]);
  	glVertex3fv(vertices[d]);
    glEnd();
}

void colorcube(void)
{

    /* map vertices to faces */
    polygon(1,0,3,2,0);
    polygon(3,7,6,2,1);
    polygon(7,3,0,4,2);
    polygon(2,6,5,1,3);
    polygon(4,5,6,7,4);
    polygon(5,4,0,1,5);
}

GLvoid swglmDraw(GLMmodel* model)
{ 	         
	GLfloat *n1, *n2, *n3; 
	GLfloat *v1, *v2, *v3;

	//get current color
	GLdouble col[4];
	glGetDoublev(GL_CURRENT_COLOR, col);

    for (unsigned int i = 0; i < model->numtriangles; i++) {
        GLMtriangle* triangle = &(model->triangles[i]);

        n1 = &model->normals[3 * triangle->nindices[0]];
        v1 = &model->vertices[3 * triangle->vindices[0]];
        n2 = &model->normals[3 * triangle->nindices[1]];
        v2 = &model->vertices[3 * triangle->vindices[1]];        
		n3 = &model->normals[3 * triangle->nindices[2]];
        v3 = &model->vertices[3 * triangle->vindices[2]];

		SwglTri(v1[0], v1[1], v1[2],
			    v2[0], v2[1], v2[2],
				v3[0], v3[1], v3[2],
				n1[0], n1[1], n1[2],
			    n2[0], n2[1], n2[2],
				n3[0], n3[1], n3[2],
				col[0], col[1], col[2],
			    col[0], col[1], col[2],
				col[0], col[1], col[2]);   
    }

}

void ReduceToUnit(double vector[3])					// Reduces A Normal Vector (3 Coordinates)
{									// To A Unit Normal Vector With A Length Of One.
	double length;							// Holds Unit Length
	// Calculates The Length Of The Vector
	length = (double)sqrt((vector[0]*vector[0]) + (vector[1]*vector[1]) + (vector[2]*vector[2]));

	if(length == 0.0f)						// Prevents Divide By 0 Error By Providing
		length = 1.0f;						// An Acceptable Value For Vectors To Close To 0.

	vector[0] /= length;						// Dividing Each Element By
	vector[1] /= length;						// The Length Results In A
	vector[2] /= length;						// Unit Normal Vector.
}

void setNormal(double v[3][3], double out[3])				// Calculates Normal For A Quad Using 3 Points
{
	double v1[3],v2[3];						// Vector 1 (x,y,z) & Vector 2 (x,y,z)
	static const int x = 0;						// Define X Coord
	static const int y = 1;						// Define Y Coord
	static const int z = 2;						// Define Z Coord

	// Finds The Vector Between 2 Points By Subtracting
	// The x,y,z Coordinates From One Point To Another.

	// Calculate The Vector From Point 1 To Point 0
	v1[x] = v[0][x] - v[1][x];					// Vector 1.x=Vertex[0].x-Vertex[1].x
	v1[y] = v[0][y] - v[1][y];					// Vector 1.y=Vertex[0].y-Vertex[1].y
	v1[z] = v[0][z] - v[1][z];					// Vector 1.z=Vertex[0].y-Vertex[1].z
	// Calculate The Vector From Point 2 To Point 1
	v2[x] = v[1][x] - v[2][x];					// Vector 2.x=Vertex[0].x-Vertex[1].x
	v2[y] = v[1][y] - v[2][y];					// Vector 2.y=Vertex[0].y-Vertex[1].y
	v2[z] = v[1][z] - v[2][z];					// Vector 2.z=Vertex[0].z-Vertex[1].z
	// Compute The Cross Product To Give Us A Surface Normal
	out[x] = v1[y]*v2[z] - v1[z]*v2[y];				// Cross Product For Y - Z
	out[y] = v1[z]*v2[x] - v1[x]*v2[z];				// Cross Product For X - Z
	out[z] = v1[x]*v2[y] - v1[y]*v2[x];				// Cross Product For X - Y

	ReduceToUnit(out);						// Normalize The Vectors
}

void Relief2Image(vector<GLfloat> src, IplImage *dst)
{
	int width = sqrt( (float) src.size() );
	int height = sqrt( (float) src.size() );

	
	for(int i=0; i < width; i++)
	{
		for(int j=0; j < height; j++)
		{
			cvSetReal2D( dst, j, i, (double) src.at( i*height + j ) );
		}
	}
}

void Image2Relief(IplImage *src, vector<GLfloat> &dst)
{
	int width = cvGetSize(src).width;
	int height = cvGetSize(src).height;

	dst.clear();
	
	for(int i=0; i < width; i++)
	{
		for(int j=0; j < height; j++)
		{
			dst.push_back( (float)cvGetReal2D( src, j, i) );
		}
	}
}

void Image2Relief(IplImage *src, vector<GLdouble> &dst)
{
	int width = cvGetSize(src).width;
	int height = cvGetSize(src).height;

	dst.clear();
	
	for(int i=0; i < width; i++)
	{
		for(int j=0; j < height; j++)
		{
			dst.push_back( cvGetReal2D( src, j, i) );
		}
	}
}

void subSample(vector<GLfloat> array1, vector<GLfloat> &array2)
{
	int width = sqrt( (float) array1.size() );
	int height = sqrt( (float) array1.size() );
	//int  height = winHeight - boundary*2;
	
	GLfloat sum;
	for(int i=0; i < width-1; i+=2)
	{
		for(int j=0; j < height-1; j+=2)
		{
			sum = array1.at( i*height + j ) + array1.at( i*height + j+1 ) + array1.at( (i+1)*height + j ) + array1.at( (i+1)*height + j+1 );
			array2.push_back(sum/4);		
		}
	}
}

void upSample(vector<GLfloat> array1, vector<GLfloat> &array2)
{
	int width = sqrt( (float) array1.size() );
	int height = sqrt( (float) array1.size() );
	/*int width = winWidth/3 - boundary*2;
	int  height = winHeight - boundary*2;*/
	float A,B,C,D;
	int x, y, index;
	float x_ratio = ( (float)(width-1) ) / (width*2 );
	float y_ratio = ( (float)(height-1) ) / (height*2);
	float x_diff, y_diff;
	
	for(int i=0; i < width*2; i++)
	{
		for(int j=0; j < height*2; j++)
		{
			x = (int)(x_ratio * i) ;
            y = (int)(y_ratio * j) ;
			x_diff = (x_ratio * i) - x ;
            y_diff = (y_ratio * j) - y ;
            index = y + x*height ;
			
			A = array1.at( index );
            B = array1.at( index+height );
            C = array1.at( index+1 );
            D = array1.at( index+height+1 );

			array2.push_back( A*(1-x_diff)*(1-y_diff) +  B*(x_diff)*(1-y_diff) + C*(y_diff)*(1-x_diff)   +  D*(x_diff*y_diff) );
		}

		/*for(int j=0; j<height; j++)
		{
			array2.push_back( array1.at( i*height + j ) );	
			array2.push_back( array1.at( i*height + j ) );		
		}*/

		/*for(int j=0; j<height; j++)
		{
			array2.push_back( array1.at( i*height + j ) );	
			array2.push_back( array1.at( i*height + j ) );		
		}*/
	}
}

void reliefSub(vector<GLfloat> src1, vector<GLfloat> src2, vector<GLfloat> &dst)
{
	int width = sqrt( (float) src1.size() );;
	int  height = sqrt( (float) src1.size() );;
	
	for(int i=0; i < width; i++)
	{
		for(int j=0; j < height; j++)
		{
			int adress = i*height + j;
			dst.push_back( src1.at( adress ) - src2.at( adress ) );	
		}
	}
}

void reliefAdd(vector<GLfloat> src1, vector<GLfloat> src2)
{
	int width = sqrt( (float) src1.size() );;
	int  height = sqrt( (float) src1.size() );;
	
	for(int i=0; i < width; i++)
	{
		for(int j=0; j < height; j++)
		{
			int adress = i*height + j;
			src1[adress] = src1.at( adress ) + src2.at( adress );	
		}
	}
}

void reliefAdd(vector<GLfloat> src1, vector<GLfloat> src2, vector<GLfloat> &dst)
{
	int width = sqrt( (float) src1.size() );;
	int  height = sqrt( (float) src1.size() );;
	
	for(int i=0; i < width; i++)
	{
		for(int j=0; j < height; j++)
		{
			int adress = i*height + j;
			dst.push_back( src1.at( adress ) + src2.at( adress ) );	
		}
	}
}

void updateLaplace(vector<GLfloat> src1, vector<GLfloat> src2, vector<GLfloat*> dst)
{
	int width = sqrt( (float) src2.size() );
	int  height = sqrt( (float) src2.size() );
	
	for(int i=0; i < width; i++)
	{
		for(int j=0; j < height; j++)
		{
			int adress = i*height + j;
			if( *dst.at(adress) == 0 )
			{
				*dst.at(adress) =  (float) ( src1.at( adress ) - src2.at( adress ) );
			}
		}
	}
}

void laplacianFilter(vector<GLfloat> src, vector<GLfloat> &dst, int aperture=7)
{
	int width = sqrt( (float) src.size() );
	int  height = sqrt( (float) src.size() );

	dst.clear();
	//vector<GLfloat> *m_src;
	//vector<GLfloat> laplace;
	for(int i=0; i < width; i++)
	{
		for(int j=0; j < height; j++)
		{
			dst.push_back( 0 );
		}
	}

	for(int i=0; i < width; i++)
	{
		for(int j=0; j < height; j++)
		{
			
			//dst.push_back( src1.at( adress ) - src2.at( adress ) );
			vector<GLfloat> sub;
			vector<GLfloat*> subLaplace;

			int ext = (aperture-1) / 2;
			if(i >= ext  && j >= ext && i+ext < width && j+ext < height)
			{
				for(int p=-ext; p <=ext; p++)
				{
					for(int q=-ext; q <=ext; q++)
					{
						//window[ aperture * (p+ext)  +  q+ext ] = cvGetReal2D(src, i+p, j+q);
						sub.push_back(  src.at( (i+p)*height + j+q )  );
						subLaplace.push_back(  &dst.at( (i+p)*height + j+q )  );
					}
				}
			}

			/*else	//boundary issues
			{
				int extU = -ext,extD = ext,extL = -ext,extR = ext;
				if( ext > j)		extU = -j;
				if( ext > i)		extL = -i;
				if( j+ext >= height)		extD = height - 1 - j;
				if( i+ext >= width)		extR =  width - 1 - i;

				//int window[ (extR - extL + 1) * (extD - extU +1) ];
				
				for(int p=extL; p <= extR; p++)
				{
					for(int q=extU; q <= extD; q++)
					{
						//window[ (extR-extL+1) * (p-extU)  +  q-extL ] = cvGetReal2D(src, i+p, j+q);
						sub.push_back(  src.at( (i+p)*height + j+q )  );
						subLaplace.push_back(  &dst.at( (i+p)*height + j+q )  );
					}
				}
			}*/
			//Remapping
			float g = src.at( i*height + j );
			for(int k=0; k < sub.size(); k++)
			{			
				sub[k] = minf( maxf( sub.at(k), g - threshold), g + threshold );
			}

			vector<GLfloat>  sub2, upSub;
			/*IplImage *img0, *img1, *imgUp;
			
			int width = sqrt( (float)sub.size() );
			img0 = cvCreateImage( cvSize(width, width), IPL_DEPTH_32F, 1);	
			Relief2Image(sub, img0);
			
			img1 = cvCreateImage( cvSize( ( width+1)/2, (width+1)/2 ), IPL_DEPTH_32F, 1);
			imgUp = cvCreateImage( cvSize( width,  width), IPL_DEPTH_32F, 1);
			cvPyrDown(img0, img1);*/
			/*for(int i=0; i < img1->width; i++)
			{
				for(int j=0; j < img1->height; j++)
				{
					double value = cvGetReal2D( img1, j, i);
				}
			}*/
			
			/*cvPyrUp(img1, imgUp);
			Image2Relief(imgUp, upSub);*/
			
			
			subSample( sub, sub2 );
			upSample( sub2, upSub );

			
			
			//Image2Relief(imgUp, upSub);
			updateLaplace( sub, upSub, subLaplace);

			/*subSample( sub2, sub3 );
			upSample( sub3, upSub );
			reliefSub( sub2, upSub, subLaplace2);*/
		}
	}


}

double compress(double x, double alpha)
{
	double c;
	c = log10(1+alpha*x) / alpha;
	return c;
}

double compress(double x, double alpha, int n)
{
	double c;
	for(int i=0; i<n; i++)
	{
		c = log(1+alpha*x) / alpha;
		if(c <= 0) break;
		x = c;
	}
	return c;
}

void BuildRelief(vector<GLfloat> &height, GLdouble *pThreadRelief, GLdouble *pThreadNormal)
{
	for(int i=0; i < height.size() / (winHeight - boundary*2) - 1; i++)
		{
			for(int j=0; j < winHeight - boundary*2 -1; j++)
			{
				int position = i*(winHeight-boundary*2) + j;
				
				GLdouble normal[3];
				GLdouble v1[3][3] =  { {i, j, height.at(position)}, {i+1, j, height.at(position+winHeight - boundary*2)}, { i, j+1, height.at(position+1)} };
				
				//glBegin(GL_TRIANGLES);
					setNormal( v1, normal );
					glNormal3dv(normal);
					memcpy( ( pThreadNormal + ( i*(winHeight - boundary*2)*2 + j*2 ) * 3 ), normal, sizeof(GLdouble)*3 );
					/*glVertex3d( i, j, height.at(position) );
					glVertex3d( i+1, j, height.at(position+winHeight - boundary*2) );
					glVertex3d( i, j+1, height.at(position+1) );*/
					memcpy( ( pThreadRelief + ( i*(winHeight - boundary*2) + j ) * 3 ), v1[0], sizeof(GLdouble)*3 );
					memcpy( ( pThreadRelief + ( (i+1)*(winHeight - boundary*2) + j ) * 3 ), v1[1], sizeof(GLdouble)*3 );
					memcpy( ( pThreadRelief + ( i*(winHeight - boundary*2) + (j+1) ) * 3 ), v1[2], sizeof(GLdouble)*3 );
					
				GLdouble v2[3][3] = { {i, j+1, height.at(position+1)}, {i+1, j, height.at(position+winHeight - boundary*2)}, {i+1, j+1, height.at(position+winHeight - boundary*2 + 1)} };
					setNormal( v2, normal);
					glNormal3dv(normal);
					memcpy( ( pThreadNormal + ( i*(winHeight - boundary*2)*2 + j*2 ) * 3 + 3), normal, sizeof(GLdouble)*3 );
					/*glVertex3d( i, j+1, height.at(position+1)  );
					glVertex3d( i+1, j, height.at(position+winHeight - boundary*2) );				
					glVertex3d( i+1, j+1, height.at(position+winHeight - boundary*2 + 1) );*/
					memcpy( ( pThreadRelief + ( i*(winHeight - boundary*2) + j+1 ) * 3 ), v2[0], sizeof(GLdouble)*3 );
					memcpy( ( pThreadRelief + ( (i+1)*(winHeight - boundary*2) + j ) * 3 ), v2[1], sizeof(GLdouble)*3 );
					memcpy( ( pThreadRelief + ( (i+1)*(winHeight - boundary*2) + (j+1) ) * 3 ), v2[2], sizeof(GLdouble)*3 );
				//glEnd();
				//glColor3f(1.0, 0.0, 0.0);
			}
		}
}

//1st Laplace
void softPath(void)
{	
	//Do not change, setting a basic transformation
	glViewport(0, 0, winWidth, winHeight);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	//glOrtho(-2.0, 2.0, -2.0, 2.0, -2.0, 2.0);
	glOrtho(0, winWidth, 0, winHeight, -2.0, 2.0);
    
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//glColor3f(1, 0, 0);
	//OpenglLine(winWidth/2, 0, 0, winWidth, winHeight, 0);


	//
	//replace the opengl function in openglPath() to sotfgl
    //


	//swClearZbuffer();



	//view transform
	glViewport(winWidth/3, 0, winWidth/3, winHeight);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	
	//swOrtho(-2.0, 2.0, -2.0, 2.0, -3.0, 3.0);
	//swFrustum(-2.0, 2.0, -2.0, 2.0, -3.0, 3.0);
	gluPerspective(60, (GLfloat)(winWidth/3)/winHeight, 0.1, 300); 

    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 4, 0, 0, 0, 0, 1, 0);

	glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);
	
	
	/*glPushMatrix();
		glTranslated(-10,15,0);
		glutSolidSphere(1,8,8);
	glPopMatrix();*/

	
	//world coordinate
	glColor3f(1, 0, 0);
	/*OpenglLine(0, 0, 0, 3, 0, 0);
	glColor3f(0, 1, 0);
	OpenglLine(0, 0, 0, 0, 3, 0);
	glColor3f(0, 0, 1);
	OpenglLine(0, 0, 0, 0, 0, 3);*/

	glPushMatrix();
		//glPolygonMode(GL_FRONT,GL_LINE);
		//multiple trackball matrix
		glMultMatrixd(TRACKM);
		
		
		GLfloat maxDepth=0,minDepth=1;
		if(relief)
		{
			disp++;
			for(int i=0; i < pyrLevel; i++)
			{
				heightPyr[i].clear();
			}
			
			for(int i=boundary; i<winWidth/3 - boundary; i++)
			{
				for(int j=boundary; j<winHeight - boundary; j++)
				{
					GLfloat depth;
					glReadPixels(i, j, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

					if( maxDepth < depth && depth !=1)
					{
						maxDepth = depth;
					}
					if( minDepth > depth)
					{
						minDepth = depth;
					}
					heightPyr[0].push_back(depth);
					
				}
			}

			

			int enhance = 50;
			double maxDepthPrime = compress(maxDepth*enhance, 5, 5);
			double minDepthPrime = compress(minDepth*enhance, 5, 5);

			for(int i=0; i<heightPyr[0].size(); i++)
			{			
				if ( heightPyr[0].at(i) == 1 )	//infinite point
				{
					heightPyr[0].at(i) = 0;
				}
				else
				{
					/*height.at(i) = compress(height.at(i)*enhance, 5, 5);
					height.at(i) = maxDepthPrime - height.at(i);
					//height.at(i) = maxDepth - height.at(i);	//transfer depth to height
					height.at(i) /= (maxDepthPrime - minDepthPrime );	//normalize to [0,1]*/

					heightPyr[0].at(i) = maxDepth - heightPyr[0].at(i);
					//height.at(i) = maxDepth - height.at(i);	//transfer depth to height
					heightPyr[0].at(i) /= (maxDepth - minDepth);	//normalize to [0,1]
				}
			}
			heightList.clear();
			vector<GLfloat> height2, upHeight, laplace;
			heightList.push_back(heightPyr[0]);
			IplImage *img0, *img1, *imgUP, *imgLa;

			int width = sqrt( (float) heightPyr[0].size() );
			//int height = sqrt( (float) height[0].size() );
			img0 = cvCreateImage( cvSize(width, width), IPL_DEPTH_32F, 1);	
			Relief2Image(heightPyr[0], img0);
			imgPyr = cvCreatePyramid(img0, pyrLevel, 0.5);
			for(int i=0; i < pyrLevel-1; i++)
			{
				subSample( heightPyr[i], heightPyr[i+1] );
			}
			//Relief2Image(height[1], img1);
			imgUP = cvCreateImage( cvSize( imgPyr[1]->width*2,  imgPyr[1]->height*2), IPL_DEPTH_32F, 1);
			imgLa = cvCreateImage( cvSize( imgPyr[1]->width*2,  imgPyr[1]->height*2), IPL_DEPTH_32F, 1);
			/*cvPyrUp(img[1], imgUP);
			cvSub(img0, imgUP, imgLa);

			Image2Relief(imgLa, laplace);*/
			
			//subSample(height, height2);
			//heightList.push_back(height2);
			//upSample(height2, upHeight);
			//reliefSub(height, upHeight, laplace);

			laplacianFilter(heightPyr[0], laplace);
			laplaceList.clear();
			laplaceList.push_back(laplace);
			//for(int i=0; i<height.size(); i++)
			//{
			//	if ( height.at(i) == 1 )	//infinite point
			//	{
			//		height.at(i) = 0;
			//	}
			//	else
			//	{
			//		height.at(i) = ( height.at(i) - minDepth ) / (maxDepth - minDepth) ;	//normalize to [0,1]
			//		height.at(i) = 1 - height.at(i);	//transfer depth to height
			//	}
			//}
			BuildRelief(laplace, pThreadRelief, pThreadNormal);
			relief = false;
		}
		//swScaled(MODELSCALE, MODELSCALE, MODELSCALE);
		glColor3f(0.6, 0.6, 0.6);
		
		glRotatef(reliefAngleX, 1, 0, 0);
		glRotatef(reliefAngleY, 0, 1, 0);
		glRotatef(reliefAngleZ, 0, 0, 1);

		glTranslated(-1.8, -1.9, 0.4);

		glScalef(0.004*scale,0.004*scale,outputHeight);
		
			for(int i=0; i < heightPyr[0].size() / (winHeight - boundary*2) - 1; i++)
			{			
				for(int j=0; j < winHeight - boundary*2 - 1; j++)
				{
					glBegin(GL_TRIANGLES);
					glNormal3dv(pThreadNormal + ( ( i*(winHeight - boundary*2)  +  j)*2 ) *3);
					glVertex3dv( pThreadRelief + ( i*(winHeight - boundary*2) + j ) * 3 );
					glVertex3dv( pThreadRelief + ( (i+1)*(winHeight - boundary*2) + j ) * 3 );
					glVertex3dv( pThreadRelief + ( i*(winHeight - boundary*2) + (j+1) ) * 3 );
					glEnd();
					
					glBegin(GL_TRIANGLES);
					glNormal3dv(pThreadNormal + ( ( i*(winHeight - boundary*2)  +  j)*2 )  *3 + 3);
					glVertex3dv( pThreadRelief + ( i*(winHeight - boundary*2) + j+1 ) * 3 );
					glVertex3dv( pThreadRelief + ( (i+1)*(winHeight - boundary*2) + j ) * 3 );
					glVertex3dv( pThreadRelief + ( (i+1)*(winHeight - boundary*2) + (j+1) ) * 3);
					glEnd();
				}
				
			}
		
		//for(int i=0; i < height.size() / (winHeight - boundary*2) - 1; i++)
		//{
		//	for(int j=0; j < winHeight - boundary*2 -1; j++)
		//	{
		//		int position = i*(winHeight-boundary*2) + j;
		//		if( height.at(position) >= 0 )
		//		{
		//			i++;
		//			i--;
		//		}
		//		if( height.at(position) >= 0.5 )
		//		{
		//			glColor3f(1.0, 1.0, 1.0);
		//		}
		//		
		//		GLdouble normal[3];
		//		GLdouble v1[3][3] =  { {i, j, height.at(position)}, {i+1, j, height.at(position+winHeight - boundary*2)}, { i, j+1, height.at(position+1)} };
		//		
		//		glBegin(GL_TRIANGLES);
		//			setNormal( v1, normal );
		//			glNormal3dv(normal);
		//			glVertex3d( i, j, height.at(position) );
		//			glVertex3d( i+1, j, height.at(position+winHeight - boundary*2) );
		//			glVertex3d( i, j+1, height.at(position+1) );
		//			
		//		GLdouble v2[3][3] = { {i, j+1, height.at(position+1)}, {i+1, j, height.at(position+winHeight - boundary*2)}, {i+1, j+1, height.at(position+winHeight - boundary*2 + 1)} };
		//			setNormal( v2, normal);
		//			glNormal3dv(normal);
		//			glVertex3d( i, j+1, height.at(position+1)  );
		//			glVertex3d( i+1, j, height.at(position+winHeight - boundary*2) );				
		//			glVertex3d( i+1, j+1, height.at(position+winHeight - boundary*2 + 1) );
		//		glEnd();
		//		//glColor3f(1.0, 0.0, 0.0);
		//	}
		//}
		//swglmDraw(MODEL);
		
	glPopMatrix();

	glPushMatrix();
		glTranslated(0, 2, 0);
		glMultMatrixd(TRACKM);

		
	glPopMatrix();
}

void equalizeHist(vector<GLfloat> src, vector<GLfloat> &dst)
{	
	int hist[1000+1];
	for(int i=0; i<= 1000; i++)
	{
		hist[i] =0;
	}
	int srcHeight = sqrt( (float)src.size() );
	for(int i=0; i< srcHeight; i++)
	{
		for(int j=0; j< srcHeight; j++)
		{
			
			hist [ (int) (src.at( i*srcHeight + j)  *1000) ] +=1;		
		}
	}
	
	double sum[1000+1];
	
	int number =  src.size() - hist[0];
	for(int i=1; i < 1000+1; i++)
	{
		if( i == 1 )
		{
			sum[i] =  hist[i] * (1000 - 1) /  number;
			//cout << "Sum " << sum[i] << std::endl;
			//sum2[i] =  cvGetReal1D(lHist1->bins, i);
		}
		else
		{
			sum[i] = sum[i-1] + (double)( hist[i] * (1000 - 1) ) / number;
			//sum2[i] = sum2[i-1] + cvGetReal1D(lHist1->bins, i);
		}
	}
	dst.clear();
	for(int i=0; i< srcHeight; i++)
	{
		for(int j=0; j< srcHeight; j++)
		{		
			if(src.at( i*srcHeight + j) == 0)
			{
				dst.push_back(0);
			}
			else
			{
				dst.push_back( sum [ (int) (src.at( i*srcHeight + j)  *1000) ]  );
			}
		}
	}
		
}

//Bas-Relief
void reliefHistogram(void)
{	
	//Do not change, setting a basic transformation
	glViewport(0, 0, winWidth, winHeight);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	//glOrtho(-2.0, 2.0, -2.0, 2.0, -2.0, 2.0);
	glOrtho(0, winWidth, 0, winHeight, -2.0, 2.0);
    
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//glColor3f(1, 0, 0);
	//OpenglLine(winWidth/2, 0, 0, winWidth, winHeight, 0);


	//
	//replace the opengl function in openglPath() to sotfgl
    //


	//swClearZbuffer();



	//view transform
	glViewport(winWidth*2/3, 0, winWidth*2/3, winHeight);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	
	//swOrtho(-2.0, 2.0, -2.0, 2.0, -3.0, 3.0);
	//swFrustum(-2.0, 2.0, -2.0, 2.0, -3.0, 3.0);
	gluPerspective(60, (GLfloat)(winWidth*2/3)/winHeight, 0.1, 300); 

    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, 4, 0, 0, 0, 0, 1, 0);

	glLightfv(GL_LIGHT2, GL_POSITION, lightPos2);
	
	/*glPushMatrix();
		glTranslated(-10,15,0);
		glutSolidSphere(1,8,8);
	glPopMatrix();*/

	
	//world coordinate
	glColor3f(1, 0, 0);
	/*OpenglLine(0, 0, 0, 3, 0, 0);
	glColor3f(0, 1, 0);
	OpenglLine(0, 0, 0, 0, 3, 0);
	glColor3f(0, 0, 1);
	OpenglLine(0, 0, 0, 0, 0, 3);*/

	glPushMatrix();
		//glPolygonMode(GL_FRONT,GL_LINE);
		//multiple trackball matrix
		glMultMatrixd(TRACKM);
		
		
		//GLfloat maxDepth=0,minDepth=1;
		if(relief2)
		{
		//	disp++;
		//	height.clear();
		//	
		//	for(int i=boundary; i<winWidth*0.5 - boundary; i++)
		//	{
		//		for(int j=boundary; j<winHeight - boundary; j++)
		//		{
		//			GLfloat depth;
		//			glReadPixels(i, j, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

		//			if( maxDepth < depth && depth !=1)
		//			{
		//				maxDepth = depth;
		//			}
		//			if( minDepth > depth)
		//			{
		//				minDepth = depth;
		//			}
		//			height.push_back(depth);
		//			
		//		}
		//	}
		//	int enhance = 50;
		//	double maxDepthPrime = compress(maxDepth*enhance, 5, 5);
		//	double minDepthPrime = compress(minDepth*enhance, 5, 5);

		//	for(int i=0; i<height.size(); i++)
		//	{			
		//		if ( height.at(i) == 1 )	//infinite point
		//		{
		//			height.at(i) = 0;
		//		}
		//		else
		//		{
		//			height.at(i) = compress(height.at(i)*enhance, 5, 5);
		//			height.at(i) = maxDepthPrime - height.at(i);
		//			//height.at(i) = maxDepth - height.at(i);	//transfer depth to height
		//			height.at(i) /= (maxDepthPrime - minDepthPrime );	//normalize to [0,1]
		//		}
		//	}
			//for(int i=0; i<height.size(); i++)
			//{
			//	if ( height.at(i) == 1 )	//infinite point
			//	{
			//		height.at(i) = 0;
			//	}
			//	else
			//	{
			//		height.at(i) = ( height.at(i) - minDepth ) / (maxDepth - minDepth) ;	//normalize to [0,1]
			//		height.at(i) = 1 - height.at(i);	//transfer depth to height
			//	}
			//}
			
			//equalizeHist(height[0], equalizeHeight);

			/*float maxHeight=0, minHeight=1000;
			for(int i=0; i < equalizeHeight.at(0).size(); i++)
			{
				float height = equalizeHeight.at(i);
				if( maxHeight < height )
				{
					maxHeight = height;
				}
				if( minHeight > height && height !=0 )
				{
					minHeight = height;
				}
			}
			float range = maxHeight - minHeight;
			for(int i=0; i < equalizeHeight.size(); i++)
			{
				equalizeHeight.at(i) /= 1000;
			}*/

			//BuildRelief(equalizeHeight, pThreadEqualizeRelief, pThreadEqualizeNormal);
			
			//heightList.push_back(height);
			vector<GLfloat> height2, upHeight, laplace[pyrLevel - 1];
			IplImage *img, *img2, *imgLa;
			/*subSample(heightList.at(0), height2);
			heightList.push_back(height2);*/
			//upSample(height2, upHeight);
			//reliefSub(heightList.at(0), upHeight, laplace);

			/*laplacianFilter(height[1], *laplace);
			laplaceList.push_back(*laplace);
			laplacianFilter(height[2], *laplace);
			laplaceList.push_back(*laplace);*/
			IplImage *tempImg;
			for(int i=1;i< pyrLevel - 1;i++)
			{
				/*tempImg = cvCreateImage( cvSize(imgPyr[i]->width, imgPyr[i]->height), IPL_DEPTH_32F, 1);
				cvGetImage(imgPyr[i], tempImg);
				Image2Relief(tempImg, height[i]);
				laplacianFilter(height[i], laplace);
				laplaceList.push_back(laplace);*/
				laplacianFilter(heightPyr[i], laplace[i]);
				laplaceList.push_back(laplace[i]);
			}
			
			//collapse the pyramid
			vector<GLfloat> compressedH;
			/*for(int i=0;i < height[pyrLevel - 1].size(); i++)
			{
				compressedH->push_back( height[pyrLevel - 1].at(i) );
			}*/
			int width = sqrt( (float) heightPyr[pyrLevel - 1].size() );
			img = cvCreateImage( cvSize( width*2,  width*2), IPL_DEPTH_32F, 1);
			imgLa = cvCreateImage( cvSize( width*2,  width*2), IPL_DEPTH_32F, 1);

			tempImg = cvCreateImage( cvSize(width, width), IPL_DEPTH_32F, 1);
			//cvGetImage(imgPyr[pyrLevel - 1], tempImg);
			Relief2Image(	heightPyr.at( pyrLevel - 1), tempImg);
			
			cvPyrUp( tempImg, img );

			tempImg = cvCreateImage( cvSize( width*2,  width*2), IPL_DEPTH_32F, 1);
			cvCopy(img, tempImg);
			Relief2Image(laplaceList.at( pyrLevel - 2), imgLa);
			cvAdd(tempImg, imgLa, img);
			//Image2Relief(img, compressedH);
			
			for(int i=pyrLevel - 2; i > 0; i--)
			{
				int width = imgPyr[i]->width;
				//img = cvCreateImage( cvSize( width,  width), IPL_DEPTH_32F, 1);
				img2 = cvCreateImage( cvSize( width*2,  width*2), IPL_DEPTH_32F, 1);
				imgLa = cvCreateImage( cvSize( width*2,  width*2), IPL_DEPTH_32F, 1);

				cvPyrUp(img, img2);
				Relief2Image(laplaceList.at(i-1), imgLa);

				img = cvCreateImage( cvSize( width*2,  width*2), IPL_DEPTH_32F, 1);
				cvAdd(img2, imgLa, img);
				/*Image2Relief(img2, compressedH);
				reliefAdd( *compressedH, laplaceList.at(i) );*/
			}
			Image2Relief(img, compressedH);
			//upHeight.clear();
		/*	height2.clear();
			upSample(heightList.at(1), height2);

			reliefAdd(height2, laplace, upHeight);
			height2.clear();
			upSample(upHeight, height2);*/
			//upHeight.clear();
			//reliefAdd(height2, laplaceList.at(0), upHeight);

			BuildRelief(compressedH, pThreadEqualizeRelief, pThreadEqualizeNormal);

			
			relief2 = false;
		}
		//swScaled(MODELSCALE, MODELSCALE, MODELSCALE);
		glColor3f(0.6, 0.6, 0.6);
		
		glRotatef(reliefAngleX, 1, 0, 0);
		glRotatef(reliefAngleY, 0, 1, 0);
		glRotatef(reliefAngleZ, 0, 0, 1);

		glTranslated(-3.8, -1.9, 0.4);

		glScalef(0.004*scale,0.004*scale,outputHeight);
		
			for(int i=0; i <heightPyr.at(0).size() / (winHeight - boundary*2) - 1; i++)
			{			
				for(int j=0; j < winHeight - boundary*2 - 1; j++)
				{
					glBegin(GL_TRIANGLES);
					glNormal3dv(pThreadEqualizeNormal + ( ( i*(winHeight - boundary*2)  +  j)*2 ) *3);
					glVertex3dv( pThreadEqualizeRelief + ( i*(winHeight - boundary*2) + j ) * 3 );
					glVertex3dv( pThreadEqualizeRelief + ( (i+1)*(winHeight - boundary*2) + j ) * 3 );
					glVertex3dv( pThreadEqualizeRelief + ( i*(winHeight - boundary*2) + (j+1) ) * 3 );
					glEnd();
					
					glBegin(GL_TRIANGLES);
					glNormal3dv(pThreadEqualizeNormal + ( ( i*(winHeight - boundary*2)  +  j)*2 )  *3 + 3);
					glVertex3dv( pThreadEqualizeRelief + ( i*(winHeight - boundary*2) + j+1 ) * 3 );
					glVertex3dv( pThreadEqualizeRelief + ( (i+1)*(winHeight - boundary*2) + j ) * 3 );
					glVertex3dv( pThreadEqualizeRelief + ( (i+1)*(winHeight - boundary*2) + (j+1) ) * 3);
					glEnd();
				}
				
			}
		
		//for(int i=0; i < height.size() / (winHeight - boundary*2) - 1; i++)
		//{
		//	for(int j=0; j < winHeight - boundary*2 -1; j++)
		//	{
		//		int position = i*(winHeight-boundary*2) + j;
		//		if( height.at(position) >= 0 )
		//		{
		//			i++;
		//			i--;
		//		}
		//		if( height.at(position) >= 0.5 )
		//		{
		//			glColor3f(1.0, 1.0, 1.0);
		//		}
		//		
		//		GLdouble normal[3];
		//		GLdouble v1[3][3] =  { {i, j, height.at(position)}, {i+1, j, height.at(position+winHeight - boundary*2)}, { i, j+1, height.at(position+1)} };
		//		
		//		glBegin(GL_TRIANGLES);
		//			setNormal( v1, normal );
		//			glNormal3dv(normal);
		//			glVertex3d( i, j, height.at(position) );
		//			glVertex3d( i+1, j, height.at(position+winHeight - boundary*2) );
		//			glVertex3d( i, j+1, height.at(position+1) );
		//			
		//		GLdouble v2[3][3] = { {i, j+1, height.at(position+1)}, {i+1, j, height.at(position+winHeight - boundary*2)}, {i+1, j+1, height.at(position+winHeight - boundary*2 + 1)} };
		//			setNormal( v2, normal);
		//			glNormal3dv(normal);
		//			glVertex3d( i, j+1, height.at(position+1)  );
		//			glVertex3d( i+1, j, height.at(position+winHeight - boundary*2) );				
		//			glVertex3d( i+1, j+1, height.at(position+winHeight - boundary*2 + 1) );
		//		glEnd();
		//		//glColor3f(1.0, 0.0, 0.0);
		//	}
		//}
		//swglmDraw(MODEL);
		
	glPopMatrix();

	glPushMatrix();
		glTranslated(0, 2, 0);
		glMultMatrixd(TRACKM);

		
	glPopMatrix();
}

//3D Scene
void openglPath(void)
{
    //view transform
	glViewport(0, 0, winWidth/3, winHeight);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	//glOrtho(-2.0, 2.0, -2.0, 2.0, -3.0, 25.0);
	//glFrustum(-2.0, 2.0, -2.0, 2.0, -3.0, 3.0);
	gluPerspective(60, (GLfloat)(winWidth*0.5)/winHeight, 0.1, 25); 
	glGetDoublev(GL_PROJECTION_MATRIX, DEBUG_M);


    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(0, 0, 4, 0, 0, 0, 0, 1, 0);
	glGetDoublev(GL_MODELVIEW_MATRIX, DEBUG_M);

	lightPos0[0] = -LIGHTP;
	lightPos0[1] = LIGHTP;
	lightPos0[2] = LIGHTP;
	
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);

	

	glTranslated(0, 0, 0);
	//world coordinate
	glColor3f(1, 0, 0);
	
	/*glPushMatrix();
		glTranslated(0,0,0);
		glutSolidSphere(1.5,9,9);
	glPopMatrix();*/
	/*OpenglLine(0, 0, 0, 3, 0, 0);
	glColor3f(0, 1, 0);
	OpenglLine(0, 0, 0, 0, 3, 0);
	glColor3f(0, 0, 1);
	OpenglLine(0, 0, 0, 0, 0, 3);*/

	
	glPushMatrix();
		//multiple trackball matrix
		glRotated(angleX,1, 0, 0);
		glRotated(angleY,0, 1, 0);
		glRotated(angleZ,0, 0, 1);

		glMultMatrixd(TRACKM);

		glScaled(MODELSCALE, MODELSCALE, MODELSCALE);
		glColor3f(1.0, 1.0, 1.0);
		glmDraw(MODEL, GLM_SMOOTH);//GLM_FLAT
		//glutSolidSphere(1, 20, 20);
	glPopMatrix();

	glPushMatrix();
		glTranslated(0, 2, 0);
		glMultMatrixd(TRACKM);

		/*glBegin(GL_TRIANGLES);
			glNormal3f(0, 0, 1);
			glColor3f(1, 0, 0);
			glVertex3f(-1, 0, 0);

			glColor3f(0, 1, 0);
			glVertex3f(1, 0, 0);

			glColor3f(0, 0, 1);
			glVertex3f(0, 1, 0);
		glEnd();		*/
	glPopMatrix();

}

/*----------------------------------------------------------------------*/
/* 
** These functions implement a simple trackball-like motion control.
*/

float lastPos[3] = {0.0F, 0.0F, 0.0F};
int curx, cury;
int startX, startY;

void trackball_ptov(int x, int y, int width, int height, float v[3])
{
    float d, a;

    /* project x,y onto a hemi-sphere centered within width, height */
    v[0] = (2.0F*x - width) / width;
    v[1] = (height - 2.0F*y) / height;
    d = (float) sqrt(v[0]*v[0] + v[1]*v[1]);
    v[2] = (float) cos((M_PI/2.0F) * ((d < 1.0F) ? d : 1.0F));
    a = 1.0F / (float) sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] *= a;
    v[1] *= a;
    v[2] *= a;
}


void mouseMotion(int x, int y)
{
    float curPos[3], dx, dy, dz;

    trackball_ptov(x, y, winWidth, winHeight, curPos);
	if(trackingMouse)
	{
		dx = curPos[0] - lastPos[0];
		dy = curPos[1] - lastPos[1];
		dz = curPos[2] - lastPos[2];

		if (dx || dy || dz) {
			angle = 90.0F * sqrt(dx*dx + dy*dy + dz*dz);

			axis[0] = lastPos[1]*curPos[2] - lastPos[2]*curPos[1];
			axis[1] = lastPos[2]*curPos[0] - lastPos[0]*curPos[2];
			axis[2] = lastPos[0]*curPos[1] - lastPos[1]*curPos[0];

			lastPos[0] = curPos[0];
			lastPos[1] = curPos[1];
			lastPos[2] = curPos[2];
		}
	} 
    glutPostRedisplay();
}

void startMotion(int x, int y)
{
    trackingMouse = true;
    redrawContinue = false;
    startX = x; startY = y;
    curx = x; cury = y;
    trackball_ptov(x, y, winWidth, winHeight, lastPos);
	trackballMove=true;
}

void stopMotion(int x, int y)
{
	trackingMouse = false;

    if (startX != x || startY != y) {
		redrawContinue = true;
    } else {
		angle = 0.0F;
		redrawContinue = false;
		trackballMove = false;
    }
}

/*----------------------------------------------------------------------*/

void displayfont()
{
    //Font
	char mss[30]="3D Scene";
	//sprintf(mss, "Score %d", Gamescore);
	glColor3f(1.0, 0.0, 0.0);  //set font color
	void * font = GLUT_BITMAP_9_BY_15;

	glWindowPos2i(10, winHeight-20);    //set font start position
	for(unsigned int i=0; i<strlen(mss); i++) {
		glutBitmapCharacter(font, mss[i]);
	}

	char mss1[30]="1st Laplace";
	
	glWindowPos2i(10+(winWidth/3), winHeight-20);    //set font start position
	for(unsigned int i=0; i<strlen(mss1); i++) {
		glutBitmapCharacter(font, mss1[i]);
	}
	
	char mss2[30]="Computing...";
	if(relief)
	{
		glColor3f(0.0, 1.0, 0.0);
		glWindowPos2i(10+(winWidth/3), winHeight-40);    //set font start position
		for(unsigned int i=0; i<strlen(mss2); i++) {
			glutBitmapCharacter(font, mss2[i]);
		}
	}
	
	char mss3[30]="Bas-Relief";
	
	glWindowPos2i(10+(winWidth*2/3), winHeight-20);    //set font start position
	for(unsigned int i=0; i<strlen(mss1); i++) {
		glutBitmapCharacter(font, mss3[i]);
	}
	//char mss3[30];
	//sprintf(mss3,"%f",scale);
	//glWindowPos2i(10+(winWidth/2), winHeight-60);    //set font start position
	//for(unsigned int i=0; i<strlen(mss1); i++) {
	//	glutBitmapCharacter(font, mss3[i]);
	//}
}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	displayfont();

    if (trackballMove) {
		glPushMatrix();
			glLoadMatrixd(TRACKM);
			glRotatef(angle, axis[0], axis[1], axis[2]);
			glGetDoublev(GL_MODELVIEW_MATRIX, TRACKM);
		glPopMatrix();	    
	}
	
	glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	openglPath();

	//we must disable the opengl's depth test, then the software depth test will work 
	//glDisable(GL_DEPTH_TEST); 
	//glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	softPath();
	//glEnable(GL_LIGHTING);
	//glEnable(GL_DEPTH_TEST); 
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
	reliefHistogram();

    glutSwapBuffers();
}

/*----------------------------------------------------------------------*/

void mouseButton(int button, int state, int x, int y)
{
	if(button==GLUT_RIGHT_BUTTON) {
		exit(0);
	}

	if(button==GLUT_LEFT_BUTTON) switch(state) 
	{
		case GLUT_DOWN:
			y=winHeight-y;
			startMotion(x, y);
			break;
		case GLUT_UP:
			y=winHeight-y;
			stopMotion(x, y);
			break;
    } 
}

void myReshape(int w, int h)
{
    winWidth = w;
    winHeight = h;

	swInitZbuffer(w/2, h);
}

void spinCube()
{
    if (redrawContinue) glutPostRedisplay();
}

void update(int i)
{
	TICK++;
	int temp=TICK%180;
	if(temp<90)
		Angle1++;
	else
		Angle1--;

	//int temp2=TICK%90;
	if(temp<90)
		Angle2+=0.5;
	else
		Angle2-=0.5;

	glutPostRedisplay();
	glutTimerFunc(33, update, ++i);
}

void myKeys(unsigned char key, int x, int y)
{
	switch(key)
	{
		case ' ':
			glPushMatrix();
				glLoadIdentity();
				glGetDoublev(GL_MODELVIEW_MATRIX, TRACKM);
			glPopMatrix();	 
			break;

		case 'a':
			MODELSCALE += 0.5;
			break;
		case 's':
			if(MODELSCALE > 0.1)
				MODELSCALE -= 0.5;
			break;

		case 'z':
			LIGHTP += 1;
			std::cout<<"LIGHTP "<<LIGHTP<<'\n';
			break;
		case 'x':
			LIGHTP -= 1;
			std::cout<<"LIGHTP "<<LIGHTP<<'\n';
			break;

        case 'Q':
        case 'q':  
			exit(0); 
			break;

        case '0':  
			DRAWTYPE=0; 
			break;
        case '1':  
			DRAWTYPE=1; 
			break;
        case '4':  
			reliefAngleY--;
			break;
        case '6':  
			reliefAngleY++;
			break;
		case '8':  
			 reliefAngleX--;
			break;
		case '2':  
			reliefAngleX++;
			break;
		case '+':  
			scale += 0.5; 
			break;
        case '-': 
			if(scale > 0.1)
			scale -= 0.5; 
			break;
		//Enter
		case 13:  
			relief = true;
			relief2 = true; 
			break;
		// Delete
		case 127 :
			reliefAngleZ++;
			break;
	}
	glutPostRedisplay();
}

void SpecialKeys(int key, int x, int y)
{
	switch(key)
	{
		case GLUT_KEY_LEFT:
			angleY--;
			break;
		case GLUT_KEY_UP:
			angleX--;
			break;
		case GLUT_KEY_PAGE_UP:
			angleZ--;
			break;
		case GLUT_KEY_RIGHT:
			angleY++;
			break;
		case GLUT_KEY_DOWN:
			angleX++;
			break;
		case GLUT_KEY_PAGE_DOWN:
			angleZ++;
			break;
		case GLUT_KEY_INSERT:
			reliefAngleZ--;
			break;
	}
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1200, 400);
    glutCreateWindow("Digital Bas-Relief from 3D Scenes");

    glutReshapeFunc(myReshape);

	vertCount = ( 400 - boundary*2 ) * ( 400 - boundary*2 );
	pThreadRelief = (GLdouble*) malloc ( sizeof(GLdouble) * vertCount *3);
	pThreadNormal = (GLdouble*) malloc ( sizeof(GLdouble) * vertCount *3 *2);
	pThreadEqualizeRelief = (GLdouble*) malloc ( sizeof(GLdouble) * vertCount *3);
	pThreadEqualizeNormal = (GLdouble*) malloc ( sizeof(GLdouble) * vertCount *3 *2);

    glutDisplayFunc(display);
    glutIdleFunc(spinCube);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
	glutKeyboardFunc(myKeys);
	glutSpecialFunc(SpecialKeys);
	glutTimerFunc(33, update, 0);

	glEnable(GL_DEPTH_TEST); 
	//glDepthRange(0.0, 5.0);
	//Font = FontCreate(wglGetCurrentDC(), "Times", 32, 0, 1);



	//Read model
	MODEL = glmReadOBJ("dragon.obj");
	glmUnitize(MODEL);
	//glmFacetNormals(MODEL);
	//glmVertexNormals(MODEL, 90);

    // Light values and coordinates
    GLfloat  ambientLight0[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    GLfloat  diffuseLight0[] = { 0.7f, 0.7f, 0.7f, 1.0f };
    GLfloat  specular0[] = { 0.7f, 0.7f, 0.7f, 1.0f };

	 GLfloat  ambientLight1[] = { 0.1f, 0.1f, 0.1f, 1.0f };
    GLfloat  diffuseLight1[] = { 0.3f, 0.3f, 0.3f, 1.0f };
    GLfloat  specular1[] = { 0.0f, 0.0f, 0.0f, 1.0f };


    GLfloat  specref[] = { 0.25f, 0.25f, 0.25f, 0.25f };
	GLfloat  shininess = 32.0f;

    // Enable lighting
    glEnable(GL_LIGHTING);

    // Setup and enable light 0
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientLight0);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuseLight0);
    glLightfv(GL_LIGHT0,GL_SPECULAR, specular0);
    glEnable(GL_LIGHT0);

	glLightfv(GL_LIGHT1,GL_AMBIENT,ambientLight1);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,diffuseLight1);
    glLightfv(GL_LIGHT1,GL_SPECULAR, specular1);

	glLightfv(GL_LIGHT2,GL_AMBIENT,ambientLight1);
    glLightfv(GL_LIGHT2,GL_DIFFUSE,diffuseLight1);
    glLightfv(GL_LIGHT2,GL_SPECULAR, specular1);
    

    // Enable color tracking
    glEnable(GL_COLOR_MATERIAL);
    // Set Material properties to follow glColor values
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	// All materials hereafter have full specular reflectivity with a high shine 
    glMaterialfv(GL_FRONT, GL_SPECULAR, specref);
    glMateriali(GL_FRONT, GL_SHININESS, shininess);

	//
	//hw3
	//
	swLightfv(GL_LIGHT0,GL_AMBIENT,ambientLight0);
    swLightfv(GL_LIGHT0,GL_DIFFUSE,diffuseLight0);
    swLightfv(GL_LIGHT0,GL_SPECULAR, specular0);

    swMaterialfv(GL_FRONT, GL_SPECULAR, specref);
    swMateriali(GL_FRONT, GL_SHININESS, shininess);


    glutMainLoop();

	free(pThreadRelief);
	free(pThreadNormal);
	free(pThreadEqualizeRelief);
	free(pThreadEqualizeNormal);

	cvWaitKey(0);
}