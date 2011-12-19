#ifndef __swgl_h__
#define __swgl_h__
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <string>
#include <list>

#include <stack>
#include "GLee.h"
#include "math3d.h"


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

//implement the following function call to archive the opengl pipeline
void swTranslated(GLdouble x, GLdouble y, GLdouble z);
void swScaled(GLdouble x, GLdouble y, GLdouble z);
void swRotated(GLdouble angle, GLdouble x, GLdouble y, GLdouble z);

void swMatrixMode(GLenum mode);
void swLoadIdentity(void);
void swLoadMatrixd(const GLdouble * m);
void swMultMatrixf(const GLdouble * m);
void swPushMatrix(void);
void swPopMatrix(void);

void swuLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ,
 	           GLdouble centerX, GLdouble centerY, GLdouble centerZ,
 	           GLdouble upX, GLdouble upY, GLdouble upZ);
void swFrustum(	GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble nearVal, GLdouble farVal);
void swuPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar);

void swViewport(GLint x, GLint y, GLsizei width, GLsizei height);

bool swTransformation(const GLdouble h[4], GLdouble w[4]);

//---------------------------------------------------------------------------
//cghw2
//---------------------------------------------------------------------------
void writepixel(int x, int y, GLdouble r, GLdouble g, GLdouble b);

//Bresenham's algorithm
bool BresenhamLine(int x1, int y1, int x2, int y2, GLdouble r, GLdouble g, GLdouble b);
bool BresenhamLine(int x1, int y1, GLdouble r1, GLdouble g1, GLdouble b1, int x2, int y2, GLdouble r2, GLdouble g2, GLdouble b2);

bool swTriangle(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble r, GLdouble g, GLdouble b);

bool swInitZbuffer(int width, int height);
bool swClearZbuffer();


//---------------------------------------------------------------------------
//cghw3
//---------------------------------------------------------------------------

void writepixelfast(int x, int y, GLdouble r, GLdouble g, GLdouble b);


bool swNormalTransformation(const GLdouble h[4], GLdouble w[4]);

void swLightfv(GLenum light, GLenum pname, const GLfloat *params);
void swMaterialfv (GLenum face, GLenum pname, const GLfloat *params);
void swMateriali (GLenum face, GLenum pname, GLint param);

//Gouraud shading: phong shading model on each vertex, then interpolate on each fragment
//vertex position:	(x1, y1, z1)   in object space
//vertex normal:	(nx1, ny1, nz1)
//vertex color:		(r1, g1, b1)
bool swTriangleG(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble nx1, GLdouble ny1, GLdouble nz1, 
			 GLdouble nx2, GLdouble ny2, GLdouble nz2, 
			 GLdouble nx3, GLdouble ny3, GLdouble nz3,
			 GLdouble r1, GLdouble g1, GLdouble b1,
			 GLdouble r2, GLdouble g2, GLdouble b2,
			 GLdouble r3, GLdouble g3, GLdouble b3);


////Phong Shading: phong shading model on each fragment
//vertex position:	(x1, y1, z1)   in object space
//vertex normal:	(nx1, ny1, nz1)
//vertex color:		(r1, g1, b1)
bool swTriangleP(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble nx1, GLdouble ny1, GLdouble nz1, 
			 GLdouble nx2, GLdouble ny2, GLdouble nz2, 
			 GLdouble nx3, GLdouble ny3, GLdouble nz3,
			 GLdouble r1, GLdouble g1, GLdouble b1,
			 GLdouble r2, GLdouble g2, GLdouble b2,
			 GLdouble r3, GLdouble g3, GLdouble b3);



#endif                  /* __swgl_h__ */
