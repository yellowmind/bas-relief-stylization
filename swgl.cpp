#include <vector>
#include "swgl.h"

#define PI 3.1415
using namespace std;
stack< vector<GLdouble> > matrixStack;
vector< vector<GLdouble> > zbuffer;
GLdouble currentZ;
//---------------------------------------------------------------------------
//cghw1
//---------------------------------------------------------------------------
GLdouble CTM_MV[16];	//Current Transformation Matrix: ModelView
GLdouble CTM_P[16];		//Current Transformation Matrix: Projection
GLdouble *CTM;

GLdouble VSM[4];    //Viewport Scale Matrix
GLdouble VTM[4];    //Viewport Translation Matrix

M3DMatrix44d MV,P,Inv,T;

typedef struct vector3Struct {
	GLdouble v[3];
} vector3;

typedef struct vector4Struct {
GLdouble v[4];
} vector4;

void normalize(GLdouble *v)
{
    GLdouble length = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] = v[0] / length;
    v[1] = v[1] / length;
    v[2] = v[2] / length;
}

vector3 cross(GLdouble *u, GLdouble *v)
{
    vector3 product;
    product.v[0] = u[1]*v[2]-u[2]*v[1];
    product.v[1] = u[2]*v[0]-u[0]*v[2];
    product.v[2] = u[0]*v[1]-u[1]*v[0];
    
    return product;
}

GLdouble dot(GLdouble *u, GLdouble *v)
{
    GLdouble product;
    product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

    return product;
}

vector4 matrixMultiply41(M3DMatrix44d a, vector4 b)
{
	vector4 product;
	product.v[0] =0;
	product.v[1] =0;
	product.v[2] =0;
	product.v[3] =0;

	for(int j=0; j<4; j++)    //Row j
	 {
		for(int k=0; k<4; k++)
		{
			product.v[j] += a[j+k*4]*b.v[k];
		}
    }

	return product;
}

bool swTransformation(const GLdouble h[4], GLdouble w[4])
{
	//p = CTM_P*CTM_MV*h
	GLdouble tempM[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	GLdouble tempV[4]={0,0,0,0};
	GLdouble v[4]={0,0,0,0};

    //CTM_P*CTM_MV
    for(int i=0;i<4;i++)         //Colum i
     {
          for(int j=0;j<4;j++)    //Row j
          {
               for(int k=0;k<4;k++)
               {
                    tempM[i*4+j] += CTM_P[j+k*4]*CTM_MV[i*4+k];
               }
          }
     }
     //
     for(int j=0;j<4;j++)    //Row j
    {
               for(int k=0;k<4;k++)
               {
                    tempV[j] += tempM[j+k*4]*h[k];
               }
    }
	

	//prespective division
	 for(int k=0;k<2;k++)
	 {
        v[k] = tempV[k] / tempV[3];
	 }
    v[2] = tempV[3];
    v[3] = 1;

	//viewport transformation
	w[0] = (v[0]+1) * VSM[0] + VTM[2];
    w[1] = (v[1]+1) * VSM[3] + VTM[1];

	currentZ = w[2];

	return true;
}

void swTranslated(GLdouble x, GLdouble y, GLdouble z)
{
	GLdouble m[16]={1,0,0,0,
					0,1,0,0,
                    0,0,1,0,
                    x,y,z,1};
    
    swMultMatrixf(m);
}

void swScaled(GLdouble x, GLdouble y, GLdouble z)
{
    GLdouble m[16]={x,0,0,0,
                    0,y,0,0,
                    0,0,z,0,
                    0,0,0,1};
    
    swMultMatrixf(m);
}

void swRotated(GLdouble angle, GLdouble x, GLdouble y, GLdouble z)
{
	GLdouble c = cos(angle*PI/180);
	GLdouble s = sin(angle*PI/180);
   
	GLdouble m[16]={x*x*(1-c)+c,y*x*(1-c)+z*s,x*z*(1-c)-y*s,0,
					x*y*(1-c)-z*s,y*y*(1-c)+c,y*z*(1-c)+x*s,0,
                    x*z*(1-c)+y*s,y*z*(1-c)-x*s,z*z*(1-c)+c,0,
                    0,0,0,1};
    swMultMatrixf(m);
}

void swMatrixMode(GLenum mode)
{
	if(mode == GL_MODELVIEW)
	{
		CTM=CTM_MV;   
	}
	else if(mode == GL_PROJECTION)
	{
		CTM=CTM_P;
	}
	else
	{
		CTM=NULL;
	}
}

void swLoadIdentity(void)
{
	for(int i=0;i<16;i++)
	{
		if(i%5 != 0)
        {
            CTM[i] = 0;
        }
        else
        {
            CTM[i] = 1;      
        }
	}
}

void swLoadMatrixf(const GLdouble * m)
{
	for(int i=0;i<16;i++)
	{
		CTM[i] = m[i];
	}
}

void swMultMatrixf(const GLdouble * m)
{
	GLdouble tempM[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
     
     for(int i=0;i<4;i++)         //Colum i
     {
          for(int j=0;j<4;j++)    //Row j
          {
               for(int k=0;k<4;k++)
               {
                    tempM[i*4+j] += CTM[j+k*4]*m[i*4+k];
               }
          }
     }
     
	for(int i=0;i<16;i++)
	{
		CTM[i] = tempM[i];
	} 
}

void swPushMatrix(void)
{
	vector<GLdouble> tempM;

	for(int i=0;i<16;i++)
	{
		tempM.push_back(CTM[i]);
	}

	matrixStack.push(tempM);
}

void swPopMatrix(void)
{
	vector<GLdouble> tempM = matrixStack.top();

	for(int i=0;i<16;i++)
	{
		CTM[i] = tempM.at(i);
	}

	matrixStack.pop();
}


void swuLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ,
 	           GLdouble centerX, GLdouble centerY, GLdouble centerZ,
 	           GLdouble upX, GLdouble upY, GLdouble upZ)
{
      GLdouble F[3] = {centerX-eyeX,centerY-eyeY,centerZ-eyeZ};
      GLdouble UP[3] = {upX, upY, upZ};
      normalize(F);
      normalize(UP);
      
      vector3 s= cross(F, UP);
      vector3 u= cross(s.v, F);
      
      GLdouble m[16]={s.v[0], u.v[0], -F[0], 0,
                        s.v[1], u.v[1], -F[1], 0,
                        s.v[2], u.v[2], -F[2], 0,
                        0, 0, 0, 1};
                        
        swMultMatrixf(m);
        swTranslated(-eyeX, -eyeY, -eyeZ);
}

void swFrustum(	GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble nearVal, GLdouble farVal)
{
     GLdouble a = (left+right) / (right-left);
     GLdouble b = (bottom+top) / (top-bottom);
     GLdouble c = -(farVal+nearVal) / (farVal-nearVal);
     GLdouble d = -(2*farVal*nearVal) / (farVal-nearVal);
     
     GLdouble m[16]={2*nearVal/(right-left), 0, 0, 0,
                        0, 2*nearVal/(top-bottom), 0, 0,
                        a, b, c, -1,
                        0, 0, d, 0,};
                        
    swMultMatrixf(m);
}

void swuPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)
{
    GLdouble radian = fovy/2*PI/180;
    GLdouble f = 1 / tan(radian);
    
    GLdouble m[16]={f/aspect, 0, 0, 0,
                    0, f, 0, 0,
                    0, 0, (zFar+zNear)/(zNear-zFar), -1,
                    0, 0, 2*zFar*zNear/(zNear-zFar), 0};
    
    swMultMatrixf(m);
}


void swViewport(GLint x, GLint y, GLsizei width, GLsizei height)
{
    VSM[0] = width/2;
    VSM[1] = 0;
    VSM[2] = 0;
    VSM[3] = height/2;
    
    VTM[0] = 1;
    VTM[1] = y;
    VTM[2] = x;
    VTM[3] = 1;

	zbuffer.clear();
	/*vector<GLdouble> row = (width, 0.0);
	zbuffer= (height, 0);*/
	for(int i=0;i<=height;i++)
	{
		vector<GLdouble> row;
		for(int j=0;j<=width;j++)
		{
			row.push_back(0);
		}
		zbuffer.push_back(row);
	}
}

//---------------------------------------------------------------------------
//cghw2
//---------------------------------------------------------------------------
GLdouble interpolate(GLdouble origin,GLdouble goal,GLdouble originY,GLdouble goalY,GLdouble currentY)
{
	GLdouble value;

	value= origin + (goal - origin)*(currentY - originY)/(goalY -originY);
	return value;
}

void writepixel(int x, int y, GLdouble r, GLdouble g, GLdouble b)
{
	GLubyte map[1]={255};

	glColor3d(r, g, b);
	glWindowPos2i(x, y);
	glBitmap(1, 1, 0, 0, 0, 0, map);
}

bool BresenhamLine(int x1, int y1, int x2, int y2, GLdouble r, GLdouble g, GLdouble b)
{
	bool steep = abs(y2 - y1) > abs(x2 - x1);
	if (steep)
	{
		swap(x1, y1);
		swap(x2, y2);
	}
	if (x1 > x2) 
	{
		swap(x1, x2);
		swap(y1, y2);
	}

	GLint deltaX = x2 - x1;
	if(deltaX == 0) writepixel(x1, y1, r, g, b);
	GLint deltaY = abs(y2 - y1);

	GLdouble error = 0;
	GLdouble deltaError =  (double) deltaY /  deltaX;

	GLint ystep;
	GLint x;
	GLint y = y1;
	
	ystep =y1 < y2 ? 1:-1;
	
	for(x=x1; x<=x2; x++)
	{
		if(steep)
		{
			writepixelfast(y, x, r, g, b);
		}
		else
		{
			writepixelfast(x, y, r, g, b);
		}
		error += deltaError;
		if(error >= 0.5)
		{
			y += ystep;
			error -= 1;
		}
	}

	return true;
}

bool swTriangle(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble r, GLdouble g, GLdouble b)
{
	vector< vector<GLdouble> > pointList;
	vector<GLdouble> point;

	point.push_back(x1);
	point.push_back(y1);
	point.push_back(z1);
	pointList.push_back(point);
	
	point.clear();
	point.push_back(x2);
	point.push_back(y2);
	point.push_back(z2);
	pointList.push_back(point);
	
	point.clear();
	point.push_back(x3);
	point.push_back(y3);
	point.push_back(z3);
	pointList.push_back(point);

	for(int i=0;i < pointList.size() - 1;i++)	//Bubble Sort the points according to the y value
	{
		for(int j=0;j < pointList.size() - i -1;j++)
		{
			if( pointList.at(j).at(1) < pointList.at(j+1).at(1) )
			{
				vector<GLdouble>  temp;
				temp = pointList.at(j);
				pointList.at(j) = pointList.at(j+1);
				pointList.at(j+1) = temp;
			}
		}
	}

	GLdouble z;
	//interpolate(origin, goal, originY, goalY, currentY)

	//draw form maximum y value point to the medium
	for(int y=pointList.at(0).at(1); y>=pointList.at(1).at(1); y--)
	{
			x1 = interpolate( pointList.at(0).at(0), pointList.at(1).at(0), pointList.at(0).at(1), pointList.at(1).at(1), y );
			z = interpolate( pointList.at(0).at(2), pointList.at(1).at(2), pointList.at(0).at(1), pointList.at(1).at(1), y );
			
			x2 = interpolate( pointList.at(0).at(0), pointList.at(2).at(0), pointList.at(0).at(1), pointList.at(2).at(1), y );
			z = interpolate( pointList.at(0).at(2), pointList.at(2).at(2), pointList.at(0).at(1), pointList.at(2).at(1), y );
			
			if( x1 - x2 == 0 ) writepixel(x1, y, r, g, b);
			else
			{
				BresenhamLine(x1, y, x2, y, r, g, b);
			}
	}
	
	//draw form minimum y value point to the medium
	for(int y=pointList.at(2).at(1)+1; y<=pointList.at(1).at(1); y++)
	{
			x1 = interpolate( pointList.at(2).at(0), pointList.at(1).at(0), pointList.at(2).at(1), pointList.at(1).at(1), y );
			z = interpolate( pointList.at(2).at(2), pointList.at(1).at(2), pointList.at(2).at(1), pointList.at(1).at(1), y );

			x2 = interpolate( pointList.at(2).at(0), pointList.at(0).at(0), pointList.at(2).at(1), pointList.at(0).at(1), y );
			z = interpolate( pointList.at(2).at(2), pointList.at(0).at(2), pointList.at(2).at(1), pointList.at(0).at(1), y );
			if( x1 - x2 == 0 ) writepixel(x1, y, r, g, b);
			else
			{
				BresenhamLine(x1, y, x2, y, r, g, b);
			}
	}

	return true;
}

bool swInitZbuffer(int width, int height)
{
	
	return true;
}

bool swClearZbuffer()
{
	zbuffer.clear();
	return true;
}

//---------------------------------------------------------------------------
//cghw3
//---------------------------------------------------------------------------
//GLdouble *ZBUFFER;
//int Z_W=-1, Z_H=-1;

GLfloat  _ambientLight[] = { 0.3f, 0.3f, 0.3f, 1.0f };
GLfloat  _diffuseLight[] = { 0.7f, 0.7f, 0.7f, 1.0f };
GLfloat  _specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat  _specref[] = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat  _shininess = 128.0f;

GLfloat _lightPos[] = { -50.f, 50.0f, 100.0f, 1.0f };

vector3 reflect(GLdouble* Incidence, GLdouble* normal)
{
	vector3 reflection;
	reflection.v[0] = 2*dot(Incidence, normal) * normal[0] - Incidence[0];
	reflection.v[1] = 2*dot(Incidence, normal) * normal[1] - Incidence[1];
	reflection.v[2] = 2*dot(Incidence, normal) * normal[2] - Incidence[2];

	return reflection;
}

void writepixelfast(int x, int y, GLdouble r, GLdouble g, GLdouble b)
{
	glBegin(GL_POINTS);
		glColor3d(r, g, b);
		glVertex2i(x, y);
	glEnd();
}

//bool swZbufferResize(int w, int h)
//{
//	//if( w!=Z_W || h!=Z_H ) {
//	//	delete [] ZBUFFER;
//	//	Z_W=w, Z_H=h;
//	//	ZBUFFER = new GLdouble[Z_W*Z_H];
//	//}
//
//	return true;
//}

//GLdouble swZbufferValue(int x, int y)
//{
//	if(x>=0 && x<Z_W && y>=0 && y<Z_H)
//		return ZBUFFER[x+y*Z_W];
//	else
//		return 1.0;
//}
//
//void swClear(GLbitfield mask)
//{
//	for(int i=0; i<Z_W*Z_H; ++i)
//		ZBUFFER[i] = 1.0;
//}

bool swNormalTransformation(const GLdouble h[4], GLdouble w[4])
{
	//apply transformation
	for(int i=0; i < 16; i++)
	{
		MV[i] = CTM_MV[i];
	}
	m3dInvertMatrix44(Inv, MV);
	m3dTransposeMatrix44(T, Inv);

	GLdouble tempV[4]={0,0,0,0};
	 for(int j=0; j<4; j++)    //Row j
    {
               for(int k=0; k<4; k++)
               {
                    tempV[j] += T[j+k*4]*h[k];
               }
    }
	for(int k=0; k<4; k++)
    {
		w[k] = tempV[k];
    }

	//Normalize
	normalize(w);

	return true;
}


void swLightfv(GLenum light, GLenum pname, const GLfloat *params)
{
	switch(pname) {
		case GL_AMBIENT:
			_ambientLight[0]=params[0];
			_ambientLight[1]=params[1];
			_ambientLight[2]=params[2];
			_ambientLight[3]=params[3];
			break;

		case GL_DIFFUSE:
			_diffuseLight[0]=params[0];
			_diffuseLight[1]=params[1];
			_diffuseLight[2]=params[2];
			_diffuseLight[3]=params[3];
			break;

		case GL_SPECULAR:
			_specular[0]=params[0];
			_specular[1]=params[1];
			_specular[2]=params[2];
			_specular[3]=params[3];
			break;

		case GL_POSITION:
			_lightPos[0]=params[0];
			_lightPos[1]=params[1];
			_lightPos[2]=params[2];
			_lightPos[3]=params[3];

			//------------------------------------------------------
			//ADD transformation to eye coordinate
			GLdouble tempV[4]={0,0,0,0};
			for(int j=0;j<4;j++)    //Row j
			{
               for(int k=0;k<4;k++)
               {
                    tempV[j] += CTM_MV[j+k*4]*_lightPos[k];
               }
			}
			_lightPos[0] = tempV[0];
			_lightPos[1] = tempV[1];
			_lightPos[2] = tempV[2];
			_lightPos[3] = tempV[3];

			break;
	}
}

void swMaterialfv (GLenum face, GLenum pname, const GLfloat *params)
{
	switch(pname) {
		case GL_SPECULAR:
			_specref[0]=params[0];
			_specref[1]=params[1];
			_specref[2]=params[2];
			_specref[3]=params[3];
			break;
	}
}

void swMateriali (GLenum face, GLenum pname, GLint param)
{
	switch(pname) {
		case GL_SHININESS:
			_shininess=param;
			break;
	}
}

bool BresenhamLine(int x1, int y1, GLdouble r1, GLdouble g1, GLdouble b1, int x2, int y2, GLdouble r2, GLdouble g2, GLdouble b2)
{
	bool steep = abs(y2 - y1) > abs(x2 - x1);
	if (steep)
	{
		swap(x1, y1);
		swap(x2, y2);
	}
	if (x1 > x2) 
	{
		swap(x1, x2);
		swap(y1, y2);
		swap(r1, r2);
		swap(g1, g2);
		swap(b1, b2);
	}

	GLint deltaX = x2 - x1;
	if(deltaX == 0) writepixelfast(x1, y1, (r1+r2)/2, (g1+g2)/2, (b1+b2)/2);
	else
	{
		GLint deltaY = abs(y2 - y1);
		GLdouble deltaR = r2 - r1;
		GLdouble deltaG = g2 - g1;
		GLdouble deltaB = b2 - b1;

		GLdouble error = 0,errorR = 0,errorG = 0,errorB = 0;
		GLdouble deltaError =  (double) deltaY /  deltaX;
		GLdouble deltaErrorR =  (double) deltaR /  deltaX;
		GLdouble deltaErrorG =  (double) deltaG /  deltaX;
		GLdouble deltaErrorB =  (double) deltaB /  deltaX;

		GLint ystep, rstep, gstep, bstep;
		GLint x;
		GLint y = y1;
		GLdouble r=r1,g=g1,b=b1;
		
		/*ystep = y1 < y2 ? 1:-1;
		rstep = r1 < r2 ? 1:-1;
		gstep = g1 < g2 ? 1:-1;
		bstep = b1 < b2 ? 1:-1;*/
		
		for(x=x1; x<=x2; x++)
		{
			if(steep)
			{
				writepixelfast(y, x, r, g, b);
			}
			else
			{
				writepixelfast(x, y, r, g, b);
			}
			error += deltaError;
			/*errorR += deltaErrorR;
			errorG += deltaErrorG;
			errorB += deltaErrorB;*/
			if(error >= 0.5)
			{
				y += ystep;
				error -= 1;
			}	
			r += deltaErrorR;
			g += deltaErrorG;
			b += deltaErrorB;
			
		}
	}

	return true;
}
//Gouraud shading
bool swTriangleG(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble nx1, GLdouble ny1, GLdouble nz1, 
			 GLdouble nx2, GLdouble ny2, GLdouble nz2, 
			 GLdouble nx3, GLdouble ny3, GLdouble nz3,
			 GLdouble r1, GLdouble g1, GLdouble b1,
			 GLdouble r2, GLdouble g2, GLdouble b2,
			 GLdouble r3, GLdouble g3, GLdouble b3)
{
	//transformation all data(vertex, normal, light vector, COP) to eye coordiante
	//normal
	GLdouble h1[4]={nx1, ny1, nz1};
	GLdouble h2[4]={nx2, ny2, nz2};
	GLdouble h3[4]={nx3, ny3, nz3};
	GLdouble n1[4],n2[4],n3[4];
	GLdouble n[4];
	
	swNormalTransformation(h1, n1);
	swNormalTransformation(h2, n2);
	swNormalTransformation(h3, n3);
	n[0]=n1[0]+n2[0]+n3[0];
	n[1]=n1[1]+n2[1]+n3[1];
	n[2]=n1[2]+n2[2]+n3[2];
	normalize(n);

	//vertex
	vector4 v1,v2,v3;
	v1.v[0] = x1;
	v1.v[1] = y1;
	v1.v[2] = z1;
	v1.v[3] = 1;
	v2.v[0] = x2;
	v2.v[1] = y2;
	v2.v[2] = z2;
	v2.v[3] = 1;
	v3.v[0] = x3;
	v3.v[1] = y3;
	v3.v[2] = z3;
	v3.v[3] = 1;

	/*for(int i=0; i < 16; i++)
	{
		P[i] = CTM_P[i];
	}
	m3dInvertMatrix44(Inv, P);*/

	v1 = matrixMultiply41(MV, v1);
	v2 = matrixMultiply41(MV, v2);
	v3 = matrixMultiply41(MV, v3);
	for(int i =0; i<3; i++)
	{
		v1.v[i] = -v1.v[i];
		v2.v[i] = -v2.v[i];
		v3.v[i] = -v3.v[i];
	}

	vector4 vv1,vv2,vv3;
	vv1 = v1;
	vv2 = v2;
	vv3 = v3;
	normalize(vv1.v);
	normalize(vv2.v);
	normalize(vv3.v);

	//light vector
	vector3 lv1,lv2,lv3;
	lv1.v[0] = _lightPos[0] + v1.v[0];
	lv1.v[1] = _lightPos[1] + v1.v[1];
	lv1.v[2] = _lightPos[2] + v1.v[2];
	normalize(lv1.v);
	lv2.v[0] = _lightPos[0] + v2.v[0];
	lv2.v[1] = _lightPos[1] + v2.v[1];
	lv2.v[2] = _lightPos[2] + v2.v[2];
	normalize(lv2.v);
	lv3.v[0] = _lightPos[0] + v3.v[0];
	lv3.v[1] = _lightPos[1] + v3.v[1];
	lv3.v[2] = _lightPos[2] + v3.v[2];
	normalize(lv3.v);

	//reflection
	vector3 rf1,rf2,rf3;
	rf1 = reflect(lv1.v, n);
	rf2 = reflect(lv2.v, n);
	rf3 = reflect(lv3.v, n);

	//modified Pong shading equation in each vertex
	GLdouble rgb1[3],rgb2[3],rgb3[3];
	//vertex1
	rgb1[0] =  r1 * ( _diffuseLight[0] *  max( dot(lv1.v, n), 0) +  _specref[0] * _specular[0] * pow(  max( dot(vv1.v, rf1.v), 0), (int)_shininess ) + _ambientLight[0] );
	rgb1[1] =  g1 * ( _diffuseLight[1] *  max( dot(lv1.v, n), 0) +  _specref[1] * _specular[1] * pow(  max( dot(vv1.v, rf1.v), 0), (int)_shininess ) + _ambientLight[1] );
	rgb1[2] =  b1 * ( _diffuseLight[2] *  max( dot(lv1.v, n), 0) +  _specref[2] * _specular[2] * pow(  max( dot(vv1.v, rf1.v), 0), (int)_shininess ) + _ambientLight[2] );
	//vertex2
	rgb2[0] =  r2 * ( _diffuseLight[0] *  max( dot(lv2.v, n), 0) +  _specref[0] * _specular[0] * pow(  max( dot(vv2.v, rf2.v), 0), (int)_shininess) + _ambientLight[0] );
	rgb2[1] =  g2 * ( _diffuseLight[1] *  max( dot(lv2.v, n), 0) +  _specref[1] * _specular[1] * pow(  max( dot(vv2.v, rf2.v), 0), (int)_shininess ) + _ambientLight[1] );
	rgb2[2] =  b2 * ( _diffuseLight[2] *  max( dot(lv2.v, n), 0) +  _specref[2] * _specular[2] * pow(  max( dot(vv2.v, rf2.v), 0), (int)_shininess) + _ambientLight[2] );
	//vertex3
	rgb3[0] =  r3 * ( _diffuseLight[0] *  max( dot(lv3.v, n), 0) +  _specref[0] * _specular[0] * pow(  max( dot(vv3.v, rf3.v), 0), (int)_shininess) + _ambientLight[0] );
	rgb3[1] =  g3 * ( _diffuseLight[1] *  max( dot(lv3.v, n), 0) +  _specref[1] * _specular[1] * pow(  max( dot(vv3.v, rf3.v), 0), (int)_shininess) + _ambientLight[1] );
	rgb3[2] =  b3 * ( _diffuseLight[2] *  max( dot(lv3.v, n), 0) +  _specref[2] * _specular[2] * pow(  max( dot(vv3.v, rf3.v), 0), (int)_shininess) + _ambientLight[2] );

	//Raterization: 
	GLdouble w1[4],w2[4],w3[4];
	vector< vector<GLdouble> > pointList;
	vector<GLdouble> point;
	h1[0] = x1;
	h1[1] = y1;
	h1[2] = z1;
	h1[3] = 1;
	h2[0] = x2;
	h2[1] = y2;
	h2[2] = z2;
	h2[3] = 1;
	h3[0] = x3;
	h3[1] = y3;
	h3[2] = z3;
	h3[3] = 1;
	//point1
	point.push_back(x1);
	point.push_back(y1);
	point.push_back(z1);
	point.push_back(rgb1[0]);
	point.push_back(rgb1[1]);
	point.push_back(rgb1[2]);

	swTransformation(h1, w1);
	point.push_back(w1[0]);
	point.push_back(w1[1]);
	pointList.push_back(point);
	point.clear();
	//point2
	point.push_back(x2);
	point.push_back(y2);
	point.push_back(z2);
	point.push_back(rgb2[0]);
	point.push_back(rgb2[1]);
	point.push_back(rgb2[2]);

	swTransformation(h2, w2);
	point.push_back(w2[0]);
	point.push_back(w2[1]);
	pointList.push_back(point);
	point.clear();
	//point3
	point.push_back(x3);
	point.push_back(y3);
	point.push_back(z3);
	point.push_back(rgb3[0]);
	point.push_back(rgb3[1]);
	point.push_back(rgb3[2]);

	swTransformation(h3, w3);
	point.push_back(w3[0]);
	point.push_back(w3[1]);
	pointList.push_back(point);

	for(int i=0;i < pointList.size() - 1;i++)	//Bubble Sort the points according to the y value
	{
		for(int j=0;j < pointList.size() - i -1;j++)
		{
			if( pointList.at(j).at(7) < pointList.at(j+1).at(7) )
			{
				vector<GLdouble>  temp;
				temp = pointList.at(j);
				pointList.at(j) = pointList.at(j+1);
				pointList.at(j+1) = temp;
			}
		}
	}

	GLdouble z;
	//interpolate(origin, goal, originY, goalY, currentY)

	//draw form maximum y value point to the medium
	for(int y=pointList.at(0).at(7)-1; y>=pointList.at(1).at(7); y--)
	{
			x1 = interpolate( pointList.at(0).at(6), pointList.at(1).at(6), pointList.at(0).at(7), pointList.at(1).at(7), y );
			x2 = interpolate( pointList.at(0).at(6), pointList.at(2).at(6), pointList.at(0).at(7), pointList.at(2).at(7), y );
			GLdouble eyeY = interpolate( pointList.at(0).at(1), pointList.at(1).at(1), pointList.at(0).at(7), pointList.at(1).at(7), y);
				r1 = interpolate( pointList.at(0).at(3), pointList.at(1).at(3), pointList.at(0).at(7), pointList.at(1).at(7), y );
				g1 = interpolate( pointList.at(0).at(4), pointList.at(1).at(4), pointList.at(0).at(7), pointList.at(1).at(7), y );
				b1 = interpolate( pointList.at(0).at(5), pointList.at(1).at(5), pointList.at(0).at(7), pointList.at(1).at(7), y );
				/*r1 = interpolate( pointList.at(0).at(3), pointList.at(1).at(3), pointList.at(0).at(1), pointList.at(1).at(1), eyeY );
				g1 = interpolate( pointList.at(0).at(4), pointList.at(1).at(4), pointList.at(0).at(1), pointList.at(1).at(1), eyeY );
				b1 = interpolate( pointList.at(0).at(5), pointList.at(1).at(5), pointList.at(0).at(1), pointList.at(1).at(1), eyeY );*/

			eyeY = interpolate( pointList.at(0).at(1), pointList.at(2).at(1), pointList.at(0).at(7), pointList.at(2).at(7), y);
				r2 = interpolate( pointList.at(0).at(3), pointList.at(2).at(3), pointList.at(0).at(7), pointList.at(2).at(7), y );
				g2 = interpolate( pointList.at(0).at(4), pointList.at(2).at(4), pointList.at(0).at(7), pointList.at(2).at(7), y );
				b2 = interpolate( pointList.at(0).at(5), pointList.at(2).at(5), pointList.at(0).at(7), pointList.at(2).at(7), y );
				/*r2 = interpolate( pointList.at(0).at(3), pointList.at(2).at(3), pointList.at(0).at(1), pointList.at(2).at(1), eyeY );
				g2 = interpolate( pointList.at(0).at(4), pointList.at(2).at(4), pointList.at(0).at(1), pointList.at(2).at(1), eyeY );
				b2 = interpolate( pointList.at(0).at(5), pointList.at(2).at(5), pointList.at(0).at(1), pointList.at(2).at(1), eyeY );*/

			if( x1 -x2 == 0 )
			{
				GLdouble r = ( r1 + r2 ) / 2;
				GLdouble g = ( g1 + g2 ) / 2;
				GLdouble b = ( b1 + b2 ) / 2;
				writepixelfast(x1, y, r, g, b);
			}
			else
			{		
				BresenhamLine(x1, y, r1, g1, b1, x2, y, r2, g2, b2);
			}
	}
	
	//draw form minimum y value point to the medium
	for(int y=pointList.at(2).at(7)+1; y<=pointList.at(1).at(7); y++)
	{
			x1 = interpolate( pointList.at(2).at(6), pointList.at(1).at(6), pointList.at(2).at(7), pointList.at(1).at(7), y );
			x2 = interpolate( pointList.at(2).at(6), pointList.at(0).at(6), pointList.at(2).at(7), pointList.at(0).at(7), y );
			GLdouble eyeY = interpolate( pointList.at(2).at(1), pointList.at(1).at(1), pointList.at(2).at(7), pointList.at(1).at(7), y);
				r1 = interpolate( pointList.at(2).at(3), pointList.at(1).at(3), pointList.at(2).at(7), pointList.at(1).at(7), y );
				g1 = interpolate( pointList.at(2).at(4), pointList.at(1).at(4), pointList.at(2).at(7), pointList.at(1).at(7), y );
				b1 = interpolate( pointList.at(2).at(5), pointList.at(1).at(5), pointList.at(2).at(7), pointList.at(1).at(7), y );
				/*r1 = interpolate( pointList.at(2).at(3), pointList.at(1).at(3), pointList.at(2).at(1), pointList.at(1).at(1), eyeY );
				g1 = interpolate( pointList.at(2).at(4), pointList.at(1).at(4), pointList.at(2).at(1), pointList.at(1).at(1), eyeY );
				b1 = interpolate( pointList.at(2).at(5), pointList.at(1).at(5), pointList.at(2).at(1), pointList.at(1).at(1), eyeY );*/

			eyeY = interpolate( pointList.at(2).at(1), pointList.at(0).at(1), pointList.at(2).at(7), pointList.at(0).at(7), y);
				r2 = interpolate( pointList.at(2).at(3), pointList.at(0).at(3), pointList.at(2).at(7), pointList.at(0).at(7), y );
				g2 = interpolate( pointList.at(2).at(4), pointList.at(0).at(4), pointList.at(2).at(7), pointList.at(0).at(7), y );
				b2 = interpolate( pointList.at(2).at(5), pointList.at(0).at(5), pointList.at(2).at(7), pointList.at(0).at(7), y );
				/*r2 = interpolate( pointList.at(2).at(3), pointList.at(0).at(3), pointList.at(2).at(1), pointList.at(0).at(1), eyeY );
				g2 = interpolate( pointList.at(2).at(4), pointList.at(0).at(4), pointList.at(2).at(1), pointList.at(0).at(1), eyeY );
				b2 = interpolate( pointList.at(2).at(5), pointList.at(0).at(5), pointList.at(2).at(1), pointList.at(0).at(1), eyeY );*/
			if( x1 -x2 == 0 )
			{
				GLdouble r = ( r1 + r2 ) / 2;
				GLdouble g = ( g1 + g2 ) / 2;
				GLdouble b = ( b1 + b2 ) / 2;
				writepixelfast(x1, y, r, g, b);
			}
			else
			{		
				BresenhamLine(x1, y, r1, g1, b1, x2, y, r2, g2, b2);
			}
	}

	return true;
}


////Phong Shading
bool swTriangleP(GLdouble x1, GLdouble y1, GLdouble z1, 
			 GLdouble x2, GLdouble y2, GLdouble z2, 
			 GLdouble x3, GLdouble y3, GLdouble z3,
			 GLdouble nx1, GLdouble ny1, GLdouble nz1, 
			 GLdouble nx2, GLdouble ny2, GLdouble nz2, 
			 GLdouble nx3, GLdouble ny3, GLdouble nz3,
			 GLdouble r1, GLdouble g1, GLdouble b1,
			 GLdouble r2, GLdouble g2, GLdouble b2,
			 GLdouble r3, GLdouble g3, GLdouble b3)
{
	//transformation all data(vertex, normal, light vector, COP) to eye coordiante
	//normal
	GLdouble h1[4]={nx1, ny1, nz1};
	GLdouble h2[4]={nx2, ny2, nz2};
	GLdouble h3[4]={nx3, ny3, nz3};
	GLdouble n1[4],n2[4],n3[4];
	GLdouble n[4];
	
	swNormalTransformation(h1, n1);
	swNormalTransformation(h2, n2);
	swNormalTransformation(h3, n3);
	n[0]=n1[0]+n2[0]+n3[0];
	n[1]=n1[1]+n2[1]+n3[1];
	n[2]=n1[2]+n2[2]+n3[2];
	normalize(n);

	//vertex
	vector4 v1,v2,v3;
	v1.v[0] = x1;
	v1.v[1] = y1;
	v1.v[2] = z1;
	v1.v[3] = 1;
	v2.v[0] = x2;
	v2.v[1] = y2;
	v2.v[2] = z2;
	v2.v[3] = 1;
	v3.v[0] = x3;
	v3.v[1] = y3;
	v3.v[2] = z3;
	v3.v[3] = 1;

	v1 = matrixMultiply41(MV, v1);
	v2 = matrixMultiply41(MV, v2);
	v3 = matrixMultiply41(MV, v3);

	vector4 vv1,vv2,vv3;
	vv1 = v1;
	vv2 = v2;
	vv3 = v3;
	normalize(vv1.v);
	normalize(vv2.v);
	normalize(vv3.v);

	//light vector
	vector3 lv1,lv2,lv3;
	lv1.v[0] = _lightPos[0] - v1.v[0];
	lv1.v[1] = _lightPos[1] - v1.v[1];
	lv1.v[2] = _lightPos[2] - v1.v[2];
	normalize(lv1.v);
	lv2.v[0] = _lightPos[0] - v2.v[0];
	lv2.v[1] = _lightPos[1] - v2.v[1];
	lv2.v[2] = _lightPos[2] - v2.v[2];
	normalize(lv2.v);
	lv3.v[0] = _lightPos[0] - v3.v[0];
	lv3.v[1] = _lightPos[1] - v3.v[1];
	lv3.v[2] = _lightPos[2] - v3.v[2];
	normalize(lv3.v);

	//reflection
	vector3 rf1,rf2,rf3;
	rf1 = reflect(lv1.v, n);
	rf2 = reflect(lv2.v, n);
	rf3 = reflect(lv3.v, n);

	//modified Pong shading equation in each vertex
	GLdouble rgb1[3],rgb2[3],rgb3[3];
	//vertex1
	rgb1[0] =  r1 * ( _diffuseLight[0] *  max( dot(lv1.v, n), 0) +  _specref[0] * _specular[0] * pow(  max( dot(vv1.v, rf1.v), 0), (int)_shininess ) + _ambientLight[0] );
	rgb1[1] =  g1 * ( _diffuseLight[1] *  max( dot(lv1.v, n), 0) +  _specref[1] * _specular[1] * pow(  max( dot(vv1.v, rf1.v), 0), (int)_shininess ) + _ambientLight[1] );
	rgb1[2] =  b1 * ( _diffuseLight[2] *  max( dot(lv1.v, n), 0) +  _specref[2] * _specular[2] * pow(  max( dot(vv1.v, rf1.v), 0), (int)_shininess ) + _ambientLight[2] );
	//vertex2
	rgb2[0] =  r2 * ( _diffuseLight[0] *  max( dot(lv2.v, n), 0) +  _specref[0] * _specular[0] * pow(  max( dot(vv2.v, rf2.v), 0), (int)_shininess) + _ambientLight[0] );
	rgb2[1] =  g2 * ( _diffuseLight[1] *  max( dot(lv2.v, n), 0) +  _specref[1] * _specular[1] * pow(  max( dot(vv2.v, rf2.v), 0), (int)_shininess ) + _ambientLight[1] );
	rgb2[2] =  b2 * ( _diffuseLight[2] *  max( dot(lv2.v, n), 0) +  _specref[2] * _specular[2] * pow(  max( dot(vv2.v, rf2.v), 0), (int)_shininess) + _ambientLight[2] );
	//vertex3
	rgb3[0] =  r3 * ( _diffuseLight[0] *  max( dot(lv3.v, n), 0) +  _specref[0] * _specular[0] * pow(  max( dot(vv3.v, rf3.v), 0), (int)_shininess) + _ambientLight[0] );
	rgb3[1] =  g3 * ( _diffuseLight[1] *  max( dot(lv3.v, n), 0) +  _specref[1] * _specular[1] * pow(  max( dot(vv3.v, rf3.v), 0), (int)_shininess) + _ambientLight[1] );
	rgb3[2] =  b3 * ( _diffuseLight[2] *  max( dot(lv3.v, n), 0) +  _specref[2] * _specular[2] * pow(  max( dot(vv3.v, rf3.v), 0), (int)_shininess) + _ambientLight[2] );

	//Raterization: 
	GLdouble w1[4],w2[4],w3[4];
	vector< vector<GLdouble> > pointList;
	vector<GLdouble> point;
	h1[0] = x1;
	h1[1] = y1;
	h1[2] = z1;
	h1[3] = 1;
	h2[0] = x2;
	h2[1] = y2;
	h2[2] = z2;
	h2[3] = 1;
	h3[0] = x3;
	h3[1] = y3;
	h3[2] = z3;
	h3[3] = 1;
	//point1
	point.push_back(x1);
	point.push_back(y1);
	point.push_back(z1);
	point.push_back(rgb1[0]);
	point.push_back(rgb1[1]);
	point.push_back(rgb1[2]);

	swTransformation(h1, w1);
	point.push_back(w1[0]);
	point.push_back(w1[1]);
	pointList.push_back(point);
	point.clear();
	//point2
	point.push_back(x2);
	point.push_back(y2);
	point.push_back(z2);
	point.push_back(rgb2[0]);
	point.push_back(rgb2[1]);
	point.push_back(rgb2[2]);

	swTransformation(h2, w2);
	point.push_back(w2[0]);
	point.push_back(w2[1]);
	pointList.push_back(point);
	point.clear();
	//point3
	point.push_back(x3);
	point.push_back(y3);
	point.push_back(z3);
	point.push_back(rgb3[0]);
	point.push_back(rgb3[1]);
	point.push_back(rgb3[2]);

	swTransformation(h3, w3);
	point.push_back(w3[0]);
	point.push_back(w3[1]);
	pointList.push_back(point);

	for(int i=0;i < pointList.size() - 1;i++)	//Bubble Sort the points according to the y value
	{
		for(int j=0;j < pointList.size() - i -1;j++)
		{
			if( pointList.at(j).at(7) < pointList.at(j+1).at(7) )
			{
				vector<GLdouble>  temp;
				temp = pointList.at(j);
				pointList.at(j) = pointList.at(j+1);
				pointList.at(j+1) = temp;
			}
		}
	}

	GLdouble z;
	//interpolate(origin, goal, originY, goalY, currentY)

	//draw form maximum y value point to the medium
	for(int y=pointList.at(0).at(7)-1; y>=pointList.at(1).at(7); y--)
	{
			x1 = interpolate( pointList.at(0).at(6), pointList.at(1).at(6), pointList.at(0).at(7), pointList.at(1).at(7), y );
			x2 = interpolate( pointList.at(0).at(6), pointList.at(2).at(6), pointList.at(0).at(7), pointList.at(2).at(7), y );
		
				r1 = interpolate( pointList.at(0).at(3), pointList.at(1).at(3), pointList.at(0).at(7), pointList.at(1).at(7), y );
				g1 = interpolate( pointList.at(0).at(4), pointList.at(1).at(4), pointList.at(0).at(7), pointList.at(1).at(7), y );
				b1 = interpolate( pointList.at(0).at(5), pointList.at(1).at(5), pointList.at(0).at(7), pointList.at(1).at(7), y );
			
				r2 = interpolate( pointList.at(0).at(3), pointList.at(2).at(3), pointList.at(0).at(7), pointList.at(2).at(7), y );
				g2 = interpolate( pointList.at(0).at(4), pointList.at(2).at(4), pointList.at(0).at(7), pointList.at(2).at(7), y );
				b2 = interpolate( pointList.at(0).at(5), pointList.at(2).at(5), pointList.at(0).at(7), pointList.at(2).at(7), y );
				

			if( x1 -x2 == 0 )
			{
				GLdouble r = ( r1 + r2 ) / 2;
				GLdouble g = ( g1 + g2 ) / 2;
				GLdouble b = ( b1 + b2 ) / 2;
				writepixelfast(x1, y, r, g, b);
			}
			else
			{		
				BresenhamLine(x1, y, r1, g1, b1, x2, y, r2, g2, b2);
			}
	}
	
	//draw form minimum y value point to the medium
	for(int y=pointList.at(2).at(7)+1; y<=pointList.at(1).at(7); y++)
	{
			x1 = interpolate( pointList.at(2).at(6), pointList.at(1).at(6), pointList.at(2).at(7), pointList.at(1).at(7), y );
			x2 = interpolate( pointList.at(2).at(6), pointList.at(0).at(6), pointList.at(2).at(7), pointList.at(0).at(7), y );

				r1 = interpolate( pointList.at(2).at(3), pointList.at(1).at(3), pointList.at(2).at(7), pointList.at(1).at(7), y );
				g1 = interpolate( pointList.at(2).at(4), pointList.at(1).at(4), pointList.at(2).at(7), pointList.at(1).at(7), y );
				b1 = interpolate( pointList.at(2).at(5), pointList.at(1).at(5), pointList.at(2).at(7), pointList.at(1).at(7), y );
			
				r2 = interpolate( pointList.at(2).at(3), pointList.at(0).at(3), pointList.at(2).at(7), pointList.at(0).at(7), y );
				g2 = interpolate( pointList.at(2).at(4), pointList.at(0).at(4), pointList.at(2).at(7), pointList.at(0).at(7), y );
				b2 = interpolate( pointList.at(2).at(5), pointList.at(0).at(5), pointList.at(2).at(7), pointList.at(0).at(7), y );

			if( x1 -x2 == 0 )
			{
				GLdouble r = ( r1 + r2 ) / 2;
				GLdouble g = ( g1 + g2 ) / 2;
				GLdouble b = ( b1 + b2 ) / 2;
				writepixelfast(x1, y, r, g, b);
			}
			else
			{		
				BresenhamLine(x1, y, r1, g1, b1, x2, y, r2, g2, b2);
			}
	} 

	return true;
}




