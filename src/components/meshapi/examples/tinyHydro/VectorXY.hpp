// 2D vector for XY geometry
// Fri Nov 21 10:50:53 PST 2014
#include <math.h>

#ifndef	_VECTORXY_H
#define	_VECTORXY_H 1

class VectorXY
{
 public:
   VectorXY(){x=0.0; y=0.0;}
   VectorXY(double xx, double yy){x=xx; y=yy;}

   double x;
   double y;
   
//----------------------------------------------

VectorXY add(const VectorXY a) const
{
   return VectorXY(x+a.x, y+a.y);
}

//----------------------------------------------

VectorXY operator+(const VectorXY a) const
{
   return VectorXY(x+a.x, y+a.y);
}

//----------------------------------------------

void accum(const VectorXY & b)
{
   x += b.x;
   y += b.y;
}

//----------------------------------------------

void operator+=(const VectorXY & b)
{
   x += b.x;
   y += b.y;
}

//----------------------------------------------

VectorXY sub(const VectorXY b) const
{
   return VectorXY(x-b.x, y-b.y);
}

//----------------------------------------------

VectorXY operator-(const VectorXY & b) const
{
   return VectorXY(x-b.x, y-b.y);
}

//----------------------------------------------

void elim(const VectorXY b)
{
   x -= b.x;
   y -= b.y;
}

//----------------------------------------------

VectorXY mult(double s) const 
{
   return VectorXY(s*x, s*y);
}

//----------------------------------------------

void scale( double s)
{
   x *= s;
   y *= s;
}

//----------------------------------------------
void operator*=(const double s)
{
   x *= s;
   y *= s;
}

//----------------------------------------------

double dot(const VectorXY & v) const
{
   return x*v.x + y*v.y;
}

//----------------------------------------------

double cross(const VectorXY & v) const
{
   return x*v.y - y*v.x;
}

//----------------------------------------------

 double mag(void) const
{
   return sqrt(x*x + y*y);
}

//----------------------------------------------

double mag2(void) const
{
   return x*x + y*y;
}

//----------------------------------------------
};

#endif
