#ifndef _OPERATORS_H_
#define _OPERATORS_H_

vector add(vector a, vector b)
{
  vector out;
  out.x = a.x + b.x;
  out.y = a.y + b.y;
  out.z = a.z + b.z;
  return out;
}

vector subtract(vector a, vector b)
{
  vector out;
  out.x = a.x - b.x;
  out.y = a.y - b.y;
  out.z = a.z - b.z;
  return out;
}

vector multiply(vector a, double b)
{
  vector out;
  out.x = a.x * b;
  out.y = a.y * b;
  out.z = a.z * b;
  return out;
}

vector divide(vector a, double b)
{
  vector out;
  out.x = a.x / b;
  out.y = a.y / b;
  out.z = a.z / b;
  return out;
}

double magnitude(vector a)
{
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

vector unit(vector a)
{
  double mag = magnitude(a);
  return vector(a.x/mag, a.y/mag, a.z/mag);
}

double dot(vector a,vector b)
{
  double dot;
  dot = (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
  return dot;
}

vector cross(vector a,vector b)
{
  vector cp;
  cp.x = a.y * b.z - b.y * a.z;
  cp.y = -(a.x * b.z - a.z * b.x);
  cp.z = a.x * b.y - a.y * b.x ;
  return cp;
}

void copy(vector a, vector b)
{
  b.x = a.x;
  b.y = a.y;
  b.z = a.z;
}

#endif
