#include "geometry.h"

void Matrix44f::multVecMatrix(const Vec3<float> &src, Vec3<float> &dst) const
{
  float a, b, c, w;
  
  a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0] + x[3][0];
  b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1] + x[3][1];
  c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2] + x[3][2];
  w = src[0] * x[0][3] + src[1] * x[1][3] + src[2] * x[2][3] + x[3][3];
  
  dst.x = a / w;
  dst.y = b / w;
  dst.z = c / w;
}

Matrix44f Matrix44f::inverse() const
{
  int i, j, k;
  Matrix44f s;
  Matrix44f t (*this);
  
  // Forward elimination
  for (i = 0; i < 3 ; i++) {
    int pivot = i;
    
    float pivotsize = t[i][i];
    
    if (pivotsize < 0)
      pivotsize = -pivotsize;
    
    for (j = i + 1; j < 4; j++) {
      float tmp = t[j][i];
      
      if (tmp < 0)
	tmp = -tmp;
      
      if (tmp > pivotsize) {
	pivot = j;
	pivotsize = tmp;
      }
    }
    
    if (pivotsize == 0) {
      // Cannot invert singular matrix
      return Matrix44f();
    }
    
    if (pivot != i) {
      for (j = 0; j < 4; j++) {
	float tmp;
        
	tmp = t[i][j];
	t[i][j] = t[pivot][j];
	t[pivot][j] = tmp;
        
	tmp = s[i][j];
	s[i][j] = s[pivot][j];
	s[pivot][j] = tmp;
      }
    }
    
    for (j = i + 1; j < 4; j++) {
      float f = t[j][i] / t[i][i];
      
      for (k = 0; k < 4; k++) {
	t[j][k] -= f * t[i][k];
	s[j][k] -= f * s[i][k];
      }
    }
  }
  
  // Backward substitution
  for (i = 3; i >= 0; --i) {
    float f;
    
    if ((f = t[i][i]) == 0) {
      // Cannot invert singular matrix
      return Matrix44f();
    }
    
    for (j = 0; j < 4; j++) {
      t[i][j] /= f;
      s[i][j] /= f;
    }
    
    for (j = 0; j < i; j++) {
      f = t[j][i];
      
      for (k = 0; k < 4; k++) {
	t[j][k] -= f * t[i][k];
	s[j][k] -= f * s[i][k];
      }
    }
  }
  
  return s;
}

  
