#include <sstream>
#include <fstream>
#include "geometry.h"

bool loadGeoFile(const char *file,
		 uint32_t &numFaces,
		 Vec3f* &verts,
		 Vec2f* &st,
		 uint32_t* &vertsIndex) {
  std::ifstream ifs;
  ifs.open(file);
  if (ifs.fail())
    return false;
  
  std::stringstream ss;
  ss << ifs.rdbuf();
  ss >> numFaces;
  uint32_t vertsIndexArraySize = 0;
  // reading face index array
  for (uint32_t i = 0; i < numFaces; ++i) {
    uint32_t tmp;
    ss >> tmp; //faceIndex[i];
    vertsIndexArraySize += tmp; //faceIndex[i];
  }
  vertsIndex = new uint32_t[vertsIndexArraySize];
  uint32_t vertsArraySize = 0;
  // reading vertex index array
  for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
    ss >> vertsIndex[i];
    if (vertsIndex[i] > vertsArraySize) vertsArraySize = vertsIndex[i];
  }
  vertsArraySize += 1;
  // reading vertices
  verts = new Vec3f[vertsArraySize];
  for (uint32_t i = 0; i < vertsArraySize; ++i) {
    ss >> verts[i].x >> verts[i].y >> verts[i].z;
  }
  // reading normals
  for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
    Vec3f normal;
    ss >> normal.x >> normal.y >> normal.z;
  }
  // reading st coordinates
  st = new Vec2f[vertsIndexArraySize];
  for (uint32_t i = 0; i < vertsIndexArraySize; ++i) {
    ss >> st[i].x >> st[i].y;
  }

  ifs.close();
  return true;
}
