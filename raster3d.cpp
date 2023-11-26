//[heaeder]
// A practical implementation of the rasterization algorithm.
//[/header]
//[compile]
// Download the raster3d.cpp, cow.h and geometry.h files to the same folder.
// Open a shell/terminal, and run the following command where the files are saved:
//
// c++ -o raster3d raster3d.cpp  -std=c++11 -O3
//
// Run with: ./raster3d. Open the file ./output.png in Photoshop or any program
// reading PPM files.
//[/compile]
//[ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//[/ignore]


#include <fstream>
#include <chrono>
#include <cassert>

#include <SDL2/SDL.h>
#include <unistd.h>

#include "geometry.h"
#include "cow.h"

static const float inchToMm = 25.4;
enum FitResolutionGate { kFill = 0, kOverscan };

const uint32_t imageWidth = 640*2;
const uint32_t imageHeight = 480*2;

static SDL_Window *sdlwin = nullptr;
static SDL_Surface *sdlscr = nullptr;

struct pixel_t {
  uint8_t r;
  uint8_t g;
  uint8_t b;
  uint8_t a;
};

template<typename T>
void incr_count(T *ptr, T amt) {
  *ptr += amt;
}

//[comment]
// Compute screen coordinates based on a physically-based camera model
// http://www.scratchapixel.com/lessons/3d-basic-rendering/3d-viewing-pinhole-camera
//[/comment]
void computeScreenCoordinates(
    const float &filmApertureWidth,
    const float &filmApertureHeight,
    const uint32_t &imageWidth,
    const uint32_t &imageHeight,
    const FitResolutionGate &fitFilm,
    const float &nearClippingPLane,
    const float &focalLength,
    float &top, float &bottom, float &left, float &right
)
{
    float filmAspectRatio = filmApertureWidth / filmApertureHeight;
    float deviceAspectRatio = imageWidth / (float)imageHeight;
    
    top = ((filmApertureHeight * inchToMm / 2) / focalLength) * nearClippingPLane;
    right = ((filmApertureWidth * inchToMm / 2) / focalLength) * nearClippingPLane;

    // field of view (horizontal)
    float fov = 2 * 180 / M_PI * atan((filmApertureWidth * inchToMm / 2) / focalLength);
    std::cerr << "Field of view " << fov << std::endl;
    
    float xscale = 1;
    float yscale = 1;
    
    switch (fitFilm) {
        default:
        case kFill:
            if (filmAspectRatio > deviceAspectRatio) {
                xscale = deviceAspectRatio / filmAspectRatio;
            }
            else {
                yscale = filmAspectRatio / deviceAspectRatio;
            }
            break;
        case kOverscan:
            if (filmAspectRatio > deviceAspectRatio) {
                yscale = filmAspectRatio / deviceAspectRatio;
            }
            else {
                xscale = deviceAspectRatio / filmAspectRatio;
            }
            break;
    }
    
    right *= xscale;
    top *= yscale;
    
    bottom = -top;
    left = -right;
}

//[comment]
// Compute vertex raster screen coordinates.
// Vertices are defined in world space. They are then converted to camera space,
// then to NDC space (in the range [-1,1]) and then to raster space.
// The z-coordinates of the vertex in raster space is set with the z-coordinate
// of the vertex in camera space.
//[/comment]
void convertToRaster(
    const Vec3f &vertexWorld,
    const Matrix44f &worldToCamera,
    const float &l,
    const float &r,
    const float &t,
    const float &b,
    const float &near,
    const uint32_t &imageWidth,
    const uint32_t &imageHeight,
    Vec3f &vertexRaster
)
{
    Vec3f vertexCamera;

    worldToCamera.multVecMatrix(vertexWorld, vertexCamera);
    
    // convert to screen space
    Vec2f vertexScreen;
    vertexScreen.x = near * vertexCamera.x / -vertexCamera.z;
    vertexScreen.y = near * vertexCamera.y / -vertexCamera.z;
    
    // now convert point from screen space to NDC space (in range [-1,1])
    Vec2f vertexNDC;
    vertexNDC.x = 2 * vertexScreen.x / (r - l) - (r + l) / (r - l);
    vertexNDC.y = 2 * vertexScreen.y / (t - b) - (t + b) / (t - b);

    // convert to raster space
    vertexRaster.x = (vertexNDC.x + 1) / 2 * imageWidth;
    // in raster space y is down so invert direction
    vertexRaster.y = (1 - vertexNDC.y) / 2 * imageHeight;
    vertexRaster.z = -vertexCamera.z;
}

float min3(const float &a, const float &b, const float &c)
{ return std::min(a, std::min(b, c)); }

float max3(const float &a, const float &b, const float &c)
{ return std::max(a, std::max(b, c)); }

float edgeFunction(const Vec3f &a, const Vec3f &b, const Vec3f &c) {
  float dx = b[0] - a[0], dy = b[1] - a[1];
  return (c[0] - a[0]) * dy - (c[1] - a[1]) * dx;
}


const Matrix44f worldToCamera = {
  0.707107, -0.331295, 0.624695, 0,
  0, 0.883452, 0.468521, 0,
  -0.707107, -0.331295, 0.624695, 0,
  -1.63871, -5.747777, -40.400412, 1};

const uint32_t ntris = (sizeof(stindices)/sizeof(stindices[0]))/3;
							 

const float nearClippingPLane = 1.0f, farClippingPLane = 1000.0f;
float filmApertureWidth = 0.980f, filmApertureHeight = 0.735f, focalLength = 20; // in mm;

int main(int argc, char **argv)
{

  sdlwin = SDL_CreateWindow("FRAMEBUFFER",
			    SDL_WINDOWPOS_UNDEFINED,
			    SDL_WINDOWPOS_UNDEFINED,
			    imageWidth,
			    imageHeight,
			    SDL_WINDOW_SHOWN);
  assert(sdlwin != nullptr);
  sdlscr = SDL_GetWindowSurface(sdlwin);
  assert(sdlscr);

  Matrix44f cameraToWorld = worldToCamera.inverse();
  
    
  // compute screen coordinates
  float t, b, l, r;
  
  computeScreenCoordinates(
			   filmApertureWidth,
			   filmApertureHeight,
			   imageWidth,
			   imageHeight,
			   kOverscan,
			   nearClippingPLane,
			   focalLength,
			   t, b, l, r);
    
    
    float *depthBuffer = new float[imageWidth * imageHeight];
    uint64_t n_tests = 0, n_passed_area = 0, n_passed_z = 0;
    int n_frames = 0;
    while(n_frames < 30) {
      n_tests  = n_passed_area = n_passed_z = 0;
      n_frames++;
      
      for (uint32_t i = 0; i < imageWidth * imageHeight; ++i) {
	depthBuffer[i] = farClippingPLane;
      }
      auto t_start = std::chrono::high_resolution_clock::now();
    
    // [comment]
    // Outer loop
    // [/comment]
    //

    SDL_LockSurface(sdlscr);
    pixel_t *pxls = reinterpret_cast<pixel_t*>(sdlscr->pixels);
    
    //#pragma omp parallel for
    for (uint32_t i = 0; i < ntris; ++i) {
        const Vec3f &v0 = vertices[nvertices[i * 3]];
        const Vec3f &v1 = vertices[nvertices[i * 3 + 1]];
        const Vec3f &v2 = vertices[nvertices[i * 3 + 2]];
        
        // [comment]
        // Convert the vertices of the triangle to raster space
        // [/comment]
        Vec3f v0Raster, v1Raster, v2Raster;
        convertToRaster(v0, worldToCamera, l, r, t, b, nearClippingPLane, imageWidth, imageHeight, v0Raster);
        convertToRaster(v1, worldToCamera, l, r, t, b, nearClippingPLane, imageWidth, imageHeight, v1Raster);
        convertToRaster(v2, worldToCamera, l, r, t, b, nearClippingPLane, imageWidth, imageHeight, v2Raster);
        
        // [comment]
        // Precompute reciprocal of vertex z-coordinate
        // [/comment]
        v0Raster.z = 1.0f / v0Raster.z,
        v1Raster.z = 1.0f / v1Raster.z,
        v2Raster.z = 1.0f / v2Raster.z;
        
        
        // [comment]
        // Prepare vertex attributes. Divde them by their vertex z-coordinate
        // (though we use a multiplication here because v.z = 1 / v.z)
        // [/comment]
        Vec2f st0 = st[stindices[i * 3]];
        Vec2f st1 = st[stindices[i * 3 + 1]];
        Vec2f st2 = st[stindices[i * 3 + 2]];

        st0 *= v0Raster.z, st1 *= v1Raster.z, st2 *= v2Raster.z;
    
        float xmin = min3(v0Raster.x, v1Raster.x, v2Raster.x);
        float ymin = min3(v0Raster.y, v1Raster.y, v2Raster.y);
        float xmax = max3(v0Raster.x, v1Raster.x, v2Raster.x);
        float ymax = max3(v0Raster.y, v1Raster.y, v2Raster.y);
        
        // the triangle is out of screen
        if (xmin > (imageWidth - 1) || xmax < 0 || ymin > (imageHeight - 1) || ymax < 0)
	  continue;

        // be careful xmin/xmax/ymin/ymax can be negative. Don't cast to uint32_t
        uint32_t x0 = std::max(int32_t(0), (int32_t)(std::floor(xmin)));
        uint32_t x1 = std::min(int32_t(imageWidth) - 1, (int32_t)(std::floor(xmax)));
        uint32_t y0 = std::max(int32_t(0), (int32_t)(std::floor(ymin)));
        uint32_t y1 = std::min(int32_t(imageHeight) - 1, (int32_t)(std::floor(ymax)));

	float recip_area = 1.0f / edgeFunction(v0Raster, v1Raster, v2Raster);

	Vec3f v0Cam, v1Cam, v2Cam;
	worldToCamera.multVecMatrix(v0, v0Cam);
	worldToCamera.multVecMatrix(v1, v1Cam);
	worldToCamera.multVecMatrix(v2, v2Cam);
	Vec3f n = (v1Cam - v0Cam).crossProduct(v2Cam - v0Cam).normalize();

	float v0Cam_xz = (v0Cam.x/-v0Cam.z);
	float v0Cam_yz = (v0Cam.y/-v0Cam.z);
	float v1Cam_xz = (v1Cam.x/-v1Cam.z);
	float v1Cam_yz = (v1Cam.y/-v1Cam.z);
	float v2Cam_xz = (v2Cam.x/-v2Cam.z);
	float v2Cam_yz = (v2Cam.y/-v2Cam.z);		
	
	
	uint32_t ylen = y1-y0, xlen = x1-x0;
	Vec3f pp_x0y0(x0 + 0.5, y0 + 0.5, 0);
	
	float l0_dy = v2Raster[1] - v1Raster[1];
	float l1_dy = v0Raster[1] - v2Raster[1];
	float l2_dy = v1Raster[1] - v0Raster[1];

	float l0_dx = v2Raster[0] - v1Raster[0];
	float l1_dx = v0Raster[0] - v2Raster[0];
	float l2_dx = v1Raster[0] - v0Raster[0];
	
	float w0_x0y0 = (pp_x0y0[0]-v1Raster[0])*l0_dy - (pp_x0y0[1]-v1Raster[1])*l0_dx;
	float w1_x0y0 = (pp_x0y0[0]-v2Raster[0])*l1_dy - (pp_x0y0[1]-v2Raster[1])*l1_dx;
	float w2_x0y0 = (pp_x0y0[0]-v0Raster[0])*l2_dy - (pp_x0y0[1]-v0Raster[1])*l2_dx;

	float w0_x0 = w0_x0y0;
	float w1_x0 = w1_x0y0;
	float w2_x0 = w2_x0y0;
	
        for (uint32_t yy = 0; yy <= ylen; yy++) {
	  uint32_t y = y0 + yy;
	  float w0_ = w0_x0, w1_ = w1_x0, w2_ = w2_x0;
	  
	  for (uint32_t xx = 0; xx <= xlen; xx++) {
	    uint32_t x = x0 + xx;
	    float w0 = w0_, w1 = w1_, w2 = w2_;
	    incr_count(&n_tests,1UL);
	    
	    if ((w0 >= 0 && w1 >= 0 && w2 >= 0)) {
	      incr_count(&n_passed_area,1UL);
	      
	      w0 *= recip_area;
	      w1 *= recip_area;
	      w2 *= recip_area;
	      //3 mults
	      
	      float z = 1.0f / (v0Raster.z * w0 + v1Raster.z * w1 + v2Raster.z * w2);
	      //1 recip, 6 mults, 2 adds
	      
	      bool failedZ = true;
	      if(depthBuffer[y * imageWidth + x] > z) {
		depthBuffer[y * imageWidth + x] = z;
		failedZ = false;
	      }
	      
	      if(not(failedZ)) {
		incr_count(&n_passed_z, 1UL);
		Vec2f st;
		float w0z = w0*z, w1z = w1*z, w2z = w2*z;
		//1 compare, 1 recip, 9 mults, 2 adds
		
		st.x = (st0.x * w0z) + (st1.x * w1z) + (st2.x * w2z);
		//1 compare, 1 recip, 12 mults, 4 adds
		
		st.y = (st0.y * w0z) + (st1.y * w1z) + (st2.y * w2z);
		//1 compare, 1 recip, 15 mults, 6 adds

		
		// [comment]
		// If you need to compute the actual position of the shaded
		// point in camera space. Proceed like with the other vertex attribute.
		// Divide the point coordinates by the vertex z-coordinate then
		// interpolate using barycentric coordinates and finally multiply
		// by sample depth.
		// [/comment]
                        
		float px = v0Cam_xz * w0 + v1Cam_xz * w1 + v2Cam_xz * w2;
		//1 compare, 1 recip, 18 mults, 8 adds
		
		float py = v0Cam_yz * w0 + v1Cam_yz * w1 + v2Cam_yz * w2;
		//1 compare, 1 recip, 21 mults, 10 adds
		
		// [comment]
		// Compute the face normal which is used for a simple facing ratio.
		// Keep in mind that we are doing all calculation in camera space.
		// Thus the view direction can be computed as the point on the object
		// in camera space minus Vec3f(0), the position of the camera in camera
		// space.
		// [/comment]
		Vec3f viewDirection(-px * z, -py * z, z); // pt is in camera space;
		//1 compare, 1 recip, 20 mults, 8 adds
		
		float t = viewDirection[0]*viewDirection[0] +
		  viewDirection[1]*viewDirection[1] +
		  viewDirection[2]*viewDirection[2];
		//1 compare, 1 recip, 23 mults, 10 adds

	      
		float t_recip = 1.0/std::sqrt(t);
		viewDirection[0] *= t_recip;
		viewDirection[1] *= t_recip;
		viewDirection[2] *= t_recip;
		//1 compare, 1 recip, 1 recip sqrt, 26 mults, 10 adds		

		float nDotView = n[0]*viewDirection[0] + n[1]*viewDirection[1] + n[2]*viewDirection[2];
		nDotView = std::max(0.0f, nDotView);
		//1 compare, 1 recip, 1 recip sqrt, 29 mults, 12 adds				
                        
		// [comment]
		// The final color is the reuslt of the faction ration multiplied by the
		// checkerboard pattern.
		// [/comment]
		const int M = 10;
		float checker = (fmod(st.x * M, 1.0) > 0.5) ^ (fmod(st.y * M, 1.0) < 0.5);
		float c = 0.3 * (1 - checker) + 0.7 * checker;
		nDotView *= c;

		pxls[y * imageWidth + x].r = nDotView * 255;
		pxls[y * imageWidth + x].g = nDotView * 255;
		pxls[y * imageWidth + x].b = nDotView * 255;		
	      }
	    }
	    
	    w0_ += l0_dy;
	    w1_ += l1_dy;
	    w2_ += l2_dy;
	  }

	  w0_x0 -= l0_dx;
	  w1_x0 -= l1_dx;
	  w2_x0 -= l2_dx;
        }
    }
    SDL_UnlockSurface(sdlscr);
    SDL_UpdateWindowSurface(sdlwin);
    
    std::cout << "n_tests        = " << n_tests << "\n";
    std::cout << "n_passed_area  = " << n_passed_area << "\n";
    std::cout << "n_passed_z     = " << n_passed_z << "\n";
    std::cout << static_cast<double>(n_passed_area)/ntris << " pixels per triange\n";
    auto t_end = std::chrono::high_resolution_clock::now();
    auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
    std::cerr << "Wall passed time:  " << passedTime << " ms" << std::endl;
    }

    delete [] depthBuffer;


    SDL_DestroyWindow(sdlwin);
    SDL_Quit();
    
    return 0;
}
