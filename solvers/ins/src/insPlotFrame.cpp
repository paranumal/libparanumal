#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
#include <png.h>
#include "ins.hpp"   // assumes your class definition is in ins.h

// Helper structures
struct Point2D {
  dfloat x, y;
};

struct Triangle {
  int v0, v1, v2;
};

struct Color {
  unsigned char r, g, b;
};

// Simple linear colormap: maps value in [minVal, maxVal] to a color from blue (min) to red (max)
static Color colormapRB(dfloat value, dfloat minVal, dfloat maxVal) {
  dfloat norm = (value - minVal) / (maxVal - minVal);
  norm = std::max((dfloat)0.0, std::min((dfloat)1.0, norm));
  unsigned char r = static_cast<unsigned char>(norm * 255);
  unsigned char g = 0;
  unsigned char b = static_cast<unsigned char>((1.0 - norm) * 255);
  return {r, g, b};
}

static Color colormap(dfloat value, dfloat minVal, dfloat maxVal) {
  // Normalize the value to [0, 1]
  dfloat norm = (value - minVal) / (maxVal - minVal);
  norm = std::max((dfloat)0.0, std::min((dfloat)1.0, norm));

  dfloat r, g, b;
  if (norm < 0.35) {
    // Lower region: blue to cyan
    r = 0;
    g = norm / 0.35;
    b = 1;
  } else if (norm < 0.66) {
    // Middle region: cyan to yellow
    r = (norm - 0.35) / (0.66 - 0.35);
    g = 1;
    b = 1 - (norm - 0.35) / (0.66 - 0.35);
  } else {
    // Upper region: yellow to red
    r = 1;
    g = 1 - (norm - 0.66) / (1 - 0.66);
    b = 0;
  }

  unsigned char R = static_cast<unsigned char>(r * 255);
  unsigned char G = static_cast<unsigned char>(g * 255);
  unsigned char B = static_cast<unsigned char>(b * 255);
  return {R, G, B};
}



// Write PNG using libpng
static bool WritePNG(const char* filename, int width, int height, unsigned char* image) {
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Error: Cannot open %s for writing.\n", filename);
    return false;
  }
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
  if (!png_ptr) { fclose(fp); return false; }
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr, nullptr);
    fclose(fp);
    return false;
  }
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    return false;
  }
  png_init_io(png_ptr, fp);
  png_set_IHDR(png_ptr, info_ptr, width, height,
               8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);
  for (int y = 0; y < height; y++) {
    png_bytep row = image + y * width * 3;
    png_write_row(png_ptr, row);
  }
  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fp);
  return true;
}

// New member function to save a PNG image instead of VTU output.
void ins_t::PlotFrame(memory<dfloat>& Q, std::string fileName, int plotNfields) {
  //--- Step 1. Interpolate plot node coordinates and Q0 for each element ---
  // (We assume here you wish to visualize the first field: Q0)
  //  std::vector<Point2D> globalPoints;

#if 1
  
  int totalPoints = mesh.Nelements * mesh.plotNp;

  memory<Point2D> globalPoints(totalPoints);
  memory<dfloat> globalScalars(totalPoints);

  // Scratch memory similar to your VTU writer.
  dlong Nscratch = std::max(mesh.Np, mesh.plotNp);
  memory<dfloat> scratch(2 * Nscratch);
  memory<dfloat> Ix(mesh.plotNp);
  memory<dfloat> Iy(mesh.plotNp);

  memory<dfloat> Iq(mesh.plotNp);

  for(dlong e = 0; e < mesh.Nelements; ++e) {
    // Compute plot node coordinates.
    mesh.PlotInterp(mesh.x + e * mesh.Np, Ix, scratch);
    mesh.PlotInterp(mesh.y + e * mesh.Np, Iy, scratch);
    // Interpolate Q field (assuming field index 0: Q0)
    //    mesh.PlotInterp(Q + e * mesh.Np * Nfields + fld * mesh.Np, Iq[fld], scratch);
    mesh.PlotInterp(Q + e * mesh.Np * plotNfields, Iq, scratch);
    
    for (int n = 0; n < mesh.plotNp; ++n) {
      globalPoints[e*mesh.plotNp + n] = {static_cast<dfloat>(Ix[n]),
	static_cast<dfloat>(Iy[n])};

#if 0
      dfloat magq = 0;

      for(int fld=0;fld<plotNfields;++fld){
	magq += Iq[fld][n]*Iq[fld][n];

      magq = sqrt(magq);
#endif 
      globalScalars[n+e*mesh.plotNp] = Iq[n];
    }
  }
  
  //--- Step 2. Build global connectivity for triangles ---
  // For 2D, your VTU writer writes triangles (cell type 5). For each element, there are mesh.plotNelements triangles.
  int totalTriangles = mesh.Nelements * mesh.plotNelements;
  memory<Triangle> globalTriangles(totalTriangles);
  
  for(dlong e = 0; e < mesh.Nelements; ++e) {
    for (int n = 0; n < mesh.plotNelements; ++n) {
      Triangle tri;
      // Each triangle uses vertices from the interpolated plot nodes.
      tri.v0 = e * mesh.plotNp + mesh.plotEToV[n * mesh.plotNverts + 0];
      tri.v1 = e * mesh.plotNp + mesh.plotEToV[n * mesh.plotNverts + 1];
      tri.v2 = e * mesh.plotNp + mesh.plotEToV[n * mesh.plotNverts + 2];
      globalTriangles[e*mesh.plotNelements+n] = tri;
    }
  }
  
  //--- Step 3. Determine bounding box of the mesh ---
  dfloat minX = globalPoints[0].x, maxX = globalPoints[0].x;
  dfloat minY = globalPoints[0].y, maxY = globalPoints[0].y;
  for (dlong i = 1; i < totalPoints; ++i) {
    minX = std::min(minX, globalPoints[i].x);
    maxX = std::max(maxX, globalPoints[i].x);
    minY = std::min(minY, globalPoints[i].y);
    maxY = std::max(maxY, globalPoints[i].y);
  }
  // Optionally add a little padding
  dfloat padding = 0.1;
  minX -= padding; maxX += padding;
  minY -= padding; maxY += padding;
  
  //--- Step 4. Set up image buffer with aspect ratio matching the mesh ---
  // Choose target resolution for the longer side.
  int targetResolution = 4800;
  dfloat meshWidth = maxX - minX;
  dfloat meshHeight = maxY - minY;
  int width, height;
  if (meshWidth >= meshHeight) {
    width = targetResolution;
    height = static_cast<int>(targetResolution * (meshHeight / meshWidth));
  } else {
    height = targetResolution;
    width = static_cast<int>(targetResolution * (meshWidth / meshHeight));
  }

  memory<unsigned char> image(width * height * 3, 255); // white background


  for(int n=0;n<width*height;++n){
    image[n*3+0] = 63;
    image[n*3+1] = 127;
    image[n*3+2] = 191;
  }
  
  // Determine scalar range for colormap.
  dfloat minVal = 0, maxVal = 1;

#if 1
  if(settings.compareSetting("VIZ COLORMAP RANGE", "TRUE")){
    settings.getSetting("VIZ COLORMAP RANGE MIN", minVal);
    settings.getSetting("VIZ COLORMAP RANGE MAX", maxVal);
  }else
    {
    minVal = globalScalars[0];
    maxVal = globalScalars[0];
    for (dlong i = 1; i < totalPoints; ++i) {
      minVal = std::min(minVal, globalScalars[i]);
      maxVal = std::max(maxVal, globalScalars[i]);
    }
  }
#endif
  
  std::cout << "min,maxVal = " << minVal << ", " << maxVal << std::endl;
  
  //--- Step 5. Compute scale factor.
  // With the image dimensions computed to match the mesh aspect ratio, the scale is:
  dfloat scale = static_cast<dfloat>(width) / meshWidth; // equals height/meshHeight
  
  //--- Step 6. Rasterize triangles ---
  for (dlong i = 0; i < totalTriangles; ++i) {
    Triangle tri = globalTriangles[i];
    Point2D aa = globalPoints[tri.v0];
    Point2D bb = globalPoints[tri.v1];
    Point2D cc = globalPoints[tri.v2];
    dfloat aVal = globalScalars[tri.v0];
    dfloat bVal = globalScalars[tri.v1];
    dfloat cVal = globalScalars[tri.v2];

    // Determine triangle bounding box in mesh coordinates.
    dfloat triMinX = std::min({aa.x, bb.x, cc.x});
    dfloat triMaxX = std::max({aa.x, bb.x, cc.x});
    dfloat triMinY = std::min({aa.y, bb.y, cc.y});
    dfloat triMaxY = std::max({aa.y, bb.y, cc.y});
    
    // Map mesh bounding box to pixel coordinates.
    int pixMinX = std::max((int)0, static_cast<int>((triMinX - minX) * scale));
    int pixMaxX = std::min(width - 1, static_cast<int>((triMaxX - minX) * scale));
    int pixMinY = std::max((int)0, static_cast<int>((triMinY - minY) * scale));
    int pixMaxY = std::min(height - 1, static_cast<int>((triMaxY - minY) * scale));
    
    // Loop over pixels in the triangle's bounding box.
    for (int py = pixMinY; py <= pixMaxY; ++py) {
      for (int px = pixMinX; px <= pixMaxX; ++px) {
        // Convert pixel center back to mesh coordinates.
        dfloat meshX = minX + (px + 0.5) / scale;
        dfloat meshY = minY + (py + 0.5) / scale;
        
        // Compute barycentric coordinates.
        dfloat denom = ((bb.y - cc.y) * (aa.x - cc.x) + (cc.x - bb.x) * (aa.y - cc.y));
        if (std::abs(denom) < 1e-12) continue; // skip degenerate triangle
        dfloat alpha = ((bb.y - cc.y) * (meshX - cc.x) + (cc.x - bb.x) * (meshY - cc.y)) / denom;
        dfloat beta  = ((cc.y - aa.y) * (meshX - cc.x) + (aa.x - cc.x) * (meshY - cc.y)) / denom;
        dfloat gamma = 1.0 - alpha - beta;
        
        if (alpha >= 0 && beta >= 0 && gamma >= 0) {
          // Interpolate scalar value.
          dfloat value = alpha * aVal + beta * bVal + gamma * cVal;
          // Map value to color.
          Color col = colormap(value, minVal, maxVal);
	  
          // Write pixel (flip y so image origin is at top-left).
	  if( py<=height-1 && px<=width-1){
	    int idx = ((height - 1 - py) * width + px) * 3;
	    image[idx]   = col.r;
	    image[idx+1] = col.g;
	    image[idx+2] = col.b;
	  }
        }
      }
    }
  }
  
  //--- Step 7. Write image to PNG ---
#if 1
  if (!WritePNG(fileName.c_str(), width, height, image.ptr())) {
    fprintf(stderr, "Error writing PNG file %s\n", fileName.c_str());
  }
#endif

#endif
  
}
