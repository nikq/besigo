#ifndef __FRAMEBUFFER_H
#define __FRAMEBUFFER_H

#include <stdio.h>
#include <stdlib.h>
#include "vectormath.h" // for vectormath
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#undef STB_IMAGE_WRITE_IMPLEMENTATION

namespace FRAMEBUFFER {

  typedef VECTORMATH::Vector3<float> Color;      // R,G,B

  class FrameBuffer {
  private:
  public:
    int width_;
    int height_;
    int sample_;
    int pass_;
    Color *film_;
    Color *raw_;
    float *weight_;
    float *dx_;
    float *dy_;

    virtual ~FrameBuffer(){ clear(); }
    FrameBuffer() : width_(0), height_(0), sample_(0), pass_(0), film_(NULL), raw_(NULL), dx_(NULL), dy_(NULL){;}

    void setup( int w, int h, int s = 1 ) {
      int i;
      float a;
      width_ = w;
      height_= h;
      sample_= s;
      film_ = (Color*) malloc( sizeof( Color ) * w * h );
      raw_  = (Color*) malloc( sizeof( Color ) * w * h );
      weight_=(float*)malloc(sizeof(float)*w*h);
      dx_   = new float[ s ];
      dy_   = new float[ s ];

      for( i = 0; i<w*h; i++ ){
        film_[i].set( 0., 0., 0. );
        weight_[i]=0.;
      }

      a = 0.;
      for(i = 0;i<sample_;i++){
        dx_[i] = (float) i / sample_;
        dy_[i] = a - floor( a );
        a      = a + (1. + sqrt(5.))/2.;
      }
    }

    void clear(void) {
      width_ = height_ = sample_ = pass_ = 0;
      if( film_ ){ free( film_ ); film_ = NULL; }
      if( raw_  ){ free( raw_ );  raw_ = NULL;  }
      if( dx_   ){ delete [] dx_; dx_ = NULL; }
      if( dy_   ){ delete [] dy_; dy_ = NULL; }
    }

    void saturate( void ) {
      int i,j;
      float lo[3] = {FLT_MAX,FLT_MAX,FLT_MAX};
      float hi[3] = {FLT_MIN,FLT_MIN,FLT_MIN};
      for(i=0;i<width_*height_;i++){
        for(j=0;j<3;j++){
          float v = film_[i].get(j);
          if( v > hi[j] ) hi[j] = v;
          if( v < lo[j] ) lo[j] = v;
        }
      }
      for(i=0;i<width_*height_;i++){
        for(j=0;j<3;j++){
          raw_[i].set(j,(film_[i].get(j)-lo[j]) / (hi[j]-lo[j]));
        }
      }
    }
    void normalize( void ) {
      int i,j;
      for(i=0;i<width_*height_;i++){
        float w = 1.;
        if( weight_[i] > 0. )
          w = weight_[i];
        for(j=0;j<3;j++){
          raw_[i].set(j,film_[i].get(j) / w);
        }
      }
    }
    void add( int x, int y, float r, float g, float b ){
      weight_[ x + y * width_ ] += 1.f;
      film_[ x + y * width_ ] = film_[ x + y * width_ ] + Color( r,g,b );
    }
    void set( int x, int y, float r, float g, float b ){
      weight_[ x + y * width_ ] = 1.f;
      film_[ x + y * width_ ] = Color( r,g,b );
    }

    int toInt(float x,float g){
      // tone mapping
      return int(pow(1-exp(-x),1/g)*255+.5);
      //return int( x*255. );
    }
    
    Color linearToneMapping( Color color, float exposure ) {
      color = color * exposure;
      return color.pow( 1.f / 2.2f );
    }
    Color Uncharted2ToneMapping( Color color, float exposure ) {
      const float A = 0.15f;
      const float B = 0.50f;
      const float C = 0.10f;
      const float D = 0.20f;
      const float E = 0.02f;
      const float F = 0.30f;
      const float W = 11.2f;
      const float gamma = 2.2f;

      float r = color.get(0) * exposure;
      float g = color.get(1) * exposure;
      float b = color.get(2) * exposure;

      r = ((r*(A*r+C*B)+D*E)/(r*(A*r+B)+D*F))-E/F;
      g = ((g*(A*g+C*B)+D*E)/(g*(A*g+B)+D*F))-E/F;
      b = ((b*(A*b+C*B)+D*E)/(b*(A*b+B)+D*F))-E/F;
      float white = ((W * (A * W + C * B) + D * E) / (W * (A * W + B) + D * F)) - E / F;
      return Color(r/white,g/white,b/white).pow(1.f/2.2f);
    }

    
    void save_ldr( const char *filename, float exposureAdjust = 1.f ) {
      unsigned char *rgb;
      int pixel;
      int x,y,i,j;

      float avg = 0.;
      for( int i=0;i<height_*width_;i++)
        avg += (raw_[i].get(0)+raw_[i].get(1)+raw_[i].get(2)) / 3.f;
      
      float scale = 0.5f / (avg / (float)(height_*width_)) * exposureAdjust;
      assert( avg == avg );
      printf("avg %f, scale %f\n",avg,scale);
      
      rgb = (unsigned char*)malloc(width_ * height_ * 3);
      for(y = 0; y < height_; y ++ ){
        for(x = 0;x < width_; x ++ ){
          i = (y*width_)+x;
          j = (((height_-y-1)*width_)+x)*3;
          Color after = Uncharted2ToneMapping( raw_[i], scale );
          //Color after = linearToneMapping( raw_[i], scale );
          pixel      = (int) (after.get(0) * 255.f);
          rgb[j ] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);
          pixel      = (int) (after.get(1) * 255.f);
          rgb[j+1] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);
          pixel      = (int) (after.get(2) * 255.f);
          rgb[j+2] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);
        }
      }
      //bmp_write( filename, width_, height_, rgb );
      stbi_write_png( filename, width_, height_, 3, rgb, width_*3 );
      free(rgb);
    }
    
    void save_bmp( const char *filename, float scale, float gamma ) {
      unsigned char *rgb;
      int pixel;
      int x,y,i,j;
      rgb = (unsigned char*)malloc(width_ * height_ * 3);
      for(y = 0; y < height_; y ++ ){
        for(x = 0;x < width_; x ++ ){
          i = (y*width_)+x;
          j = (((height_-y-1)*width_)+x)*3;
          pixel      = (int) toInt( raw_[i].get(0) * scale, gamma );
          rgb[j ] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);

          pixel      = (int) toInt( raw_[i].get(1) * scale, gamma );
          rgb[j+1] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);

          pixel      = (int) toInt( raw_[i].get(2) * scale, gamma );
          rgb[j+2] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);
        }
      }
      bmp_write( filename, width_, height_, rgb );
      free(rgb);
    }
    void save_png( const char *filename, float scale, float gamma ) {
      unsigned char *rgb;
      int pixel;
      int x,y,i,j;
      rgb = (unsigned char*)malloc(width_ * height_ * 3);
      for(y = 0; y < height_; y ++ ){
        for(x = 0;x < width_; x ++ ){
          i = (y*width_)+x;
          j = (((height_-y-1)*width_)+x)*3;
          pixel      = (int) toInt( raw_[i].get(0) * scale, gamma );
          rgb[j ] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);

          pixel      = (int) toInt( raw_[i].get(1) * scale, gamma );
          rgb[j+1] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);

          pixel      = (int) toInt( raw_[i].get(2) * scale, gamma );
          rgb[j+2] = (pixel>255) ? 255: ((pixel<0) ? 0:pixel);
        }
      }
      stbi_write_png( filename, width_, height_, 3, rgb, 3 * width_ );
      free(rgb);
    }
    void save_hdr( const char *filename, float scale, float gamma ) {
      float *rgb;
      int x,y,i,j;
      float avg = 0.;
      for( int i=0;i<height_*width_;i++)
        avg += (raw_[i].get(0)+raw_[i].get(1)+raw_[i].get(2)) / 3.f;
      
      float s = 0.5f / (avg / (float)(height_*width_));
      rgb = (float*)malloc(width_ * height_ * 3 * sizeof(float));
      for(y = 0; y < height_; y ++ ){
        for(x = 0;x < width_; x ++ ){
          i = (y*width_)+x;
          j = i*3;
          rgb[j+0] = raw_[i].get(0) * s;
          rgb[j+1] = raw_[i].get(1) * s;
          rgb[j+2] = raw_[i].get(2) * s;
        }
      }
      hdr_write( filename, width_, height_, rgb );
      free(rgb);
    }


    typedef struct {
      unsigned char  id[2];
      unsigned long  filesize;
      unsigned short reserve1;
      unsigned short reserve2;
      unsigned long  offset; // 26:12,54:40,122:108
      unsigned long  headsize; // 12,40,108
      unsigned long  width;
      unsigned long  height;
      unsigned short plane;
      unsigned short bpp;
      unsigned long  method;
      unsigned long  datasize;
      unsigned long  x_dpm;
      unsigned long  y_dpm;
      unsigned long  palnum;
      unsigned long  imp_pal;
    }BMP_header;

    inline float maxf( float a, float b ){ return (a>b) ? a:b ; }

    int hdr_write( const char *name,int width,int height,float *rgb){
      FILE *fp;
      fp = fopen( name, "wb" );
      if(!fp)
        return -1;
      fprintf( fp, "#?RADIANCE\n");
      fprintf( fp, "# Made with RLR\n");
      fprintf( fp, "FORMAT=32-bit_rle_rgbe\n");
      fprintf( fp, "EXPOSURE=          1.0000000000000\n");
      fprintf( fp, "\n");
      fprintf( fp, "-Y %d +X %d\n",height,width);

      //unsigned char magic4[4] = {0x02, 0x02, (width>>8)&0xFF, width&0xFF};
      unsigned char line[ 0x7FFF*4 ];

      //for(int y=0;y<height;y++){
      for(int y=height-1;y>=0;y--){

        //fwrite( magic4, 1, 4, fp );
        // RGBE
        for(int x=0;x<width;x++){
          float r = rgb[ (y*width + x)*3 + 0 ];
          float g = rgb[ (y*width + x)*3 + 1 ];
          float b = rgb[ (y*width + x)*3 + 2 ];
          float m = maxf(maxf(r,g),b);

          int e;
          int iR, iG, iB, iE;

          if(m <= 1e-32) {
            iR = iG = iB = iE = 0;
          } else {
            m = frexp(m, &e) * 255.9999 / m;

            iR = (int)(r * m);
            iG = (int)(g * m);
            iB = (int)(b * m);
            iE = (int)(e + 128);
          }
          line[ x * 4 + 0 ] = iR;
          line[ x * 4 + 1 ] = iG;
          line[ x * 4 + 2 ] = iB;
          line[ x * 4 + 3 ] = iE;
        }

        fwrite( line, 1, width*4, fp );
      }

      fclose(fp);
      return 0;
    }

    int bmp_write( const char *name,int width,int height,unsigned char *rgb) {
      BMP_header bh;
      FILE *fp;
      int x,y,k;
      int r,g,b;
      int mod;

      bh.id[0] = 'B';
      bh.id[1] = 'M';
      bh.filesize = 54 + (width*height*3);
      bh.offset   = 54;
      bh.headsize = 40;
      bh.width    = width;
      bh.height   = height;
      bh.plane    = 1;
      bh.bpp      = 24;
      bh.method   = 0;
      bh.datasize = 0;
      bh.x_dpm    = 3779;
      bh.y_dpm    = 3779;
      bh.palnum   = 0;
      bh.imp_pal  = 0;

      fp = fopen(name,"wb");
      if(!fp) return -1;
      fwrite( bh.id, 1,2 ,fp);
      fwrite( &(bh.filesize), 4,1,fp);
      fwrite( &(bh.reserve1), 2,1,fp);
      fwrite( &(bh.reserve2), 2,1,fp);
      fwrite( &(bh.offset  ), 4,1,fp);
      fwrite( &(bh.headsize), 4,1,fp);
      fwrite( &(bh.width   ), 4,1,fp);
      fwrite( &(bh.height  ), 4,1,fp);
      fwrite( &(bh.plane   ), 2,1,fp);
      fwrite( &(bh.bpp     ), 2,1,fp);
      fwrite( &(bh.method  ), 4,1,fp);
      fwrite( &(bh.datasize), 4,1,fp);
      fwrite( &(bh.x_dpm   ), 4,1,fp);
      fwrite( &(bh.y_dpm   ), 4,1,fp);
      fwrite( &(bh.palnum  ), 4,1,fp);
      fwrite( &(bh.imp_pal ), 4,1,fp);

      mod = (width*3)%4;
      if(mod)
        mod = 4-mod;
      //printf("mod:%d\n",mod);
      for(y=0;y<height;y++){
        for(x=0;x<width;x++){
          k = ((height-y-1)*width + x)*3;

          r = rgb[k  ]; // R
          g = rgb[k+1]; // G
          b = rgb[k+2]; // B
          fputc(b,fp);
          fputc(g,fp);
          fputc(r,fp);
        }
        if(mod>0){
          for(x=0;x<mod;x++){
            fputc(0,fp);
          }
        }
      }
      fclose(fp);
      return 0;
    }

  };

} // FRAMEBUFFER

#endif
