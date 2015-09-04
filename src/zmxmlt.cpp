
// GAJAI!!!!!!!!
// BESIGO!!!!!!

#include <stdio.h>
#include <algorithm>
#include <vector>
#include <stack>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "vectormath.h"
#include "framebuffer.h"
#include "mtseq.h"
#include "color3.h"
#include "xor128.h"
#include "model.h"
#include "kelemenmlt.h"
#include "lensnuma.h"
#include "stb_image.c"

#include "include/embree2/rtcore.h" // EMBREE.
#include "include/embree2/rtcore_ray.h" // EMBREE.
#include <xmmintrin.h>
#include <pmmintrin.h>


// #define CHECK // チェックモード。低解像度、7スレッド
// #define DEBUG // デバッグモード。低解像度、1スレッド



// 上のパストレで生成したパスを保存しておく
struct PathSample {
  int x, y;

  double wavelength;
  int    wavelength_index;
  double F;
  double weight;
  
  PathSample(const int x_ = 0, const int y_ = 0, const double &wavelength_ = 0., const double F_ = 0., const double weight_ = 1.0) :
  x(x_), y(y_), wavelength(wavelength_), F(F_), weight(weight_) {}
  
  double luminance( void ){
    // return SPECTRUM::getXYZ( wavelength, F ).b_;
    COLOR3::Color3 rgb = SPECTRUM::getXYZ( wavelength, F ).get( COLOR3::RGB );
    return std::max( std::max( rgb.a_, rgb.b_ ), rgb.c_ );
  }

  FRAMEBUFFER::Color color( void ) {
    COLOR3::Color3 rgb = SPECTRUM::getXYZ( wavelength, F ).get( COLOR3::RGB );
    return FRAMEBUFFER::Color( rgb.a_, rgb.b_, rgb.c_ );
  }
};


class Texture {
public:
  int w_;
  int h_;
  std::vector<COLOR3::Color3> rgb_;
  std::vector<std::vector<float>> pdf_;
  std::vector<std::vector<float>> cdf_;
  
  Texture():w_(0),h_(0){;}
  bool isValid(){ return !w_ && !h_; }
  void load( const char *fn ){
    int n;
    float * data = stbi_loadf(fn, &w_, &h_, &n, 0);
    if(!data){
      w_ = h_ = 1;
      rgb_.resize(1);
      rgb_[0].set(COLOR3::RGB,0,0,0);
    }else{
      rgb_.resize( w_ * h_ * 3 );
      for(int i=0;i<w_*h_;i++)
        rgb_[i] = COLOR3::Color3( COLOR3::RGB, data[i*3],data[i*3+1],data[i*3+2] ).get( COLOR3::XYZ );
      free( data );
    }
  }

  const COLOR3::Color3 & pixelNN( int x, int y ) const {
    if( x < 0 ) x = 0;
    if( x >=w_) x = w_-1;
    if( y < 0 ) y = 0;
    if( y >=h_) y = h_-1;
    return rgb_[x+y*w_];
  }
  const COLOR3::Color3 pixel( float x, float y ) const {
    int xl = (int)(x-0.5);
    int yl = (int)(y-0.5);
    int xh = xl + 1;
    int yh = yl + 1;
    float xm = x - xl;
    float ym = y - yl;
    COLOR3::Color3 ll = pixelNN( xl, yl );
    COLOR3::Color3 lh = pixelNN( xl, yh );
    COLOR3::Color3 hl = pixelNN( xh, yl );
    COLOR3::Color3 hh = pixelNN( xh, yh );
    COLOR3::Color3 l(
      ll.space_,
      ll.a_ * (1.f - ym) + lh.a_ * ym,
      ll.b_ * (1.f - ym) + lh.b_ * ym,
      ll.c_ * (1.f - ym) + lh.c_ * ym);
    COLOR3::Color3 h(
      hl.space_,
      hl.a_ * (1.f - ym) + hh.a_ * ym,
      hl.b_ * (1.f - ym) + hh.b_ * ym,
      hl.c_ * (1.f - ym) + hh.c_ * ym);
    COLOR3::Color3 p(
      l.space_,
      l.a_ * (1.f - xm) + h.a_ * xm,
      l.b_ * (1.f - xm) + h.b_ * xm,
      l.c_ * (1.f - xm) + h.c_ * xm);
    return p;
  }

  const COLOR3::Color3 &sphere( VECTORMATH::Vector3<double>& v ){
    double t = (atan2( v.x , v.z ) / M_PI) + 1.;
    double p = acos( v.y ) / M_PI;
    return pixelNN( w_ * (1. - t / 2.), h_ * p );
  }

  const COLOR3::Color3 &tex2D( double u, double v ){
    return pixelNN( w_ * u, h_ * v );
  }
};

class Scene {
public:
  
  MODEL::MQO model_;        // mesh model. contains materials, verticies, indicies.
  Texture  background_;
  RTCScene scene_;
  
  void init( const char *model, const char *bgfn ){
    scene_ = rtcNewScene( RTC_SCENE_STATIC | RTC_SCENE_HIGH_QUALITY, RTC_INTERSECT1 );
    model_.load( model, scene_ );
    rtcCommit( scene_ );
    background_.load( bgfn );
  }

  // 背景に逃げたレイの扱い.
  double radiance_bg( KelemenMLT& mlt, VECTORMATH::Vector& pos, VECTORMATH::Vector& dir, double lambda ){
    const COLOR3::Color3 &bg( background_.sphere( dir ) );
    return SPECTRUM::getLuminance( bg, lambda );
  }
  
  // 直接照明サンプリング.
  bool isNaN( volatile double b ){
    return !(b == b);
  }
  double solidAngle( VECTORMATH::Vector& p, VECTORMATH::Vector& v1, VECTORMATH::Vector& v2, VECTORMATH::Vector& v3 ) {
    // http://be.nucl.ap.titech.ac.jp/~koba/cgi-bin/moin.cgi/%E7%AB%8B%E4%BD%93%E8%A7%92
    VECTORMATH::Vector e1 = v1 - p; e1.normalize();
    VECTORMATH::Vector e2 = v2 - p; e2.normalize();
    VECTORMATH::Vector e3 = v3 - p; e3.normalize();
    double cosA = e1.dot( e2 ); double c2a = cosA*cosA;
    double cosB = e2.dot( e3 ); double c2b = cosB*cosB;
    double cosC = e3.dot( e1 ); double c2c = cosC*cosC;
    double subx = 1. - c2a - c2b - c2c + 2. * cosA * cosB * cosC;
    if( subx < 0. )
      return 0.;
    double x = sqrt( subx ) * (1. + cosA + cosB + cosC);
    double y = (cosA + cosB + cosC + c2a + c2b + c2c + cosA*cosB + cosA*cosC + cosB*cosC - cosA*cosB*cosC );
    // if( fabs(x)<1e-6 && fabs(y)<1e-6 ) return 0.;
    double r = atan2( x, y );
    // assert( !isNaN(r));
    return r;
  }
  double direct_radiance( KelemenMLT& mlt, VECTORMATH::Vector& pos, VECTORMATH::Vector& normal, double lambda ){
    VECTORMATH::Vector v1,v2,v3;
    VECTORMATH::Vector emitterVert = model_.emitterDir( mlt,v1,v2,v3 );
    VECTORMATH::Vector dir = emitterVert - pos;
    dir.normalize();
    if( normal.dot( dir ) < 0. )
      return 0.;
    
    double radiance = 0.;
    RTCRay ray;
    {
      ray.org[0] = pos.x;
      ray.org[1] = pos.y;
      ray.org[2] = pos.z;
      ray.dir[0] = dir.x;
      ray.dir[1] = dir.y;
      ray.dir[2] = dir.z;
      ray.tnear  = 1e-6; // 近すぎるのは交差しない.
      ray.tfar   = 1e8;
      ray.geomID = RTC_INVALID_GEOMETRY_ID;
      ray.primID = RTC_INVALID_GEOMETRY_ID;
      ray.instID = RTC_INVALID_GEOMETRY_ID;
      ray.mask   = 0xFFFFFFFF;
      ray.time   = 0.f;
    }
    rtcIntersect( scene_, ray );
    if( ray.geomID == RTC_INVALID_GEOMETRY_ID )
      return radiance_bg( mlt, pos, dir, lambda );

    VECTORMATH::Vector light_normal = model_.normal( ray );
    double dot0 = normal.dot( dir );
    double dot1 = light_normal.dot( dir * -1. );
    if( dot1 < 0. )
      return 0.;
    
    double G = dot0 * dot1 / ray.tfar;
    double r = solidAngle( pos, v1, v2, v3 );
    return model_.material(ray)->emitter(lambda) * G * (r / M_PI) / M_PI;
  }
  
  VECTORMATH::Vector reflect( const VECTORMATH::Vector& I, const VECTORMATH::Vector& N ){
    return I - N * 2. * N.dot( I );
  }
  bool refract( const VECTORMATH::Vector& I, VECTORMATH::Vector& N, double eta, VECTORMATH::Vector& result ){
    double cosi  = - I.dot( N ); // dot(-i, n);
    double cost2 = 1.0 - eta * eta * (1.0 - cosi * cosi);
    if( cost2 < 0. ) return false; // 全反射.
    result = I * eta + (N * (eta * cosi - sqrtf( cost2 )));
    return true;
  }
  
  // pathtrace!
  double radiance( int depth, KelemenMLT& mlt, NUMA::Vector& pos, NUMA::Vector& dir, double lambda, double *dist = NULL ){
    
    if( depth > 32 )
      return 0.;
    
    // do intersect.
    dir.normalize();
    RTCRay ray;
    {
      ray.org[0] = pos.x;
      ray.org[1] = pos.y;
      ray.org[2] = pos.z;
      ray.dir[0] = dir.x;
      ray.dir[1] = dir.y;
      ray.dir[2] = dir.z;
      ray.tnear  = 1e-6; // 近すぎるのは交差しない.
      ray.tfar   = 1e8;
      ray.geomID = RTC_INVALID_GEOMETRY_ID;
      ray.primID = RTC_INVALID_GEOMETRY_ID;
      ray.instID = RTC_INVALID_GEOMETRY_ID;
      ray.mask   = 0xFFFFFFFF;
      ray.time   = 0.f;
    }
    rtcIntersect( scene_, ray );
    
    if( ray.geomID == RTC_INVALID_GEOMETRY_ID ){
      // no hit.
      return radiance_bg( mlt, pos, dir, lambda );
    }
    
    // hit to the scene!
    if( dist )
      *dist = ray.tfar;
    
    VECTORMATH::Vector hit (
      ray.org[0] + ray.dir[0] * ray.tfar,
      ray.org[1] + ray.dir[1] * ray.tfar,
      ray.org[2] + ray.dir[2] * ray.tfar );
    
    VECTORMATH::Vector geometry_normal = model_.normal( ray ); // geometry_normal( ray.Ng[0], ray.Ng[1], ray.Ng[2] ); // プリミティブの法線.
    // VECTORMATH::Vector geometry_normal( ray.Ng[0], ray.Ng[1], ray.Ng[2] ); // プリミティブの法線.
    geometry_normal.normalize();
    VECTORMATH::Vector ray_normal;
    bool               into = false;
    
    if( geometry_normal.dot( dir ) < 0.f ){ // レイは物体に入り込む方向.
      ray_normal = geometry_normal *  1.;
      into = true;
    }else{
      ray_normal = geometry_normal * -1.; // レイは物体から出る方向.
    }
    VECTORMATH::Vector up(0,1,0);
    VECTORMATH::Vector binormal = up.cross( ray_normal );
    VECTORMATH::Vector tangent  = ray_normal.cross( binormal );
    
    MODEL::Material *mat = model_.material( ray );
    /* COLOR3::Color3 normalRGB( COLOR3::RGB,
                              into ? 1. : 0.,
                              0., 0. );
    return SPECTRUM::getReflectance( normalRGB, lambda );  debug out. */
    
    double russianRandom      = mlt.NextSample();
    double materialColor      = SPECTRUM::getReflectance( mat->color_, lambda ); // この波長における色.
    double russianProbability = materialColor;
    double emission = mat->emitter( lambda );
    if( depth > 5 ){
      if( russianRandom >= russianProbability )
        return emission; //ロシアンルーレット
    }else{
      russianProbability = 1.f;
    }
    
    double materialRandom = mlt.NextSample();
    if( materialRandom < mat->refraction_ ) {

      // 屈折.
      double matIor  = mat->ior( lambda );
      double iorNow  = into ? 1. : matIor;
      double iorNext = into ? matIor : 1.;
      VECTORMATH::Vector refractDir;
      VECTORMATH::Vector reflectDir;
      bool refracted = refract( dir, ray_normal, iorNow/iorNext, refractDir );
      reflectDir = reflect( dir, ray_normal );

      if( !refracted ){

        // 全反射してしまった.
        return emission + radiance( depth+1, mlt, hit, reflectDir, lambda ) / russianProbability;
        //return SPECTRUM::getReflectance( COLOR3::Color3( COLOR3::RGB,1,0,0 ), lambda );

      }else{

        double a  = iorNow - iorNext;
        double b  = iorNow + iorNext;
        double R0 = (a*a)/(b*b);
        double c  = 1.0 + dir.dot( ray_normal );
        assert( 0. <= c && c <= 1. );
        double Re = R0 + (1.0 - R0) * pow( c, 5. ); // 反射からの寄与.
        double nnt2 = (iorNow/iorNext)*(iorNow/iorNext);
        double Tr = (1. - Re)*nnt2; // 屈折からの寄与.
        double probability = 0.25 + 0.5 * Re;

        if( mlt.NextSample() < probability ){
          //return SPECTRUM::getReflectance( COLOR3::Color3( COLOR3::RGB,0,1,0 ), lambda );
          // 反射側を追う. 色は乗らない.
          return emission + radiance( depth+1, mlt, hit, reflectDir, lambda ) * Re / (probability * russianProbability);
        }else{
          //return SPECTRUM::getReflectance( COLOR3::Color3( COLOR3::RGB,0,0,1 ), lambda );
          // 屈折側を追う.
          double dist;
          double scale = 100.;
          double rad = radiance( depth+1, mlt, hit, refractDir, lambda , &dist );
          double transmittance = 1. - (dist / scale);
          if( transmittance < 0. )
            transmittance = 0.;
          return emission + materialColor * rad * transmittance * Tr / ((1.-probability) * russianProbability);
        }
      }
    } /*  else{     // DISNEY BRDFへのトライ.
      //return SPECTRUM::getReflectance( COLOR3::Color3( COLOR3::RGB,0,0,1 ), lambda );
      VECTORMATH::Vector nextDir;
      do{
        nextDir.set( mlt.NextSample(), mlt.NextSample(), mlt.NextSample() );
      }while( nextDir.length2() > 1.);
      if( nextDir.dot( ray_normal ) < 0. )
        nextDir = nextDir * -1.;

      assert( ray_normal.dot( nextDir ) > 0. );
      assert( ray_normal.dot( dir * -1. ) > 0. );
      
      double b = mat->DisneyBRDF( nextDir, dir * -1., ray_normal, tangent, binormal, lambda );
      return radiance( depth+1, mlt, hit, nextDir, lambda) * b; } */
    else if( materialRandom < mat->diffuse_ ){
      // diffuse反射

      // NEE.
      const int shadowRays = 1;
      double directLight = 0.;
      for (int i = 0; i < shadowRays; i ++ )
        directLight += direct_radiance( mlt, hit, ray_normal, lambda ) / shadowRays;

      double phi= 2. * M_PI * mlt.NextSample();
      double r  = mlt.NextSample();

      double x  = cos( phi ) * sqrt( 1.0 - r );
      double y  = sin( phi ) * sqrt( 1.0 - r );
      double z  = sqrt( r );

      VECTORMATH::Vector basis2 = ray_normal;
      VECTORMATH::Vector basis1( 0.f, 0.f, 0.f );
      if( -0.6 < ray_normal.x && ray_normal.x < 0.6 )
        basis1.x = 1.;
      else if( -0.6 < ray_normal.y && ray_normal.y < 0.6 )
        basis1.y = 1.;
      else if( -0.6 < ray_normal.z && ray_normal.z < 0.6 )
        basis1.z = 1.;
      else
        basis1.x = 1.;

      VECTORMATH::Vector basis0 = basis1.cross( basis2 );
      basis0.normalize();
      basis1 = basis2.cross( basis0 );
      basis1.normalize();

      VECTORMATH::Vector dx = basis0 * x;
      VECTORMATH::Vector dy = basis1 * y;
      VECTORMATH::Vector dz = basis2 * z;
      VECTORMATH::Vector nextDir = dx + dy + dz;

      return emission + (materialColor * (directLight + radiance( depth+1, mlt, hit, nextDir, lambda ))) / russianProbability;

    }else if( materialRandom < mat->reflection_ ){

      // specular反射
      VECTORMATH::Vector nextDir = reflect( dir, ray_normal );
      return emission + radiance( depth+1, mlt, hit, nextDir, lambda ) / russianProbability;
      
    }
    // 吸収.
    
    return emission;
  }
  double pathtrace( KelemenMLT& mlt, NUMA::Vector& cpos, NUMA::Vector& cdir, double lambda ){
    VECTORMATH::Vector pos,dir;
    model_.createCameraRay( cpos, cdir, pos, dir ); // カメラローカルからワールド座標にマップする.
    return radiance( 0, mlt, pos, dir, lambda );
  }
};

// MLTのために新しいパスをサンプリングする関数。
PathSample generate_new_path(
  KelemenMLT &mlt,
  NUMA::Lens &lens,
  Scene& scene,
  const int width,
  const int height )
{
  double wavelength = 380. + 400. * mlt.NextSample();
  double weight = width * height * 4.;
  double fx = mlt.NextSample() - 0.5;
  double fy = mlt.NextSample() - 0.5;

  NUMA::Vector pos,dir;
  double tryout = 0.;
  bool   result = false;
  tryout += 1.;
  result = lens.createRaySample(
    mlt,
    lens.getImageSurfaceZ(),
    (fy * lens.getImageSurfaceR())*height/width,
    fx * lens.getImageSurfaceR(),
    wavelength, pos, dir);
  double light = result ? scene.pathtrace( mlt, pos, dir, wavelength ) : 0.;
  return PathSample( (fx+0.5)*(width-1), (fy+0.5)*(height-1), wavelength, light / tryout, 1.0 / (1.0 / weight));
}


class KelemenMLTRender {
  
public:

  KelemenMLT mlt_;
  double b_;
  double p_large_;
  int accept_;
  int reject_;
  int width_;
  int height_;

  std::vector<FRAMEBUFFER::Color> image_;
  PathSample old_path_;
  
  void init_mlt( int seed, const int width, const int height, NUMA::Lens& lens, Scene& scene ) {
    printf("%d\n",omp_get_thread_num());
    mlt_.xor128_.setSeed( omp_get_thread_num() + seed );
    
    width_ = width;
    height_= height;
    
    image_.resize(width * height);
    for(int i=0;i<width_*height_;i++)
      image_[i].set(0,0,0);

    mlt_.large_step = 1;
    while( 1 ) {
      mlt_.InitUsedRandCoords();
      PathSample sample = generate_new_path( mlt_, lens, scene, width, height );
      mlt_.global_time ++;
      while (!mlt_.primary_samples_stack.empty()) // スタック空にする
        mlt_.primary_samples_stack.pop();
      double lum = sample.luminance();
      if( lum > 0. ){
        // 有効なパスが見つかった！
        b_       = lum;
        p_large_ = 0.1;
        accept_ = 0;
        reject_ = 0;
        old_path_ = sample;
        break;
      }
    }
    printf("seed found\n");
    return;
  }

  void step( NUMA::Lens& lens, Scene& scene ){
    
    mlt_.large_step = mlt_.rand01() < p_large_;
    mlt_.InitUsedRandCoords();
    
    PathSample new_path = generate_new_path(mlt_,lens,scene,width_,height_);

    double old_lum = old_path_.luminance();
    double new_lum = new_path.luminance();
    double a = std::min(1.0, new_lum / old_lum );
    double new_path_weight = (a + mlt_.large_step) / (new_lum / b_ + p_large_);
    double old_path_weight = (1.0 - a) / (old_lum / b_ + p_large_);
    
    image_[new_path.y  * width_ + new_path.x ] = image_[new_path.y  * width_ + new_path.x ] + new_path.color() * new_path.weight  * new_path_weight;
    image_[old_path_.y * width_ + old_path_.x] = image_[old_path_.y * width_ + old_path_.x] + old_path_.color()* old_path_.weight * old_path_weight;

    if (mlt_.rand01() < a) { // 受理
      accept_ ++;
      old_path_ = new_path;
      if (mlt_.large_step)
        mlt_.large_step_time = mlt_.global_time;
      mlt_.global_time ++;
      while (!mlt_.primary_samples_stack.empty()) // スタック空にする.
        mlt_.primary_samples_stack.pop();
    } else { // 棄却
      reject_ ++;
      int idx = mlt_.used_rand_coords - 1;
      while (!mlt_.primary_samples_stack.empty()) {
        mlt_.primary_samples[idx --] = mlt_.primary_samples_stack.top();
        mlt_.primary_samples_stack.pop();
      }
    }
  }
  
  void accum( std::vector<FRAMEBUFFER::Color> & image ){
    for(int i=0;i<width_*height_;i++)
      image[i] = image[i] + image_[i];
  }
};

void tonemap( std::vector<FRAMEBUFFER::Color> &accum, std::vector<unsigned char> &rgbx ){
  double scale = 255. ;
  rgbx.resize( accum.size() * 4 );
  for(int i=0;i<accum.size();i++){
    int ir = accum[i].x * scale;
    int ig = accum[i].y * scale;
    int ib = accum[i].z * scale;
    rgbx[i*4]= ir > 255 ? 255:ir;
    rgbx[i*4+1]= ig > 255 ? 255:ig;
    rgbx[i*4+2]= ib > 255 ? 255:ib;
    rgbx[i*4+3]= 0;
  }
}


char * modelName = "";
char * lensName  = "";
char * bgName    = "";
bool isZemax;
bool isLensText;
float focusAdjust;
float evAdjust;
int   imageSize;
float filmSize;
float anamoRatio;
float irisScale;
float timelimit;
bool  saveHdr;
bool  lensDump;
bool autoFocus;

void parseArgs( int argc, char **argv )
{
  isZemax = false;
  isLensText = false;

  imageSize = 2048;
  evAdjust = 1.;
  focusAdjust = 0.;
  filmSize = 0.;
  irisScale = 1.;
  timelimit = 900.;
  saveHdr = false;
  lensDump = false;
  autoFocus = true;
  
  for(int i=1;i<argc;i++){
    if( !strcmp( argv[i], "/saveHdr") )
      saveHdr = true;
    if( !strcmp( argv[i], "/lensDump") )
      lensDump = true;
    if( !strcmp( argv[i], "/manualFocus"))
      autoFocus = false;
    
    if( !strcmp( argv[i], "/model" ) ){
      modelName = argv[i+1];
      i++;
    }
    if( !strcmp( argv[i], "/lens" ) ){
      lensName = argv[i+1];
      isLensText = true;
      isZemax = false;
      i++;
    }
    if( !strcmp( argv[i], "/zmx" ) ){
      lensName = argv[i+1];
      isLensText = false;
      isZemax = true;
      i++;
    }
    if( !strcmp( argv[i], "/bg" ) ){
      bgName = argv[i+1];
      i++;
    }
    if( !strcmp( argv[i], "/focus" ) ){
      focusAdjust = atof( argv[i+1] );
      i++;
      printf("focus %f\n",focusAdjust);
    }
    if( !strcmp( argv[i], "/exposure" ) ){
      evAdjust = atof( argv[i+1] );
      i++;
      printf("exposure %f\n",evAdjust);
    }
    if( !strcmp( argv[i], "/filmSize" ) ){
      filmSize = atof( argv[i+1] );
      i++;
      printf("filmSize %f\n",filmSize);
    }
    if( !strcmp( argv[i], "/imageSize" ) ){
      imageSize = atoi( argv[i+1] );
      i++;
      printf("imageSize %d\n",imageSize);
    }
    if( !strcmp( argv[i], "/anamo" ) ){
      anamoRatio = atof( argv[i+1] );
      i++;
      printf("anamo %f\n",anamoRatio);
    }
    if( !strcmp( argv[i], "/iris" ) ){
      irisScale = atof( argv[i+1] );
      i++;
      printf("iris %f\n",irisScale);
    }
    if( !strcmp( argv[i], "/timelimit" ) ){
      timelimit = atof( argv[i+1] );
      i++;
      printf("timelimit %f\n",timelimit);
    }
  }
}

int main(int argc, char **argv) {

  printf("init embree.\n");
  rtcInit( NULL );
  _MM_SET_FLUSH_ZERO_MODE    ( _MM_FLUSH_ZERO_ON     );
  _MM_SET_DENORMALS_ZERO_MODE( _MM_DENORMALS_ZERO_ON );
  printf("done.\n");
  SPECTRUM::initXyzPdf();
  
  Scene scene;
  KelemenMLTRender render;
  NUMA::Lens lens;

  parseArgs( argc, argv );
  
  scene.init( modelName, bgName );
  if( isLensText )
    NUMA::LOADER::LENS::load( lensName, lens );
  else if( isZemax )
    NUMA::LOADER::ZEMAX::load( lensName, lens );
  lens.setIrisScale( irisScale );
  
  if( autoFocus ){
    double bf = NUMA::autofocus( lens, scene.model_.cameraFocus_ );
    double backfrange = lens.getImageSurfaceZ();
    printf("backfrange %f adjust %f, autofocus %f\n",backfrange,focusAdjust,bf);
    lens.setImageSurfaceZ( bf + focusAdjust );
  }else{
    double backfrange = lens.getImageSurfaceZ();
    lens.setImageSurfaceZ( backfrange + focusAdjust );
  }
  
  if( filmSize > 0. )
    lens.setImageSurfaceR( filmSize );

  if( lensDump )
  {
    char fn[129];
    sprintf(fn,"%d_0001_xz.png",time(NULL));
    lens.dump_XZ(fn,435.8,lens.getImageSurfaceZ());
  }
  lens.dump();
  printf("dumped\n");
  
  int width  = imageSize;
  int num_parallel = omp_get_max_threads();
#ifdef CHECK
  num_parallel = omp_get_max_threads() - 1;
#endif
#ifdef DEBUG
  num_parallel = 1;
#endif

  FRAMEBUFFER::FrameBuffer fb;
  int height = width * 3/4;
  fb.setup(width,height);

  printf("render in %d parallel\n",num_parallel);
  std::vector<KelemenMLTRender> renders(num_parallel);
  std::vector<FRAMEBUFFER::Color> accum( width * height );
  double div = 0.;

  srand(time(NULL));
#pragma omp parallel for
  for(int i=0;i<num_parallel;i++)
    renders[i].init_mlt(rand(),width,height,lens,scene);
  printf("mlt init\n");
  
  int s = 0;
  double steps = 0.;

  time_t start_time;
  time_t previous_time;
  time( &start_time );
  time( &previous_time );

  for(s=1;s;s++){
    
#if defined(CHECK) || defined(DEBUG)
#define MUTATIONS 0x100000
#else
#define MUTATIONS 0x400000
#endif
#pragma omp parallel for
    for(int i=0;i<num_parallel;i++){
      for(int step=0;step<MUTATIONS;step++)
        renders[i].step( lens, scene );
      printf("%d accept, %d reject (%f %% accept)\n",renders[i].accept_,renders[i].reject_,
             renders[i].accept_ * 100. / (renders[i].reject_ + renders[i].accept_) );
    }
    steps += (double)MUTATIONS;

    time_t current_time;
    time( &current_time );
    double singlestep = difftime( current_time, start_time ) / (double)s;
    bool willover = false;
    if( timelimit > 0. && difftime( current_time, start_time ) + singlestep >= timelimit ){
      printf("time will be up!\n");
      willover = true;
    }
    
#if defined(CHECK) || defined(DEBUG)
    {
#else
    if( difftime( current_time, previous_time ) >= 30.0 || willover ) {
#endif
      // save
      previous_time = current_time;

      for(int j=0;j<width*height;j++)
        accum[j].set(0,0,0);
      for(int i=0;i<num_parallel;i++)
        renders[i].accum( accum );
      for(int j=0;j<width*height;j++)
        accum[j] = accum[j] * (1./ steps);

      // レンズのために反転する.
      for(int j=0;j<height;j++){
        for(int i=0;i<width;i++){
          fb.set( width - i - 1, height - j - 1, accum[j*width+i].x, accum[j*width+i].y,accum[j*width+i].z );
        }
      }

      fb.normalize();
      char fn[129];
      if( saveHdr ){
        sprintf(fn,"%d_%08d.hdr",start_time,s);
        fb.save_hdr(fn,1.,1.);
      }
      sprintf(fn,"%d_%08d.png",start_time,s);
      fb.save_ldr(fn,evAdjust);
      printf("%d / %f\n",s,difftime(current_time,start_time));
      if( willover )
        break;
    }
  }
  rtcExit();
}

