/*
 * copyright(c) 2015 Hajime UCHIMURA.
 * do not remove these text.
 * please contact before business use.
 */

#ifndef __LENSNUMA_H
#define __LENSNUMA_H

#include "kelemenmlt.h"
#include "vectormath.h"

namespace NUMA { // レンズは沼です.

  typedef VECTORMATH::Vector3<double> Vector;

  class Intersection{
  public:
    bool   hit_;
    Vector from_;
    Vector dir_;
    Vector point_;
    Vector normal_;
    
    Intersection() : hit_(false) {;}
    void set( bool hit, const Vector& from, const Vector& dir, const Vector& point, const Vector& norm ){
      hit_   = hit;
      from_  = from;
      dir_   = dir.normal();
      point_ = point;
      normal_= norm.normal();
    }
  };

  void reflect( const Vector& I, const Vector& N, Vector& result ){
    result = I - N * 2. * N.dot( I );
  }
  
  bool refract( const Vector& I, Vector& N, double eta, Vector& result ){
    if( I.dot( N ) > 0. )
      N = N * -1.;
    double cosi  = - I.dot( N ); // dot(-i, n);
    double cost2 = 1.0 - eta * eta * (1.0 - cosi * cosi);
    if( cost2 < 0. ) return false; // 全反射.
    result = I * eta + (N * (eta * cosi - sqrtf( cost2 )));
    return true;
  }
  
  class Surface{
  public:
    
    typedef enum {
      NONE,     // 空間.
      APERTURE_HEXAGONAL, // 六角絞り
      APERTURE_CIRCLE,    // 円絞り. xyアスペクトで楕円も.
      STANDARD, // 球面レンズ
      CYLINDER_X, // シリンドリカルレンズX.
      CYLINDER_Y, // シリンドリカルレンズY. anamo.
      CYLINDER_Z, // シリンドリカルレンズZ.
    } TYPE;
    
    TYPE  type_;
    double center_;
    double radius_;
    double diameter_;
    double thickness_; // 次の面までの距離.
    double iris_x_,iris_y_; // 円絞り楕円率.
    double ior_;
    double abbe_vd_;
    double reflection_; // 反射率.
    bool   is_coated_;
    double coat_thickness_; // 275nm = 550nmの半波長.
    double coat_ior_;       // MgF2で1.38
    double roughness_;
    
    Surface(){
      init();
    }

    void init(){
      type_    = NONE;
      center_  = 0.;
      radius_  = 0.;
      diameter_= 0.;
      iris_x_  = 1.;
      iris_y_  = 1.;
      thickness_ = 0.;
      roughness_ = 0.; // 面の荒れ.
      
      ior_     = 1.; // at D light.
      abbe_vd_ = 1.; // at D light.
      reflection_ = 0.1;
      is_coated_ = false;
    }

    double ior( double lambda ) const {
      
      //return ior_ + (lambda / 1000.) / ( 100. * abbe_vd_ );
      /* 
      if( abbe_vd_ <= 0. )
        return ior_;
      
      const double A1 = -8.39627218911159000000;
      const double B1 =  0.07307694336347550000;
      const double A2 =  5.43767332268365000000;
      const double B2 = -0.04894067114512750000;
      
      double A  = B1 + A1 * abbe_vd_;
      double B  = B2 + A2 * abbe_vd_;
      return ior_ + B + A * (lambda / 1000.);
       */
      // コーシーの式.
      double B = (ior_ - 1.) / abbe_vd_ * 0.52345;
      double A = ior_ - B / 0.34522792;
      double C = lambda / 1000.;
      double ret = A + B / (C*C);
      assert( ret == ret );
      return ret;
    }

    // coat_thickness must be in nm.
    // lambda is in nm.
    double single_coat_reflect( double lambda, double ior_in, double coat_ior, double coat_thickness, const Vector& dir, const Vector& norm ){
      double cosTheta1 = dir.dot( norm );
      double cosTheta2 = (ior_in / coat_ior) * sqrtf( (coat_ior*coat_ior)/(ior_in*ior_in) - (1. - cosTheta1 * cosTheta1 ) );
      double distance_diff = 2. * coat_thickness * (coat_ior / ior_in) * cosTheta2;
      double m = fmod( distance_diff, lambda ) / lambda; // ズレ幅0-1,
      return cos( m * M_PI ); // 0.5のとき0になるような値.
    }

    double reflection( double lambda, double ior_now, double ior_next, const Vector& dir, const Vector& norm, double& Re, double& Tr ){

      if( is_coated_ ) // シングルコート.
        return single_coat_reflect( lambda, ior_now, coat_ior_, coat_thickness_, dir, norm );
      // コーティング無し
      double a  = ior_now - ior_next;
      double b  = ior_now + ior_next;
      double R0 = (a*a)/(b*b);
      double t  = dir.dot( norm );
      double c  = (t < 0.) ? 1.0 + t : 1.0 - t;
      //assert( 0. <= c && c <= 1. );
      Re = R0 + (1.0 - R0) * pow( c, 5. ); // 反射からの寄与.
      double nnt2 = (ior_now/ior_next)*(ior_now/ior_next);
      Tr = (1. - Re)*nnt2; // 屈折からの寄与.
      return Re;
    }

    double inline maxf( double a, double b ){ return (a>b) ? a:b; }
    bool intersect( const Vector& origin, const Vector& dir, Intersection& hit, double irisScale = 1. ){

      if( radius_ == 0. ){
        double t = (center_ - origin.x) / dir.x;
        Vector p = origin + dir * t;
        double r = sqrtf( p.z * p.z + p.y * p.y ); // 面上の高さ.

        if( r < diameter_ ){
          if( type_ == APERTURE_CIRCLE ){
            double py = p.y * iris_y_ * irisScale;
            double pz = p.z * iris_x_ * irisScale;
            double d  = ( py*py + pz*pz ) - diameter_ * diameter_; // circular iris.
            if( d < 0. ){
              hit.set( true, origin, dir, p, Vector( dir.x > 0. ? -1. : 1., 0., 0. ) );
              return true;
            }
          }else if( type_ == APERTURE_HEXAGONAL ){
            double iris = diameter_ / (iris_x_ * irisScale);
            ////// hexagonal iris
            double ax = p.y > 0. ? p.y : -p.y;
            double ay = p.z > 0. ? p.z : -p.z;
            double d  = maxf( ax - iris, maxf( ax + ay * 0.57735, ay * 1.1547 ) - iris );
            if( d < 0. ){
              hit.set( true, origin, dir, p, Vector( dir.x > 0. ? -1. : 1., 0., 0. ) );
              return true;
            }
          }else{
            hit.set( true, origin, dir, p, Vector( dir.x > 0. ? -1. : 1., 0., 0. ) );
            return true;
          }
        }
        // 絞られた.
        return false;
      }

      if( type_ == CYLINDER_Y ){
        
        Vector center( center_, 0., 0. );
        
        // Y軸シリンドリカルレンズ. X^2+Z^2 = radius^2
        double px = origin.x - center_;
        double a = dir.x * dir.x + dir.z * dir.z;
        double b = 2. * px * dir.x + 2. * origin.z * dir.z;
        double c = px * px + origin.z * origin.z - radius_ * radius_;
        double d = b * b - 4. * a * c ;
        if( d < 0. )
          return false;
        
        double s  = sqrt( d );
        double t1 = ( -b + s ) / (2. * a);
        double t2 = ( -b - s ) / (2. * a);
        if( t1 <= 0. && t2 <= 0. )
          return false;
        
        if( t2 >= 0. && t1 >= 0. ){
          Vector p1 = origin + dir * t1;
          Vector p2 = origin + dir * t2;
          Vector n1 = p1 - center; n1.y = 0.; n1.normalize();
          Vector n2 = p2 - center; n2.y = 0.; n2.normalize();
          if( radius_ < 0. ){
            // 右(x+)向きに凸レンズ.
            // n.xが正の方を採用.
            if( n1.x > 0. && (p1.y < diameter_ && p1.y > -diameter_) ){
              hit.set( true, origin, dir, p1, n1 );
              return true;
            }else if( n2.x > 0. && (p2.y < diameter_ && p2.y > -diameter_) ){
              hit.set( true, origin, dir, p2, n2 );
              return true;
            }
            return false;
          }else{
            // 左(x-)向きに凸レンズ.
            // n.xが負の方を採用.
            if( n1.x < 0. && (p1.y < diameter_ && p1.y > -diameter_) ){
              hit.set( true, origin, dir, p1, n1 );
              return true;
            }else if( n2.x < 0. && (p2.y < diameter_ && p2.y > -diameter_)){
              hit.set( true, origin, dir, p2, n2 );
              return true;
            }
            return false;
          }
          
        }else{
          
          double t;
          t = (t1 > t2) ? t1 : t2; // どちらかが負なら大きい方が距離.
          
          Vector p = origin + dir * t;
          if( p.y > diameter_ || p.y <= -diameter_ )
            return false;
          Vector n;
          n.x = p.x - center.x;
          n.z = p.z - center.z;
          n.y = 0.;
          n.normalize();
          if( radius_ < 0. && n.x > 0. ){
            hit.set( true, origin, dir, p, n );
            return true;
          }else if( radius_ > 0. && n.x < 0. ){
            hit.set( true, origin, dir, p, n );
            return true;
          }
          
          return false;
        }
        
        return false;
        
      }else if( type_ == STANDARD ){
        
        // 球面レンズ.
        
        Vector center( center_, 0., 0. );
        Vector oc = origin - center;

        double a = dir.dot( dir );
        double b = 2. * dir.dot( oc );
        double c = oc.dot( oc ) - radius_ * radius_;
        double d2 = b*b - 4. * a * c;

        if( d2 < 0. )
          return false;

        double d  = sqrtf( d2 );
        double t1 = (-b - d) / (2. * a);
        double t2 = (-b + d) / (2. * a);
        if( t1 <= 0. && t2 <= 0. )
          return false;

        double t;
        if( t2 >= 0. && t1 >= 0. ){
          // 両方正なので両方テストする必要がある.
          Vector p1 = origin + dir * t1;
          Vector n1 = p1 - center;
          Vector p2 = origin + dir * t2;
          Vector n2 = p2 - center;
          if( radius_ < 0. ){
            // 右(x+)向きに凸レンズ.
            // n.xが正の方を採用.
            if( n1.x > 0. ){
              double r = sqrtf( p1.y * p1.y + p1.z * p1.z );
              if( r < diameter_ ){
                hit.set( true, origin, dir, p1, n1 );
                return true;
              }
            }else if( n2.x > 0. ){
              double r = sqrtf( p2.y * p2.y + p2.z * p2.z );
              if( r < diameter_ ){
                hit.set( true, origin, dir, p2, n2 );
                return true;
              }
            }
            return false;
          }else{
            // 左(x-)向きに凸レンズ.
            // n.xが負の方を採用.
            if( n1.x < 0. ){
              double r = sqrtf( p1.y * p1.y + p1.z * p1.z );
              if( r < diameter_ ){
                hit.set( true, origin, dir, p1, n1 );
                return true;
              }
            }else if( n2.x < 0. ){
              double r = sqrtf( p2.y * p2.y + p2.z * p2.z );
              if( r < diameter_ ){
                hit.set( true, origin, dir, p2, n2 );
                return true;
              }
            }
            return false;
          }
        }else{
          t = (t1 > t2) ? t1 : t2; // どちらかが負なら大きい方が距離.
          Vector p = origin + dir * t;
          Vector n = p - center;
          if( radius_ < 0. ){
            // 右(x+)向きに凸レンズ.
            // n.xが正の方を採用.
            if( n.x > 0. ){
              double r = sqrtf( p.y * p.y + p.z * p.z );
              if( r < diameter_ ){
                hit.set( true, origin, dir, p, n );
                return true;
              }
            }
            return false;
          }else{
            // 左(x-)向きに凸レンズ.
            // n.xが負の方を採用.
            if( n.x < 0. ){
              double r = sqrtf( p.y * p.y + p.z * p.z );
              if( r < diameter_ ){
                hit.set( true, origin, dir, p, n );
                return true;
              }
            }
            return false;
          }
        }
      }
      return false;
    }
  };
  
  typedef std::vector<Surface> SurfaceSet;
  
  class Lens {
  public:
    
    SurfaceSet surfaces_;
    double imageSurfaceZ_;
    double imageSurfaceR_; // 像面高さ.
    double irisScale_;
    
    Lens() {
      surfaces_.clear();
      imageSurfaceZ_ = 100.;
      imageSurfaceR_ = 100.;
      irisScale_ = 1.;
    }

    /* 近軸厚肉焦点距離
    double singleFocus( Surface& front, Surface& rear, double lambda, double& bf ){
      // 1/f = (n - 1)( 1/R1 - 1/R2 + tc*(n-1)/R1R2n );
      double ior = front.ior( lambda );
      double c1  = (front.radius_ == 0.) ? 0. : 1./front.radius_; // 平面の時は0.
      double c2  = (rear.radius_  == 0.) ? 0. : 1./rear.radius_;
      double c3  = (front.radius_ == 0. || rear.radius_ == 0.) ? 0. : front.thickness_ * (ior-1.)/(ior*front.radius_*rear.radius_);
      double f   = 1.f/((ior-1.)*(c1-c2+c3));
      bf = f * ((front.radius_==0.) ? 1. : ( 1. - front.thickness_ * (ior-1.)/(ior*front.radius_)));
    }*/
    
    double getImageSurfaceR(void){
      return imageSurfaceR_;
    }
    void setImageSurfaceR(double r){
      imageSurfaceR_ = r;
    }
    double getImageSurfaceZ(void){
      return imageSurfaceZ_;
    }
    void setImageSurfaceZ(double z){
      imageSurfaceZ_ = z;
    }
    void setIrisScale( double i ){
      irisScale_ = i;
    }
    
    // create a ray sample from lens cylinders.
    bool createRaySample( KelemenMLT &mlt, double fx, double fy, double fz, double lambda, Vector& origin, Vector& dir ){
      Intersection hit;

      double ior_now = 1.; // 空気.
      Vector p( fx, fy, fz );
      Vector lookat;
      Vector d( -1., 0., 0. );

      double r1, r2;
      do{
        r1 = mlt.NextSample() * 2. - 1.;
        r2 = mlt.NextSample() * 2. - 1.;
      }while( r1 * r1 + r2 * r2 > 1.);
      //printf("p %f,%f,%f\nr %f %f\n",p.x,p.y,p.z,r1,r2);

      Surface& last( surfaces_.back() );
      lookat.set( last.center_ - last.radius_,
                  last.diameter_ * r1,
                  last.diameter_ * r2);
      //printf("l %f,%f,%f\n",lookat.x,lookat.y,lookat.z);
      d = lookat - p;
      d.normalize();
      // printf("d %f,%f,%f\n",d.x,d.y,d.z);

      int index = surfaces_.size() - 1 ;
      int delta = -1; // 外にでる方向.
      
      while( 1 ) {

        // printf("tracing surface %d with %f,%f,%f\n",index,d.x,d.y,d.z);
        // 次の面がもうない.
        if( index < 0 || index >= surfaces_.size() ){
          // printf("no next surface\n");
          break;
        }
        
        bool result = surfaces_[ index ].intersect( p, d, hit, irisScale_ );
        
        if( result ){
          
          // 面にあたった.
          if( surfaces_[ index ].type_ == Surface::APERTURE_CIRCLE ||
              surfaces_[ index ].type_ == Surface::APERTURE_HEXAGONAL ){
            
            // 絞りだ.
            p = hit.point_;
            // dは変化なし.
            // ior変化なし.
            index = index + delta;
            
          }else{
            // 右に進むなら今のガラスの屈折率に、
            // 左に進むなら一個前のガラスの屈折率に
            // 左がもうないなら空気の屈折率に.
            double nextIor;
            if( d.x > 0. ){
              nextIor = surfaces_[ index ].ior( lambda );
            }else if( index > 0 ){
              nextIor = surfaces_[ index - 1 ].ior( lambda );
            }else{
              nextIor = 1.;
            }
            
            Vector normal = hit.normal_;
            double roughness = surfaces_[index].roughness_;
            if( roughness > 0. ){
              if( mlt.NextSample() < roughness ){
                Vector n2 ( mlt.NextSample(), mlt.NextSample(), mlt.NextSample() );
                n2.normalize();
                if( n2.dot( normal ) < 0. )
                  n2 = n2 * -1.;
                normal = n2;
                /*
                double t = (mlt.NextSample() * 2. - 1.) * M_PI; // theta
                double x = normal.x;
                double z = normal.z;
                double s = sin( t );
                double c = cos( t );
                normal.x = x * c - z * s;
                normal.z = x * s + z * c;
                  */
                //normal.normalize();
              }
            }

            double Re, Tr;
            surfaces_[ index ].reflection( lambda, ior_now, nextIor, hit.dir_,  normal, Re, Tr );
            if( mlt.NextSample() < Re ){
              
              // 反射しちゃった.
              Vector nextDir;
              reflect( d, normal, nextDir );
              d = nextDir.normal();
              // iorは変化なし.
              p = hit.point_;
              delta = -delta;
              index = index + delta;
              
            }else{
              
              // 屈折する.
              Vector nextDir;
              
              if( !refract( d, normal, ior_now / nextIor, nextDir ) ){
                // 全反射してしまった.
                reflect( d, normal, nextDir );
                delta = -delta;
              }
              
              d = nextDir.normal() ;
              p = hit.point_;
              
              ior_now = nextIor;
              index = index + delta;
            }
          }
        }else{
          // printf("not intersecting\n");
          break;
        }
      }
      if( index == -1 ){
        // 外に出た!
        origin = p * 0.001; // mm>m換算する.
        dir    = d;
        return true;
      }
      // 外に出なかった…
      return false;
    }
    
    bool traceRay( double fx, double fy, double fz, double x, double y, double lambda, Vector& origin, Vector& dir ){
      Intersection hit;
      double ior_now = 1.; // 空気.
      Vector p( fx, fy, fz );
      Vector lookat;
      Vector d( -1., 0., 0. );
      Surface& last( surfaces_.back() );
      lookat.set( last.center_ - last.radius_, last.diameter_ * x, last.diameter_ * y );
      //printf("p %f,%f,%f\n",fx,fy,fz);
      //printf("l %f,%f,%f\n",lookat.x,lookat.y,lookat.z);
      d = lookat - p;
      d.normalize();
      int index = surfaces_.size() - 1 ;
      int delta = -1; // 外にでる方向.
      while( 1 ) {
        if( index < 0 || index >= surfaces_.size() ){
          break;
        }
        bool result = surfaces_[ index ].intersect( p, d, hit, irisScale_ );
        //printf("%d:%d\n",index,result);
        if( result ){
          // 面にあたった.
          if( surfaces_[ index ].type_ == Surface::APERTURE_CIRCLE ||
              surfaces_[ index ].type_ == Surface::APERTURE_HEXAGONAL ){
            p = hit.point_;
            index = index + delta;
          }else{
            // 右に進むなら今のガラスの屈折率に、
            // 左に進むなら一個前のガラスの屈折率に
            // 左がもうないなら空気の屈折率に.
            double nextIor;
            if( d.x > 0. ){
              nextIor = surfaces_[ index ].ior( lambda );
            }else if( index > 0 ){
              nextIor = surfaces_[ index - 1 ].ior( lambda );
            }else{
              nextIor = 1.;
            }
            // 屈折する.
            Vector nextDir;
            refract( d, hit.normal_, ior_now / nextIor, nextDir );
            d = nextDir.normal() ;
            p = hit.point_;
            
            ior_now = nextIor;
            index = index + delta;
          }
        }else{
          //printf("not intersecting\n");
          break;
        }
      }
      if( index == -1 ){
        origin = p * 0.001; // mm>m換算する.
        dir    = d;
        return true;
      }
      return false;
    }

    void drawLine(FRAMEBUFFER::FrameBuffer& fb, double x1, double y1, double x2, double y2, double scale, FRAMEBUFFER::Color& color ) {
      if( x1 > x2 ){
        double t;
        t = x1; x1 = x2; x2 = t;
        t = y1; y1 = y2; y2 = t;
      }
      int w = fb.width_/2;
      int h = fb.height_/2;
      int ix1 = w + x1 * scale ;
      int ix2 = w + x2 * scale ;
      int iy1 = h + y1 * scale ;
      int iy2 = h + y2 * scale ;
      int dx = abs( ix2-ix1 );
      int dy = abs( iy2-iy1 );
      if( dx > dy ){
        for(int x=ix1;x<ix2;x+=(ix1<ix2)?1:-1){
          int y = iy1 + (iy2-iy1)*(x-ix1)/(ix2-ix1);
          if( 0 <= x && x < fb.width_ && 0 <= y && y < fb.height_ ){
            fb.film_  [ y * fb.width_ + x ]  = fb.film_[ y * fb.width_ + x ] + color;
            fb.weight_[ y * fb.width_ + x ] += 1.;
          }
        }
      }else{
        for(int y=iy1;y<iy2;y+=(iy1<iy2)?1:-1){
          int x = ix1 + (ix2-ix1)*(y-iy1)/(iy2-iy1);
          if( 0 <= x && x < fb.width_ && 0 <= y && y < fb.height_ ){
            fb.film_  [ y * fb.width_ + x ]  = fb.film_[ y * fb.width_ + x ] + color;
            fb.weight_[ y * fb.width_ + x ] += 1.;
          }
        }
      }
    }
    
    void dump( void ){
      for(int i=0;i<surfaces_.size();i++){
        printf("%d type %d %f %f %f %f %f %f %f %f %f\n",i,
               surfaces_[i].type_,
               surfaces_[i].diameter_,
               surfaces_[i].radius_,
               surfaces_[i].center_,surfaces_[i].ior_,surfaces_[i].abbe_vd_,
               surfaces_[i].roughness_,
               surfaces_[i].reflection_,
               surfaces_[i].iris_x_,surfaces_[i].iris_y_);
      }
    }
    void dump_XZ( char *fn, double lambda, double backfocus ){
      FRAMEBUFFER::FrameBuffer fb;
      RANDOM::MTSequence mts;
      mts.init_genrand( time(NULL) );
      
      fb.setup( 4096, 2048 );
      //fb.setup( 1024, 256 );
      double scale = (fb.width_/2) / backfocus;
      
      for(int y=0;y<fb.height_;y++) {
        fb.film_  [ y * fb.width_ + (fb.width_/2) ].set( 1.f,0.f,0.f );
        fb.weight_[ y * fb.width_ + (fb.width_/2) ] = 1.;
      }
      for(int x=0;x<fb.width_;x++) {
        fb.film_  [ (fb.height_/2) * fb.width_ + x ].set( 0.,1.,0. );
        fb.weight_[ (fb.height_/2) * fb.width_ + x ] = 1.;
      }
      for(int i=0;i<surfaces_.size();i++){
        double ior = surfaces_[i].ior_;
        double vd  = surfaces_[i].abbe_vd_;
        double B = (ior - 1.) / vd * 0.52345;
        double A = ior - B / 0.34522792;
        if( surfaces_[i].radius_ != 0. ){
          int idiam = (int)( surfaces_[i].diameter_ * scale );
          for(int iy = (fb.height_/2) - idiam; iy < (fb.height_/2) + idiam; iy++ ) {
            if( iy < 0 || iy >= fb.height_ )
              continue;
            double y = (iy - (fb.height_/2)) / scale ; // x = sqrt( r^2 - y^2 )
            double x = sqrt( surfaces_[i].radius_ * surfaces_[i].radius_ - y * y );
            if( surfaces_[i].radius_ > 0. )x = -x;
            int ix = (fb.width_/2) + ( x + surfaces_[i].center_ ) * scale ;
            if( ix < 0 || ix > fb.width_ )
              continue;
            fb.film_  [ iy * fb.width_ + ix ].set( 0., 1., 1. );
            fb.weight_[ iy * fb.width_ + ix ] = 1.;
          }
        }else{
          int idiam = (int)( surfaces_[i].diameter_ * scale );
          for(int iy = (fb.height_/2) - idiam; iy < (fb.height_/2) + idiam; iy++ ) {
            if( iy < 0 || iy >= fb.height_ )
              continue;
            double y = (iy - (fb.height_/2)) / scale ; // x = sqrt( r^2 - y^2 )
            double x = 0.f;
            int ix = (fb.width_/2) + ( x + surfaces_[i].center_ ) * scale ;
            if( ix < 0 || ix > fb.width_ )
              continue;
            fb.film_  [ iy * fb.width_ + ix ].set( 1., 1., 0. );
            fb.weight_[ iy * fb.width_ + ix ] = 1.;
          }
        }
      }
      Vector p;
      Vector d;
      Intersection hit;

      double ior_now = 1.; // 空気.
      
      int samples = 100;
      for(int i = 0 ; i < samples; i ++ ){
        Vector p( backfocus, 0., 0. );
        if( i < samples/2 )
          p.z = imageSurfaceR_ / 2.;
        
        Vector lookat;
        Vector d( -1., 0., 0. );

        double r1, r2;
        do{
          r1 = 0.;// mts.genrand_real1() * 2. - 1.;
          r2 = mts.genrand_real1() * 2. - 1.;
        }while( r1 * r1 + r2 * r2 > 1.);
        
        Surface& last( surfaces_.back() );
        lookat.set( last.center_ - last.radius_,
                    last.diameter_ * r1,
                    last.diameter_ * r2);
        d = lookat - p;
        d.normalize();
        
        int index = surfaces_.size() - 1 ;
        int delta = -1;
        int count = 0;
        while( 1 ) {
          if( index < 0 || index >= surfaces_.size() ){
            break;
          }
          bool result = surfaces_[ index ].intersect( p, d, hit, irisScale_ );
          if( result ){
            drawLine( fb, p.x, p.z,
                      hit.point_.x, hit.point_.z,
                      scale,
                      ( ior_now <= 1.0 ) ? FRAMEBUFFER::Color(1.,1.,1.):FRAMEBUFFER::Color(0.,0.,1.));
            if( surfaces_[ index ].type_ == Surface::APERTURE_CIRCLE ||
                surfaces_[ index ].type_ == Surface::APERTURE_HEXAGONAL ){
              p = hit.point_;
              index += delta;
            }else{
              double nextIor;
              if( d.x > 0. ){
                nextIor = surfaces_[ index ].ior( lambda );
              }else if( index > 0 ){
                nextIor = surfaces_[ index - 1 ].ior( lambda );
              }else{
                nextIor = 1.;
              }

              double Re, Tr;
              surfaces_[ index ].reflection( lambda, ior_now, nextIor, hit.dir_,  hit.normal_, Re, Tr );
              if( mts.genrand_realF() < Re ){
                Vector nextDir;
                reflect( d, hit.normal_, nextDir );
                d = nextDir.normal();
                p = hit.point_;
                delta = -delta;
                index += delta;
              }else{
                Vector nextDir;
                refract( d, hit.normal_, ior_now / nextIor, nextDir );
                d = nextDir.normal() ;
                p = hit.point_;
                ior_now = nextIor;
                index += delta;
              }
            }
          }else{
            break;
          }
        }
        if( index == -1 ){
          drawLine( fb,
                    p.x , p.z,
                    p.x + d.x * 300.,
                    p.z + d.z * 300.,
                    scale, FRAMEBUFFER::Color( 1.,0.,1.) );
        }
      }
      fb.normalize();
      fb.save_png(fn,1.,1.);
    }
    
    void dump_XY( char *fn, double lambda, double backfocus ){
      FRAMEBUFFER::FrameBuffer fb;
      RANDOM::MTSequence mts;
      mts.init_genrand( time(NULL) );
      
      fb.setup( 4096, 2048 );
      double scale = (fb.width_/2) / backfocus;
      
      for(int y=0;y<fb.height_;y++) {
        fb.film_  [ y * fb.width_ + (fb.width_/2) ].set( 1.f,0.f,0.f );
        fb.weight_[ y * fb.width_ + (fb.width_/2) ] = 1.;
      }
      for(int x=0;x<fb.width_;x++) {
        fb.film_  [ (fb.height_/2) * fb.width_ + x ].set( 0.,1.,0. );
        fb.weight_[ (fb.height_/2) * fb.width_ + x ] = 1.;
      }
      for(int i=0;i<surfaces_.size();i++){
        double ior = surfaces_[i].ior_;
        double vd  = surfaces_[i].abbe_vd_;
        double B = (ior - 1.) / vd * 0.52345;
        double A = ior - B / 0.34522792;
        if( surfaces_[i].radius_ != 0. ){
          int idiam = (int)( surfaces_[i].diameter_ * scale );
          for(int iy = (fb.height_/2) - idiam; iy < (fb.height_/2) + idiam; iy++ ) {
            if( iy < 0 || iy >= fb.height_ )
              continue;
            double y = (iy - (fb.height_/2)) / scale ; // x = sqrt( r^2 - y^2 )
            double x = sqrt( surfaces_[i].radius_ * surfaces_[i].radius_ - y * y );
            if( surfaces_[i].type_ == Surface::CYLINDER_Y )
              x = 0.;
            
            if( surfaces_[i].radius_ > 0. )x = -x;
            int ix = (fb.width_/2) + ( x + surfaces_[i].center_ ) * scale ;
            if( ix < 0 || ix > fb.width_ )
              continue;
            fb.film_  [ iy * fb.width_ + ix ].set( 0., 1., 1. );
            fb.weight_[ iy * fb.width_ + ix ] = 1.;
          }
        }else{
          int idiam = (int)( surfaces_[i].diameter_ * scale );
          for(int iy = (fb.height_/2) - idiam; iy < (fb.height_/2) + idiam; iy++ ) {
            if( iy < 0 || iy >= fb.height_ )
              continue;
            double y = (iy - (fb.height_/2)) / scale ; // x = sqrt( r^2 - y^2 )
            double x = 0.f;
            int ix = (fb.width_/2) + ( x + surfaces_[i].center_ ) * scale ;
            if( ix < 0 || ix > fb.width_ )
              continue;
            fb.film_  [ iy * fb.width_ + ix ].set( 1., 1., 0. );
            fb.weight_[ iy * fb.width_ + ix ] = 1.;
          }
        }
      }
      Vector p;
      Vector d;
      Intersection hit;

      double ior_now = 1.; // 空気.
      
      int samples = 100;
      for(int i = 0 ; i < samples; i ++ ){
        Vector p( backfocus, 0., 0. );
        if( i < samples/2 )
          p.y = imageSurfaceR_ / 2.;
        
        Vector lookat;
        Vector d( -1., 0., 0. );

        double r1, r2;
        do{
          r1 = mts.genrand_real1() * 2. - 1.;
          r2 = 0.; // mts.genrand_real1() * 2. - 1.;
        }while( r1 * r1 + r2 * r2 > 1.);
        
        Surface& last( surfaces_.back() );
        lookat.set( last.center_ - last.radius_,
                    last.diameter_ * r1,
                    last.diameter_ * r2);
        d = lookat - p;
        d.normalize();
        
        int index = surfaces_.size() - 1 ;
        int count = 0;
        while( 1 ) {
          if( index < 0 || index >= surfaces_.size() ){
            break;
          }
          bool result = surfaces_[ index ].intersect( p, d, hit, irisScale_ );
          if( result ){
            drawLine( fb, p.x, p.y,
                      hit.point_.x, hit.point_.y,
                      scale,
                      ( ior_now <= 1.0 ) ? FRAMEBUFFER::Color(1.,1.,1.):FRAMEBUFFER::Color(0.,0.,1.));
            if( surfaces_[ index ].type_ == Surface::APERTURE_CIRCLE ||
                surfaces_[ index ].type_ == Surface::APERTURE_HEXAGONAL ){
              p = hit.point_;
              index = ( d.x > 0. ) ? index + 1 : index - 1;
            }else{
              double nextIor;
              if( d.x > 0. ){
                nextIor = surfaces_[ index ].ior( lambda );
              }else if( index > 0 ){
                nextIor = surfaces_[ index - 1 ].ior( lambda );
              }else{
                nextIor = 1.;
              }

              double Re, Tr;
              surfaces_[ index ].reflection( lambda, ior_now, nextIor, hit.dir_,  hit.normal_, Re, Tr );
              if( mts.genrand_realF() < Re ){
                Vector nextDir;
                reflect( d, hit.normal_, nextDir );
                d = nextDir.normal();
                p = hit.point_;
                index = ( d.x > 0. ) ? index + 1 : index - 1;
              }else{
                Vector nextDir;
                refract( d, hit.normal_, ior_now / nextIor, nextDir );
                d = nextDir.normal() ;
                p = hit.point_;
                ior_now = nextIor;
                index = ( d.x > 0. ) ? index + 1 : index - 1;
              }
            }
          }else{
            break;
          }
        }
        if( index == -1 ){
          drawLine( fb,
                    p.x , p.y,
                    p.x + d.x * 300.,
                    p.y + d.y * 300.,
                    scale, FRAMEBUFFER::Color( 0.,1.,1.) );
        }
      }
      fb.normalize();
      fb.save_png(fn,1.,1.);
    }
    
  };

  namespace LOADER{
    namespace LENS {

      
      char *tokenize( char *s, char *token ) {
        char *p = (char *)s;
        while( *p && (*p == ' ' || *p == '\t') )
          p++;
        if( *p == 0 )
          return NULL;
        while(*p != ' ' && *p != '\t' && *p != '\r' && *p != '\n' && *p) {
          //printf("%c\n",*p);
          *token = *p;
          token++; p++;
        }
        *token = 0;
        if( *p == '\r' || *p == '\n' || *p == 0 )
          return NULL;
        return p;
      }
      
      void load( const char *filename, Lens& lens ){
        
        FILE *fp = fopen( filename, "rb" );
        if(!fp)return;
        double sumz = 0.;
        float R,curv,disz,nd,vd;
        float bf = 0.;
        float ir = 0.;//image sensor.
        while( !feof(fp) ){
          char line[1024],*p,token[256];
          
          fgets(line,1024,fp);
          if( feof( fp ) )
            break;
          
          p = line;
          p = tokenize( p, token );
          if( token[0] == 'S' || token[0] == 'I' || token[0] == 'Y' ){
            bool radiusInvert = false;
            Surface surface;
            char t;
            
            t = token[0];
            if( strcmp(token,"S") == 0 )
              surface.type_ = Surface::STANDARD;
            else if( strcmp(token,"SZ") == 0 ){
              surface.type_ = Surface::STANDARD;
              radiusInvert = true;
            }else if( strcmp(token,"Y") == 0 )
              surface.type_ = Surface::CYLINDER_Y;
            else if( strcmp(token,"I6") == 0 )
              surface.type_ = Surface::APERTURE_HEXAGONAL;
            else if( strcmp(token,"I") == 0 )
              surface.type_ = Surface::APERTURE_CIRCLE;
            
            p = tokenize( p, token );
            surface.diameter_  = atof( token ) / 2.f;
            p = tokenize( p, token );
            surface.radius_  = atof( token );
            if( radiusInvert && surface.radius_ != 0. )
              surface.radius_ = 1. / surface.radius_;
            
            surface.center_  = sumz + surface.radius_;
            p = tokenize( p, token ); // 面厚み.
            surface.thickness_ = atof( token );
            sumz += surface.thickness_;
            // surface.iris_x_ = 2.f;
            //  surface.iris_y_ = 1.f;
            if( surface.type_ == Surface::APERTURE_CIRCLE ||
                surface.type_ == Surface::APERTURE_HEXAGONAL ){
              // 楕円率もセットする?
              if( p ){
                p = tokenize( p, token );
                surface.iris_x_ = atof( token );
              }
              if( p ){
                p = tokenize( p, token );
                surface.iris_y_ = atof( token );
              }
            }else{
              // レンズなので、あっべ数等を設定する.
              p = tokenize( p, token );
              surface.ior_     = atof( token );
              surface.abbe_vd_ = 1.;
              if( p ){
                p = tokenize( p, token );
                surface.abbe_vd_ = atof( token );
              }
              if( p ){
                p = tokenize( p, token );
                surface.reflection_ = atof( token );
              }
              if( p ){
                p = tokenize( p, token );
                surface.roughness_ = atof( token );
              }
            }
            lens.surfaces_.push_back( surface ); // レンズフラッシュ
          }else if( token[0] == 'B' ){
            // backfocus.
            p = tokenize( p, token );
            bf = atof( token );
          }else if( token[0] == 'C' ){
            // imagesensor
            p = tokenize( p, token );
            ir = atof( token );
          }
        }
        lens.imageSurfaceZ_ = sumz + bf;
        lens.imageSurfaceR_ = ir;
        fclose( fp );
      }
    }
    namespace ZEMAX{
      
      char *tokenize( char *s, char *token ) {
        char *p = (char *)s;
        while( *p && (*p == ' ' || *p == '\t') )
          p++;
        if( *p == 0 )
          return NULL;
        while(*p != ' ' && *p != '\t' && *p != '\r' && *p != '\n' && *p) {
          //printf("%c\n",*p);
          *token = *p;
          token++; p++;
        }
        *token = 0;
        if( *p == '\r' || *p == '\n' || *p == 0 )
          return NULL;
        return p;
      }
      
      void load( const char *filename, Lens& lens ){
        
        FILE *fp = fopen( filename, "rb" );
        if(!fp)return;
        int surfaceIndex = -1;
        Surface surface;

        bool isAperture = false;
        double sumz = 0.;
        while( !feof(fp) ){
          char line[1024],*p,token[256];
          
          fgets(line,1024,fp);
          p = line;
          p = tokenize( p, token );

          if( strcmp( token, "SURF" ) == 0 ){
            if( surfaceIndex >= 0 ){
              if( surface.diameter_ > 0. ){
                surface.center_ = surface.center_ + surface.radius_; // 中心位置を調整しておく.
                surface.type_   = isAperture ? Surface::APERTURE_HEXAGONAL : Surface::STANDARD; // STOPは絞り.
                lens.surfaces_.push_back( surface ); // レンズフラッシュ
              }
            }
            p = tokenize( p, token );
            surfaceIndex = atoi( token );
            surface.init();
            isAperture = false;
          }
          
          if( strcmp( token, "TYPE" ) == 0 ){
            p = tokenize( p, token );
            if( strcmp( token, "STANDARD" ) == 0 )
              surface.type_ = Surface::STANDARD;
            else{
              printf("unknown surface: %s\n",token);
            }
          }
          if( strcmp( token, "STOP") == 0 )
            isAperture = true;
          
          if( strcmp( token, "CURV" ) == 0 ){
            p = tokenize( p, token );
            double curve = atof( token );
            if( curve != 0. )
              surface.radius_ = 1. / curve;
            else
              surface.radius_ = 0.;
          }
          
          if( strcmp( token, "DISZ" ) == 0 ){
            p = tokenize( p, token );
            double disz = atof( token );
            if( strcmp( token , "INFINITY" ) == 0 )
              disz = 0.;
            surface.thickness_= disz;
            surface.center_ = sumz / 2.;
            sumz += disz;
          }
          if( strcmp( token, "DIAM" ) == 0 ){
            p = tokenize( p, token );
            surface.diameter_ = atof( token );
          }
          
          if( strcmp( token, "GLAS" ) == 0 ){
            p = tokenize( p, token ); // name
            p = tokenize( p, token ); // nazo1
            p = tokenize( p, token ); // nazo2
            
            p = tokenize( p, token ); // ior
            surface.ior_ = atof( token );
            p = tokenize( p, token ); // abbe
            surface.abbe_vd_ = atof( token ); // fitting was done in inverted.
            if( surface.abbe_vd_ <= 0. )
              surface.abbe_vd_ =1.;
          }
        }
        // lens.surfaces_.push_back( surface ); // レンズフラッシュ
        lens.imageSurfaceZ_ = surface.center_;
        fclose(fp);
        
      }
      
    }
  }

  double focusDistance( Lens& lens, double fx ){
    Vector origin, dir;
    bool hit = false;
    double r = 0.01;
    double l = 500.;
    while(!hit){
      printf("%f,%f\n",r,l);
      hit = lens.traceRay( fx, 0., 0., r, 0.0, l, origin, dir );
      r = 0.1 * rand()/(double)RAND_MAX; // 通らないならパスを変えてみる.
    }
    double t = origin.y / (-dir.y);
    t = (t<0.)?FLT_MAX:t;
    if(!hit) return -1;
    return t;
  }
  double autofocus( Lens& lens, double dist ){
    Surface& last( lens.surfaces_.back() );
    double z = last.center_ - last.radius_;
    double lo = z;
    double hi = z+1000.;
    printf("lo %f,hi %f\n",lo,hi);
    int iter = 100;
    double fx,t ;
    do{
      fx = (lo+hi)/2.;
      t = focusDistance( lens, fx );
      if( t < 0. )
        break;
      printf("%f - %f : fx %f , target %f t %f\n",lo,hi,fx,dist,t);
      if( t > dist ){ // tを小さくしたい.
        lo = fx;
      }else{
        hi = fx;
      }
    }while(fabs(dist-t)>1e-8 && --iter);
    if( t < 0. )
      return 0.;
    return fx;
  }
}

#endif
