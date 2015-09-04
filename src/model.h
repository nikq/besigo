/*
 *
 * MQO file loader.
 *
 */
#ifndef __MODEL_H
#define __MODEL_H

#include <stdio.h>
#include <string.h>
#include <vector>
#include "include/embree2/rtcore.h"
#include "include/embree2/rtcore_ray.h"
#include "color3.h"
#include "kelemenmlt.h"

namespace MODEL {
  
  struct Material{

    float metallic_;// 0 1 0
    float subsurface_;// 0 1 0
    float specular_;// 0 1 .5
    float roughness_;// 0 1 .5
    float specularTint_;// 0 1 0
    float anisotropic_;// 0 1 0
    float sheen_;// 0 1 0
    float sheenTint_;// 0 1 .5
    float clearcoat_;// 0 1 0
    float clearcoatGloss_;// 0 1 1

    COLOR3::Color3 color_; // diffuse.
    float diffuse_;    // 拡散反射率.
    float refraction_; // 透過率.
    float reflection_; // 反射率.
    float emitter_;    // 発光.
    float emitter_color_; // 発光色.(色温度)
    float ior_;
    float abbe_;

    Material(){
      color_.set(COLOR3::XYZ,0.950450003147125,1.,1.08891701698303); // D65 white.
      ior_  = 1.;
      abbe_ = 1.;
    }
    void normalize(void){
      // 正規化する.
      // mia_materialを参考するに, 透明度はdiffuseを奪う.
      // 反射率はdiffuseと透明度を奪う.
      printf("\nbefore    : D %f, R %f, T %f\n",diffuse_,reflection_,refraction_);
      if( (diffuse_ + refraction_) > 1.f )
        diffuse_ -= (refraction_ + diffuse_ - 1.f);
      if( (reflection_ + diffuse_ + refraction_) > 1.f ){
        float remain = reflection_ + refraction_ + diffuse_ - 1.f;
        diffuse_    -= remain / 2.f;
        refraction_ -= remain / 2.f;
      }
      // pdfからCDFにする.
      // diffuse, reflection, refraction, absorb.
      printf("normalized: D %f, R %f, T %f\n",diffuse_,reflection_,refraction_);
      diffuse_ += refraction_;
      reflection_ += diffuse_;
      //reflection_ += diffuse_;
      //refraction_ += reflection_; // ( = reflection_ + diffuse_ )
      printf("CDF       : D %f, R %f, T %f\n",diffuse_,reflection_,refraction_);
      
      
      // "water" shader(3) col(0.820 0.890 1.000 0.000) dif(0.010) amb(0.600) emi(0.000) spc(0.000) power(5.00) refract(1.330)
      // default values.
      metallic_ = 0;
      subsurface_ = 0;
      specular_ = 0.5;
      specularTint_ = 0.;
      anisotropic_ = 0;
      sheen_ = 0;
      sheenTint_ = 0.5;
      clearcoat_ = 0;
      clearcoatGloss_ = 1.;
    }
    double ior( double lambda ) const {
      // コーシーの式.
      double B = (ior_ - 1.) / abbe_ * 0.52345;
      double A = ior_ - B / 0.34522792;
      double C = lambda / 1000.;
      return A + B / (C*C);
    }

    double emitter( double lambda ) const {
      return SPECTRUM::normalizedPlanck( emitter_color_, lambda ) * emitter_;
    }
    
    const float PI = 3.14159265358979323846;
    inline float sqr(float x) { return x*x; }
    inline float clamp( float a, float l, float h ){
      return (a<l) ? l : ((a>h) ? h : a);
    }
    inline float mix( float a, float b, float c ){
      return a*(1-c)+b*c;
    }
    inline float max( float a, float b ){
      return (a>b) ? a : b;
    }
    inline float SchlickFresnel(float u) {
      float m = clamp(1-u, 0, 1);
      float m2 = m*m;
      return m2*m2*m; // pow(m,5)
    }
    inline float GTR1(float NdotH, float a){
      if (a >= 1) return 1/PI;
      float a2 = a*a;
      float t = 1 + (a2-1)*NdotH*NdotH;
      return (a2-1) / (PI*log(a2)*t);
    }
    inline float GTR2(float NdotH, float a){
      float a2 = a*a;
      float t = 1 + (a2-1)*NdotH*NdotH;
      return a2 / (PI * t*t);
    }
    inline float GTR2_aniso(float NdotH, float HdotX, float HdotY, float ax, float ay){
      return 1 / ( PI * ax*ay * sqr( sqr(HdotX/ax) + sqr(HdotY/ay) + NdotH*NdotH ));
    }
    inline float smithG_GGX(float Ndotv, float alphaG) {
      float a = alphaG*alphaG;
      float b = Ndotv*Ndotv;
      return 1/(Ndotv + sqrt(a + b - a*b));
    }
    float DisneyBRDF(
      VECTORMATH::Vector& L,
      VECTORMATH::Vector& V,
      VECTORMATH::Vector& N,
      VECTORMATH::Vector& X,
      VECTORMATH::Vector& Y,
      float Lambda ) {
      
      float NdotL = N.dot(L);
      float NdotV = N.dot(V);
      if (NdotL < 0 || NdotV < 0) return 0.;

      VECTORMATH::Vector H = L+V;
      H.normalize();
      float NdotH = N.dot(H);
      float LdotH = L.dot(H);

      float Cdlin = SPECTRUM::getReflectance( color_, Lambda );
      float Cdlum = Cdlin; // luminance approx.

      float Ctint = 1.; // normalize lum. to isolate hue+sat
      float Cspec0 = mix(specular_*.08*mix(1., Ctint, specularTint_), Cdlin, metallic_);
      float Csheen = mix(1, Ctint, sheenTint_);

      // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing
      // and mix in diffuse retro-reflection based on roughness
      float FL = SchlickFresnel(NdotL);
      float FV = SchlickFresnel(NdotV);
      float Fd90 = 0.5 + 2 * LdotH*LdotH * roughness_;
      float Fd = mix(1, Fd90, FL) * mix(1, Fd90, FV);

      // Based on Hanrahan-Krueger brdf approximation of isotropic bssrdf
      // 1.25 scale is used to (roughly) preserve albedo
      // Fss90 used to "flatten" retroreflection based on roughness
      float Fss90 = LdotH*LdotH*roughness_;
      float Fss = mix(1, Fss90, FL) * mix(1, Fss90, FV);
      float ss = 1.25 * (Fss * (1 / (NdotL + NdotV) - .5) + .5);

      // specular
      float aspect = sqrt(1-anisotropic_*.9);
      float ax = max(.001, sqr(roughness_)/aspect);
      float ay = max(.001, sqr(roughness_)*aspect);
      float Ds = GTR2_aniso(NdotH, H.dot(X), H.dot(Y), ax, ay);
      float FH = SchlickFresnel(LdotH);
      float Fs = mix(Cspec0, 1, FH);
      float roughg = sqr(roughness_*.5+.5);
      float Gs = smithG_GGX(NdotL, roughg) * smithG_GGX(NdotV, roughg);

      // sheen
      float Fsheen = FH * sheen_ * Csheen;

      // clearcoat (ior = 1.5 -> F0 = 0.04)
      float Dr = GTR1(NdotH, mix(.1,.001,clearcoatGloss_));
      float Fr = mix(.04, 1, FH);
      float Gr = smithG_GGX(NdotL, .25) * smithG_GGX(NdotV, .25);

      return ((1/PI) * mix(Fd, ss, subsurface_)*Cdlin + Fsheen)
        * (1-metallic_)
          + Gs*Fs*Ds + .25*clearcoat_*Gr*Fr*Dr;
    }

  };

  class MQO {
    
  public:
    
    class Vertex{
    public:
      Vertex(){;}
      Vertex( float _x, float _y, float _z ) : x(_x), y(_y), z(_z), a(0.f) {;}
      float x, y, z, a ;
    };
    
    typedef std::vector< Vertex       > Vertices;
    typedef std::vector< unsigned int > Indices;
    typedef std::vector< Material > Materials;
    typedef struct {
      std::vector< VECTORMATH::Vector > normal_;
      std::vector< int > material_;
      std::vector< int > index_;
    } Object;

    class Triangle{
    public:
      Triangle(){;}
      VECTORMATH::Vector v_[3];
    };
    typedef std::vector<Triangle> Triangles;
    
  public:

    MQO(){ ; }
    virtual ~MQO(){ ; }
    
    VECTORMATH::Vector cameraPos_;
    VECTORMATH::Vector cameraDirX_;
    VECTORMATH::Vector cameraDirY_;
    VECTORMATH::Vector cameraDirZ_;
    std::vector<Object> objects_;
    Materials           materials_;
    Triangles           emitters_; // for NEE
    double              cameraFocus_;
    
    void createCameraRay(
      VECTORMATH::Vector &cameraOutPos,
      VECTORMATH::Vector &cameraOutDir,
      VECTORMATH::Vector &worldRayPos,
      VECTORMATH::Vector &worldRayDir ){
      VECTORMATH::Vector cameraLocalPos = cameraOutPos - VECTORMATH::Vector(0.,0.,0.); // カメラローカルにおける座標.

      // レンズシミュレータは-xが出力方向なので、ベクトルの向きを読み替える.
      // レンズローカル座標は 右が+z 上が+y 奥が-x
      // ワールド座標は       右が+x 上が+y 奥が+z
      worldRayPos = cameraPos_ +
        cameraDirX_ * cameraOutPos.z +
        cameraDirY_ * cameraOutPos.y +
        cameraDirZ_ * -cameraOutPos.x;
      worldRayDir = 
        cameraDirX_ * cameraOutDir.z +
        cameraDirY_ * cameraOutDir.y +
        cameraDirZ_ * -cameraOutDir.x;
    }
    
    
    // 受け付ける文法："hagehoge"  / hoge( hoge hoge hoge )
    //トークンセパレータ：\x20 || \t
    char *tokenize( char *s, char *token ) {
      char *p = (char *)s;
      while( *p && (*p == ' ' || *p == '\t') )
        p++;
      if( *p == 0 )
        return NULL;
      if( *p == '\"' ) {
        p++;
        while(*p != '\"' && *p){
          *token = *p; token++; p++;
        }
        p++;
        *token = 0;
      } else {
        while(*p != ' ' && *p != '\t' && *p != '\r' && *p != '\n' && *p) {
          *token = *p;
          if(*p == '(') {
            token++; p++;
            while( *p != ')' && *p ){
              *token = *p; token++; p++;
            }
            *token = *p;
          }
          token++; p++;
        }
        *token = 0;
      }
      if( *p == '\r' || *p == '\n' || *p == 0 )
        return NULL;
      return p;
    }
    
    void calcSmoothNormals(
      Vertices& v,
      std::vector<VECTORMATH::Vector>& n,
      std::vector<int>& index ){
      for(int i=0;i<v.size();i++)
        n[i].set(0,0,0);
      for(int i=0;i<index.size();i+=3){
        VECTORMATH::Vector v1(
          v[ index[ i ] ].x,
          v[ index[ i ] ].y,
          v[ index[ i ] ].z );
        VECTORMATH::Vector v2(
          v[ index[ i +1] ].x,
          v[ index[ i +1] ].y,
          v[ index[ i +1] ].z );
        VECTORMATH::Vector v3(
          v[ index[ i +2] ].x,
          v[ index[ i +2] ].y,
          v[ index[ i +2] ].z );
        VECTORMATH::Vector e1 = v2 - v1;
        VECTORMATH::Vector e2 = v3 - v1;
        VECTORMATH::Vector nv = e2.cross( e1 );
        nv.normalize();
        n[index[i  ]] = n[index[i  ]] + nv;
        n[index[i+1]] = n[index[i+1]] + nv;
        n[index[i+2]] = n[index[i+2]] + nv;
      }
      for(int i=0;i<n.size();i++)
        n[i].normalize();
    }
    
    bool load( const char *fn, RTCScene& scene ) {
      
      FILE *fp;
      char line[1024];
      char *p,*p2;
      char token[256],t2[256];

      int nmat;
      int nvert;
      int nface;
      int i,a = 0,b = 0,c = 0,d = 0,m = 0;
      int depth;
      int totalface = 0, currface = 0;

      double x,y,z;

      fp = fopen(fn,"rb");
      if(fp == NULL)
        return false;

      fgets(line,1024,fp);
      if(strcmp(line,"Metasequoia Document\r\n") != 0)
        return false;

      fgets(line,1024,fp);
      //if(strcmp(line,"Format Text Ver 1.1\r\n") != 0) return false;
      
      materials_.clear();
      materials_.push_back( Material() );
      
      VECTORMATH::Vector lookat;
      
      //ロード開始
      while(1) {
        fgets(line,1024,fp);
        //printf("\"%s\"\n",line);
        if(strncmp(line,"Eof",3) == 0 || feof(fp))
          break;

        p = line;
        p = tokenize(p,token);

        if(strcmp(token,"Scene") == 0) {
          VECTORMATH::Vector pos, look, clook, cpos;
          while(1) {
            fgets(line,1024,fp);
            p = line;
            p = tokenize(p,token);
            if( token[0] == '}' ) {
              break;
            }
          }
        }

        if(strncmp(token,"Material",6) == 0) {
          materials_.clear();
          p = tokenize(p,token);
          nmat = atoi(token);
          printf("[MQO] %d materials\n",nmat);
          for(i=0;i<nmat;i++) {
            //マテリアルの読み込み
            fgets(line,1024,fp);
            p = line;
            p = tokenize(p,token);
            printf("[MAT] %2d:%12s :",i,token);
            //終了
            if(token[0] == '}')
              break;

            // 構文サンプル

            Material m;
            // "mat1"  shader(3) vcol(1) col(1.000 1.000 1.000 1.000) dif(1.000) amb(0.250) emi(0.000) spc(0.000) power(0.00) tex("00tex_master.BMP")
            // "floor" shader(3)         col(1.000 1.000 1.000 1.000) dif(0.800) amb(0.600) emi(0.000) spc(0.000) power(5.00)

            while( 1 ) {
              p = tokenize( p, token );
              if( strncmp( token, "col(" , 4) == 0 ){
                double r,g,b;
                p2 = tokenize(token+4,t2);
                r = atof( t2 );
                p2 = tokenize(p2,t2);
                g = atof( t2 );
                p2 = tokenize(p2,t2);
                b = atof( t2 );
                m.color_ = COLOR3::Color3( COLOR3::sRGB, r, g, b );
                // alpha goes to Refract
                p2 = tokenize(p2,t2);
                m.refraction_ = 1.0f - atof( t2 );
              }
              if( strncmp( token, "dif(" , 4) == 0 ){
                p2 = tokenize( token+4, t2 );
                m.diffuse_ = atof( t2 );
              }
              if( strncmp( token, "amb(" , 4) == 0 ){
                p2 = tokenize( token+4, t2 );
                m.roughness_ = atof( t2 ); // 読み替え.
              }
              if( strncmp( token, "emi(" , 4) == 0 ){
                p2 = tokenize( token+4, t2 );
                m.emitter_ = atof( t2 );
              }
              if( strncmp( token, "power(" , 6) == 0 ){
                p2 = tokenize( token+6, t2 );
                m.emitter_color_ = atof( t2 ) * 1000.;
              }
              if( strncmp( token, "spc(" , 4) == 0 ){
                p2 = tokenize( token+4, t2 );
                m.reflection_ = atof( t2 ); // 読み替え.
              }
              if( strncmp( token, "refract(" , 8) == 0 ){
                p2 = tokenize( token+8, t2 );
                m.ior_  = atof( t2 ); // 読み替え.
                m.abbe_ = 100 - (m.ior_ - 1.0) * 150 ; // 超適当にアッベ数を決める.
              }
              if( strncmp( token, "tex(" , 4) == 0 ){
                p2 = tokenize( token+4, t2 );
                // m.tex_.load( t2 );
              }
              if( !p )
                break;
            }
            m.normalize();
            printf("C(%5.3f,%5.3f,%5.3f),", m.color_.a_, m.color_.b_, m.color_.c_);
            printf("E %5.3f, D %5.3f,S %5.3f,R %5.3f IOR %f vd %f\n", m.emitter_, m.diffuse_, m.reflection_, m.refraction_ , m.ior_, m.abbe_ );
            materials_.push_back( m );
          }
        }

        if(strncmp(token,"Object",6) == 0) {

          int      materialIndex;
          Vertices vertex;
          Indices  index;
          std::vector<int> material;
          
          //object
          char objname[1024];
          p = tokenize(p,objname);
          printf("[OBJ] %-20s ",objname);
          depth = 1;
          currface = totalface;
          
          while(1) {
            fgets(line,1024,fp);
            p = line;
            p = tokenize(p,token);
            if(strcmp(token,"vertex")==0) {
              //vertex
              p = tokenize(p,token);
              nvert = atoi(token);
              printf("%5d,",nvert);
              depth++;
              for(i=0;i<nvert;i++) {
                fgets(line,1024,fp);
                p = line;
                p = tokenize(p,token);
                x = atof(token);
                p = tokenize(p,token);
                y = atof(token);
                p = tokenize(p,token);
                z = atof(token);
                //printf("%d : % 8.3f % 8.3f % 8.3f\n",i,x,y,z);
                vertex.push_back( Vertex( x,y,z ) );
              }
            } else if(strcmp(token,"face")==0) {
              //face
              p = tokenize(p,token);
              nface = atoi(token);
              depth++;
              printf("%5d,",nface);
              for(i=0;i<nface;i++) {
                fgets(line,1024,fp);
                p = line;
                p = tokenize( p, token );
                if( token[0] == '3' ) {
                  //p = tokenize( p, token );
                  int m, a, b, c;
                  while( 1 ) {
                    p = tokenize( p, token );
                    if( token[0] == 'V' ) {
                      p2 = tokenize( token+2, t2 );
                      a  = atoi( t2 );
                      p2 = tokenize( p2, t2 );
                      b  = atoi( t2 );
                      p2 = tokenize( p2, t2 );
                      c  = atoi( t2 );
                    }
                    if( token[0] == 'M' ) {
                      p2 = tokenize( token+2, t2 );
                      m  = atoi( t2 );
                    }
                    if( !strncmp(token,"UV(",3)  ){
                      float uv[6];
                      p2 = tokenize( token+3, t2 );
                      uv[0] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[1] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[2] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[3] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[4] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[5] = atof(t2);
                    }
                    if( !p )
                      break;
                  }
                  
                  materialIndex = m; // 代表.
                  index.push_back( a );
                  index.push_back( b );
                  index.push_back( c );
                  material.push_back(m);
                  
                  if( materials_[m].emitter_ > 0. ){
                    // エミッター登録.
                    Triangle t;
                    t.v_[0].set(vertex[a].x,vertex[a].y,vertex[a].z);
                    t.v_[1].set(vertex[b].x,vertex[b].y,vertex[b].z);
                    t.v_[2].set(vertex[c].x,vertex[c].y,vertex[c].z);
                    emitters_.push_back(t);
                  }
                  
                  totalface++;
                } else if(token[0] == '4') {
                  int m, a, b, c, d;
                  while( 1 ) {
                    p  = tokenize( p, token );
                    if( token[0] == 'V' ){
                      p2 = tokenize( token+2, t2 );
                      a  = atoi( t2 );
                      p2 = tokenize( p2, t2 );
                      b  = atoi( t2 );
                      p2 = tokenize( p2, t2 );
                      c  = atoi( t2 );
                      p2 = tokenize( p2, t2 );
                      d  = atoi( t2 );
                    }
                    if( token[0] == 'M' ){
                      p2 = tokenize( token+2, t2 );
                      m  = atoi( t2 );
                    }
                    if( !strncmp(token,"UV(",3)  ){
                      float uv[8];
                      p2 = tokenize( token+3, t2 );
                      uv[0] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[1] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[2] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[3] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[4] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[5] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[6] = atof(t2);
                      p2    = tokenize( p2, t2 );
                      uv[7] = atof(t2);
                    }
                    if( !p )
                      break;
                  }
                  
                  materialIndex = m;
                  index.push_back(a);
                  index.push_back(b);
                  index.push_back(c);
                  index.push_back(a);
                  index.push_back(c);
                  index.push_back(d);
                  material.push_back(m);
                  material.push_back(m);
                  if( materials_[m].emitter_ > 0. ){
                    // エミッター登録.
                    Triangle t1,t2;
                    t1.v_[0].set(vertex[a].x,vertex[a].y,vertex[a].z);
                    t1.v_[1].set(vertex[b].x,vertex[b].y,vertex[b].z);
                    t1.v_[2].set(vertex[c].x,vertex[c].y,vertex[c].z);
                    t2.v_[0].set(vertex[a].x,vertex[a].y,vertex[a].z);
                    t2.v_[1].set(vertex[c].x,vertex[c].y,vertex[c].z);
                    t2.v_[2].set(vertex[d].x,vertex[d].y,vertex[d].z);
                    emitters_.push_back(t1);
                    emitters_.push_back(t2);
                  }
                  
                  totalface++;
                  totalface++;
                }
              }
              
            } else if(strcmp(token,"BVertex")==0) {
              printf("バイナリ形式はサポートされていません.\n");
              break;
            } else if(token[0] == '}') {
              //オブジェクト終端のチェック
              depth--;
              if(depth == 0){
                printf("%6d triangles\n", totalface - currface);
                break;
              }
            }
          }

          if( strcmp( objname, "CAMERA_POS" ) == 0 ){
            // カメラの位置を固定するロケータ
            float cx = 0.f;
            float cy = 0.f;
            float cz = 0.f;
            for( int i=0;i<vertex.size();i++){
              cx += vertex[i].x;
              cy += vertex[i].y;
              cz += vertex[i].z;
            }
            cameraPos_.set( cx/(float)vertex.size(), cy/(float)vertex.size(), cz/(float)vertex.size() );
            printf("camera position locator found : %f,%f,%f\n",cameraPos_.x,cameraPos_.y,cameraPos_.z);
            unsigned geomID = rtcNewTriangleMesh( scene, RTC_GEOMETRY_STATIC, 0,0 );
            Object obj;
            objects_.push_back( obj );
          }else if( strcmp( objname, "CAMERA_DIR" )==0 ){
            // カメラの視線を固定するロケータ
            float cx = 0.f;
            float cy = 0.f;
            float cz = 0.f;
            for( int i=0;i<vertex.size();i++){
              cx += vertex[i].x;
              cy += vertex[i].y;
              cz += vertex[i].z;
            }
            lookat = VECTORMATH::Vector( cx/(float)vertex.size(), cy/(float)vertex.size(), cz/(float)vertex.size() );
            unsigned geomID = rtcNewTriangleMesh( scene, RTC_GEOMETRY_STATIC, 0,0 );
            Object obj;
            objects_.push_back( obj );
          }else{
            if( vertex.size() > 0 && index.size() > 0 ) {

              int faces = totalface - currface;
              // 登録すべきジオメトリがある場合はembreeに登録.
              unsigned geomID = rtcNewTriangleMesh( scene, RTC_GEOMETRY_STATIC, faces, vertex.size() );
              printf("regist %d faces, %d verts > geomID %d\n",faces,vertex.size(),geomID);

              // fill.
              {
                float* v = (float*)rtcMapBuffer( scene, geomID, RTC_VERTEX_BUFFER );
                for(int p = 0;p<vertex.size();p++){
                  v[p*4+0] = vertex[p].x;
                  v[p*4+1] = vertex[p].y;
                  v[p*4+2] = vertex[p].z;
                  v[p*4+3] = vertex[p].a;
                }
                rtcUnmapBuffer(scene,geomID,RTC_VERTEX_BUFFER);

                int  * i = (int  *)rtcMapBuffer( scene, geomID, RTC_INDEX_BUFFER );
                for(int p = 0;p<index.size();p++)
                  i[p] = index[p];
                rtcUnmapBuffer(scene,geomID,RTC_INDEX_BUFFER);
                Object obj;
                obj.index_.resize( index.size() );
                for(int p = 0;p<index.size();p++)
                  obj.index_[p] = index[p];
                obj.normal_.resize( vertex.size() );
                obj.material_.resize( index.size() / 3 );
                for(int p=0;p<index.size()/3;p++)
                  obj.material_[p] = material[p];
                calcSmoothNormals( vertex, obj.normal_, obj.index_ );
                objects_.push_back( obj );
              }
            }
          }
        }
      }
      {
        cameraDirZ_ = lookat - cameraPos_;
        cameraFocus_ = cameraDirZ_.length();
        cameraDirZ_.normalize();
        VECTORMATH::Vector up ( 0., 1., 0. );
        cameraDirX_ = cameraDirZ_.cross( up );
        cameraDirY_ = cameraDirX_.cross( cameraDirZ_ );
        printf("camera lookat locator found : %f,%f,%f\n",lookat.x, lookat.y, lookat.z );
        printf("camera pose X : %f,%f,%f\n", cameraDirX_.x, cameraDirX_.y, cameraDirX_.z );
        printf("camera pose Y : %f,%f,%f\n", cameraDirY_.x, cameraDirY_.y, cameraDirY_.z );
        printf("camera pose Z : %f,%f,%f\n", cameraDirZ_.x, cameraDirZ_.y, cameraDirZ_.z );
      }

      printf("[MQO] %d faces,Done.\n",totalface);
      fclose(fp);
      return true;
    }
    
    VECTORMATH::Vector normal( const RTCRay & ray ){
      //printf("geomID %d, primID %d\n",ray.geomID,ray.primID);
      int i1 = objects_[ ray.geomID ].index_[ ray.primID * 3     ];
      int i2 = objects_[ ray.geomID ].index_[ ray.primID * 3 + 1 ];
      int i3 = objects_[ ray.geomID ].index_[ ray.primID * 3 + 2 ];
      VECTORMATH::Vector n1 = objects_[ray.geomID].normal_[i1];
      VECTORMATH::Vector n2 = objects_[ray.geomID].normal_[i2];
      VECTORMATH::Vector n3 = objects_[ray.geomID].normal_[i3];
      VECTORMATH::Vector n = n1 * ( 1 - ray.u - ray.v ) + n2 * ray.u + n3 * ray.v;
      return n.normal();
    }
    Material* material( const RTCRay& ray ){
      return &(materials_[objects_[ ray.geomID ].material_[ ray.primID ]]);
    }

    VECTORMATH::Vector emitterDir(
      KelemenMLT& mlt,
      VECTORMATH::Vector& v1,
      VECTORMATH::Vector& v2,
      VECTORMATH::Vector& v3 ){
      int ri = emitters_.size() * mlt.NextSample();
      double u = mlt.NextSample();
      double v = mlt.NextSample();
      v1 = emitters_[ri].v_[0];
      v2 = emitters_[ri].v_[1];
      v3 = emitters_[ri].v_[2];
      return emitters_[ri].v_[0] * (1. - u - v ) + emitters_[ri].v_[1] * u + emitters_[ri].v_[2] * v;
    }
  };
}

#endif
