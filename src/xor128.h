/*
 * copyright(c) 2015 Hajime UCHIMURA.
 * do not remove these text.
 * please contact before business use.
 */

#ifndef __XOR128_H
#define __XOR128_H

namespace RANDOM {

  // Xorshift
  // from Wikipedia
  class xor128{
  public:
    unsigned int x,y,z,w;
    xor128(){
      x = 123456789;
      y = 362436069;
      z = 521288629;
      w = 88675123;
    }
    inline unsigned int step(void) {
      unsigned int t;
      t = x ^ (x << 11);
      x = y; y = z; z = w;
      return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }
    void setSeed( unsigned u ){
      x ^= u;
    }
    inline double rand01() { return (double)step()/(UINT_MAX); }
  };

}

#endif
