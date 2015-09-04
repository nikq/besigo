/*
 * mtseq.h
 *
 * implements Mersenne Twister.
 */


#ifndef __MT_SEQUENCE
#define __MT_SEQUENCE

namespace RANDOM {

  /* Period parameters */
#define MT_N   624
#define MT_M   397
#define MTMIX(u,v) (((u) & 0x80000000UL)|((v) & 0x7fffffffUL))
#define TWIST(u,v) ((MTMIX(u,v)>>1)^((v)&1UL?0x9908b0dfUL:0UL))

  class MTSequence {
  private:
    unsigned long *mt_state, *mt_next, mt_left;

  public:
    virtual ~MTSequence(){
      if(mt_state)
        delete[] mt_state;
      mt_state = NULL;
    }
    MTSequence(){
      mt_state = new unsigned long [MT_N];
      mt_left  = 1;          // initial state
      init_genrand( 12345 ); // initial seed
    }

    void init_genrand(unsigned long s) {
      int j;
      mt_state[0] = s & 0xffffffffUL;
      for (j=1; j<MT_N; j++) {
        mt_state[j] = (1812433253UL * (mt_state[j-1] ^ (mt_state[j-1] >> 30)) + j);
        mt_state[j] &= 0xffffffffUL;  /* for >32 bit machines */
      }
      mt_left  = 1;
    }
    void   next_state(void);
    long   genrand_int31(void);
    double genrand_real1(void);
    float  genrand_realF(void);
  };

  inline void MTSequence::next_state(void)
    {
      unsigned long *p = mt_state;
      int j;

      mt_left = MT_N;
      mt_next = mt_state;

      for (j=MT_N - MT_M+1; --j; p++)
        *p = p[MT_M] ^ TWIST(p[0], p[1]);

      for (j=MT_M; --j; p++)
        *p = p[MT_M - MT_N] ^ TWIST(p[0], p[1]);

      *p = p[MT_M - MT_N] ^ TWIST(p[0], mt_state[0]);
    }

  /* generates a random number on [0,0x7fffffff]-interval */
  inline long MTSequence::genrand_int31(void)
    {
      unsigned long y;

      if (--mt_left == 0)
        next_state();
      y = *mt_next++;

      /* Tempering */
      y ^= (y >> 11);
      y ^= (y << 7 ) & 0x9d2c5680UL;
      y ^= (y << 15) & 0xefc60000UL;
      y ^= (y >> 18);

      return (long)(y>>1);
    }

  /* generates a random number on [0,1]-real-interval */
  inline double MTSequence::genrand_real1(void)
    {
      unsigned long y;

      if (--mt_left == 0)
        next_state();
      y = *mt_next++;
      /* Tempering */
      y ^= (y >> 11);
      y ^= (y << 7 ) & 0x9d2c5680UL;
      y ^= (y << 15) & 0xefc60000UL;
      y ^= (y >> 18);
      return (double)y * (1.0/4294967295.0);
    }

  /* generates a random number on [0,1]-real-interval */
  inline float MTSequence::genrand_realF(void)
    {
      unsigned long y;

      if (--mt_left == 0)
        next_state();
      y = *mt_next++;
      /* Tempering */
      y ^= (y >> 11);
      y ^= (y << 7 ) & 0x9d2c5680UL;
      y ^= (y << 15) & 0xefc60000UL;
      y ^= (y >> 18);
      return (float)y / 4294967295.0f;
    }

#undef MT_N
#undef MT_M
#undef MTMIX
#undef TWIST

}

#endif
