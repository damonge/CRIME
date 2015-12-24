#include "common.h"

Rng *init_rng(unsigned long seed)
{
  Rng *rng=my_malloc(sizeof(Rng));

  rng->mt[0]=seed & 0xffffffffUL;
  for (rng->mti=1;rng->mti<RNG_NRAN;rng->mti++) {
    rng->mt[rng->mti]=
      (1812433253UL*(rng->mt[rng->mti-1] ^ (rng->mt[rng->mti-1] >> 30))+rng->mti); 
    rng->mt[rng->mti]&=0xffffffffUL;
  }

  rng->calc_gauss=1;

  return rng;
}

void end_rng(Rng *rng)
{
  free(rng);
}

unsigned long rand_ulong(Rng *rng)
{
  unsigned long y;
  static unsigned long mag01[2]={0x0UL,RNG_MATRIX_A};

  if (rng->mti>=RNG_NRAN) {
    int kk;

    for (kk=0;kk<RNG_NRAN-RNG_MRAN;kk++) {
      y=(rng->mt[kk]&RNG_UPPER_MASK)|(rng->mt[kk+1]&RNG_LOWER_MASK);
      rng->mt[kk]=rng->mt[kk+RNG_MRAN] ^ (y>>1) ^ mag01[y & 0x1UL];
    }
    for(;kk<RNG_NRAN-1;kk++) {
      y=(rng->mt[kk]&RNG_UPPER_MASK)|(rng->mt[kk+1]&RNG_LOWER_MASK);
      rng->mt[kk]=rng->mt[kk+(RNG_MRAN-RNG_NRAN)] ^ (y>>1) ^ mag01[y & 0x1UL];
    }
    y=(rng->mt[RNG_NRAN-1]&RNG_UPPER_MASK)|(rng->mt[0]&RNG_LOWER_MASK);
    rng->mt[RNG_NRAN-1]=rng->mt[RNG_MRAN-1] ^ (y>>1) ^ mag01[y&0x1UL];

    rng->mti=0;
  }
  
  y=rng->mt[rng->mti++];

  y^=(y>>11);
  y^=(y<<7)&0x9d2c5680UL;
  y^=(y<<15)&0xefc60000UL;
  y^=(y>>18);

  return y;
}

double rand_real01(Rng *rng)
{
  return rand_ulong(rng)*(1.0/4294967296.0); 
}

double rand_gauss(Rng *rng)
{
  double x;

  if(rng->calc_gauss) {
    rng->phi=2*M_PI*rand_real01(rng);
    rng->u=sqrt(-2*log(1-rand_real01(rng)));
    x=rng->u*cos(rng->phi);
    rng->calc_gauss=0;
  }
  else {
    x=rng->u*sin(rng->phi);
    rng->calc_gauss=1;
  }

  return x;
}
