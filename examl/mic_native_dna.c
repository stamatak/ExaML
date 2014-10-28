#include <immintrin.h>
#include <string.h>
#include <math.h>

#include "axml.h"
#include "mic_native.h"

static const int states = 4;
static const int statesSquare = 16;
static const int span = 4 * 4;
static const int maxStateValue = 16;

/* Common functions */

void updateModel_MIC(pInfo* part)
{
  double
    *EV               = part->EV,
    *tipVector        = part->tipVector,
    *aEV              = part->mic_EV,
    *aTipVector       = part->mic_tipVector;

  const int
    states = part->states,
    span = 4 * states,
    maxState = getUndetermined(part->dataType) + 1;

  int
    k, l;
  #pragma ivdep
  for (l = 0; l < 4 * states * states; ++l)
  {
    aEV[l] = EV[(l / span) * states + (l % states)];
  }

  for(int k = 0; k < maxState; k++)
  {
    #pragma ivdep
    for(int l = 0; l < states; l++)
    {
	aTipVector[k*span + l] = aTipVector[k*span + states + l] = aTipVector[k*span + 2*states + l] = aTipVector[k*span + 3*states + l] = tipVector[k*states + l];
    }
  }
}

/* DNA */

void makeP_DNA_MIC(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right,
               boolean saveMem, int maxCat)
{
  int
    i,
    j,
    k,
    span = states * numberOfCategories;

  /* assign some space for pre-computing and later re-using functions */

  double lz1[4] __attribute__((align(BYTE_ALIGNMENT)));
  double lz2[4] __attribute__((align(BYTE_ALIGNMENT)));
  double d1[4] __attribute__((align(BYTE_ALIGNMENT)));
  double d2[4] __attribute__((align(BYTE_ALIGNMENT)));


  /* multiply branch lengths with eigenvalues */
  for(i = 1; i < states; i++)
    {
      lz1[i] = EIGN[i] * z1;
      lz2[i] = EIGN[i] * z2;
    }


  /* loop over the number of rate categories, this will be 4 for the GAMMA model and
     variable for the CAT model */

  for(i = 0; i < numberOfCategories; i++)
    {
      /* exponentiate the rate multiplied by the branch */

      for(j = 1; j < states; j++)
	{
	  d1[j] = EXP(rptr[i] * lz1[j]);
	  d2[j] = EXP(rptr[i] * lz2[j]);
	}

      /* now fill the P matrices for the two branch length values */

      for(j = 0; j < states; j++)
	{
	  /* left and right are pre-allocated arrays */

	  left[i * states + j] = 1.0;
	  right[i * states + j] = 1.0;

	  for(k = 1; k < states; k++)
	    {
	      left[k * span + i * states + j]  = d1[k] * EI[states * j + k];
	      right[k * span + i * states + j] = d2[k] * EI[states * j + k];
	    }
	}
    }


  /* if memory saving is enabled and we are using CAT we need to do one additional P matrix
     calculation for a rate of 1.0 to compute the entries of a column/tree site comprising only gaps */


  if(saveMem)
    {
      i = maxCat;

      for(j = 1; j < states; j++)
	{
	  d1[j] = EXP (lz1[j]);
	  d2[j] = EXP (lz2[j]);
	}

      for(j = 0; j < states; j++)
	{
	  left[statesSquare * i  + states * j] = 1.0;
	  right[statesSquare * i + states * j] = 1.0;

	  for(k = 1; k < states; k++)
	    {
	      left[statesSquare * i + states * j + k]  = d1[k] * EI[states * j + k];
	      right[statesSquare * i + states * j + k] = d2[k] * EI[states * j + k];
	    }
	}
    }
}

void precomputeTips_DNA_MIC(int tipCase, double *tipVector, double *left, double *right,
                  double *umpLeft, double *umpRight,
                  int numberOfCategories)
{
  /* no precomputation needed if both children are inner nodes */
  if (tipCase == INNER_INNER)
    return;

  const int
    span 	= states * 4,
    umpSize 	= span * 16;

  for(int k = 0; k < umpSize; ++k)
  {
      umpLeft[k] = 0.0;
      umpRight[k] = 0.0;
  }

  for(int i = 0; i < maxStateValue; ++i)
  {
    for(int l = 0; l < states; ++l)
    {
	#pragma ivdep
	#pragma vector aligned
	for(int k = 0; k < span; ++k)
	{
	    umpLeft[span * i + k] +=  tipVector[i * states + l] * left[l * span + k];
	    if (tipCase == TIP_TIP)
	      umpRight[span * i + k] +=  tipVector[i * states + l] * right[l * span + k];
	}
    }
  }
}

inline void mic_fma4x16(const double* inv, double* outv, double* mulv)
{
    __mmask8 k1 = _mm512_int2mask(0x0F);
    __mmask8 k2 = _mm512_int2mask(0xF0);

    __m512d acc1 = _mm512_setzero_pd();
    __m512d acc2 = _mm512_setzero_pd();

    __m512d t;

    for(int k = 0; k < 4; k++)
    {
        t = _mm512_mask_extload_pd(t, k1, &inv[0 + k], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        t = _mm512_mask_extload_pd(t, k2, &inv[4 + k], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

        __m512d m = _mm512_load_pd(&mulv[k * 16]);
        acc1 = _mm512_fmadd_pd(t, m, acc1);

        t = _mm512_mask_extload_pd(t, k1, &inv[8 + k], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        t = _mm512_mask_extload_pd(t, k2, &inv[12 + k], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

        m = _mm512_load_pd(&mulv[k * 16 + 8]);
        acc2 = _mm512_fmadd_pd(t, m, acc2);
    }

    _mm512_store_pd(&outv[0], acc1);
    _mm512_store_pd(&outv[8], acc2);
}

void newviewGTRGAMMA_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement,
                  double *umpLeft, double *umpRight)
{
    __m512d minlikelihood_MIC = _mm512_set1_pd(minlikelihood);
    __m512d twotothe256_MIC = _mm512_set1_pd(twotothe256);
    __m512i absMask_MIC = _mm512_set1_epi64(0x7fffffffffffffffULL);

    int addScale = 0;

  /* we assume that P-matrix and eigenvectors are in correct layout already */
  double
    *aEV = extEV,
    *aRight = right,
    *aLeft = left,
    *umpX1 = umpLeft,
    *umpX2 = umpRight;
  
  switch(tipCase)
  {
    case TIP_TIP:
      {
	#pragma noprefetch umpX1,umpX2
	for (int i = 0; i < n; i++)
        {
            _mm_prefetch((const char *)&x3[span*(i+8)], _MM_HINT_ET1);
            _mm_prefetch((const char *)&x3[span*(i+8) + 8], _MM_HINT_ET1);

            _mm_prefetch((const char *)&x3[span*(i+1)], _MM_HINT_ET0);
            _mm_prefetch((const char *)&x3[span*(i+1) + 8], _MM_HINT_ET0);

            const double *uX1 = &umpX1[16 * tipX1[i]];
            const double *uX2 = &umpX2[16 * tipX2[i]];

            double uX[16] __attribute__((align(64)));

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < 16; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
            }

            double* v3 = &x3[i * 16];

            mic_fma4x16(uX, v3, aEV);
        } // sites loop
      }
      break;
    case TIP_INNER:
      {
        #pragma noprefetch umpX1
	for (int i = 0; i < n; i++)
        {
            _mm_prefetch((const char *)&x2[span*(i+16)], _MM_HINT_T1);
            _mm_prefetch((const char *)&x2[span*(i+16) + 8], _MM_HINT_T1);
            _mm_prefetch((const char *)&x3[span*(i+16)], _MM_HINT_ET1);
            _mm_prefetch((const char *)&x3[span*(i+16) + 8], _MM_HINT_ET1);

            _mm_prefetch((const char *)&x2[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char *)&x2[span*(i+1) + 8], _MM_HINT_T0);
            _mm_prefetch((const char *)&x3[span*(i+1)], _MM_HINT_ET0);
            _mm_prefetch((const char *)&x3[span*(i+1) + 8], _MM_HINT_ET0);

            /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */
            double* uX1 = &umpX1[span * tipX1[i]];
            double uX2[16] __attribute__((align(64)));
            double uX[16] __attribute__((align(64)));

            const double* v2 = &(x2[16 * i]);

            mic_fma4x16(v2, uX2, aRight);

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < 16; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
            }

            double* v3 = &(x3[span * i]);

            mic_fma4x16(uX, v3, aEV);

            __m512d t1 = _mm512_load_pd(&v3[0]);
	    t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax1 = _mm512_reduce_gmax_pd(t1);
            __m512d t2 = _mm512_load_pd(&v3[8]);
	    t2 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t2), absMask_MIC));
            double vmax2 = _mm512_reduce_gmax_pd(t2);

            if(vmax1 < minlikelihood && vmax2 < minlikelihood)
            {
	      /*	t1 = _mm512_mul_pd(t1, twotothe256_MIC);
        	_mm512_store_pd(&v3[0], t1);
        	t2 = _mm512_mul_pd(t2, twotothe256_MIC);
        	_mm512_store_pd(&v3[8], t2);*/
	     
#pragma vector aligned nontemporal
	      for(int l = 0; l < span; l++)
		v3[l] *= twotothe256;

                addScale += wgt[i];
            }
        } // site loop
      }
      break;
    case INNER_INNER:
    {
        for (int i = 0; i < n; i++)
        {
            _mm_prefetch((const char *) &x1[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x1[span*(i+8) + 8], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2[span*(i+8) + 8], _MM_HINT_T1);
            _mm_prefetch((const char *) &x3[span*(i+8)], _MM_HINT_ET1);
            _mm_prefetch((const char *) &x3[span*(i+8) + 8], _MM_HINT_ET1);

            _mm_prefetch((const char *) &x1[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x1[span*(i+1) + 8], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2[span*(i+1) + 8], _MM_HINT_T0);
            _mm_prefetch((const char *) &x3[span*(i+1)], _MM_HINT_ET0);
            _mm_prefetch((const char *) &x3[span*(i+1) + 8], _MM_HINT_ET0);

            double uX1[16] __attribute__((align(64)));
            double uX2[16] __attribute__((align(64)));
            double uX[16] __attribute__((align(64)));

            const double* v1 = &(x1[span * i]);
            const double* v2 = &(x2[span * i]);

            mic_fma4x16(v1, uX1, aLeft);
            mic_fma4x16(v2, uX2, aRight);

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < 16; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
            }

            double* v3 =  &(x3[span * i]);

            mic_fma4x16(uX, v3, aEV);

            __m512d t1 = _mm512_load_pd(&v3[0]);
	    t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax1 = _mm512_reduce_gmax_pd(t1);
            __m512d t2 = _mm512_load_pd(&v3[8]);
	    t2 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t2), absMask_MIC));
            double vmax2 = _mm512_reduce_gmax_pd(t2);

            if(vmax1 < minlikelihood && vmax2 < minlikelihood)
            {
	      /* t1 = _mm512_mul_pd(t1, twotothe256_MIC);
        	_mm512_store_pd(&v3[0], t1);
        	t2 = _mm512_mul_pd(t2, twotothe256_MIC);
        	_mm512_store_pd(&v3[8], t2);
	      */
	      
#pragma vector aligned nontemporal
	      for(int l = 0; l < span; l++)
		v3[l] *= twotothe256;

	      addScale += wgt[i];
            }
        }
    } break;
    default:
//      assert(0);
      break;
  }

  *scalerIncrement = addScale;

}

double evaluateGAMMA_MIC(int *wgt, double *x1_start, double *x2_start, double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable)
{
    double sum = 0.0;

    /* the left node is a tip */
    if(tipX1)
    {
	double
	  *aTipVec = tipVector;

        /* loop over the sites of this partition */
        for (int i = 0; i < n; i++)
        {
          /* access pre-computed tip vector values via a lookup table */
          const double *x1 = &(aTipVec[16 * tipX1[i]]);
          /* access the other(inner) node at the other end of the branch */
          const double *x2 = &(x2_start[span * i]);

          double term = 0.;

          #pragma ivdep
          #pragma vector aligned
          for(int j = 0; j < 16; j++)
              term += x1[j] * x2[j] * diagptable[j];

          term = log(0.25 * fabs(term));

          sum +=  wgt[i] * term;
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            _mm_prefetch((const char *) &x1_start[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x1_start[span*(i+8) + 8], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2_start[span*(i+8)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2_start[span*(i+8) + 8], _MM_HINT_T1);

            _mm_prefetch((const char *) &x1_start[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x1_start[span*(i+1) + 8], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2_start[span*(i+1)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2_start[span*(i+1) + 8], _MM_HINT_T0);

          const double *x1 = &(x1_start[span * i]);
          const double *x2 = &(x2_start[span * i]);

          double term = 0.;

          #pragma ivdep
          #pragma vector aligned
          for(int j = 0; j < 16; j++)
              term += x1[j] * x2[j]  * diagptable[j];

          term = log(0.25 * fabs(term));

          sum +=  wgt[i] * term;
        }
    }

    return sum;
}

void sumGAMMA_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
    const double
      *aTipVec = tipVector;

    switch(tipCase)
    {
      case TIP_TIP:
      {
        #pragma unroll(8)
        for(int i = 0; i < n; i++)
        {
            const double *left  = &(aTipVec[16 * tipX1[i]]);
            const double *right = &(aTipVec[16 * tipX2[i]]);

            #pragma ivdep
            #pragma vector aligned nontemporal
            for(int l = 0; l < 16; l++)
            {
                sumtable[i * span + l] = left[l] * right[l];
            }
        }
      } break;
      case TIP_INNER:
      {
	#pragma unroll(8)
	for(int i = 0; i < n; i++)
        {
          _mm_prefetch((const char *) &x2_start[span*(i+32)], _MM_HINT_T1);
          _mm_prefetch((const char *) &x2_start[span*(i+32) + 8], _MM_HINT_T1);

          _mm_prefetch((const char *) &x2_start[span*(i+8)], _MM_HINT_T0);
          _mm_prefetch((const char *) &x2_start[span*(i+8) + 8], _MM_HINT_T0);

          const double *left = &(aTipVec[16 * tipX1[i]]);
          const double *right = &(x2_start[span * i]);

          #pragma ivdep
          #pragma vector aligned nontemporal
          for(int l = 0; l < 16; l++)
          {
              sumtable[i * span + l] = left[l] * right[l];
          }
        }
      } break;
      case INNER_INNER:
      {
	#pragma unroll(8)
        for(int i = 0; i < n; i++)
        {
            _mm_prefetch((const char *) &x1_start[span*(i+16)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x1_start[span*(i+16) + 8], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2_start[span*(i+16)], _MM_HINT_T1);
            _mm_prefetch((const char *) &x2_start[span*(i+16) + 8], _MM_HINT_T1);

            _mm_prefetch((const char *) &x1_start[span*(i+4)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x1_start[span*(i+4) + 8], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2_start[span*(i+4)], _MM_HINT_T0);
            _mm_prefetch((const char *) &x2_start[span*(i+4) + 8], _MM_HINT_T0);

            const double *left  = &(x1_start[span * i]);
            const double *right = &(x2_start[span * i]);

            #pragma ivdep
            #pragma vector aligned nontemporal
            for(int l = 0; l < 16; l++)
            {
                sumtable[i * span + l] = left[l] * right[l];
            }
        }
      } break;
  //    default:
  //      assert(0);
    }
}

void coreGTRGAMMA_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wgt)
{
    double diagptable0[16] __attribute__((align(64)));
    double diagptable1[16] __attribute__((align(64)));
    double diagptable2[16] __attribute__((align(64)));
    double diagptable01[16] __attribute__((align(64)));
    double diagptable02[16] __attribute__((align(64)));

    /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */

    for(int i = 0; i < 4; i++)
    {
        const double ki = gammaRates[i];
        const double kisqr = ki * ki;

        diagptable0[i*4] = 1.;
        diagptable1[i*4] = 0.;
        diagptable2[i*4] = 0.;

        for(int l = 1; l < 4; l++)
        {
          diagptable0[i * 4 + l]  = exp(EIGN[l] * ki * lz);
          diagptable1[i * 4 + l] = EIGN[l] * ki;
          diagptable2[i * 4 + l] = EIGN[l] * EIGN[l] * kisqr;
        }
    }

    #pragma ivdep
    for(int i = 0; i < 16; i++)
    {
        diagptable01[i] = diagptable0[i] * diagptable1[i];
        diagptable02[i] = diagptable0[i] * diagptable2[i];
    }

    /* loop over sites in this partition */

    const int aligned_width = upper % 8 == 0 ? upper / 8 : upper / 8 + 1;

    double dlnLBuf[8] __attribute__((align(64)));
    double d2lnLBuf[8] __attribute__((align(64)));
    for (int j = 0; j < 8; ++j)
    {
        dlnLBuf[j] = 0.;
        d2lnLBuf[j] = 0.;
    }

    __mmask16 k1 = _mm512_int2mask(0x000000FF);

    for (int i = 0; i < aligned_width; i++)
    {
        _mm_prefetch((const char *) &sumtable[i * span * 8], _MM_HINT_T0);
        _mm_prefetch((const char *) &sumtable[i * span * 8 + 8], _MM_HINT_T0);

        /* access the array with pre-computed values */
        const double *sum = &sumtable[i * span * 8];

        /* initial per-site likelihood and 1st and 2nd derivatives */

        double invBuf[8] __attribute__((align(64)));
        double d1Buf[8] __attribute__((align(64)));
        double d2Buf[8] __attribute__((align(64)));

        __m512d invVec;
        __m512d d1Vec;
        __m512d d2Vec;
        int mask = 0x01;

        #pragma noprefetch sum
        #pragma unroll(8)
        for(int j = 0; j < 8; j++)
        {
            _mm_prefetch((const char *) &sum[span*(j+8)], _MM_HINT_T1);
            _mm_prefetch((const char *) &sum[span*(j+8) + 8], _MM_HINT_T1);

            _mm_prefetch((const char *) &sum[span*(j+1)], _MM_HINT_T0);
            _mm_prefetch((const char *) &sum[span*(j+1) + 8], _MM_HINT_T0);

            __m512d d0_1 = _mm512_load_pd(&diagptable0[0]);
            __m512d d0_2 = _mm512_load_pd(&diagptable0[8]);

            __m512d d01_1 = _mm512_load_pd(&diagptable01[0]);
            __m512d d01_2 = _mm512_load_pd(&diagptable01[8]);

            __m512d d02_1 = _mm512_load_pd(&diagptable02[0]);
            __m512d d02_2 = _mm512_load_pd(&diagptable02[8]);

            __m512d s_1 = _mm512_load_pd(&sum[j*16]);
            __m512d s_2 = _mm512_load_pd(&sum[j*16 + 8]);
            __m512d inv_1 = _mm512_mul_pd(d0_1, s_1);
            __m512d d1_1 = _mm512_mul_pd(d01_1, s_1);
            __m512d d2_1 = _mm512_mul_pd(d02_1, s_1);

            __m512d inv_2 = _mm512_fmadd_pd(d0_2, s_2, inv_1);
            __m512d d1_2 = _mm512_fmadd_pd(d01_2, s_2, d1_1);
            __m512d d2_2 = _mm512_fmadd_pd(d02_2, s_2, d2_1);

            __mmask8 k1 = _mm512_int2mask(mask);
            mask <<= 1;

            // reduce
            inv_2 = _mm512_add_pd (inv_2, _mm512_swizzle_pd(inv_2, _MM_SWIZ_REG_CDAB));
            inv_2 = _mm512_add_pd (inv_2, _mm512_swizzle_pd(inv_2, _MM_SWIZ_REG_BADC));
            inv_2 = _mm512_add_pd (inv_2, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(inv_2), _MM_PERM_BADC)));
            invVec = _mm512_mask_mov_pd(invVec, k1, inv_2);

            d1_2 = _mm512_add_pd (d1_2, _mm512_swizzle_pd(d1_2, _MM_SWIZ_REG_CDAB));
            d1_2 = _mm512_add_pd (d1_2, _mm512_swizzle_pd(d1_2, _MM_SWIZ_REG_BADC));
            d1_2 = _mm512_add_pd (d1_2, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d1_2), _MM_PERM_BADC)));
            d1Vec = _mm512_mask_mov_pd(d1Vec, k1, d1_2);

            d2_2 = _mm512_add_pd (d2_2, _mm512_swizzle_pd(d2_2, _MM_SWIZ_REG_CDAB));
            d2_2 = _mm512_add_pd (d2_2, _mm512_swizzle_pd(d2_2, _MM_SWIZ_REG_BADC));
            d2_2 = _mm512_add_pd (d2_2, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d2_2), _MM_PERM_BADC)));
            d2Vec = _mm512_mask_mov_pd(d2Vec, k1, d2_2);
        }

        _mm512_store_pd(&invBuf[0], invVec);
        _mm512_store_pd(&d1Buf[0], d1Vec);
        _mm512_store_pd(&d2Buf[0], d2Vec);

        #pragma ivdep
        #pragma vector aligned
        for (int j = 0; j < 8; ++j)
        {
            const double inv_Li = 1.0 / invBuf[j];

            const double d1 = d1Buf[j] * inv_Li;
            const double d2 = d2Buf[j] * inv_Li;

            dlnLBuf[j] += wgt[i * 8 + j] * d1;
            d2lnLBuf[j] += wgt[i * 8 + j] * (d2 - d1 * d1);
        }
    } // site loop

    double dlnLdlz = 0.;
    double d2lnLdlz2 = 0.;
    for (int j = 0; j < 8; ++j)
    {
        dlnLdlz += dlnLBuf[j];
        d2lnLdlz2 += d2lnLBuf[j];
    }

    *ext_dlnLdlz   = dlnLdlz;
    *ext_d2lnLdlz2 = d2lnLdlz2;
}
