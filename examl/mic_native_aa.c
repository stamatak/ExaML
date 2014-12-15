#include <immintrin.h>
#include <string.h>
#include <math.h>

#include "axml.h"
#include "mic_native.h"

static const int states = 20;
static const int statesSquare = 20 * 20;
static const int span = 20 * 4;
static const int maxStateValue = 23;

void makeP_PROT_MIC(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right,
               boolean saveMem, int maxCat)
{
  int
    i,
    j,
    k,
    span = states * numberOfCategories;

  /* assign some space for pre-computing and later re-using functions */

  double lz1[20] __attribute__((align(BYTE_ALIGNMENT)));
  double lz2[20] __attribute__((align(BYTE_ALIGNMENT)));
  double d1[20] __attribute__((align(BYTE_ALIGNMENT)));
  double d2[20] __attribute__((align(BYTE_ALIGNMENT)));


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

void precomputeTips_PROT_MIC(int tipCase, double *tipVector, double *left, double *right,
                  double *umpLeft, double *umpRight,
                  int numberOfCategories)
{
  /* no precomputation needed if both children are inner nodes */
  if (tipCase == INNER_INNER)
    return;

  const int
    span 	= states * numberOfCategories,
    umpSize 	= span * maxStateValue;

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

inline void mic_fma4x80(const double* inv, double* outv, double* mulv)
{
    __mmask8 k1 = _mm512_int2mask(0x0F);
    __mmask8 k2 = _mm512_int2mask(0xF0);
    for(int l = 0; l < 80; l += 40)
    {
        __m512d t = _mm512_setzero_pd();

        t = _mm512_extload_pd(&inv[l], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        __m512d m = _mm512_load_pd(&mulv[l]);
        __m512d acc = _mm512_load_pd(&outv[l]);
        __m512d r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l], r);

        m = _mm512_load_pd(&mulv[l + 8]);
        acc = _mm512_load_pd(&outv[l + 8]);
        r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l + 8], r);

        t = _mm512_mask_extload_pd(t, k1, &inv[l], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        t = _mm512_mask_extload_pd(t, k2, &inv[l+20], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);

        m = _mm512_load_pd(&mulv[l + 16]);
        acc = _mm512_load_pd(&outv[l + 16]);
        r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l + 16], r);

        t = _mm512_extload_pd(&inv[l+20], _MM_UPCONV_PD_NONE, _MM_BROADCAST_1X8, _MM_HINT_NONE);
        m = _mm512_load_pd(&mulv[l + 24]);
        acc = _mm512_load_pd(&outv[l + 24]);
        r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l + 24], r);

        m = _mm512_load_pd(&mulv[l + 32]);
        acc = _mm512_load_pd(&outv[l + 32]);
        r = _mm512_fmadd_pd(t, m, acc);
        _mm512_store_pd(&outv[l + 32], r);
    }
}


void newviewGTRGAMMAPROT_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement,
                  double *umpLeft, double *umpRight)
{
  __m512d minlikelihood_MIC = _mm512_set1_pd(minlikelihood);
  __m512d twotothe256_MIC = _mm512_set1_pd(twotothe256);
  __m512i absMask_MIC = _mm512_set1_epi64(0x7fffffffffffffffULL);

  int addScale = 0;

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
        /* multiply all possible tip state vectors with the respective P-matrices
        */

        for (int i = 0; i < n; i++)
        {
            const double *uX1 = &umpX1[span * tipX1[i]];
            const double *uX2 = &umpX2[span * tipX2[i]];

            double uX[span] __attribute__((align(BYTE_ALIGNMENT)));
            double* v3 = &x3[i * span];

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
                for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

        } // sites loop
      }
      break;
    case TIP_INNER:
      {
        for (int i = 0; i < n; i++)
        {
            #pragma unroll(10)
            for (int j = 0; j < span; j += 8)
            {
                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T1);
//                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T0);
            }

            /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */
            double* uX1 = &umpX1[span * tipX1[i]];
            double uX2[span] __attribute__((align(BYTE_ALIGNMENT)));
            double uX[span] __attribute__((align(BYTE_ALIGNMENT)));

            double* v3 = &(x3[span * i]);

            const double* v2 = &(x2[span * i]);

            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX2[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
		#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aRight[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&v2[k], uX2, &aRight[k * span]);
            }

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
		#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax = _mm512_reduce_gmax_pd(t1);
            for (int l = 8; l < span; l += 8)
            {
                __m512d t = _mm512_load_pd(&v3[l]);
                t = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t), absMask_MIC));
                double vmax2 = _mm512_reduce_gmax_pd(t);
                vmax = MAX(vmax, vmax2);
            }

            if (vmax < minlikelihood)
            {
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
      /* same as above, without pre-computations */

        for (int i = 0; i < n; i++)
        {

            #pragma unroll(10)
            for (int j = 0; j < span; j += 8)
            {
                _mm_prefetch((const char *)&x1[span*(i+1) + j], _MM_HINT_T1);
                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T1);
//                _mm_prefetch((const char *)&x1[span*(i+1) + j], _MM_HINT_T0);
//                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T0);
            }


            double uX1[span] __attribute__((align(BYTE_ALIGNMENT)));
            double uX2[span] __attribute__((align(BYTE_ALIGNMENT)));
            double uX[span] __attribute__((align(BYTE_ALIGNMENT)));

            double* v3 = &(x3[span * i]);

            const double* v1 = &(x1[span * i]);
            const double* v2 = &(x2[span * i]);

            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX1[l] = 0.;
                uX2[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
		#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aRight[span*(k+1) + j], _MM_HINT_T0);
                    _mm_prefetch((const char *)&aLeft[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&v1[k], uX1, &aLeft[k * span]);
                mic_fma4x80(&v2[k], uX2, &aRight[k * span]);
            }

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
		#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax = _mm512_reduce_gmax_pd(t1);
            for (int l = 8; l < span; l += 8)
            {
                __m512d t = _mm512_load_pd(&v3[l]);
                t = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t), absMask_MIC));
                double vmax2 = _mm512_reduce_gmax_pd(t);
                vmax = MAX(vmax, vmax2);
            }

            if (vmax < minlikelihood)
            {
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



double evaluateGAMMAPROT_MIC(int *wgt, double *x1_start, double *x2_start, double *tipVector,
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
	    const double *x1 = &(aTipVec[span * tipX1[i]]);
	    /* access the other(inner) node at the other end of the branch */
	    const double *x2 = &(x2_start[span * i]);

	    #pragma unroll(10)
	    for (int k = 0; k < span; k += 8)
	    {
		    _mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
		    _mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
	    }

	    double term = 0.;

	    #pragma ivdep
	    #pragma vector aligned
	    #pragma noprefetch x2
	    for(int j = 0; j < span; j++) {
	      term += x1[j] * x2[j] * diagptable[j];
	    }

	    term = log(0.25 * fabs(term));

	    sum += wgt[i] * term;
	  }
      }
    else
      {
	for (int i = 0; i < n; i++)
	  {
	    #pragma unroll(10)
	    for (int k = 0; k < span; k += 8)
	    {
	      _mm_prefetch((const char *) &x1_start[span*(i+2) + k], _MM_HINT_T1);
	      _mm_prefetch((const char *) &x1_start[span*(i+1) + k], _MM_HINT_T0);

	      _mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
	      _mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
	    }

	    const double *x1 = &(x1_start[span * i]);
	    const double *x2 = &(x2_start[span * i]);

	    double term = 0.;

	    #pragma ivdep
	    #pragma vector aligned
	    #pragma noprefetch x1 x2
	    for(int j = 0; j < span; j++)
	      term += x1[j] * x2[j] * diagptable[j];

	    term = log(0.25 * fabs(term));

	    sum += wgt[i] * term;
	  }
      }

    return sum;
}

void sumGAMMAPROT_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
    double
      *aTipVec = tipVector;

    switch(tipCase)
    {
      case TIP_TIP:
      {
        for(int i = 0; i < n; i++)
	  {
	    const double *left  = &(aTipVec[span * tipX1[i]]);
	    const double *right = &(aTipVec[span * tipX2[i]]);

	    #pragma ivdep
	    #pragma vector aligned nontemporal
	    for(int l = 0; l < span; l++)
	      {
		sumtable[i * span + l] = left[l] * right[l];
	      }
	  }
      } break;
      case TIP_INNER:
      {
        for(int i = 0; i < n; i++)
        {
	  #pragma unroll(10)
	  for (int k = 0; k < span; k += 8)
	    {
	      _mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
	      _mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
	    }

          const double *left = &(aTipVec[span * tipX1[i]]);
          const double *right = &(x2_start[span * i]);

          #pragma ivdep
          #pragma vector aligned nontemporal
		  #pragma noprefetch right
          for(int l = 0; l < span; l++)
          {
              sumtable[i * span + l] = left[l] * right[l];
          }
        }
      } break;
      case INNER_INNER:
      {
        for(int i = 0; i < n; i++)
        {
	    #pragma unroll(10)
	    for (int k = 0; k < span; k += 8)
	      {
		_mm_prefetch((const char *) &x1_start[span*(i+2) + k], _MM_HINT_T1);
		_mm_prefetch((const char *) &x1_start[span*(i+1) + k], _MM_HINT_T0);

		_mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
		_mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
	      }

            const double *left  = &(x1_start[span * i]);
            const double *right = &(x2_start[span * i]);

            #pragma ivdep
            #pragma vector aligned nontemporal
	    #pragma noprefetch left right
            for(int l = 0; l < span; l++)
	      {
		sumtable[i * span + l] = left[l] * right[l];
	      }
        }
      } break;
  //    default:
  //      assert(0);
    }
}

void coreGTRGAMMAPROT_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wgt)
{
    static const int states = 20;
    static const int span = 20 * 4;

    double diagptable0[span] __attribute__((align(BYTE_ALIGNMENT)));
    double diagptable1[span] __attribute__((align(BYTE_ALIGNMENT)));
    double diagptable2[span] __attribute__((align(BYTE_ALIGNMENT)));
    double diagptable01[span] __attribute__((align(BYTE_ALIGNMENT)));
    double diagptable02[span] __attribute__((align(BYTE_ALIGNMENT)));

    /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */

    for(int i = 0; i < 4; i++)
    {
        const double ki = gammaRates[i];
        const double kisqr = ki * ki;

        diagptable0[i*states] = 1.;
        diagptable1[i*states] = 0.;
        diagptable2[i*states] = 0.;

        for(int l = 1; l < states; l++)
        {
          diagptable0[i * states + l]  = exp(EIGN[l] * ki * lz);
          diagptable1[i * states + l] = EIGN[l] * ki;
          diagptable2[i * states + l] = EIGN[l] * EIGN[l] * kisqr;
        }
    }

    #pragma ivdep
    for(int i = 0; i < span; i++)
    {
        diagptable01[i] = diagptable0[i] * diagptable1[i];
        diagptable02[i] = diagptable0[i] * diagptable2[i];
    }

    /* loop over sites in this partition */

    const int aligned_width = upper % 8 == 0 ? upper / 8 : upper / 8 + 1;

    double dlnLBuf[8] __attribute__((align(BYTE_ALIGNMENT)));
    double d2lnLBuf[8] __attribute__((align(BYTE_ALIGNMENT)));
    for (int j = 0; j < 8; ++j)
    {
        dlnLBuf[j] = 0.;
        d2lnLBuf[j] = 0.;
    }

    __mmask16 k1 = _mm512_int2mask(0x000000FF);

    for (int i = 0; i < aligned_width; i++)
    {
        /* access the array with pre-computed values */
        const double *sum = &sumtable[i * span * 8];

        /* initial per-site likelihood and 1st and 2nd derivatives */

        double invBuf[8] __attribute__((align(BYTE_ALIGNMENT)));
        double d1Buf[8] __attribute__((align(BYTE_ALIGNMENT)));
        double d2Buf[8] __attribute__((align(BYTE_ALIGNMENT)));

        __m512d invVec;
        __m512d d1Vec;
        __m512d d2Vec;
        int mask = 0x01;

        #pragma noprefetch sum
        #pragma unroll(8)
        for(int j = 0; j < 8; j++)
        {
	    #pragma unroll(10)
	    for (int k = 0; k < span; k += 8)
	      {
		_mm_prefetch((const char *) &sum[span*(j+2) + k], _MM_HINT_T1);
		_mm_prefetch((const char *) &sum[span*(j+1) + k], _MM_HINT_T0);
	      }

            __m512d inv_1 = _mm512_setzero_pd();
            __m512d d1_1 = _mm512_setzero_pd();
            __m512d d2_1 = _mm512_setzero_pd();

            for (int offset = 0; offset < span; offset += 8)
	      {
		__m512d d0_1 = _mm512_load_pd(&diagptable0[offset]);
		__m512d d01_1 = _mm512_load_pd(&diagptable01[offset]);
		__m512d d02_1 = _mm512_load_pd(&diagptable02[offset]);
		__m512d s_1 = _mm512_load_pd(&sum[j*span + offset]);

		inv_1 = _mm512_fmadd_pd(d0_1, s_1, inv_1);
		d1_1 = _mm512_fmadd_pd(d01_1, s_1, d1_1);
		d2_1 = _mm512_fmadd_pd(d02_1, s_1, d2_1);
	      }

            __mmask8 k1 = _mm512_int2mask(mask);
            mask <<= 1;

            // reduce
            inv_1 = _mm512_add_pd (inv_1, _mm512_swizzle_pd(inv_1, _MM_SWIZ_REG_CDAB));
            inv_1 = _mm512_add_pd (inv_1, _mm512_swizzle_pd(inv_1, _MM_SWIZ_REG_BADC));
            inv_1 = _mm512_add_pd (inv_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(inv_1), _MM_PERM_BADC)));
            invVec = _mm512_mask_mov_pd(invVec, k1, inv_1);

            d1_1 = _mm512_add_pd (d1_1, _mm512_swizzle_pd(d1_1, _MM_SWIZ_REG_CDAB));
            d1_1 = _mm512_add_pd (d1_1, _mm512_swizzle_pd(d1_1, _MM_SWIZ_REG_BADC));
            d1_1 = _mm512_add_pd (d1_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d1_1), _MM_PERM_BADC)));
            d1Vec = _mm512_mask_mov_pd(d1Vec, k1, d1_1);

            d2_1 = _mm512_add_pd (d2_1, _mm512_swizzle_pd(d2_1, _MM_SWIZ_REG_CDAB));
            d2_1 = _mm512_add_pd (d2_1, _mm512_swizzle_pd(d2_1, _MM_SWIZ_REG_BADC));
            d2_1 = _mm512_add_pd (d2_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d2_1), _MM_PERM_BADC)));
            d2Vec = _mm512_mask_mov_pd(d2Vec, k1, d2_1);
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


/****
 *       PROTEIN - LG4
 */
void updateModel_LG4_MIC(pInfo* part)
{
  double
    **EV              = part->EV_LG4,
    **tipVector       = part->tipVector_LG4,
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
      aEV[l] = EV[(l % span) / states][(l / span) * states + (l % states)];
    }

  for(int k = 0; k < maxState; k++)
    {
      for(int j = 0; j < 4; j++)
      {
	for(int l = 0; l < states; l++)
	  {
	    aTipVector[k*span + j*states + l] = tipVector[j][k*states + l];
	  }
      }
    }
}

void makeP_PROT_LG4_MIC(double z1, double z2, double *rptr, double *EI[4],  double *EIGN[4], int numberOfCategories, double *left, double *right)
{
  int
    i,
    j,
    k;

  double
    d1[64],
    d2[64];

  for(i = 0; i < numberOfCategories; i++)
    {
      for(j = 1; j < states; j++)
	{
	  d1[j] = EXP (rptr[i] * EIGN[i][j] * z1);
	  d2[j] = EXP (rptr[i] * EIGN[i][j] * z2);
	}

      for(j = 0; j < states; j++)
	{
	  left[i * states + j] = 1.0;
	  right[i * states + j] = 1.0;

	  for(k = 1; k < states; k++)
	    {
	      left[k * span + i * states + j]  = d1[k] * EI[i][states * j + k];
	      right[k * span + i * states + j] = d2[k] * EI[i][states * j + k];
	    }
	}
    }
}

void precomputeTips_PROT_LG4_MIC(int tipCase, double *tipVector[4], double *left, double *right,
                  double *umpLeft, double *umpRight,
                  int numberOfCategories)
{
  /* no precomputation needed if both children are inner nodes */
  if (tipCase == INNER_INNER)
    return;

  const int
    span 	= states * numberOfCategories,
    umpSize 	= span * maxStateValue;

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
	    umpLeft[span * i + k] +=  tipVector[k/20][i * states + l] * left[l * span + k];
	    if (tipCase == TIP_TIP)
	      umpRight[span * i + k] +=  tipVector[k/20][i * states + l] * right[l * span + k];
	}
    }
  }
}

void newviewGTRGAMMAPROT_LG4_MIC(int tipCase,
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
        for (int i = 0; i < n; i++)
        {
            const double *uX1 = &umpX1[span * tipX1[i]];
            const double *uX2 = &umpX2[span * tipX2[i]];

            double uX[span] __attribute__((align(BYTE_ALIGNMENT)));
            double* v3 = &x3[i * span];

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
                for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

        } // sites loop
      }
      break;
    case TIP_INNER:
      {
        for (int i = 0; i < n; i++)
        {
            #pragma unroll(10)
            for (int j = 0; j < span; j += 8)
            {
                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T1);
//                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T0);
            }

            /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */
            double* uX1 = &umpX1[span * tipX1[i]];
            double uX2[span] __attribute__((align(BYTE_ALIGNMENT)));
            double uX[span] __attribute__((align(BYTE_ALIGNMENT)));

            double* v3 = &(x3[span * i]);

            const double* v2 = &(x2[span * i]);

            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX2[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
				#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aRight[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&v2[k], uX2, &aRight[k * span]);
            }

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
		#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax = _mm512_reduce_gmax_pd(t1);
            for (int l = 8; l < span; l += 8)
            {
                __m512d t = _mm512_load_pd(&v3[l]);
                t = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t), absMask_MIC));
                double vmax2 = _mm512_reduce_gmax_pd(t);
                vmax = MAX(vmax, vmax2);
            }

            if (vmax < minlikelihood)
            {
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

            #pragma unroll(10)
            for (int j = 0; j < span; j += 8)
            {
                _mm_prefetch((const char *)&x1[span*(i+1) + j], _MM_HINT_T1);
                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T1);
//                _mm_prefetch((const char *)&x1[span*(i+1) + j], _MM_HINT_T0);
//                _mm_prefetch((const char *)&x2[span*(i+1) + j], _MM_HINT_T0);
            }


            double uX1[span] __attribute__((align(BYTE_ALIGNMENT)));
            double uX2[span] __attribute__((align(BYTE_ALIGNMENT)));
            double uX[span] __attribute__((align(BYTE_ALIGNMENT)));

            double* v3 = &(x3[span * i]);

            const double* v1 = &(x1[span * i]);
            const double* v2 = &(x2[span * i]);

            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX1[l] = 0.;
                uX2[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
		#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aRight[span*(k+1) + j], _MM_HINT_T0);
                    _mm_prefetch((const char *)&aLeft[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&v1[k], uX1, &aLeft[k * span]);
                mic_fma4x80(&v2[k], uX2, &aRight[k * span]);
            }

            #pragma ivdep
            #pragma vector aligned
            for(int l = 0; l < span; ++l)
            {
                uX[l] = uX1[l] * uX2[l];
                v3[l] = 0.;
            }

            for(int k = 0; k < states; ++k)
            {
		#pragma unroll(10)
            	for (int j = 0; j < span; j += 8)
                {
                    _mm_prefetch((const char *)&aEV[span*(k+1) + j], _MM_HINT_T0);
                }

                mic_fma4x80(&uX[k], v3, &aEV[k * span]);
            }

            __m512d t1 = _mm512_load_pd(&v3[0]);
            t1 = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t1), absMask_MIC));
            double vmax = _mm512_reduce_gmax_pd(t1);
            for (int l = 8; l < span; l += 8)
            {
                __m512d t = _mm512_load_pd(&v3[l]);
                t = _mm512_castsi512_pd(_mm512_and_epi64(_mm512_castpd_si512(t), absMask_MIC));
                double vmax2 = _mm512_reduce_gmax_pd(t);
                vmax = MAX(vmax, vmax2);
            }

            if (vmax < minlikelihood)
            {
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



double evaluateGAMMAPROT_LG4_MIC(int *wgt, double *x1_start, double *x2_start, double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable, double *weights)
{
    double wtable[span] __attribute__((align(BYTE_ALIGNMENT)));

    /* pre-multiply diagptable entries with the corresponding weights */
    for(int j = 0; j < 4; j++)
      for(int k = 0; k < states; k++)
	{
	  wtable[j * states + k] = diagptable[j * states + k] * weights[j];
	}

    double sum = 0.0;

    /* the left node is a tip */
    if(tipX1)
    {
        /* loop over the sites of this partition */
        for (int i = 0; i < n; i++)
        {
          const double
    	    *aTipVec = tipVector;

	  /* access pre-computed tip vector values via a lookup table */
	  const double *x1 = &(aTipVec[span * tipX1[i]]);
	  /* access the other(inner) node at the other end of the branch */
	  const double *x2 = &(x2_start[span * i]);

	  #pragma unroll(10)
	  for (int k = 0; k < span; k += 8)
	    {
	      _mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
	      _mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
	    }

	  double term = 0.;

	  #pragma ivdep
	  #pragma vector aligned
	  #pragma noprefetch x2
	  for(int j = 0; j < span; j++)
	    {
	      term += x1[j] * x2[j] * wtable[j];
	    }

	  term = log(fabs(term));

	  sum += wgt[i] * term;
        }
    }
    else
    {
      for (int i = 0; i < n; i++)
	{
	  #pragma unroll(10)
	  for (int k = 0; k < span; k += 8)
	    {
	      _mm_prefetch((const char *) &x1_start[span*(i+2) + k], _MM_HINT_T1);
	      _mm_prefetch((const char *) &x1_start[span*(i+1) + k], _MM_HINT_T0);

	      _mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
	      _mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
	    }

	  const double *x1 = &(x1_start[span * i]);
	  const double *x2 = &(x2_start[span * i]);

	  double term = 0.;

	  #pragma ivdep
	  #pragma vector aligned
	  #pragma noprefetch x1 x2
	  for(int j = 0; j < span; j++)
	    term += x1[j] * x2[j] * wtable[j];

	  term = log(fabs(term));

	  sum += wgt[i] * term;
	}
    }

    return sum;
}

void sumGAMMAPROT_LG4_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n)
{
    const double
      *aTipVec = tipVector;

    switch(tipCase)
    {
      case TIP_TIP:
      {
        for(int i = 0; i < n; i++)
        {
            const double *left  = &(aTipVec[span * tipX1[i]]);
            const double *right = &(aTipVec[span * tipX2[i]]);

            #pragma ivdep
            #pragma vector aligned nontemporal
            for(int l = 0; l < span; l++)
            {
                sumtable[i * span + l] = left[l] * right[l];
            }
        }
      } break;
      case TIP_INNER:
      {
        for(int i = 0; i < n; i++)
        {
	  #pragma unroll(10)
	  for (int k = 0; k < span; k += 8)
	    {
	      _mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
	      _mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
	    }

          const double *left = &(aTipVec[span * tipX1[i]]);
          const double *right = &(x2_start[span * i]);

          #pragma ivdep
          #pragma vector aligned nontemporal
	  #pragma noprefetch right
          for(int l = 0; l < span; l++)
          {
              sumtable[i * span + l] = left[l] * right[l];
          }
        }
      } break;
      case INNER_INNER:
      {
        for(int i = 0; i < n; i++)
	  {
	      #pragma unroll(10)
	      for (int k = 0; k < span; k += 8)
	      {
		_mm_prefetch((const char *) &x1_start[span*(i+2) + k], _MM_HINT_T1);
		_mm_prefetch((const char *) &x1_start[span*(i+1) + k], _MM_HINT_T0);

		_mm_prefetch((const char *) &x2_start[span*(i+2) + k], _MM_HINT_T1);
		_mm_prefetch((const char *) &x2_start[span*(i+1) + k], _MM_HINT_T0);
	      }

	      const double *left  = &(x1_start[span * i]);
	      const double *right = &(x2_start[span * i]);

	      #pragma ivdep
	      #pragma vector aligned nontemporal
	      #pragma noprefetch left right
	      for(int l = 0; l < span; l++)
	      {
		  sumtable[i * span + l] = left[l] * right[l];
	      }
	  }
      } break;
  //    default:
  //      assert(0);
    }
}

void coreGTRGAMMAPROT_LG4_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN[4], double *gammaRates,
    double lz, int *wgt, double *weights)
{
    double diagptable0[span] __attribute__((align(BYTE_ALIGNMENT)));
    double diagptable1[span] __attribute__((align(BYTE_ALIGNMENT)));
    double diagptable2[span] __attribute__((align(BYTE_ALIGNMENT)));
    double diagptable01[span] __attribute__((align(BYTE_ALIGNMENT)));
    double diagptable02[span] __attribute__((align(BYTE_ALIGNMENT)));

    /* pre-compute the derivatives of the P matrix for all discrete GAMMA rates */

    for(int i = 0; i < 4; i++)
    {
        const double ki = gammaRates[i];
        const double kisqr = ki * ki;

        diagptable0[i*states] = 1. * weights[i];
        diagptable1[i*states] = 0.;
        diagptable2[i*states] = 0.;

        for(int l = 1; l < states; l++)
        {
          diagptable0[i * states + l]  = exp(EIGN[i][l] * ki * lz) * weights[i];
          diagptable1[i * states + l] = EIGN[i][l] * ki;
          diagptable2[i * states + l] = EIGN[i][l] * EIGN[i][l] * kisqr;
        }
    }

    #pragma ivdep
    for(int i = 0; i < span; i++)
    {
        diagptable01[i] = diagptable0[i] * diagptable1[i];
        diagptable02[i] = diagptable0[i] * diagptable2[i];
    }

    /* loop over sites in this partition */

    const int aligned_width = upper % 8 == 0 ? upper / 8 : upper / 8 + 1;

    double dlnLdlz = 0.;
    double d2lnLdlz2 = 0.;

    __mmask16 k1 = _mm512_int2mask(0x000000FF);

    for (int i = 0; i < aligned_width; i++)
    {
        /* access the array with pre-computed values */
        const double *sum = &sumtable[i * span * 8];

        /* initial per-site likelihood and 1st and 2nd derivatives */

        double invBuf[8] __attribute__((align(BYTE_ALIGNMENT)));
        double d1Buf[8] __attribute__((align(BYTE_ALIGNMENT)));
        double d2Buf[8] __attribute__((align(BYTE_ALIGNMENT)));

        __m512d invVec;
        __m512d d1Vec;
        __m512d d2Vec;
        int mask = 0x01;

        #pragma noprefetch sum
        #pragma unroll(8)
        for(int j = 0; j < 8; j++)
        {
	    #pragma unroll(10)
	    for (int k = 0; k < span; k += 8)
	    {
		    _mm_prefetch((const char *) &sum[span*(j+2) + k], _MM_HINT_T1);
		    _mm_prefetch((const char *) &sum[span*(j+1) + k], _MM_HINT_T0);
	    }

            __m512d inv_1 = _mm512_setzero_pd();
            __m512d d1_1 = _mm512_setzero_pd();
            __m512d d2_1 = _mm512_setzero_pd();

            for (int offset = 0; offset < span; offset += 8)
            {
                __m512d d0_1 = _mm512_load_pd(&diagptable0[offset]);
                __m512d d01_1 = _mm512_load_pd(&diagptable01[offset]);
                __m512d d02_1 = _mm512_load_pd(&diagptable02[offset]);
                __m512d s_1 = _mm512_load_pd(&sum[j*span + offset]);

                inv_1 = _mm512_fmadd_pd(d0_1, s_1, inv_1);
                d1_1 = _mm512_fmadd_pd(d01_1, s_1, d1_1);
                d2_1 = _mm512_fmadd_pd(d02_1, s_1, d2_1);
            }

            __mmask8 k1 = _mm512_int2mask(mask);
            mask <<= 1;

            // reduce
            inv_1 = _mm512_add_pd (inv_1, _mm512_swizzle_pd(inv_1, _MM_SWIZ_REG_CDAB));
            inv_1 = _mm512_add_pd (inv_1, _mm512_swizzle_pd(inv_1, _MM_SWIZ_REG_BADC));
            inv_1 = _mm512_add_pd (inv_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(inv_1), _MM_PERM_BADC)));
            invVec = _mm512_mask_mov_pd(invVec, k1, inv_1);

            d1_1 = _mm512_add_pd (d1_1, _mm512_swizzle_pd(d1_1, _MM_SWIZ_REG_CDAB));
            d1_1 = _mm512_add_pd (d1_1, _mm512_swizzle_pd(d1_1, _MM_SWIZ_REG_BADC));
            d1_1 = _mm512_add_pd (d1_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d1_1), _MM_PERM_BADC)));
            d1Vec = _mm512_mask_mov_pd(d1Vec, k1, d1_1);

            d2_1 = _mm512_add_pd (d2_1, _mm512_swizzle_pd(d2_1, _MM_SWIZ_REG_CDAB));
            d2_1 = _mm512_add_pd (d2_1, _mm512_swizzle_pd(d2_1, _MM_SWIZ_REG_BADC));
            d2_1 = _mm512_add_pd (d2_1, _mm512_castsi512_pd(_mm512_permute4f128_epi32(_mm512_castpd_si512(d2_1), _MM_PERM_BADC)));
            d2Vec = _mm512_mask_mov_pd(d2Vec, k1, d2_1);
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

            dlnLdlz += wgt[i * 8 + j] * d1;
            d2lnLdlz2 += wgt[i * 8 + j] * (d2 - d1 * d1);
        }
    } // site loop

    *ext_dlnLdlz   = dlnLdlz;
    *ext_d2lnLdlz2 = d2lnLdlz2;
}

