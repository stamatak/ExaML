#ifndef MIC_NATIVE_H_
#define MIC_NATIVE_H_

//#define VECTOR_PADDING 8
//#define GET_PADDED_WIDTH(w) w % VECTOR_PADDING == 0 ? w : w + (VECTOR_PADDING - (w % VECTOR_PADDING))

// general functions
void updateModel_MIC(pInfo* part);

// DNA data

void makeP_DNA_MIC(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories,
               double *left, double *right, boolean saveMem, int maxCat);

void precomputeTips_DNA_MIC(int tipCase, double *tipVector, double *left, double *right,
                  double *umpLeft, double *umpRight,
                  int numberOfCategories);

void newviewGTRGAMMA_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement,
                  double *umpLeft, double *umpRight);

double evaluateGAMMA_MIC(int *wptr,
                 double *x1_start, double *x2_start,
                 double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable);

void sumGAMMA_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

void coreGTRGAMMA_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr);

void sumcoreGTRGAMMA_MIC(int tipCase, double *x1_start, double *x2_start, double *tipVector,
        unsigned char *tipX1, unsigned char *tipX2, const int n,
        volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wgt);

// protein data - single matrix

void makeP_PROT_MIC(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories,
               double *left, double *right, boolean saveMem, int maxCat);

void precomputeTips_PROT_MIC(int tipCase, double *tipVector, double *left, double *right,
                  double *umpLeft, double *umpRight,
                  int numberOfCategories);

void newviewGTRGAMMAPROT_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement,
                  double *umpLeft, double *umpRight);

double evaluateGAMMAPROT_MIC(int *wptr,
                 double *x1_start, double *x2_start,
                 double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable);

void sumGAMMAPROT_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

void coreGTRGAMMAPROT_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr);


// protein data - LG4

void updateModel_LG4_MIC(pInfo* part);

void makeP_PROT_LG4_MIC(double z1, double z2, double *rptr, double *EI[4],  double *EIGN[4], int numberOfCategories, double *left, double *right);

void precomputeTips_PROT_LG4_MIC(int tipCase, double *tipVector[4], double *left, double *right,
                  double *umpLeft, double *umpRight,
                  int numberOfCategories);

void newviewGTRGAMMAPROT_LG4_MIC(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement,
                  double *umpLeft, double *umpRight);

double evaluateGAMMAPROT_LG4_MIC(int *wptr,
                 double *x1_start, double *x2_start,
                 double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable);

void sumGAMMAPROT_LG4_MIC(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

void coreGTRGAMMAPROT_LG4_MIC(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN[4], double *gammaRates, double lz, int *wrptr);



#endif /* MIC_NATIVE_H_ */
