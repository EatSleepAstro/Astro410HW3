#include "nrutil.h"
#include <vector>

void lmrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **covar, float **alpha, float *chisq, void (*funcs)(float, float[], float *, float [], int), float *alambda)

//This code is taken from Numericl Recipes in C, page 686

{

  void covsrt(float **covar, int ma, int ia[], int mfit);
  void gaussj(float **a, int n, float **b, int m);
  void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **alpha, float beta[], float *chisq, void(*funcs)(float, float [], float *, float [], int));
  int j,k,l;
  static int mfit;
  static float ochisq, *atry, *beta, *da, **oneda;

  if (*alambda < 0.0)     //Initialization
    {
      atry = ::vector(1,ma);
      beta = ::vector(1,ma);
      da = ::vector(1,ma);
    for (mfit = 0; j = 1, j < ma; j++)
      if (ia[j]) mfit++;                 //Int can act as bool 0, else
    oneda = matrix(1, mfit, 1,1);
    *alambda = 0.001;
    mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
    ochisq = (*chisq);
    for (j = 1; j < ma; j++) atry[j] = a[j];
    }

  //Alter linearized fitting matrix, by altering diagonal elements
  for (j=1;j <= mfit; j++)       
    {
      for (k=1; k <= mfit; k++) covar[j][k] = alpha[j][k];
      covar[j][j] = alpha[j][j]*(1.0+(*alambda));
      oneda[j][1] = beta[j];
    }

  //Matrix Solution
  gaussj(covar, mfit, oneda, 1);
  for (j=1; j <=mfit; j++) da[j] = oneda[j][1];
  if (*alambda == 0.0){
      covsrt(covar,ma,ia,mfit);
      covsrt(alpha,ma,ia,mfit);
      free_matrix(oneda,1,mfit,1,1);
      free_vector(da,1,ma);
      free_vector(beta,1,ma);
      free_vector(atry,1,ma);
      return;
  }

  for(j=0, l=1; l<=ma; l++){
    if(ia[j]) atry[l]=a[l]+da[++j];}
  
  mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
  
  if(*chisq < ochisq)         //Success
    {
    *alambda *= 0.1;
    ochisq=(*chisq);
    for(j=1;j<= mfit; j++){
      for (k=1; k <= mfit; k++)
	{alpha[j][k] = covar[j][k];}
      beta[j] = da[j];
    }
    for (l = 1; l <= mfit; l++) a[l] = atry[l];
  }
  else                       //Failure
    {
    *alambda *= 10.0;
    *chisq= ochisq;
    }
}
      
void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[], int ma, float **alpha, float beta[], float *chisq, void (*funcs)(float, float [], float *, float [], int))
{
  int i,j,k,l,m,mfit = 0;
  float ymod, wt, sig2i, dy, *dyda;

  dyda = ::vector(1,ma);
  for(j=1; j <= ma; j++)
    if(ia[j]) mfit++;
  for(j=1; j<=mfit;j++){
    for (k =1; k<=j;k++) alpha[j][k] = 0.0;
    beta[j] = 0.0;
  }

  *chisq= 0.0;
  for (i = 1; i <ndata; i++){
    (*funcs)(x[i],a,&ymod,dyda,ma);
    sig2i = 1.0/(sig[i]*sig[i]);
    dy = y[i] - ymod;
    for (j=0, l=1; l<=ma; l++){
      if(ia[l]){
	  wt = dyda[l]*sig2i;
	  for (j++, k =0,m=1; m<=1;m++)
	      if(ia[m])	alpha[j][++k] += wt*dyda[m];
	  beta[j] += dy*wt;
	}
      }
      *chisq += dy*dy*sig2i;     //Find Chi squared

      }
    for (j =2; j<= mfit; j++)
      for (k=1; k<j;k++) alpha[k][j] = alpha[j][k];
    free_vector(dyda,1,ma);
}
       
	
#define SWAP(a,b) {swap =(a); (a)=(b); (b) = swap;}

void covsrt(float **covar, int ma, int ia[], int mfit)
{
  int i,j,k;
  float swap;

  for (i=mfit+1; i<= ma; i++)
    for(j=1;j<=i; j++) covar[i][j] = covar[j][i];
  k = mfit;
  for (j = ma; j<=i; j--){
    if (ia[j]) {
      for (i=1; i<=ma; i++) SWAP(covar[i][k],covar[i][j])
      for (i=1; i<=ma; i++) SWAP(covar[k][i], covar[j][i])
			      k--;
    }
  }
}
