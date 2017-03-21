/*
 ** fit.c -- solution to Homework2
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "../gsl2nr/util.c"
#include "../gsl2nr/mrqmin.c"

#define MAX_N 100 /* maximum file size */

#define ACC 1.0e-3 /* iterative accuracy for chisq */

#define SQ(x) ((x)*(x))

void lorentz(float x,float a[],float *y,float dyda[],int na)
{
    /* a[1] = Lorentz lineshape halfwidth, a[2] = resonant frequency */

    float denom = SQ(x - a[2]) + SQ(a[1]);

    assert(denom);
    denom = 1.0/denom;
    *y = M_1_PI*a[1]*denom;
    if (dyda) {
        dyda[1] = *y*(1.0/a[1] - 2.0*a[1]*denom);
        dyda[2] = *y*2.0*(x - a[2])*denom;
        }
    }

void gauss(float x,float a[],float *y,float dyda[],int na)
{
    /* a[1] = Gaussian halfwidth, a[2] = center */

    const double fact = sqrt(M_LN2*M_1_PI);

    float ia,ia2;
    float nu2 = SQ(x - a[2]);

    assert(a[1]);
    ia = 1/a[1];
    ia2 = SQ(ia);

    *y = ia*fact*exp(-M_LN2*nu2*ia2);
    if (dyda) {
        dyda[1] = *y*ia*(2.0*M_LN2*nu2*ia2 - 1.0);
        dyda[2] = *y*2.0*M_LN2*(x - a[2])*ia2;
        }
    }

void fit(float nu[],float phi[],float sig[],int n,float a[],int ia[],int ma,
         void (*func)(float,float *,float *,float *,int),char *label)
{
    float **covar,**alpha,chisq,alamda,chiold,lamold;
    int i;

    /* allocate working matrices */

    covar = matrix(1,2,1,2);
    alpha = matrix(1,2,1,2);

    /* initialize */

    alamda = -1.0;
    chiold = HUGE_VAL;
    lamold = 0.001;

    /* iterate to best fit */

    do {
        mrqmin(nu-1,phi-1,sig-1,n,a-1,ia-1,ma,covar,alpha,&chisq,func,&alamda);
        (void) printf("%s:",label);
        for (i=0;i<ma;i++)
            (void) printf(" a[%i] = %g",i,a[i]);
        (void) printf(" chisq = %g alamda = %g\n",chisq,alamda);
        if (alamda < lamold && (chiold - chisq)/chisq < ACC)
            break;
        chiold = chisq;
        lamold = alamda;
        } while (1);

    /* one last call for error estimates */

    alamda = 0.0;
    mrqmin(NULL,NULL,NULL,0,a-1,ia-1,ma,covar,alpha,NULL,NULL,&alamda);
    (void) printf("Error estimates:");
    for (i=0;i<ma;i++)
        (void) printf(" a[%i]: %g",i,sqrt(covar[i+1][i+1]));
    (void) printf("\n");

    /* release resources */

    free_matrix(alpha,1,2,1,2);
    free_matrix(covar,1,2,1,2);
    }

int read_data(char *filename,float *x,float *y,float *e,int *n)
{
    FILE *fp;

    if ((fp = fopen(filename,"r")) == NULL) {
        (void) fprintf(stderr,"Unable to open \"%s\" for reading\n",filename);
        return 1;
        }

    *n = 0;
    while (*n < MAX_N && (fscanf(fp,"%f%f%f",&x[*n],&y[*n],&e[*n]) == 3))
        ++*n;

    (void) fclose(fp);

    (void) printf("%i data point%s read in%s\n",*n,*n==1?"":"s",
                  *n==MAX_N?" (arrays full)":"");

    return 0;
    }

int main(int argc,char *argv[])
{
    float nu[MAX_N],phi[MAX_N],sig[MAX_N],a[2];
    int ia[2] = {1,1};
    int n;

    if (argc != 2) {
        (void) fprintf(stderr,"Usage: %s dat-file\n",argv[0]);
        return 1;
        }

    if (read_data(argv[1],nu,phi,sig,&n))
        return 1;

    /* fit the Lorentzian */

    a[0] = 20.0; /* starting guess for alpha */
    a[1] = 50.0; /* starting guess for nu0 */
    (void) printf("\nLorentz definitions: a[0] = alpha, a[1] = nu0\n");
    fit(nu,phi,sig,n,a,ia,2,lorentz,"Lorentz");

    /* fit the Gaussian */

    a[0] = 20.0; /* starting guess for alpha */
    a[1] = 50.0; /* starting guess for nu0 */
    (void) printf("\nGauss definitions: a[0] = alpha, a[1] = nu0\n");
    fit(nu,phi,sig,n,a,ia,2,gauss,"Gauss");

    return 0;
    }
