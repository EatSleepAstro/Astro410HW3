void gauss(float x,float a[],float *y,float dyda[],int na)
{
    /* a[1] = Gaussian halfwidth, a[2] = center */

    const double fact = sqrt(M_LN2*M_1_PI);

    float ia,ia2;
    float nu2 = (x - a[2])*(x - a[2]);

    assert(a[1]);
    ia = 1/a[1];
    ia2 = (ia)*(ia);

    *y = ia*fact*exp(-M_LN2*nu2*ia2);
    if (dyda) {
        dyda[1] = *y*ia*(2.0*M_LN2*nu2*ia2 - 1.0);
        dyda[2] = *y*2.0*M_LN2*(x - a[2])*ia2;
        }
    }
