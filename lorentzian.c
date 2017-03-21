void lorentzian(float x,float a[],float *y,float dyda[],int na)
{
    /* a[1] = Lorentz lineshape halfwidth, a[2] = resonant frequency */

    float denom = (x - a[2])*(x - a[2]) + (a[1])*(a[1]);

    assert(denom);
    denom = 1.0/denom;
    *y = M_1_PI*a[1]*denom;
    if (dyda) {
        dyda[1] = *y*(1.0/a[1] - 2.0*a[1]*denom);
        dyda[2] = *y*2.0*(x - a[2])*denom;
        }
}
