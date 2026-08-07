// Prob3++ on the same scan.  SetMNS takes Dm2_32; SetDensityConversion
// carries 0.5 times the 0.9924922 ratio to this library's constant.
#include <cstdio>
#include <cmath>
#include "BargerPropagator.h"
int main()
{
    BargerPropagator b;
    b.SetWarningSuppression(true);
    b.SetDensityConversion(0.5*1.514423e-4/1.52588e-4);
    for (int i = 0; i < 150; i++)
    {
        double E = pow(10.0, log10(0.6) + (log10(20.0)-log10(0.6))*i/149.0);
        b.SetMNS(0.310, 2.240e-2, 0.582, 7.39e-5, 2.525e-3-7.39e-5,
                 217.0*M_PI/180.0, E, true, 1);
        b.propagateLinear(1, 1299.999945, 3.0);
        printf("%.10g %.15g\n", E, b.GetProb(2, 1));
    }
    return 0;
}
