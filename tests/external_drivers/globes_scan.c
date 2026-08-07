/* GLoBES on the constant-density energy scan of the comparison figure.
 * P(nu_mu -> nu_e) at L = 1300 km through 3 g/cm^3, 150 energies from
 * 0.6 to 20 GeV, standard settings.  Density carries the 0.992093 mass
 * defect and no electron fraction, GLoBES applying its own 0.5. */
#include <stdio.h>
#include <math.h>
#include <globes/globes.h>
int main(int argc, char *argv[])
{
    glb_params p;
    int i;
    glbInit(argv[0]);
    p = glbAllocParams();
    glbDefineParams(p, asin(sqrt(0.310)), asin(sqrt(2.240e-2)),
                    asin(sqrt(0.582)), 217.0*M_PI/180.0, 7.39e-5, 2.525e-3);
    glbSetDensityParams(p, 1.0, GLB_ALL);
    glbSetOscillationParameters(p);
    for (i = 0; i < 150; i++)
    {
        double E = pow(10.0, log10(0.6) + (log10(20.0)-log10(0.6))*i/149.0);
        printf("%.10g %.15g\n", E,
               glbConstantDensityProbability(2, 1, +1, E, 1300.0,
                                             3.0*0.992093));
    }
    glbFreeParams(p);
    return 0;
}
