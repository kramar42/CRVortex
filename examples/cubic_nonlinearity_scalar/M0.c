#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
    const int n = 1000;
    FILE  * fp;

    fp = fopen("sequence.txt", "w");

    int i;

    double delta = 10.0 / n;
    double lambda = 1.0;

    double a[n], b[n], c[n], diag[n];

    double Nn = 0.0, Np = 0.0;

    double d[n], psi[n], psip[n];
    double S = 0.0;

    int fl = 0;

    for (i = 0; i < n; i++)
    {
        psi[i] = pow(2.71, -i * delta * (i * delta));
    }

    for (i = 0; i < n; i++)
    {
        a[i] =  -1.0 / (delta * delta) + 1.0 / (2.0 * i * delta * delta);
        c[i] = lambda + 2 / (delta * delta); // + m * m / (i * i * delta * delta);
        b[i] =  -1.0 / (delta * delta) - 1.0 / (2.0 * i * delta * delta);
    }

    a[0] = 0.0;
    b[0] = -2 / (delta * delta);

    while(((S - 1) * (S - 1)) > 0.000001)
    {
        Np = 0.0;
        Nn = 0.0;
        fl++;

        for (i = 0; i < n; i++)
        {
            diag[i] = c[i];
            d[i] = psi[i] * psi[i] * psi[i];
            psip[i] = psi[i];
        }

        for (i = 1; i < n; i++)
        {
            diag[i] = diag[i] - b[i-1] * a[i] / diag[i-1];
            d[i] = d[i] - d[i-1] * a[i] / diag[i-1];
        }

        psi[n-1] = d[n-1] / diag[n-1];

        for (i = n-2; i>= 0; i--)
        {
            d[i] = d[i] - d[i+1] * b[i] / diag[i+1];
            psi[i] = d[i] / diag[i];
        }

        for (i = 0; i < n; i++)
        {
            Np += (psip[i] * psip[i] * i * delta * delta);
            Nn += (psi[i] * psi[i] * i * delta * delta);
        }

        S = pow((Np / Nn), 0.75);

        for (i = 0; i < n; i++)
        {
            psi[i] *= S;
        }

        printf("%lf\n", S);
    }

    Nn = 2 * 3.1415 * Nn;

    printf("%lf\n", Nn);

    for (i = 0; i < n; i++)
    {
        fprintf(fp, "%lf %lf\n", (i * delta), psi[i]);
    }

    fclose(fp);
    return 0;
}

