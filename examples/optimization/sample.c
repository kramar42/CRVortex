#include <immintrin.h>

void find_stat(__m256**, __m256**, int, double, double);

int main(void)
{
	/* __m128i a = _mm_set_epi32(4, 3, 2, 1); */


    int i;
    int deg = 10 - 3;
    int cnum = 1 << deg;
    /* double **RePsi, **ImPsi; */
    double size = 60;

    double dx = size / cnum;

	__m256** RePsi, **ImPsi;
    RePsi = (__m256 **) malloc (cnum * sizeof(__m256 *));
    ImPsi = (__m256 **) malloc (cnum * sizeof(__m256 *));

    for (i = 0; i < cnum; i++)
    {
        RePsi[i] = (__m256 *) malloc (cnum * sizeof(__m256));
        ImPsi[i] = (__m256 *) malloc (cnum * sizeof(__m256));
    }

    find_stat(RePsi, ImPsi, cnum, dx, size);

    for (i = 0; i < cnum; i++)
    {
        free(RePsi[i]);
        free(ImPsi[i]);
    }

    free(RePsi);
    free(ImPsi);

	return 0;
}

void find_stat(__m256** Rdat, __m256** Idat, int cnum, double dx, double size)
{
    const int n = 10 * cnum;
	const

    int i, j;

    float _delta = size / n / 2;
	__m256 delta = _mm256_broadcast_ss(&_delta);

    double lambda = 0.85;

    __m256 a[n], b[n], c[n], diag[n];

    double Nn = 0.0, Np = 0.0;

    double d[n], psip[n];
	__m256 psi[n];
    double S = 0.0;

    int fl = 0;

	float _incr[8] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f};
	float _eights[8] = {8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};

	__m256 incr = _mm256_load_ps(_incr);
	__m256 eights = _mm256_load_ps(_eights);
	__m256 counter = _mm256_xor_ps(incr, incr);
	counter = _mm256_add_ps(counter, incr);

    for (i = 0; i < n; i++)
    {
		psi[i] = _mm256_mul_ps(delta, _mm256_mul_ps(counter, delta));
		counter = _mm256_add_ps(counter, incr);
        /* psi[i] = pow(2.71, -i * delta * (i * delta)); */
    }

    for (i = 0; i < n; i++)
    {
		a[i] = _mm256_mul_ps(delta, delta);
        /* a[i] = -1.0 / (delta * delta) + 1.0 / (2.0 * i * delta * delta); */
        /* c[i] = lambda + 2 / (delta * delta); */
        /* b[i] = -1.0 / (delta * delta) - 1.0 / (2.0 * i * delta * delta); */
    }

    /* a[0] = 0.0; */
    /* b[0] = -2 / (delta * delta); */

    /* while (((S - 1) * (S - 1)) > 0.000000001) */
    /* { */
    /*     Np = 0.0; */
    /*     Nn = 0.0; */
    /*     fl++; */

    /*     for (i = 0; i < n; i++) */
    /*     { */
    /*         diag[i] = c[i]; */

	/* d[i] = -2 * psi[i] * psi[i] * psi[i] + 4 * psi[i] * psi[i] * psi[i] * psi[i] * psi[i]; */

    /*         psip[i] = psi[i]; */
    /*     } */

    /*     for (i = 1; i < n; i++) */
    /*     { */
    /*         diag[i]=diag[i]-b[i-1]*a[i]/diag[i-1]; */
    /*         d[i]=d[i]-d[i-1]*a[i]/diag[i-1]; */
    /*     } */

    /*     psi[n-1]=d[n-1]/diag[n-1]; */

    /*     for (i=n-2; i>=0; i--) */
    /*     { */
    /*         d[i]=d[i]-d[i+1]*b[i]/diag[i+1]; */
    /*         psi[i]=d[i]/diag[i]; */
    /*     } */

    /*     for (i=0; i<n; i++) */
    /*     { */
    /*         Np+=(psip[i]*psip[i]*i*delta*delta); */
    /*         Nn+=(psi[i]*psi[i]*i*delta*delta); */
    /*     } */

    /*     S=pow((Np/Nn),0.6); */

    /*     for (i=0; i<n; i++) */
    /*     { */
    /*         psi[i]*=S; */
    /*     } */
    /* } */

	/*
    printf("%lf\n",psi[0]);

    int r;
    double buf;

    for (i=0;i<cnum;i++)
	{
        for (j=0;j<cnum;j++)
        {
            buf=sqrt((i-cnum/2)*(i-cnum/2)+(j-cnum/2)*(j-cnum/2)*1.0)*20;

            r=int_double(buf);

            if ((i==cnum/2) && (j==cnum/2))
			{
                printf("%d\n",r);
			}
            if ((r*delta*2)>size)
			{
                Rdat[i][j]=0.0;
			}
            else
			{
                Rdat[i][j]=psi[r];
			}
            Idat[i][j]=Rdat[i][j];
        }
	}
	*/
}
