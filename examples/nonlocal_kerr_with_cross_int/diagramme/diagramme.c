#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
    const int n = 1000;

    FILE *fp_1;
    FILE *fp_2;
    FILE *fp_t1;
    FILE *fp_t2;
    FILE *fp_dev;
    FILE *fp_dev2;
    FILE *fp_diag;

    fp_1=fopen("sequence1.txt","w");
    fp_2=fopen("sequence2.txt","w");
    fp_t1=fopen("sequencet1.txt","w");
    fp_t2=fopen("sequencet2.txt","w");
    fp_dev=fopen("sequencedev.txt","w");
    fp_dev2=fopen("sequencedev2.txt","w");
    fp_diag=fopen("diagramme.txt","w");

    int i;

    double delta=20.0/n;
    double lambda_1=0.5;
    double lambda_2=0.5;
    double alfa=1.0;
    double sigma=0.75;

    double m_1=1.0;
    double m_2=1.0;

    double bot_1[n],mid_1[n],top_1[n],diag_1[n];
    double bot_2[n],mid_2[n],top_2[n],diag_2[n];
    double bot_3[n],mid_3[n],top_3[n],diag_3[n];
    double bot_4[n],mid_4[n],top_4[n],diag_4[n];

    double norm_1,norm_p_1;
    double norm_2=1.0,norm_p_2=2.0;
    double norm_3,norm_p_3;

    double col_1[n], psi_1[n],psi_p_1[n];
    double col_2[n], psi_2[n],psi_p_2[n];
    double col_3[n], theta1[n], theta_p1[n];
    double col_4[n], theta2[n];

    double S1=0.0,S2=0.0;

    int fl=0;

    double prev=2.0;

    for(sigma=0.999; sigma<1; sigma+=0.0001)
    {
        prev=2.0;
        lambda_1=1.0;
        lambda_2=1.0;

        do
        {
            prev=norm_2;
            S1=0.0;
            S2=0.0;

            for (i=0; i<n; i++)
            {
                psi_1[i]=pow((i*delta),m_1)*pow(2.71, -i*delta*(i*delta));
                psi_2[i]=pow((i*delta),m_2)*pow(2.71, -i*delta*(i*delta));
                theta1[i]=pow(2.71, (1-i*delta*i*delta));
                theta2[i]=pow(2.71, (1-i*delta*i*delta));
            }

            for (i=1; i<n; i++)
            {
                bot_1[i]=-1.0/(delta*delta)+1.0/(2.0*i*delta*delta);
                mid_1[i]=lambda_1+2/(delta*delta)+m_1*m_1/(i*i*delta*delta);
                top_1[i]=-1.0/(delta*delta)-1.0/(2.0*i*delta*delta);

                bot_2[i]=-1.0/(delta*delta)+1.0/(2.0*i*delta*delta);
                mid_2[i]=lambda_2+2/(delta*delta)+m_2*m_2/(i*i*delta*delta);
                top_2[i]=-1.0/(delta*delta)-1.0/(2.0*i*delta*delta);

                bot_3[i]=-1.0/(delta*delta)+1.0/(2.0*i*delta*delta);
                top_3[i]=-1.0/(delta*delta)-1.0/(2.0*i*delta*delta);
                mid_3[i]=alfa+2/(delta*delta);

                bot_4[i]=-1.0/(delta*delta)+1.0/(2.0*i*delta*delta);
                top_4[i]=-1.0/(delta*delta)-1.0/(2.0*i*delta*delta);
                mid_4[i]=alfa+2/(delta*delta);
            }

            bot_1[0]=0.0;
            top_1[0]=-2/(delta*delta);

            bot_2[0]=0.0;
            top_2[0]=-2/(delta*delta);

            bot_3[0]=0.0;
            top_3[0]=-2/(delta*delta);
            mid_3[0]=alfa+2/(delta*delta);

            bot_4[0]=0.0;
            top_4[0]=-2/(delta*delta);
            mid_4[0]=alfa+2/(delta*delta);

            while(((S1-1)*(S1-1))>0.00000000002)
            {
                norm_p_1=0.0;
                norm_1=0.0;
                norm_p_2=0.0;
                norm_2=0.0;
                norm_p_3=0.0;
                norm_3=0.0;

                fl++;

                for (i=0; i<n; i++)
                {
                    diag_1[i]=mid_1[i];
                    diag_2[i]=mid_2[i];

                    col_1[i]=psi_1[i]*theta1[i];
                    col_2[i]=psi_2[i]*theta2[i];

                    psi_p_1[i]=psi_1[i];
                    psi_p_2[i]=psi_2[i];
                }

                for (i=2; i<n; i++)
                {
                    diag_1[i]=diag_1[i]-top_1[i-1]*bot_1[i]/diag_1[i-1];
                    col_1[i]=col_1[i]-col_1[i-1]*bot_1[i]/diag_1[i-1];

                    diag_2[i]=diag_2[i]-top_2[i-1]*bot_2[i]/diag_2[i-1];
                    col_2[i]=col_2[i]-col_2[i-1]*bot_2[i]/diag_2[i-1];
                }

                psi_1[n-1]=col_1[n-1]/diag_1[n-1];
                psi_2[n-1]=col_2[n-1]/diag_2[n-1];

                for (i=n-2; i>=1; i--)
                {
                    col_1[i]=col_1[i]-col_1[i+1]*top_1[i]/diag_1[i+1];
                    psi_1[i]=col_1[i]/diag_1[i];

                    col_2[i]=col_2[i]-col_2[i+1]*top_2[i]/diag_2[i+1];
                    psi_2[i]=col_2[i]/diag_2[i];
                }

                psi_1[0]=0.0;

                for (i=1; i<n; i++)
                {
                    norm_p_1+=(psi_p_1[i]*psi_p_1[i]*i*delta*delta);
                    norm_1+=(psi_1[i]*psi_1[i]*i*delta*delta);

                    norm_p_2+=(psi_p_2[i]*psi_p_2[i]*i*delta*delta);
                    norm_2+=(psi_2[i]*psi_2[i]*i*delta*delta);
                }

                S1=pow((norm_p_1/norm_1),0.75);
                S2=pow((norm_p_2/norm_2),0.75);

                for (i=0; i<n; i++)
                {
                    psi_1[i]*=S1;
                    psi_2[i]*=S2;
                }

                for (i=0; i<n; i++)
                {
                    diag_3[i]=mid_3[i];
                    col_3[i]=psi_1[i]*psi_1[i]+sigma*psi_2[i]*psi_2[i];

                    theta_p1[i]=theta1[i];

                    diag_4[i]=mid_4[i];
                    col_4[i]=sigma*psi_1[i]*psi_1[i]+psi_2[i]*psi_2[i];
                }

                for (i=1; i<n; i++)
                {
                    diag_3[i]=diag_3[i]-top_3[i-1]*bot_3[i]/diag_3[i-1];
                    col_3[i]=col_3[i]-col_3[i-1]*bot_3[i]/diag_3[i-1];

                    diag_4[i]=diag_4[i]-top_4[i-1]*bot_4[i]/diag_4[i-1];
                    col_4[i]=col_4[i]-col_4[i-1]*bot_4[i]/diag_4[i-1];
                }

                theta1[n-1]=col_3[n-1]/diag_3[n-1];
                theta2[n-1]=col_4[n-1]/diag_4[n-1];

                for (i=n-2; i>=0; i--)
                {
                    col_3[i]=col_3[i]-col_3[i+1]*top_3[i]/diag_3[i+1];
                    theta1[i]=col_3[i]/diag_3[i];

                    col_4[i]=col_4[i]-col_4[i+1]*top_4[i]/diag_4[i+1];
                    theta2[i]=col_4[i]/diag_4[i];
                }

                for (i=1; i<n; i++)
                {
                    norm_p_3+=(theta_p1[i]*theta_p1[i]*i*delta*delta);
                    norm_3+=(theta1[i]*theta1[i]*i*delta*delta);
                }
            }

            lambda_2-=0.00001;
            lambda_1+=0.00001;
        }
        while(((prev-norm_2)*(prev-norm_2))>0.000001);
        printf("%lf %lf\n",sigma,((lambda_1-1.0)*2));
        fprintf(fp_diag,"%lf %lf\n",sigma,((lambda_1-1.0)*2));
    }

    fclose(fp_diag);
    fclose(fp_1);
    fclose(fp_2);
    fclose(fp_t1);
    fclose(fp_t2);
    fclose(fp_dev);
    fclose(fp_dev2);
    return 0;
}

