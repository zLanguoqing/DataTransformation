#include "stdio.h"
#include "DataTransformation.h"
int main()
{
    int i,j;
    double pr[64],pi[64],fr[64],fi[64];
    for(i = 0; i <= 63; i++)
    {
        pr[i] = exp(-0.1 *(i+0.5));
        pi[i] = 0.0;
    }
    printf("\n");
    for(i = 0; i <= 15; i++)
    {
        for(j = 0;j <=3;j++)
        {
            printf("%13.5e  ",pr[4*i+j]);

        }
        printf("\n");
    }
    printf("\n");
    FFT(pr,pi,64,6,fr,fi,0,1);
    for (i = 0; i <= 15; i++)
    {
        for (j = 0; j <= 3; j++)
        {
            printf("%13.5e  ", fr[4 * i + j]);
        }
        printf("\n");
    }
    printf("\n");
    for (i = 0; i <= 15; i++)
    {
        for (j = 0; j <= 3; j++)
        {
            printf("%13.5e  ", fi[4 * i + j]);
        }
        printf("\n");
    }
    printf("\n");
    for (i = 0; i <= 15; i++)
    {
        for (j = 0; j <= 3; j++)
        {
            printf("%13.5e  ", pr[4 * i + j]);
        }
        printf("\n");
    }
    printf("\n");
    for (i = 0; i <= 15; i++)
    {
        for (j = 0; j <= 3; j++)
        {
            printf("%13.5e  ", pi[4 * i + j]);
        }
        printf("\n");
    }
    printf("\n");
    FFT(fr,fi,64,6,pr,pi,1,1);
        for (i = 0; i <= 15; i++)
    {
        for (j = 0; j <= 3; j++)
        {
            printf("%13.5e  ", fr[4 * i + j]);
        }
        printf("\n");
    }
    printf("\n");
    int ii;
    double p1[8],x[8];
    for(ii = 0; ii <=7; ii++)
        p1[ii] = ii+1;
    
    Walsh(p1,8,3,x);
    printf("\n");
    for(ii = 0; ii<=7; ii++)    
        printf("x(%2d) = %13.5e \n",ii,x[ii]);
    printf("\n");

    return 0;
}