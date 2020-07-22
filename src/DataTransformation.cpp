#include "DataTransformation.h"
/* 
 快速傅里叶变换
 参数说明
 pr 当l为0，存放n个采样输入的实部，返回离散傅里叶变换的模；
    当l=1，存放傅里叶变换的n个实部，返回傅里叶变换的模长
 pi 当l为0，存放n个采样输入的虚部，返回离散傅里叶变换的幅角；
    当l=1，存放傅里叶变换的n个虚部，返回傅里叶变换的幅角，单位度

 n 采样点数  k 满足n=2^k
 fr 当l=0 返回傅里叶变换的n个实部，当l=1 返回逆傅里叶变换的n个实部
 fi 当l=0 返回傅里叶变换的n个虚部，当l=1 返回逆傅里叶变换的n个虚部
 l 当l =0 表示要求本函数计算傅里叶变换，当 l=1 表示要求本函数计算逆傅里叶变换
 il 当il =0 表示不要求本函数计算傅里叶变换或者逆傅里叶变换的模长与幅角
    当il =1 表示要求本函数计算傅里叶变换或者逆傅里叶变换的模长与幅角
*/
void FFT(double pr[], double pi[], int n, int k, double fr[], double fi[], int l, int il)
{
    int it, m, is, i, j, nv, l0;
    double p, q, s, vr, vi, poddr, poddi;
    for (it = 0; it <= n - 1; it++)
    {
        m = it;
        is = 0;
        for (i = 0; i <= k - 1; i++)
        {
            j = m / 2;
            is = 2 * is + (m - 2 * j);
            m = j;
        }
        fr[it] = pr[is];
        fi[it] = pi[is];
    }
    pr[0] = 1.0;
    pi[0] = 0.0;
    p = 6.283185306 / (1.0 * n);
    pr[1] = cos(p);
    pi[1] = -sin(p);
    if (l != 0)
        pi[1] = -pi[1];
    for (i = 2; i <= n - 1; i++)
    {
        p = pr[i - 1] * pr[1];
        q = pi[i - 1] * pi[1];
        s = (pr[i - 1] + pi[i - 1]) * (pr[1] + pi[1]);
        pr[i] = p - q;
        pi[i] = s - p - q;
    }
    for (it = 0; it <= n - 2; it = it + 2)
    {
        vr = fr[it];
        vi = fi[it];
        fr[it] = vr + fr[it + 1];
        fi[it] = vi + fi[it + 1];
        fr[it + 1] = vr - fr[it + 1];
        fi[it + 1] = vi - fi[it + 1];
    }
    m = n / 2;
    nv = 2;
    for (l0 = k - 2; l0 >= 0; l0--)
    {
        m = m / 2;
        nv = 2 * nv;
        for (it = 0; it <= (m - 1) * nv; it = it + nv)
            for (j = 0; j <= (nv / 2) - 1; j++)
            {
                p = pr[m * j] * fr[it + j + nv / 2];
                q = pi[m * j] * fi[it + j + nv / 2];
                s = pr[m * j] + pi[m * j];
                s = s * (fr[it + j + nv / 2] + fi[it + j + nv / 2]);
                poddr = p - q;
                poddi = s - p - q;
                fr[it + j + nv / 2] = fr[it + j] - poddr;
                fi[it + j + nv / 2] = fi[it + j] - poddi;
                fr[it + j] = fr[it + j] + poddr;
                fi[it + j] = fi[it + j] + poddi;
            }
    }
    if (l != 0)
        for (i = 0; i <= n - 1; i++)
        {
            fr[i] = fr[i] / (1.0 * n);
            fi[i] = fi[i] / (1.0 * n);
        }
    if (il != 0)
        for (i = 0; i <= n - 1; i++)
        {
            pr[i] = sqrt(fr[i] * fr[i] + fi[i] * fi[i]);
            if (fabs(fr[i]) < 0.000001 * fabs(fi[i]))
            {
                if ((fi[i] * fr[i]) > 0)
                    pi[i] = 90.0;
                else
                    pi[i] = -90.0;
            }
            else
            {
                pi[i] = atan(fi[i] / fr[i]) * 360.0 / 6.283185306;
            }
        }
    return;
}

/*
快速Walsh变换
p 存放长度为n=2^k的给定输入序列
n 输入序列的长度
k 满足n=2^k
x 返回输入序列pi的沃什变换序列
*/
void Walsh(double p[],int n,int k,double x[])
{
    int m,l,it,ii,i,j,is;
    double q;
    m = 1;
    l = n;
    it =2;
    x[0] = 1;
    ii = n/2;
    x[ii] = 2;
    for(i = 1; i <= k-1;i++)
    {
        m = m+ m;
        l = l /2;
        it = it +it;
        for(j = 0; j <= m-1; j++)
        {
            x[j*l+l/2] = it+1-x[j*l];

        }
    }
    for(i = 0; i <= n-1;i++)
    {
        ii = (int)(x[i] -1);
        x[i] = p[ii];
    }
    l=1;
    for(i =1; i <= k; i++)
    {
        m = n / (2*l)-1;
        for(j = 0; j <=m; j++)
        {
            it = 2*l*j;
            for(is = 0; is <= l-1; is++)
            {
                q = x[it +is] +x[it +is+l];
                x[it+is+l] = x[it+is]-x[it+is+l];
                x[it+is] =q;
            }
        }
        l = 2*l;
    }
    return;
}