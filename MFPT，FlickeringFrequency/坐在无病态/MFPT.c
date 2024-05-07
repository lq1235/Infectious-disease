
#include <math.h>
#include <stdlib.h>
#include "stdio.h"
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define TINY 1.0e-20;

#define Dim 7
#define La 200
#define Lb 200
#define XMAX 1.0
#define DC 0.0000017

double a1;
FILE *fp1, *fp2;

void mnewt(int ntrial, double x[], double tolx, double tolf);
double gasdev(long * idum);
void force(double x[], double fvec[]);
double ran1(long *idum);
double dl(double xi[], double xf[]);

main() 
{	
    long   r = 1, *point;
    point = &r;
    char st1[20], st2[20];
	int ntrial = 51;                 
    int i, j,  filenum, num;               
	double h_sqrt, gx[Dim], aux[Dim], xh[Dim], fvec[Dim], fxh[Dim],noise[Dim];   
	double x[Dim],  xi[Dim],xf[Dim];
	double ti, tau, tau2, var_t; 	
	double tolx = 1.0e-4, tolf = 1.0e-4;  
	//可以修改的参数 
	double   h =0.002; 
	double  r0 = 0.1; 
	double  Tn =1000500000;
	double x1[Dim] = {0}, x2[Dim] = {0}; 
	double diff[Dim] = {DC, DC, DC, DC, DC, DC};
	double x_max[Dim] = {XMAX, XMAX, XMAX, XMAX, XMAX, XMAX,XMAX}, x_min[Dim] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0};
	
	sprintf(st1, "MFPT_h=%0.3f_r=%0.2f.txt", h, r0); 
	fp1 = fopen(st1, "w+");

	for (filenum = 1; filenum <= 16; filenum++) 
	{   
	    tau = 0;   
		printf("filenum:%d\n", filenum);  
	
	 a1=0.0021- (filenum-1.0)*0.00005;
	
	       	x2[0] = 1;x2[1] = 0;x2[2] = 1;x2[3] = 0, x2[4] = 0;x2[5] = 0,x2[6] = 0; //无病态 
	       	x1[0] = 0.1483;x1[1] = 0.2027;x1[2] = 0.1891;x1[3] = 0.0721;x1[4] = 0.1902;x1[5] = 0.3268,x1[6] = 0.0748; //有病态 
			mnewt(ntrial,x1,tolx,tolf);    
			mnewt(ntrial,x2,tolx,tolf);    
		
		printf("x1=%f,%f,%f,%f,%f,%f,%f\n", x1[0], x1[1], x1[2], x1[3], x1[4], x1[5],x1[6]); 
		printf("x2=%f,%f,%f,%f,%f,%f,%f\n", x2[0], x2[1], x2[2], x2[3], x2[4], x2[5],x2[6]); 
		
		for (i = 0; i < Dim; i++)
		{
			xi[i] = x2[i]; 
		}	
		for (i = 0; i < Dim; i++)
		{
			xf[i] = x1[i]; 	
			gx[i] = sqrt(2 * diff[i]);
		}
		sprintf(st2, "FPT_h=%0.3f_r0=%0.2f_a1=%0.5f.txt", h, r0,a1); 
		fp2 = fopen(st2, "w+"); 
	
		for (num = 1; num <= 100; num++) 
		{
			for (i = 0; i < Dim; i++) 
			{
				x[i] = xi[i];        
			}
	
			h_sqrt = sqrt(h);
			for (j = 1; j <= Tn; j++) 
			{ 
				force(x, fvec);
				for (i = 0; i < Dim; i++) 
				{
					noise[i] = h_sqrt * gasdev(point);
					aux[i] = h * fvec[i] + noise[i] * gx[i];
					xh[i] = x[i] + aux[i];
				}
				force(xh, fxh);
				for (i = 0; i < Dim; i++) 
				{
					x[i] = x[i] + 0.5 * (aux[i] + h * fxh[i] + noise[i] * gx[i]);                         
				
					if (x[i] < x_min[i])                                                                 
						x[i] = 2 * x_min[i] - x[i]; 
					else if (x[i] > x_max[i])
						x[i] = 2 * x_max[i] - x[i];
					else
						x[i] = x[i];
				}
			
				if (dl(x, xf) <= r0) 
				{                             
					ti = j * h;               
					fprintf(fp2, "%d	%f\n", num, ti); 
                    printf("%d\n",num);   
					break;                     
				}
			}
			tau = tau + ti;    
			tau2 = tau2 + ti * ti; 
			
		}
       
		tau = tau / (num -1); 
		tau2 = tau2 / (num - 1); 
		var_t = (tau2 - tau * tau) / (tau * tau); 
		printf("a1=%0.5f,tau=%f,vart=%f\n", a1, tau, var_t);
		fprintf(fp1, "%0.5f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f\n", a1, tau, var_t, xi[0], xi[1], xi[2], xi[3], xi[4],xi[5], xi[6], xf[0],xf[1], xf[2], xf[3], xf[4], xf[5],xf[6]);
	}                                                                                                                               
	fclose(fp1);
	fclose(fp2);
}

void force(double xp[], double fvec[]) 
{

   	  double N = 10000;
	  double K = 100000;
	  double  w = 0.04;
	  double   r = 0.002;
	  double  q = 0.0064;

	fvec[0] = q*xp[1] - a1*(K/N)*xp[3];                                                       
	fvec[1] = r*(1-xp[0]-xp[1]) - q*xp[1] ;              
	fvec[2]=q*xp[5] + w*xp[0]*xp[3]/(xp[0]+xp[1]) - 2*a1*(K/N)*(xp[2]*xp[3])/xp[0];    
	fvec[3]=2*a1*(K/N)*(xp[2]*xp[3])/xp[0] + q*xp[6] - r*xp[3] - w*xp[3] - a1*(xp[3]+(K/N)*xp[3]*xp[3]/xp[0]); 
	fvec[4]=a1*(xp[3]+(K/N)*xp[3]*xp[3]/xp[0]) - 2*r*xp[4];
	fvec[5]=r*xp[3] + w*xp[1]*xp[3]/(xp[0]+xp[1]) + 2*q*(1-xp[2]-xp[3]-xp[4]-xp[5]-xp[6]) - q*xp[5] - a1*(K/N)*xp[3]*xp[5]/xp[0] + w*xp[0]*xp[6]/(xp[0]+xp[1]);
	fvec[6]=2*r*xp[4] + a1*(K/N)*xp[3]*xp[5]/xp[0] - q*xp[6] - r*xp[6] - w*xp[6] ;	  
	 
	
	
}

double fmax(double a, double b) 
{
	return (a > b) ? a : b;
}

double fmin(double a, double b) 
{
	return (a < b) ? a : b;
}

double dl(double xi[], double xf[]) 
{ 
	int i;
	double l = 0;
	for (i = 0; i < Dim; i++)
		l = l + (xf[i] - xi[i]) * (xf[i] - xi[i]);
	l = sqrt(l);

	return l;
}

double dot(double x[], double y[]) 
{
	int i = 0;
	double sum = 0;
	for (i = 0; i < Dim; i++)
		sum = sum + x[i] * y[i];
	return sum;
}

//Newton-Raphson methold solve the root of nonlinear equations.
void mnewt(int ntrial, double x[], double tolx, double tolf) 
{
	void lubksb(double a[Dim][Dim], int n, int indx[], double b[]);
	void ludcmp(double a[Dim][Dim], int n, int indx[], double * d);
	void fdjac(double x[], double fvec[], double df[Dim][Dim], void (*vecfunc)(double [], double []));
	int k, i, indxx[Dim] = {0};
	double errx, errf, d, a = 0;
	double pp[Dim], fvec[Dim], df[Dim][Dim] = {0};

	for (k = 0; k < ntrial; k++) 
	{
		//printf("x=%f,%f,%f,%f,%f,%f\n",x[0],x[1],x[2],x[3],x[4],x[5]);
		force(x, fvec);
		//printf("fvec1=%f,%f,%f,%f,%f,%f\n",fvec[0],fvec[1],fvec[2],fvec[3],fvec[4],fvec[5]);
		fdjac(x, fvec, df, force);
		//printf("fvec2=%f,%f,%f,%f,%f,%f\n",fvec[0],fvec[1],fvec[2],fvec[3],fvec[4],fvec[5]);
		errf = 0.0;
		for (i = 0; i < Dim; i++)
			errf += fabs(fvec[i]);
		if (errf <= tolf)
			return;
		for (i = 0; i < Dim; i++)
			pp[i] = -fvec[i];
		ludcmp(df, Dim, indxx, &d);
		lubksb(df, Dim, indxx, pp);
		errx = 0.0;
		for (i = 0; i < Dim; i++) 
		{
			errx += fabs(pp[i]);
			x[i] += pp[i];
		}
		//printf("x=%f,%f,%f,%f,%f,%f\n",x[0],x[1],x[2],x[3],x[4],x[5]);
		//getchar();
		if (errx <= tolx)  
			return;
	}
	return;
}

void ludcmp(double a[Dim][Dim], int n, int indx[], double *d) 
{
	int i, imax = 0, j, k;
	double big, dum, sum, temp;
	double vv[Dim];
	*d = 1.0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(a[i][j])) > big)
				big = temp;
		//if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		if (big == 0.0) {
			big = TINY
		}//{printf("Singular matrix in routine ludcmp");exit(0);}
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ( (dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0)
			a[j][j] = TINY;
		if (j != n - 1) {
			dum = 1.0 / (a[j][j]);
			for (i = j + 1; i < n; i++)
				a[i][j] *= dum;
		}
	}
}

void lubksb(double a[Dim][Dim], int n, int indx[], double b[]) 
{
	int i, ii = 0, ip, j;
	double sum;
	for (i = 0; i < n; i++) 
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii != 0)
			for (j = ii - 1; j < i; j++)
				sum -= a[i][j] * b[j];
		else if (sum != 0.0)
			ii = i + 1;
		b[i] = sum;
	}
	//printf("aa");
	for (i = n - 1; i >= 0; i--) 
	{
		sum = b[i];
		for (j = i + 1; j < n; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}

void fdjac(double x[], double fvec[], double df[Dim][Dim], void (*vecfunc)(double [], double [])) 
{
	int i, j;
	double h, temp, ff[Dim] = {1.0}, eps = 1.0e-6;
	for (j = 0; j < Dim; j++)
	{
		temp = x[j];
		h = eps * fabs(temp);
		if (h == 0.0)
			h = eps;
		x[j] = temp + h;
		h = x[j] - temp;
		(*vecfunc)(x, ff);
		x[j] = temp;
		for (i = 0; i < Dim; i++) 
		{
			df[i][j] = (ff[i] - fvec[i]) / h;   
		}
	}
}

double gasdev(long *idum) 
{

	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;

	if (*idum < 0)
		iset = 0;
	if  (iset == 0) 
	{
		do 
		{
			v1 = 2.0 * ran1(idum) - 1.0;
			v2 = 2.0 * ran1(idum) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} 
		while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	} else {
		iset = 0;
		return gset;
	}
}

double ran1(long *idum) 
{
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) 
	{
		if (-(*idum) < 1)
			*idum = 1;
		else *idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--) 
		{
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0)
				*idum += IM;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0)
		*idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX)
		return RNMX;
	else
		return temp;
}
