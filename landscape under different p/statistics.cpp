
#include <stdio.h> 
#include <math.h>
#include <stdlib.h> 

#define TINY 1.0e-20
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define Dim 7               
#define La 200
#define Lb 200
#define XMAX 1.0
#define D 0.0000017

void rk2(double h,double x[],double fvec[],double xmin[],double xmax[],long *point);
void force(double xp[],double fvec[]);
double ran1(long *idum);
double gasdev(long *idum); 

FILE *fp2;
double a1;

main()
{  
    long r=1,*point;point=&r;
	char st2[20];
	int index1,index2,A,B,filenum;
	double p[La][Lb]={0};
	double xc[2],xi[2],xf[2]; 
	double f_cyc[La][Lb]={0},f_cdk1[La][Lb]={0},fvec[Dim]={0};
	double iter=1.0,Tnf;
	double x_min[Dim]={0,0,0,0,0,0,0},x_max[Dim]={XMAX,0.25,XMAX,0.09,0.27,0.35,0.16}; //x_max[Dim]={1,1,1,1,1,1,1}; //  
 
	char a_1[20]="p"; 
	double Tn=1e12;
	double Tni=1e5;
	double h=0.005,x[Dim]={0.5,0.1,0.8,0.04,0.04,0.1,0.01};//1.586221,0,1.849857,0,1.849857,0
	index1 = 2;
	index2 = 3;

	Tnf = Tn+Tni;

	{	

	 a1 =0.00170;  //parameter 
	  printf("filenum:%d\n",filenum);
       for(A=0;A<La;A++)                           
        {
	        for(B=0;B<Lb;B++)
	        {
		        p[A][B]= 0.0;
				f_cyc[A][B] = 0.0;
				f_cdk1[A][B] = 0.0;                  
			}
		}

	for(iter=1.0;iter<=Tnf;iter=iter+1)                                      
	{
		rk2(h,x,fvec,x_min,x_max,point);   
			if(fmod(iter,10000000)==0)  //function of fmod is remainder 
	           printf("%e	%f\n",iter,iter*h);   
		xc[0]=x[index1-1]; xi[0]=x_min[index1-1]; xf[0]=x_max[index1-1];
		xc[1]=x[index2-1]; xi[1]=x_min[index2-1]; xf[1]=x_max[index2-1]; 

		if(iter>=Tni)               
		{
			A = (int)((xc[0]-xi[0])*La/(xf[0]-xi[0])); 
			B = (int)((xc[1]-xi[1])*Lb/(xf[1]-xi[1]));
			p[A][B] = p[A][B] + 1;
			f_cyc[A][B] = f_cyc[A][B] + fvec[index1-1];    
			f_cdk1[A][B] = f_cdk1[A][B] + fvec[index2-1];
		}	
	}
		sprintf(st2,"pp%d%d_%s=%0.5f.txt",index1,index2,a_1,a1);         
		fp2=fopen(st2,"w+");

		for(A=0;A<La;A++)
		{
			for(B=0;B<La;B++)
			{
					if(p[A][B]==0) 
					fprintf(fp2,"%d	%d	0	0.0	0.0\n",A,B);
					else 
					fprintf(fp2,"%d	%d	%f	%e	%e\n",A,B,p[A][B],f_cyc[A][B]/p[A][B],f_cdk1[A][B]/p[A][B]); 
				
			}         
		}
	fclose(fp2);
 	}
}
void rk2(double h,double x[],double q[],double xmin[],double xmax[],long *point)
{
	int i;
	double gauss[Dim],x1[Dim],q1[Dim];
		force(x,q);    
		for(i=0;i<Dim;i++)
		{    
        	gauss[i]=gasdev(point);

		}

		for(i=0;i<Dim;i++)
		{	

			
			x[i]=x[i]+q[i]*h+sqrt(h)*sqrt(2*D)*gauss[i]; 
			
			if(x[i]<xmin[i]) x[i]=2*xmin[i]-x[i];
			else if(x[i]>xmax[i]) x[i]=2*xmax[i]-x[i]; 
			else x[i]=x[i];
		}	
}
void force(double xp[],double fvec[])        
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
double gasdev(long *idum) 
{
	 
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	if (*idum < 0) iset=0;
	if  (iset == 0) 
	{
		do {
			v1=2.0*ran1(idum)-1.0; 
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} 
		while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} 
	else 
	{
		iset=0;
		return gset;
	}
}
double ran1(long *idum)   
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;
	if (*idum <= 0 || !iy) 
	{
		if (-(*idum) < 1) 
		*idum=1;
		else 
		*idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) 
		{
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) 
			*idum += IM;
			if (j < NTAB)
			 iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX)
	 return RNMX;
	else 
	 return temp;
}
