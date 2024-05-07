#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M 41                                                 
#define Dim 7
#define Min_num 2
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define DC 0.0000017                      
FILE *fp;
double diff[Dim][Dim],invd[Dim][Dim],invdx[Dim][Dim][Dim];
double e_eff,lambda=100.0;
double vx[Min_num][Dim]={0.18085, 0.19504, 0.26094, 0.07900, 0.16748, 0.32534, 0.06021 ,               	1, 0, 1, 0, 0, 0, 0       	};  
	
double saddle[Dim]={0,0,0,0,0,0}, diag[Dim][Dim]={1.0,0,0,0,0,0,0,1.0,0,0,0,0,0,0,1.0,0,0,0,0,0,0,1.0,0,0,0,0,0,0,1.0,0,0,0,0,0,0,1.0};
double x_max[Dim]={1,0.25,1,0.09,0.27,0.35,0.16},x_min[Dim]={0.0,0.0,0.0,0.0,0.0,0.0,0.0};
double a1=0.00158;

double veff(double x[]);
void path_mc(double xi[],double xf[],double path[M][Dim]);
double veff_min(double xi[]);
void diffx(double x[]);
void force(double xp[],double fvec[]);
void fdjac(double x[], double fvec[], double df[Dim][Dim],void (*vecfunc)(double [], double []));
double dl(double inv[Dim][Dim],double xi[],double xf[]);
double dot(double x[],double y[]);
double dshj(double xi[],double xf[]);
double dt(double xi[], double xf[]);
double xmax(double a,double b);
double xmin(double a,double b);
double time(double x[M][Dim]);
double shj_lambda(double x[M][Dim]);
double shj(double x[M][Dim]);

void path_mc(double xi[],double xf[],double x_path[M][Dim]);
double veff_min(double x[]);

main()
{
    int n,ndiff,fixpoint,i,j,k,nk=1,num,grid[Dim];
    double  itgl[Min_num],transition_time,temp=0,pb=0,phi,l=0,t=10;
    double c=0,e_old,amplitude=1,ds=0,deltat=0,smin=0,e1=0,e2=0,x1[Dim],x2[Dim];
    double ps[Min_num][M][Dim],psr[Min_num][M][Dim],x_path[M][Dim]={0};

    char st[20];
    num=21;
    
    
    e_eff = -xmin(veff(vx[0]),veff(vx[1])); 
    //e_eff=0.0;
   
	printf("e1=%f,e2=%f,e_eff=%f\n",veff(vx[0]),veff(vx[1]),e_eff);
	
 
    
	    fixpoint=1;
		//define the staring and ending position
		
		for(i=0;i<Dim;i++)
		{
			ps[fixpoint-1][0][i]=vx[1][i];           
			ps[fixpoint-1][M-1][i]=vx[0][i];         
		}
                path_mc(ps[fixpoint-1][0],ps[fixpoint-1][M-1],x_path);

		for(j=1;j<M-1;j++)
			for(i=0;i<Dim;i++)
				ps[fixpoint-1][j][i]=x_path[j][i];
				//printf("hehe\n");
				itgl[fixpoint-1]=shj(ps[fixpoint-1]);
  
	smin=itgl[0];
            printf("phi=%f,fixpoint=%d\n",itgl[fixpoint-1],fixpoint);
            sprintf(st,"pathbaDC%0.2f_a1_%0.1f.txt",DC,a1);           
	    printf("num=%d\n",num);
            if((fp=fopen(st,"w"))==NULL)
            {
                printf("Cannot open file. \n");
                exit(0);
            }
            //向文件中写入数据 
            for(i=0;i<M;i++)
            {

			fprintf(fp,"%f	%f	%f	%f	%f	%f	%f	%f\n",1.0*i,ps[fixpoint-1][i][0],ps[fixpoint-1][i][1],ps[fixpoint-1][i][2],ps[fixpoint-1][i][3],ps[fixpoint-1][i][4],ps[fixpoint-1][i][5],ps[fixpoint-1][i][6]);  //向文件中写入数据 

	        }

            fclose(fp);
            if((fp=fopen("shj.txt","a+"))==NULL)  
            {
                printf("Cannot open file. \n");
                exit(0);
            }
		transition_time=time(ps[fixpoint-1]);
		fprintf(fp,"num	M	a1	lambda	e_eff	x0	x1	x2	x3	x4	x5	x6	Shj	time\n");
		fprintf(fp,"%d	%d	%0.1f	%0.1e	%0.2f	%0.2f	%0.2f	%0.2f	%0.2f	%0.2f	%0.2f	%0.2f	%f	%f\n",num,M,a1,lambda,e_eff,ps[fixpoint-1][M-1][0],ps[fixpoint-1][M-1][1],ps[fixpoint-1][M-1][2],ps[fixpoint-1][M-1][3],ps[fixpoint-1][M-1][4],ps[fixpoint-1][M-1][5],ps[fixpoint-1][M-1][6],itgl[fixpoint-1],transition_time);
        fclose(fp);
}
double veff(double x[])
{
    int i,j,k;
    void force(double x[],double fvec[]);
	void fdjac(double x[], double fvec[], double df[Dim][Dim],void (*vecfunc)(double [], double []));
	double m=0.01,bb=0,cc=0,vf;
	double fvec[Dim],df[Dim][Dim]={0};
	diffx(x);
	force(x,fvec);
	fdjac(x,fvec,df,force);
    for(i=0;i<Dim;i++)
		for(j=0;j<Dim;j++)
		{
			bb=bb+invd[i][j]*fvec[i]*fvec[j];
			for(k=0;k<Dim;k++)
				cc=cc+diff[j][k]*invd[i][j]*df[i][k]+diff[j][k]*fvec[i]*invdx[k][i][j];
		}
   
     vf=0.25*bb+0.5*cc;  
	  
	return vf;
} 

void path_mc(double xi[],double xf[],double path[M][Dim])
{
	int ii,i,j,k,n,m,pj,mk,acceptance_num=0,cycle=1,times=5000;
	int inttau,tau,pos_start,pos_end,pos_ex,n_vari=3,Steps=300,si;//+++++
	long ra=1,*point;
	double itgl=0,itglr=0,l=0,r,t=1.0e-2,t_init,tfactor1=0.9,rho=0.05,tfactor2=0.7;
	double S,S1,delta,delta_x,avg_old,avg_new,fluc_old,fluc_new,dE;
	double psr[M][Dim],amplitude={0},ran1(long *idum),xc[Dim];
	point=&ra;
	for(i=0;i<Dim;i++){
		psr[0][i]=path[0][i]=xi[i];
		psr[M-1][i]=path[M-1][i]=xf[i];
	}
	for(i=1;i<=M-2;i++){
		for(n=0;n<Dim;n++){
			path[i][n]=psr[i][n]=psr[0][n]+(psr[M-1][n]-psr[0][n])*i/(M-1.0);
		}
	}
    itgl=shj(path);
    printf("shj_straight=%f\n",itgl);
        printf("i	amplitude	t	acceptance_ratio	shj	penalty\n");
		t_init=t;
        amplitude=0.01;//可以修改 
	for(i=1;i<=Steps;i++)
	{
		acceptance_num=0;
		for(si=1;si<=times;si++)
		{
			S = 0;
			S1 = 0;
			do
			{
				pos_start=rand()%M;
				pos_end=rand()%M;
				if (pos_start > pos_end) {
				pos_ex = pos_end;
				pos_end = pos_start;
				pos_start = pos_ex;
				}
			} while(pos_end-pos_start<2);

				for (tau=pos_start+1; tau<=pos_end; tau++)
				{
					S = S + dshj(path[tau-1],path[tau]);
					if (tau<pos_end)
					{
						for(n=0;n<Dim;n++)
						{
							inttau = rand();
							delta = amplitude*rand()/(double)RAND_MAX;
							delta_x = delta*(x_max[n]-x_min[n])/(double)M;
							//delta_x = delta/(double)M;
							if (inttau%n_vari == 0) psr[tau][n] = path[tau][n];
							if (inttau%n_vari == 1 && path[tau][n]+delta_x <= x_max[n]) psr[tau][n] = path[tau][n] + delta_x;
							if (inttau%n_vari == 2 && path[tau][n]-delta_x >= x_min[n]) psr[tau][n] = path[tau][n] - delta_x;
						}
					}
					for(m=0;m<Dim;m++) xc[m]=(psr[tau-1][m]+psr[tau][m])/2.0;
					if(e_eff+veff(xc)>=0) {
						S1 = S1 + dshj(psr[tau-1],psr[tau]);
					}
					else {
						S1 = S1 + 10000.0;
					}
				}
				avg_old = 0;
				avg_new = 0;
				for (k=pos_start+1; k<=pos_end; k++)
				{
					avg_old = avg_old + dl(diag,path[k-1],path[k]);
					avg_new = avg_new + dl(diag,psr[k-1],psr[k]);
				}
				avg_old = avg_old/(pos_end-pos_start);
				avg_new = avg_new/(pos_end-pos_start);
				fluc_old=0;
				fluc_new=0;
				for (k=pos_start+1; k<=pos_end; k++)
				{
					fluc_old = fluc_old + pow(dl(diag,path[k-1],path[k])-avg_old,2);
					fluc_new = fluc_new + pow(dl(diag,psr[k-1],psr[k])-avg_new,2);
				}
				S = S + lambda*fluc_old;
				S1 = S1 + lambda*fluc_new;
				if (S>S1)
				{
					for(k=pos_start+1;k<pos_end;k++)
					{
						for (n=0;n<Dim;n++) path[k][n] = psr[k][n];
					}
					S = S1;
					acceptance_num++;
				}
				for(k=pos_start+1;k<pos_end;k++)
					{
						for (n=0;n<Dim;n++) psr[k][n] = path[k][n];
					}
		}
		itgl=shj(path);
		rho=acceptance_num*1.0/times;
		printf("%d	%f	%0.2e	%0.3f	%f	%f\n",i,amplitude,t,rho,itgl,shj_lambda(path)-itgl);
		t=t*tfactor1;
	}
}

double veff_min(double xi[])
{
	int ii,i,j,k,n,m,pj,mk,acceptance_num=0,times=1000;
	int inttau,n_vari=2,Steps=100,si;
	double v0,v,l=0,r,t=1.0e-5,t_init,tfactor1=0.9,rho=0.05,tfactor2=0.7;
	double S,S1,delta,delta_x,x[Dim],x_new[Dim];
	double psr[M][Dim],amplitude={0},ran1(long *idum),xc[Dim];
	
	for (n=0;n<Dim;n++) x[n] = x_new[n] = xi[n];
    v0=veff(xi);
	t_init=t;
    amplitude=0.0001;
	for(i=1;i<=Steps;i++)
	{
		acceptance_num=0;
		for(si=0;si<=times;si++)
		{
			S = veff(x);
			for(n=0;n<Dim;n++)
				{
					inttau = rand();
					delta = amplitude*rand()/(double)RAND_MAX;
					delta_x = delta*(x_max[n]-x_min[n]);
					if (inttau%n_vari == 0 && x[n]+delta_x <= x_max[n]) x_new[n] = x[n] + delta_x;
					if (inttau%n_vari == 1 && x[n]-delta_x >= x_min[n]) x_new[n] = x[n] - delta_x;
				}
			S1 = veff(x_new);
			if (S>S1)
			{
				for (n=0;n<Dim;n++) x[n] = x_new[n];
				S = S1;
				acceptance_num++;
			}
		}
		v=veff(x);
		rho=acceptance_num*1.0/times;
		t=t*tfactor1; 	
	}
	for(j=0;j<Dim;j++) printf("x[%d]=%f	",j,x[j]);
	for(j=0;j<Dim;j++) xi[j] = x[j];
	printf("\n");
	printf("v0=%f	vf=%f\n",v0,v);
	return v;
}

void diffx(double x[])
{
	int i,j,k;
	
	diff[0][0]=DC;
	diff[0][1]=0;
	diff[0][2]=0;
	diff[0][3]=0;
	diff[0][4]=0;
	diff[0][5]=0;
    diff[0][6]=0;
	
	diff[1][0]=0;
    diff[1][1]=DC;
    diff[1][2]=0;
    diff[1][3]=0;
    diff[1][4]=0;
    diff[1][5]=0;
    diff[1][6]=0;
    
    diff[2][0]=0;
    diff[2][1]=0;
    diff[2][2]=DC;
    diff[2][3]=0;
    diff[2][4]=0;
    diff[2][5]=0;
    diff[2][6]=0;
    
    diff[3][0]=0;
    diff[3][1]=0;
    diff[3][2]=0;
    diff[3][3]=DC;
    diff[3][4]=0;
    diff[3][5]=0;
    diff[3][6]=0;
    
    diff[4][0]=0;
    diff[4][1]=0;
    diff[4][2]=0;
    diff[4][3]=0;
    diff[4][4]=DC;
    diff[4][5]=0;
    diff[4][6]=0;
    
    diff[5][0]=0;
    diff[5][1]=0;
    diff[5][2]=0;
    diff[5][3]=0;
    diff[5][4]=0;
    diff[5][5]=DC;
    diff[5][6]=0;
    
    
    diff[6][0]=0;
    diff[6][1]=0;
    diff[6][2]=0;
    diff[6][3]=0;
    diff[6][4]=0;
    diff[6][5]=0;
    diff[6][6]=DC;
    
    // inverse of diffusion coefficient
	invd[0][0]=1.0/diff[0][0];
	invd[0][1]=0;
	invd[0][2]=0;
	invd[0][3]=0;
	invd[0][4]=0;
	invd[0][5]=0;
	invd[0][6]=0;
	
	
	invd[1][0]=0;
	invd[1][1]=1.0/diff[1][1];
	invd[1][2]=0;
	invd[1][3]=0;
	invd[1][4]=0;
	invd[1][5]=0;
	invd[1][6]=0;
	
	invd[2][0]=0;
	invd[2][1]=0;
	invd[2][2]=1.0/diff[2][2];
	invd[2][3]=0;
	invd[2][4]=0;
	invd[2][5]=0;
	invd[2][6]=0;
	
	
	invd[3][0]=0;
	invd[3][1]=0;
	invd[3][2]=0;
	invd[3][3]=1.0/diff[3][3];
	invd[3][4]=0;
	invd[3][5]=0;
	invd[3][6]=0;
	
	
	invd[4][0]=0;
	invd[4][1]=0;
	invd[4][2]=0;
	invd[4][3]=0;
	invd[4][4]=1.0/diff[4][4];
	invd[4][5]=0;
	invd[4][6]=0;
	  
	
	invd[5][0]=0;
	invd[5][1]=0;
	invd[5][2]=0;
	invd[5][3]=0;
	invd[5][4]=0;
	invd[5][5]=1.0/diff[5][5];
	invd[5][6]=0;
	
	invd[6][0]=0;
	invd[6][1]=0;
	invd[6][2]=0;
	invd[6][3]=0;
	invd[6][4]=0;
	invd[6][5]=0;
	invd[6][6]=1.0/diff[6][6];
    // d D ^-1 / dx 
    for(i=0;i++;i<=Dim){
    	for(j=0;j++;j<=Dim){
    		for(k=0;k++;k<=Dim){
    			invdx[i][j][k]=0;
    		}
    	}
    }

}

void force(double xp[],double fvec[])
{
   	  double N = 10000;
	  double K = 100000;
	  double  w = 0.04;
	  double   r = 0.002;
	  double  q = 0.0064;
	 // double p=0.0018; 
	fvec[0] = q*xp[1] - a1*(K/N)*xp[3];                                                       
	fvec[1] = r*(1-xp[0]-xp[1]) - q*xp[1] ;              
	fvec[2]=q*xp[5] + w*xp[0]*xp[3]/(xp[0]+xp[1]) - 2*a1*(K/N)*(xp[2]*xp[3])/xp[0];    
	fvec[3]=2*a1*(K/N)*(xp[2]*xp[3])/xp[0] + q*xp[6] - r*xp[3] - w*xp[3] - a1*(xp[3]+(K/N)*xp[3]*xp[3]/xp[0]); 
	fvec[4]=a1*(xp[3]+(K/N)*xp[3]*xp[3]/xp[0]) - 2*r*xp[4];
	fvec[5]=r*xp[3] + w*xp[1]*xp[3]/(xp[0]+xp[1]) + 2*q*(1-xp[2]-xp[3]-xp[4]-xp[5]-xp[6]) - q*xp[5] - a1*(K/N)*xp[3]*xp[5]/xp[0] + w*xp[0]*xp[6]/(xp[0]+xp[1]);
	fvec[6]=2*r*xp[4] + a1*(K/N)*xp[3]*xp[5]/xp[0] - q*xp[6] - r*xp[6] - w*xp[6] ;	 
}

void fdjac(double x[], double fvec[], double df[Dim][Dim],void (*vecfunc)(double [], double []))
{
	int i,j;
	double h,temp,ff[Dim]={1.0},eps=1.0e-6;

	for (j=0;j<Dim;j++) {
		temp=x[j];
		h=eps*fabs(temp);
		if (h == 0.0) h=eps;
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(x,ff);
		x[j]=temp;
		for (i=0;i<Dim;i++){
			df[i][j]=(ff[i]-fvec[i])/h;
		}
	}
}

double dl(double inv[Dim][Dim],double xi[],double xf[])
{
    int i,j;
    double l=0;
    for(i=0;i<Dim;i++)
		for(j=0;j<Dim;j++)
			l=l+inv[i][j]*(xf[i]-xi[i])*(xf[j]-xi[j]);
    l=sqrt(l);
    return l;
}

double dot(double x[],double y[])
{
    int i=0;
    double sum=0;
    for(i=0;i<Dim;i++)
        sum=sum+x[i]*y[i];
    return sum;
}

double dshj(double xi[],double xf[])
{
	int i,j;
    double xc[Dim],ds=0,vi,fdl=0,dx[Dim],fd[Dim];
    double check_rest(double x[]),veff(double xi[]),dot(double x[],double y[]);
    double fvec[Dim];
	for(i=0;i<Dim;i++) xc[i]=(xf[i]+xi[i])/2;
	force(xc,fvec);
	vi=veff(xc);
	for(i=0;i<Dim;i++) dx[i]=xf[i]-xi[i];
	for(i=0;i<Dim;i++)
		for(j=0;j<Dim;j++)
			fdl=fdl+invd[i][j]*fvec[j]*dx[i];//——————————————————————散度力 
//	ds=dl(invd,xi,xf)*sqrt(vi+e_eff)-fdl/2.0;
		//ds=dl(invd,xi,xf)*sqrt(vi+e_eff);//-fdl/2.0;
		ds=dl(invd,xi,xf)*sqrt(vi+e_eff);//-fdl/2.0;  //可以修改的地方 
	
    return ds;
}

double dt(double xi[], double xf[])
{
	int i;
	double xc[Dim],vi,dt=0;
	double veff(double xi[]);
	for(i=0;i<Dim;i++) xc[i]=(xf[i]+xi[i])/2;
	//for(i=0;i<Dim;i++) xc[i]=xi[i];
	vi=veff(xc);
	dt=dl(diag,xi,xf)/sqrt(4*(vi+e_eff));
	return dt;
}

double xmax(double a,double b)
{
	return (a>b)?a:b;
}

double xmin(double a,double b)
{
	return (a<b)?a:b;
}

double time(double x[M][Dim])
{
	int i;
	double sum_t=0;
	double dt(double xi[],double xf[]);
	for(i=0;i<=M-2;i++)
		sum_t=sum_t+dt(x[i],x[i+1]);
	return sum_t;
}

double shj_lambda(double x[M][Dim])
{
	int i;
	double s=0,lave=0;
	double dl(double inv[Dim][Dim],double xi[],double xf[]),dshj(double xi[],double xf[]);
	lave=0;	
	for(i=0;i<=M-2;i++)
	{
		diffx(x[i]);
		lave=lave+dl(diag,x[i],x[i+1]);
	}
	lave=lave/(M-1);
	for(i=0;i<=M-2;i++)
	{
		s=s+dshj(x[i],x[i+1])+lambda*pow(dl(diag,x[i],x[i+1])-lave,2);
	}
	return s;
}

double shj(double x[M][Dim])
{
	int i;
	double s=0,lave=0;
	double dl(double inv[Dim][Dim],double xi[],double xf[]),dshj(double xi[],double xf[]);

	for(i=0;i<=M-2;i++)
	{
		s=s+dshj(x[i],x[i+1]);
	}
	return s;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
#undef M
#undef Dim
#undef Gn
#undef Min_num

