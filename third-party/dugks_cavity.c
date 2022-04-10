#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define L 1.0
#define M0 30
#define N0 30
#define M (M0+2)
#define N (N0+2)
#define M1 (M+1)
#define N1 (N+1)
#define Qx 3
#define Qy 3
#define imin 1
#define imax (N-2)
#define jmin 1
#define jmax (M-2)
#define Re 1000.0
#define CFL 0.1
#define sgn(x) ((x)==0.0?0.0:((x)<0.0?-1.0:1.0))
#define fmin(x,y) ((x)<(y)?(x):(y))
#define u_wall 0.1
#define a 2.0

void gks_ini(void);
void SlopeY(int i);
void SlopeX(int j);
void InterpX(int j);
void InterpY(int i);
void boundary(void);
void Evol(void);
double feq(int kx, int ky, double Ux, double Uy, double RHO);
void  datadeal(void);
double dx,dy, dt,RT,w,wb,tau,nu;
double f[M][N][Qx][Qy]; // f~ at the cell center
double f_plus[M][N][Qx][Qy]; // f_bar^+ at the cell center
double rho[M][N],ux[M][N],uy[M][N]; //cell center
double xc[N],yc[M],ti[M1];
double xf_face[N1][Qx][Qy], yf_face[M1][Qx][Qy];      // cell interface
//double minmod(double a, double b, double c);
double ex[Qx]={0., 1., -1.};
double ey[Qy]={0., 1., -1.};
int re[Qx]={0,2,1};
double tpx[Qx]={2.0/3.0, 1.0/6.0, 1.0/6.0};
double tpy[Qy]={2.0/3.0, 1.0/6.0, 1.0/6.0};

void main()
{
  int m,readdata,mmax;
  double err, u_old,t;

  RT=1.0/3;
  for(m=0;m<Qx;m++) {ex[m]=sqrt(3*RT)*ex[m];ey[m]=sqrt(3*RT)*ey[m];}
  nu=L*u_wall/Re;

  dx=L/(imax-imin+1);
  for(m=imin;m<=imax+1;m++)
  {
    t=(m-imin)*dx;
    ti[m]=(tanh(a*(t-0.5))/tanh(0.5*a)+1.0)/2;
  }
  for(m=imin;m<=imax;m++) xc[m]=yc[m]=0.5*(ti[m]+ti[m+1]);
  xc[imin-1]=-xc[imin];xc[imax+1]=2-xc[imax];
  yc[jmin-1]=-yc[jmin];yc[jmax+1]=2-yc[jmax];

  dy=dx=xc[imin]-xc[imin-1];
  dt=CFL*dx/sqrt(6*RT);
  tau=nu/RT;
  w=1.5*dt/(2*tau+dt);
  wb=dt/(2*tau+0.5*dt); //for interface

  printf("w=%e \n",w);
  gks_ini();
  u_old=ux[M/2][N/2];
  m=0;
AA:
  printf("input mmax:\n");
  scanf("%d",&mmax);
  mmax+=m;
  printf("dt=%lf mmax=%d   \n", dt, mmax);
  err=1.0;
  while((m<mmax) && err>1.0e-6)
     {
       m++;
       Evol();
       if(m%1000==0)
       {
         err=fabs(ux[M/2][N/2]-u_old)/u_wall;
         u_old=ux[M/2][N/2];
         printf("err=%e  u=%e  m=%d\n", err, u_old, m);
       }
     }
  datadeal();
  printf("Continue? (yes=1 no=0)\n");
  scanf("%d",&readdata);
  if(readdata) goto AA;
}

void gks_ini()
{
  int i,j, kx, ky;
  for(j=jmin;j<=jmax;j++) for(i=imin;i<=imax;i++)
  {
    ux[j][i]=uy[j][i]=0.0; rho[j][i]=1.0;
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)  f[j][i][kx][ky]=feq(kx,ky,ux[j][i],uy[j][i],rho[j][i]);
  }
}

double feq(int kx, int ky, double Ux, double Uy, double RHO)
{
  double uv,eu,x;

  eu=(ex[kx]*Ux+ey[ky]*Uy)/RT;
  uv=(Ux*Ux+Uy*Uy)/RT;
  x=tpx[kx]*tpy[ky]*RHO*(1.0+eu+0.5*(eu*eu-uv));
//  x=tpx[kx]*tpy[ky]*RHO*(1.0+3*eu+4.5*eu*eu-1.5*uv);
//  x=x*(1.0+tau*(ex[kx]-Ux)*force/RT);
  return x;
}
/*
double minmod(double a, double b, double c)
{
  double sa, sb, sc;
  sa=sgn(a); sb=sgn(b);sc=sgn(c);
  if(sa==sb && sb==sc)
     return sa*fmin(fabs(a),fmin(fabs(b),fabs(c)));
  else return 0;
}
*/

void boundary()
{
  int i,j, kx, ky;
  double rho_b,epsL, epsR;

//left & right walls
  epsL=(xc[imin]-xc[imin-1])/(xc[imin+1]-xc[imin]);
  epsR=(xc[imax+1]-xc[imax])/(xc[imax]-xc[imax-1]);
  for(j=jmin;j<=jmax;j++) for(ky=0;ky<Qy;ky++) for(kx=0;kx<Qx;kx++)
  {
    f_plus[j][imin-1][kx][ky]=(1+epsL)*f_plus[j][imin][kx][ky]-epsL*f_plus[j][imin+1][kx][ky];
    f_plus[j][imax+1][kx][ky]=(1+epsR)*f_plus[j][imax][kx][ky]-epsR*f_plus[j][imax-1][kx][ky];
  }

  for(j=jmin;j<=jmax;j++) //left wall
  {
    for(kx=0;kx<Qx;kx++) //bounce back
    {
      if(ex[kx]>0) for(ky=0;ky<Qy;ky++)
      f_plus[j][imin-1][kx][ky]=f_plus[j][imin-1][re[kx]][re[ky]]+f_plus[j][imin][re[kx]][re[ky]]-f_plus[j][imin][kx][ky];
    }
  }
  for(j=jmin;j<=jmax;j++) //right wall
  {
    for(kx=0;kx<Qx;kx++) //bounce back
    {
      if(ex[kx]<0) for(ky=0;ky<Qy;ky++)
      f_plus[j][imax+1][kx][ky]=f_plus[j][imax+1][re[kx]][re[ky]]+f_plus[j][imax][re[kx]][re[ky]]-f_plus[j][imax][kx][ky];
    }
  }

// top & bottom
  epsL=(yc[jmin]-yc[jmin-1])/(yc[jmin+1]-yc[jmin]);
  epsR=(yc[jmax+1]-yc[jmax])/(yc[jmax]-yc[jmax-1]);
  for(i=imin-1;i<=imax+1;i++) for(ky=0;ky<Qy;ky++) for(kx=0;kx<Qx;kx++)
  {
    f_plus[jmin-1][i][kx][ky]=(1+epsL)*f_plus[jmin][i][kx][ky]-epsL*f_plus[jmin+1][i][kx][ky];
    f_plus[jmax+1][i][kx][ky]=(1+epsR)*f_plus[jmax][i][kx][ky]-epsR*f_plus[jmax-1][i][kx][ky];
  }

  for(i=imin-1;i<=imax+1;i++) //bottom wall
  {
    for(ky=0;ky<Qy;ky++) //bounce back
    {
      if(ey[ky]>0) for(kx=0;kx<Qx;kx++)
      f_plus[jmin-1][i][kx][ky]=f_plus[jmin-1][i][re[kx]][re[ky]]+f_plus[jmin][i][re[kx]][re[ky]]-f_plus[jmin][i][kx][ky];
    }
  }

  for(i=imin-1;i<=imax+1;i++)
  {
    rho_b=0.0;
    for(ky=0;ky<Qy;ky++) //bounce back
    {
      if(ey[ky]==0) for(kx=0;kx<Qx;kx++) rho_b+=0.5*(f_plus[jmax+1][i][kx][ky]+f_plus[jmax][i][kx][ky]);
      else if(ey[ky]>0) for(kx=0;kx<Qx;kx++) rho_b+=(f_plus[jmax+1][i][kx][ky]+f_plus[jmax][i][kx][ky]);
    }
    for(ky=0;ky<Qy;ky++)
    {
      if(ey[ky]<0) for(kx=0;kx<Qx;kx++)
      f_plus[jmax+1][i][kx][ky]=f_plus[jmax+1][i][re[kx]][re[ky]]+f_plus[jmax][i][re[kx]][re[ky]]-f_plus[jmax][i][kx][ky]
                                +4*rho_b*tpx[kx]*tpy[ky]*(ex[kx]*u_wall)/RT;
    }
  }
}

void InterpX(int j)   // f at cell interface: X-direction
{
  int i, kx, ky, iL,jL,jR;
  double x, y, fc, dfx, dfy, ux_face, uy_face, rho_face;
  double hR,hL,AL,AR,AC;

  jL=j-1;  jR=j+1;
  hL=yc[j]-yc[jL];
  hR=yc[jR]-yc[j];
  AC=(hR-hL)/(hL*hR); AL=hR/(hL*(hR+hL)); AR=hL/(hR*(hR+hL));
  for(i=imin;i<=imax+1;i++) // inner nodes
  {
   iL=i-1;
   dx=xc[i]-xc[iL];
	for(ky=0;ky<Qy;ky++) for(kx=0;kx<Qx;kx++)
    {
      fc=0.5*(f_plus[j][i][kx][ky]+f_plus[j][iL][kx][ky]);
      dfx=(f_plus[j][i][kx][ky]-f_plus[j][iL][kx][ky])/dx;
      dfy=0.5*(AC*(f_plus[j][iL][kx][ky]+f_plus[j][i][kx][ky])
               + AR*(f_plus[jR][iL][kx][ky]+f_plus[jR][i][kx][ky])
               - AL*(f_plus[jL][iL][kx][ky]+f_plus[jL][i][kx][ky]));
      x=0.5*ex[kx]*dt; y=0.5*ey[ky]*dt;//half time step
      xf_face[i][kx][ky]=fc-x*dfx-y*dfy;
    }
  }

//the original f at interface
  for(i=imin;i<=imax+1;i++)
  {
    ux_face=uy_face=rho_face=0.0;
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
    {
      rho_face+=xf_face[i][kx][ky];
      ux_face+=ex[kx]*xf_face[i][kx][ky];
      uy_face+=ey[ky]*xf_face[i][kx][ky];
    }
    ux_face/=rho_face;    uy_face/=rho_face;
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
      xf_face[i][kx][ky]=(1.0-0.5*wb)*xf_face[i][kx][ky]+0.5*wb*feq(kx,ky, ux_face, uy_face,rho_face);
  }
}

void InterpY(int i)   // f at cell interface
{
  int j, kx, ky, iL, iR, jL;
  double fc, x, y, ux_face, uy_face, rho_face, dfx, dfy;
  double hR,hL,AL,AR,AC;

  iL=i-1; iR=i+1;
  hL=xc[i]-xc[iL];hR=xc[iR]-xc[i];
  AC=(hR-hL)/(hL*hR); AL=hR/(hL*(hR+hL)); AR=hL/(hR*(hR+hL));

// y-direction: no-slip
  for(j=jmin;j<=jmax+1;j++)
  {
    jL=j-1;
    dy=yc[j]-yc[jL];
    for(ky=0;ky<Qy;ky++) for(kx=0;kx<Qx;kx++)
    {
      fc=0.5*(f_plus[j][i][kx][ky]+f_plus[jL][i][kx][ky]);
      dfy=(f_plus[j][i][kx][ky]-f_plus[jL][i][kx][ky])/dy;
      dfx=0.5*(AC*(f_plus[jL][i][kx][ky]+f_plus[j][i][kx][ky])
               + AR*(f_plus[jL][iR][kx][ky]+f_plus[j][iR][kx][ky])
               - AL*(f_plus[jL][iL][kx][ky]+f_plus[j][iL][kx][ky]));
      y=0.5*ey[ky]*dt; x=0.5*ex[kx]*dt;//half time step
      yf_face[j][kx][ky]=fc-x*dfx-y*dfy;
    }
  }

//origional DFs
  for(j=jmin;j<=jmax+1;j++)
  {
    ux_face=uy_face=rho_face=0.0;
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
    {
      rho_face+=yf_face[j][kx][ky];
      ux_face+=ex[kx]*yf_face[j][kx][ky];
      uy_face+=ey[ky]*yf_face[j][kx][ky];
    }
    ux_face/=rho_face;    uy_face/=rho_face;
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
      yf_face[j][kx][ky]=(1.0-0.5*wb)*yf_face[j][kx][ky]+0.5*wb*feq(kx,ky, ux_face, uy_face,rho_face);
  }
}

void Evol()
{
  int i,j, kx, ky;
  double FM;

  // f_plus in each cell
  for(j=jmin;j<=jmax;j++) for(i=imin;i<=imax;i++) for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
  {
    FM=feq(kx, ky, ux[j][i], uy[j][i], rho[j][i]);
    f_plus[j][i][kx][ky]=f[j][i][kx][ky]-w*(f[j][i][kx][ky]-FM);
  }

  boundary();

  //update f: X-direction
  for(j=jmin;j<=jmax;j++)
  {
	 InterpX(j);
    for(i=imin;i<=imax;i++)
    {
      dx=(xc[i+1]-xc[i-1])*0.5;
      for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
      {
       f[j][i][kx][ky]=(4.0*f_plus[j][i][kx][ky]-f[j][i][kx][ky])/3.0
		  +ex[kx]*dt/dx*(xf_face[i][kx][ky]-xf_face[i+1][kx][ky]);
      }
    }
  }

  //update f: Y-direction
  for(i=imin;i<=imax;i++)
  {
	InterpY(i);
    for(j=jmin;j<=jmax;j++)
    {
      dy=(yc[j+1]-yc[j-1])*0.5;
      for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
      {
       f[j][i][kx][ky]=f[j][i][kx][ky]
		  +ey[ky]*dt/dy*(yf_face[j][kx][ky]-yf_face[j+1][kx][ky]);
      }
    }
  }

  //update macroscopic variables in each cell
  for(j=jmin;j<=jmax;j++) for(i=imin;i<=imax;i++)
  {
    rho[j][i]=ux[j][i]=uy[j][i]=0.0;
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
    {
      rho[j][i]+=f[j][i][kx][ky];
      ux[j][i]+=ex[kx]*f[j][i][kx][ky];
      uy[j][i]+=ey[ky]*f[j][i][kx][ky];
    }
    ux[j][i]/=rho[j][i];  uy[j][i]/=rho[j][i];
  }
}

void  datadeal()
   {
      int i,j;
      FILE *fp;

     fp=fopen("xc.dat","w");
     for (i=imin; i<=imax; i++) fprintf(fp,"%e ",xc[i]);
     fclose(fp);

     fp=fopen("yc.dat","w");
     for(j=jmin;j<=jmax;j++)   fprintf(fp,"%e ",yc[j]);
     fclose(fp);

     fp=fopen("rho.dat","w");
     for(j=jmin;j<=jmax;j++)	{
         for (i=imin; i<=imax; i++) fprintf(fp,"%e ",rho[j][i]);
         fprintf(fp,"\n");
      }
     fclose(fp);
     fp=fopen("u.dat","w");
     for(j=jmin;j<=jmax;j++)	{
         for (i=imin; i<=imax; i++) fprintf(fp,"%e ",ux[j][i]);
         fprintf(fp,"\n");
      }
     fclose(fp);
     fp=fopen("v.dat","w");
     for(j=jmin;j<=jmax;j++)	{
         for (i=imin; i<=imax; i++) fprintf(fp,"%e ",uy[j][i]);
         fprintf(fp,"\n");
      }
     fclose(fp);
}
