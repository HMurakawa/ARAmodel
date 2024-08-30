#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define Nx 400 
#define Nt 100000
#define Ntheta 70 // Ntheta%4!=0 Ntheta>=4*Nr
#define Nr 16

#define L1 2.80 //domain size
#define Radius 1.0 
double tend = 1000.0;

double dt = 1.0/(double)Nt;
double tstep = 1.0;//1.0/(double)Nt;
double dx = 2.0*(double)L1/((double)Nx);
double dtheta = 2.0*M_PI/Ntheta;
double dr = Radius/Nr;
double tol_stationary = 1.0e-7;


#define diffusioncoeficient 1.0



#define distance 0.62
#define init_radius_u .040
#define init_radius_v .120


#define a11 1.0
#define a12 1.0
#define a22 1.0
#define a21 a12

#define DUMPON 1

double F0_11 = 5000.0;
double Fr11 = -4000.0;
double Fa11 = 500.0;
double R0_11 = .14;
double R1_11 = 0.55;
double R2_11 = 0.88;


double F0_12 = 50.0;
double F0_13 = 0.0;
double F0_21 = 50.0;
double F0_22 = 25.0;
double F0_23 = .0;
double F0_31 = 0.0;
double F0_32 = .0;
double F0_33 = 100.0;

double Fr12 = 0.0;
double Fr13 = 0.0;
double Fr21 = 0.0;
double Fr22 = 0.0;
double Fr23 = 0.0;
double Fr31 = 0.0;
double Fr32 = 0.0;
double Fr33 = 0.0;

double Fa12 = (0.0);
double Fa13 = (0.0);
double Fa21 = 0.0;
double Fa22 = 0.0;
double Fa23 = 0.0;
double Fa31 = 0.0;
double Fa32 = 0.0;
double Fa33 = 0.0;

double R0_12 = .3;
double R0_13 = .14;
double R0_21 = .3;
double R0_22 = .14;
double R0_23 = .14;
double R0_31 = .14;
double R0_32 = .14;
double R0_33 = .14;

double R1_12 = 0.55;
double R1_13 = 0.55;
double R1_21 = .55;
double R1_22 = .55;
double R1_23 = .55;
double R1_31 = 0.55;
double R1_32 = 0.55;
double R1_33 = 0.55;

double R2_12 = 0.88;
double R2_13 = 0.88;
double R2_21 = .88;
double R2_22 = .88;
double R2_23 = .88;
double R2_31 = .88;
double R2_32 = .88;
double R2_33 = .88;



#define c10 .0
#define c11 .0
#define c12 .0
#define c20 .0
#define c21 .0
#define c22 .0


//5 20 5 20 0 0 0 0
//20 10 5 10 0 0 0 0
//20 19 5 19


double cutoff=1.0e-2;//1.0/Nt; //for visualization

void init(double u[],double v[],double w[],double z[]);
void init2(double u[],double v[],double w[],double z[]);
void init3(double u[],double v[],double w[],double z[]);
void calc_Ju();
//void calc_Jv();
void SOLVE(void);
double xiu(double u, double v, double w);
double xiv(double u, double v, double w);
double xiw(double u, double v, double w);

double C1, C2, C3, C4;
double calc_Hstab();
double calc_Hstab2();

int num_spots();
int count_numspots_sqare(int N);
double dist(double u[], double v[]);
double calc_ip(double u[], double v[], double w[]);
double dist_min(int num);
void angles(double xtmp[][2], int N, double distmin);
double ave_mindist(int num);
double med_mindist(int num);
int count_dist1(int N, double ave);
int count_distsqrt2(int N, double ave);
int count_bonds(int N, double ave);
double count_surrounding(int N, double ave);
void min_arearatio_v_around_spot(int N, double rad, double *min, double *ave);
double norm_l1_3(double u[], double v[], double w[]);

double u[Nx*Nx], v[Nx*Nx], w[Nx*Nx], z0[Nx*Nx], neuron[Nx*Nx], z[Nx*Nx], Jusx[Nx*Nx], Jusy[Nx*Nx], Jvsx[Nx*Nx], Jvsy[Nx*Nx], Juwx[Nx*Nx], Juwy[Nx*Nx], Jvwx[Nx*Nx], Jvwy[Nx*Nx], 
Jwsx[Nx*Nx], Jwsy[Nx*Nx], Jwwx[Nx*Nx], Jwwy[Nx*Nx];
double u_old[Nx*Nx], v_old[Nx*Nx], w_old[Nx*Nx], tmp1[Nx*Nx], tmp2[Nx*Nx], tmp3[Nx*Nx], tmp4[Nx*Nx];
double u_dif[Nx*Nx], v_dif[Nx*Nx], w_dif[Nx*Nx];
double spot[Nx*Nx][2];

FILE *fptr, *fpinit;
char fname[32];
int counterdata = 0;


void main(void)
{
  int     i,j,n;
  long    kaisu;
  double  x, y, t;
  
  double vol_v, vol_v_init;

  kaisu = tstep / dt;
  dt = tstep / kaisu;


  t = 0.0;
  
  init2(u,v,w,z);

  for(i=0; i<Nx*Nx; i++)
    z0[i]=0.0;
    

  double aerr, err, vol;
  while (t <= tend) {
	
	for(i=0; i<Nx*Nx; i++){
		u_old[i] = u[i]; 
		v_old[i] = v[i]; 
		w_old[i] = w[i];
	}		
		
    for (n = 0; n < kaisu; ++n){    
      
	  if(n%100==0)
	      calc_Ju();
      
      SOLVE();

      t += dt;
    }
	
	for(i=0; i<Nx*Nx; i++){
		u_dif[i] = u[i]-u_old[i]; 
		v_dif[i] = v[i]-v_old[i]; 
		w_dif[i] = w[i]-w_old[i]; 
	}		
  }

  for (j=0; j<Nx; j++){
      y=((j+0.5)*dx-L1);
    for (i=0; i<Nx; i++){
      x=((i)*dx-L1);
		printf("%f %f %f\n", x, y, u[i+Nx*j]);
	}
	printf("\n");
  }


}	    	




double omega11(double r){
	if(r<R0_11)
		return F0_11;
	if(r<R1_11)
		return Fr11;
	if(r<R2_11)
		return 0.0;//Fr*(R2-r)/(R2-R1);
	if(r<=1)
		return Fa11;//4.0*Fa*(1.0-r)*(r-R2)/(1.0-R2)/(1.0-R2);
}

double omega12(double r){
	if(r<R0_12)
		return F0_12;
	if(r<R1_12)
		return Fr12;
	if(r<R2_12)
		return 0.0;//Fr*(R2-r)/(R2-R1);
	if(r<=1)
		return Fa12;//4.0*Fa*(1.0-r)*(r-R2)/(1.0-R2)/(1.0-R2);
}

double omega13(double r){
	if(r<R0_13)
		return F0_13;
	if(r<R1_13)
		return Fr13;
	if(r<R2_13)
		return 0.0;//Fr*(R2-r)/(R2-R1);
	if(r<=1)
		return Fa13;//4.0*Fa*(1.0-r)*(r-R2)/(1.0-R2)/(1.0-R2);
}

double omega21(double r){
	if(r<R0_21)
		return F0_21;
	if(r<R1_21)
		return Fr21;
	if(r<R2_21)
		return 0.0;//Fr*(R2-r)/(R2-R1);
	if(r<=1)
		return Fa21;//4.0*Fa*(1.0-r)*(r-R2)/(1.0-R2)/(1.0-R2);
}

double omega22(double r){
	if(r<R0_22)
		return F0_22;
	if(r<R1_22)
		return Fr22;
	if(r<R2_22)
		return 0.0;//Fr*(R2-r)/(R2-R1);
	if(r<=1)
		return Fa22;//4.0*Fa*(1.0-r)*(r-R2)/(1.0-R2)/(1.0-R2);
}

double omega23(double r){
	if(r<R0_23)
		return F0_23;
	if(r<R1_23)
		return Fr23;
	if(r<R2_23)
		return 0.0;//Fr*(R2-r)/(R2-R1);
	if(r<=1)
		return Fa23;//4.0*Fa*(1.0-r)*(r-R2)/(1.0-R2)/(1.0-R2);
}

double omega31(double r){
	if(r<R0_31)
		return F0_31;
	if(r<R1_31)
		return Fr31;
	if(r<R2_31)
		return 0.0;//Fr*(R2-r)/(R2-R1);
	if(r<=1)
		return Fa31;//4.0*Fa*(1.0-r)*(r-R2)/(1.0-R2)/(1.0-R2);
}

double omega32(double r){
	if(r<R0_32)
		return F0_32;
	if(r<R1_32)
		return Fr32;
	if(r<R2_32)
		return 0.0;//Fr*(R2-r)/(R2-R1);
	if(r<=1)
		return Fa32;//4.0*Fa*(1.0-r)*(r-R2)/(1.0-R2)/(1.0-R2);
}

double omega33(double r){
	if(r<R0_33)
		return F0_33;
	if(r<R1_33)
		return Fr33;
	if(r<R2_33)
		return 0.0;//Fr*(R2-r)/(R2-R1);
	if(r<=1)
		return Fa33;//4.0*Fa*(1.0-r)*(r-R2)/(1.0-R2)/(1.0-R2);
}


int num_spots(){
	int i,j, num=0;
	
    for(j=1; j<Nx-1; j++){
		for(i=1; i<Nx-1; i++){
			if(u[(i)+Nx*(j)]>0.1 && u[(i)+Nx*(j)]>u[(i+1)+Nx*(j)] && u[(i)+Nx*(j)]>u[(i-1)+Nx*(j)] && u[(i)+Nx*(j)]>u[(i)+Nx*(j+1)] && u[(i)+Nx*(j)]>u[(i)+Nx*(j-1)]){
				spot[num][0] = (double)i*dx -L1;
				spot[num][1] = (double)j*dx -L1;
				num++;
			}
		}
	}
	return num;
}



double norm_l1_3(double u[], double v[], double w[]){
	int i;
	double l1 = 0.0;
	for(i=0; i<Nx*Nx; i++){
		l1 += fabs(u[i]);
		l1 += fabs(v[i]);
		l1 += fabs(w[i]);
	}
	l1 /= (double)(Nx*Nx);	
	return l1;
}


double dist(double u[], double v[]){
  double sum = 0;
  int i;
  for(i=0; i<2; i++)
    sum += (u[i]-v[i])*(u[i]-v[i]);
  return sqrt(sum);
}

double calc_ip(double u[], double v[], double w[]){
  double sum = 0;
  int i;
  for(i=0; i<2; i++)
    sum += (v[i]-u[i])*(w[i]-u[i]);
  return sum;
}

void bubble_sort(double a[], int N){
  int i, j;
  double tmp;
  
  for(i=0; i<N-1; i++){
    for(j=1; j<N-i; j++){
      if(a[j-1]<a[j]){
        tmp = a[j-1];
        a[j-1] = a[j];
        a[j] = tmp;
      }
    }
  }
}

double med_mindist(int num){
	int i,j;
	double dij, min[num], med=0;
	
	for(i=0; i<num; i++){
		min[i]=L1*2;
	    for(j=0; j<num; j++){
			if(i!=j){
				dij = dist(spot[i],spot[j]);
				if(min[i]>dij)
					min[i]=dij;
			}
		}
	}
	bubble_sort(min,num);
	return min[num/2];
}


double dist_min(int num){
	int i,j;
	double dij, min=L1*2;
	
	for(i=0; i<num; i++){
	    for(j=i+1; j<num; j++){
			dij = dist(spot[i],spot[j]);
			if(min>dij)
				min=dij;
		}
	}
	return min;
}

double ave_mindist(int num){
	int i,j;
	double dij, min, ave=0;
	
	for(i=0; i<num; i++){
		min=L1*2;
	    for(j=0; j<num; j++){
			if(i!=j){
				dij = dist(spot[i],spot[j]);
				if(min>dij)
					min=dij;
			}
		}
		ave += min;
	}
	return ave/num;
}

double count_surrounding(int N, double ave){
	int i,j, counter=0, num=0;
	double dij;
	
    for(i=0; i<N; i++){
	  if(spot[i][0]>-L1*0.6 && spot[i][0]<L1*0.6 && spot[i][1]>-L1*0.6 && spot[i][1]<L1*0.6){
//		  num=0;
          counter++;
		  for(j=0; j<N; j++){
			  dij = dist(spot[i],spot[j]);
			  if(j!=i && dij < (((ave)*(1.41421356)*1.1)))
				  num++;
		  }
        }
//		printf("%d\n", num);
      }
      return (double)num/(double)counter;
}


int count_numspots_sqare(int N){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<N; i++){
	  if(spot[i][0]>-L1*0.6 && spot[i][0]<L1*0.6 && spot[i][1]>-L1*0.6 && spot[i][1]<L1*0.6){
          counter++;
	  }
      }
      return counter;
}

int count_bonds(int N, double ave){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<N; i++){
	  if(spot[i][0]>-L1*0.6 && spot[i][0]<L1*0.6 && spot[i][1]>-L1*0.6 && spot[i][1]<L1*0.6){
      for(j=0; j<N; j++){
          dij = dist(spot[i],spot[j]);
          if(dij > ((ave)*0.8) && dij < (((ave)*(1.41421356)*1.1)))
          counter++;
        }
	  }
      }
      return counter/2;
}


int count_dist1(int N, double ave){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<N; i++){
	  if(spot[i][0]>-L1*0.6 && spot[i][0]<L1*0.6 && spot[i][1]>-L1*0.6 && spot[i][1]<L1*0.6){
      for(j=0; j<N; j++){
          dij = dist(spot[i],spot[j]);
          if(dij > ((ave)*0.9) && dij < ((ave)*1.1))
          counter++;
        }
	  }
      }
      return counter/2;
}

int count_distsqrt2(int N, double ave){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<N; i++){
	  if(spot[i][0]>-L1*0.6 && spot[i][0]<L1*0.6 && spot[i][1]>-L1*0.6 && spot[i][1]<L1*0.6){
      for(j=0; j<N; j++){
          dij = dist(spot[i],spot[j]);
          if(dij > (ave)*(1.41421356)*0.9 && dij < ((ave)*(1.41421356)*1.1))
          counter++;
        }
	  }
      }
      return counter/2;
}


void min_arearatio_v_around_spot(int N, double rad, double *min, double *ave){
	int i,j,n,m, count=0;
	double area[4], tmp;
	double p[2];
	double tol = 0.01;

	*min = 1.0;
	*ave = 0.0;

    for(n=0; n<N; n++){
		for(m=0; m<4; m++)
			area[m] = 0.0;
		
    for(j=1; j<Nx-1; j++){
		for(i=1; i<Nx-1; i++){
			p[0] = (double)i*dx -L1;
			p[1] = (double)j*dx -L1;
				if( dist(p,spot[n])<rad && v[(i)+Nx*(j)]>tol){
					if(p[0]<spot[n][0])  area[0]+=1.0;
					if(p[0]>spot[n][0])  area[1]+=1.0;
					if(p[1]<spot[n][1])  area[2]+=1.0;
					if(p[1]>spot[n][1])  area[3]+=1.0;
				}
			}
		}

		bubble_sort(area,4);
		if(area[0]!=0){
			tmp = area[3]/area[0];
			count++;
			*ave += tmp;
			if(*min>tmp) *min = tmp;
		}

	}
	*ave /= count;
}

void angles(double xtmp[][2], int N, double distmin){
  int i, j, k, m;
  double dij, dkj, dik, ip, V;

/*
  FILE *fp1;
  char fname1[32];

  sprintf(fname1,"angles.txt");
  fp1 = fopen(fname1, "w");
  */
    for(i=0; i<N; i++){
      for(j=i+1; j<N; j++){
          dij = dist(xtmp[i],xtmp[j]);
          if(dij > ((distmin)*0.8) && dij < ((distmin)*1.6))
//             fprintf(fp1, "%.15f\n", dij); 
             printf("dist: %.2f\n", dij); 
        }
      }

//  fclose(fp1);

////////

/*
  if(sigma_const<0.0)
    sprintf(fname1,"angles10.0.3_b%.1lf_sf_a%.1lf_t%.1f.txt", beta, alpha, t);
  else
    sprintf(fname1,"angles10.0.3_b%.1lf_s%.2lf_a%.1lf_t%.1f.txt", beta, sigma_const, alpha, t);
  fp1 = fopen(fname1, "w");
*/

    for(i=0; i<N; i++){
       for(j=0; j<N; j++){
        if(j!=i){
          dij = dist(xtmp[i],xtmp[j]);
          if(dij > ((distmin)*0.8) && dij < ((distmin)*1.6)){
            for(k=j+1; k<N; k++){
              if(k!=i){              
                dik = dist(xtmp[i],xtmp[k]);
                if(dij > ((distmin)*0.8) && dik < ((distmin)*1.6)){
                  ip = calc_ip(xtmp[i],xtmp[j],xtmp[k]);
//                  fprintf(fp1, "%.15f\n", acos(ip/dij/dik));
                  printf("angle: %.2f\n", acos(ip/dij/dik)*360.0/M_PI); 
//                  fprintf(fp1, "%d, %d, %d\n", i,j,k); 
                }
              }
            }
          
          }
        }
      }  
    }

//  fclose(fp1);



}

 
 

double g_outside(double u, double v, double w){
	return 1.0-u-v-w;
}




double xiu(double u, double v, double w){
	return u+v+w;
}

double xiv(double u, double v, double w){
	return u+v+w;
}

double xiw(double u, double v, double w){
	return u+v+w;
}


void calc_Ju(){
  int i, j, it, ir, ixx, iyy;
  double theta, r, x, y, xx, yy, S11, S12, S13, S21, S22, S23, S31, S32, S33, etax, etay;
  
  for (j=0; j<Nx; j++){
      y=((j)*dx-L1);
    for (i=0; i<Nx; i++){
      x=((i+0.5)*dx-L1);
      
      //
      Jusy[i+Nx*j] = 0.0;
       for(it=0; it<Ntheta; it++){
				theta = it*dtheta+M_PI/2.0;
				etax = cos(theta); 	etay = sin(theta);
				for(ir=1; ir<=Nr; ir++){
					r = (ir-0.5)*dr;
					S11 = omega11(r)*r*dtheta*dr;
					
					xx = x+etax*r;
					if(xx<-L1) xx += 2*L1;
					else if(xx>L1) xx -= 2*L1;
					if(etax>0)
  					ixx = (xx+L1)/dx;
  				else
  				  ixx = -(int)((-xx-L1)/dx);
					yy = y+etay*r;
					if(yy<-L1) yy += 2*L1;
					else if(yy>L1) yy -= 2*L1;
					if(etay>0)
  					iyy = (yy+L1)/dx;
  			  else
  			    iyy = -(int)((-yy-L1)/dx);
					
					Jusy[i+Nx*j] += ( S11*u[ixx+Nx*iyy]+S12*v[ixx+Nx*iyy]+S13*w[ixx+Nx*iyy] )*etay; 
				}
			}
			
			if(j==0){
  			Jusy[i+Nx*j] *= (g_outside(u[i+Nx*(Nx-1)],v[i+Nx*(Nx-1)],w[i+Nx*(Nx-1)])+g_outside(u[i+Nx*(j)],v[i+Nx*(j)],w[i+Nx*(j)]))/2.0;
	  	}
			else{
			  Jusy[i+Nx*j] *= (g_outside(u[i+Nx*(j-1)],v[i+Nx*(j-1)],w[i+Nx*(j-1)])+g_outside(u[i+Nx*(j)],v[i+Nx*(j)],w[i+Nx*(j)]))/2.0;
			}
			
    }
  }

  for (j=0; j<Nx; j++){
      y=((j+0.5)*dx-L1);
    for (i=0; i<Nx; i++){
      x=((i)*dx-L1);
      
      Juwx[i+Nx*j] = 0.0;

       for(it=0; it<Ntheta; it++){
				theta = it*dtheta;
				etax = cos(theta); 	etay = sin(theta);
				for(ir=1; ir<=Nr; ir++){
					r = (ir-0.5)*dr;
					S11 = omega11(r)*r*dtheta*dr;
					
					xx = x+etax*r;
					if(xx<-L1) xx += 2*L1;
					else if(xx>L1) xx -= 2*L1;
					if(etax>0)
  					ixx = (xx+L1)/dx;
  				else
  				  ixx = -(int)((-xx-L1)/dx);
					yy = y+etay*r;
					if(yy<-L1) yy += 2*L1;
					else if(yy>L1) yy -= 2*L1;
					if(etay>0)
  					iyy = (yy+L1)/dx;
  			  else
  			    iyy = -(int)((-yy-L1)/dx);
					
					Juwx[i+Nx*j] += ( S11*u[ixx+Nx*iyy]+S12*v[ixx+Nx*iyy]+S13*w[ixx+Nx*iyy] )*etax; 
				}
			}
			
			if(i==0){
  			Juwx[i+Nx*j] *= (g_outside(u[(Nx-1)+Nx*j],v[(Nx-1)+Nx*j],w[(Nx-1)+Nx*j])+g_outside(u[(i)+Nx*j],v[(i)+Nx*j],w[(i)+Nx*j]))/2.0; 
	  	}
	  	else{
				Juwx[i+Nx*j] *= (g_outside(u[(i-1)+Nx*j],v[(i-1)+Nx*j],w[(i-1)+Nx*j])+g_outside(u[(i)+Nx*j],v[(i)+Nx*j],w[(i)+Nx*j]))/2.0; 
			} 
    }
  }

}


void SOLVE(void){
  int i,j; 
  double g, lam=dt/dx;


  for(i=0; i<Nx*Nx; i++){
			tmp1[i] = u[i];
	}
	

  for(j=0; j<Nx; j++){
    for(i=0; i<Nx-1; i++){
      g = ( -(xiu(tmp1[i+1+Nx*j],tmp2[i+1+Nx*j],tmp3[i+1+Nx*j])-xiu(tmp1[i+Nx*j],tmp2[i+Nx*j],tmp3[i+Nx*j]))/dx + Juwx[i+1+Nx*j] )*lam;
      if(g>=0.0)
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j)];
      else
        u[i+Nx*j] -= g*tmp1[(i+1)+Nx*(j)];
    }
    i=Nx-1;
      g = ( -(xiu(tmp1[Nx*j],tmp2[Nx*j],tmp3[Nx*j])-xiu(tmp1[i+Nx*j],tmp2[i+Nx*j],tmp3[i+Nx*j]))/dx + Juwx[Nx*j] )*lam;
      if(g>=0.0)
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j)];
      else
        u[i+Nx*j] -= g*tmp1[Nx*(j)];
  }

  for(j=0; j<Nx-1; j++){
    for(i=0; i<Nx; i++){
      g = ( -(xiu(tmp1[i+Nx*(j+1)],tmp2[i+Nx*(j+1)],tmp3[i+Nx*(j+1)])-xiu(tmp1[i+Nx*(j)],tmp2[i+Nx*(j)],tmp3[i+Nx*(j)]))/dx + Jusy[i+Nx*(j+1)] )*lam;
      if(g>=0.0)
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j)];
      else
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j+1)];
    }
  }
  j=Nx-1; 
    for(i=0; i<Nx; i++){
      g = ( -(xiu(tmp1[i],tmp2[i],tmp3[i])-xiu(tmp1[i+Nx*(j)],tmp2[i+Nx*(j)],tmp3[i+Nx*(j)]))/dx + Jusy[i] )*lam;
      if(g>=0.0)
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j)];
      else
        u[i+Nx*j] -= g*tmp1[(i)];
    }


  for(j=0; j<Nx; j++){
		i=0;
      g = ( (xiu(tmp1[i+Nx*j],tmp2[i+Nx*j],tmp3[i+Nx*j])-xiu(tmp1[Nx-1+Nx*j],tmp2[Nx-1+Nx*j],tmp3[Nx-1+Nx*j]))/dx - Juwx[i+Nx*j] )*lam;
      if(g>=0.0)
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j)];
      else
        u[i+Nx*j] -= g*tmp1[(Nx-1)+Nx*(j)];		
    for(i=1; i<Nx; i++){
      g = ( (xiu(tmp1[i+Nx*j],tmp2[i+Nx*j],tmp3[i+Nx*j])-xiu(tmp1[i-1+Nx*j],tmp2[i-1+Nx*j],tmp3[i-1+Nx*j]))/dx - Juwx[i+Nx*j] )*lam;
      if(g>=0.0)
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j)];
      else
        u[i+Nx*j] -= g*tmp1[(i-1)+Nx*(j)];
    }
  }

  j=0;
    for(i=0; i<Nx; i++){
      g = ( (xiu(tmp1[i+Nx*(j)],tmp2[i+Nx*(j)],tmp3[i+Nx*(j)])-xiu(tmp1[i+Nx*(Nx-1)],tmp2[i+Nx*(Nx-1)],tmp3[i+Nx*(Nx-1)]))/dx - Jusy[i+Nx*j] )*lam;
      if(g>=0.0)
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j)];
      else
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(Nx-1)];
    }
  for(j=1; j<Nx; j++){
    for(i=0; i<Nx; i++){
      g = ( (xiu(tmp1[i+Nx*(j)],tmp2[i+Nx*(j)],tmp3[i+Nx*(j)])-xiu(tmp1[i+Nx*(j-1)],tmp2[i+Nx*(j-1)],tmp3[i+Nx*(j-1)]))/dx - Jusy[i+Nx*j] )*lam;
      if(g>=0.0)
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j)];
      else
        u[i+Nx*j] -= g*tmp1[(i)+Nx*(j-1)];
    }
  }


}





double positive(double x){
  if(x>0.0)
    return x;
  else 
    return 0.0;
}



double calc_Hstab(){
	return F0_11*R0_11*R0_11*R0_11/3.0+C1*R0_11*R0_11/2.0  +  Fr11*R1_11*R1_11*R1_11/3.0+C2*R1_11*R1_11/2.0  -  Fr11*R0_11*R0_11*R0_11/3.0-C2*R0_11*R0_11/2.0  +  C3*R2_11*R2_11/2.0  -  C3*R1_11*R1_11/2.0  +  Fa11/3.0+C4/2.0  -  Fa11*R2_11*R2_11*R2_11/3.0-C4*R2_11*R2_11/2.0;
}

double calc_Hstab2(){
	return Fr11*R1_11*R1_11*R1_11/3.0+C2*R1_11*R1_11/2.0   +  C3*R2_11*R2_11/2.0  -  C3*R1_11*R1_11/2.0  +  Fa11/3.0+C4/2.0  -  Fa11*R2_11*R2_11*R2_11/3.0-C4*R2_11*R2_11/2.0;
}


/*
void init(double u[],double v[]){
  int i, j, k;
  double x, y, tmpx, tmpy, tmp;
  
  for (i = 0; i < Nx; ++i) {
    for(j =0; j<Nx; j++){
      u[i+Nx*j] = 0.0;
      v[i+Nx*j] = 0.0;
    }
  }  

  srand(12);
   for (i = 0; i < Nx; ++i) {
    for (j = 0; j < Nx; ++j) {
		tmpx = 0.5+0.01*(2*(double)rand()/RAND_MAX-1);
		tmpy = 0.5+0.01*(2*(double)rand()/RAND_MAX-1);
      u[i+Nx*j] += tmpx;
      v[i+Nx*j] += tmpy;

    }
  }
}
*/

/*
void init(double u[],double v[]){
  int i, j, k;
  double x, y, tmpx, tmpy, tmp, totalu=0.0, totalv=0.0;
  
  	for (i = 0; i < Nx; ++i) {
    	for (j = 0; j < Nx; ++j) {
			u[(i)+Nx*(j)] = 0;
			v[(i)+Nx*(j)] = 0;
		}
	}

  	fpinit = fopen("init_triangle_restrict1.dat", "r");
//	fpinit = fopen("init_square_restrict4.dat", "r");

	for (i = 0; i < Nx; ++i) {
    	for (j = 0; j < Nx; ++j) {
     	  x=((i+0.5)*dx-L1);
    	  y=((j+0.5)*dx-L1);

//	      fprintf(fpinit,"%e %e %e\n", x, y, u[(i)+Nx*(j)]);
	      if(fscanf(fpinit,"%lf", &tmp1[(i)+Nx*(j)])<0){
			  printf("err init\n");
			  exit(1);
		  };

	    }
  	}

	fclose(fpinit);

	
  	for (i = 0; i < Nx; ++i) {
    	for (j = 0; j < Nx; ++j) {
     		x=((i+0.5)*dx-L1);
    		y=((j+0.5)*dx-L1);
			if(x<-0.8*L1 || x>0.8*L1 || y>0.8*L1 || y<-0.8*L1)
				tmp1[(i)+Nx*(j)] = 0;

		}
	}


  	for (i = 0; i < Nx; ++i) {
    	for (j = 0; j < Nx; ++j) {
     		x=((i+0.5)*dx-L1);
    		y=((j+0.5)*dx-L1);
			if(y>0.15*L1 || (x>0.15*L1 && y>-0.15*L1)){
				u[(i)+Nx*(j)] = tmp1[(i)+Nx*(j)];
				v[(i)+Nx*(j)] = 0;
			}
			else if(x<-0.5*L1 && y<-0.5*L1){
				u[(i)+Nx*(j)] = 0;
				v[(i)+Nx*(j)] = 0;
			}
			else{
				u[(i)+Nx*(j)] = 0;
				v[(i)+Nx*(j)] = tmp1[(i)+Nx*(j)];
			}
		}
	}

}

*/




double X[128][2];

void init(double u[],double v[],double w[],double z[]){
  int i, j, k, l, m, n, nn, IN=5;
  int NN = (3*IN*IN-3*IN+1);
  double random_intensity =0.02;
  double random_intensityr =0.02;
//  double random_intensity =0.4;

  double x, y, tmpx, tmpy, tmp, totalu=0.0, totalv=0.0;
  
  double ix=0.0, iy=0.0;

  srand(1);

  n=0;
  nn=IN;
  for(i=0; i<IN; i++){
		ix = i*0.5*distance;
		iy = -i*sqrt(3.0)/2.0*distance;
		for(j=0; j<nn; j++){
			double randshiftx=0.0, randshifty=0.0;
//			if(n==4 || n==5 || n==8 || n==9 || n==10 || n==13 || n==14){
			  randshiftx = distance*((2*(double)rand()/RAND_MAX-1))*random_intensity; 
			  randshifty = distance*((2*(double)rand()/RAND_MAX-1))*random_intensity; 
//			}
			X[n][1] = ix+j*0.5*distance+ randshiftx;
			X[n][0] = iy+j*sqrt(3.0)/2.0*distance+ randshifty;
			n++;
		}
		nn++;
	}
	nn-=2;
	ix += distance;
  for(i=0; i<IN-1; i++){
		for(j=0; j<nn; j++){
			double randshiftx=0.0, randshifty=0.0;
//			if(n==4 || n==5 || n==8 || n==9 || n==10 || n==13 || n==14){
			  randshiftx = ((2*(double)rand()/RAND_MAX-1))*distance*random_intensity; 
			  randshifty = ((2*(double)rand()/RAND_MAX-1))*distance*random_intensity; 
//			}
			X[n][1] = ix+j*0.5*distance + randshiftx;
			X[n][0] = iy+j*sqrt(3.0)/2.0*distance + randshifty;
			n++;
		}
		nn--;
		ix += distance;
	}
	
    double centerx = 0.0;
    for(i=0; i<NN; i++)
      centerx += X[i][1];
    centerx /= (NN);
    for(i=0; i<NN; i++)
      X[i][1] += -centerx;

    double centery = 0.0;
    for(i=0; i<NN; i++)
      centery += X[i][0];
    centery /= (NN);
    for(i=0; i<NN; i++)
      X[i][0] += -centery;
  

  
  
  for (i = 0; i < Nx; ++i) {
    for(j =0; j<Nx; j++){
      u[i+Nx*j] = 0.0;
      v[i+Nx*j] = 0.0;
      w[i+Nx*j] = 0.0;
      z[i+Nx*j] = 0.0;
    }
  }  


  for(nn=0; nn<NN; nn++){
//if(nn!=22){
		tmpx = X[nn][0];
		tmpy = X[nn][1];

    double randtmpx=((2*(double)rand()/RAND_MAX-1))*distance*random_intensityr, 
    randtmpy=((2*(double)rand()/RAND_MAX-1))*distance*random_intensityr;

    for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Nx; ++j) {
        x=((i+0.5)*dx-L1);
        y=((j+0.5)*dx-L1);
        double rr=sqrt((tmpx-x)*(tmpx-x)+(tmpy-y)*(tmpy-y));
        double theta=atan2((y-tmpy)/rr,(x-tmpx)/rr);
        
        if( rr<init_radius_v ){
          v[i+Nx*j] = .1;//3;
        }
        rr=sqrt((x-tmpx+randtmpx)*(x-tmpx+randtmpx)+(y-tmpy+randtmpy)*(y-tmpy+randtmpy));
        if( rr<(init_radius_u)){     
	        u[i+Nx*j] = .2;
  	      v[i+Nx*j] = .0;
        }
        
//}
        
      }
    }
	}
	
  

   for (i = 0; i < Nx; ++i) {
    for (j = 0; j < Nx; ++j) {
      x=((i+0.5)*dx-L1);
      y=((j+0.5)*dx-L1);
//      if(x*x+y*y <= IN*distance*IN*distance/1.250 && u[i+Nx*j]<=cutoff && v[i+Nx*j]<=cutoff)
      if( u[i+Nx*j]<=cutoff && v[i+Nx*j]<=cutoff)
        w[i+Nx*j] = .040;
    }
  }
  
}




void init2(double u[],double v[],double w[],double z[]){
  int i, j, k, l, m, n, nn, IN=3;
  int NN = (3*IN*IN-3*IN+1);
  double random_intensity =0.08;
  double random_intensityr =0.06;
//  double random_intensity =0.4;

  double x, y, tmpx, tmpy, tmp, totalu=0.0, totalv=0.0;
  
  double ix=0.0, iy=0.0;

  srand(1);

	
  
  for (i = 0; i < Nx; ++i) {
    for(j =0; j<Nx; j++){
     x=((i+0.5)*dx-L1);
      y=((j+0.5)*dx-L1);
      u[i+Nx*j] = 0.0;
      v[i+Nx*j] = 0.0;
      w[i+Nx*j] = 0.0;
      z[i+Nx*j] = 0.0;
    }
  }  

double radi=0.2;
  for (i = 0; i < Nx; ++i) {
    for(j =0; j<Nx; j++){
     x=((i+0.5)*dx-L1);
      y=((j+0.5)*dx-L1);
	if(x*x+y*y<radi*radi)
      u[i+Nx*j] = 0.1;
    }
  }  

  
}


void init3(double u[],double v[],double w[],double z[]){
  int i, j, k, l, m, n, nn, IN=4, JN=8;
  int NN = (2*IN*JN);
  double random_intensity =0.02;
  double random_intensityr =0.02;


  double x, y, tmpx, tmpy, tmp, totalu=0.0, totalv=0.0;
  
  double ix=0.0, iy=0.0;

  srand(1);

  n=0;
  for(i=0; i<IN; i++){
		ix = 0.0;
		iy = -i*distance;
		for(j=0; j<JN; j++){
			double randshiftx=0.0, randshifty=0.0;
			X[n][0] = ix+j*distance+ randshiftx;
			X[n][1] = iy+ randshifty;
			n++;
		}
	}

	double iiy = iy;
  int sign=1;
  for(i=0; i<IN; i++){
		sign *= -1;
		ix += sign*0.5*distance;
		iy = iiy -(i+1)*sqrt(3.0)/2.0*distance;
		for(j=0; j<JN; j++){
			double randshiftx=0.0, randshifty=0.0;
			X[n][0] = ix+j*distance+ randshiftx;
			X[n][1] = iy+ randshifty;
			n++;
		}
	}


    double centerx = 0.0;
    for(i=0; i<NN; i++)
      centerx += X[i][1];
    centerx /= (NN);
    for(i=0; i<NN; i++)
      X[i][1] += -centerx;

    double centery = 0.0;
    for(i=0; i<NN; i++)
      centery += X[i][0];
    centery /= (NN);
    for(i=0; i<NN; i++)
      X[i][0] += -centery;
  

  
  
  for (i = 0; i < Nx; ++i) {
    for(j =0; j<Nx; j++){
      u[i+Nx*j] = 0.0;
      v[i+Nx*j] = 0.0;
      w[i+Nx*j] = 0.0;
      z[i+Nx*j] = 0.0;
    }
  }  


  for(nn=0; nn<NN; nn++){
		tmpx = X[nn][0];
		tmpy = X[nn][1];

    double randtmpx=((2*(double)rand()/RAND_MAX-1))*distance*random_intensityr, 
    randtmpy=((2*(double)rand()/RAND_MAX-1))*distance*random_intensityr;

    for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Nx; ++j) {
        x=((i+0.5)*dx-L1);
        y=((j+0.5)*dx-L1);
        double rr=sqrt((tmpx-x)*(tmpx-x)+(tmpy-y)*(tmpy-y));
        double theta=atan2((y-tmpy)/rr,(x-tmpx)/rr);
        
        rr=sqrt((x-tmpx+randtmpx)*(x-tmpx+randtmpx)+(y-tmpy+randtmpy)*(y-tmpy+randtmpy));
        if( rr<(init_radius_u)){     
	        u[i+Nx*j] = .2;
  	    //  v[i+Nx*j] = .0;
        }
        
       
      }
    }
	}
	
  

  
}



