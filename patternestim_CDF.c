#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define Nnumerinteg 500

#define sqr(x)  ((x) * (x))
#define min(a, b) ((a) < (b) ? (a) : (b))

double dist(double u[], double v[]);
double calc_ip(double u[], double v[], double w[]);
double dist_min(int num);
void angles(double xtmp[][2], int N, double distmin);
int index_angles_4(double spot[][2], int N, double lspot[][2], int lN, double distmin, double alln4angles[]);
int index_angles_6(double spot[][2], int N, double lspot[][2], int lN, double distmin, double alln4angles[]);
int index_angles_n(double spot[][2], int N, double lspot[][2], int lN, double distmin, double alln4angles[]);
double ave_mindist(int num);
double med_mindist(double spot[][2], int num);
int count_dist1(int N, double ave);
int count_distsqrt2(int N, double ave);
int count_bonds(int N, double ave);
double count_surrounding(int N, double ave);

void bubble_sort(double a[], int N);
void bubble_sort_up(double a[], int N);

int get_l_spots(double spot[][2], int N, double meddist, double lspot[][2]);
double l_med_mindist(int num, int lnum);
double l_dist_min(int num, int lnum);
double l_ave_mindist(int num, int lnum);
double l_count_surrounding(int N, int lN, double ave); 
int l_count_bonds(int N, int lN, double ave);
int l_count_dist1(int N, int lN, double ave);
int l_count_distsqrt2(int N, int lN, double ave);
int init(double dist_x, double dist_y, double p_intensity, int randindex, double spot[][2]);


double spot[10000][2], lspot[10000][2], tmps[10000][2];
double spot_G[10000][2], lspot_G[10000][2]; //Generated

double alln4angles[60000], alln6angles[60000], allnangles[60000];
int Ialln4angles = 0, Ialln6angles = 0, Iallnangles = 0;

double alln4angles_G[60000], alln6angles_G[60000], allnangles_G[60000];
int Ialln4angles_G = 0, Ialln6angles_G = 0, Iallnangles_G = 0;

double CDFdata4[Nnumerinteg+1], CDFgenerated4[Nnumerinteg+1], CDFdata6[Nnumerinteg+1], CDFgenerated6[Nnumerinteg+1], CDFdatan[Nnumerinteg+1], CDFgeneratedn[Nnumerinteg+1];

double print_CDF(double a[], double b[]);

double CDF_calc1(double x, double a[], int I);
void CDF_calc2(double a[], int Ia, double cdf[]);
double l2err_CDF(double a[], double b[]);

FILE *fp, *fptr;
char fname[32];



#define Nx 200 
#define Nt 2000
#define Ntheta 70 // Ntheta%4!=0 Ntheta>=4*Nr
#define Nr 20
#define L1 2.80 //domain size
#define Radius 1.0 
double dx = 2.0*(double)L1/((double)Nx);
double u[Nx*Nx], v[Nx*Nx], w[Nx*Nx];
double Fa11, R2_11, F0_22, F0_12, F0_21, x, y;
double min_arearatio_v, ave_arearatio_v;

#define distance 0.61
#define init_radius_u .040
#define init_radius_v .120

int num_spots();
void min_arearatio_v_around_spot(int N, double rad, double *min, double *ave);

int main(void){
    int numspots=0, numspots_G, i, j;
    int numlspots=0, numlspots_G;
    double meddist, meddist_G;
    double err_4, err_6, err_n;

  printf("# Fa_11, R2_11, n4: min(minD,index), minD minP minerr, n6: min(minD,index), minD minP minerr, n: min(minD,index), minD minP minerr\n");

  int cptmp = 500;
	for(int cp22=cptmp; cp22<=cptmp+19; cp22++){
		Fa11 = cp22*1;
    for(int cp=0; cp<=0; cp++){
  	  R2_11 = 0.64;//0.55+0.01*cp;

      sprintf(fname,"data/Fa_11=%.0f_R2_11=%.3f.dat", Fa11, R2_11); 

    	fptr = fopen(fname, "r");
    	if(fptr == NULL) {
		    exit(1);
    	}
  
      for (i = 0; i < Nx; ++i) {
        for (j = 0; j < Nx; ++j) {
          x=((i+0.5)*dx-L1);
          y=((j+0.5)*dx-L1);
          int ste = fscanf(fptr, "%lf %lf %lf", &u[i+Nx*j], &v[i+Nx*j], &w[i+Nx*j]);
        }
      }

      fclose(fptr);

	numspots = num_spots();
	min_arearatio_v_around_spot(numspots, init_radius_v, &min_arearatio_v, &ave_arearatio_v);

  meddist = med_mindist(spot,numspots);
  numlspots = get_l_spots(spot,numspots,meddist,lspot);

  Ialln4angles = index_angles_4(spot, numspots, lspot, numlspots, meddist, alln4angles);
  Ialln6angles = index_angles_6(spot, numspots, lspot, numlspots, meddist, alln6angles);
  Iallnangles = index_angles_n(spot, numspots, lspot, numlspots, meddist, allnangles);

  
  int ND, NP;
  double Dy, P_intensity;
  ND = 1000;
  NP = 400;

  double minerr_4=1.0, minD_4, minP_4;
  double minerr_6=1.0, minD_6, minP_6;
  double minerr_n=1.0, minD_n, minP_n;

  double Dmin=0.5, Dmax=1.0, dD=(Dmax-Dmin)/ND;
  double Pmin=0.0, Pmax=0.4, dP=(Pmax-Pmin)/NP;

  CDF_calc2(alln4angles, Ialln4angles, CDFdata4);
  CDF_calc2(alln6angles, Ialln6angles, CDFdata6);
  CDF_calc2(allnangles, Iallnangles, CDFdatan);


double index=1.0;
	if(numspots<=38 || numspots>=100){
    index = -2.0;
	}
	else if(ave_arearatio_v<0.8){
		index = -1.0;
  }
if(index>0)
for(int nD=0; nD<ND; nD++){
  Dy = Dmin+dD*nD;

for(int nP=0; nP<NP; nP++){
  P_intensity = Pmin+dP*nP;

  numspots_G=init(1.0,Dy,P_intensity,nD,spot_G);

  meddist_G = med_mindist(spot_G,numspots_G);
  numlspots_G = get_l_spots(spot_G,numspots_G,meddist_G,lspot_G);

  Ialln4angles_G = index_angles_4(spot_G, numspots_G, lspot_G, numlspots_G, meddist_G, alln4angles_G);
  Ialln6angles_G = index_angles_6(spot_G, numspots_G, lspot_G, numlspots_G, meddist_G, alln6angles_G);
  Iallnangles_G = index_angles_n(spot_G, numspots_G, lspot_G, numlspots_G, meddist_G, allnangles_G);

  CDF_calc2(alln4angles_G, Ialln4angles_G, CDFgenerated4);
  CDF_calc2(alln6angles_G, Ialln6angles_G, CDFgenerated6);
  CDF_calc2(allnangles_G, Iallnangles_G, CDFgeneratedn);
  
  err_4 = l2err_CDF(CDFdata4, CDFgenerated4); 
  err_6 = l2err_CDF(CDFdata6, CDFgenerated6); 
  err_n = l2err_CDF(CDFdatan, CDFgeneratedn); 

  if(minerr_4>err_4){
    minerr_4 = err_4;
    minD_4 = Dy; minP_4 = P_intensity;
  }
  if(minerr_6>err_6){
    minerr_6 = err_6;
    minD_6 = Dy; minP_6 = P_intensity;
  }
  if(minerr_n>err_n){
    minerr_n = err_n;
    minD_n = Dy; minP_n = P_intensity;
  }

}

}
  
  printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", Fa11, R2_11, min(minD_4,index), minD_4, minP_4, minerr_4, min(minD_6,index), minD_6, minP_6, minerr_6, min(minD_n,index), minD_n, minP_n, minerr_n);

    }
  }
  
    return 0;
}


int get_l_spots(double spot[][2], int N, double meddist, double lspot[][2]){
    int n, i=0;
    for(n=0; n<N; n++){

            	int jj, numl=0;
	double dij;
		  for(jj=0; jj<N; jj++){
			  dij = dist(spot[n],spot[jj]);
			  if(dij > 0.01*(meddist) && dij < (((meddist)*(1.41421356)*1.2)))
				  numl++;
		  }
      if(numl>=6){

            lspot[i][0] = spot[n][0];
            lspot[i][1] = spot[n][1];
            i++;
        }
    }
    return i;
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


void bubble_sort_up(double a[], int N){
  int i, j;
  double tmp;
  
  for(i=0; i<N-1; i++){
    for(j=1; j<N-i; j++){
      if(a[j-1]>a[j]){
        tmp = a[j-1];
        a[j-1] = a[j];
        a[j] = tmp;
      }
    }
  }
}

double med_mindist(double spot[][2], int num){
	int i,j;
	double dij, min[num], med=0;
	
	for(i=0; i<num; i++){
		min[i]=dist(spot[0],spot[1]);
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

double l_med_mindist(int num, int lnum){
	int i,j;
	double dij, min[num], med=0;
	
	for(i=0; i<lnum; i++){
		min[i]=dist(lspot[0],spot[1]);
	    for(j=0; j<num; j++){
				dij = dist(lspot[i],spot[j]);
				if(dij > 0.01 && min[i]>dij)
					min[i]=dij;
		}
	}
	bubble_sort(min,lnum);
	return min[lnum/2];
}

double dist_min(int num){
	int i,j;
	double dij, min=dist(spot[0],spot[1]);
	
	for(i=0; i<num; i++){
	    for(j=i+1; j<num; j++){
			dij = dist(spot[i],spot[j]);
			if(min>dij)
				min=dij;
		}
	}
	return min;
}

double l_dist_min(int num, int lnum){
	int i,j;
	double dij, min=dist(lspot[0],spot[1]);
	
	for(i=0; i<lnum; i++){
	    for(j=0; j<num; j++){
			dij = dist(lspot[i],spot[j]);
			if(dij > 0.01 && min>dij)
				min=dij;
		}
	}
	return min;
}

double ave_mindist(int num){
	int i,j;
	double dij, min, ave=0;
	
	for(i=0; i<num; i++){
		min=dist(spot[0],spot[1]);
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

double l_ave_mindist(int num, int lnum){
	int i,j;
	double dij, min, ave=0;
	
	for(i=0; i<lnum; i++){
		min=dist(lspot[0],lspot[1]);
	    for(j=0; j<num; j++){
				dij = dist(lspot[i],lspot[j]);
				if(dij > 0.01 && min>dij)
					min=dij;
		}
		ave += min;
	}
	return ave/lnum;
}

double count_surrounding(int N, double ave){
	int i,j, counter=0, num=0;
	double dij;
	
    for(i=0; i<N; i++){
          counter++;
		  for(j=0; j<N; j++){
			  dij = dist(spot[i],spot[j]);
			  if(j!=i && dij < (((ave)*(1.41421356)*1.1)))
				  num++;
		  }
      }
      return (double)num/(double)N;
}

double l_count_surrounding(int N, int lN, double ave){
	int i,j, counter=0, num=0;
	double dij;
	
    for(i=0; i<lN; i++){
          counter++;
		  for(j=0; j<N; j++){
			  dij = dist(lspot[i],spot[j]);
			  if(dij > 0.01*ave && dij < (((ave)*(1.41421356)*1.1)))
				  num++;
		  }
      }
      return (double)num/(double)lN;
}

int count_bonds(int N, double ave){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<N; i++){
      for(j=i+1; j<N; j++){
          dij = dist(spot[i],spot[j]);
          if(dij > ((ave)*0.8) && dij < (((ave)*(1.41421356)*1.1)))
          counter++;
        }
      }
      return counter;
}

int l_count_bonds(int N, int lN, double ave){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<lN; i++){
      for(j=0; j<N; j++){
          dij = dist(lspot[i],spot[j]);
          if(dij > ((ave)*0.8) && dij < (((ave)*(1.41421356)*1.1)))
          counter++;
        }
      }
      return counter;
}


int count_dist1(int N, double ave){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<N; i++){
      for(j=i+1; j<N; j++){
          dij = dist(spot[i],spot[j]);
          if(dij > ((ave)*0.9) && dij < ((ave)*1.1))
          counter++;
        }
      }
      return counter;
}

int l_count_dist1(int N, int lN, double ave){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<lN; i++){
      for(j=0; j<N; j++){
          dij = dist(lspot[i],spot[j]);
          if(dij > ((ave)*0.9) && dij < ((ave)*1.1))
          counter++;
        }
      }
      return counter;
}

int count_distsqrt2(int N, double ave){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<N; i++){
      for(j=i+1; j<N; j++){
          dij = dist(spot[i],spot[j]);
          if(dij > (ave)*(1.41421356)*0.9 && dij < ((ave)*(1.41421356)*1.1))
          counter++;
        }
      }
      return counter;
}

int l_count_distsqrt2(int N, int lN, double ave){
	int i,j, counter=0;
	double dij;
	
    for(i=0; i<lN; i++){
      for(j=0; j<N; j++){
          dij = dist(lspot[i],spot[j]);
          if(dij > (ave)*(1.41421356)*0.9 && dij < ((ave)*(1.41421356)*1.1))
          counter++;
        }
      }
      return counter;
}



void angles(double xtmp[][2], int N, double distmin){
  int i, j, k, m;
  double dij, dkj, dik, ip, V;

    for(i=0; i<N; i++){
      for(j=i+1; j<N; j++){
          dij = dist(xtmp[i],xtmp[j]);
        }
      }

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
                  printf("angle: %.2f\n", acos(ip/dij/dik)*180.0/M_PI); 
                }
              }
            }
          
          }
        }
      }  
    }

}


double  get_n4spots(double n4[][2], int I, double spot[][2], int N, double lspot[][2], int lN, double meddist){
  int i,j,k,m;
  double dij, dijm, dkj, dik, ip, tmp0, tmp1, tmph[2], angm, ang;
  int count=0;

   for(j=0; j<N; j++){
      dij = dist(lspot[I],spot[j]);
      if(dij > ((meddist)*0.1) && dij < ((meddist)*2.)){
        tmps[count][0] = spot[j][0];
        tmps[count][1] = spot[j][1];
        count++;
      }
   }
  for(i=0; i<count-1; i++){
    for(j=1; j<count-i; j++){
      dijm = dist(lspot[I],tmps[j-1]);
      dij  = dist(lspot[I],tmps[j]);
      if(dijm>dij){
        tmp0 = tmps[j-1][0];  tmp1 = tmps[j-1][1];
        tmps[j-1][0] = tmps[j][0]; tmps[j-1][1] = tmps[j][1];
        tmps[j][0] = tmp0;  tmps[j][1] = tmp1;
      }
    }
  }

  for(i=0; i<4; i++){
    n4[i][0] = tmps[i][0];
    n4[i][1] = tmps[i][1];    
  }
  tmph[0]=lspot[I][0]+meddist;  tmph[1]=lspot[I][1];
  for(i=0; i<4-1; i++){
    for(j=1; j<4-i; j++){

      angm = atan2(n4[j-1][1]-lspot[I][1], n4[j-1][0]-lspot[I][0]);
      ang = atan2(n4[j][1]-lspot[I][1], n4[j][0]-lspot[I][0]);

      if(angm>ang){
        tmp0 = n4[j-1][0];  tmp1 = n4[j-1][1];
        n4[j-1][0] = n4[j][0]; n4[j-1][1] = n4[j][1];
        n4[j][0] = tmp0;  n4[j][1] = tmp1;
      }
    }
  }

  double meand = 0.0;
  for(i=0; i<4; i++){
    meand += dist(lspot[I],n4[i]);
  }
  meand /= 4.0;

  return meand;
}
 
int index_angles_4(double spot[][2], int N, double lspot[][2], int lN, double distmin, double alln4angles[]){
  int i, j, k, m, I, counter=0;
  double dij, dijp, dik, ip, V, ang, tmp0, ave=0.0, ave01=0.0, ave23=0.0, ave45=0.0;
  double n4[7][2];
  double angle[6];

  int Ialln4angles=0;

    for(I=0; I<lN; I++){
        get_n4spots(n4,I,spot,N,lspot,lN,distmin);
        n4[4][0] = n4[0][0]; n4[4][1] = n4[0][1]; 

        for(i=0; i<4; i++){
          dij  = dist(lspot[I],n4[i]);
          dijp  = dist(lspot[I],n4[i+1]);
          ip = calc_ip(lspot[I],n4[i],n4[i+1]);
          angle[i] = acos(ip/dij/dijp)*180.0/M_PI;

          if(angle[i]>0.01 && angle[i]<180){
            alln4angles[Ialln4angles] = angle[i];
            Ialln4angles++;
          }
        }

        //sorting
        for(i=0; i<4-1; i++){
          for(j=1; j<4-i; j++){
            if(angle[j-1]>angle[j]){
              tmp0 = angle[j-1];
              angle[j-1] = angle[j];
              angle[j] = tmp0;
            }
          }
        }

      ave01 += (angle[0]+angle[1])/2.0;
      ave23 += (angle[2]+angle[3])/2.0;
      ave45 += (angle[4]+angle[5])/2.0;
    }
    ave01 /= lN;
    ave23 /= lN;
    ave45 /= lN;
  bubble_sort_up(alln4angles,Ialln4angles);

  return Ialln4angles;
}


void  get_n6spots(double n6[][2], int I, double spot[][2], int N, double lspot[][2], int lN, double meddist){
  int i,j,k,m;
  double dij, dijm, dkj, dik, ip, tmp0, tmp1, tmph[2], angm, ang;
  int count=0;

   for(j=0; j<N; j++){
      dij = dist(lspot[I],spot[j]);
      if(dij > ((meddist)*0.1) && dij < ((meddist)*2.)){
        tmps[count][0] = spot[j][0];
        tmps[count][1] = spot[j][1];
        count++;
      }
   }
  for(i=0; i<count-1; i++){
    for(j=1; j<count-i; j++){
      dijm = dist(lspot[I],tmps[j-1]);
      dij  = dist(lspot[I],tmps[j]);
      if(dijm>dij){
        tmp0 = tmps[j-1][0];  tmp1 = tmps[j-1][1];
        tmps[j-1][0] = tmps[j][0]; tmps[j-1][1] = tmps[j][1];
        tmps[j][0] = tmp0;  tmps[j][1] = tmp1;
      }
    }
  }

  for(i=0; i<6; i++){
    n6[i][0] = tmps[i][0];
    n6[i][1] = tmps[i][1];    
  }
// bubble sort in angle from positive horizontal line
  tmph[0]=lspot[I][0]+meddist;  tmph[1]=lspot[I][1];
  for(i=0; i<6-1; i++){
    for(j=1; j<6-i; j++){

      angm = atan2(n6[j-1][1]-lspot[I][1], n6[j-1][0]-lspot[I][0]);
      ang = atan2(n6[j][1]-lspot[I][1], n6[j][0]-lspot[I][0]);

      if(angm>ang){
        tmp0 = n6[j-1][0];  tmp1 = n6[j-1][1];
        n6[j-1][0] = n6[j][0]; n6[j-1][1] = n6[j][1];
        n6[j][0] = tmp0;  n6[j][1] = tmp1;
      }
    }
  }

}
 


int index_angles_6(double spot[][2], int N, double lspot[][2], int lN, double distmin, double alln6angles[]){
  int i, j, k, m, I, counter=0;
  double dij, dijp, dik, ip, V, ang, tmp0, ave=0.0, ave01=0.0, ave23=0.0, ave45=0.0;
  double n6[7][2];
  double angle[6];

  int Ialln6angles = 0;

    for(I=0; I<lN; I++){
        get_n6spots(n6,I,spot,N,lspot,lN,distmin);
        n6[6][0] = n6[0][0]; n6[6][1] = n6[0][1]; 

        for(i=0; i<6; i++){
          dij  = dist(lspot[I],n6[i]);
          dijp  = dist(lspot[I],n6[i+1]);
          ip = calc_ip(lspot[I],n6[i],n6[i+1]);
          angle[i] = acos(ip/dij/dijp)*180.0/M_PI;

          if(angle[i]>0.01 && angle[i]<180){
            alln6angles[Ialln6angles] = angle[i];
            Ialln6angles++;
          }
        }

        //sorting
        for(i=0; i<6-1; i++){
          for(j=1; j<6-i; j++){
            if(angle[j-1]>angle[j]){
              tmp0 = angle[j-1];
              angle[j-1] = angle[j];
              angle[j] = tmp0;
            }
          }
        }

      ave01 += (angle[0]+angle[1])/2.0;
      ave23 += (angle[2]+angle[3])/2.0;
      ave45 += (angle[4]+angle[5])/2.0;
    }
    ave01 /= lN;
    ave23 /= lN;
    ave45 /= lN;
  bubble_sort_up(alln6angles,Ialln6angles);

  return Ialln6angles;

}




int get_nspots(double n[][2], int I, double spot[][2], int N, double lspot[][2], int lN, double meddist){
  int i,j,k,m;
  double dij, dijm, dkj, dik, ip, tmp0, tmp1, tmph[2], angm, ang;
  int count=0;
  double mean_n4_d;

  double tmp4[4][2];
  mean_n4_d = get_n4spots(tmp4,I,spot,N,lspot,lN,meddist);

   for(j=0; j<N; j++){
      dij = dist(lspot[I],spot[j]);
      if(dij > ((mean_n4_d)*0.1) && dij < ((mean_n4_d)*(1.41421356)*1.1)){
        n[count][0] = spot[j][0];
        n[count][1] = spot[j][1];
        count++;
      }
   }   

// bubble sort in angle from positive horizontal line
  tmph[0]=lspot[I][0]+mean_n4_d;  tmph[1]=lspot[I][1];
  for(i=0; i<count-1; i++){
    for(j=1; j<count-i; j++){

      angm = atan2(n[j-1][1]-lspot[I][1], n[j-1][0]-lspot[I][0]);
      ang = atan2(n[j][1]-lspot[I][1], n[j][0]-lspot[I][0]);

      if(angm>ang){
        tmp0 = n[j-1][0];  tmp1 = n[j-1][1];
        n[j-1][0] = n[j][0]; n[j-1][1] = n[j][1];
        n[j][0] = tmp0;  n[j][1] = tmp1;
      }
    }
  }
  return count;
}
 


int index_angles_n(double spot[][2], int N, double lspot[][2], int lN, double distmin, double allnangles[]){
  int i, j, k, m, I, counter=0;
  double dij, dijp, dik, ip, V, ang, tmp0, ave=0.0, ave01=0.0, ave23=0.0, ave45=0.0;
  double nspots[32][2];
  int numnspots;
  double angle[32];

  int Iallnangles = 0;
    for(I=0; I<lN; I++){
        numnspots = get_nspots(nspots,I,spot,N,lspot,lN,distmin);
        nspots[numnspots][0] = nspots[0][0]; nspots[numnspots][1] = nspots[0][1]; 

        for(i=0; i<numnspots; i++){
          dij  = dist(lspot[I],nspots[i]);
          dijp  = dist(lspot[I],nspots[i+1]);
          ip = calc_ip(lspot[I],nspots[i],nspots[i+1]);
          angle[i] = acos(ip/dij/dijp)*180.0/M_PI;

          if(angle[i]>0.01 && angle[i]<180){
            allnangles[Iallnangles] = angle[i];
            Iallnangles++;
          }
        }

        //sorting
        for(i=0; i<numnspots-1; i++){
          for(j=1; j<numnspots-i; j++){
            if(angle[j-1]>angle[j]){
              tmp0 = angle[j-1];
              angle[j-1] = angle[j];
              angle[j] = tmp0;
            }
          }
        }

      ave01 += (angle[0]+angle[1])/2.0;
      ave23 += (angle[2]+angle[3])/2.0;
      ave45 += (angle[4]+angle[5])/2.0;
    }
    ave01 /= lN;
    ave23 /= lN;
    ave45 /= lN;
  bubble_sort_up(allnangles,Iallnangles);


  return Iallnangles;
}


double CDFTri(double x){
  if(x<60.0)
    return 0.0;
  return 1.0;
}

double CDFQuad(double x){
  if(x<90.0)
    return 0.0;
  return 1.0;
}

double CDF(double x, double a[], int I){
  int i, tmpi=I-1; 
  if(x<a[0])
    tmpi = 0;
  for(i=0; i<I-1; i++){
    if(a[i]<=x && x<=a[i+1]){
      tmpi = i;
    }
  }
      return (double)(tmpi+1)/I;

  return 1.0;
}

double CDF_calc1(double x, double a[], int I){
  int i, tmpi=I-1; 
  if(x<a[0])
    tmpi = 0;
  for(i=0; i<I-1; i++){
    if(a[i]<=x && x<=a[i+1]){
      tmpi = i;
    }
  }
      return (double)(tmpi+1)/I;

}

void CDF_calc2(double a[], int Ia, double cdf[]){
  int i;
  double maxang=180.0, dx, x;
  dx = maxang/Nnumerinteg;
  
  for(i=0; i<=Nnumerinteg; i++){
    x = i*dx;
    cdf[i] = CDF_calc1(x,a,Ia);
  }
}



double print_CDF(double a[], double b[]){
  int i;
  double err, maxang=180.0, dx, x;
  dx = maxang/Nnumerinteg;
  
  err = 0.0;
  for(i=0; i<=Nnumerinteg; i++){
    x = i*dx;
    printf("%f %f %f\n", x, a[i], b[i]);
  }

  return 0.0;
}

double print_CDF_old(double a[], int Ia, double b[], int Ib, double c[], int Ic){
  int i;
  double err, maxang=180.0, dx, x;
  dx = maxang/Nnumerinteg;
  
  err = 0.0;
  for(i=0; i<Nnumerinteg-1; i++){
    x = i*dx;
    printf("%f %f %f %f %f %f\n", x, CDF(x,a,Ia), CDF(x,b,Ib), CDF(x,c,Ic), CDFTri(x), CDFQuad(x));
  }

  return 0.0;
}


double l2err_CDF(double a[], double b[]){
  int i;
  double err, maxang=M_PI, dx, x;
  dx = maxang/Nnumerinteg;
  
  err = 0.0;
  for(i=0; i<Nnumerinteg; i++){
    x = i*dx;
    err += sqr( a[i]-b[i] );
  }
  err *= dx/maxang;
  err = sqrt(err);

  return err;
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


int init(double dist_x, double dist_y, double p_intensity, int randindex, double spot[][2]){
  int i, j, k, l, m, n, nn, IN=32, JN=32;
  int NN = (IN*JN);

  double x, y, tmpx, tmpy, tmp, totalu=0.0, totalv=0.0;
  double pr, ptheta;
  
  double ix=0.0, iy=0.0;

  srand(randindex);

  n=0;

	double iiy = iy;
  int sign=1;
  for(i=0; i<IN; i++){
		sign *= -1;
		ix += sign*0.5*dist_x;
		iy = iiy +(i)*dist_y;
		for(j=0; j<JN; j++){
      pr = dist_x*p_intensity*(2*(double)rand()/RAND_MAX-1);
      ptheta = (2*(double)rand()/RAND_MAX-1)*M_PI;
			spot[n][0] = ix+j*dist_x  + pr*cos(ptheta);
			spot[n][1] = iy  + pr*sin(ptheta);
			n++;
		}
	}

    double centerx = 0.0;
    for(i=0; i<NN; i++)
      centerx += spot[i][1];
    centerx /= (NN);
    for(i=0; i<NN; i++)
      spot[i][1] += -centerx;

    double centery = 0.0;
    for(i=0; i<NN; i++)
      centery += spot[i][0];
    centery /= (NN);
    for(i=0; i<NN; i++)
      spot[i][0] += -centery;

  return NN;  
}
