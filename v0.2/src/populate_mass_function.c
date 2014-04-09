#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define LINELENGTH 1024


#include "populate_mass_function.h"

void get_cubic_splines(double *, double *, int , double *);
double spline_inter(double *, double *, int , double *, double );
int read_mass_function(char *, double **, double **, long *,double,long *);
static int compare_double(const void * , const void * );



long populate_mass_function(char *filename, double Mmin, double Lbox, float **halo_masses){
	fprintf(stderr,"populate_mass_function.c V4.0\n");		
	long Npoints,Nshort,Nhalos,i,j;
	double *mass_arr,*mass_inv_arr, *dens_arr, *dens_inv_arr, *y2, d_rand;
	double dens_max;
	long t;


	//Reading cumulative mass function
	fprintf(stderr,"Reading Mass Function... \n");		
	if (read_mass_function(filename,&mass_arr,&dens_arr,&Npoints,Mmin,&Nshort)!=1){
		fprintf(stderr,"ERROR: could not read %s\n",filename);
		return -1;	
	}
	Nshort+=2; //avoid border effects
	fprintf(stderr,"... read\n");		
	t = time(NULL);
	srand(t); //needed? or initiallised in main.c?
	fprintf(stderr,"seed: %ld\n",t);	


	//prepare splines for Cummulative density as a funtion of M: n(M)
	y2 = (double *) calloc(Npoints,sizeof(double));

	get_cubic_splines(mass_arr,dens_arr,Npoints,y2);

	dens_max = spline_inter(mass_arr,dens_arr,Npoints,y2,Mmin);
	Nhalos = (long)(dens_max*Lbox*Lbox*Lbox+0.5);
	

	fprintf(stderr,"Number of halos: %ld\n",Nhalos);	



	fprintf(stderr,"Inverting Mass Function... \n");	
	//prepare splines for inverse function M(n)
	mass_inv_arr = (double *) calloc(Nshort,sizeof(double));
	dens_inv_arr = (double *) calloc(Nshort,sizeof(double)); 
	for (i=0;i<Nshort;i++){
		j = Npoints-1-i;
		mass_inv_arr[i]=mass_arr[j];
		dens_inv_arr[i]=dens_arr[j];
	}
	fprintf(stderr,"...done!\n");
	*halo_masses = (float *) calloc(Nhalos,sizeof(float));
	get_cubic_splines(dens_inv_arr,mass_inv_arr,Nshort,y2);



	//Generate masses (parallelisable)
	fprintf(stderr,"Generating halo masses...\n");	
	#pragma omp parallel for private(i,d_rand) shared(dens_max,halo_masses,Nhalos,dens_inv_arr,mass_inv_arr,Nshort,y2) default(none) //Not sure it wont slow down the non serial access of memory ?
	for (i=0;i<Nhalos;i++){
		d_rand = (double) rand()*dens_max/RAND_MAX;
		(*halo_masses)[i] = (float) spline_inter(dens_inv_arr,mass_inv_arr,Nshort,y2,d_rand);
	}
	fprintf(stderr,"...sorting them...\n");	
	
	qsort(*halo_masses, Nhalos, sizeof(*halo_masses[0]), compare_double);
	

	fprintf(stderr,"...done\n");

	return Nhalos;
}


int read_mass_function(char *fname, double **x, double **y, long *N, double xmin, long *Nred){
	long Nlines=0, N_gtmin=0;
	FILE *f;
	char line[LINELENGTH];
	double X,Y;

 	if ((f = fopen(fname, "r"))==NULL) {
                perror("fopen:");
                fprintf(stderr,"Error while reading the file %s\n",fname);
                return -1;
        }
	while(fgets(line,LINELENGTH,f)!=NULL) {
                if (line[0] != '#') {
			Nlines++;	
		}	
        }
        fclose(f);

	(*x) = (double *) calloc (Nlines,sizeof(double));
	(*y) = (double *) calloc (Nlines,sizeof(double));

	Nlines=0;
	f = fopen (fname,"r");
	while(fgets(line,LINELENGTH,f)!=NULL) {
                if (line[0] != '#') {
			sscanf(line,"%lf %lf",&X,&Y);	
			(*x)[Nlines]=X;
			(*y)[Nlines]=Y;
			Nlines++;
			if (X>xmin)
				N_gtmin++;
		}	
        }
        fclose(f);

	*N = Nlines;
	*Nred = N_gtmin;
	return 1;
}



void get_cubic_splines(double *x, double *y, int N, double *y2){
        int i;
        double *u,s,p;

        u = (double *) calloc(N,sizeof(double));

        y2[0]=0.0; //Natural spline
        u[0]=0.0;
        //fprintf(stderr,"%f\n",y2[0]);


        for (i=1;i<N-1;i++){
                s = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
                p = s * y2[i-1] + 2.0;
                y2[i] = (s - 1.0) / p;
                u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
                u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-s*u[i-1])/p;
                //fprintf(stderr,"%f %f %f\n",s,p,y2[i]);
        }

        y2[N-1]=0.0; //Natural spline

        for (i=N-2;i>=0;i--){
                y2[i] = y2[i]*y2[i+1]+u[i];
                //fprintf(stderr,"%f\n",y2[i]);

        }
        free(u);
}



double spline_inter(double *x, double *y, int N, double *y2, double xi){

        double h,b,a;
        int i_low=0,i_high=N-1,i_mid;

        while(i_high-i_low>1){
                i_mid=(i_high+i_low)/2;
                if (x[i_mid]>xi)
                        i_high = i_mid;
                else
                        i_low = i_mid;
        }

        h = x[i_high] - x[i_low];
        a = (x[i_high]-xi)/h;
        b = (xi-x[i_low])/h;
        return (a*y[i_low]+b*y[i_high]+((a*a*a-a)*y2[i_low]+(b*b*b-b)*y2[i_high])*h*h/6.0);
}

static int compare_double(const void * a, const void * b)
{
  if (*(double*)a > *(double*)b) return -1;
  else if (*(double*)a < *(double*)b) return 1;
  else return 0;  
}
