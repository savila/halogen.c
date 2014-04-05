#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "read_snapshot.h"
#include "populate_mass_function.h"
#include "place_halos.h"


#define LUNIT 1.
#define MUNIT 1.0e10
#define SWP 0
#define LGADGET 0
#define DGADGET 0
#define NMIN 20
#define MMAX 1.0e16 //unused

#define OVD (200.)
#define rho_crit (27.755e10)

#define Nalpha 1
int Nlin = 64;
double alpha[Nalpha]={2.0};
double Malpha[Nalpha];



int write_halogen_cat(char *, float *, float *, float *, double *, long);

int main(int argc, char **argv){
	if (argc!=5){
		fprintf(stderr,"usage: %s gadget_file gadget_format mass_function_file output_file\n",argv[0]);	
		return -1;
	}
	int format;
	float Lbox, mpart, *x, *y, *z, *hx, *hy, *hz;	
	char inname[256], outname[256], mf_file[256];
	long Npart, Nhalos;
	double *HaloMass;

	strcpy(inname,argv[1]);
	format = atoi(argv[2]);
	strcpy(mf_file,argv[3]);
	strcpy(outname,argv[4]);

	fprintf(stderr,"Reading Gadget files...\n");
	if (read_snapshot(inname, LUNIT, MUNIT, format,SWP, LGADGET, DGADGET, &x, &y, &z, &Npart, &mpart, &Lbox)==0)
		fprintf(stderr,"Gadget file(s) correctly read!\n");
	else {
		fprintf(stderr,"error: Something went wrong reading the gadget file %s\n",inname);
		return -1;
	}
	fprintf(stderr,"Check: Npart=%ld, mpart=%e, Lbox=%f\n",Npart,mpart,Lbox);
	fprintf(stderr,"x[0]= %f, y[0]= %f, z[0]= %f\n",x[0],y[0],z[0]);
	fprintf(stderr,"      ...\n");
	fprintf(stderr,"x[%ld]= %f, y[%ld]= %f, z[%ld]= %f\n",Npart-1,x[Npart-1],Npart-1,y[Npart-1],Npart-1,z[Npart-1]);


	Nhalos = populate_mass_function(mf_file,NMIN*mpart,MMAX,Lbox,&HaloMass);
	if (Nhalos<0)
		fprintf(stderr,"error: Couldnt create HaloMass array\n");	

	hx = (float *) calloc(Nhalos,sizeof(float));
	hy = (float *) calloc(Nhalos,sizeof(float));
	hz = (float *) calloc(Nhalos,sizeof(float));


	Malpha[Nalpha-1]=NMIN*mpart*0.99;
	if (place_halos(Nhalos, HaloMass, Nlin, Npart, x, y, z, Lbox, rho_crit*OVD,mpart, alpha, Malpha, Nalpha, hx, hy, hz)==0)
		fprintf(stderr,"Halos placed!\n");
	else
		fprintf(stderr,"error in placing the halos\n");

	
	write_halogen_cat(outname,hx,hy,hz,HaloMass,Nhalos);
	return 0;	
}

int write_halogen_cat(char *filename, float *x, float *y, float *z, double *M, long N){
	FILE *f;
	long i;
	if ((f=fopen(filename,"w") )== NULL){
		fprintf(stderr,"Couldnt open output file %s\n",filename);
		return -1;
	}
	for(i=0;i<N;i++){
		fprintf(f,"%f %f %f %e\n",x[i],y[i],z[i],M[i]);
	}
	fclose(f);
	return 0;
}
