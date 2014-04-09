#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#include "read_snapshot.h"
#include "populate_mass_function.h"
#include "place_halos.h"


float LUNIT;
float MUNIT;
int SWP;
int LGADGET;
int DGADGET;
float Mmin;


float OVD;
char rho_ref[8];
#define rho_crit (27.755e10)

int Nalpha=0;
int Nlin;
double *alpha;
double *Malpha;

#define LINELENGTH 256

#define NParam 14
#define NParamOpt 5
char ParameterList[NParam][32] = {"Snapshot","GadgetFormat","MassFunctionFile","OutputFile","NCellsLin","alphaFile","rho_ref","Overdensity","MinNumPartPerHalo","GadL_Unit","GadM_Unit","GadSwap","GadDouble","GadLong"};
int ParameterSet[NParam];
int NParametersSet = 0;


char Snapshot[LINELENGTH];
char OutputFile[LINELENGTH],MassFunctionFile[LINELENGTH], alphaFile[LINELENGTH];
int format;



int read_input_file(char *);

int write_halogen_cat(char *, float *, float *, float *, float *, float *, long);

int main(int argc, char **argv){

	fprintf(stderr,"\n*******************************************************************\n");
	fprintf(stderr,"**                                                               **\n");
	fprintf(stderr,"**            =           HALOGEN V0.3          =                **\n");
	fprintf(stderr,"**                                                               **\n");
	fprintf(stderr,"**                                                               **\n");
	fprintf(stderr,"**                                            let there be dark  **\n");
	fprintf(stderr,"*******************************************************************\n\n");

	if (argc!=2){
		fprintf(stderr,"usage: %s inputfile\n",argv[0]);	
		return -1;
	}

	float Lbox, mpart, *x, *y, *z, *hx, *hy, *hz, *hR, om_m;	
	char inname[256];
	long Npart, Nhalos;
	float *HaloMass, rho;

	strcpy(inname,argv[1]);

	fprintf(stderr,"Reading input file...\n");
	if (read_input_file(inname)<0)
		return -1;
	fprintf(stderr,"... file read correctly!\n");


	fprintf(stderr,"Reading Gadget files...\n");
	if (read_snapshot(Snapshot, format, &x, &y, &z, &Npart, &mpart, &Lbox, &om_m)==0)
		fprintf(stderr,"Gadget file(s) correctly read!\n");
	else {
		fprintf(stderr,"error: Something went wrong reading the gadget file %s\n",inname);
		return -1;
	}
	fprintf(stderr,"Check: Npart=%ld, mpart=%e, Lbox=%f\n",Npart,mpart,Lbox);
	fprintf(stderr,"x[0]= %f, y[0]= %f, z[0]= %f\n",x[0],y[0],z[0]);
	fprintf(stderr,"      ...\n");
	fprintf(stderr,"x[%ld]= %f, y[%ld]= %f, z[%ld]= %f\n",Npart-1,x[Npart-1],Npart-1,y[Npart-1],Npart-1,z[Npart-1]);


	Nhalos = populate_mass_function(MassFunctionFile,Mmin*mpart,Lbox,&HaloMass);
	if (Nhalos<0)
		fprintf(stderr,"error: Couldnt create HaloMass array\n");	

	hx = (float *) calloc(Nhalos,sizeof(float));
	hy = (float *) calloc(Nhalos,sizeof(float));
	hz = (float *) calloc(Nhalos,sizeof(float));
	hR = (float *) calloc(Nhalos,sizeof(float));

	if (strcmp(rho_ref,"crit")==0)
		rho = OVD*rho_crit;
	else 
		rho = OVD*rho_crit*om_m;

	if (place_halos(Nhalos, HaloMass, Nlin, Npart, x, y, z, Lbox, rho,-1,mpart, alpha, Malpha, Nalpha, hx, hy, hz,hR)==0)
		fprintf(stderr,"Halos placed!\n");
	else
		fprintf(stderr,"error in placing the halos\n");
	
	write_halogen_cat(OutputFile,hx,hy,hz,HaloMass,hR,Nhalos);
	return 0;	
}

int write_halogen_cat(char *filename, float *x, float *y, float *z, float *M, float *R,long N){
	FILE *f;
	long i;
	if ((f=fopen(filename,"w") )== NULL){
		fprintf(stderr,"Couldnt open output file %s\n",filename);
		return -1;
	}
	for(i=0;i<N;i++){
		fprintf(f,"%f %f %f %e %f\n",x[i],y[i],z[i],M[i],R[i]);
	}
	fclose(f);
	return 0;
}


int read_input_file(char *name){
	FILE *f;
	char line[LINELENGTH],word[LINELENGTH];
	int i;
	float x,y;
	for (i=0;i<NParam;i++) 
		ParameterSet[i]=0;	

	if ((f = fopen(name,"r"))==NULL){
		fprintf(stderr,"Could not open input file %s\n",name);
		return -1;
	}
        while(fgets(line,LINELENGTH,f)!=NULL) {	
                if (line[0] != '#') {				//Omit comment lines
			if (sscanf(line,"%s",word)!=EOF){ 	//Check it is not an empty line
			  for (i=0;i<NParam;i++)		//Find the parameter specified
				if (strcmp(ParameterList[i],word)==0)
					break;
			  switch (i){	
				case 0:
					sscanf(line,"%s",Snapshot);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;
				case 1:
					sscanf(line,"%d",&format);
					if (format!=1 && format !=2){
						fprintf(stderr,"ERROR: invalid format: %d. Please select 1 or 2\n",format);
						return -1;
					}	
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;
					
				case 2:
					sscanf(line,"%s",MassFunctionFile);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 3:
					sscanf(line,"%s",OutputFile);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 4:
					sscanf(line,"%d",&Nlin);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;
				case 5:
					sscanf(line,"%s",alphaFile);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 6:
					sscanf(line,"%s",rho_ref);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					if ((strcmp(rho_ref,"crit")!=0) && (strcmp(rho_ref,"matter")!=0)){
							fprintf(stderr,"ERROR: Not valid option for %s: %s.\nPlease select \"crit\" or \"matter\". \n",word,rho_ref);
							return -1;
					}
					ParameterSet[i]++;
					break;

				case 7:
					sscanf(line,"%f",&OVD);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 8: 
					sscanf(line,"%f",&Mmin);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 9:
					sscanf(line,"%f",&LUNIT);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 10:
					sscanf(line,"%f",&MUNIT);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 11:
					sscanf(line,"%d",&SWP);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 12:
					sscanf(line,"%d",&DGADGET);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 13:
					sscanf(line,"%d",&LGADGET);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				default:
					fprintf(stderr,"WARNING: Unknown parameter %s\n",word);
					break;
			}//switch
		  }// if (non-empty line)
		}//if (no comment line)
	}//while	
	fclose(f);

	if(NParametersSet!=NParam){
		for (i=0;i<NParam;i++){
			if (ParameterSet[i]==0)
				fprintf(stderr,"ERROR: Parameter %s not set in input file\n",ParameterList[i]);
		}
		return -1;
	}
	
	if ((f = fopen(alphaFile,"r"))==NULL){
		fprintf(stderr,"Could not open input file %s\n",alphaFile);
		return -1;
	}
        while(fgets(line,LINELENGTH,f)!=NULL)
                if (line[0] != '#')
			Nalpha++;	
	fclose(f);
	
	alpha = (double *) calloc(Nalpha,sizeof(double));
	Malpha = (double *) calloc(Nalpha,sizeof(double));
	
	if ((f = fopen(alphaFile,"r"))==NULL){
		fprintf(stderr,"Could not open alpha file %s\n",alphaFile);
		return -1;
	}

        while(fgets(line,LINELENGTH,f)!=NULL){
                if (line[0] != '#'){
			fgets(line,LINELENGTH,f);
			sscanf(line,"%f %f",&x,&y);
			alpha[i]=x;
			Malpha[i]=y;
		}
	}
	fclose(f);
	return 0;
}

