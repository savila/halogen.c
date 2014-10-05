/*=============================================================================
 *                              LIBRARIES
 *=============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "read_snapshot.h"
#include "populate_mass_function.h"

#include "place_halos.h"


/*=============================================================================
 *                              COMMON DEFINES
 *=============================================================================*/

#define rho_crit (27.755e10)
#define LINELENGTH 256
#define NParam 18

char ParameterList[NParam][32] = {"Snapshot","GadgetFormat","MassFunctionFile",
		"OutputFile","NCellsLin","alphaFile","rho_ref","Overdensity","Mmin",
		"GadL_Unit","GadM_Unit","GadSwap","GadDouble","GadLong","Seed","recalc_frac",
		"nthreads","Beta"};


int ParameterSet[NParam];
int NParametersSet = 0;
long seed;
int nthreads;

char Snapshot[LINELENGTH];
char OutputFile[LINELENGTH],MassFunctionFile[LINELENGTH], alphaFile[LINELENGTH];
int format;
float recalc_frac;

int Nlin;
float Mmin;
float OVD;
char rho_ref[8];

int Nalpha=0;
double *alpha,*betavec;
double *Malpha;

float LUNIT, MUNIT;
int SWP, LGADGET, DGADGET;



/*=============================================================================
 *                              PROTOTYPES
 *=============================================================================*/

int read_input_file(char *);

int write_halogen_cat(char *, float *, float *, float *, float *, float *, float *, float *, float *, long);




/*=============================================================================
 *                              MAIN()
 *=============================================================================*/

int main(int argc, char **argv){

	fprintf(stderr,"\n*******************************************************************\n");
	fprintf(stderr,"**                                                               **\n");
	fprintf(stderr,"**            =         HALOGEN V0.7.3          =                **\n");
	fprintf(stderr,"**                                                               **\n");
	fprintf(stderr,"**                                                               **\n");
	fprintf(stderr,"**                                            let there be dark  **\n");
	fprintf(stderr,"*******************************************************************\n\n");

	if (argc!=2){
		fprintf(stderr,"usage: %s inputfile\n",argv[0]);	
		return -1;
	}

	float Lbox, mpart, *x, *y, *z, *vx,*vy,*vz,*hx, *hy, *hz, *hvx,*hvy,*hvz,*hR, om_m;	
	char inname[256];
	long Npart, Nhalos, **ListOfParticles, *NPartPerCell;
	float *HaloMass, rho;
	double *MassLeft;
	float *dens;
	long i;

	strcpy(inname,argv[1]);

#ifdef VERB
	fprintf(stderr,"#def VERB\n");
#endif 
#ifdef DEBUG
	fprintf(stderr,"#def DEBUG \n");
#endif
#ifdef ULTRADEBUG
	fprintf(stderr,"#def ULTRADEBUG \n");
#endif
#ifdef ONLYBIG
	fprintf(stderr,"#def ONLYBIG\n");
#endif
#ifdef NO_EXCLUSION
	fprintf(stderr,"#def NO_EXCLUSION\n");
#endif
#ifdef NO_MASS_CONSERVATION
	fprintf(stderr,"#def NO_MASS_CONSERVATION\n");
#endif
#ifdef MASS_OF_PARTS
	fprintf(stderr,"#def MASS_OF_PARTS\n");
#endif
#ifdef RANKED
	fprintf(stderr,"#def RANKED\n");
#endif
#ifdef NDENS
	fprintf(stderr,"#def NDENS\n");
#endif
#ifdef NGP
	fprintf(stderr,"#def NGP\n");
#endif


	fprintf(stderr,"\nReading input file...\n");
	if (read_input_file(inname)<0)
		return -1;
	fprintf(stderr,"... file read correctly!\n");


	fprintf(stderr,"Reading Gadget file(s)...\n");
	if (read_snapshot(Snapshot, format, LUNIT, MUNIT, SWP, LGADGET, DGADGET,Nlin,nthreads,&x, &y, &z, &vx, &vy, &vz, &Npart, &mpart, &Lbox, &om_m,&ListOfParticles,&NPartPerCell,&MassLeft,&dens)==0)
		fprintf(stderr,"Gadget file(s) correctly read!\n");
	else {
		fprintf(stderr,"error: Something went wrong reading the gadget file %s\n",inname);
		return -1;
	}
	
	#ifdef VERB
	fprintf(stderr,"\n\tCheck: Npart=%ld, mpart=%e, Lbox=%f\n",Npart,mpart,Lbox);
	fprintf(stderr,"\tx[0]= %f, y[0]= %f, z[0]= %f\n",(x)[0],(y)[0],(z)[0]);
	fprintf(stderr,"\t      ...\n");
	fprintf(stderr,"\tx[%ld]= %f, y[%ld]= %f, z[%ld]= %f\n\n",Npart-1,(x)[Npart-1],Npart-1,(y)[Npart-1],Npart-1,(z)[Npart-1]);
	#endif
	
	//TEMPORARY OUTPUT FOR TESTING AGAINST 0.6.7
	for(i=0;i<10;i++){
		fprintf(stderr,"ListOfParticles[%ld][0] = %ld, NPartPerCell[%ld] = %ld\n",i,ListOfParticles[i][0],i,NPartPerCell[i]);
	}

	if (seed<0){
		seed = time(NULL);
		fprintf(stderr,"Seed used: %ld\n",seed);
	}

	//Generate the halo masses from the mass function
	fprintf(stderr,"Generating Halo Masses...\n");
	Nhalos = populate_mass_function(MassFunctionFile,Mmin,Lbox,&HaloMass,seed,nthreads);
	if (Nhalos<0)
		fprintf(stderr,"error: Couldnt create HaloMass array\n");	
	fprintf(stderr,"...Halo Masses Generated\n");

	// TEMPORARY OUTPUT FOR TESTING AGAINST 0.6.7
	for(i=0;i<10;i++){
		fprintf(stderr,"HaloMass[%ld] = %f\n",i,HaloMass[i]);
	}
	for(i=Nhalos-10;i<Nhalos;i++){
		fprintf(stderr,"HaloMass[%ld] = %f\n",i,HaloMass[i]);
	}

	//Allocalte memory for the halo XYZR, and MassLeft vector
	hx = (float *) calloc(Nhalos,sizeof(float));
	hy = (float *) calloc(Nhalos,sizeof(float));
	hz = (float *) calloc(Nhalos,sizeof(float));
	hvx = (float *) calloc(Nhalos,sizeof(float));
	hvy = (float *) calloc(Nhalos,sizeof(float));
	hvz = (float *) calloc(Nhalos,sizeof(float));
	hR = (float *) calloc(Nhalos,sizeof(float));
	
	//MassLeft = (double *) calloc(Nlin*Nlin*Nlin,sizeof(double));
	
	//density at the boundary of a halo
	if (strcmp(rho_ref,"crit")==0)
		rho = OVD*rho_crit;
	else 
		rho = OVD*rho_crit*om_m;

	//place the halos
	fprintf(stderr,"Placing halos down...\n");
	
	// Check the M-alpha vector against produced halos.
	fprintf(stderr,"MALPHA, HALOMASS: %e, %e", Malpha[Nalpha-1],HaloMass[Nhalos-1]);
	if( Malpha[Nalpha-1] > HaloMass[Nhalos-1]){
		if (Malpha[Nalpha-1] < 1.2*HaloMass[Nhalos-1]) {
			Malpha[Nalpha-1] = HaloMass[Nhalos-1]/1.01;
		}
		else{
			fprintf(stderr,"ERROR: LOWEST M(alpha) LARGER THAN LOWEST MASS: %e < %e",Malpha[Nalpha-1],HaloMass[Nhalos-1]);
			exit(0);
		}
	}
	fprintf(stderr,"MALPHA, HALOMASS: %e, %e", Malpha[Nalpha-1],HaloMass[Nhalos-1]);

	fprintf(stderr,"\n");
	if (place_halos(Nhalos,HaloMass, Nlin, Npart, x, y, z, vx,vy,vz,Lbox, rho,seed,mpart, nthreads,alpha, betavec, Malpha, Nalpha,recalc_frac,hx, hy, hz, hvx,hvy,hvz, hR,ListOfParticles,NPartPerCell,MassLeft,dens)==0){
		fprintf(stderr,"...halos placed correctly\n");
	}


	/*
	int place_halos(long Nend, float *HaloMass, long Nlin, long NTotPart, float *PartX,
                float *PartY, float *PartZ, float *PartVX, float *PartVY, float *PartVZ,
                float L, float rho_ref, long seed, float mp, int nthreads, double *alpha, double *beta, double *Malpha,
                long Nalpha,int Nrecalc, float *HaloX, float *HaloY, float *HaloZ, float *HaloVX,
                float *HaloVY, float *HaloVZ,float *HaloR,long **ListOfPart,
                long *NPartPerCell,double *RemMass, float *dens){
*/


	else {
		fprintf(stderr,"Problem placing halos\n");
		return -1;
	}
	fprintf(stderr,"\n");

	for(i=0;i<10;i++){
		fprintf(stderr,"(hx, hy, hz)[%ld] = (%f, %f, %f)",i,hx[i],hy[i],hz[i]);
	}
	for(i=Nhalos-10;i<Nhalos;i++){
		fprintf(stderr,"(hx, hy, hz)[%ld] = (%f, %f, %f)",i,hx[i],hy[i],hz[i]);
	}
	//writting output	
	fprintf(stderr,"Writing Halo catalogue...\n");
	write_halogen_cat(OutputFile,hx,hy,hz,hvx,hvy,hvz,HaloMass,hR,Nhalos);
	fprintf(stderr,"...halo catalogue written in %s\n",OutputFile);
	
	free(hx);free(hy);free(hz);free(hR);
	free(alpha); free(Malpha);


	fprintf(stderr,"\n*******************************************************************\n");
	fprintf(stderr,"**                        ... and there were dark matter haloes  **\n");
	fprintf(stderr,"*******************************************************************\n\n");

	return 0;	
}


/*=============================================================================
 *                              I/O
 *=============================================================================*/

int write_halogen_cat(char *filename, float *x, float *y, float *z, float *vx, float *vy, float *vz, float *M, float *R,long N){
	FILE *f;
	long i;
	if ((f=fopen(filename,"w") )== NULL){
		fprintf(stderr,"Couldnt open output file %s\n",filename);
		return -1;
	}
	for(i=0;i<N;i++){
		#ifndef FITTING
		fprintf(f,"%f %f %f %f %f %f %e %f\n",x[i],y[i],z[i],vx[i],vy[i],vz[i],M[i],R[i]);
		#else
		fprintf(f,"%f %f %f %f %f %f %e %f\n",x[i],y[i],z[i],vx[i],vy[i],vz[i],M[i],R[i]);
		#endif
	}
	fclose(f);
	return 0;
}


int read_input_file(char *name){
	FILE *f;
	char line[LINELENGTH],word[LINELENGTH];
	int i;
	float x,y;
	float beta;
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
					sscanf(line,"%s %s",word,Snapshot);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;	
					#ifdef VERB 
					fprintf(stderr,"\t%s: %s\n",word,Snapshot);
					#endif
					ParameterSet[i]++;
					break;
				case 1:
					sscanf(line,"%s %d",word,&format);
					if (format!=1 && format !=2){
						fprintf(stderr,"ERROR: invalid parameter %s: %d. Please select 1 or 2\n",word,format);
						return -1;
					}	
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;	
					#ifdef VERB 
					fprintf(stderr,"\t%s: %d\n",word,format);
					#endif
					ParameterSet[i]++;
					break;
					
				case 2:
					sscanf(line,"%s %s",word,MassFunctionFile);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %s\n",word,MassFunctionFile);
					#endif
					ParameterSet[i]++;
					break;

				case 3:
					sscanf(line,"%s %s",word,OutputFile);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %s\n",word,OutputFile);
					#endif
					ParameterSet[i]++;
					break;

				case 4:
					sscanf(line,"%s %d",word,&Nlin);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %d\n",word,Nlin);
					#endif
					ParameterSet[i]++;
					break;
				case 5:
					sscanf(line,"%s %s",word,alphaFile);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %s\n",word,alphaFile);
					#endif
					ParameterSet[i]++;
					break;

				case 6:
					sscanf(line,"%s %s",word,rho_ref);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %s\n",word,rho_ref);
					#endif
					if ((strcmp(rho_ref,"crit")!=0) && (strcmp(rho_ref,"matter")!=0)){
							fprintf(stderr,"ERROR: Not valid option for %s: %s.\nPlease select \"crit\" or \"matter\". \n",word,rho_ref);
							return -1;
					}
					ParameterSet[i]++;
					break;

				case 7:
					sscanf(line,"%s %f",word,&OVD);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %f\n",word,OVD);
					#endif
					ParameterSet[i]++;
					break;

				case 8: 
					sscanf(line,"%s %f",word,&Mmin);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %f.\n",word,Mmin);
					#endif
					ParameterSet[i]++;
					break;

				case 9:
					sscanf(line,"%s %f",word,&LUNIT);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					ParameterSet[i]++;
					break;

				case 10:
					sscanf(line,"%s %f",word,&MUNIT);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %f\n",word,MUNIT);
					#endif
					ParameterSet[i]++;
					break;

				case 11:
					sscanf(line,"%s %d",word, &SWP);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %d\n",word,SWP);
					#endif
					ParameterSet[i]++;
					break;

				case 12:
					sscanf(line,"%s %d",word,&DGADGET);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB 
					fprintf(stderr,"\t%s: %d\n",word,DGADGET);
					#endif
					ParameterSet[i]++;
					break;

				case 13:
					sscanf(line,"%s %d",word,&LGADGET);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %d\n",word,LGADGET);
					#endif
					ParameterSet[i]++;
					break;
				case 14:
					sscanf(line,"%s %ld",word,&seed);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %ld\n",word,seed);
					#endif
					ParameterSet[i]++;
					break;
				case 15:
					sscanf(line,"%s %f",word,&recalc_frac);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %f\n",word,recalc_frac);
					#endif
					ParameterSet[i]++;
					break;
				case 16:
					sscanf(line,"%s %d",word,&nthreads);
					if (ParameterSet[i]>0)
						fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
					else
						NParametersSet++;
					#ifdef VERB
					fprintf(stderr,"\t%s: %d\n",word,nthreads);
					#endif
					ParameterSet[i]++;
					break;
                                case 17:
                                        sscanf(line,"%s %f",word,&beta);
                                        if (ParameterSet[i]>0)
                                                fprintf(stderr,"WARNING: Parameter %s set more than once\n",word);
                                        else
                                                NParametersSet++;
                                        #ifdef VERB
                                        fprintf(stderr,"\t%s: %d\n",word,nthreads);
                                        #endif
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
	betavec = (double *) calloc(Nalpha,sizeof(double));
	Malpha = (double *) calloc(Nalpha,sizeof(double));
	
	if ((f = fopen(alphaFile,"r"))==NULL){
		fprintf(stderr,"Could not open alpha file %s\n",alphaFile);
		return -1;
	}

	i=0;
        while(fgets(line,LINELENGTH,f)!=NULL){
                if (line[0] != '#'){
			//fgets(line,LINELENGTH,f);
			sscanf(line,"%f %f",&x,&y);
			alpha[i]=(double) x;
			if (beta==-1.0)
				betavec[i]=x;
			else
				betavec[i]=beta;
			Malpha[i]= (double) y;
			i++;
		}
	}
	fclose(f);
	return 0;
}

