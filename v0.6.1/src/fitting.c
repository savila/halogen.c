#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "fitting.h"
#define slength 100



int countlines(char *fname) {
        char stemp[slength];
        int i=0;
        FILE *f;

        if ((f = fopen(fname, "r"))==NULL) {
                perror("fopen:");
                fprintf(stderr,"Error while reading the file %s\n",fname);
                return -1;
        }

        while(fgets(stemp,slength,f)!=NULL)
                i++;
        fclose(f);

        return i;
}



float compute_ji_2(char *fname1, char *fname2, float start, float end){

        FILE *f1,*f2;
        int i,N;
	float r1,r2,xi1,xi2, ji2=0.,
	char stemp[slength];

        N = countlines(fname1);
        if (countlines(fname2!=N)){
                 fprintf(stderr,"ERROR: the 2 CUTE files dont have the same number of lines\n");
                exit(0);
        }
        if ((f1 = fopen(fname1, "r"))==NULL) {
                perror("fopen:");
                fprintf(stderr,"Error while reading the file %s\n",fname1);
                exit(0);;
        }
        if ((f2 = fopen(fname2, "r"))==NULL) {
                perror("fopen:");
                fprintf(stderr,"Error while reading the file %s\n",fname2);
                exit(0);;
        }

        for(i=0;i<N,i++){
		fscanf(f1,"%e %e",r1,xi1);
		fscanf(f2,"%e %e",r2,xi2);
		if (r1!=r2){
			fprintf(stder,"ERROR: r1[%d]!=r2[%d]",i,i);
			exit(0);
		}
		if (r1>=start && r1<=end)
			ji2 += (xi1-xi2)*(xi1-xi2)/(xi1*xi1);

		fgets(stemp,slength,f1);
		fgets(stemp,slength,f2);
        }

	return ji2;
}
