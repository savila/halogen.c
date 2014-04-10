#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> 
#include <omp.h> 

#include "place_halos.h"
 

#define MAXTRIALS (20)




//long select_cell_rnd(long *, long *, long *); 
long select_cell(); 
long select_heaviest_cell(long *, long *, long *); 
long select_part(long );
int check_HaloR_in_mesh(long,float *, float *, float * , float *,long,long,long);
int check_HaloR_in_cell(long ,float *, float *, float * , float *,long ,long,long);
void ComputeCumulative(double);


long **ListOfPart, **ListOfHalos, *NPartPerCell, *NHalosPerCell;
double *CumulativeProb; 
double *MassLeft;
float Lbox,TotProb;
long NCells,NTotCells;
double mpart;




float square(float a){
	return a*a;
}


float R_from_mass(float Mass,float rho) {
	return  (float) pow((3./(4.*rho*M_PI)*Mass),(1./3.));
}

long check_limit(long i, long N){
	if (i==N)
		return 0;
	if (i<0 || i>N){
		fprintf(stderr,"particle assigned to cell %ld\nExiting...",i);
		exit(0);
	}
	return i;
		
}


int place_halos(long NHalosTot, float *HaloMass, long Nlin, long NTotPart, float *PartX, float *PartY, float *PartZ, float L, float rho_ref, long seed, float mp, double *alpha, double *Malpha,long Nalpha,float *HaloX, float *HaloY, float *HaloZ, float *HaloR){

fprintf(stderr,"\tThis is place_halos.c v9.2\n");


//Initiallising -------------------------------------------------
	long i,j,k,lin_ijk,check, icell, Nmin;
	long *count,trials;
	long ihalo,ilong, ipart, Nhalos,i_alpha;
	double invL = 1./L, ProbDiff;
	float Mcell,Mhalo,Mchange,exp; 
	float R;
	time_t t0;

	#ifdef VERB
	time_t t1,t2,t3,t4,t5;
	float diff;
	#endif

	NCells = Nlin;
	Lbox = L;
	
	t0=time(NULL);
	NTotCells = NCells*NCells*NCells;
	
	//Allocate memory for the arrays 
	NPartPerCell = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	NHalosPerCell = (long *) calloc(NCells*NCells*NCells,sizeof(long));
	MassLeft = (double *) calloc(NCells*NCells*NCells,sizeof(double));
	count = (long *) calloc(NCells*NCells*NCells,sizeof(long));


	CumulativeProb = (double *) calloc(NTotCells, sizeof(double));
  	if(CumulativeProb == NULL) {
    		fprintf(stderr,"\tplace_halos(): could not allocate %ld  array for CumulativeProb[]\nABORTING",NTotCells);
    		exit(-1);
	}
	fprintf(stderr,"\tUsing OMP with %d threads\n",omp_get_max_threads());

        //Initiallise random numbers
	#ifdef VERB
        fprintf(stderr,"\tinput seed: %ld.    time0: %ld. Used: ",seed,t0);
	#endif

        if (seed>=0){
                srand(seed);
		#ifdef VERB
        	fprintf(stderr,"%ld \n",seed);
		#endif
	}
        else {
                srand(t0);
		#ifdef VERB
        	fprintf(stderr,"%ld \n",t0);
		#endif
	}

	mpart = (double) mp;

	Nmin = (long)ceil(HaloMass[NHalosTot-1]/mpart);

	
	#ifdef VERB
	fprintf(stderr,"\n\tParticles and Halos placed in %ld^3 cells\n",NCells);
	fprintf(stderr,"\tBOX = %f  lcell =%f   rho_ref = %e  invL %f\n",L,L/NCells,rho_ref,invL);
	fprintf(stderr,"\tNhalos = %ld NPart = %ld\n",NHalosTot, NTotPart);
	#endif
	#ifdef DEBUG
	fprintf(stderr,"\n\tRAND_MAX=%d\n",RAND_MAX);
	fprintf(stderr,"\tX[0] = %f Y[0] = %f Z[0] = %f\n",PartX[0],PartY[0],PartZ[0]);
	fprintf(stderr,"\tX[1] = %f Y[1] = %f Z[1] = %f\n",PartX[1],PartY[1],PartZ[1]);
	fprintf(stderr,"\tM[0] = %e \n",HaloMass[0]);
	fprintf(stderr,"\tM[1] = %e \n",HaloMass[1]);
	fprintf(stderr,"\tM[%ld] = %e \n",NHalosTot-1,HaloMass[NHalosTot-1]);
	fprintf(stderr,"\n\tMinimmum mass= %e. Minimum part per halo = %ld. mpart %e\n",HaloMass[NHalosTot-1],Nmin,mpart);
	#endif	

	
	#ifdef VERB
	t1=time(NULL);
 	diff = difftime(t1,t0);
	fprintf(stderr,"\ttime of initialisation %f\n",diff);
	#endif
// ------------------------------------------------- Initiallised



#ifdef VERB
	fprintf(stderr,"\tAssigning particles to grid ...\n");
#endif


//Assign particles to grid ------------------------------------
	//count particles per cell
	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (long) (invL * PartX[ilong]*NCells);
		j = (long) (invL * PartY[ilong]*NCells);
		k = (long) (invL * PartZ[ilong]*NCells);
		if (i<0 || i>=NCells || j<0 || j>=NCells || k<0 || k>=NCells){	
			fprintf(stderr,"\tWARNING: Particle %ld at [%f,%f,%f] seems to be out of the right box interval [0.,%f)",ilong,PartX[ilong],PartY[ilong],PartZ[ilong],L);	
			i=check_limit(i,NCells);
			j=check_limit(j,NCells);
			k=check_limit(k,NCells);
			fprintf(stderr,", placed at cell [%ld,%ld,%ld]\n",i,j,k);
		}
		lin_ijk = k+j*NCells+i*NCells*NCells;
		NPartPerCell[lin_ijk]++;
#ifdef DEBUG
		if(ilong<10 || ilong > NTotPart -10 || ilong==243666)
			fprintf(stderr,"\tipart=%ld  cell: %ld=[%ld,%ld,%ld] Parts in cell=%ld, Pos= [%f,%f,%f]\n",ilong,lin_ijk,i,j,k,NPartPerCell[lin_ijk],PartX[ilong],PartY[ilong],PartZ[ilong]);
#endif
	}
#ifdef VERB
	fprintf(stderr,"\t... particles counted ...\n");
	t2=time(NULL);
 	diff = difftime(t2,t1);
	fprintf(stderr,"\ttime counting %f\n",diff);
#endif
	//Alloc Enough Memory
	ListOfPart = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	ListOfHalos = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
	for (i=0;i<NCells;i++){
	for (j=0;j<NCells;j++){
	for (k=0;k<NCells;k++){
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk] = (long *) calloc(NPartPerCell[lin_ijk],sizeof(long));
		Nhalos = (long) (NPartPerCell[lin_ijk]/Nmin);
		ListOfHalos[lin_ijk] = (long *) calloc(Nhalos,sizeof(long));
		MassLeft[lin_ijk] = (double) NPartPerCell[lin_ijk]*mpart; 
#ifdef ULTRADEBUG
		if (lin_ijk<10 || lin_ijk > (NCells*NCells*NCells) - 10){
			fprintf(stderr,"\tAllocated %ld (longs) in ListOfPart(%ld=[%ld,%ld,%ld])\n",NPartPerCell[lin_ijk],lin_ijk,i,j,k);
			fprintf(stderr,"\tAllocated %ld (longs) in ListOfHalos(%ld=[%ld,%ld,%ld])\n",Nhalos,lin_ijk,i,j,k);
		}
#endif		
	}	
	}
	}

#ifdef VERB
	fprintf(stderr,"\t... memory allocated ...\n");
	t3=time(NULL);
 	diff = difftime(t3,t2);
	fprintf(stderr,"\ttime allocating %f\n",diff);
#endif

	for (ilong=0;ilong<NTotPart;ilong++) {
		i = (long) (invL * PartX[ilong]*NCells);
		j = (long) (invL * PartY[ilong]*NCells);
		k = (long) (invL * PartZ[ilong]*NCells);
		i=check_limit(i,NCells);
		j=check_limit(j,NCells);
		k=check_limit(k,NCells);
		lin_ijk = k+j*NCells+i*NCells*NCells;
		ListOfPart[lin_ijk][count[lin_ijk]] = ilong;
		count[lin_ijk]++;
	}

#ifdef DEBUG
        fprintf(stderr,"\tMass_cell[0]=%e",MassLeft[0]);
#endif

	Mhalo = HaloMass[0];
	i_alpha = 0;
	while(Mhalo<Malpha[i_alpha]) {
		i_alpha++;
		if (i_alpha==Nalpha){
			fprintf(stderr,"\tERROR: No alpha low enough found\n");
			exit(0);
		}
	}	
	Mchange = Malpha[i_alpha];
	exp = alpha[i_alpha];
	ComputeCumulative(exp);
#ifdef VERB
        fprintf(stderr,"\tNumber of alphas: %ld\n",Nalpha);
        fprintf(stderr,"\tUsing alpha_%ld=%f for M>%e\n",i_alpha,exp,Mchange);
#endif



//----------------------------------- Particles assigned to grid


#ifdef VERB
	fprintf(stderr,"\t ...done\n\n");
	t4=time(NULL);
 	diff = difftime(t4,t3);
	fprintf(stderr,"\ttime of the actual assignment %f\n",diff);
#endif
#ifdef DEBUG
	fprintf(stderr,"\t Mass Function\n");
	for (ihalo=0;ihalo<15;ihalo++){
		fprintf(stderr,"\thalo %ld: ",ihalo);
		fprintf(stderr,"M=%e\n",HaloMass[ihalo]);
	}
#endif

#ifdef VERB
	fprintf(stderr,"\n\tPlacing Halos...\n\n");
#endif

	for (ihalo=0;ihalo<NHalosTot;ihalo++){

		#ifdef DEBUG
		fprintf(stderr,"\n\t- Halo %ld ",ihalo);
		#endif
		
		do {		
		  #ifndef RANKED				
		  lin_ijk = select_cell();
		 
		  k=lin_ijk%(NCells);
		  j=((lin_ijk-k)/NCells)%NCells;
	  	  i=(lin_ijk-k-j*NCells)/(NCells*NCells);
		  #else
		  lin_ijk=select_heaviest_cell(&i,&j,&k);		  
		  #endif


		  trials=0;
		  do {
			ipart = select_part(lin_ijk);		
               		HaloX[ihalo] = PartX[ipart];
               		HaloY[ihalo] = PartY[ipart];
               		HaloZ[ihalo] = PartZ[ipart];
			R=R_from_mass(HaloMass[ihalo],rho_ref);
			HaloR[ihalo]= R;
			#ifdef NO_EXCLUSION
			check = 1;
			#else
			check = check_HaloR_in_mesh(ihalo,HaloX,HaloY,HaloZ,HaloR,i,j,k);
			#endif
			if (check==0){
				#ifdef DEBUG
				fprintf(stderr,"Refused part : %ld\n",ipart);
				#endif
				trials++;
			}
			if (trials == MAXTRIALS){
				#ifdef DEBUG
				fprintf(stderr,"MAXTRIALS=%d reached, selecting another cell\n",MAXTRIALS);
				#endif
				break;
			}
		  } while (check==0);//PART excluded
	        } while(check==0); //if reached MAXTRIALS, select another cell
		
                Mcell=MassLeft[lin_ijk];
		Mhalo= HaloMass[ihalo];

		while (Mhalo < Mchange){
			i_alpha++;		
			Mchange = Malpha[i_alpha];
			exp = alpha[i_alpha];
			ComputeCumulative(exp);
		#ifdef VERB
        		fprintf(stderr,"\n\tUsing alpha_%ld=%f for M>%e\n",i_alpha,exp,Mchange);
		#endif
		}
		
		#ifndef NO_MASS_CONSERVATION 
                if (Mcell>HaloMass[ihalo])
			MassLeft[lin_ijk] -= Mhalo; 
                else
                        MassLeft[lin_ijk] = 0.;
		#endif

		ProbDiff = pow(MassLeft[lin_ijk]/mpart,exp)-pow(Mcell/mpart,exp);

		#ifdef DEBUG
		fprintf(stderr,"\n \tassigned to cell %ld=[%ld,%ld,%ld]\n\t Before: Mcell=%e, TotProb=%e. ",lin_ijk,i,j,k,Mcell,TotProb);
		#endif

                #pragma omp parallel for private(icell) shared(CumulativeProb,ProbDiff,NTotCells,lin_ijk) default(none)
                for(icell=lin_ijk;icell<NTotCells;icell++){
                        CumulativeProb[icell]+=ProbDiff;
                }
                TotProb+=ProbDiff;

		#ifdef DEBUG
		fprintf(stderr," After: Mcell=%e, TotProb=%e.   ProbDiff=%e, Mhalo=%e\n",MassLeft[lin_ijk],TotProb,ProbDiff,Mhalo);
		#endif
	

		#ifdef DEBUG
		fprintf(stderr,"\thalo %ld assigned to particle %ld at [%f,%f,%f]. R= %f, M= %e\n",ihalo,ipart,HaloX[ihalo],HaloY[ihalo],HaloZ[ihalo],R,Mhalo);
		#endif

		ListOfHalos[lin_ijk][NHalosPerCell[lin_ijk]]=ihalo;
		NHalosPerCell[lin_ijk]++;
	}

#ifdef VERB
	t5=time(NULL);
 	diff = difftime(t5,t4);
	fprintf(stderr,"\ttime placing %f\n",diff);
 	diff = difftime(t5,t0);
	fprintf(stderr,"\ttotal time in .c %f\n",diff);
	fprintf(stderr,"\n\tPlacement done!!!\n");
#endif
	free(count); free(NPartPerCell); free(ListOfPart);
	free(CumulativeProb);
	return 0;
}



void ComputeCumulative(double alpha){
	long i;
	TotProb = 0.;
	for(i=0;i<NTotCells;i++){
		CumulativeProb[i] = TotProb + pow(MassLeft[i]/mpart,alpha);
		TotProb = CumulativeProb[i];
	}
}



long select_part(long ijk){
	long i_rnd,ipart;

	i_rnd = (long) (NPartPerCell[ijk] * ((double)rand()/(RAND_MAX+1.0)));
        ipart = ListOfPart[ijk][i_rnd];

	return ipart;
}

long select_heaviest_cell(long *x, long *y, long *z) {
	long i,j,k,lin_ijk, out_ijk=-1;
	float max=0.0;	
	for (i=0;i<NCells;i++){
        for (j=0;j<NCells;j++){
        for (k=0;k<NCells;k++){
		lin_ijk = k+j*NCells+i*NCells*NCells;
		if (max<MassLeft[lin_ijk]){
			max=MassLeft[lin_ijk];
			*x=i;
			*y=j;
			*z=k;
			out_ijk = lin_ijk;	
		}
	}
	}
	}
	return out_ijk;
}

/*
long select_cell_rnd(long *x, long *y, long *z) {
	(*x) = rand()%NCells;
	(*y) = rand()%NCells;
	(*z) = rand()%NCells;
	return (*z)+(*y)*NCells+(*x)*NCells*NCells;
}
*/

long select_cell() {
        long   i_low, i_up, i_mid;
        double d_rand;

        d_rand = TotProb * ((double)rand()/(RAND_MAX+1.0));

        i_low = 0;
        i_up  = NTotCells-1;
        while(i_low < i_up) {
                i_mid = (i_low+i_up)/2;
                if(d_rand < CumulativeProb[i_mid])
                        i_up  = i_mid;
                else
                        i_low = i_mid + 1;
         }
        return(i_low);
}





int check_HaloR_in_cell(long ipart,float *PartX, float *PartY, float *PartZ, float *PartR, long i,long j, long k){
        long jpart,jj;
        double X=PartX[ipart],Y=PartY[ipart],Z=PartZ[ipart];
	#ifndef ONLYBIG
	double R=PartR[ipart];
	#endif 


//#ifdef _PERIODIC
                if (i==-1){
                        i = NCells -1;
                        X = X + Lbox;
                }
                else if (i==NCells){
                        i = 0;
                        X = X - Lbox;
                }

                if (j==-1){
                        j = NCells -1;
                        Y = Y + Lbox;
                }
                else if (j==NCells){
                        j = 0;
                        Y = Y - Lbox;
                }
                if (k==-1){
                        k = NCells -1;
                        Z = Z + Lbox;
                }
                else if (k==NCells){
                        k = 0;
                        Z = Z - Lbox;
                }
//#endif
  if (i>=0 && i<NCells && j>=0 && j<NCells && k>=0 && k<NCells){
#ifdef ULTRADEBUG
       fprintf(stderr,"Checking cell [%ld,%ld,%ld] for halo %ld",i,j,k,ipart);
#endif
        long lin_ijk = k+j*NCells+i*NCells*NCells;
#ifdef ULTRADEBUG
	fprintf(stderr,"= %ld.",lin_ijk);
#endif
	for (jj=0; jj<NHalosPerCell[lin_ijk]; jj++){
#ifdef ULTRADEBUG
		fprintf(stderr,"jj=%ld/%ld ",jj,NHalosPerCell[lin_ijk]);
#endif
		jpart=ListOfHalos[lin_ijk][jj];
#ifdef ULTRADEBUG
		fprintf(stderr,"jpart=%ld ",jpart);
#endif
		#ifdef ONLYBIG
		if ((square(X-PartX[jpart])+square(Y-PartY[jpart])+square(Z-PartZ[jpart]))<square(PartR[jpart])) 
		#else
		if ((square(X-PartX[jpart])+square(Y-PartY[jpart])+square(Z-PartZ[jpart]))<square(R+PartR[jpart])) 
		#endif
		{
#ifdef DEBUG
                        fprintf(stderr,"\nChecking cell [%ld,%ld,%ld] for halo %ld",i,j,k,ipart);
                        fprintf(stderr," lin_ijk= %ld.",lin_ijk);
                        fprintf(stderr,"jj=%ld/%ld ",jj,NHalosPerCell[lin_ijk]);
                        fprintf(stderr,"jpart=%ld ",jpart);
                        fprintf(stderr,"refused!\n");
#endif
			return 0;
		}
	}
	return 1;
  }
  else {
	fprintf(stderr,"WARNING: Computing distances outside the box\n");
	return 1;
  }
}




int check_HaloR_in_mesh(long ihalo,float *X, float *Y, float *Z , float *R,long i,long j,long k){
	int l,m,n;
	for (l=i-1;l<=i+1;l++){
	for (m=j-1;m<=j+1;m++){
	for (n=k-1;n<=k+1;n++){
		if (check_HaloR_in_cell(ihalo,X,Y,Z,R,l,m,n)==0)
			return 0;

	}
	}
	}

	return 1;
}
