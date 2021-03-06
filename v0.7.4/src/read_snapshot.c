/*==========================================================================
 * read_gadget:
 * ------------
 * initially coded by Alexander Knebe, modified by Santiago Avila
 *
 * routine to read the data of interest from one gadget snapshot
 *
 *=========================================================================*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "read_snapshot.h"
/*=============================================================================
 *                              COMMON DEFINES
 *=============================================================================*/

double       GADGET_LUNIT;
double       GADGET_MUNIT;
int          SWAPBYTES;
int          FORMAT;
int          LGADGET;
int          DGADGET;
unsigned int blklen;


#define MAXSTRING          2048 
#define GADGET_SKIP        ReadUInt(icfile,&blklen,SWAPBYTES);
#define SIZEOFGADGETHEADER 256
#define MZERO             (1e-10)

#define TRUE         1
#define FALSE        0

#define V_FACT (0.7228)

int NTHREADS;
/*=============================================================================
 *                                PROTOTYPES
 *=============================================================================*/

int ReadFloat          (FILE *fptr,float *n, int swap);
int ReadDouble         (FILE *fptr,double *n,int swap);
int ReadUInt           (FILE *fptr,unsigned int *n,int swap);
int ReadInt            (FILE *fptr,int *n,int swap);
int ReadChars          (FILE *fptr,char *s,int n);






/*=============================================================================
 *                                STRUCTURES
 *=============================================================================*/
struct info_gadget
{
  int      no_gadget_files;
  int      i_gadget_file;
  long   *(np[6]);
  long     nall;
  
  struct io_gadget_header
 {
  int      np[6];
  double   massarr[6];
  double   expansion;
  double   redshift;
  int      flagsfr;
  int      flagfeedback;
  int      nall[6];
  int      flagcooling;
  int      NumFiles;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     unused[SIZEOFGADGETHEADER- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
 } header; 
  
} gadget;


/*=============================================================================
 *                                PROTOTYPES
 *=============================================================================*/
int read_gadget(FILE *icfile,float ***out_x,float ***out_y, float ***out_z,float ***out_vx,float ***out_vy, float ***out_vz,int ifile);
long get_pid(int i);

int wrap_int(int ii,int n_grid);

float *interp_dens(int ,float *,float *, float *,double *);


double *pos_2_cic(int ,float *,float *, float *);

/*=============================================================================
 *                                   MAIN
 *=============================================================================*/
int read_snapshot(char *infile_name, int format, float lunit, float munit, int swp, int glong, int gdouble, int NCells, int nthreads, float **out_x, float **out_y, float **out_z, float **out_vx, float **out_vy, float **out_vz,long *out_Np, float *out_mp, float *out_L, float *out_omega_0,long ***ListOfPart,long **NPartPerCell,double **MassLeft,float **dens){
  char    gadget_file[MAXSTRING];
  int     no_gadget_files, i_gadget_file;
  long  ipart=0,ii,i,j,k,lin_ijk;
  FILE   *icfile;
  long   *PartPerFile;
  double invL;
  long *count;

  float **file_x, **file_y, **file_z, **file_vx,**file_vy,**file_vz; 

  FORMAT       = format;
  GADGET_LUNIT = lunit;
  GADGET_MUNIT = munit;
  SWAPBYTES    = swp;
  LGADGET      = glong;
  DGADGET      = gdouble;

  if (nthreads<1){
	NTHREADS = omp_get_max_threads();
  }else{
	NTHREADS = nthreads;
  }
  fprintf(stderr,"\tUsing OMP with %d threads\n",NTHREADS);
  
  /*==================================
   * are there multiple GADGET files?
   *==================================*/
  if((icfile = fopen(infile_name,"rb")) == NULL)
   {
    /* maybe there are multiple GADGET files ... count them! */
    no_gadget_files = 0;
    i_gadget_file   = 0;
    sprintf(gadget_file,"%s.%d",infile_name,i_gadget_file);
    while((icfile = fopen(gadget_file,"rb")) != NULL)
     {
      fclose(icfile);
      no_gadget_files++;
      i_gadget_file++;
      sprintf(gadget_file,"%s.%d",infile_name,i_gadget_file);
     }

    if(no_gadget_files > 1)
     {
      fprintf(stderr,"\n\treading GADGET data from %d files:\n",no_gadget_files);
      gadget.no_gadget_files  = no_gadget_files;
      if(!((PartPerFile)=(long *) calloc(no_gadget_files, sizeof(long))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_x)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_y)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_z)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
       #ifndef FITTING
      if(!((file_vx)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vy)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vz)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
       #endif
      
      /* allocate temporary storage for no. of particles arrays */
      gadget.np[0]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[1]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[2]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[3]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[4]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[5]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
   

    
      /* read multi-GADGET files one by one */
      #pragma omp parallel for num_threads(NTHREADS) private(i_gadget_file,icfile,gadget_file) shared(no_gadget_files,gadget,infile_name,stderr,file_x,file_y,file_z,file_vx,file_vy,file_vz,PartPerFile) default(none)
      for(i_gadget_file=0; i_gadget_file<no_gadget_files; i_gadget_file++)
       {
        sprintf(gadget_file,"%s.%d",infile_name,i_gadget_file);

//        fprintf(stderr,"\n\t===================================================================\n");
  //      fprintf(stderr,"\t=> reading %s\n\n",gadget_file);
        icfile = fopen(gadget_file,"rb");
        
        /* tell read_gadget() which file we are using at the moment */
        gadget.i_gadget_file = i_gadget_file;
        
        /* read files... */
        PartPerFile[i_gadget_file] = read_gadget(icfile,&file_x,&file_y,&file_z,&file_vx,&file_vy,&file_vz,i_gadget_file);
        fclose(icfile);
        fprintf(stderr,"\tDone with %s\n\n",gadget_file);
       } 
      
      /* free temporary storage again */
      free(gadget.np[0]);
      free(gadget.np[1]);
      free(gadget.np[2]);
      free(gadget.np[3]);
      free(gadget.np[4]);
      free(gadget.np[5]);
        fprintf(stderr,"\tEnd of parallel reading\n");
     }
    else
     {
      /* there are no multi-GADGET files */
      fprintf(stderr,"\n\ninput: could not open file  %s\n",infile_name);
      fprintf(stderr,"Remember: if using multiple files, do not include the dot at the end of the filename.\n");
      exit(0);
     }
   }
  /*===============================
   * there is only one GADGET file
   *===============================*/
  else
   {
    gadget.no_gadget_files  = 1;
    i_gadget_file    = 0;
    no_gadget_files  = 1;
      if(!((PartPerFile)=(long *) calloc(no_gadget_files, sizeof(long))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_x)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_y)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_z)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
       #ifndef FITTING
      if(!((file_vx)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vy)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
      if(!((file_vz)=(float **) calloc(no_gadget_files, sizeof(float *))))
      {
        fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
        exit(1);
      }
	#endif
    
    /* allocate temporary storage for no. of particles arrays */
    gadget.np[0]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[1]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[2]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[3]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[4]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[5]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));

    PartPerFile[i_gadget_file] = read_gadget(icfile,&file_x,&file_y,&file_z,&file_vx,&file_vy,&file_vz,0);
    fclose(icfile);
    
    /* remove temporary storage again */
    free(gadget.np[0]);
    free(gadget.np[1]);
    free(gadget.np[2]);
    free(gadget.np[3]);
    free(gadget.np[4]);
    free(gadget.np[5]);
   }
	fprintf(stderr,"\tAllocating...\n");
	//Distribute Particles
        (*NPartPerCell) = (long *) calloc(NCells*NCells*NCells,sizeof(long));
        if( *NPartPerCell == NULL) {
                fprintf(stderr,"\tplace_halos(): could not allocate %d array for NPartPerCell[]\nABORTING",NCells*NCells*NCells);
                exit(-1);
        }
	fprintf(stderr,"\t...NPartPerCell allocaled, continue...\n");
        (*ListOfPart) = (long **) calloc(NCells*NCells*NCells,sizeof(long *));
        if( *ListOfPart == NULL) {
                fprintf(stderr,"\tplace_halos(): could not allocate %d array for NPartPerCell[]\nABORTING",NCells*NCells*NCells);
                exit(-1);
        }
	fprintf(stderr,"\t...ListOFPART allocaled, continue...\n");

	count = (long *) calloc(NCells*NCells*NCells,sizeof(long));
        if( count == NULL) {
                fprintf(stderr,"\tplace_halos(): could not allocate %d array for NPartPerCell[]\nABORTING",NCells*NCells*NCells);
                exit(-1);
        }
        fprintf(stderr,"\t Counters allocated \n\n");
	ipart=0;
	invL = 1.0/gadget.header.BoxSize;
	//Counting Particles in cells
      	for(i_gadget_file=0; i_gadget_file<no_gadget_files; i_gadget_file++)
        for (ii=0;ii<PartPerFile[i_gadget_file];ii++) {

                if (file_x[i_gadget_file][ii]==gadget.header.BoxSize)
                        file_x[i_gadget_file][ii]=0.;
                if (file_y[i_gadget_file][ii]==gadget.header.BoxSize)
                        file_y[i_gadget_file][ii]=0.;
                if (file_z[i_gadget_file][ii]==gadget.header.BoxSize)
                        file_z[i_gadget_file][ii]=0.;
                i = (long) (invL * file_x[i_gadget_file][ii]*NCells);
                j = (long) (invL * file_y[i_gadget_file][ii]*NCells);
                k = (long) (invL * file_z[i_gadget_file][ii]*NCells);
                if (i<0 || i>=NCells || j<0 || j>=NCells || k<0 || k>=NCells){
                        fprintf(stderr,"\tERROR: Particle %ld at [%f,%f,%f] seems to be out of the right box interval [0.,%f)",ipart,file_x[i_gadget_file][ii],file_y[i_gadget_file][ii],file_z[i_gadget_file][ii],gadget.header.BoxSize);
                }

                lin_ijk = k+j*NCells+i*NCells*NCells;
                (*NPartPerCell)[lin_ijk]++;
		ipart++;
        }
        fprintf(stderr,"\t Particles Counted \n\n");

	//Alloccing memory for particle
        for (i=0;i<NCells;i++){
        for (j=0;j<NCells;j++){
        for (k=0;k<NCells;k++){
                lin_ijk = k+j*NCells+i*NCells*NCells;
		//fprintf(stderr,"lin_ijk=%d\n",lin_ijk);
                (*ListOfPart)[lin_ijk] = (long *) calloc((*NPartPerCell)[lin_ijk],sizeof(long));
        }
        }
        }
        fprintf(stderr,"\t Grid allocated \n\n");
	(*out_x) = (float *) calloc(ipart,sizeof(float));
	(*out_y) = (float *) calloc(ipart,sizeof(float));
	(*out_z) = (float *) calloc(ipart,sizeof(float));
	 #ifndef FITTING
	(*out_vx) = (float *) calloc(ipart,sizeof(float));
	(*out_vy) = (float *) calloc(ipart,sizeof(float));
	(*out_vz) = (float *) calloc(ipart,sizeof(float));
	 #endif
        fprintf(stderr,"\t vector allocated \n\n");
	ipart=0;
	//Distributing Particles
        //#pragma omp parallel for private(i_gadget_file,ipart,i,k,j,lin_ijk,ii) shared(PartPerFile,invL,NCells,ListOfPart,count,out_x,out_z,out_y,out_vx,out_vz,out_vy,file_x,file_y,file_z,file_vx,file_vy,file_vz,stderr,gadget,no_gadget_files) default(none)
      	for(i_gadget_file=0; i_gadget_file<no_gadget_files; i_gadget_file++)
        for (ii=0;ii<PartPerFile[i_gadget_file];ii++) {
                i = (long) (invL * file_x[i_gadget_file][ii]*NCells);
                j = (long) (invL * file_y[i_gadget_file][ii]*NCells);
                k = (long) (invL * file_z[i_gadget_file][ii]*NCells);
                if (i<0 || i>=NCells || j<0 || j>=NCells || k<0 || k>=NCells){
                        fprintf(stderr,"\tERROR: Particle %ld at [%f,%f,%f] seems to be out of the right box interval [0.,%f)",ipart,file_x[i_gadget_file][ii],file_y[i_gadget_file][ii],file_z[i_gadget_file][ii],gadget.header.BoxSize);
                }
		(*out_x)[ipart] = file_x[i_gadget_file][ii];
		(*out_y)[ipart] = file_y[i_gadget_file][ii];
		(*out_z)[ipart] = file_z[i_gadget_file][ii];
		 #ifndef FITTING
		(*out_vx)[ipart] = V_FACT * file_vx[i_gadget_file][ii];
		(*out_vy)[ipart] = V_FACT * file_vy[i_gadget_file][ii];
		(*out_vz)[ipart] = V_FACT * file_vz[i_gadget_file][ii];
		 #endif
                lin_ijk = k+j*NCells+i*NCells*NCells;
                (*ListOfPart)[lin_ijk][count[lin_ijk]] = ipart;
                count[lin_ijk]++;
		ipart++;
        }

        fprintf(stderr,"\t Particles Distributed \n\n");
 	
	(*MassLeft) =  pos_2_cic(NCells,(*out_x),(*out_y),(*out_z));
	(*dens) = (float *) interp_dens(NCells,(*out_x),(*out_y),(*out_z),(*MassLeft));
 
	*out_mp  = GADGET_MUNIT*gadget.header.massarr[1];
	*out_Np  = gadget.nall; 
	*out_L   = gadget.header.BoxSize;


	return 0;
  

  
}


/*=============================================================================
 *                                READ_GADGET
 *=============================================================================*/
int read_gadget(FILE *icfile, float ***out_x,float ***out_y,float ***out_z, float ***out_vx,float ***out_vy,float ***out_vz,int ifile)
{

  int            i;
  long           no_part;
  int            massflag;
  char           DATA[MAXSTRING];
  float          fdummy[3];
  double         ddummy[3];
  double         x_fac;
//  long           pid;

  
  /*================= read in GADGET IO header =================*/
  if(FORMAT == 2)
   {
    GADGET_SKIP;
    fread(DATA,sizeof(char),blklen,icfile);
    DATA[4] = '\0';
  if (ifile==0)
    fprintf(stderr,"\treading... %s\n",DATA);
    
    GADGET_SKIP;
   }
  
  GADGET_SKIP;
  
  ReadInt(icfile,&(gadget.header.np[0]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.np[1]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.np[2]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.np[3]),SWAPBYTES);    /* number of particles in current file */
  ReadInt(icfile,&(gadget.header.np[4]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.np[5]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[0]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[1]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[2]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[3]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[4]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.massarr[5]),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.expansion),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.redshift),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.flagsfr),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.flagfeedback),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[0]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[1]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[2]),SWAPBYTES);  /* total number of particles in simulation */
  ReadInt(icfile,&(gadget.header.nall[3]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[4]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.nall[5]),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.flagcooling),SWAPBYTES);
  ReadInt(icfile,&(gadget.header.NumFiles),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.BoxSize),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.Omega0),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.OmegaLambda),SWAPBYTES);
  ReadDouble(icfile,&(gadget.header.HubbleParam),SWAPBYTES);
  ReadChars(icfile,&(gadget.header.unused[0]),SIZEOFGADGETHEADER- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8);
  
  GADGET_SKIP;
  /*================= read in GADGET IO header =================*/

  
  /* keep track of no. of particles in each GADGET file */
  gadget.np[0][gadget.i_gadget_file] = gadget.header.np[0];
  gadget.np[1][gadget.i_gadget_file] = gadget.header.np[1];
  gadget.np[2][gadget.i_gadget_file] = gadget.header.np[2];
  gadget.np[3][gadget.i_gadget_file] = gadget.header.np[3];
  gadget.np[4][gadget.i_gadget_file] = gadget.header.np[4];
  gadget.np[5][gadget.i_gadget_file] = gadget.header.np[5];
  
  /* conversion factors to Mpc/h, km/sec, Msun/h */
  x_fac  = GADGET_LUNIT;
  //m_fac  = GADGET_MUNIT;
  
  /* count total no. of particles in current file (and set massflag) */
  massflag    = 0;
  no_part     = 0;
  gadget.nall = 0;
  for(i=0;i<6;i++) 
   {
    no_part     += gadget.header.np[i];
    gadget.nall += gadget.header.nall[i];
    if(gadget.header.massarr[i] < MZERO && gadget.header.np[i] > 0)
      massflag=1;  
   }  
  
  #ifdef VERB
  if (ifile==0){
  fprintf(stderr,"\n\texpansion factor: %lf\n",             gadget.header.expansion);
  fprintf(stderr,"\tredshift:         %lf\n",             gadget.header.redshift);
  fprintf(stderr,"\tboxsize:          %lf (%lf Mpc/h)\n", gadget.header.BoxSize,gadget.header.BoxSize*GADGET_LUNIT);
  fprintf(stderr,"\tomega0:           %lf\n",             gadget.header.Omega0);
  fprintf(stderr,"\tlambda0:          %lf\n",             gadget.header.OmegaLambda);
  fprintf(stderr,"\tHubbleParam:      %lf\n\n",           gadget.header.HubbleParam);
  }
  #endif


  #ifdef DEBUG
  fprintf(stderr,"\n\tgas:    np[0]=%9d\t nall[0]=%9d\t massarr[0]=%g\n",gadget.header.np[0],gadget.header.nall[0],gadget.header.massarr[0]); 
  fprintf(stderr,"\thalo:   np[1]=%9d\t nall[1]=%9d\t massarr[1]=%g\n",gadget.header.np[1],gadget.header.nall[1],gadget.header.massarr[1]); 
  fprintf(stderr,"\tdisk:   np[2]=%9d\t nall[2]=%9d\t massarr[2]=%g\n",gadget.header.np[2],gadget.header.nall[2],gadget.header.massarr[2]); 
  fprintf(stderr,"\tbulge:  np[3]=%9d\t nall[3]=%9d\t massarr[3]=%g\n",gadget.header.np[3],gadget.header.nall[3],gadget.header.massarr[3]); 
  fprintf(stderr,"\tstars:  np[4]=%9d\t nall[4]=%9d\t massarr[4]=%g\n",gadget.header.np[4],gadget.header.nall[4],gadget.header.massarr[4]); 
  fprintf(stderr,"\tbndry:  np[5]=%9d\t nall[5]=%9d\t massarr[5]=%g\n",gadget.header.np[5],gadget.header.nall[5],gadget.header.massarr[5]); 
  fprintf(stderr,"\n\t-> reading %ld/%ld particles from  GADGET file #%d/%d...\n\n", no_part,gadget.nall, gadget.i_gadget_file+1, gadget.no_gadget_files);
  #endif


    #ifdef DEBUG
    fprintf(stderr,"\t-> allocating %f GB of RAM for particles\n\n",(float)(no_part*6*sizeof(float))/1024./1024./1024.);
    #endif



    if(!((*out_x)[ifile]=(float *) calloc(no_part, sizeof(float))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }

    if(!((*out_y)[ifile]=(float *) calloc(no_part, sizeof(float))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
    if(!((*out_z)[ifile]=(float *) calloc(no_part, sizeof(float))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
    #ifndef FITTING 
    if(!((*out_vx)[ifile]=(float *) calloc(no_part, sizeof(float))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
    if(!((*out_vy)[ifile]=(float *) calloc(no_part, sizeof(float))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
    if(!((*out_vz)[ifile]=(float *) calloc(no_part, sizeof(float))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
     #endif

   
  /*================= read in GADGET particles =================*/
  if(FORMAT == 2)
   {
    GADGET_SKIP;
    fread(DATA,sizeof(char),blklen,icfile);
    DATA[4] = '\0';
  if (ifile==0)
    fprintf(stderr,"\treading... %s",DATA);
    //GADGET_SKIP;
    
    GADGET_SKIP;
   }
  else
   {
  if (ifile==0)
    fprintf(stderr,"\treading ");
   }
  
  GADGET_SKIP;
  #ifdef DEBUG
  fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
  #endif

  fprintf(stderr,"Thread %d/%d\n",omp_get_thread_num(),omp_get_num_threads());
  //#pragma omp barrier
  


  for(i=0;i<no_part;i++)
   {    
    /* read */
     if(DGADGET)
      {
       ReadDouble(icfile,&(ddummy[0]),SWAPBYTES);
       ReadDouble(icfile,&(ddummy[1]),SWAPBYTES);
       ReadDouble(icfile,&(ddummy[2]),SWAPBYTES);
      }
     else
      {
       ReadFloat(icfile,&(fdummy[0]),SWAPBYTES);
       ReadFloat(icfile,&(fdummy[1]),SWAPBYTES);
       ReadFloat(icfile,&(fdummy[2]),SWAPBYTES);
      }
    
   
    //if (i==0)
//	fprintf(stderr,"Thread %d/%d: nstart = %ld\n",omp_get_thread_num(),omp_get_num_threads(),pid);

 
    /* storage and conversion to comoving physical units */
    (*out_x)[ifile][i] = fdummy[0] * x_fac;
    (*out_y)[ifile][i] = fdummy[1] * x_fac;
    (*out_z)[ifile][i] = fdummy[2] * x_fac;
   }

     fprintf(stderr,"\tpositions done: ifile %d\n",ifile);

  #ifdef DEBUG
//  fprintf(stderr,"Pos[X]=%12.6g Pos[Y]=%12.6g Pos[Z]=%12.6g ... ",(*out_x)[no_part-1],(*out_y)[no_part-1],(*out_z)[no_part-1]);
  #endif 
  GADGET_SKIP;
  #ifdef DEBUG
  fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
  #endif
  /*================= read in GADGET particles =================*/
  /*================= read in GADGET velocities =================*/  
  #ifndef FITTING
  if(FORMAT == 2)
   {
    GADGET_SKIP;
    fread(DATA,sizeof(char),blklen,icfile);
    DATA[4] = '\0';
    fprintf(stderr,"reading... %s",DATA);
    GADGET_SKIP;
   }
  else
   {
  if (ifile==0)
    fprintf(stderr,"reading ");
   }

  GADGET_SKIP;
  if (ifile==0)
  fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);

  for(i=0;i<no_part;i++)
   {
    /* read */
    if(DGADGET)
     {
      ReadDouble(icfile,&(ddummy[0]),SWAPBYTES);
      ReadDouble(icfile,&(ddummy[1]),SWAPBYTES);
      ReadDouble(icfile,&(ddummy[2]),SWAPBYTES);
     }
    else
     {
      ReadFloat(icfile,&(fdummy[0]),SWAPBYTES);
      ReadFloat(icfile,&(fdummy[1]),SWAPBYTES);
      ReadFloat(icfile,&(fdummy[2]),SWAPBYTES);
      ddummy[0] = fdummy[0];
      ddummy[1] = fdummy[1];
      ddummy[2] = fdummy[2];
     }

     (*out_vx)[ifile][i] = fdummy[0];
     (*out_vy)[ifile][i] = fdummy[1];
     (*out_vz)[ifile][i] = fdummy[2];
   }
     fprintf(stderr,"\tvelocities  done: ifile %d\n",ifile);
   #endif
  /*================= read in GADGET velocities =================*/


 
  /* massflag == 1 indicates that massarr[i] = 0 which shouldnt be the case for HALOGEN */
  if(massflag==1) 
   {
	fprintf(stderr,"ERROR: HALOGEN does not expect to encounter varying masses\n");
	exit(0);
   }
   return no_part;
  
}

/*=============================================================================
 *                        get proper position in *out_x[] array
 *=============================================================================*/
long get_pid(int i)
{
  long pid;
  long itype, ifile;
  
  pid = 0;
  for(ifile=0; ifile<gadget.i_gadget_file; ifile++)
    for(itype=0; itype<6; itype++)
      pid += gadget.np[itype][ifile];
  pid += i;
  
  return(pid);
}

/*
 Read a possibly byte swapped integer
 */
int ReadInt(FILE *fptr,int *n,int swap)
{
  unsigned char *cptr,tmp;

  if(sizeof(int) != 4)
   {
    fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
    exit(0);
   }

  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped unsigned integer
 */
int ReadUInt(FILE *fptr,unsigned int *n,int swap)
{
  unsigned char *cptr,tmp;

  if(sizeof(int) != 4)
   {
    fprintf(stderr,"ReadUInt: sizeof(int)=%ld and not 4\n",sizeof(int));
    exit(0);
   }

  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped double precision number
 Assume IEEE
 */
int ReadDouble(FILE *fptr,double *n,int swap)
{
  unsigned char *cptr,tmp;

  if(sizeof(double) != 8)
   {
    fprintf(stderr,"ReadDouble: sizeof(double)=%ld and not 8\n",sizeof(double));
    exit(0);
   }

  if (fread(n,8,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[7];
    cptr[7] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[6];
    cptr[6] = tmp;
    tmp     = cptr[2];
    cptr[2] = cptr[5];
    cptr[5] = tmp;
    tmp     = cptr[3];
    cptr[3] = cptr[4];
    cptr[4] = tmp;
  }

  return(TRUE);
}

/*
 Read a possibly byte swapped floating point number
 Assume IEEE format
 */
int ReadFloat(FILE *fptr,float *n, int swap)
{
  unsigned char *cptr,tmp;

  if(sizeof(float) != 4)
   {
    fprintf(stderr,"ReadFloat: sizeof(float)=%ld and not 4\n",sizeof(float));
    exit(0);
   }

  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap)
   {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
   }
  return(TRUE);
}


/*
 Read an array of n characters
 NOTE: the difference to ReadString() is that we do not '\0'-terminate the array
 */
int ReadChars(FILE *fptr,char *s,int n)
{
  int i,c;

  if(sizeof(char) != 1)
   {
    fprintf(stderr,"ReadChars: sizeof(char)=%ld and not 1\n",sizeof(char));
    exit(0);
   }

  s[0] = '\0';
  for (i=0;i<n;i++) {
    c = fgetc(fptr);
    if (c == EOF)
      return(FALSE);
    s[i] = c;
  }
  return(TRUE);
}

//Module by D.Alonso
double *pos_2_cic(int n_grid,float *X,float *Y, float *Z)
{
  //////
  // Returns a grid[n_grid*n_grid*n_grid] with    
  // the density field calculated using the CIC
  // (cloud-in-cell) technique from a set of np
  // particle positions **pos.
  fprintf(stderr,"Calculating CIC...\n");
  int np = gadget.nall;
  double l_box=gadget.header.BoxSize;
  double mpart = GADGET_MUNIT*gadget.header.massarr[1];
  double *grid;

  grid =(double *)calloc(n_grid*n_grid*n_grid,sizeof(double));
  if(grid==NULL) {
	fprintf(stderr,"Couldnt allocate grid\n");
  }

  double t1,t2;
  t1=omp_get_wtime();
  fprintf(stderr,"Using %d threads\n",omp_get_max_threads());

#pragma omp parallel num_threads(NTHREADS) default(none)      \
  shared(X,Y,Z,n_grid,np,grid,l_box,mpart)
  {
    int ii;
    double agrid=l_box/n_grid;
    double ivcell=1./(agrid*agrid*agrid);
#pragma omp for  nowait
    for(ii=0;ii<np;ii++) {
      int i0x,i0y,i0z;
      int i1x,i1y,i1z;
      double a1x,a1y,a1z;
      double a0x,a0y,a0z;
      double xp=X[ii];
      double yp=Y[ii];
      double zp=Z[ii];
      i0x=(int)(xp/agrid);
      i0y=(int)(yp/agrid);
      i0z=(int)(zp/agrid);

      a1x=xp-(i0x+0.5)*agrid;
      a1y=yp-(i0y+0.5)*agrid;
      a1z=zp-(i0z+0.5)*agrid;

      if(a1x<0) {
        a1x=-a1x;
        i1x=wrap_int(i0x-1,n_grid);
      }
      else
        i1x=wrap_int(i0x+1,n_grid);
      a0x=agrid-a1x;

      if(a1y<0) {
        a1y=-a1y;
        i1y=wrap_int(i0y-1,n_grid);
      }
      else
        i1y=wrap_int(i0y+1,n_grid);
      a0y=agrid-a1y;

      if(a1z<0) {
        a1z=-a1z;
        i1z=wrap_int(i0z-1,n_grid);
      }
      else
        i1z=wrap_int(i0z+1,n_grid);
      a0z=agrid-a1z;

#pragma omp atomic
      grid[i0z+n_grid*(i0y+n_grid*i0x)]+=
        a0x*a0y*a0z*ivcell*mpart;
#pragma omp atomic
      grid[i1z+n_grid*(i0y+n_grid*i0x)]+=
        a1x*a0y*a0z*ivcell*mpart;
#pragma omp atomic
      grid[i0z+n_grid*(i1y+n_grid*i0x)]+=
        a0x*a1y*a0z*ivcell*mpart;
#pragma omp atomic
      grid[i1z+n_grid*(i1y+n_grid*i0x)]+=
        a1x*a1y*a0z*ivcell*mpart;
#pragma omp atomic
      grid[i0z+n_grid*(i0y+n_grid*i1x)]+=
        a0x*a0y*a1z*ivcell*mpart;
#pragma omp atomic
      grid[i1z+n_grid*(i0y+n_grid*i1x)]+=
        a1x*a0y*a1z*ivcell*mpart;
#pragma omp atomic
      grid[i0z+n_grid*(i1y+n_grid*i1x)]+=
        a0x*a1y*a1z*ivcell*mpart;
#pragma omp atomic
      grid[i1z+n_grid*(i1y+n_grid*i1x)]+=
        a1x*a1y*a1z*ivcell*mpart;
    } //end omp for
  } //end omp parallel
  t2=omp_get_wtime();
  fprintf(stderr,"Took %.3lf seconds\n",t2-t1);


  return grid;
}


float *interp_dens(int n_grid,float *X,float *Y, float *Z,double *grid)
{


  fprintf(stderr,"Interpolating back to particles...\n");
  int np = gadget.nall;
  double l_box=gadget.header.BoxSize;
  float *dens;
  dens = (float *)calloc(np,sizeof(float));
  if(dens==NULL) {
	fprintf(stderr,"Couldnt allocate *dens\n");
  }

  double t1,t2;
  t1=omp_get_wtime();
  fprintf(stderr,"Using %d threads\n",omp_get_max_threads());

#pragma omp parallel num_threads(NTHREADS) default(none)      \
  shared(X,Y,Z,n_grid,np,grid,l_box,dens)
  {
    int ii;
    double agrid=l_box/n_grid,d_temp;
    double ivcell=1./(agrid*agrid*agrid);
#pragma omp for nowait
    for(ii=0;ii<np;ii++) {
      int i0x,i0y,i0z;
      int i1x,i1y,i1z;
      double a1x,a1y,a1z;
      double a0x,a0y,a0z;
      double xp=X[ii];
      double yp=Y[ii];
      double zp=Z[ii];
      i0x=(int)(xp/agrid);
      i0y=(int)(yp/agrid);
      i0z=(int)(zp/agrid);

      a1x=xp-(i0x+0.5)*agrid;
      a1y=yp-(i0y+0.5)*agrid;
      a1z=zp-(i0z+0.5)*agrid;

      if(a1x<0) {
        a1x=-a1x;
        i1x=wrap_int(i0x-1,n_grid);
      }
      else
        i1x=wrap_int(i0x+1,n_grid);
      a0x=agrid-a1x;

      if(a1y<0) {
        a1y=-a1y;
        i1y=wrap_int(i0y-1,n_grid);
      }
      else
        i1y=wrap_int(i0y+1,n_grid);
      a0y=agrid-a1y;

      if(a1z<0) {
        a1z=-a1z;
        i1z=wrap_int(i0z-1,n_grid);
      }
      else
        i1z=wrap_int(i0z+1,n_grid);
      a0z=agrid-a1z;

//#pragma omp atomic read
      d_temp = grid[i0z+n_grid*(i0y+n_grid*i0x)]*ivcell*ivcell*a0x*a0y*a0z;
//#pragma omp atomic
      d_temp += grid[i1z+n_grid*(i0y+n_grid*i0x)]*ivcell*ivcell*a1x*a0y*a0z;
//#pragma omp atomic
      d_temp += grid[i0z+n_grid*(i1y+n_grid*i0x)]*ivcell*ivcell*a0x*a1y*a0z;
//#pragma omp atomic
      d_temp += grid[i1z+n_grid*(i1y+n_grid*i0x)]*ivcell*ivcell*a1x*a1y*a0z;
//#pragma omp atomic
      d_temp += grid[i0z+n_grid*(i0y+n_grid*i1x)]*ivcell*ivcell*a0x*a0y*a1z;
//#pragma omp atomic
      d_temp += grid[i1z+n_grid*(i0y+n_grid*i1x)]*ivcell*ivcell*a1x*a0y*a1z;
//#pragma omp atomic
      d_temp += grid[i0z+n_grid*(i1y+n_grid*i1x)]*ivcell*ivcell*a0x*a1y*a1z;
//#pragma omp atomic
      d_temp += grid[i1z+n_grid*(i1y+n_grid*i1x)]*ivcell*ivcell*a1x*a1y*a1z;
	
      dens[ii]= (float) d_temp;
    } //end omp for
  } //end omp parallel

  t2=omp_get_wtime();
  fprintf(stderr,"Took %.3lf seconds\n",t2-t1);


  return dens;
}


int wrap_int(int ii,int n_grid)
{
  if(ii<0)
    return wrap_int(ii+n_grid,n_grid);
  else
    return ii%n_grid;
}
