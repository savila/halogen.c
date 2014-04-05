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

long         IDmin =  4294967296;
long         IDmax = -4294967296; //simply for debugging purposes


#define MAXSTRING          2048 
#define GADGET_SKIP        ReadUInt(icfile,&blklen,SWAPBYTES);
#define SIZEOFGADGETHEADER 256
#define MZERO             (1e-10)
#define X                  0
#define Y                  1
#define Z                  2

#define TRUE         1
#define FALSE        0




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
void read_gadget(FILE *icfile,float **out_x,float **out_y, float **out_z);
long get_pid(int i);



/*=============================================================================
 *                                   MAIN
 *=============================================================================*/
int read_snapshot(char *infile_name, int format, float **out_x, float **out_y, float **out_z, long *out_Np, float *out_mp, float *out_L, float *out_omega_0){
  char    gadget_file[MAXSTRING];
  int     no_gadget_files, i_gadget_file;
  FILE   *icfile;

  
  FORMAT       = format;
  
  #ifdef _KPC
  	GADGET_LUNIT = 1000.;
  #else
	GADGET_LUNIT=1.;
  #endif

  GADGET_MUNIT = 1.0e10;

  #ifdef _SWAPBYTES
  	SWAPBYTES    = 1;
  #else
	SWAPBYTES    =0;
  #endif


  #ifdef _LONGGADGET
  	LGADGET      = 1;
  #else
	LGADGET =0;
  #endif

  #ifdef _DOUBLEGADGET
  	DGADGET      = 1;
  #else
	DGADGET      = 0;
  #endif
  
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
      fprintf(stderr,"\nreading GADGET data from %d files:\n",no_gadget_files);
      gadget.no_gadget_files  = no_gadget_files;
      
      /* allocate temporary storage for no. of particles arrays */
      gadget.np[0]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[1]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[2]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[3]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[4]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[5]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
            
      /* read multi-GADGET files one by one */
      for(i_gadget_file=0; i_gadget_file<no_gadget_files; i_gadget_file++)
       {
        sprintf(gadget_file,"%s.%d",infile_name,i_gadget_file);
        fprintf(stderr,"\n===================================================================\n");
        fprintf(stderr,"=> reading %s\n\n",gadget_file);
        icfile = fopen(gadget_file,"rb");
        
        /* tell read_gadget() which file we are using at the moment */
        gadget.i_gadget_file = i_gadget_file;
        
        /* read files... */
        read_gadget(icfile,out_x,out_y,out_z);
        fclose(icfile);
       } 
      
      /* free temporary storage again */
      free(gadget.np[0]);
      free(gadget.np[1]);
      free(gadget.np[2]);
      free(gadget.np[3]);
      free(gadget.np[4]);
      free(gadget.np[5]);
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
    gadget.i_gadget_file    = 0;
    
    /* allocate temporary storage for no. of particles arrays */
    gadget.np[0]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[1]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[2]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[3]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[4]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    gadget.np[5]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
    
    read_gadget(icfile,out_x,out_y,out_z);
    fclose(icfile);
    
    /* remove temporary storage again */
    free(gadget.np[0]);
    free(gadget.np[1]);
    free(gadget.np[2]);
    free(gadget.np[3]);
    free(gadget.np[4]);
    free(gadget.np[5]);
   }
  
  
  /* check the range of the IDs of "halo" particles */
  #ifdef _DEBUG
  fprintf(stderr,"\nquick ID check:\n");
  fprintf(stderr,"   IDmin = %ld\n",IDmin);
  fprintf(stderr,"   IDmax = %ld\n",IDmax);
  fprintf(stderr,"   IDmax-IDmin = %ld  vs.  nall[1] = %d\n\n",IDmax-IDmin,gadget.header.nall[1]);
  #endif 
  

	*out_mp  = GADGET_MUNIT*gadget.header.massarr[1];
	*out_Np  = gadget.nall; 
	*out_L   = gadget.header.BoxSize;


	return 0;
  

  
}


/*=============================================================================
 *                                READ_GADGET
 *=============================================================================*/
void read_gadget(FILE *icfile, float **out_x,float **out_y,float **out_z)
{

  int            i;
  long           no_part;
  int            massflag;
  char           DATA[MAXSTRING];
  float          fdummy[3];
  double         ddummy[3];
  double         x_fac, v_fac, m_fac;
  long           pid;

  
  /*================= read in GADGET IO header =================*/
  if(FORMAT == 2)
   {
    GADGET_SKIP;
    fread(DATA,sizeof(char),blklen,icfile);
    DATA[4] = '\0';
    fprintf(stderr,"reading... %s\n",DATA);
    
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
  v_fac  = sqrt(gadget.header.expansion);
  m_fac  = GADGET_MUNIT;
  
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
  
  #ifdef _VERB
  fprintf(stderr,"expansion factor: %lf\n",             gadget.header.expansion);
  fprintf(stderr,"redshift:         %lf\n",             gadget.header.redshift);
  fprintf(stderr,"boxsize:          %lf (%lf Mpc/h)\n", gadget.header.BoxSize,gadget.header.BoxSize*GADGET_LUNIT);
  fprintf(stderr,"omega0:           %lf\n",             gadget.header.Omega0);
  fprintf(stderr,"lambda0:          %lf\n",             gadget.header.OmegaLambda);
  fprintf(stderr,"HubbleParam:      %lf\n\n",           gadget.header.HubbleParam);
  
  fprintf(stderr,"gas:    np[0]=%9d\t nall[0]=%9d\t massarr[0]=%g\n",gadget.header.np[0],gadget.header.nall[0],gadget.header.massarr[0]); 
  fprintf(stderr,"halo:   np[1]=%9d\t nall[1]=%9d\t massarr[1]=%g\n",gadget.header.np[1],gadget.header.nall[1],gadget.header.massarr[1]); 
  fprintf(stderr,"disk:   np[2]=%9d\t nall[2]=%9d\t massarr[2]=%g\n",gadget.header.np[2],gadget.header.nall[2],gadget.header.massarr[2]); 
  fprintf(stderr,"bulge:  np[3]=%9d\t nall[3]=%9d\t massarr[3]=%g\n",gadget.header.np[3],gadget.header.nall[3],gadget.header.massarr[3]); 
  fprintf(stderr,"stars:  np[4]=%9d\t nall[4]=%9d\t massarr[4]=%g\n",gadget.header.np[4],gadget.header.nall[4],gadget.header.massarr[4]); 
  fprintf(stderr,"bndry:  np[5]=%9d\t nall[5]=%9d\t massarr[5]=%g\n",gadget.header.np[5],gadget.header.nall[5],gadget.header.massarr[5]); 
  
  fprintf(stderr,"\n-> reading %ld/%ld particles from  GADGET file #%d/%d...\n\n", no_part,gadget.nall, gadget.i_gadget_file+1, gadget.no_gadget_files);
  #endif
  /* allocate particle array (only once when reading the first file, of course!) */
  if(gadget.i_gadget_file == 0)
   {
    fprintf(stderr,"-> allocating %f GB of RAM for particles\n\n",(float)(gadget.nall*3*sizeof(float))/1024./1024./1024.);
    if(!((*out_x)=(float *) calloc(gadget.nall, sizeof(float))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
    if(!((*out_y)=(float *) calloc(gadget.nall, sizeof(float))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }
    if(!((*out_z)=(float *) calloc(gadget.nall, sizeof(float))))
     {
      fprintf(stderr,"\nfailed to allocate memory for GADGET data\n");
      exit(1);
     }

   }
   
  /*================= read in GADGET particles =================*/
  if(FORMAT == 2)
   {
    GADGET_SKIP;
    fread(DATA,sizeof(char),blklen,icfile);
    DATA[4] = '\0';
    fprintf(stderr,"reading... %s",DATA);
    //GADGET_SKIP;
    
    GADGET_SKIP;
   }
  else
   {
    fprintf(stderr,"reading ");
   }
  
  GADGET_SKIP;
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
     
    /* get proper position in  array */
    pid = get_pid(i);
    
    /* storage and conversion to comoving physical units */
    (*out_x)[pid] = ddummy[0] * x_fac;
    (*out_y)[pid] = ddummy[1] * x_fac;
    (*out_z)[pid] = ddummy[2] * x_fac;
   }
  fprintf(stderr,"Pos[X]=%12.6g Pos[Y]=%12.6g Pos[Z]=%12.6g ... ",*out_x[no_part-1],*out_y[no_part-1],*out_z[no_part-1]);
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
  /*================= read in GADGET particles =================*/
  

 
  /* massflag == 1 indicates that massarr[i] = 0 which shouldnt be the case for HALOGEN */
  if(massflag==1) 
   {
	fprintf(stderr,"ERROR: HALOGEN does not expect to encounter variable masses\n If needed, please, contact us to implement it\n");
	exit(0);
   }

  
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

