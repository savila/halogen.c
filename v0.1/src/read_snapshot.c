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

#include "utility.h"
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
//#define GADGET_SKIP        {ReadUInt(icfile,&blklen,SWAPBYTES);fprintf(stderr,"GADGET_SKIP: block=%d\n",blklen);}
#define GADGET_SKIP        ReadUInt(icfile,&blklen,SWAPBYTES);
#define SIZEOFGADGETHEADER 256
#define MZERO             (1e-10)
#define X                  0
#define Y                  1
#define Z                  2


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

struct particle_data //should get rid of this struct 
{
  double     Pos[3];       /* particle position   */  
  double     Vel[3];       /* particle velocity   */  
  double     Mass;         /* particle mass       */
  double     u;            /* gas internal energy */
  long       ID;           /* particle IDs        */
} *Part;


/*=============================================================================
 *                                PROTOTYPES
 *=============================================================================*/
void read_gadget(FILE *icfile);
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
        read_gadget(icfile);
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
      fprintf(stderr,"\n\ninput: could not open file with IC's  %s\n",infile_name);
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
    
    read_gadget(icfile);
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
  fprintf(stderr,"\nquick ID check:\n");
  fprintf(stderr,"   IDmin = %ld\n",IDmin);
  fprintf(stderr,"   IDmax = %ld\n",IDmax);
  fprintf(stderr,"   IDmax-IDmin = %ld  vs.  nall[1] = %d\n\n",IDmax-IDmin,gadget.header.nall[1]);
  //exit(0);
  

	*out_mp  = GADGET_MUNIT*gadget.header.massarr[1];
	*out_Np  = gadget.nall; 
	*out_L   = gadget.header.BoxSize;
	(*out_x) = (float *) calloc(*out_Np,sizeof(float)); 
	(*out_y) = (float *) calloc(*out_Np,sizeof(float)); 
	(*out_z) = (float *) calloc(*out_Np,sizeof(float)); 

	long i;
	for (i=0; i<gadget.nall; i++){
		(*out_x)[i]=Part[i].Pos[X];
		(*out_y)[i]=Part[i].Pos[Y];
		(*out_z)[i]=Part[i].Pos[Z];
	}


	return 0;
  

  
}


/*=============================================================================
 *                                READ_GADGET
 *=============================================================================*/
void read_gadget(FILE *icfile)
{

  
  double         tot_mass[6];
  
  int            i,j,k;
  long           no_part;
  int            massflag;
  char           DATA[MAXSTRING];
  float          fdummy[3];
  double         ddummy[3];
  double         x_fac, v_fac, m_fac;
  long           pid, ldummy;
  unsigned int   idummy;
  
  /*================= read in GADGET IO header =================*/
  if(FORMAT == 2)
   {
    GADGET_SKIP;
    fread(DATA,sizeof(char),blklen,icfile);
    DATA[4] = '\0';
    fprintf(stderr,"reading... %s\n",DATA);
    //GADGET_SKIP;
    
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
  //exit(0);
  
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
  
  /* be verbose */
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
  
  fprintf(stderr,"\n-> reading %ld particles from  GADGET file #%d/%d...\n\n", no_part, gadget.i_gadget_file+1, gadget.no_gadget_files);
  
  /* allocate particle array (only once when reading the first file, of course!) */
  if(gadget.i_gadget_file == 0)
   {
    fprintf(stderr,"-> allocating %f GB of RAM for particles\n\n",(float)(gadget.nall*sizeof(struct particle_data))/1024./1024./1024.);
    if(!(Part=(struct particle_data *) calloc(gadget.nall, sizeof(struct particle_data))))
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
     
    /* get proper position in Part[] array */
    pid = get_pid(i);
    
    /* storage and conversion to comoving physical units */
    Part[pid].Pos[0] = ddummy[0] * x_fac;
    Part[pid].Pos[1] = ddummy[1] * x_fac;
    Part[pid].Pos[2] = ddummy[2] * x_fac;      
   }
  fprintf(stderr,"Pos[X]=%12.6g Pos[Y]=%12.6g Pos[Z]=%12.6g ... ",Part[no_part-1].Pos[X],Part[no_part-1].Pos[Y],Part[no_part-1].Pos[Z]);
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
  /*================= read in GADGET particles =================*/
  
  
  
  /*================= read in GADGET velocities =================*/
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
    
    /* get proper position in Part[] array */
    pid = get_pid(i);
    
    /* storage and conversion to comoving physical units */
    Part[pid].Vel[0] = ddummy[0] * v_fac;
    Part[pid].Vel[1] = ddummy[1] * v_fac;
    Part[pid].Vel[2] = ddummy[2] * v_fac; 
   }
  fprintf(stderr,"Vel[X]=%12.6g Vel[Y]=%12.6g Vel[Z]=%12.6g ... ",Part[no_part-1].Vel[X],Part[no_part-1].Vel[Y],Part[no_part-1].Vel[Z]);
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
  /*================= read in GADGET velocities =================*/
  
  
  /*================= read in GADGET id's =================*/
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
    /* get proper position in Part[] array */
    pid = get_pid(i);

    if(LGADGET)
     {
      ReadLong(icfile,&ldummy,SWAPBYTES);
      Part[pid].ID = ldummy;
     }
    else
     {
      ReadUInt(icfile,&idummy,SWAPBYTES);
      Part[pid].ID = (long) idummy;
     }
    
    /* check the ID range of the "halo" particles */
    if(gadget.header.np[0] <= i && i < gadget.header.np[0]+gadget.header.np[1])
     {
      if(Part[pid].ID > IDmax) IDmax = Part[pid].ID; 
      if(Part[pid].ID < IDmin) IDmin = Part[pid].ID; 
     }
   }
  
  fprintf(stderr,"ID=%12ld ...  ",Part[no_part-1].ID);
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
  /*================= read in GADGET id's =================*/
  
  
  k = 0;
  /* massflag == 1 indicates that massarr[i] = 0 and hence need to read in particle masses */
  if(massflag==1) 
   {
    /*================= read in GADGET individual particle masses =================*/
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
    
    for(i=0;i<6;i++)
     {
      tot_mass[i] = 0.;
      if (gadget.header.np[i] > 0 && gadget.header.massarr[i] < MZERO  ) 
       {
        
        fprintf(stderr,"  %d    ",i);
        
        for(j=0; j<gadget.header.np[i]; j++)
         {
          /* read */
          if(DGADGET)
           {
            ReadDouble(icfile,&(ddummy[0]),SWAPBYTES);
           }
          else
           {
            ReadFloat(icfile,&(fdummy[0]),SWAPBYTES);
            ddummy[0] = fdummy[0];
           }

          /* get proper position in Part[] array */
          pid = get_pid(k);
          
          /* store */
          Part[pid].Mass  = ddummy[0];
          tot_mass[i]    += ddummy[0];
          k++;
         }
       }
      else
       {
        /* simply copy appropriate massarr[i] to particles */
        for(j=0; j<gadget.header.np[i]; j++) 
         {
          /* get proper position in Part[] array */
          pid = get_pid(k);
                   
          /* store */
          Part[pid].Mass = gadget.header.massarr[i];
          k++;
         }
        tot_mass[i] = gadget.header.np[i]*gadget.header.massarr[i];
       }
     }
    
    GADGET_SKIP;
    fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
    /*================= read in GADGET individual particle masses =================*/
   }
  
  /* simply copy appropriate massarr[i] to particles */
  else 
   {
    k=0;
    for(i=0;i<6;i++)
     {
      for(j=0;j<gadget.header.np[i];j++) 
       {
        /* get proper position in Part[] array */
        pid = get_pid(k);
             
        /* store */
        Part[pid].Mass = gadget.header.massarr[i];
        k++;
       }
      tot_mass[i] = gadget.header.np[i]*gadget.header.massarr[i];
     }
   }
  
  /*============ convert masses to Msun/h and set particle type ============*/
  k=0;
  
  // 1. gas (no fiddling with Part[].u, please!)
  i=0;
  for(j=0; j<gadget.header.np[i]; j++)
   {
    pid = get_pid(k);
    Part[pid].Mass *= m_fac;
    k++;
   }
  
  // 2. all other species
  for(i=1; i<6; i++)
   {
    for(j=0; j<gadget.header.np[i]; j++)
     {
      /* get proper position in Part[] array */
      pid = get_pid(k);
      Part[pid].Mass *= m_fac;
      Part[pid].u     = -i;
      k++;
     }
   }

  /*================= read in GADGET gas particle energies =================*/
  if(gadget.header.np[0] > 0) 
   {      
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
     
     for(i=0; i<gadget.header.np[0]; i++)
      {
       /* store */
       if(DGADGET)
        {
         ReadDouble(icfile,&(ddummy[0]),SWAPBYTES);
        }
       else
        {
         ReadFloat(icfile,&(fdummy[0]),SWAPBYTES);
         ddummy[0] = fdummy[0];
        }
       
       /* get proper position in Part[] array */
       pid = get_pid(i);
              
       /* store additional gas particle property */
       Part[pid].u = ddummy[0];         
      }
     
     GADGET_SKIP;
     fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
   } 
  /*================= read in GADGET gas particle energies =================*/
  
  
  /* be verbose */
  fprintf(stderr,"\n");
  if(gadget.header.np[0] > 0) fprintf(stderr,"    gas:    tot_mass[0]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[0]*GADGET_MUNIT,tot_mass[0]/(double)gadget.header.np[0]*GADGET_MUNIT);
  if(gadget.header.np[1] > 0) fprintf(stderr,"    halo:   tot_mass[1]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[1]*GADGET_MUNIT,tot_mass[1]/(double)gadget.header.np[1]*GADGET_MUNIT);
  if(gadget.header.np[2] > 0) fprintf(stderr,"    disk:   tot_mass[2]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[2]*GADGET_MUNIT,tot_mass[2]/(double)gadget.header.np[2]*GADGET_MUNIT);
  if(gadget.header.np[3] > 0) fprintf(stderr,"    bulge:  tot_mass[3]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[3]*GADGET_MUNIT,tot_mass[3]/(double)gadget.header.np[3]*GADGET_MUNIT);
  if(gadget.header.np[4] > 0) fprintf(stderr,"    stars:  tot_mass[4]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[4]*GADGET_MUNIT,tot_mass[4]/(double)gadget.header.np[4]*GADGET_MUNIT);
  if(gadget.header.np[5] > 0) fprintf(stderr,"    bndry:  tot_mass[5]=%16.8g Msun/h (%16.8g Msun/h per particle)\n",tot_mass[5]*GADGET_MUNIT,tot_mass[5]/(double)gadget.header.np[5]*GADGET_MUNIT);
  
  fprintf(stderr,"===================================================================\n");
}

/*=============================================================================
 *                        get proper position in Part[] array
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

