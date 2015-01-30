include ../Makefile.defs

EXEC   = 2LPT_HALOGEN

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o \
	 read_snapshot.o place_halos.o populate_mass_function.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile


#OPT   +=  -DPRODUCEGAS   # Set this to automatically produce gas particles 
                         # for a single DM species in the input file by interleaved by a half a grid spacing


#OPT   +=  -DMULTICOMPONENTGLASSFILE  # set this if the initial glass file contains multiple components

#OPT   +=  -DDIFFERENT_TRANSFER_FUNC  # set this if you want to implement a transfer function that depends on
                                     # particle type

OPT   +=  -DNO64BITID    # switch this on if you want normal 32-bit IDs
#OPT   +=  -DCORRECT_CIC  # only switch this on if particles start from a glass (as opposed to grid)

#OPT += -DONLY_ZA # swith this on if you want ZA initial conditions (2LPT otherwise)



OPTIONS =  $(OPT)


GSL_INCL =  -I/home/savila/libs/include/gsl
GSL_LIBS= -L/home/savila/libs/lib
FFTW_INCL= -I/home/savila/libs/include
FFTW_LIBS= -L/home/savila/libs/lib

CC       =  /home/savila/libs/bin/mpicc  #-g -Wall -fbounds-check    # sets the C-compiler (default)
MPICHLIB =  -L/home/dmonge/mpi/lib



OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)


FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas -fopenmp

CFLAGS =   $(OPTIONS) $(DEFS) $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 

halogen:	halogen.c\
		read_snapshot.o\
		read_snapshot.h\
		place_halos.o\
		place_halos.h\
		populate_mass_function.o\
		populate_mass_function.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -o halogen read_snapshot.o place_halos.o  populate_mass_function.o halogen.c $(MLIB)


.PHONY : clean
clean:
	rm -f $(OBJS) $(EXEC)


read_snapshot.o:        read_snapshot.c\
                read_snapshot.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -c read_snapshot.c 

place_halos.o:          place_halos.c\
                        place_halos.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -c place_halos.c 

populate_mass_function.o:               populate_mass_function.c\
                                        populate_mass_function.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -c populate_mass_function.c 

