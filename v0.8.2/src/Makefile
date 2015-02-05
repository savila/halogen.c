include ../Makefile.defs

EXEC   = 2LPT-HALOGEN

OBJS   = main.o power.o allvars.o save.o read_param.o  read_glass.o  \
         nrsrc/nrutil.o nrsrc/qromb.o nrsrc/polint.o nrsrc/trapzd.o \
	 read_snapshot.o place_halos.o populate_mass_function.o

INCL   = allvars.h proto.h  nrsrc/nrutil.h  Makefile


OPTIONS =  $(OPT)


GSL_INCL =  -I/home/savila/libs/include
GSL_LIBS= -L/home/savila/libs/lib
FFTW_INCL= -I/home/savila/libs/include
FFTW_LIBS= -L/home/savila/libs/lib

CCMPI       =  /home/savila/libs/bin/mpicc  #-g -Wall -fbounds-check    # sets the C-compiler (default)
MPICHLIB =  -L/home/dmonge/mpi/lib



OPTIMIZE =   -O3 -Wall    # optimization and warning flags (default)


FFTW_LIB =  $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw

LIBS   =   -lm  $(MPICHLIB)  $(FFTW_LIB)  $(GSL_LIBS)  -lgsl -lgslcblas -fopenmp

CFLAGS =   $(OPTIONS) $(DEFS) $(OPTIMIZE)  $(FFTW_INCL) $(GSL_INCL)

$(EXEC): $(OBJS) 
	$(CCMPI) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 

halogen:	halogen.c\
		read_snapshot.o\
		read_snapshot.h\
		place_halos.o\
		place_halos.h\
		populate_mass_function.o\
		populate_mass_function.h
	$(CC) $(OPTIMIZE) $(LIBS) $(DEFS) -o halogen read_snapshot.o place_halos.o  populate_mass_function.o halogen.c

fit:            fit.c\
                halogen.c\
                read_snapshot.o\
                read_snapshot.h\
                place_halos.o\
                place_halos.h\
                populate_mass_function.o\
                populate_mass_function.h\
                correlate.o\
                correlate.h
	$(CC) $(OPT) $(LIBS) $(DEFS) $(GSL_INCL) -o fit cute.src/common.o cute.src/define.o cute.src/io.o cute.src/neighbors.o cute.src/correlator.o correlate.o read_snapshot.o place_halos.o  populate_mass_function.o fit.c 



clean:
	rm -f $(OBJS) $(EXEC) halogen fit


read_snapshot.o:        read_snapshot.c\
                read_snapshot.h
	$(CC) $(OPT) $(LIBS) $(DEFS) -c read_snapshot.c 

place_halos.o:          place_halos.c\
                        place_halos.h
	$(CC) $(OPT) $(LIBS) $(DEFS) -c place_halos.c 

populate_mass_function.o:               populate_mass_function.c\
                                        populate_mass_function.h
	$(CC) $(OPT) $(LIBS) $(DEFS) -c populate_mass_function.c 

correlate.o:		cute\
			correlate.c\
			correlate.h
	$(CC) $(OPT) $(LIBS) $(DEFS) $(GSL_INCL) -c correlate.c

cute:                   
	$(MAKE) -C cute.src cute

