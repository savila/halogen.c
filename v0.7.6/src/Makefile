CC		=	gcc
#OPT		=	-O2 -fopenmp 
OPT		=	-O2 -Wall -fopenmp -g
MLIB		=	-lm
#GSLIB           =       -I/usr/local/Cellar/gsl/1.15/include -L/usr/local/Cellar/gsl/1.15/lib -lgslcblas -lgsl
GSLIB    = -I/home/savila/libs/include -L/home/savila/libs/lib -lgslcblas -lgsl
include ../Makefile.defs
	
test:		test_reading.c\
                read_snapshot.o\
                read_snapshot.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -o test read_snapshot.o test_reading.c $(MLIB)


halogen:	halogen.c\
		read_snapshot.o\
		read_snapshot.h\
		place_halos.o\
		place_halos.h\
		populate_mass_function.o\
		populate_mass_function.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -o halogen read_snapshot.o place_halos.o  populate_mass_function.o halogen.c $(MLIB)

fit:		fit.c\
		halogen.c\
                read_snapshot.o\
                read_snapshot.h\
                place_halos.o\
                place_halos.h\
                populate_mass_function.o\
                populate_mass_function.h\
		correlate.o\
		correlate.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -o fit cute.src/common.o cute.src/define.o cute.src/io.o cute.src/neighbors.o cute.src/correlator.o correlate.o read_snapshot.o place_halos.o  populate_mass_function.o fit.c $(MLIB)

correlate.o:		cute\
			correlate.c\
			correlate.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -c correlate.c

cute:			
			$(MAKE) -C cute.src cute

read_snapshot.o:	read_snapshot.c\
		read_snapshot.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -c read_snapshot.c 

place_halos.o:		place_halos.c\
			place_halos.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -c place_halos.c 

populate_mass_function.o:		populate_mass_function.c\
					populate_mass_function.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -c populate_mass_function.c 

clean:
	rm -f *.o *~ halogen
	rm -f cute.src/*.o
