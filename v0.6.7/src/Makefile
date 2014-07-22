CC		=	gcc
#OPT		=	-O2 -fopenmp 
OPT		=	-O2 -Wall -fopenmp -g
MLIB		=	-lm

include ../Makefile.defs	

halogen:	halogen.c\
		read_snapshot.o\
		read_snapshot.h\
		place_halos.o\
		place_halos.h\
		populate_mass_function.o\
		populate_mass_function.h
	$(CC) $(OPT) $(GSLIB) $(DEFS) -o halogen read_snapshot.o place_halos.o  populate_mass_function.o halogen.c $(MLIB)
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
