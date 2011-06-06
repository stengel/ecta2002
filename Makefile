# the file "Depend" must exist for make to run

.SUFFIXES:
.SUFFIXES: .c .o .h

CC = gcc
# -lm  links the math library in
CFLAGS = -ansi -Wall -O3 

# needed for gmp library
LIBDIR     = /usr/local/lib

%.o : %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

DEPFILE = Depend

COMOBJ  = alloc.o col.o mp.o rat.o

INLEMOBJ  = lemke.o inlemke.o 

TREEOBJ = leaves.o prior.o rataux.o sfnf.o treegen.o treedef.o

GMPOBJ = gmpwrap.o glemke.o

METHOBJ = normform.o seqform.o rsf.o main.o gambit.o interface.o

ALLOBJ = $(COMOBJ) $(INLEMOBJ) $(TREEOBJ) $(METHOBJ) $(GMPOBJ)

# local variables for single substitutions

BINTREE = $(COMOBJ) $(TREEOBJ) $(METHOBJ) lemke.o

GMP = $(COMOBJ) $(TREEOBJ) $(METHOBJ) $(GMPOBJ) 

gbin: $(GMP)
	$(CC) $(CFLAGS) -L$(LIBDIR) -lgmp $(GMP)  -o gbin

bintr: $(BINTREE) 	
	$(CC) $(CFLAGS) -lm $(BINTREE)  -o bintr

INLEMKE = $(COMOBJ) $(INLEMOBJ)
inlemke: $(INLEMKE)
	$(CC) $(CFLAGS) $(INLEMKE) -o inlemke 

depend:: 
	gcc -MM $(ALLOBJ:.o=.c) > $(DEPFILE)

.PHONY : clean 

clean:
	-rm -f *.o core *.exe $(DEPFILE); touch $(DEPFILE); make depend

include $(DEPFILE)
