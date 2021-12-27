IDIR=./include
SDIR=./src
CC=mpicc
CFLAGS=-O3 -I$(IDIR)

ODIR=./obj

LIBS=-lm

_DEPS = ut.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ut.o blur.mpi.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) 

blur: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o


