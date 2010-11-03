SHELL = /bin/sh
EXEC =  `pwd`
#EXEC = .
FC =  pathf90 
#FC =   mpif90 
CMD =	$(EXEC)/gfdl_ocean.united
FFLAGS = -Ofast -byteswapio
FFLAGS = -Ofast 
FFLAGS = -O2 -OPT:Ofast -fno-math-errno -DDEC_ALPHA
#FFLAGS =
SRC1 = MAIN.f ocean.f newwind.f	

BJ1 = ${SRC1:.f=.o}

.SUFFIXES: .f .o

.f.o:
	$(FC) $(FFLAGS) -c $<

all: $(CMD)

$(CMD):       $(BJ1)
	$(FC) $(FFLAGS) -o $@ $(BJ1)

$(OBJ):      
	$(FC) $(FFLAGS) -c $<

clean:
	-rm -f $(OBJ) 

allclean: clean all
