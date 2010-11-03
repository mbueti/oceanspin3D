# Biju Thomas ... on 05/17/2007
SHELL = /bin/sh
SRC = MAIN.f ocean.f
OBJ = ${SRC:.f=.o}
FC  = ifort -w90 -w95 -cm -nowarn 
LDFLAGS =   
CMD =   ./gfdl_ocean.united
FFLAGS = -ftz -stack_temps -safe_cray_ptr -O2 -convert big_endian
.SUFFIXES: .o .f .f90

.f.o:
	$(FC) $(FFLAGS) -c $*.f

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build field.
#
all:	$(CMD)
$(CMD):	$(OBJ)
	$(FC) $(FFLAGS) -o $(@) $(OBJ) 

clean:
	rm -f $(OBJ) 
