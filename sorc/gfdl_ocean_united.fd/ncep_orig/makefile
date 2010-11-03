CF = ncepxlf
FFLAGS = -O3 -qstrict -qnosave  -qsmp=auto
LIBS=           -L/nwprod/lib -lw3_d 

LDFLAGS=       

.f.o: comblk.h
	$(CF) -c $(FFLAGS) $<

gfdl_ocean_united: ocean.o MAIN.o 
	$(CF) $(LDFLAGS) $(FFLAGS) $(LIBS) -o gfdl_ocean_united ocean.o MAIN.o $(LIBS)

