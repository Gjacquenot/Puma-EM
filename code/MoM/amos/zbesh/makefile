include ../../../../makefile.inc

OBJECTS_lib = d1mach.o zbesh_interface.o zbinu.o zbunk.o zmlt.o zseri.o zuchk.o zuni2.o zunk2.o dgamln.o xerror.o zacon.o zbesh.o zbknu.o zdiv.o zlog.o zrati.o zshch.o zunhj.o zunik.o zuoik.o i1mach.o zabs.o zairy.o zbesj.o zbuni.o zexp.o zmlri.o zs1s2.o zsqrt.o zuni1.o zunk1.o zwrsk.o zacai.o zasyi.o zkscl.o 

main: libAMOS


libAMOS: $(OBJECTS_lib)
	ar -r libAMOS.a $(OBJECTS_lib)

%.o : %.cpp
	$(CC) $(INCLUDE_PATH) $(CFLAGS) $< -o $@

%.o : %.f
	$(F77) $(F_FLAGS) $< -o $@

clean: 
	rm -f *.o libAMOS.a *~ main
