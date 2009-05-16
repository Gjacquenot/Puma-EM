include ../../../../makefile.inc

OBJECTS_lib = d1mach.o zbesh_interface.o zbinu.o zbunk.o zmlt.o zseri.o zuchk.o zuni2.o zunk2.o dgamln.o xerror.o zacon.o zbesh.o zbknu.o zdiv.o zlog.o zrati.o zshch.o zunhj.o zunik.o zuoik.o i1mach.o zabs.o zairy.o zbesj.o zbuni.o zexp.o zmlri.o zs1s2.o zsqrt.o zuni1.o zunk1.o zwrsk.o zacai.o zasyi.o zkscl.o 

main: libAMOS


libAMOS: $(OBJECTS_lib)
	ar -r libAMOS.a $(OBJECTS_lib)

zbesh_interface.o: zbesh_interface.cpp
	$(CC) $(CFLAGS) zbesh_interface.cpp
zabs.o: zabs.f
	$(F77) $(F_FLAGS) zabs.f
zacai.o: zacai.f
	$(F77) $(F_FLAGS) zacai.f
zacon.o: zacon.f
	$(F77) $(F_FLAGS) zacon.f
zairy.o: zairy.f
	$(F77) $(F_FLAGS) zairy.f
zasyi.o: zasyi.f
	$(F77) $(F_FLAGS) zasyi.f
zbesh.o: zbesh.f
	$(F77) $(F_FLAGS) zbesh.f
zbesj.o: zbesj.f
	$(F77) $(F_FLAGS) zbesj.f
zbinu.o: zbinu.f
	$(F77) $(F_FLAGS) zbinu.f
zbknu.o: zbknu.f
	$(F77) $(F_FLAGS) zbknu.f
zbuni.o: zbuni.f
	$(F77) $(F_FLAGS) zbuni.f
zbunk.o: zbunk.f
	$(F77) $(F_FLAGS) zbunk.f
zdiv.o: zdiv.f
	$(F77) $(F_FLAGS) zdiv.f
zexp.o: zexp.f
	$(F77) $(F_FLAGS) zexp.f
zkscl.o: zkscl.f
	$(F77) $(F_FLAGS) zkscl.f
zlog.o: zlog.f
	$(F77) $(F_FLAGS) zlog.f
zmlri.o: zmlri.f
	$(F77) $(F_FLAGS) zmlri.f
zmlt.o: zmlt.f
	$(F77) $(F_FLAGS) zmlt.f
zrati.o: zrati.f
	$(F77) $(F_FLAGS) zrati.f
zsqrt.o: zsqrt.f
	$(F77) $(F_FLAGS) zsqrt.f
zs1s2.o: zs1s2.f
	$(F77) $(F_FLAGS) zs1s2.f
zseri.o: zseri.f
	$(F77) $(F_FLAGS) zseri.f
zshch.o: zshch.f
	$(F77) $(F_FLAGS) zshch.f
zuchk.o: zuchk.f
	$(F77) $(F_FLAGS) zuchk.f
zunhj.o: zunhj.f
	$(F77) $(F_FLAGS) zunhj.f
zuni1.o: zuni1.f
	$(F77) $(F_FLAGS) zuni1.f
zuni2.o: zuni2.f
	$(F77) $(F_FLAGS) zuni2.f
zunik.o: zunik.f
	$(F77) $(F_FLAGS) zunik.f
zunk1.o: zunk1.f
	$(F77) $(F_FLAGS) zunk1.f
zunk2.o: zunk2.f
	$(F77) $(F_FLAGS) zunk2.f
zuoik.o: zuoik.f
	$(F77) $(F_FLAGS) zuoik.f
zwrsk.o: zwrsk.f
	$(F77) $(F_FLAGS) zwrsk.f
dgamln.o: dgamln.f
	$(F77) $(F_FLAGS) dgamln.f
gamln.o: gamln.f
	$(F77) $(F_FLAGS) gamln.f
i1mach.o: i1mach.f
	$(F77) $(F_FLAGS) i1mach.f
d1mach.o: d1mach.f
	$(F77) $(F_FLAGS) d1mach.f
xerror.o: xerror.f
	$(F77) $(F_FLAGS) xerror.f


clean: 
	rm -f *.o libAMOS.a *~ main
