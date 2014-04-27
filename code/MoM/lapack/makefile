include ../../../makefile.inc

OBJECTS_lib = zgels_interface.o ztrsm_interface.o zgels.o zunmlq.o xerbla.o dlamch.o ilaenv.o lsame.o zgelqf.o zgelss.o  zgeqr2.o zgeqrf.o zlange.o zlascl.o zlaset.o  zunmqr.o zlacpy.o zbdsqr.o dlascl.o zgelq2.o zlarfg.o zunmbr.o zlarfb.o zlarft.o ztrmv.o ieeeck.o ztrmm.o zungbr.o zswap.o zunm2r.o zlassq.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasq5.o dlasq6.o dcopy.o zunml2.o zdrot.o dlasrt.o dlas2.o dlasv2.o zgemv.o zscal.o zgebrd.o zlabrd.o zgebd2.o zunglq.o zladiv.o dladiv.o dlaset.o dlartg.o dznrm2.o zungqr.o zung2r.o zlarf.o zgerc.o dlapy3.o dlabad.o zcopy.o zdscal.o zlacgv.o zungl2.o zgemm.o zdrscl.o ztrsm.o zlasr.o i1mach.o

main: libZGELS # zgels.so

libZGELS: $(OBJECTS_lib)
	ar -r libZGELS.a $(OBJECTS_lib)

%.o : %.cpp
	$(CC) $(INCLUDE_PATH) $(CFLAGS) $< -o $@

%.o : %.f
	$(F77) $(F_FLAGS) $< -o $@

zgels.so: $(OBJECTS_lib)
	#f2py zgels.f dlamch.f -m zgels -h zgels.pyf --overwrite-signature
	f2py -c zgels.pyf $(OBJECTS_lib)

clean: 
	rm -f zgels.so *.o *.a *~ *.pyc

