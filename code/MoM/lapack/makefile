include ../../../makefile.inc

OBJECTS_lib = zgels_interface.o ztrsm_interface.o zgels.o zunmlq.o xerbla.o dlamch.o ilaenv.o lsame.o zgelqf.o zgelss.o  zgeqr2.o zgeqrf.o zlange.o zlascl.o zlaset.o  zunmqr.o zlacpy.o zbdsqr.o dlascl.o zgelq2.o zlarfg.o zunmbr.o zlarfb.o zlarft.o ztrmv.o ieeeck.o ztrmm.o zungbr.o zswap.o zunm2r.o zlassq.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasq5.o dlasq6.o dcopy.o zunml2.o zdrot.o dlasrt.o dlas2.o dlasv2.o zgemv.o zscal.o zgebrd.o zlabrd.o zgebd2.o zunglq.o zladiv.o dladiv.o dlaset.o dlartg.o dznrm2.o zungqr.o zung2r.o zlarf.o zgerc.o dlapy3.o dlabad.o zcopy.o zdscal.o zlacgv.o zungl2.o zgemm.o zdrscl.o ztrsm.o zlasr.o i1mach.o


main: libZGELS # zgels.so

libZGELS: $(OBJECTS_lib)
	ar -r libZGELS.a $(OBJECTS_lib)

zgels_interface.o: zgels_interface.cpp
	$(CC) $(INCLUDE_PATH) $(CFLAGS) zgels_interface.cpp
ztrsm_interface.o: ztrsm_interface.cpp
	$(CC) $(INCLUDE_PATH) $(CFLAGS) ztrsm_interface.cpp
zgels.o: zgels.f
	$(F77) $(F_FLAGS) zgels.f
zunmlq.o: zunmlq.f
	$(F77) $(F_FLAGS) zunmlq.f
xerbla.o: xerbla.f
	$(F77) $(F_FLAGS) xerbla.f
dlamch.o: dlamch.f
	$(F77) $(F_FLAGS) dlamch.f
ilaenv.o: ilaenv.f
	$(F77) $(F_FLAGS) ilaenv.f
lsame.o: lsame.f
	$(F77) $(F_FLAGS) lsame.f
zgelqf.o: zgelqf.f
	$(F77) $(F_FLAGS) zgelqf.f
zgelss.o: zgelss.f
	$(F77) $(F_FLAGS) zgelss.f
zgeqr2.o: zgeqr2.f
	$(F77) $(F_FLAGS) zgeqr2.f
zgeqrf.o: zgeqrf.f
	$(F77) $(F_FLAGS) zgeqrf.f
zlange.o: zlange.f
	$(F77) $(F_FLAGS) zlange.f
zlascl.o: zlascl.f
	$(F77) $(F_FLAGS) zlascl.f
zlaset.o: zlaset.f
	$(F77) $(F_FLAGS) zlaset.f
zunmqr.o: zunmqr.f
	$(F77) $(F_FLAGS) zunmqr.f
zlacpy.o: zlacpy.f
	$(F77) $(F_FLAGS) zlacpy.f
zbdsqr.o: zbdsqr.f
	$(F77) $(F_FLAGS) zbdsqr.f
dlascl.o: dlascl.f
	$(F77) $(F_FLAGS) dlascl.f
zgelq2.o: zgelq2.f
	$(F77) $(F_FLAGS) zgelq2.f
zlarfg.o: zlarfg.f
	$(F77) $(F_FLAGS) zlarfg.f
zunmbr.o: zunmbr.f
	$(F77) $(F_FLAGS) zunmbr.f
zlarfb.o: zlarfb.f
	$(F77) $(F_FLAGS) zlarfb.f
zlarft.o: zlarft.f
	$(F77) $(F_FLAGS) zlarft.f
ztrmv.o: ztrmv.f
	$(F77) $(F_FLAGS) ztrmv.f
ieeeck.o: ieeeck.f
	$(F77) $(F_FLAGS) ieeeck.f
ztrmm.o: ztrmm.f
	$(F77) $(F_FLAGS) ztrmm.f
zungbr.o: zungbr.f
	$(F77) $(F_FLAGS) zungbr.f
zswap.o: zswap.f
	$(F77) $(F_FLAGS) zswap.f
zunm2r.o: zunm2r.f
	$(F77) $(F_FLAGS) zunm2r.f
zlassq.o: zlassq.f
	$(F77) $(F_FLAGS) zlassq.f
dlasq1.o: dlasq1.f
	$(F77) $(F_FLAGS) dlasq1.f
dlasq2.o: dlasq2.f
	$(F77) $(F_FLAGS) dlasq2.f
dlasq3.o: dlasq3.f
	$(F77) $(F_FLAGS) dlasq3.f
dlasq4.o: dlasq4.f
	$(F77) $(F_FLAGS) dlasq4.f
dlasq5.o: dlasq5.f
	$(F77) $(F_FLAGS) dlasq5.f
dlasq6.o: dlasq6.f
	$(F77) $(F_FLAGS) dlasq6.f
dcopy.o: dcopy.f
	$(F77) $(F_FLAGS) dcopy.f
zunml2.o: zunml2.f
	$(F77) $(F_FLAGS) zunml2.f
zdrot.o: zdrot.f
	$(F77) $(F_FLAGS) zdrot.f
zrotg.o: zrotg.f
	$(F77) $(F_FLAGS) zrotg.f
dlasrt.o: dlasrt.f
	$(F77) $(F_FLAGS) dlasrt.f
dlas2.o: dlas2.f
	$(F77) $(F_FLAGS) dlas2.f
dlasv2.o: dlasv2.f
	$(F77) $(F_FLAGS) dlasv2.f
zgemv.o: zgemv.f
	$(F77) $(F_FLAGS) zgemv.f
zscal.o: zscal.f
	$(F77) $(F_FLAGS) zscal.f
zgebrd.o: zgebrd.f
	$(F77) $(F_FLAGS) zgebrd.f
zlabrd.o: zlabrd.f
	$(F77) $(F_FLAGS) zlabrd.f
zgebd2.o: zgebd2.f
	$(F77) $(F_FLAGS) zgebd2.f
zunglq.o: zunglq.f
	$(F77) $(F_FLAGS) zunglq.f
zladiv.o: zladiv.f
	$(F77) $(F_FLAGS) zladiv.f
dladiv.o: dladiv.f
	$(F77) $(F_FLAGS) dladiv.f
dlaset.o: dlaset.f
	$(F77) $(F_FLAGS) dlaset.f
dlartg.o: dlartg.f
	$(F77) $(F_FLAGS) dlartg.f
dznrm2.o: dznrm2.f
	$(F77) $(F_FLAGS) dznrm2.f
zungqr.o: zungqr.f
	$(F77) $(F_FLAGS) zungqr.f
zung2r.o: zung2r.f
	$(F77) $(F_FLAGS) zung2r.f
zlarf.o: zlarf.f
	$(F77) $(F_FLAGS) zlarf.f
zgerc.o: zgerc.f
	$(F77) $(F_FLAGS) zgerc.f
dlapy3.o: dlapy3.f
	$(F77) $(F_FLAGS) dlapy3.f
dlabad.o: dlabad.f
	$(F77) $(F_FLAGS) dlabad.f
zcopy.o: zcopy.f
	$(F77) $(F_FLAGS) zcopy.f
zdscal.o: zdscal.f
	$(F77) $(F_FLAGS) zdscal.f
zlacgv.o: zlacgv.f
	$(F77) $(F_FLAGS) zlacgv.f
zungl2.o: zungl2.f
	$(F77) $(F_FLAGS) zungl2.f
zgemm.o: zgemm.f
	$(F77) $(F_FLAGS) zgemm.f
zdrscl.o: zdrscl.f
	$(F77) $(F_FLAGS) zdrscl.f
ztrsm.o: ztrsm.f
	$(F77) $(F_FLAGS) ztrsm.f
zlasr.o: zlasr.f
	$(F77) $(F_FLAGS) zlasr.f
i1mach.o: i1mach.f
	$(F77) $(F_FLAGS) i1mach.f

zgels.so: $(OBJECTS_lib)
	#f2py zgels.f dlamch.f -m zgels -h zgels.pyf --overwrite-signature
	f2py -c zgels.pyf $(OBJECTS_lib)

clean: 
	rm -f zgels.so *.o *.a *~ *.pyc

