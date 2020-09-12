include ../../../makefile.inc

ZSOURCESF = ZGMRESREVCOM.f zaxpy.f zrotg.f zdotc.f dcabs1.f ztrsv.f xerbla.f zgemv.f zscal.f dznrm2.f ZSTOPTEST2.f zcopy.f lsame.f

OBJECTS_lib = iterative.o

iterative: $(ZOURCESF)
	#f2py2.4 $(ZSOURCESF) -m zgmresrevcom -h zgmresrevcom.pyf --overwrite-signature
	f2py2.4 -c zgmresrevcom.pyf $(ZSOURCESF)


libITERATIVE: $(OBJECTS_lib)
	ar -r libITERATIVE.a $(OBJECTS_lib)
iterative.o: iterative.cpp
	$(MPICC) $(INCLUDE_PATH) $(CFLAGS) iterative.cpp

clean:
	rm -f *.o libITERATIVE.a *~
