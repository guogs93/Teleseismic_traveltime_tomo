#
#   Makefile for program aktsurf
#
SHELL=/bin/zsh
FF = gfortran
AKF = libsun.o libtau.o ellip.o sphdist.o
TAUL = libtau.o libsun.o
all: aktsurf remodl setbrn
#
#  Create program aktsurf
#
aktsurf: aktsurf.o $(AKF)
	$(FF) -o aktsurf aktsurf.o $(AKF)
aktsurf.o: aktsurf.f
	$(FF) -c aktsurf.f
#
#  Create program remodlv
#
remodl: remodl.o $(TAUL) emtvelin.o
	$(FF) -o remodl remodl.o emtvelin.o $(TAUL)
remodl.o: remodl.f
	$(FF) -c remodl.f
#
#  Create program setbrn
#
setbrn: setbrn.o $(TAUL)
	$(FF) -o setbrn setbrn.o $(TAUL)
setbrn.o: setbrn.f
	$(FF) -c setbrn.f
#
#  Create other object files
#
emtvelin.o: emtvelin.f
	$(FF) -c emtvelin.f

libsun.o: libsun.f
	$(FF) -c libsun.f

libtau.o: libtau.f
	$(FF) -c libtau.f

ellip.o: ellip.f
	$(FF) -c ellip.f

sphdist.o: sphdist.f
	$(FF) -c sphdist.f
#
# Remove object files
#
clean:
	rm *.o
