
nvfortran -c -O3 -fpic -stdpar=gpu -shared \
	-o liblwcuf.so ../kinds.f90 ../sim.F90 ../liblwcuf.CUF