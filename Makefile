FC = ifort

CMPLFLAGS = -O3 #-xW -tpp7

LIB_CARTER = -llapack -lblas  #-ffree-line-length-none

LIBRARIES = $(LIB_CARTER)

SOURCE = calcdress18.f90 van_der_waals_t3.f90 rfchannelsav2.f90 rbrbpotential.f90 Hypf2b.SubRdressRb.f90 gensub2.f90 resanalyzersub.f90

OBJECT=$(SOURCE) 

Rmatrix_calc: $(OBJECT)
	$(FC) $(CMPLFLAGS) -o Rmatrix_calc $(OBJECT) $(LIBRARIES)

.f90.o:              
	$(FC) $(CMPLFLAGS) -c $< 
.f.o: 
	$(FC) $(CMPLFLAGS)  -c $<

clean:
	rm $(DEL)
