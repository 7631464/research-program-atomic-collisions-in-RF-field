This program uses multichannel quantum defect theory to calculate the scattering between two Rb87 identical bosons in a single color radio frequency field.

Please read the Makefile to see how to compile the program. 
The entire code includes 7 source files: 

Calcdress[16-18]: 16 to 18 are the version numbers. This is the main code to solve the scattering equation. It use R-matrix method or frame transformation method to calculate the short range K matrix and then use QDT to calculate the physical S matrix. The program can also calculate multichannel molecular states. 16 use the frame transformation method. 17 scan the collision energy. 18 use the coupled channel R-matrix method. 
rfchannelsav2: subroutine that calculates the asymptotic Eigen channels with rf, using Floquet method. One floquet block is defined as follows: mf= even, n=0; mf=odd, n=1.
van_der_waals_t3: subroutines that generate the B-spline basis and calculate the reference wave functions using R-matrix method.
Hypf.SubdressRb.f90 : subroutines that calculate the hyperfine states of Rb87 and all asymptotic channels of Rb87+Rb87 collision without rf. 
contains the suboutine to calculate the hyperfine states and thresholds of two atoms.
rbrbpotential.f90: The singlet and triplet potential of Rb-Rb.  
gensub2.f90 : Subroutines of calculating the CG coefficients and 3j,6j ,9j coefficients.
resanalyzersub.f90: subroutines that analyze the resonance properties, i.e. position , width, decay width.


parameterrb87wL.inp: This is the input file that specifies parameters in each calculation,
	mfab: the Toal m_F=m_fa+m_fb. It defines a certain collision manifold without rf.
	thresh: the label of the collision channel, for example thresh=1 is the lowest collision channel in a certain mfab manifold without rf. 
	mfab and thresh together defines a specific collision channel without rf. With the inclusion of rf, the program with choose a dressed channel with the most component in (mfab, thresh) state. 
	iflag: switch between coupled channel calculation(iflag=1) and frame transformation (iflag=0).
	iflagBorwL: switch between scanning static B field(iflagBorwL=0) and scanning rf frequency(iflagBorwL=1), or scanning both(iflagBorwL=2);
	iflagmol: The program will calculate scattering physics (iflagmol=0) or molecular states (iflagmol=1)
	Einput: the collision energy in atomic units;
	c6,c8,c10: the long range parameters of the interatomic potential, in frame transformation calculation c8 and c10 are chosen to be zero;
	m1, m2: the mass of the two atoms in atomic units.
	Rmatch: the matching distance in atomic units. The short range solution are matched to the reference wave functions at this distance.
	Binitial, Bfinal, Numpoints: specify how to scan the static B field in units of Gauss. When iflagBorwL=1, Binitial is used.
	wLinital, wLfinal, NumwLpoints: specify how to scan the rf frequency in units of MHz. When iflagBorwL=0, wLinitial is used.
	Bx: the B field amplitude of the rf field in units of Gauss;
	wL: not used;
	nblock: the number of Floquet blocks used in the calculation. nblock=3 is usually enough to converge.
	

alist.dat: QDT parameter A as a function sqrt(E);
glist.dat: QDT parameter G as a function sqrt(E);
etalist.dat: QDT parameter eta as a function sqrt(E);
gamalist.dat: QDT parameter gamma as a function sqrt(-E);

