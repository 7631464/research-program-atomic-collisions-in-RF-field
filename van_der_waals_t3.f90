!==============================================================================================

subroutine phiphase(rinit,r0,l,energy,c6,c8,c10,ttanphi)
!This subroutine calculates the zero energy phi phases
!uses the rrmatrix and the wkbboundcond subroutines
!Input: rinit, r0, l, energy,c6,c8,c10
!output: ttanphi, phi
!Everything is calculated in van der Waals length scale units.
!rinit is the initial point of the grind, in Ruzic's converntion rinit=R_x
!r0 is the final point of the grid. Typically you need a few VdW length scale units, i.e r0~2-10
!l is the angular momentum and energy is typically setted equal to zero.

Implicit double precision (a-h,o-z)
integer*4, allocatable ::		lobound(:),upbound(:)
integer*4 				i,j,k,left,right,NumBeta,lpts,order,np,matdim,ii,ij
real*8 					rinit,r0,energy,field,beta,l,momentum,energ,ml,Rampl,Sampl,cosgamma,singamma
real*8 					lomom,phib,pi,flast,glast,regeta,irregeta,irregscat,logderf(2),logderg(2)
real*8					c6,c8,c10,b6,sc8,sc10,fbc,gbc,fpbc,gpbc,r0x,v,x,vm,dervflast,dervglast,ttanphi,djm,bjm,sc6
real*8					phi,wros
double precision, allocatable :: 	ub(:,:),pts(:),xleg(:),wleg(:), weight(:),surfaceL(:), surface(:),derfsol(:),dergsol(:)
double precision, allocatable :: 	sectors(:),acoef(:),bcoef(:),dat(:,:),uf(:),ug(:),grid(:),fsol(:),gsol(:),duf(:),dug(:)
double precision, allocatable ::	derbs(:,:),kesav(:,:),pesav(:,:),oversav(:,:)
Real*8					BJ(0:250),DJ(0:250),BY(0:250),DY(0:250)

pi=dacos(-1.d0)


!========================================
!NUMERICAL PARAMETERS FOR THE RMATRIX AND  
! GAULEG AND MAKEBASIS SUBROUTINES
!========================================
    left=2									!PARAMTERS THAT DEFINE THE NUMBER OF REACTION ZONES, IE LEFT=RIGHT=2
    right=2									!MEANS TWO REACTION ZONES ONE ON THE L.H.S. AND THE OTHER ON R.H.S.
    order=12									!THE ORDER+1 OF THE B-SPLINES
    np=1000									!NUMBER OF BASIS FUNCTIONS
    lpts=64									!NUMBER OF LEGENDRE POINTS
    matdim=np+order+left/2+right/2-3						!THE OVERALL MATRIX DIMENSIONALITY NEEDED IN RMATRIX SUBROUNTINE
!========================================


!===============================================
!MAKE B-SPLINE BASIS FUNCTIONS AND LEGENDRE MESH
!===============================================
allocate(lobound(matdim))
allocate(upbound(matdim))
allocate(ub((order+1)*lpts,matdim*2),pts((np-1)*lpts),xleg(lpts),wleg(lpts), weight((np-1)*lpts))
allocate(surface(2*matdim),surfaceL(2*matdim),derbs((order+1)*lpts,matdim))
Allocate(sectors(matdim))

 call GAULEG(-1.d0,1.d0,xleg,wleg,lpts)												
 !print*,'CALL MAKEBASIS ROUTINE'
 call makebasis(ub,lobound,upbound,matdim,order,left,right,np,r0,rinit,lpts,xleg,&						 
 		&wleg,weight,pts,surfaceL,surface,sectors,6)
isav=0

!------------------------------
!CALCULATE THE SURFACE TERMS
  do j=1,matdim
      If (j.eq.1) then
	  surfaceL(1)=1.d0
	  surface(1)=0.d0
      elseif (j.eq.matdim) then
	  surfaceL(matdim)=0.d0
	  surface(matdim)=1.d0
      else
	  surfaceL(j)=0.d0
	  surface(j)=0.d0
      endif
  enddo
!------------------------------

!Calculation of the the WKB boundary conditions F_hat and G_hat with the phase phi set it to zero
ttanphi=0.d0
call wkbboundcond(l,rinit,c8,c10,energy,ttanphi,fbc,gbc,fpbc,gpbc)
!end of calculation

logderf(1)=fpbc
logderf(2)=fbc
logderg(1)=gpbc
logderg(2)=gbc

allocate(uf(matdim),ug(matdim),duf(matdim),dug(matdim),grid(lpts*matdim),fsol(lpts*matdim),gsol(lpts*matdim)&
&,derfsol(lpts*matdim),dergsol(lpts*matdim),kesav(matdim,matdim),pesav(matdim,matdim),oversav(matdim,matdim))
!    write(*,*) 'Call R-matrix Subroutine'
    call rrmatrix(rinit,r0,energy,field,l,c6,c8,c10,beta,logderf,logderg,left,right,&
    &lpts,np,order,matdim,ub,derbs,lobound,upbound,pts,weight,surface,surfaceL,&
    &Rampl,Sampl,cosgamma,singamma,flast,dervflast,glast,dervglast,uf,ug,duf,dug,kesav,pesav,oversav,isav)

!calculate the zero energy functions \chi_{+},\chi_{-}
v=l/2.d0+0.25d0
x=0.5d0/ r0**2.d0

call JYV(V,X,VM,BJ,DJ,BY,DY)
!The Wronskian of \chi_{-} and \chi_{+} functions, ie W(X_{-},X_{+})
wros=dsqrt(r0)*bj(int(v))*((0.5d0*by(int(v))/dsqrt(r0)+dsqrt(r0)*(-1.d0/ r0**3.d0)*dy(int(v))))&
&-dsqrt(r0)*by(int(v))*((0.5d0*bj(int(v))/dsqrt(r0)+dsqrt(r0)*(-1.d0/ r0**3.d0)*dj(int(v))))
!-------------------------------------------------------
!end of \chi calculations

!calculation of the phase phi that goes in the definition of the WKB boundary conditions of F_hat and G_hat 
!functions
 ttanphi=-(dsqrt(r0)*bj(int(v))*dervglast-(0.5d0*bj(int(v))/dsqrt(r0)+dsqrt(r0)*(-1.d0/ r0**3.d0)*dj(int(v)))*glast)&
 &/(dsqrt(r0)*bj(int(v))*dervflast-(0.5d0*bj(int(v))/dsqrt(r0)+dsqrt(r0)*(-1.d0/ r0**3.d0)*dj(int(v)))*flast)
! End of the calculation

deallocate(lobound,upbound,ub,derbs,pts,xleg,wleg)
deallocate(weight,surface,surfaceL,sectors)
deallocate(uf,ug,duf,dug,grid,fsol,gsol,derfsol,dergsol)

return
end

subroutine wkbboundcond(l,rinit,c8,c10,energy,ttanphi,fbc,gbc,fpbc,gpbc)
!*****************************
!COMMENTS
!*****************************
!In order to avoid the direct computation of the phi phase
!the reference wavefunctions f_hat and g_hat from Ruzic's paper here are scaled by cosphi.
!Namely, my F_hat and G_hat are equal to F_hat=f_hat/cosphi and G_hat=g_hat/cosphi
!The reason of using this definition is that tanphi is a multivalue function.
!Thus in order to compute the phi phase one has to specify the quadrant of the computed phi phase.
!Any mistake in the quadrant of the phase phi might yield wrong qdt parameters. Therefore with the definition that I used
!I avoid any type of such mistakes.
!In addition, the equations in Ruzic's paper for the computation of the qdt parameters G, taneta and tangamma are unchanged.
!However, the qdt parameter A has to be calculated by a different equation which elliminates the overall scale factor 1/cosphi from the
!reference functions F_hat and G_hat.
!for this puprose I used the relation f/g=A (f_hat/(g_hat+gqdt*f_hat)), which is derived by eq.9 in Ruzic's paper.
Implicit double precision (a-h,o-z)
integer*4 			i,j,N
real*8				l, rinit,r0x,expo, c8,c10,energy,fbc,gbc,fpbc,gpbc
real*8				phaseint, initvalue,pi,sum0,k,derk,integrand,r,phi,ttanphi
double precision, allocatable::	xleg(:),wleg(:)

pi=dacos(-1.d0)
!phase integral set it to zero
phaseint=0.d0

!Local momentum and its derivative with Langer correction
  k=energy+1.d0/rinit**6.d0+c8/rinit**8.d0+c10/rinit**10.d0-((l+1.d0/2.d0)**2.d0)/rinit**2.d0
if (k>0d0) then
k=dsqrt(k)
else
write(*,*)"k<0",k
stop
end if
  derk=-(0.5d0/k)*(6.d0/rinit**7.d0+8.d0*c8/rinit**9.d0+10.d0*c10/rinit**11.d0-2.d0*((l+1.d0/2.d0)**2.d0)/rinit**3.d0)
!Local momentum and its derivative without langer correction
!    k=dsqrt(energy+1.d0/rinit**6.d0+c8/rinit**8.d0+c10/rinit**10.d0-l*(l+1)/rinit**2.d0)
!    derk=-(0.5d0/k)*(6.d0/rinit**7.d0+8.d0*c8/rinit**9.d0+10.d0*c10/rinit**11.d0-2.d0*l*(l+1)/rinit**3.d0)
!   
  fbc=dsqrt(1.d0/(k))*(dsin(phaseint)+dcos(phaseint)*ttanphi)
  gbc=-dsqrt(1.d0/(k))*(dcos(phaseint)-dsin(phaseint)*ttanphi)
  fpbc=-0.5d0*(derk/k)*fbc-k*gbc
  gpbc=-0.5d0*(derk/k)*gbc+k*fbc
return
end

!================================================================================================

subroutine solutregirrg(matdim,np,order,pts,lpts,lobound,upbound,ub,derbs,uf,ug,dug,duf,grid,fsol,gsol,derfsol,dergsol)
Implicit double precision (a-h,o-z)
integer*4 			j,jp,ip,order,ir,i,jo,iimin,iimax,lobound(matdim),upbound(matdim),count0,lpts,matdim,np,ij
real*8 				sum0,sum1,ub((order+1)*lpts,matdim*2),grid(lpts*matdim),fsol(lpts*matdim),gsol(lpts*matdim)
real*8				derbs((order+1)*lpts,matdim*2)
real*8				uf(matdim),ug(matdim),duf(matdim),dug(matdim),sum2,sum3,derfsol(lpts*matdim),dergsol(lpts*matdim)
real*8 				rad,pts((np-1)*lpts),summa(lpts*matdim)
double precision, allocatable:: allsplin(:,:),alldersplin(:,:)

allocate(allsplin(lpts*matdim,matdim),alldersplin(lpts*matdim,matdim))

do ij=1,lpts*matdim
summa(ij)=0.d0
enddo

do j=1,matdim
 do jp=1,matdim 
 iiMIN=MAX(lpts*(lobound(j)-1)+1,lpts*(lobound(jp)-1)+1)
 iiMAX=MIN((upbound(j))*lpts,(upbound(jp))*lpts)
 count0=1+(jp-1)*lpts     
     do ir=iimin,iimax
	    i=ir-(lobound(j)-1)*lpts
  
	    rad=pts(ir)

 		      if (ir.ge.(jp-1)*lpts.and. ir.le.jp*(lpts)) then
		      grid(count0)=rad
		      endif
		     
		      if (jp.le.order+1) then
		      allsplin(i,jp)=ub(i,jp)
		      alldersplin(i,jp)=derbs(i,jp)
		      else		   
		      allsplin(i+(jp-order-1)*lpts,jp)=ub(i,jp)
		      alldersplin(i+(jp-order-1)*lpts,jp)=derbs(i,jp)
		      endif
	  count0=1+count0
     enddo
  enddo
enddo

      do i=1,lpts*(np-1)
	  sum0=0.d0
	  sum1=0.d0
	  sum2=0.d0
	  sum3=0.d0
	      do j=1,matdim
		  sum0=sum0+allsplin(i,j)*uf(j)
		  sum1=sum1+allsplin(i,j)*ug(j)
		  sum2=sum2+alldersplin(i,j)*duf(j)
		  sum3=sum3+alldersplin(i,j)*dug(j)
	      enddo
	  fsol(i)=sum0
	  gsol(i)=sum1
	  derfsol(i)=sum2
	  dergsol(i)=sum3
      enddo      
deallocate(allsplin)
return
end

!==============================================================================================

Subroutine rrmatrix(rinit,r0,energy,field,l,c6,c8,c10,beta,logderf,logderg,left,right,&
&lpts,np,order,matdim,ub,derbs,lobound,upbound,pts,weight,surface,surfaceL,&
&Rampl,Sampl,cosgamma,singamma,flast,dervflast,glast,dervglast,uf,ug,duf,dug,kesav,pesav,oversav,isav)
implicit double precision (a-h,o-z)
character :: 			V,N,U
integer*4, allocatable ::	ipv(:)
integer*4 			lobound(matdim),upbound(matdim)
integer*4 			matdim,order,left,right,np,lpts,i,rj,ri,iimax,iimin,jp,ir,j,ip,k,info,jl,count0,count1,kl,ku,jo
INTEGER*4 			INFO1,info2,info4,info5,info3
real*8 				ub((order+1)*lpts,matdim*2),surface(matdim),surfaceL(matdim),sectors(matdim),t1,t2,derbs((order+1)*lpts,matdim)
real*8 				pts((np-1)*lpts),weight((np-1)*lpts),uf(matdim),ug(matdim),sc8,sc10,dug(matdim),duf(matdim)
Real*8				rad,u1,u2,du1,du2,pot,sum1,sum0,sum2,sum3,det,w11,w12,w22,theta,uL1,uL2,uR1,uR2,gammall,dervglast
real*8 				gammarr,lambdarr,lambdall,eigen1,eigen2,field, beta,rinit,r0,energy,l,momentum,singamma, cosgamma,dervflast
real*8				c6,c8,c10
real*8 				logderf(2),logderg(2),gomega(2,2),detgomega,localmom,dervlocalmom,tanphi,pi,tanphi2,Rampl,Sampl,flast,glast
double precision, allocatable:: KE(:,:),PE(:,:),OVERLAP(:,:),gama(:,:),lambda(:,:)
double precision, allocatable:: evalr(:),b(:),work(:),work1(:),work2(:),gcr(:),gcl(:),omega(:,:)
double precision, allocatable::	gammarc(:),gammacc(:,:),gammacl(:),gammalc(:),gammacr(:),invgammacc(:,:)
double precision, allocatable:: logder(:,:),cg(:),cf(:),vecmat(:,:),work4(:),work5(:),AB(:,:),bs(:,:)
double precision:: kesav(matdim,matdim),pesav(matdim,matdim),oversav(matdim,matdim)
double precision, external::vatsup_unret,vxssgp_unret
pi=dacos(-1.d0)

call CPU_TIME(t1)



allocate(KE(matdim,matdim),PE(matdim,matdim),OVERLAP(matdim,matdim))
allocate(gama(matdim,matdim),lambda(matdim,matdim))

if (isav==0) then
do j=1,matdim									!CALCULATION OF THE POTENTIAL AND KINETIC TERM OF THE HAMILTONIAN
    do jp=1,matdim								!ON A B-SPLINE BASIS EMPLOYING THE VARIATIONAL PRINCIPAL
	OVERLAP(j,jp)=0.d0							!OVERLAP CORRESPONDS TO THE OVERLAP OF THE WAVEFUNCTION
	KE(j,jp)=0.d0
	PE(j,jp)=0.d0
	sum0=0.d0
	sum1=0.d0
	sum2=0.d0
	
	iiMIN=MAX(lpts*(lobound(j)-1)+1,lpts*(lobound(jp)-1)+1)
	iiMAX=MIN((upbound(j))*lpts,(upbound(jp))*lpts)
	count0=1+(jp-1)*lpts     
	    do ir=iimin,iimax
		    i=ir-(lobound(j)-1)*lpts
		    ip=ir-(lobound(jp)-1)*lpts

		    rad=pts(ir)
		    u1=ub(i,j)
		    u2=ub(ip,jp)
		    du1=ub(i,j+matdim)
		    du2=ub(ip,jp+matdim)
		    derbs(i,j)=du1
		    pot=(l*(l+1.d0))/rad**2.d0-(c6/rad**6.d0+c8/rad**8.d0+c10/rad**10.d0)

		    sum0=sum0+weight(ir)*u1*u2
		    sum1=sum1+weight(ir)*u1*u2*(-pot)
		    sum2=sum2+weight(ir)*du1*du2
	    enddo
	count1=1+count1    
	KE(j,jp)=-sum2
	PE(j,jp)=sum1
	OVERLAP(j,jp)=sum0
	gama(j,jp)=KE(j,jp)+PE(j,jp)+energy*overlap(j,jp)
	lambda(j,jp)=surface(j)*surface(jp)+surfaceL(j)*surfaceL(jp)
    enddo
enddo
KESAV=ke
PESAV=pe
oversav=overlap
deallocate(KE,PE,OVERLAP)
isav=1

elseif(isav==1) then

KE=KESAV
PE=PESAV
OVERLAP=oversav
gama=ke+pe+energy*overlap
do j=1,matdim
do jp=1,matdim
lambda(j,jp)=surface(j)*surface(jp)+surfaceL(j)*surfaceL(jp)
end do
end do
deallocate(KE,PE,OVERLAP)
else
write(*,*)" isav is incorrect"
stop
end if




    !====================================================
    !THE PARTIONING OF THE GAMMA AND LAMBDA MATRICES INTO
    !OPEN AND CLOSED "CHANNELS"
    !====================================================
	gammall=gama(1,1)									
	gammarr=gama(matdim,matdim)								
	lambdall=lambda(1,1)
	lambdarr=lambda(matdim,matdim)
    deallocate(lambda)

    allocate(gammarc(matdim-2),gammacr(matdim-2),gammalc(matdim-2),gammacl(matdim-2)&
    &,invgammacc(matdim-2,matdim-2))

	do j=1,matdim-2
	    gammalc(j)=gama(1,j+1)
	    gammacl(j)=gama(j+1,1)
	enddo
	
	do j=1,matdim-2
	    gammarc(j)=gama(matdim,j+1)
	    gammacr(j)=gama(j+1,matdim)
	enddo


	do j=1,matdim-2
	    do jp=1,matdim-2
		invgammacc(j,jp)=gama(1+j,jp+1)
	    enddo
	enddo
      
    deallocate(gama)
    !=============================================

    !=============================================
    !BANDED LU SOLVER FOR THE INVERSION OF GAMMACC
    !MATRIX
    !=============================================
	kl=order
	ku=order
    allocate(work1(matdim-2),work2(5*matdim),ab(2*kl+ku+1,matdim-2))	

	do j=1, matdim-2									!BANDED STORAGE OF THE GAMMACC MATRIX INTO AB 
	    do i=max(1,j-ku),min(matdim-2,j+kl)
		AB(kl+ku+1+i-j,j) = invgammacc(i,j)
	    enddo
	enddo

    deallocate(invgammacc)
    allocate(bs(matdim-2,matdim-2),ipv(matdim-2))
	
	do j=1,matdim-2
	    do i=1,matdim-2
		if (i.eq.j) then
		    bs(i,j)=1.d0									!THE IDENTITY MATRIX
		else
		    bs(i,j)=0.d0
		endif
	    enddo
	enddo
	

    call DGBSV( matdim-2, KL, KU, matdim-2,ab, 2*kl+ku+1,ipv,bs,matdim-2,INFO3)  			!LINEAR SOLVER SUBROUTINE
	if (abs(info3).gt. 1.d-16) then
	  print*, 'ROUTINE DGBSV: INFO=/=0',info3
	  stop
	ENDIF
    deallocate(ab)

    Allocate(gammacc(matdim-2,matdim-2))
	do i=1,matdim-2
	  do j=1,matdim-2
	    gammacc(i,j)=bs(i,j)									!STORING THE INVERTED GAMMACC MATRIX
	  enddo
	enddo
    deallocate(bs,ipv)
    !==============================================================

    
    !===============================================================
    !CALCULATION OF THE OMEGA MATRIX
    !================================================================
    allocate(gcr(matdim-2),gcl(matdim-2),omega(2,2))

      do j=1,matdim-2
	  sum0=0.d0
	  sum1=0.d0
	      do jp=1,matdim-2
		  sum0=sum0+gammacc(j,jp)*gammacl(jp)
		  sum1=sum1+gammacc(j,jp)*gammacr(jp)
	      enddo
	  gcl(j)=sum0
	  gcr(j)=sum1
      enddo
      
    deallocate(gammacc)

      sum0=0.d0
      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
	  do j=1,matdim-2
		sum0=sum0+gammalc(j)*gcl(j)
		sum1=sum1+gammalc(j)*gcr(j)
		sum2=sum2+gammarc(j)*gcl(j)
		sum3=sum3+gammarc(j)*gcr(j)
	  enddo
      omega(1,1)=gammall-sum0
      omega(1,2)=-sum1
      omega(2,1)=-sum2
      omega(2,2)=gammarr-sum3
    !=======================================================

    !========================================================
    !THE EIGENVALUES AND THE EIGENVECTORS OF THE OMEGA MATRIX
    !========================================================

    Allocate(evalr(2),work(8*matdim), vecmat(2,2))

    eigen1=0.5d0*(omega(1, 1) + omega(2, 2) -&							!THE EIGENVALUES OF THE OMEGA MATRIX
	  &dsqrt(omega(1, 1)**2.d0 + 4.d0*omega(1, 2)*omega(2, 1) -&
	  &2.d0*omega(1, 1)*omega(2, 2) + omega(2, 2)**2.d0))
    eigen2=0.5d0*(omega(1, 1) + omega(2, 2) +&
	&dsqrt(omega(1, 1)**2.d0 + 4.d0*omega(1, 2)*omega(2, 1) -&
	&2.d0*omega(1, 1)*omega(2, 2) + omega(2, 2)**2.d0)) 

    evalr(1)=eigen1
    evalr(2)=eigen2
    vecmat(1,1)=1.d0 										!THE EIGENVECTORS OF THE OMEGA MATRIX
    vecmat(1,2)=1.d0
    vecmat(2,1)=((-omega(1,1) + omega(2,2) - dsqrt(omega(1,1)**2.d0+4.d0*omega(1,2)*omega(1,2)&
		&-2.d0*omega(1,1)*omega(2,2)+omega(2,2)**2.d0))/(2.d0*omega(1,2)))

    vecmat(2,2)=((-omega(1,1) + omega(2,2) + dsqrt(omega(1,1)**2.d0+4.d0*omega(1,2)*omega(1,2)&
		&-2.d0*omega(1,1)*omega(2,2)+omega(2,2)**2.d0))/(2.d0*omega(1,2)))

    ul1=vecmat(1,1)/dsqrt(vecmat(2,1)**2.d0+1.d0)						!TWO LINEAR INDEPENDENT AND NORMALISED VECTORS 
    ul2=vecmat(1,2)/dsqrt(vecmat(2,2)**2.d0+1.d0)						!FOR THE LHS OF THE GRID (UL1,UL2) AND TWO FOR  
    ur1=vecmat(2,1)/dsqrt((vecmat(2,1)**2.d0+1.d0))						!THE RHS, IE (UR1, UR2)
    ur2=vecmat(2,2)/dsqrt(vecmat(2,2)**(2.d0)+1.d0)

    !====================================================================================================
    
    allocate(logder(2,2),work4(2),work5(10),cf(2),cg(2))

    det=(evalr(1)-evalr(2))*uL1*uL2
    logder(1,1)=uL2/det
    logder(1,2)=-evalr(2)*uL2/det
    logder(2,1)=-uL1/det
    logder(2,2)=evalr(1)*uL1/det

    cf=matmul(logder,logderf)
    cg=matmul(logder,logderg)


!=================================================================
!Linear combination of the two eigenvectors of Omega matrix
!Calculate the CC channels coefficients for the f and g solutions
!================================================================
    uf(1)=cf(1)*ul1+cf(2)*ul2
    ug(1)=cg(1)*ul1+cg(2)*ul2
    
    uf(matdim)=cf(1)*ur1+cf(2)*ur2
    ug(matdim)=cg(1)*ur1+cg(2)*ur2
    
    duf(1)=evalr(1)*cf(1)*ul1+evalr(2)*cf(2)*ul2
    dug(1)=evalr(1)*cg(1)*ul1+evalr(2)*cg(2)*ul2
    
    duf(matdim)=-evalr(1)*cf(1)*ur1-evalr(2)*cf(2)*ur2
    dug(matdim)=-evalr(1)*cg(1)*ur1-evalr(2)*cg(2)*ur2
    do j=1,matdim-2
	uf(j+1)=-gcr(j)*(cf(1)*ur1+cf(2)*ur2)-gcl(j)*(cf(1)*ul1+cf(2)*ul2)
	ug(j+1)=-gcr(j)*(cg(1)*ur1+cg(2)*ur2)-gcl(j)*(cg(1)*ul1+cg(2)*ul2)
	duf(j+1)=-gcr(j)*(duf(matdim))-gcl(j)*(duf(1))
	dug(j+1)=-gcr(j)*(dug(matdim))-gcl(j)*(dug(1))
    enddo
     
      flast=uf(matdim)									!THE REGULAR SOLUTION EVALUATED AT R0 (R.H.S. GRID POINT)
      dervflast=-evalr(1)*cf(1)*ur1-evalr(2)*cf(2)*ur2
      glast=ug(matdim)									!THE IRREGULAR SOLUTION EVALUATED AT R0 (R.H.S. GRID POINT)	
      dervglast=-evalr(1)*cg(1)*ur1-evalr(2)*cg(2)*ur2
!====================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++
!Calculate R & S amplitudes, cosgamma, singamma
!++++++++++++++++++++++++++++++++++++++++++++++++++++
  !====================================================
  !WITH LANGER CORRECTION
  !====================================================
      localmom=dsqrt(0.5d0*energy-2.d0*(0.5d0*((l+0.5d0)**2.d0)/r0**2.d0-0.5d0*(1.d0-beta)/r0**1.d0-field*r0/8.d0))	!WKB LOCAL MOMENTUM CALCULATION
      dervlocalmom=- (-(field/8.d0) - ( ((l+0.5d0)**2.d0))/r0**3.d0 + (0.5d0*(1.d0 - beta))/r0**2.d0)/localmom		!DERIVATIVE OF THE LOCAL MOMENTUM
  !=======================================================
  !WITHOUT THE LANGER CORRECTION
  !=======================================================
!	localmom=dsqrt(0.5d0*energy-2.d0*(0.5d0*(l*(l+1.d0))/r0**2.d0-0.5d0*(1.d0-beta)/r0**1.d0-field*r0/8.d0))	!WKB LOCAL MOMENTUM CALCULATION
!	dervlocalmom=- (-(field/8.d0) - ( (l*(l+1.d0)))/r0**3.d0 + (0.5d0*(1.d0 - beta))/r0**2.d0)/localmom		!DERIVATIVE OF THE LOCAL MOMENTUM
  !====================================================

  tanphi=localmom/(0.5d0*dervlocalmom/localmom+(-evalr(1)*cf(1)*ur1-evalr(2)*cf(2)*ur2)/(cf(1)*ur1+cf(2)*ur2))
  tanphi2=localmom/(0.5d0*dervlocalmom/localmom+(-evalr(1)*cg(1)*ur1-evalr(2)*cg(2)*ur2)/(cg(1)*ur1+cg(2)*ur2))

  Rampl=1.d0/(dsqrt(2.d0/(pi*localmom))*dabs((tanphi/dsqrt(1.d0+tanphi**2.d0))/(cf(1)*ur1+cf(2)*ur2)))
  Sampl=1.d0/(dsqrt(2.d0/(pi*localmom))*dabs((tanphi2/dsqrt(1.d0+tanphi2**2.d0))/(cg(1)*ur1+cg(2)*ur2)))
  singamma=1.d0/(Rampl*Sampl)
  cosgamma=singamma/tanphi+(Rampl/sampl)*((cg(1)*ur1+cg(2)*ur2)/(cf(1)*ur1+cf(2)*ur2))
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++

deallocate(gammarc,gammacr,gammalc,gammacl,work1,work2&
	    &,work,evalr,gcr,gcl,omega,logder,work4,work5,cf,cg,vecmat)

call CPU_TIME(t2)
!write(*,*) 'CPU TIME', (t2-t1)/60,'MIN'

return
end 

!c23456789
!c	u   stores the splines and their first derivatives
!c	lobound,upbound   contain the 1st and last sector where the spline is non-zero	
!c	matdim=np+order+left/2+right/2-3
!c	order+1   is the order of the spline, don't ask me why
!c	left,right:  0 means u->0
!c		     1 means u'->0
!c		     2 means u is unrestricted
!c	np  is the number of knot points, or one more than the number of sectors	
!c	r0  is the box size
!c	lpts is the number of Legendre Points in each sector
!c	xleg,wleg are Legendre integration information
!c	weight is an output vector giving the Gauss-Legendre integration weightings
!c	pts is an output vector giving the values of r at which the splines were computed
!c	
!c	The storage is such that u(i,j) means the jth spline evaluated at
!c	the point pts(i+lpts*(lobound(i)-1))
!c
!c
!c
!==============================================================================================
!==============================================================================================

	subroutine makebasis(u,lobound,upbound,matdim,order,left,right,np,r0,rinit,lpts,xleg,&
	&wleg,weight,pts,surfaceL,surface,sectors,igridchoice)

	implicit double precision (a-h,o-z)

	integer*4 matdim,order,left,right,np,lpts,count0,cmax,cmin,i,j,k,kmin,kmax,li,lj,jmin,jmax,l
 	integer*4 lobound(matdim),upbound(matdim),Derv0,Derv1
!chg	double precision lobound(matdim),upbound(matdim),Derv0,Derv1
 	integer*4,allocatable :: t1(:)
!chg	double precision,allocatable :: t1(:)
	real*8 u((order+1)*lpts,matdim*2),r0,rinit,weight((np-1)*lpts),pts((np-1)*lpts)
	real*8 scale0,zero,xleg(lpts),wleg(lpts),surface(2*matdim),sectors(matdim)
	real*8 int,bspline,surfaceL(2*matdim),expo
	double precision, allocatable :: points(:),t(:)

	Derv0=0
	Derv1=1
	allocate(points(np))
if (igridchoice>1) then
	expo=dfloat(igridchoice)
	points=(/((i*(r0**(1.d0/expo)-rinit**(1.d0/expo))/dfloat(np-1) +rinit**(1.d0/expo)),i=0,np-1,1)/)
	points=points**expo
! 	points=(/((i*(dsqrt(r0)-dsqrt(rinit))/dfloat(np-1) +dsqrt(rinit)),i=0,np-1,1)/)
! 	points=points**2.d0
else if (igridchoice==1) then
points=(/((i*(r0-rinit)/dfloat(np-1) +rinit),i=0,np-1,1)/)

else
write(*,*) "igridchoice is incorrect"
stop
end if


	sectors=points!(/(points(i+1)-points(i), i=1,np-1,1)/)

	allocate(t(np+2*order))

	t(1:order)=points(1)
	t(order+1:order+np)=points
        t(order+np+1:2*order+np)=points(np)

	allocate(t1(np+2*order))

	t1(1:order)=1	
	t1(order+1:order+np)=(/(i,i=1,np,1)/)
	t1(order+np+1:2*order+np)=np
!c	write(85,9875) (t1(i),t(i),i=1,2*order+np)

   	lobound(1:matdim)=t1(2-left/2:1+matdim-left/2)
   	upbound(1:matdim)=t1(3+order-left/2:2+matdim+order-left/2)-1

!chg 	lobound(1:matdim)=t(2-left/2:1+matdim-left/2)
!chg 	upbound=t(3+order-left/2:2+matdim+order-left/2)-1
!c	write(86,*) lobound
!c	write(87,*) upbound

	do k=0,np-2  
	 scale0=(points(k+2)-points(k+1))/2.d0;
	 zero=(points(k+2)+points(k+1))/2.d0; 
	 pts(1+k*lpts:(k+1)*lpts)=zero+scale0*xleg;
	 weight(1+k*lpts:(k+1)*lpts)=scale0*wleg;   
!c	 write(88,*) k,np
!c	 write(88,9876) pts(1+k*lpts:(k+1)*lpts)
!  9876	 format(6(1x,1pd12.5))
!  9875	 format(3(1x,i5,1x,1pd12.5))
	enddo

	cmax=matdim
	cmin=1
!open(unit=17,file="bsplinebasis.dat")
	if (left==1) then
	 count0=1
	 cmin=2
	 do k=1+(lobound(count0)-1)*lpts,upbound(count0)*lpts
	  l=k-(lobound(count0)-1)*lpts
	  u(l,count0)=bspline(t,order+1,Derv0,np,1,pts(k))+bspline(t,order+1,Derv0,np,2,pts(k))
	  u(l,count0+matdim)=bspline(t,order+1,Derv1,np,1,pts(k))+bspline(t,order+1,Derv1,np,2,pts(k))
!c	  u(l,count0+2*matdim)=bspline(t,order+1,2,np,1,pts(k))+bspline(t,order+1,2,np,2,pts(k))
!write(17,*) pts(k), u(l,count0)
	 enddo
	endif

	if (right==1) then
	 cmax=matdim-1
	 count0=matdim
	 do k=1+(lobound(count0)-1)*lpts,upbound(count0)*lpts
	  l=k-(lobound(count0)-1)*lpts
	  u(l,count0)=bspline(t,order+1,0,np,np+order-1,pts(k))+bspline(t,order+1,0,np,np+order-2,pts(k))
	  u(l,count0+matdim)=bspline(t,order+1,1,np,np+order-1,pts(k))+bspline(t,order+1,1,np,np+order-2,pts(k))
!c	  u(l,count0+2*matdim)=bspline(t,order+1,2,np,np+order-1,pts(k))+bspline(t,order+1,2,np,np+order-2,pts(k))
!write(17,*) pts(k), u(l,count0)
	 enddo
	endif

!c	write(6,*) 'cmax=',cmax,cmin
!write(*,*) "cmin=",cmin,"cmax=",cmax


	do count0=cmin,cmax
	 i=count0+1-left/2
	 do k=1+(lobound(count0)-1)*lpts,upbound(count0)*lpts
	  l=k-(lobound(count0)-1)*lpts
!c	   write(83,*) count0,i,k,l,'=count0,i,k,l'
	   u(l,count0)=bspline(t,order+1,Derv0,np,i,pts(k))
!c	   write(83,*) lobound(count0),upbound(count0),pts(k),u(l,count0)
	  u(l,count0+matdim)=bspline(t,order+1,Derv1,np,i,pts(k))
!c	   write(83,*) lobound(count0),upbound(count0),pts(k),u(l,count0+matdim)
!c	  u(l,count0+2*matdim)=bspline(t,order+1,2,np,i,pts(k))
!write(17,*) pts(k), u(l,count0)
	 enddo

	  surface(count0)=0.d0
	  surface(count0+matdim)=0.d0
	  	  surfaceL(count0)=0.d0
	  surfaceL(count0+matdim)=0.d0
	  
	  if(count0 .ge. cmax-1) then	

		surface(count0) = bspline(t,order+1,Derv0,np,i,r0-1.d-14)
		surface(count0+matdim) = bspline(t,order+1,Derv1,np,i,r0-1.d-14)
	  endif
	  if(count0 .le. 2) then
		surfaceL(count0) = bspline(t,order+1,Derv0,np,i,rinit+1.d-14)
!		write(6,*) 'surfaceL(count0)=',surfaceL(count0),count0
		surfaceL(count0+matdim) = bspline(t,order+1,Derv1,np,i,rinit+1.d-14)
        endif
	enddo
	
!     write(6,*) 'done calculating in makebasis ',surface(cmax),cmax,r0
	deallocate(t1,t)
	deallocate(points)


	end

         
        double precision function bspline(t,k,jderiv,xnp,n,x)
           
        integer*4 k,jderiv,xnp,n,q,qmin,jcmin,jcmax,j
	double precision bcoef(xnp+k-1)
        double precision dr(k),dl(k),aj(k),t(*),x,EPS
	  data EPS/1.d-13/
        integer :: i=1
          
        bcoef=0.d0
        bcoef(n)=1.d0
        bspline=0.d0
           
!       write(66,*) 'in bspline now!',t(xnp+2*(k-1))-x,x

        if (x>=t(1) .and. x<=t(xnp+k) .and. jderiv<k) then
           
         if (i>2*k+xnp-2) then
          i=1
         endif
         if (t(i)<=x) then
          q=i 
         else
          q=1
         endif  
           
         do while (q<xnp+k .and. t(q)<=x)
          q=q+1
         enddo

         i=q-1

         if (k+n-1>0) then
          if (k==1 .and. jderiv==0) then
           bspline=bcoef(i)
          endif
        
          jcmin=max(1,1+k-i)
          jcmax=min(k,k+n-1)
        
          if (i-k+1<0) then
           dl(1:k-i-1)=x-t(1)
           dl(k-i:k-1)=x-t(1:i)
          else
           dl=x-t(i-k+2:i)
          endif
!   print*,'dl',dl
!   read(*,*)
          if (n-i+1<0) then
           dr(1:jcmax)=t(i+1:i+jcmax)-x
           dr(jcmax+1:k-1)=dr(jcmax)
          else
           dr(1:k-1)=t(i+1:i+k-1)-x
          endif
! print*,'dr',dr
! read(*,*)
          aj(k)=0
          aj(jcmin:jcmax)=bcoef(i-k+jcmin:i-k+jcmax)
        
          if (jderiv>0) then
           do j=1,jderiv
            aj(1:k-j)=dfloat(k-j)*(aj(2:k-j+1)-aj(1:k-j))/(dl(j:k-1)+dr(1:k-j))
           enddo
          endif
        
          if (jderiv==k-1) then
           bspline=aj(1)
          else
           do j=jderiv+1,k-1
            aj(1:k-j)=(aj(2:k-j+1)*dl(j:k-1)+aj(1:k-j)*dr(1:k-j))/(dl(j:k-1)+dr(1:k-j))
!             Print*, (dl(j:k-1)+dr(1:k-j)),1.d0/(dl(j:k-1)+dr(1:k-j)),k,j
!  !           read(*,*)
           enddo
           bspline=aj(1)
          endif
         endif
        endif
          
        end

!     upotp1=(upot(2)-upot(1))/(radmesh(2)-radmesh(1))
!     upotpn=(upot(numR)-upot(numR-1))/(radmesh(numR)-radmesh(numR-1))

!     upotp1=(upot(2)-upot(1))/(radmesh(2)-radmesh(1))
!     upotpn=(upot(numR)-upot(numR-1))/(radmesh(numR)-radmesh(numR-1))

!     call spline(radmesh,upot,numR,upotp1,upotpn,y2)
!     call splint(radmesh,upot,y2,numR,rk,eng)


      subroutine spline(x,y,n,yp1,ypn,y2)
      implicit real*8(a-h,o-z)
      dimension x(n),y(n),y2(n)
      parameter (nmax=500)
      dimension u(nmax)

      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+&
     &1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig&
     &*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end


      subroutine splint(xa,ya,y2a,n,x,y)
      implicit real*8(a-h,o-z)
      dimension xa(n),y2a(n),ya(n)

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.d0) THEN 
      PRINT*,  'bad xa input in splint'
      GOTO 2
      ENDIF
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))*(h**&
     &2.d0)/6.d0
2 CONTINUE
      return
      end
! ! c******************************************************************
! !       subroutine GetGaussFactors(File,Points,x,w)
! ! 
! ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ! c
! ! c     This subroutine retrieves the points and weights for
! ! c      Gaussian Quadrature from a given file
! ! c
! ! c     Variables:
! ! c      File		name of file containing points and 
! ! c			 weights
! ! c      Points		number of points to be used for 
! ! c			quadrature
! ! c      x		array of nodes
! ! c      w		array of weights
! ! c
! ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! !       
! !       integer Points
! !       double precision x(Points),w(Points)
! !       character*64 File
! !       
! !       integer i,tempPoints
! !       
! !       open(unit=7,file=File(1:index(File,' ')-1))
! ! 
! !        do i = 1,18
! !         read(7,*)
! !        enddo
! !        read(7,*) tempPoints
! !        do while (tempPoints .ne. Points)
! !         do i = 1,tempPoints
! !          read(7,*)
! !         enddo
! !         read(7,*) tempPoints
! !        enddo
! !  
! !        do i = 1,Points
! !         read(7,*) x(i),w(i)
! !        enddo
! !       
! !       close(unit=7)
! ! 
! !       return
! !       end
! ! 
!       double precision function BSpline(Order,Deriv,xNumPoints,xPoints,n,t,b,x)
! 
!       integer Order,Deriv,xNumPoints,n
!       double precision xPoints(*),t(*),b(xNumPoints+Order),x
! 
!       integer i
!       double precision bvalue
! 
!       b = 0.0d0
!       b(n) = 1.0d0
! 
!       BSpline = bvalue(t,b,n,Order+1,x,Deriv)
! 
!       return
!       end
! 
!       subroutine CalcBasisFuncs(Left,Right,Order,xPoints,LegPoints,xLeg,MatrixDim,xBounds,xNumPoints,Deriv,u)
! 
!       integer Left,Right,Order,LegPoints,MatrixDim,xBounds(*),xNumPoints,Deriv
!       double precision xPoints(*),xLeg(*)
!       double precision u(LegPoints,xNumPoints,MatrixDim)
! 
!       integer i,k,l,Count
!       integer, allocatable :: t(:)
!       double precision BSpline
!       double precision ax,bx,xIntScale,xScaledZero
!       double precision, allocatable :: x(:,:),tx(:),b(:)
! 
!       allocate(x(LegPoints,xNumPoints))
!       allocate(t(xNumPoints+2*Order),tx(xNumPoints+2*Order),b(xNumPoints+Order))
! 
!       do i = 1,Order
!        t(i) = 1
!        tx(i) = xPoints(1)
!       enddo
!       do i = 1,xNumPoints
!        t(i+Order) = i
!        tx(i+Order) = xPoints(i)
!       enddo
!       do i = 1,Order
!        t(i+Order+xNumPoints) = xNumPoints
!        tx(i+Order+xNumPoints) = xPoints(xNumPoints)
!       enddo
! 
!       select case (Left)
!        case (0:1)
!         select case (Right)
!          case (0:1)
!           do i = 2,xNumPoints+2*Order-1
!            xBounds(i-1) = t(i)
!           enddo
!          case (2)
!           do i = 2,xNumPoints+2*Order
!            xBounds(i-1) = t(i)
!           enddo
!         end select
!        case (2)
!         select case (Right)
!          case (0:1)
!           do i = 1,xNumPoints+2*Order-1
!            xBounds(i) = t(i)
!           enddo
!          case (2)
!           do i = 1,xNumPoints+2*Order
!            xBounds(i) = t(i)
!           enddo
!         end select
!       end select
! 
!       deallocate(t)
! 
!       u = 0.0d0
! 
!       do k = 1,xNumPoints-1
!        ax = xPoints(k)
!        bx = xPoints(k+1)
!        xIntScale = 0.5d0*(bx-ax)
!        xScaledZero = 0.5d0*(bx+ax)
!        do l = 1,LegPoints
!         x(l,k) = xIntScale*xLeg(l)+xScaledZero
!        enddo
!       enddo
! 
!       Count = 1
!       select case (Left)
!        case (0)
!         do k = xBounds(Count),xBounds(Count+Order+1)-1
!          do l = 1,LegPoints
!           u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,2,tx,b,x(l,k))
!          enddo
!         enddo
!         Count = Count + 1
!        case (1)
!         do k = xBounds(Count),xBounds(Count+Order+1)-1
!          do l = 1,LegPoints
!           u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,1,tx,b,x(l,k))+BSpline(Order,Deriv,xNumPoints,xPoints,2,tx,b,x(l,k))
!          enddo
!         enddo
!         Count = Count + 1
!        case(2)
!         do i = 1,2
!          do k = xBounds(Count),xBounds(Count+Order+1)-1
!           do l = 1,LegPoints
!            u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(l,k))
!           enddo
!          enddo
!          Count = Count + 1
!         enddo
!       end select
! 
!       do i = 3,xNumPoints+Order-3
!        do k = xBounds(Count),xBounds(Count+Order+1)-1
!         do l = 1,LegPoints
!          u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(l,k))
!         enddo
!        enddo
!        Count = Count + 1
!       enddo
! 
!       select case (Right)
!        case (0)
!         do k = xBounds(Count),xBounds(Count+Order+1)-1
!          do l = 1,LegPoints
!           u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,tx,b,x(l,k))
!          enddo
!         enddo
!        case (1)
!         do k = xBounds(Count),xBounds(Count+Order+1)-1
!          do l = 1,LegPoints
!           u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-2,tx,b,x(l,k))+&
! 			&BSpline(Order,Deriv,xNumPoints,xPoints,xNumPoints+Order-1,tx,b,x(l,k))
!          enddo
!         enddo
!        case(2)
!         do i = xNumPoints+Order-2,xNumPoints+Order-1
!          do k = xBounds(Count),xBounds(Count+Order+1)-1
!           do l = 1,LegPoints
!            u(l,k,Count) = BSpline(Order,Deriv,xNumPoints,xPoints,i,tx,b,x(l,k))
!           enddo
!          enddo
!          Count = Count + 1
!         enddo
!       end select
! 
!       deallocate(x,tx,b)
! 
!       return
!       end
!       
!            double precision function bvalue ( t, bcoef, n, k, x, jderiv )
! ! c  from  * a practical guide to splines *  by c. de boor    
! ! calls  interv
! ! c
! ! calculates value at  x  of  jderiv-th derivative of spline from b-repr.
! ! c  the spline is taken to be continuous from the right, EXCEPT at the
! ! c  rightmost knot, where it is taken to be continuous from the left.
! ! c
! ! c******  i n p u t ******
! ! c  t, bcoef, n, k......forms the b-representation of the spline  f  to
! ! c        be evaluated. specifically,
! ! c  t.....knot sequence, of length  n+k, assumed nondecreasing.
! ! c  bcoef.....b-coefficient sequence, of length  n .
! ! c  n.....length of  bcoef  and dimension of spline(k,t),
! ! c        a s s u m e d  positive .
! ! c  k.....order of the spline .
! ! c
! ! c  w a r n i n g . . .   the restriction  k .le. kmax (=20)  is imposed
! ! c        arbitrarily by the dimension statement for  aj, dl, dr  below,
! ! c        but is  n o w h e r e  c h e c k e d  for.
! ! c
! ! c  x.....the point at which to evaluate .
! ! c  jderiv.....integer giving the order of the derivative to be evaluated
! ! c        a s s u m e d  to be zero or positive.
! ! c
! ! c******  o u t p u t  ******
! ! c  bvalue.....the value of the (jderiv)-th derivative of  f  at  x .
! ! c
! ! c******  m e t h o d  ******
! ! c     The nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo-
! ! c  cated with the aid of  interv . The  k  b-coeffs of  f  relevant for
! ! c  this interval are then obtained from  bcoef (or taken to be zero if
! ! c  not explicitly available) and are then differenced  jderiv  times to
! ! c  obtain the b-coeffs of  (d**jderiv)f  relevant for that interval.
! ! c  Precisely, with  j = jderiv, we have from x.(12) of the text that
! ! c
! ! c     (d**j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) )
! ! c
! ! c  where
! ! c                   / bcoef(.),                     ,  j .eq. 0
! ! c                   /
! ! c    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1)
! ! c                   / ----------------------------- ,  j .gt. 0
! ! c                   /    (t(.+k-j) - t(.))/(k-j)
! ! c
! ! c     Then, we use repeatedly the fact that
! ! c
! ! c    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) )
! ! c  with
! ! c                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1)
! ! c    a(.,x)  =    ---------------------------------------
! ! c                 (x - t(.))      + (t(.+m-1) - x)
! ! c
! ! c  to write  (d**j)f(x)  eventually as a linear combination of b-splines
! ! c  of order  1 , and the coefficient for  b(i,1,t)(x)  must then be the
! ! c  desired number  (d**j)f(x). (see x.(17)-(19) of text).
! ! c
!       integer jderiv,k,n,   i,ilo,imk,j,jc,jcmin,jcmax,jj,kmax,kmj,km1&
!       &                     ,mflag,nmi,jdrvp1
!       parameter (kmax = 20)
! ! C     double precision bcoef(n),t(1),x,   aj(20),dl(20),dr(20),fkmj
! !       double precision bcoef(*),t(*),x,   aj(kmax),dl(kmax),dr(kmax),fkmj
! ! c      dimension t(n+k)
! ! c former fortran standard made it impossible to specify the length of  t
! ! c  precisely without the introduction of otherwise superfluous addition-
! ! c  al arguments.
!       bvalue = 0.0d0
!       if (jderiv .ge. k)                go to 99
! ! c
! ! c  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
! ! c      t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
! ! c      outside the support of  the spline  f , hence  bvalue = 0.
! ! c      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
! ! c      at  t(n+k) where it is leftcontinuous.)
!       call interv ( t, n+k, x, i, mflag )
!       if (mflag .ne. 0)                 go to 99
! ! c  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i).
!       km1 = k - 1
!       if (km1 .gt. 0)                   go to 1
!       bvalue = bcoef(i)
!                                         go to 99
! ! c
! ! c  *** store the k b-spline coefficients relevant for the knot interval
! ! c     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
! ! c     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
! ! c     from input to zero. set any t.s not obtainable equal to t(1) or
! ! c     to t(n+k) appropriately.
!     1 jcmin = 1
!       imk = i - k
!       if (imk .ge. 0)                   go to 8
!       jcmin = 1 - imk
!       do j=1,i
!          dl(j) = x - t(i+1-j)
!       enddo
!       do j=i,km1
!          aj(k-j) = 0.0d0
!          dl(j) = dl(i)
!       enddo
!                                         go to 10
!     8 do j=1,km1
!          dl(j) = x - t(i+1-j)
!       enddo
! ! c
!    10 jcmax = k
!       nmi = n - i
!       if (nmi .ge. 0)                   go to 18
!       jcmax = k + nmi
!       do j=1,jcmax
!          dr(j) = t(i+j) - x
!       enddo
!       do j=jcmax,km1
!          aj(j+1) = 0.0d0
!          dr(j) = dr(jcmax)
!       enddo
!                                         go to 20
!    18 do j=1,km1
!          dr(j) = t(i+j) - x
!       enddo
! ! c
!    20 do jc=jcmin,jcmax
!          aj(jc) = bcoef(imk + jc)
!       enddo
! ! c
! ! c               *** difference the coefficients  jderiv  times.
!       if (jderiv .eq. 0)                go to 30
!       do j=1,jderiv
!          kmj = k-j
!          fkmj = float(kmj)
!          ilo = kmj
!          do jj=1,kmj
!             aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
!             ilo = ilo - 1
!          enddo
!       enddo
! ! c
! ! c  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
! ! c     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
!    30 if (jderiv .eq. km1)              go to 39
!       jdrvp1 = jderiv + 1     
!       do j=jdrvp1,km1
!          kmj = k-j
!          ilo = kmj
!          do jj=1,kmj
!             aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
!             ilo = ilo - 1
!          enddo
!       enddo
!    39 bvalue = aj(1)
! ! c
!    99                                   return
!       end
!        subroutine interv ( xt, lxt, x, left, mflag )
! ! c  from  * a practical guide to splines *  by C. de Boor    
! ! computes  left = max( i :  xt(i) .lt. xt(lxt) .and.  xt(i) .le. x )  .
! ! c
! ! c******  i n p u t  ******
! ! c  xt.....a double precision sequence, of length  lxt , assumed to be nondecreasing
! ! c  lxt.....number of terms in the sequence  xt .
! ! c  x.....the point whose location with respect to the sequence  xt  is
! ! c        to be determined.
! ! c
! ! c******  o u t p u t  ******
! ! c  left, mflag.....both integers, whose value is
! ! c
! ! c   1     -1      if               x .lt.  xt(1)
! ! c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
! ! c   i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
! ! c   i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
! ! c
! ! c        In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
! ! c        indicates that  x  lies outside the CLOSED interval
! ! c        xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
! ! c        intervals is due to the decision to make all pp functions cont-
! ! c        inuous from the right, but, by returning  mflag = 0  even if
! ! C        x = xt(lxt), there is the option of having the computed pp function
! ! c        continuous from the left at  xt(lxt) .
! ! c
! ! c******  m e t h o d  ******
! ! c  The program is designed to be efficient in the common situation that
! ! c  it is called repeatedly, with  x  taken from an increasing or decrea-
! ! c  sing sequence. This will happen, e.g., when a pp function is to be
! ! c  graphed. The first guess for  left  is therefore taken to be the val-
! ! c  ue returned at the previous call and stored in the  l o c a l  varia-
! ! c  ble  ilo . A first check ascertains that  ilo .lt. lxt (this is nec-
! ! c  essary since the present call may have nothing to do with the previ-
! ! c  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
! ! c  ilo  and are done after just three comparisons.
! ! c     Otherwise, we repeatedly double the difference  istep = ihi - ilo
! ! c  while also moving  ilo  and  ihi  in the direction of  x , until
! ! c                      xt(ilo) .le. x .lt. xt(ihi) ,
! ! c  after which we use bisection to get, in addition, ilo+1 = ihi .
! ! c  left = ilo  is then returned.
! ! c
!       integer left,lxt,mflag,   ihi,ilo,istep,middle
!       double precision x,xt(lxt)
!       data ilo /1/
!       save ilo  
!       ihi = ilo + 1
!       if (ihi .lt. lxt)                 go to 20
!          if (x .ge. xt(lxt))            go to 110
!          if (lxt .le. 1)                go to 90
!          ilo = lxt - 1
!          ihi = lxt
! ! c
!    20 if (x .ge. xt(ihi))               go to 40
!       if (x .ge. xt(ilo))               go to 100
! ! c
! ! c              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
!       istep = 1
!    31    ihi = ilo
!          ilo = ihi - istep
!          if (ilo .le. 1)                go to 35
!          if (x .ge. xt(ilo))            go to 50
!          istep = istep*2
!                                         go to 31
!    35 ilo = 1
!       if (x .lt. xt(1))                 go to 90
!                                         go to 50
! ! c              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
!    40 istep = 1
!    41    ilo = ihi
!          ihi = ilo + istep
!          if (ihi .ge. lxt)              go to 45
!          if (x .lt. xt(ihi))            go to 50
!          istep = istep*2
!                                         go to 41
!    45 if (x .ge. xt(lxt))               go to 110
!       ihi = lxt
! ! c
! ! c           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
!    50 middle = (ilo + ihi)/2
!       if (middle .eq. ilo)              go to 100
! ! c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
!       if (x .lt. xt(middle))            go to 53
!          ilo = middle
!                                         go to 50
!    53    ihi = middle
!                                         go to 50
! ! c**** set output and return.
!    90 mflag = -1
!       left = 1
!                                         return
!   100 mflag = 0
!       left = ilo
!                                         return
!   110 mflag = 1
! 	  if (x .eq. xt(lxt)) mflag = 0
!       left = lxt
!   111 if (left .eq. 1)                  return
! 	  left = left - 1
! 	  if (xt(left) .lt. xt(lxt))        return
! 										go to 111
!       end
! 
!       
        SUBROUTINE GAULEG(X1,X2,X,W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X1,X2,X(N),W(N),pi
      PARAMETER (EPS=3.D-14)
      pi=dacos(-1.d0)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=COS(pi*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
      
      subroutine main(A,n1,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !parameter (n=60)
      real*8 A(n1,n1),C(n1,n1),indx(n1),g!(n1,n1)
      do i=1,n1
      do j=1,n1
	C(i,j)=0.d0
	enddo
	C(i,i)=1.d0
      enddo

      call ludcmp(A,n1,n1,indx,g)

      do j=1,n1
      call lubksb(A,n1,n1,indx,C(1,j))
      enddo
      return
      end


      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  integer:: II,I,LL,J
	  real*8 A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUMa=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUMa=SUMa-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUMa.NE.0.d0) THEN
          II=I
        ENDIF
        B(I)=SUMa
12    CONTINUE
      DO 14 I=N,1,-1
        SUMa=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUMa=SUMa-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUMa/A(I,I)
14    CONTINUE
      RETURN
      END



      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  PARAMETER (NMAX=1000,TINY=1.d-20)
	  integer:: I,J,K
      real*8 A(NP,NP),INDX(N),VV(NMAX)
      D=1.d0
      DO 12 I=1,N
        AAMAX=0.d0
        DO 11 J=1,N
          IF (dABS(A(I,J)).GT.AAMAX) AAMAX=dABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.d0) then
         print*, 'singular matrix'
         stop
        endif 
!PAUSE 'Singular matrix.'
        VV(I)=1.d0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUMa=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUMa=SUMa-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUMa
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.d0
        DO 16 I=J,N
          SUMa=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUMa=SUMa-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUMa
          ENDIF
          DUM=VV(I)*dABS(SUMa)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.d0)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.d0)A(N,N)=TINY
      RETURN
      END

 SUBROUTINE AIRYA(X,AI,BI,AD,BD)
! C
! C       ======================================================
! C       Purpose: Compute Airy functions and their derivatives
! C       Input:   x  --- Argument of Airy function
! C       Output:  AI --- Ai(x)
! C                BI --- Bi(x)
! C                AD --- Ai'(x)
! C                BD --- Bi'(x)
! C       Routine called:
! C                AJYIK for computing Jv(x), Yv(x), Iv(x) and
! C                Kv(x) with v=1/3 and 2/3
! C       ======================================================
! C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PIR=0.318309886183891D0
        C1=0.355028053887817D0
        C2=0.258819403792807D0
        SR3=1.732050807568877D0
        Z=XA**1.5/1.5D0
        XQ=DSQRT(XA)
        CALL AJYIK(Z,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2)
        IF (X.EQ.0.0D0) THEN
           AI=C1
           BI=SR3*C1
           AD=-C2
           BD=SR3*C2
        ELSE IF (X.GT.0.0D0) THEN
           AI=PIR*XQ/SR3*VK1
           BI=XQ*(PIR*VK1+2.0D0/SR3*VI1)
           AD=-XA/SR3*PIR*VK2
           BD=XA*(PIR*VK2+2.0D0/SR3*VI2)
        ELSE
           AI=0.5D0*XQ*(VJ1-VY1/SR3)
           BI=-0.5D0*XQ*(VJ1/SR3+VY1)
           AD=0.5D0*XA*(VJ2+VY2/SR3)
           BD=0.5D0*XA*(VJ2/SR3-VY2)
        ENDIF
        RETURN
        END


        SUBROUTINE AJYIK(X,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2)
! C
! C       =======================================================
! C       Purpose: Compute Bessel functions Jv(x) and Yv(x),
! C                and modified Bessel functions Iv(x) and
! C                Kv(x), and their derivatives with v=1/3,2/3
! C       Input :  x --- Argument of Jv(x),Yv(x),Iv(x) and
! C                      Kv(x) ( x \F2 0 )
! C       Output:  VJ1 --- J1/3(x)
! C                VJ2 --- J2/3(x)
! C                VY1 --- Y1/3(x)
! C                VY2 --- Y2/3(x)
! C                VI1 --- I1/3(x)
! C                VI2 --- I2/3(x)
! C                VK1 --- K1/3(x)
! C                VK2 --- K2/3(x)
! C       =======================================================
! C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0D0) THEN
           VJ1=0.0D0
           VJ2=0.0D0
           VY1=-1.0D+300
           VY2=1.0D+300
           VI1=0.0D0
           VI2=0.0D0
           VK1=-1.0D+300
           VK2=-1.0D+300
           RETURN
        ENDIF
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        GP1=.892979511569249D0
        GP2=.902745292950934D0
        GN1=1.3541179394264D0
        GN2=2.678938534707747D0
        VV0=0.444444444444444D0
        UU0=1.1547005383793D0
        X2=X*X
        K0=12
        IF (X.GE.35.0) K0=10
        IF (X.GE.50.0) K0=8
        IF (X.LE.12.0) THEN
           DO 25 L=1,2
              VL=L/3.0D0
              VJL=1.0D0
              R=1.0D0
              DO 15 K=1,40
                 R=-0.25D0*R*X2/(K*(K+VL))
                 VJL=VJL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 20
15            CONTINUE
20            A0=(0.5D0*X)**VL
              IF (L.EQ.1) VJ1=A0/GP1*VJL
              IF (L.EQ.2) VJ2=A0/GP2*VJL
25         CONTINUE
        ELSE
           DO 40 L=1,2
              VV=VV0*L*L
              PX=1.0D0
              RP=1.0D0
              DO 30 K=1,K0
                 RP=-0.78125D-2*RP*(VV-(4.0*K-3.0)**2.0)*(VV-&
     &              (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*X2)
30               PX=PX+RP
              QX=1.0D0
              RQ=1.0D0
              DO 35 K=1,K0
                 RQ=-0.78125D-2*RQ*(VV-(4.0*K-1.0)**2.0)*(VV-&
     &              (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*X2)
35               QX=QX+RQ
              QX=0.125D0*(VV-1.0)*QX/X
              XK=X-(0.5D0*L/3.0D0+0.25D0)*PI
              A0=DSQRT(RP2/X)
              CK=DCOS(XK)
              SK=DSIN(XK)
              IF (L.EQ.1) THEN
                 VJ1=A0*(PX*CK-QX*SK)
                 VY1=A0*(PX*SK+QX*CK)
              ELSE IF (L.EQ.2) THEN
                 VJ2=A0*(PX*CK-QX*SK)
                 VY2=A0*(PX*SK+QX*CK)
              ENDIF
40         CONTINUE
        ENDIF
        IF (X.LE.12.0D0) THEN
           DO 55 L=1,2
              VL=L/3.0D0
              VJL=1.0D0
              R=1.0D0
              DO 45 K=1,40
                 R=-0.25D0*R*X2/(K*(K-VL))
                 VJL=VJL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 50
45            CONTINUE
50            B0=(2.0D0/X)**VL
              IF (L.EQ.1) UJ1=B0*VJL/GN1
              IF (L.EQ.2) UJ2=B0*VJL/GN2
55         CONTINUE
           PV1=PI/3.0D0
           PV2=PI/1.5D0
           VY1=UU0*(VJ1*DCOS(PV1)-UJ1)
           VY2=UU0*(VJ2*DCOS(PV2)-UJ2)
        ENDIF
        IF (X.LE.18.0) THEN
           DO 70 L=1,2
              VL=L/3.0D0
              VIL=1.0D0
              R=1.0D0
              DO 60 K=1,40
                 R=0.25D0*R*X2/(K*(K+VL))
                 VIL=VIL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 65
60            CONTINUE
65            A0=(0.5D0*X)**VL
              IF (L.EQ.1) VI1=A0/GP1*VIL
              IF (L.EQ.2) VI2=A0/GP2*VIL
70         CONTINUE
        ELSE
           C0=DEXP(X)/DSQRT(2.0D0*PI*X)
           DO 80 L=1,2
              VV=VV0*L*L
              VSL=1.0D0
              R=1.0D0
              DO 75 K=1,K0
                 R=-0.125D0*R*(VV-(2.0D0*K-1.0D0)**2.0)/(K*X)
75               VSL=VSL+R
              IF (L.EQ.1) VI1=C0*VSL
              IF (L.EQ.2) VI2=C0*VSL
80         CONTINUE
        ENDIF
        IF (X.LE.9.0D0) THEN
           DO 95 L=1,2
              VL=L/3.0D0
               IF (L.EQ.1) GN=GN1
               IF (L.EQ.2) GN=GN2
               A0=(2.0D0/X)**VL/GN
               SUM=1.0D0
               R=1.0D0
               DO 85 K=1,60
                  R=0.25D0*R*X2/(K*(K-VL))
                  SUM=SUM+R
                  IF (DABS(R).LT.1.0D-15) GO TO 90
85             CONTINUE
90            IF (L.EQ.1) VK1=0.5D0*UU0*PI*(SUM*A0-VI1)
              IF (L.EQ.2) VK2=0.5D0*UU0*PI*(SUM*A0-VI2)
95         CONTINUE
        ELSE
           C0=DEXP(-X)*DSQRT(0.5D0*PI/X)
           DO 105 L=1,2
              VV=VV0*L*L
              SUM=1.0D0
              R=1.0D0
              DO 100 K=1,K0
                 R=0.125D0*R*(VV-(2.0*K-1.0)**2.0)/(K*X)
100              SUM=SUM+R
              IF (L.EQ.1) VK1=C0*SUM
              IF (L.EQ.2) VK2=C0*SUM
105        CONTINUE
        ENDIF
        RETURN
        END

	SUBROUTINE JYV(V,X,VM,BJ,DJ,BY,DY)
! C
! C       =======================================================
! C       Purpose: Compute Bessel functions Jv(x) and Yv(x)
! C                and their derivatives
! C       Input :  x --- Argument of Jv(x) and Yv(x)
! C                v --- Order of Jv(x) and Yv(x)
! C                      ( v = n+v0, 0 ף v0 < 1, n = 0,1,2,... )
! C       Output:  BJ(n) --- Jn+v0(x)
! C                DJ(n) --- Jn+v0'(x)
! C                BY(n) --- Yn+v0(x)
! C                DY(n) --- Yn+v0'(x)
! C                VM --- Highest order computed
! C       Routines called:
! C            (1) GAMMA for computing gamma function
! C            (2) MSTA1 and MSTA2 for computing the starting
! C                point for backward recurrence
! C       =======================================================
! C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(0:*),DJ(0:*),BY(0:*),DY(0:*)
        EL=.5772156649015329D0
        PI=dacos(-1.d0)		!3.141592653589793D0
        RP2=.63661977236758D0
        X2=X*X
        N=INT(V)
        V0=V-N
        IF (X.LT.1.0D-100) THEN
           DO 10 K=0,N
              BJ(K)=0.0D0
              DJ(K)=0.0D0
              BY(K)=-1.0D+300
10            DY(K)=1.0D+300
           IF (V0.EQ.0.0) THEN
              BJ(0)=1.0D0
              DJ(1)=0.5D0
           ELSE
              DJ(0)=1.0D+300
           ENDIF
           VM=V  
           RETURN
        ENDIF
        IF (X.LE.12.0) THEN
           DO 25 L=0,1
              VL=V0+L
              BJVL=1.0D0
              R=1.0D0
              DO 15 K=1,40
                 R=-0.25D0*R*X2/(K*(K+VL))
                 BJVL=BJVL+R
                 IF (DABS(R).LT.DABS(BJVL)*1.0D-15) GO TO 20
15            CONTINUE
20            VG=1.0D0+VL
              CALL GAMMA(VG,GA)
              A=(0.5D0*X)**VL/GA
              IF (L.EQ.0) BJV0=BJVL*A
              IF (L.EQ.1) BJV1=BJVL*A
25         CONTINUE
        ELSE
           K0=11
           IF (X.GE.35.0) K0=10
           IF (X.GE.50.0) K0=8
           DO 40 J=0,1
              VV=4.0D0*(J+V0)*(J+V0)
              PX=1.0D0
              RP=1.0D0
              DO 30 K=1,K0
                 RP=-0.78125D-2*RP*(VV-(4.0*K-3.0)**2.0)*(VV-&
     &              (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*X2)
30               PX=PX+RP
              QX=1.0D0
              RQ=1.0D0
              DO 35 K=1,K0
                 RQ=-0.78125D-2*RQ*(VV-(4.0*K-1.0)**2.0)*(VV-&
     &              (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*X2)
35               QX=QX+RQ
              QX=0.125D0*(VV-1.0)*QX/X
              XK=X-(0.5D0*(J+V0)+0.25D0)*PI
              A0=DSQRT(RP2/X)
              CK=DCOS(XK)
              SK=DSIN(XK)
              IF (J.EQ.0) THEN
                 BJV0=A0*(PX*CK-QX*SK)
                 BYV0=A0*(PX*SK+QX*CK)
              ELSE IF (J.EQ.1) THEN
                 BJV1=A0*(PX*CK-QX*SK)
                 BYV1=A0*(PX*SK+QX*CK)
              ENDIF
40         CONTINUE
        ENDIF
        BJ(0)=BJV0
        BJ(1)=BJV1
        DJ(0)=V0/X*BJ(0)-BJ(1)
        DJ(1)=-(1.0D0+V0)/X*BJ(1)+BJ(0)
        IF (N.GE.2.AND.N.LE.INT(0.9*X)) THEN
           F0=BJV0
           F1=BJV1
           DO 45 K=2,N
              F=2.0D0*(K+V0-1.0D0)/X*F1-F0
              BJ(K)=F
              F0=F1
45            F1=F
        ELSE IF (N.GE.2) THEN
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              N=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F2=0.0D0
           F1=1.0D-100
           DO 50 K=M,0,-1
              F=2.0D0*(V0+K+1.0D0)/X*F1-F2
              IF (K.LE.N) BJ(K)=F
              F2=F1
50            F1=F
           IF (DABS(BJV0).GT.DABS(BJV1)) THEN
               CS=BJV0/F
           ELSE
               CS=BJV1/F2
           ENDIF
           DO 55 K=0,N
55            BJ(K)=CS*BJ(K)
        ENDIF
        DO 60 K=2,N
60         DJ(K)=-(K+V0)/X*BJ(K)+BJ(K-1)
        IF (X.LE.12.0D0) THEN
           IF (V0.NE.0.0) THEN
              DO 75 L=0,1
                 VL=V0+L
                 BJVL=1.0D0
                 R=1.0D0
                 DO 65 K=1,40
                    R=-0.25D0*R*X2/(K*(K-VL))
                    BJVL=BJVL+R
                    IF (DABS(R).LT.DABS(BJVL)*1.0D-15) GO TO 70
65               CONTINUE
70               VG=1.0D0-VL
                 CALL GAMMA(VG,GB)
                 B=(2.0D0/X)**VL/GB
                 IF (L.EQ.0) BJU0=BJVL*B
                 IF (L.EQ.1) BJU1=BJVL*B
75            CONTINUE
              PV0=PI*V0
              PV1=PI*(1.0D0+V0)
              BYV0=(BJV0*DCOS(PV0)-BJU0)/DSIN(PV0)
              BYV1=(BJV1*DCOS(PV1)-BJU1)/DSIN(PV1)
           ELSE
              EC=DLOG(X/2.0D0)+EL
              CS0=0.0D0
              W0=0.0D0
              R0=1.0D0
              DO 80 K=1,30
                 W0=W0+1.0D0/K
                 R0=-0.25D0*R0/(K*K)*X2
80               CS0=CS0+R0*W0
              BYV0=RP2*(EC*BJV0-CS0)
              CS1=1.0D0
              W1=0.0D0
              R1=1.0D0
              DO 85 K=1,30
                 W1=W1+1.0D0/K
                 R1=-0.25D0*R1/(K*(K+1))*X2
85               CS1=CS1+R1*(2.0D0*W1+1.0D0/(K+1.0D0))
              BYV1=RP2*(EC*BJV1-1.0D0/X-0.25D0*X*CS1)
           ENDIF
        ENDIF
        BY(0)=BYV0
        BY(1)=BYV1
        DO 90 K=2,N
           BYVK=2.0D0*(V0+K-1.0D0)/X*BYV1-BYV0
           BY(K)=BYVK
           BYV0=BYV1
90         BYV1=BYVK
        DY(0)=V0/X*BY(0)-BY(1)
        DO 95 K=1,N
95         DY(K)=-(K+V0)/X*BY(K)+BY(K-1)
        VM=N+V0
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
! C
! C       ==================================================
! C       Purpose: Compute gamma function ג(x)
! C       Input :  x  --- Argument of ג(x)
! C                       ( x is not equal to 0,-1,-2,תתת)
! C       Output:  GA --- ג(x)
! C       ==================================================
! C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,&
     &          -0.6558780715202538D0, -0.420026350340952D-1,&
     &          0.1665386113822915D0,-.421977345555443D-1,&
     &          -.96219715278770D-2, .72189432466630D-2,&
     &          -.11651675918591D-2, -.2152416741149D-3,&
     &          .1280502823882D-3, -.201348547807D-4,&
     &          -.12504934821D-5, .11330272320D-5,&
     &          -.2056338417D-6, .61160950D-8,&
     &          .50020075D-8, -.11812746D-8,&
     &          .1043427D-9, .77823D-11,&
     &          -.36968D-11, .51D-12,&
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)
! C
! C       ===================================================
! C       Purpose: Determine the starting point for backward  
! C                recurrence such that the magnitude of    
! C                Jn(x) at that point is about 10^(-MP)
! C       Input :  x     --- Argument of Jn(x)
! C                MP    --- Value of magnitude
! C       Output:  MSTA1 --- Starting point   
! C       ===================================================
! C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)
! C
! C       ===================================================
! C       Purpose: Determine the starting point for backward
! C                recurrence such that all Jn(x) has MP
! C                significant digits
! C       Input :  x  --- Argument of Jn(x)
! C                n  --- Order of Jn(x)
! C                MP --- Significant digit
! C       Output:  MSTA2 --- Starting point
! C       ===================================================
! C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
