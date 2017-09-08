program main

implicit double precision (a-h,o-z)

integer*4, allocatable ::		lobound(:),upbound(:)
integer*4 				i,j,k,left,right,NumBeta,lpts,order,np,matdim,ii,ij,jj
real*8 					rinit,r0,energy,field,beta,l,momentum,energ,ml,Rampl,Sampl,cosgamma,singamma
real*8 					lomom,phib,pi,flast,glast,regeta,irregeta,irregscat,logderf(2),logderg(2)
real*8					c6,c8,c10,b6,sc8,sc10,fbc,gbc,fpbc,gpbc,r0x,v,x,vm,dervflast,dervglast,ttanphi,djm,bjm,sc6
real*8					phi,ttanphi0, kmom, spherj,sphern,derspherj,dersphern,energyHF,wro1,wro2,taneta,gfunc,dergfunc
real*8					kclosed,tangamma
double precision, allocatable ::        aqdt(:),gqdt(:),eta(:)
double precision, allocatable :: 	ub(:,:),pts(:),xleg(:),wleg(:), weight(:),surfaceL(:), surface(:),derfsol(:),dergsol(:)
double precision, allocatable :: 	sectors(:),acoef(:),bcoef(:),dat(:,:),uf(:),ug(:),grid(:),fsol(:),gsol(:),duf(:),dug(:)
double precision, allocatable ::	derbs(:,:),kesav(:,:),pesav(:,:),oversav(:,:),kesav2(:,:),pesav2(:,:),oversav2(:,:)
double precision, allocatable ::    ksrmat2(:,:)
Real*8					BJ(0:250),DJ(0:250),BY(0:250),DY(0:250)
double precision:: energylist(400),energylisttemp(400),energylistall(5,400)
double precision:: largef , dervlargef,mu,bfield
double precision:: mfab
double precision:: estepdd
integer:: intmfab
integer:: iindex(10),jindex(10),indextest(60)
double precision:: alphasav(100,6),alphasym2bsav(100,6)
integer*4, allocatable ::ipv(:)
double precision, allocatable::work(:),work2(:)
double precision, allocatable ::ke(:,:),pe0(:,:),pe1(:,:), overlap(:,:)
double precision:: detlist(10), etriallist(10)
integer:: fdsteparray(5)
double precision:: fdcoeffarray(5)
double precision,allocatable:: Ksrmat(:,:),Fmatrix(:,:),logdervF(:),Fpmatrix(:,:),invImatrix(:,:),&
& Jmatrix(:,:),evalksr(:),Rmatrix(:,:),gama(:)
double precision,allocatable::fvector(:),fpvector(:),gvector(:),gpvector(:)
real*8 Thresholds(400), VQcal(0:1,400,400),thresholds2(50),vqcal2(0:1,50,50)
double precision:: boundvec(400)
!double precision:: kmomlist(0:1200), tanetalist(0:700),etalist(0:700),gqdtlist(0:700),aqdtlist(0:700),gamalist(0:1200),tangamalist(0:1200)
double precision:: kmomlist1(796),kmomlist2(1651),etalist(796,3),glist(796,3),alist(796,3),gamalist(1651,3)
double precision::muk,MHz
complex(8):: asarray(400),phasearray(400),trlnsphys(5)
complex(8),allocatable::evalsb(:),evals(:),evalsf(:),vl(:,:),vr(:,:),sphyssav(:,:),evalsb2(:),evalsf2(:),evalsall(:,:)
double precision,allocatable::invkqq(:,:)
double precision, allocatable:: dktilde(:,:),dkinv(:,:),dk(:,:)
double precision, allocatable:: dinvmed(:,:),crosssection(:,:)
double precision, allocatable:: evalqphys(:), RWORK(:)
complex(8), allocatable:: cwork(:)
complex(8), allocatable::sphys(:,:),sphysf2(:,:),sphysb2(:,:),cinvmed(:,:),cmed(:,:),qphys(:,:),sphysb(:,:),sphysf(:,:),dsphys(:,:),dsphyspre(:,:)
complex(8), allocatable:: sphysdev(:,:,:)
complex(8)::cj
complex(8)::trlnsphysb,trlnsphysf,dtrlnsde,dtrlnsdepre
double precision:: bx,wL
double precision:: evectorssav(400,400),evectorssav1(400,400)
double precision,allocatable:: ardeparray(:,:),aideparray(:,:)
integer:: index(20)
double precision:: alpha1(100,6),alphasym2b(100,6)
double precision:: Pmatrix(400,400)
double precision:: timedelay(5001)
double precision:: wltable(11)
double precision, allocatable::evectemp(:)
integer,allocatable::iorder(:)
double precision::blarge(20),bsmall(20)
integer:: iblarge(20), ibsmall(20)
integer:: numhypfarray(-4:4)
double precision,allocatable::bdepanalyze(:,:)
integer:: nblock
common /nblock/ nblock

cj=(0d0,1d0)
pi=dacos(-1.d0)
call setup
!  Initializing Single and Triplet Potentials      
   VSinglet = vxssgp_unret(0.d0,-0.5430d-4,4700.d0)
   VTriplet = vatsup_unret(0.d0,-1.7505d-4,4700.d0)
call cpu_time(totaltime1)
!  plot singlet and triplet potentials
open(14,file='runtime.dat')
open(100,file='potrb.dat')
   do ir = 1,1001
      write(100,*) 1.0d-1*ir,vxssgp_unret(1.0d-1*ir,-0.5430d-4,4700.d0),vatsup_unret(1.0d-1*ir,-1.7505d-4,4700.d0)
   enddo  

read(*,*)
read(*,*) mfab,ithbare,iflag,iflagBorwL,iflagmol
read(*,*)

read(*,*)
read(*,*) Einput
read(*,*)

read(*,*)
read(*,*) sc6,sc8,sc10
read(*,*)

read(*,*)
read(*,*) dm1,dm2
read(*,*)

read(*,*)
read(*,*) r0
read(*,*)

read(*,*)
read(*,*)Binitial,Bfinal,numpoints
read(*,*)

read(*,*)
read(*,*)wLinitial,wLfinal,numwLpoints
read(*,*)

read(*,*)
read(*,*) bx,ww,nblock

if (nblock>11) then
write(*,*) 'Too many floquet blocks'
stop
end if


!==============================================================================
!physical constants
!==============================================================================


Gauss = 4.254382547308656d-10
ceau2uk=3.15771d11
muK = 3.166815355d-12
MHz = 1.519829846006321d-10


!==============================================================================
!constants of the system Rb87+Rb87
!==============================================================================

mu=dm1*dm2/(dm1+dm2)!158425.7361670824d0/2d0

b6=(2d0*mu*sc6)**(0.25d0)
c6= 1.d0
c8=(sc8/sc6)/b6**2.d0
c10=(sc10/sc6)/b6**4.d0

ebeta=1d0/(2d0*mu*b6**2)
write(*,*)"iflag=",iflag
if (iflag==1) then
write(*,*) "Exact CC calculation"
else
write(*,*) "Frame Transformation"
end if

write(*,*)"b6=",b6,"ebeta=",ebeta

!------------------------------------------------------------------------------


!================================================================================
!NUMERICAL PARAMETERS FOR THE RMATRIX AND
! GAULEG AND MAKEBASIS SUBROUTINES
! numerical and physical adjustable paramters
!================================================================================
left=0									!PARAMTERS THAT DEFINE THE NUMBER OF REACTION ZONES, IE LEFT=RIGHT=2
right=2									!MEANS TWO REACTION ZONES ONE ON THE L.H.S. AND THE OTHER ON R.H.S.
order=8									!THE ORDER OF THE B-SPLINES
np=2000									!NUMBER OF BASIS FUNCTIONS
lpts=64									!NUMBER OF LEGENDRE POINTS
matdim=np+order+left/2+right/2-3						!THE OVERALL MATRIX DIMENSIONALITY NEEDED IN RMATRIX SUBROUNTINE
!========================================

rinit=1.d-1									!THE GRID POINT ON THE L.H.S. OF THE INTERPARTICLE DISTANCE
!r0=35.d0/143.9d0									!THE GRID POINT ON THE R.H.S. OF THE INTERPARTICLE DISTANCE
					
l=0.d0									!THE ORBITAL ANGULAR MOMENTUM

rinitp=1d-4

!write(*,*)Binitial,Bfinal,numpoints
energy=Einput/ebeta
r0=r0/b6
Bstep=(Bfinal-Binitial)/dfloat(numpoints-1)
wLstep=(wLfinal-wLinitial)/dfloat(numwLpoints-1)
bx=bx*Gauss
wL=wL*MHz
write(*,*) Einput,energy
write(*,*) c6,c8,c10
Bfield=Binitial*gauss
!call rfchannel(intMFab,Bfield,bx,wL,NumStates,Thresholds,VQcal,evectorssav)

!================================================================================================
!------------------------------------------------------------------------------
!Import long range QDT PARAMTERS
!------------------------------------------------------------------------------
!================================================================================================

open(16,file="alist.dat",action="read")
open(17,file="glist.dat",action="read")
open(18,file="etalist.dat",action="read")
open(19,file="gamalist.dat",action="read")
etalist=0d0
glist=0d0
alist=0d0
gamalist=0d0

do i=2,796
read(16,*) Alist(i,1),Alist(i,2),Alist(i,3)
read(17,*) glist(i,1),glist(i,2),glist(i,3)
read(18,*) etalist(i,1),etalist(i,2),etalist(i,3)
end do

do i=2,1651
read(19,*) gamalist(i,1),gamalist(i,2),gamalist(i,3)
enddo



!===============================================
!MAKE B-SPLINE BASIS FUNCTIONS AND LEGENDRE MESH
!===============================================
!allocate matrices in constructing Bspline basis
allocate(lobound(matdim))
allocate(upbound(matdim))
allocate(ub((order+1)*lpts,matdim*2),pts((np-1)*lpts),xleg(lpts),wleg(lpts), weight((np-1)*lpts))
allocate(surface(2*matdim),surfaceL(2*matdim),derbs((order+1)*lpts,matdim))
Allocate(sectors(matdim))
allocate(uf(matdim),ug(matdim),duf(matdim),dug(matdim),grid(lpts*matdim),fsol(lpts*matdim),gsol(lpts*matdim)&
&,derfsol(lpts*matdim),dergsol(lpts*matdim))


allocate(ke(matdim,matdim),pe0(matdim,matdim),pe1(matdim,matdim),overlap(matdim,matdim))

call GAULEG(-1.d0,1.d0,xleg,wleg,lpts)
!------------------------------------------------------------------------------------------------------------------
isav=0
!go to 305
!------------------------------------------------------------------------------------
!call makebasis to make B-spline basis for F calculation
!------------------------------------------------------------------------------------
call cpu_time(tmake1)
call makebasis(ub,lobound,upbound,matdim,order,left,right,np,r0,rinitp,lpts,xleg,&
&wleg,weight,pts,surfaceL,surface,sectors,6)
call cpu_time(tmake2)
write(*,*) "time in makebasisi=",tmake2-tmake1

write(*,*) "matdim=",matdim

!------------------------------
!CALCULATE THE SURFACE TERMS
do j=1,matdim
If (j.eq.1) then
surfaceL(1)=0.d0
surface(1)=0.d0
elseif (j.eq.matdim) then
surfaceL(matdim)=0.d0
surface(matdim)=1.d0
else
surfaceL(j)=0.d0
surface(j)=0.d0

endif
enddo
!-----------------------------


call cpu_time(tmake1)
call basicmatrix(rinitp,r0,ebeta,energy,l,c6,c8,c10,b6,left,right,&
&lpts,np,order,matdim,ub,lobound,upbound,pts,weight,&
& ke,pe0,pe1,overlap)
call cpu_time(tmake2)

write(*,*) "time in make basic matrix=",tmake2-tmake1

!initiate the B-spline basis of calculating fprime and gprime
allocate(fvector(1),fpvector(1),gvector(1),gpvector(1))
energylist=0d0
call fgreference(energylist, l, 1,rinit,r0,c6,c8,c10,fvector,fpvector,gvector,gpvector)
deallocate (fvector,fpvector,gvector,gpvector)

305 continue
!------------------------------------------------------------------------------------------------------------------
! calculate asymptotic channels and thresholds
!------------------------------------------------------------------------------------------------------------------
open(42,file='Bdep.dat')
!open(124,file='ksrmat.dat')
open(45,file='phaseshift.dat')
open(23,file='boundenergy.dat')
open(678,file='boundvector.dat')
open(125,file='thresholds.dat')
open(43, file='crosssection.dat')
open(44,file='evector.dat')
open(46,file='sphys.dat')
open(896,file='trlnsphys.dat')
open(897,file='evalqphys.dat')
open(898,file='reqphys.dat')
open(899,file='imqphys.dat')
open(791,file='timedelay.dat')
open(39,file='channelinfo.dat')
open(40,file='molinfo.dat')
open(675,file='ardeparray.dat')
open(676,file='aideparray.dat')
open(677,file='resinfo.dat')
call setup

!------------------------------------------------------------------------------------------------------------------
! Start to iterate ( either scan B field or rf frequency wL ) 
!------------------------------------------------------------------------------------------------------------------

if (iflagBorwL==0) then
numwLpoints=1
else if (iflagBorwL==1) then
numpoints=1
end if
numhypfarray=(/1,2,5,6,8,6,5,2,1/)
intmfab=int(mfab)
ilower=0
do i=-4,intmfab-1
ilower=ilower+numhypfarray(i)
end do
if (mod(nblock,2)==1) then
istart=36*(nblock-1)/2+ilower+ithbare
else
istart=36*(nblock-2)/2+ilower+ithbare
end if
write(*,*) 'istart=',istart
isav=0

wltable=(/ 1087.2951328d0, 1087.31510707d0, 1087.33575235d0, 1087.35715662d0,1087.37934189d0, 1087.40237416d0, 1087.42634144d0, 1087.45133171d0,1087.47742198d0, 1087.50475525d0, 1087.53349653d0 /)

bxinitial=18d0*gauss
bxfinal=28d0*gauss
numbx=1
do ibx=1,numbx
!bx=sqrt(bxinitial**2+(bxfinal**2-bxinitial**2)*(ibx-1)/dfloat(numbx))
evalqphyspre=-1d3
!wLinitial=wltable(ibx)-0.0005d0!1084.77d0-0.000856559d0*(bx/gauss)**2-0.02d0
!1087.15d0+0.0096114d0+0.0004011d0*(bx/gauss)**2-0.01d0!985.3121478d0 - 0.00007145277d0*(bx/gauss)**2-0.01d0
!wLfinal=wltable(ibx)+0.0005d0!1084.77d0 - 0.000856559d0* (bx/gauss)**2+0.02d0
!1087.15d0+0.0096114d0+0.0004011d0*(bx/gauss)**2+0.1d0!985.3121478d0 - 0.00007145277d0*(bx/gauss)**2+0.01d0
!if (bx/gauss<20d0) then
!numwLpoints=2001
!wLstep=(wLfinal-wLinitial)/dfloat(numwLpoints-1)
!else
!numwLpoints=2001
!wLstep=(wLfinal-wLinitial)/dfloat(numwLpoints-1)
!end if

allocate(ardeparray(numpoints,numwLpoints),aideparray(numpoints,numwLpoints))
!allocate(bdepanalyze(numpoints,2))
do ib=1,numpoints                       !Scan Magneic field 
isav=0
Bfield=(Binitial+(ib-1)*bstep)*gauss
diffB=bstep/10d0*gauss
!write(*,*) 'mfab=',mfab
!call Hyperfine2b(int(mfab),Bfield,NumSym2bStates,Thresholds2,VQcal2,indexsav,alphasav,alphasym2bsav)
!write(*,*) Bfield/gauss, Thresholds2(1:numsym2bstates)
!call rfchannel(int(MFab),Bfield,bx,wL,NumStates,Thresholds,VQcal,evectorssav,indextest)

do iwL=1,numwLpoints                    !Scan photon frequency

wL=(wLinitial+(iwL-1)*wLstep)*MHz


!call rfchannel(Bfield+diffB,bx,wL,nphoton,NumStates,Thresholds,VQcal,evectorssav1,indextest,isav)
!isav=0
call rfchannel(Bfield,bx,wL,nblock,NumStates,Thresholds,VQcal,evectorssav,indextest,isav)
!call Hyperfine2b(intmfab,bfield,NumStates,Thresholds,VQcal,index,alpha1,alphasym2b)
!write(125,'(400e20.10)') wL/(MHz),(thresholds(i),i=1,numstates)
!write(*,*) "numstates=", numstates
!do i=1,numstates
!write(792,'(400e20.10)') (evectorssav(i,j),j=1,numstates)
!end do
!cycle
!write(*,*) 'good', iwL
!write(211,'(400e20.10)') thresholds(50:70)
!write(212,'(400e20.10)') (evectorssav(101,i),i=50,70)
!cycle
!do j=91,95
!write(499,'(400e20.10)') bfield/gauss, dfloat(j),(evectorssav(i,j),i=1,numstates)
!end do 

!write(66,'(400e20.10)')wL/(MHz),(Thresholds(i),i=1,numstates)
!stop
!if (evectorssav(29,1)<0d0) then
!evectorssav(1:numstates,1)=evectorssav(1:numstates,1)*(-1d0)
!end if

!if (evectorssav1(29,1)<0d0) then
!evectorssav1(1:numstates,1)=evectorssav1(1:numstates,1)*(-1d0)
!end if

!if (evectorssav(23,2)<0d0) then
!evectorssav(1:numstates,2)=evectorssav(1:numstates,2)*(-1d0)
!end if

!if (evectorssav1(23,2)<0d0) then
!evectorssav1(1:numstates,2)=evectorssav1(1:numstates,2)*(-1d0)
!end if

!if (evectorssav(23,3)<0d0) then
!evectorssav(1:numstates,3)=evectorssav(1:numstates,3)*(-1d0)
!end if

!if (evectorssav1(23,3)<0d0) then
!evectorssav1(1:numstates,3)=evectorssav1(1:numstates,3)*(-1d0)
!end if

!do i=1,3
!do j=1,3

!total1=0d0
!total2=0d0
!do m=1,numstates
!total1=total1+evectorssav1(m,i)*evectorssav(m,j)
!total2=total2+evectorssav(m,i)*evectorssav1(m,j)
!end do
!
!Pmatrix(i,j)=(total1-total2)/2d0/diffB*gauss

!end do
!end do

!write(399,'(400e20.10)') Bfield/gauss,Pmatrix(1,2),Pmatrix(1,3),Pmatrix(2,3)


456 continue
!write(*,*)indextest

!write(467,*)"ithresh=",ithresh

!Select the state with the largest component of |1111> from the middle floquet
!block as entrance channel the states of |1111>+2hv is 36+29 state before
!diagonalization 
!the |11>+|2-1> should be in 53rd state if FB=3
!the |11>+|11> should be in 36+29 state if FB=3
!if (ib==1 .and. iwL==1) then
dlargestcomp=0d0
do i=1, numstates
if (abs(evectorssav(istart,i))>dlargestcomp) then
dlargestcomp=abs(evectorssav(istart,i))
itag=i
end if
end do
!end if
!write(*,*) 'so good 437'
!itag=53
!write(908,*) wL/(MHz), evectorssav(36+29,itag),itag
!cycle
do ithresh=itag,itag
!write(*,*) 'so good 439'
!================================================================================================

if (iflag==1) then
!-----------------------------------------------------------------------
!CALCULATE KSR MATRIX FROM THE EXACT CC SOLUTION APPROACH
! allocate shorrange matrices depending on energy/Bfield
allocate(Ksrmat(NumStates,NumStates),ksrmat2(NumStates,NumStates),Fmatrix(NumStates,NumStates),&
logdervF(NumStates),evalksr(numstates),fpmatrix(NumStates,NumStates))
allocate(invImatrix(NumStates,NumStates),Jmatrix(NumStates,NumStates),Rmatrix(numstates,numstates))
allocate(fvector(NumStates),fpvector(NumStates),gvector(NumStates),gpvector(NumStates))
allocate(ipv(NumStates),work(5*NumStates))

write(*,*) "numstates=",numstates

do i=1,numstates
energylist(i)=energy+(Thresholds(ithresh)-Thresholds(i))/ebeta
end do
write(*,*) "thresholds="
!write(*,*) thresholds(1:numstates)
write(*,*) "energylist="
!write(*,*) energylist(1:numstates)



!------------------------------------------------------------------------------------------------------------------------

call cpu_time(tmake1)
call rrmatrixall(rinitp,r0,ebeta,energy,ithresh,l,c6,c8,c10,beta,left,right,&
&lpts,np,order,matdim,ub,lobound,upbound,pts,weight,surface,surfaceL,ke,pe0,pe1,overlap,NumStates,&
&Thresholds,VQcal,Fmatrix,logdervF,Fpmatrix)
call cpu_time(tmake2)
write(*,*) "time in rrmatrixall=",tmake2-tmake1
write(*,*)"================================================"
write(*,*) "print logdervF"
write(*,*) logdervF(1:numstates)
write(*,*)"================================================"
!stop
!-------------------------------------indextest(i)-----------------------------------------------------------------------------------
!set up the long range reference function \hat{f} and \hat{g} calculation
!------------------------------------------------------------------------------------------------------------------------

call fgreference(energylist, l, numstates,rinit,r0,c6,c8,c10,fvector,fpvector,gvector,gpvector)
!------------------------------------------------------------------------------------------------------------------------
write(*,*) fvector(numstates),fpvector(numstates),gvector(numstates),gpvector(numstates)
!stop
do i=1,NumStates
wrogf=gvector(i)*fpvector(i)-gpvector(i)*fvector(i)
do jb=1,NumStates
wro1=gvector(i)*Fpmatrix(i,jb)-gpvector(i)*Fmatrix(i,jb)
wro2=fvector(i)*Fpmatrix(i,jb)-fpvector(i)*Fmatrix(i,jb)
invImatrix(i,jb)=wro1/wrogf
Jmatrix(i,jb)=wro2/wrogf
end do
end do


call DGETRF( numstates, numstates, invImatrix,numstates, IPV, INFO )
call DGETRI( numstates, invImatrix, numstates, IPV, WORK, 5*numstates, INFO3 )
ksrmat=matmul(Jmatrix,invImatrix)
ksrmat2=ksrmat
open(55,file='ksrmat.dat')
!write(*,*) 'indextest',indextest
!write(*,*)'numstates=',numstates
!do i=1,numstates
!do j=1,numstates
!mm=indextest(i)
!nn=indextest(j)
!ksrmat2(i,j)=ksrmat(mm,nn)
!end do
!end do
!write(*,*)ksrmat2
!do i=1,numstates
!write(55,'(200e20.10)') (ksrmat2(i,j),j=1,numstates)
!end do
!stop


call DSYEV( 'V', 'U', numstates, ksrmat2, numstates, evalksr, WORK, 5*numstates, INFO )
write(41,'(100e20.10)') datan(evalksr)/pi

deallocate(ipv,work)
deallocate(ksrmat2,Fmatrix,logdervF,evalksr,fpmatrix)
deallocate(invImatrix,Jmatrix,Rmatrix)
deallocate(fvector,fpvector,gvector,gpvector)

!-----------------------------------------------------------------------
!CALCULATE KSR MATRIX FROM THE FRAME TRANSFORMATION APPROACH
else if (iflag==0) then
allocate(Ksrmat(NumStates,NumStates))

do i=1,numstates
energylist(i)=energy+(Thresholds(ithresh)-Thresholds(i))/ebeta
end do

!tanuSqdt= -0.14491954479344132      tanuTqdt= -0.25339662366532090
!uSqdt=  -4.5810401986873771E-002 uTqdt=  -7.8995894446919207E-002

!make the quantum defects field dependent by interpolation using two points
!(9.1G,0.98),(1007G,1) for singlet and (9.1G,1.03),(1007G,1.2) for triplet 
delu0=1.231d0!0.74d0!1.231d0!(1d0-0.98d0)/(1007d0-9.1d0)*(Bfield/gauss-9.1d0)+0.98d0
delu1=0.9905d0!0.74d0!0.9905d0!(1.2d0-1.03d0)/(1007d0-9.1d0)*(Bfield/gauss-9.1d0)+1.03d0
uqdt0=-4.646258483579639d-2*delu0!-4.5810401986873771d-2*delu0!(0.98d0)
uqdt1=-7.883272773919514d-2*delu1!-7.8995894446919207d-2*delu1!(1.03d0)
tanuSqdt= -0.14491954479344132d0
tanuTqdt= -0.25339662366532090d0

do i=1,numstates
do j=1,numstates

ksrmat(i,j)=VQcal(0,i,j)*tan(pi*uqdt0)+VQcal(1,i,j)*tan(pi*uqdt1)

end do
end do

else
write(*,*) 'iflag is wrong'
stop
end if
!================================================================================================

if (iflagmol==0) go to 301
!================================================================================================
!calculate the molecular states by solving det(ksrmat+cot(gama))=0
Do NumOpen = NumStates,1,-1
If(Thresholds(NumOpen)<(energy)*ebeta+thresholds(ithresh)) Exit
Enddo

write(*,*) "numopen=",numopen

NumClose=NumStates-numopen


allocate(ipv(NumClose),work(5*NumClose))
allocate(gama(NumClose))
allocate(invkqq(NumClose,NumClose))
allocate(evalksr(NumClose))
write(*,*) "NumClose=",NumClose


!do i=1,51
!estep=1d-10
!etrial=0d0-estep*dfloat(i-1)
!determinant1=determin(Etrial,ithresh, numstates,ksrmat,Thresholds,ebeta,gamalist,kmomlist)
!write(43,*) etrial/MHz,determinant1
!end do


estepdd=1d-10
itrial=0
detlist=0d0
etriallist=0d0
do i=1,100
!write(*,*) i,estepdd
tol=1d-8
Etrial1=-estepdd*dfloat(i-1)
Etrial2=-estepdd*dfloat(i)
!write(*,*) Etrial1,Etrial2
determinant1=determin(Etrial1, ithresh, numstates,ksrmat,Thresholds,ebeta,gamalist,kmomlist,bfield,wL,evectorssav,0)
determinant2=determin(Etrial2, ithresh, numstates,ksrmat,Thresholds,ebeta,gamalist,kmomlist,bfield,wL,evectorssav,0)
!write(*,*) determinant1,determinant2

if (determinant1*determinant2*1d8<0d0) then
!write(*,*) 'itrial=',itrial
itrial=itrial+1
do numiter=1,200
!determinant1=determin(Etrial1,1,  numstates,ksrmat,Thresholds,ebeta,gamalist,kmomlist)
!determinant2=determin(Etrial2,1,  numstates,ksrmat,Thresholds,ebeta,gamalist,kmomlist)
Etrial3=Etrial1+(Etrial2-Etrial1)/2d0

determinant3=determin(Etrial3, ithresh,numstates,ksrmat,Thresholds,ebeta,gamalist,kmomlist,bfield,wL,evectorssav,0)
if (abs(determinant3)<tol .or. abs(Etrial2-Etrial1)/abs(Etrial2)<1d-10) exit

if (determinant1*determinant3>0d0) then
Etrial1=Etrial3
determinant1=determinant3
else
Etrial2=Etrial3
determinant2=determinant3
end if

if (numiter==200) write(*,*) "maximum interation reached"
end do
!write(*,*) etrial3
determinant=determinant3
!write(*,*)"start to write results"
!write(*,'(3e20.10,I4)') Bfield/Gauss,etrial3/MHz,determinant,numiter
!write(56,'(3e20.10,I4)') Bfield/Gauss,etrial3/MHz,determinant,numiter
detlist(itrial)=determinant
etriallist(itrial)=etrial3
if (abs(Etrial3/MHz)<1d-1) then

dummy=determin(Etrial3, ithresh,numstates,ksrmat,Thresholds,ebeta,gamalist,kmomlist,bfield,wL,evectorssav,1)
end if

end if

end do
write(23,'(100e20.10)')Bfield/gauss, wL/(MHz), (etriallist(i)/MHz,i=itrial,1,-1)
deallocate(gama,invkqq)
deallocate(ipv,work,evalksr)
deallocate(ksrmat)
write(*,*) 'good', iwl
!stop
cycle
!go to 1006
!================================================================================================




!================================================================================================
!Matching the short range Kmatrix Ksr to long range QDT parameters to get physical observables
301 continue
!deallocate(ipv,work)
Do NumOpen = NumStates,1,-1
If(Thresholds(NumOpen)<(energy)*ebeta+thresholds(ithresh)) Exit
Enddo

!write(*,*) "numopen=",numopen
NumClose=NumStates-numopen

write(*,*) 'ithresh=',ithresh,'numopen=',numopen,'numclose=',numclose

allocate(evectemp(numstates),iorder(numstates))
if (ib==1 .and. iwL==1) then
do j=1,ithresh
evectemp=evectorssav(1:numstates,j)
call sortorder(evectemp,iorder,numstates)
write(39,'(10I15)') j, iorder(1:5)
write(39,'(I15,10e15.5)') j, evectemp(1:5)**2  
enddo
endif
deallocate(evectemp,iorder)


if (ib==2) then
do j=1,ithresh
!write(44,'(400e20.10)') Bfield/gauss,dfloat(j),(evectorssav(i,j),i=1,numstates)
end do 
write(44,*)"--------------------------------------------------------"
end if
!write(*,*) "SO FAR OK"
allocate(gama(numclose))
allocate(invkqq(NumClose,NumClose))
allocate(eta(numopen),aqdt(numopen),gqdt(numopen))
allocate(dktilde(numopen,numopen),dkinv(numopen,numopen),dk(numopen,numopen),sphys(numopen,numopen),sphysb(numopen,numopen),sphysf(numopen,numopen))
!allocate(sphysf2(numopen,numopen),sphysb2(numopen,numopen))
allocate(sphysdev(5,numopen,numopen))
allocate(evalsall(5,numopen),evalsb(numopen),evals(numopen),evalsf(numopen),vr(numopen,numopen),vl(numopen,numopen),evalsf2(numopen),evalsb2(numopen))
allocate(dinvmed(numopen,numopen),cinvmed(numopen,numopen),cmed(numopen,numopen))
allocate(cwork(5*numopen))
allocate(qphys(numopen,numopen))
allocate(dsphys(numopen,numopen),dsphyspre(numopen,numopen))
allocate(work2(5*numopen),evalqphys(numopen),RWORK(3*numopen-2))
allocate(crosssection(numopen,numopen))
!go to 112
sphysdev=(0d0,0d0)
energylisttemp=energylist
fdsteparray=(/ -2,-1,0,1,2 /)
fdcoeffarray=(/ 1d0/12d0,-2d0/3d0,0d0,2d0/3d0,-1d0/12d0 /)
ediv=1d3
!write(*,*) 'so far so good'
if (energylist(numopen)-energy/ediv<0d0 .or. energylist(numopen+1)+energy/ediv>0d0) ediv=10d0*energy/abs(energylist(numopen+1))
dsphys=(0d0,0d0)
dtrlnsde=(0d0,0d0)
isph=1
do while (isph<=20)
sphysdev=(0d0,0d0)
isphflag=0
!dsphyspre=dsphys
do iter=1,5
do i=1,400
energylistall(iter,i)=energylisttemp(i)+dfloat(fdsteparray(iter))*energy/ediv
end do 
end do
dtrlnsdepre=dtrlnsde

do iter=1,5
allocate(ipv(NumClose),work(5*NumClose))
energylist=energylistall(iter,:)
!write(*,'(400e22.12)') energylist(1:numopen)

gama=0d0
eta=0d0
aqdt=0d0
gqdt=0d0

do j=1,numopen
eta(j)=csplint(sqrt(energylist(j)),etalist,796)
aqdt(j)=csplint(sqrt(energylist(j)),alist,796)
gqdt(j)=csplint(sqrt(energylist(j)),glist,796)
end do

!write(*,*) "eta=",eta,aqdt,gqdt
do j=numopen+1,numstates
gama(j-numopen)=csplint(sqrt(abs(energylist(j))),gamalist,1651)
end do

!write(*,*) gamalist(1:5)
!write(*,*) "tangama=",tan(gama)
!write(*,*) "eta=",eta
!write(*,*) "a=",aqdt
!write(*,*) "g=",gqdt

invkqq=ksrmat(numopen+1:numstates,numopen+1:numstates)
do j=1,numclose
invkqq(j,j)=invkqq(j,j)+1d0/tan(gama(j))
end do
!write(*,*) "numclose=",numclose,size(invkqq,1),size(invkqq,2)
!write(*,*) "invkqq",invkqq

call DGETRF( numclose, numclose, invkqq, numclose, IPV, INFO )

if (abs(info)>1d-8) then
write(*,*) "info =" , info
stop
end if
!write(*,*)'so far so good 1' 
call DGETRI( numclose, invkqq, numclose, IPV, WORK, 5*numclose, INFO3)
if (abs(info3)>1d-8) then
write(*,*) "info3 =" , info3
stop
end if
!write(*,*)"good invert invkqq"

dkinv=0d0
do ni=1,numopen
do nj=1,numopen

do i=1,numclose
do j=1,numclose

dkinv(ni,nj)=dkinv(ni,nj)+ksrmat(ni,i+numopen)*invkqq(i,j)*ksrmat(j+numopen,nj)

enddo
enddo

enddo
enddo


dktilde=ksrmat(1:numopen,1:numopen)-dkinv


dinvmed=0d0
do i=1,numopen
do j=1,numopen
dinvmed(i,j)=gqdt(i)*dktilde(i,j)
end do
end do

do i=1,numopen
dinvmed(i,i)=dinvmed(i,i)+1d0
end do


deallocate(ipv,work)
allocate(ipv(numopen),work(5*numopen))
call DGETRF( numopen, numopen, dinvmed, numopen, IPV, INFO )

if (abs(info)>1d-8) then
write(*,*) "info =" , info
stop
end if

call DGETRI( numopen, dinvmed, numopen, IPV, WORK, 5*numopen, INFO3)
if (abs(info3)>1d-8) then
write(*,*) "info3 =" , info3
stop
end if
!write(*,*) "invert dinvmed"
dk=0d0

do i=1,numopen
do j=1,numopen

do m=1,numopen
dk(i,j)=dk(i,j)+sqrt(aqdt(i))*dktilde(i,m)*dinvmed(m,j)*sqrt(aqdt(j))
end do

end do
end do
!write(*,*)'so far so good 2' 
cmed=cj*dk
cinvmed=-cj*dk
do i=1,numopen
cinvmed(i,i)=cinvmed(i,i)+(1d0,0d0)
cmed(i,i)=cmed(i,i)+(1d0,0d0)
end do
call zgetrf(numopen,numopen,cinvmed,numopen,ipv,info)
if (abs(info)>1d-8) then
write(*,*) "info =" , info
stop
end if

call zGETRI( numopen, cinvmed, numopen, IPV, cWORK, 5*numopen, INFO3)
if (abs(info3)>1d-8) then
write(*,*) "info3 =" , info3
stop
end if
!write(*,*) "invert cinvmed"

!sphysdev(iter,:,:)=(0d0,0d0)
do i=1,numopen
do j=1,numopen

do m=1,numopen
sphysdev(iter,i,j)=sphysdev(iter,i,j)+exp(cj*eta(i))*cmed(i,m)*cinvmed(m,j)*exp(cj*eta(j))
end do

end do 
end do

deallocate(ipv,work)
end do

!dsphys=(sphysf-sphysb)/(energy/ediv*2d0)
do iter=1,5
sphysb=sphysdev(iter,:,:)
call   ZGEEV( 'N', 'N', numopen, sphysb, numopen, evalsb, VL, numopen, VR, numopen, cWORK, 5*numopen, RWORK, INFO )
if (abs(info)>1d-8) then
write(*,*) 'info=',info
stop
end if
evalsall(iter,:)=evalsb
end do 

trlnsphys=(0d0,0d0)
do i=1,numopen
if (abs(abs(evalsb(i))-1d0)>1d-6) then 
write(*,*) 'evals is not unitary', i
stop
end if
end do 
do iter=1,5
do i=1,numopen
trlnsphys(iter)=trlnsphys(iter)+zlog(evalsall(iter,i))
end do 
end do
dtrlnsde=(0d0,0d0) 
do iter=1,5
dtrlnsde=dtrlnsde+(trlnsphys(iter)*fdcoeffarray(iter))/(energy/ediv)
enddo
!write(*,*) trlnsphys(1:4),dtrlnsde!trlnsphysf,dtrlnsde
!stop

if (abs(dtrlnsde-dtrlnsdepre)/abs(dtrlnsde)>1d-6) then
isph=isph+1
ediv=ediv*2d0
cycle
else
exit
end if

end do
!This is the end of iter part of iteration 
!---------------------------------------------------------------------------------

!write(*,*) dtrlnsde


!write(*,*) 'good,isphs=', isph
!stop

sphys=sphysdev(3,:,:)

go to 908

write(*,*) 'numclose=',numclose, 'numopen=',numopen

allocate(ipv(NumClose),work(5*NumClose))

energylist=energylisttemp

gama=0d0
eta=0d0
aqdt=0d0
gqdt=0d0

do j=1,numopen
eta(j)=csplint(sqrt(energylist(j)),etalist,796)
aqdt(j)=csplint(sqrt(energylist(j)),alist,796)
gqdt(j)=csplint(sqrt(energylist(j)),glist,796)
end do

!write(*,*) "eta=",eta,aqdt,gqdt
do j=numopen+1,numstates
gama(j-numopen)=csplint(sqrt(abs(energylist(j))),gamalist,1651)
end do



!write(*,*) gamalist(1:5)
!write(*,*) "tangama=",tan(gama)
!write(*,*) "eta=",eta
!write(*,*) "a=",aqdt
!write(*,*) "g=",gqdt

invkqq=ksrmat(numopen+1:numstates,numopen+1:numstates)
do j=1,numclose
invkqq(j,j)=invkqq(j,j)+1d0/tan(gama(j))
end do
!write(*,*) "numclose=",numclose,size(invkqq,1),size(invkqq,2)
!write(*,*) "invkqq",invkqq

call DGETRF( numclose, numclose, invkqq, numclose, IPV, INFO )

if (abs(info)>1d-8) then
write(*,*) "info =" , info
stop
end if

call DGETRI( numclose, invkqq, numclose, IPV, WORK, 5*numclose, INFO3)
if (abs(info3)>1d-8) then
write(*,*) "info3 =" , info3
stop
end if
!write(*,*)"good invert invkqq"
!write(*,*) numopen,numclose
dkinv=0d0
do ni=1,numopen
do nj=1,numopen
dummy1=0d0
do i=1,numclose
do j=1,numclose
dummy1=dummy1+ksrmat(ni,i+numopen)*invkqq(i,j)*ksrmat(j+numopen,nj)
!write(*,*) dummy1
!dkinv(ni,nj)=dkinv(ni,nj)+ksrmat(ni,i+numopen)*invkqq(i,j)*ksrmat(j+numopen,nj)

enddo
enddo
dkinv(ni,nj)=dummy1
enddo
enddo
!write(*,*) ((dkinv(i,j),i=1,27),j=1,27)
!write(*,*)'good in dkinv'
!stop

dktilde=ksrmat(1:numopen,1:numopen)-dkinv


dinvmed=0d0
do i=1,numopen
do j=1,numopen
dinvmed(i,j)=gqdt(i)*dktilde(i,j)
end do
end do

do i=1,numopen
dinvmed(i,i)=dinvmed(i,i)+1d0
end do


deallocate(ipv,work)
allocate(ipv(numopen),work(5*numopen))
call DGETRF( numopen, numopen, dinvmed, numopen, IPV, INFO )

if (abs(info)>1d-8) then
write(*,*) "info =" , info
stop
end if

call DGETRI( numopen, dinvmed, numopen, IPV, WORK, 5*numopen, INFO3)
if (abs(info3)>1d-8) then
write(*,*) "info3 =" , info3
stop
end if
!write(*,*) "invert dinvmed"
dk=0d0

do i=1,numopen
do j=1,numopen

do m=1,numopen
dk(i,j)=dk(i,j)+sqrt(aqdt(i))*dktilde(i,m)*dinvmed(m,j)*sqrt(aqdt(j))
end do

end do
end do

cmed=cj*dk
cinvmed=-cj*dk
do i=1,numopen
cinvmed(i,i)=cinvmed(i,i)+(1d0,0d0)
cmed(i,i)=cmed(i,i)+(1d0,0d0)
end do
call zgetrf(numopen,numopen,cinvmed,numopen,ipv,info)
if (abs(info)>1d-8) then
write(*,*) "info =" , info
stop
end if

call zGETRI( numopen, cinvmed, numopen, IPV, cWORK, 5*numopen, INFO3)
if (abs(info3)>1d-8) then
write(*,*) "info3 =" , info3
stop
end if
!write(*,*) "invert cinvmed"

sphys=(0d0,0d0)
do i=1,numopen
do j=1,numopen

do m=1,numopen
sphys(i,j)=sphys(i,j)+exp(cj*eta(i))*cmed(i,m)*cinvmed(m,j)*exp(cj*eta(j))
end do

end do
end do

deallocate(ipv,work)


!write(*,*) 'sphys is good'


908 continue







!---------------------------------------------------------------------------------
!construct the time delay matrix
!---------------------------------------------------------------------------------
!dsphys=(sphysf-sphysb)/(energy/100d0*2d0)
!qphys=cj*matmul(sphys,conjg(transpose(dsphys)))

!diagonalize time delay matrix Qphys
!call ZHEEV( 'V', 'U', numopen,qphys, numopen,evalqphys, cWORK, 5*numopen, RWORK,INFO )
!if (abs(info)>1d-8) then
!write(*,*) 'info=',info
!stop
!end if
write(897,'(400e22.12)') wL/(MHz),real(dtrlnsde/(0d0,1d0)), aimag(dtrlnsde/(0d0,1d0))
!write(896,'(400e22.12)') wL/(MHz),real(trlnsphys/(0d0,1d0)),aimag(trlnsphys/(0d0,1d0))
!write(897,'(400e20.10)') evalqphys(1:numopen)
!do i=1,numopen
!write(898,'(400e20.10)') (real(qphys(i,j)),j=1,numopen)
!end do
!do i=1,numopen
!write(899,'(400e20.10)') (aimag(qphys(i,j)),j=1,numopen)
!end do
!write(898,*)'---------------------------'
!write(899,*)'---------------------------'
!write(*,*) "sphys=",sphys

!ascurrent=-real(tan(zlog(sphys(ithresh,ithresh))/(0.d0,2.d0))/sqrt(2d0*mu*Einput))
!write(*,*) "a=",ascurrent
do i=1,ithresh
asarray(i)=-tan(zlog(sphys(i,i))/(0.d0,2.d0))/sqrt(2d0*mu*(Einput+thresholds(ithresh)-thresholds(i)))
phasearray(i)=zlog(sphys(i,i))/(0.d0,2.d0)
end do 

do i=1,numopen
do j=1,numopen
if (i==j) then
crosssection(i,j)=pi/(2d0*mu*(Einput+thresholds(ithresh)-thresholds(i)))*abs(sphys(i,j)-1d0)**2
else 
crosssection(i,j)=pi/(2d0*mu*(Einput+thresholds(ithresh)-thresholds(i)))*abs(sphys(i,j))**2
end if
end do 
end do
 
!write(42,'(400e20.10)') wL/(MHz),Bfield/(4.254382547308656d-10),(real(asarray(i)),i=1,numopen),(aimag(asarray(i)),i=1,numopen),dfloat(ithresh)
!write(42,'(10e20.10)') wL/(MHz),Bfield/(4.254382547308656d-10),-real(tan(zlog(sphys(ithresh,ithresh))/(0.d0,2.d0))/sqrt(2d0*mu*Einput)),&
!-aimag(tan(zlog(sphys(ithresh,ithresh))/(0.d0,2.d0))/sqrt(2d0*mu*Einput)),dfloat(ithresh)
!if (ibx==1 .or. ibx==10) then

!ardeparray(Bfield,frequency)

ardeparray(ib,iwL)=real(asarray(ithresh))
aideparray(ib,iwL)=aimag(asarray(ithresh))
!bdepanalyze(ib,1)=Bfield/gauss
!bdepanalyze(ib,2)=real(asarray(ithresh))
bdepcurrent=real(asarray(ithresh))
if (iwL>2) then
if(bdepcurrent< bdeppre1 .and. bdeppre1>bdeppre2) then
awllarge=bdeppre1
iwLlarge=iwL-1 
write(1234,*) wL/MHz, bdepcurrent,iwL-1
end if
if(bdepcurrent> bdeppre1 .and. bdeppre1<bdeppre2) then 
awlsmall=bdeppre1
iwLsmall=iwL-1
write(1234,*) wL/MHz, bdepcurrent,iwL-1
end if
end if 
bdeppre2=bdeppre1
bdeppre1=bdepcurrent

write(42,'(400e25.15)') bx/gauss,wL/(MHz),Bfield/gauss,real(asarray(ithresh)),aimag(asarray(ithresh)),dfloat(ithresh),evectorssav(istart,ithresh),dfloat(isph)!,evectorssav(101,ithresh), dfloat(isph)!,real(asarray(ithresh)-aspre)/wLstep
!write(45,'(400e22.12)') bx/gauss,wL/(MHz), Bfield/(4.254382547308656d-10), real(phasearray(1:ithresh)) , aimag(phasearray(1:ithresh))
!write(46,'(400e22.12)') bx/gauss,wL/(MHz), Bfield/(4.254382547308656d-10), real(sphys(ithresh,ithresh)) , aimag(sphys(ithresh,ithresh))
!end if
do i=ithresh,ithresh
!write(43,'(400e20.10)') (crosssection(i,j),j=1,numopen)
end do 

!if (evalqphys(numopen)<evalqphyspre .and. evalqphys(numopen)>1d0) then
!write(791,'(400e22.12)') bx/gauss, (bx/gauss)**2,wL/(MHz),evalqphyspre2,evalqphyspre,evalqphys(numopen),real(asarray(ithresh)),aimag(asarray(ithresh))
!igoflag=1
!end if

timedelay(iwL)=real(dtrlnsde/(0d0,1d0))

evalqphyspre2=evalqphyspre
evalqphyspre=evalqphys(numopen)
aspre=asarray(ithresh)

deallocate(invkqq)
deallocate(gama)
deallocate(ksrmat)
deallocate(eta,aqdt,gqdt)
deallocate(dktilde,dkinv,dk,sphys)
deallocate(dinvmed,cinvmed,cmed)
deallocate(cwork)
deallocate(qphys)
deallocate(sphysb,sphysf,dsphys,dsphyspre)
deallocate(sphysdev,evalsall)
deallocate(evals,evalsb,evalsf,vl,vr,evalsb2,evalsf2)
deallocate(work2,evalqphys,RWORK)
deallocate(crosssection)


!if (igoflag==1) go to 107
!if (ascurrent<-1d0 .and. aspre>1d0) then
!if ( abs(496.6d0-Bfield/(4.254382547308656d-10))<1.5d0) then
!write(15,*) jj, Bfield/(4.254382547308656d-10)
!deallocate(Ksrmat,energylist)
!deallocate(gama)
!exit
!end if
enddo     !end iteration of ithresh
enddo     !end iteration of iwL
enddo     !end iteration of ib

do i=1,numpoints
write(675,'(1001e20.10)') (ardeparray(i,j),j=1,numwLpoints)
write(676,'(1001e20.10)') (aideparray(i,j),j=1,numwLpoints)
end do
go to 67
!=========================================================================================================
! find resonance information. First find resonance position by selecting the extremum,
! that is, find the point where (a(i)>a(i+1) .and. a(i)>a(i-1)) .or. (a(i)<a(i+1) .and. a(i)<a(i-1))
!=========================================================================================================
!ardeparray(Bfield,frequency)
blarge=0d0
bsmall=0d0
iblarge=0
ibsmall=0
istep1=0
istep2=0

abg=ardeparray(1,1)
do iwL=1,numwLpoints
do i=2,numpoints-1

if  (ardeparray(i,iwL)>ardeparray(i+1,iwL) .and. ardeparray(i,iwL)>ardeparray(i-1,iwL)) then
istep1=istep1+1
blarge(istep1)=ardeparray(i,iwL)
iblarge(istep1)=i
end if

if (ardeparray(i,iwL)<ardeparray(i+1,iwL) .and. ardeparray(i,iwL)<ardeparray(i-1,iwL)) then
istep2=istep2+1
bsmall(istep2)=ardeparray(i,iwL)
ibsmall(istep2)=i
endif

end do   !end iteraction of finding res on B, the number of res is istep1==step2
write(*,*) 'istep1',istep1,'istep2',istep2
numres=istep1


w0pos=wL/MHz
do j=1,numres

b0pos=Binitial+(dfloat(iblarge(j)+ibsmall(j))/2d0-1d0)*bstep
tunability=blarge(j)-bsmall(j)
gamab=abs(dfloat(iblarge(j)-ibsmall(j)))*bstep
deltab=tunability/abg*(gamab/2d0)

!write(677,'(100e20.10)') bx/gauss,w0pos,b0pos,deltab,gamab,tunability
enddo     ! end of writing resinfo for one wL

enddo   ! end iteraction of scaning wL for finding resinfo
!=========================================================================================================
67 continue
!call resanalyzer(bdepanalyze,numpoints)
w0pos=wLinitial+(dfloat(iwLlarge+iwLsmall)/2d0-1d0)*wLstep
tunability=awLlarge-awLsmall
gamaw=abs(dfloat(iwLlarge-iwLsmall))*wLstep
deltaw=tunability/abg*(gamaw/2d0)
write(677,'(10e22.12)') (bx/gauss)**2, w0pos, tunability, gamaw
107 continue
igoflag=0
deallocate(ardeparray,aideparray)
!write(791,'(5000e20.10)') bx/gauss, (timedelay(i), i=1,numwLpoints)
!deallocate(bdepanalyze)
end do    !end iterationof ibx
!=======================================
deallocate(uf,ug,duf,dug,grid,fsol,gsol&
&,derfsol,dergsol)
deallocate(lobound)
deallocate(upbound)
deallocate(ub,pts,xleg,wleg, weight)
deallocate(surface,surfaceL,derbs)
deallocate(sectors)

call cpu_time(totaltime2)
write(14,*) 'total runing time=',totaltime2-totaltime1

end program

!================================================================================================


Subroutine basicmatrix(rinit,r0,ebeta,energy,l,c6,c8,c10,beta,left,right,&
&lpts,np,order,matdim,ub,lobound,upbound,pts,weight,&
& ke,pe0,pe1,overlap)
implicit double precision (a-h,o-z)

integer*4 			lobound(matdim),upbound(matdim)
integer*4 			matdim,order,left,right,np,lpts,i,rj,ri,iimax,iimin,jp,ir,j,ip,k,info,jl,count0,count1,kl,ku,jo

real*8 				ub((order+1)*lpts,matdim*2),surface(matdim),surfaceL(matdim),sectors(matdim),t1,t2
real*8 				pts((np-1)*lpts),weight((np-1)*lpts),uf(matdim),ug(matdim),sc8,sc10,dug(matdim),duf(matdim)
Real*8				rad,u1,u2,du1,du2,pot,sum1,sum0,sum2,sum3,det,w11,w12,w22,theta,uL1,uL2,uR1,uR2,gammall,dervglast
real*8 				beta,rinit,r0,energy,l

double precision:: ke(matdim,matdim),pe0(matdim,matdim),pe1(matdim,matdim),overlap(matdim,matdim)
double precision, external::vatsup_unret,vxssgp_unret
pi=dacos(-1.d0)

do j=1,matdim									!CALCULATION OF THE POTENTIAL AND KINETIC TERM OF THE HAMILTONIAN
do jp=1,matdim								!ON A B-SPLINE BASIS EMPLOYING THE VARIATIONAL PRINCIPAL
OVERLAP(j,jp)=0.d0							!OVERLAP CORRESPONDS TO THE OVERLAP OF THE WAVEFUNCTION
KE(j,jp)=0.d0
PE0(j,jp)=0.d0
PE1(j,jp)=0.d0
sum0=0.d0
sum1=0.d0
sum2=0.d0
sum3=0d0
sum4=0d0

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

pot1=(l*(l+1.d0))/rad**2.d0+vatsup_unret(rad*beta,-1.7505d-4*1.055d0,4700.d0)/ebeta

pot0=(l*(l+1.d0))/rad**2.d0+vxssgp_unret(rad*beta,-0.5430d-4*1.4d0,4700.d0)/ebeta

sum0=sum0+weight(ir)*u1*u2
sum1=sum1+weight(ir)*u1*u2*(pot0)
sum3=sum3+weight(ir)*u1*u2*(pot1)
sum2=sum2+weight(ir)*du1*du2
enddo
count1=1+count1
KE(j,jp)=sum2
PE0(j,jp)=sum1
PE1(j,jp)=sum3
OVERLAP(j,jp)=sum0
enddo
enddo

end subroutine

!================================================================================================

subroutine rrmatrixall(rinitp,r0,ebeta,energy,ithresh,l,c6,c8,c10,beta,left,right,&
&lpts,np,order,matdim,ub,lobound,upbound,pts,weight,surface,surfaceL,ke,pe0,pe1,overlap,NumStates,&
&Thresholds,VQcal,Fmatrix,logdervF,Fpmatrix)


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
real*8              Thresholds(400), VQcal(0:1,400,400)
double precision, allocatable:: gama(:,:),lambda(:,:),gamaall(:,:),lambdaall(:,:)
double precision, allocatable:: evalr(:),b(:),work(:),work1(:),work2(:),gcr(:),gcl(:),omega(:,:)
double precision, allocatable::	gammaoc(:,:),gammacc(:,:),gammaco(:,:),invgammacc(:,:),gammatemp(:,:),gammaoo(:,:)
double precision, allocatable:: logder(:,:),cg(:),cf(:),vecmat(:,:),work4(:),work5(:),AB(:,:),bs(:,:)
double precision:: ke(matdim,matdim),pe0(matdim,matdim),pe1(matdim,matdim),overlap(matdim,matdim)
double precision:: Fmatrix(NumStates,NumStates),logdervF(NumStates),Fpmatrix(NumStates,NumStates)


!allocate(gama(matdim,matdim),lambda(matdim,matdim))
!allocate(gamaall(numstates*matdim,numstates*matdim),lambdaall(numstates*matdim,numstates*matdim))

!gamaall=0d0
!lambdaall=0d0

!=============================================
!write(*,*) 'pe0='
!write(*,*) (pe0(i,i),i=1000,1010)
!write(*,*) 'pe1='
!write(*,*) (pe1(i,i),i=1000,1010)
!write(*,*) 'ke='
!write(*,*) (ke(i,i),i=1000,1010)
!write(*,*) 'overlap='
!write(*,*) (overlap(i,i),i=1000,1010)

call cpu_time(tmake1)
write(*,*)'start rrmatrixall calculation'
allocate(gammaoo(numstates,numstates),gammaco(numstates*(matdim-1),numstates),&
gammaoc(numstates,numstates*(matdim-1)))!,invgammacc(numstates*(matdim-1),numstates*(matdim-1)))

!allocate(ipv(numstates*(matdim-1)),work2(5*numstates*(matdim-1)))
allocate(omega(numstates,numstates),gammatemp(numstates,numstates),evalr(numstates),work(5*numstates))

write(*,*)'good allocation'
do ic=1,numstates
do jc=1,numstates

if (ic==jc) then
gammaoo(ic,jc)=-ke(matdim,matdim)-VQcal(0,ic,jc)*pe0(matdim,matdim)-VQcal(1,ic,jc)*pe1(matdim,matdim)&
&+(energy+(thresholds(ithresh)-thresholds(ic))/ebeta)*overlap(matdim,matdim)
gammaoo(ic,jc)=gammaoo(ic,jc)/overlap(matdim,matdim)
else
gammaoo(ic,jc)=-VQcal(0,ic,jc)*pe0(matdim,matdim)-VQcal(1,ic,jc)*pe1(matdim,matdim)
gammaoo(ic,jc)=gammaoo(ic,jc)/overlap(matdim,matdim)
end if

end do
end do

write(*,*)'good in gamaoo'
!=============================================
!write(*,*) (gammaoo(1,i),i=1,numstates)


kl=order*numstates
ku=order*numstates
allocate(work1(numstates*(matdim-1)),work2(5*numstates*(matdim-1)),ab(2*kl+ku+1,numstates*(matdim-1)))
write(*,*) "good0"


call cpu_time(time1)
go to 101
istep=0
jstep=0
do ii=1,matdim-1
do ic=1,numstates

istep=istep+1
jstep=0
do jj=1,matdim-1
do jc=1,numstates
jstep=jstep+1


if (ic==jc) then
invgammacc(istep,jstep)=-ke(ii,jj)-VQcal(0,ic,jc)*pe0(ii,jj)-VQcal(1,ic,jc)*pe1(ii,jj)&
&+(energy+(thresholds(ithresh)-thresholds(ic))/ebeta)*overlap(ii,jj)
invgammacc(istep,jstep)=invgammacc(istep,jstep)/sqrt(overlap(ii,ii)*overlap(jj,jj))
else
invgammacc(istep,jstep)=-VQcal(0,ic,jc)*pe0(ii,jj)-VQcal(1,ic,jc)*pe1(ii,jj)
invgammacc(istep,jstep)=invgammacc(istep,jstep)/sqrt(overlap(ii,ii)*overlap(jj,jj))
end if

end do
end do

end do
end do
write(*,*)'good in invgamacc'
101 continue


write(*,*)'start good in invgamacc'
istep=0
jstep=0
do jj=1,matdim-1
do jc=1,numstates

jstep=jstep+1
istep=0
do ii=1,matdim-1
do ic=1,numstates
istep=istep+1

if (istep>=max(1,jstep-ku).and. istep<=min(numstates*(matdim-1),jstep+kl)) then

if (ic==jc) then
AB(kl+ku+1+istep-jstep,jstep) = (-ke(ii,jj)-VQcal(0,ic,jc)*pe0(ii,jj)-VQcal(1,ic,jc)*pe1(ii,jj)&
&+(energy+(thresholds(ithresh)-thresholds(ic))/ebeta)*overlap(ii,jj))/sqrt(overlap(ii,ii)*overlap(jj,jj))
else

AB(kl+ku+1+istep-jstep,jstep) = (-VQcal(0,ic,jc)*pe0(ii,jj)-VQcal(1,ic,jc)*pe1(ii,jj))/sqrt(overlap(ii,ii)*overlap(jj,jj))
end if

end if

end do
end do

end do
end do

!write(*,*) 'AB='
!write(*,*) AB(1,1:10)
!stop
call cpu_time(time2)
write(*,*) "time in construct gamacc=",time2-time1
!=============================================









call cpu_time(time1)
istep=0
jstep=0
do ii=matdim,matdim
do ic=1,numstates
istep=istep+1
jstep=0
do jj=1,matdim-1
do jc=1,numstates
jstep=jstep+1


if (ic==jc) then
gammaoc(istep,jstep)=-ke(ii,jj)-VQcal(0,ic,jc)*pe0(ii,jj)-VQcal(1,ic,jc)*pe1(ii,jj)&
&+(energy+(thresholds(ithresh)-thresholds(ic))/ebeta)*overlap(ii,jj)
gammaoc(istep,jstep)=gammaoc(istep,jstep)/sqrt(overlap(ii,ii)*overlap(jj,jj))
else
gammaoc(istep,jstep)=-VQcal(0,ic,jc)*pe0(ii,jj)-VQcal(1,ic,jc)*pe1(ii,jj)
gammaoc(istep,jstep)=gammaoc(istep,jstep)/sqrt(overlap(ii,ii)*overlap(jj,jj))
end if

end do
end do

end do
end do
call cpu_time(time2)
write(*,*) "time in construct gamaoc=",time2-time1





istep=0
jstep=0
do ii=1,matdim-1
do ic=1,numstates

istep=istep+1
jstep=0
do jj=matdim,matdim
do jc=1,numstates
jstep=jstep+1


if (ic==jc) then
gammaco(istep,jstep)=-ke(ii,jj)-VQcal(0,ic,jc)*pe0(ii,jj)-VQcal(1,ic,jc)*pe1(ii,jj)&
&+(energy+(thresholds(ithresh)-thresholds(ic))/ebeta)*overlap(ii,jj)
gammaco(istep,jstep)=gammaco(istep,jstep)/sqrt(overlap(ii,ii)*overlap(jj,jj))
else
gammaco(istep,jstep)=-VQcal(0,ic,jc)*pe0(ii,jj)-VQcal(1,ic,jc)*pe1(ii,jj)
gammaco(istep,jstep)=gammaco(istep,jstep)/sqrt(overlap(ii,ii)*overlap(jj,jj))
end if

end do
end do

end do
end do



allocate(bs(numstates*(matdim-1),numstates),ipv(numstates*(matdim-1)))


write(*,*) size(ab),size(ab,1),size(ab,2)

bs=gammaco

call cpu_time(tmake1)
call DGBSV(numstates*(matdim-1),KL,KU,numstates,ab,2*kl+ku+1,ipv,bs,numstates*(matdim-1),INFO3)  !LINEAR SOLVER SUBROUTINE

call cpu_time(tmake2)
write(*,*) "time in construct and inverse matrix=",tmake2-tmake1
!stop

if (abs(info3).gt. 1.d-16) then
print*, 'ROUTINE DGBSV: INFO=/=0',info3
stop
ENDIF
deallocate(ab)

call cpu_time(tmake1)
do j=1,numstates
do jp=1,numstates

sum=0d0
do mm=1,numstates*(matdim-1)

sum=sum+gammaoc(j,mm)*bs(mm,jp)

end do
gammatemp(j,jp)=sum

end do
end do
call cpu_time(tmake2)
write(*,*) "time in construct omega=",tmake2-tmake1


omega=gammaoo-gammatemp

call DSYEV( 'V', 'U', numstates, omega, numstates, evalr, WORK, 5*numstates, INFO )

if (abs(info)>1d-6) then
write(*,*) "info(dsyev)=",info
stop
end if


Fmatrix=omega
logdervF=evalr*overlap(matdim,matdim)
do i=1,numstates
Fpmatrix(:,i)=-Fmatrix(:,i)*logdervF(i)
end do

end subroutine



!------------------------------------------------------------------------------------------------------------------------

subroutine fgreference(energylist,l,numstates,rinit,r0,c6,c8,c10,fvector,fpvector,gvector,gpvector)

implicit double precision (a-h,o-z)
double precision:: l,logderf(2),logderg(2)
integer*4, allocatable ::		lobound(:),upbound(:)
double precision, allocatable :: 	ub(:,:),pts(:),xleg(:),wleg(:), weight(:),surfaceL(:), surface(:),derfsol(:),dergsol(:)
double precision, allocatable :: 	sectors(:),acoef(:),bcoef(:),dat(:,:),uf(:),ug(:),grid(:),fsol(:),gsol(:),duf(:),dug(:)
double precision, allocatable ::	derbs(:,:),kesav(:,:),pesav(:,:),oversav(:,:),kesav2(:,:),pesav2(:,:),oversav2(:,:)


double precision:: fvector(numstates),fpvector(numstates),gvector(numstates),gpvector(numstates),energylist(numstates)
integer::order,left, right

SAVE

pi=dacos(-1.d0)


left=2									!PARAMTERS THAT DEFINE THE NUMBER OF REACTION ZONES, IE LEFT=RIGHT=2
right=2									!MEANS TWO REACTION ZONES ONE ON THE L.H.S. AND THE OTHER ON R.H.S.
order=8									!THE ORDER OF THE B-SPLINES
np=2000									!NUMBER OF BASIS FUNCTIONS
lpts=64
matdim=np+order+left/2+right/2-3						!THE OVERALL MATRIX DIMENSIONALITY NEEDED IN RMATRIX SUBROUNTINE
!========================================

if (abs(energylist(1))<1d-10) then

allocate(lobound(matdim))
allocate(upbound(matdim))
allocate(ub((order+1)*lpts,matdim*2),pts((np-1)*lpts),xleg(lpts),wleg(lpts), weight((np-1)*lpts))
allocate(surface(2*matdim),surfaceL(2*matdim),derbs((order+1)*lpts,matdim))
Allocate(sectors(matdim))
allocate(uf(matdim),ug(matdim),duf(matdim),dug(matdim),grid(lpts*matdim),fsol(lpts*matdim),gsol(lpts*matdim)&
&,derfsol(lpts*matdim),dergsol(lpts*matdim))


allocate(kesav(matdim,matdim),pesav(matdim,matdim),oversav(matdim,matdim))


call GAULEG(-1.d0,1.d0,xleg,wleg,lpts)

isav=0
call cpu_time(tmake1)
call makebasis(ub,lobound,upbound,matdim,order,left,right,np,r0,rinit,lpts,xleg,&
&wleg,weight,pts,surfaceL,surface,sectors,6)
call cpu_time(tmake2)
write(*,*) "time in makebasisi=",tmake2-tmake1

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

ii=32
print*, 'COMPUTE THE ZERO ENERGY PHASE TANPHI'
call phiphase(rinit,0.2d0*dble(ii),l,0.d0,c6,c8,c10,ttanphi)
write(*,*) ttanphi
!============================================
return
end if

do i=1,numstates
energy=energylist(i)


call wkbboundcond(l,rinit,c8,c10,energy,ttanphi,fbc,gbc,fpbc,gpbc)

logderf(1)=fpbc							!THE DERIVATIVE OF THE REGULAR FUNCTION F
logderf(2)=fbc							!THE REGULAR FUNCTION F
logderg(1)=gpbc							!THE DERIVATIVE OF THE IRREGULAR FUNCTION G
logderg(2)=gbc							!THE IRREGULAR FUNCTION G

call rrmatrix(rinit,r0,energy,field,l,c6,c8,c10,beta,logderf,logderg,left,right,&
&lpts,np,order,matdim,ub,derbs,lobound,upbound,pts,weight,surface,surfaceL,&
&Rampl,Sampl,cosgamma,singamma,flast,dervflast,glast,dervglast,uf,ug,duf,dug,kesav,pesav,oversav,isav)

fvector(i)=flast
fpvector(i)=dervflast
gvector(i)=glast
gpvector(i)=dervglast

end do
!write(*,*)fvector
!write(*,*)gvector
!write(*,*)fpvector

call cpu_time(tmake2)
write(*,*) "time in fhat and ghat=",tmake2-tmake1

return
end subroutine


!------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------
! the function determin calculate the determinant det|cos(gama)*Ksr+sin(gama)|
double precision function determin(Etrial,ithresh,numstates,ksrmat,Thresholds,ebeta,gamalist,kmomlist,bfield,wL,evectorssav,ivecflag)

implicit double precision (a-h,o-z)
double precision:: ksrmat(numstates,numstates),Thresholds(400),gamalist(1651,3),kmomlist(0:1200),energylist2(400)
double precision:: boundvec(400),evectorssav(400,400),bvector(numstates)
double precision:: MHz,bfield,gauss
double precision,allocatable::work(:), invkqq(:,:), evalksr(:),gama(:)
double precision,allocatable::evectemp(:)
integer,allocatable::iorder(:)
integer:: istateinfo(5)
save detrel
common /nblock/ nblock

allocate(work(50*numstates))
allocate(gama(numstates-ithresh+1))
allocate(invkqq(numstates-ithresh+1,numstates-ithresh+1))
allocate(evalksr(numstates-ithresh+1))
MHz=1.519829846006321d-10
gauss=4.254382547308656d-10
pi=acos(-1d0)
do i=1,numstates-ithresh+1
energylist2(i)=(Etrial+Thresholds(ithresh)-Thresholds(ithresh+i-1))/ebeta
end do




do j=1,numstates-ithresh+1
gama(j)=csplint(sqrt(abs(energylist2(j))),gamalist,1651)
end do




do i=1,numstates-ithresh+1
do j=1,numstates-ithresh+1
invkqq(i,j)=sin(gama(i))*ksrmat(ithresh+i-1,ithresh+j-1)
end do
end do

do j=1,numstates-ithresh+1
invkqq(j,j)=invkqq(j,j)+cos(gama(j))
end do

call DSYEV( 'V', 'U', numstates-ithresh+1, invkqq, numstates-ithresh+1, evalksr, WORK, 50*numstates, INFO )
if (abs(info)>1d-8) then
write(*,*) "info=/=0", info
stop
end if

if (ivecflag==1) then
!write(*,*) 'evalksr='
!write(*,'(400e20.10)') evalksr
!write(*,*)'ithresh=',ithresh,'numstates=',numstates
evalsmallest=1d5
do i=1,numstates-ithresh+1
if (abs(evalksr(i))<evalsmallest) then 
evalsmallest=abs(evalksr(i))
itag=i
end if
end do
!write(*,*) 'invkqq='
do i=1,numstates-ithresh+1
!write(*,'(400e20.8)') (invkqq(i,j),j=1,numstates-ithresh+1)
end do  
write(*,*) 'itag=',itag,'numstates=',numstates






bvector=0d0
do i=1,numstates
do j=1,numstates-ithresh+1

bvector(i)=bvector(i)+evectorssav(i,j+ithresh-1)*invkqq(j,itag)

end do 
end do 


allocate(evectemp(numstates),iorder(numstates))
evectemp=bvector
call sortorder(evectemp,iorder,numstates)
write(40,'(10e20.8)') bfield/gauss, wL/(MHz), etrial/MHz
write(40,'(10I15)')  iorder(1:5)
do i=1,5
call stateinterpreter(iorder(i),nblock,istateinfo)
write(40,'(10I15)') iorder(i), istateinfo(1:5)
end do
write(40,'(10e15.5)')  evectemp(1:5)**2
write(40,*)'=================================================='
deallocate(evectemp,iorder)




write(678,'(400e20.8)') itag, bfield/gauss,wL/(MHz),etrial/MHz, evalksr(itag),(invkqq(j,itag),j=1,numstates-ithresh+1)
!write(399,'(400e20.10)') bfield/gauss,wL/(MHz),etrial/MHz,evalksr(1:numstates-ithresh+1)

end if

if (abs(etrial)<1d-15) then
determinant=1d0
do i=1,numstates-ithresh+1
determinant=determinant*evalksr(i)
end do
determin=1d0
detrel=determinant
else

determinant=1d0/detrel
do i=1,numstates-ithresh+1
determinant=determinant*evalksr(i)
end do
determin=determinant
end if

deallocate(gama,invkqq)
deallocate(work,evalksr)

end function

subroutine sortorder(b,iordera,n)
implicit double precision (a-h,o-z)
double precision:: a(n),b(n)
integer::iordera(n)
a=b
do i=1,5
amax=0d0
do j=1,n
if (a(j)**2>amax) then
amax=a(j)**2
jtag=j
endif
enddo
b(i)=a(jtag)
a(jtag)=0d0
iordera(i)=jtag 
enddo
end subroutine




double precision function csplint(x,xmesh,n)
implicit double precision (a-h,o-z)
double precision::x
double precision::xmesh(n,3)

! 1st column of xmesh are x points
! 2nd column of xmesh are f(x) points
! 3rd column of xmesh are f''(x) points

do i=1,n-1
if (x>=xmesh(i,1) .and. x<xmesh(i+1,1)) then
aa=(xmesh(i+1,1)-x)/(xmesh(i+1,1)-xmesh(i,1))
bb=1-aa
cc=(aa**3-aa)*(xmesh(i+1,1)-xmesh(i,1))**2/6d0
dd=(bb**3-bb)*(xmesh(i+1,1)-xmesh(i,1))**2/6d0
csplint=aa*xmesh(i,2)+bb*xmesh(i+1,2)+cc*xmesh(i,3)+dd*xmesh(i+1,3)
return
exit
end if
enddo
if (i==n) then
csplint=(xmesh(n,2)-xmesh(n-10,2))/(xmesh(n,1)-xmesh(n-10,1))*(x-xmesh(n,1))+xmesh(n,2)
return
end if
end function





subroutine fppgen(xmesh,n)
implicit double precision (a-h,o-z)
double precision:: xmesh(n,3)
integer:: ipv(n-2)
double precision::bs(n-2,1)
double precision,allocatable::ab(:,:)

kl=1
ku=1
allocate(ab(2*kl+ku+1,n-2))
do i=1,n-2
do j=1,n-2

if (j-i==-1) then
AB(KL+KU+1+i-j,j)=(xmesh(i+1,1)-xmesh(i,1))/6d0
elseif (j==i) then
AB(KL+KU+1+i-j,j)=(xmesh(i+2,1)-xmesh(i,1))/3d0
elseif (j-i==1) then
AB(KL+KU+1+i-j,j)=(xmesh(i+2,1)-xmesh(i+1,1))/6d0
end if

end do
end do

do i=1,n-2
bs(i,1)=(xmesh(i+2,2)-xmesh(i+1,2))/(xmesh(i+2,1)-xmesh(i+1,1))-(xmesh(i+1,2)-xmesh(i,2))/(xmesh(i+1,1)-xmesh(i,1))
enddo

call DGBSV(n-2,KL,KU,1,ab,2*kl+ku+1,ipv,bs,n-2,INFO3)
if (info3/=0) then
write(*,*) 'info3=/=0'
stop
end if

xmesh(1,3)=0d0
xmesh(n,3)=0d0
do i=2,n-1
xmesh(i,3)=bs(i-1,1)
end do

end subroutine



