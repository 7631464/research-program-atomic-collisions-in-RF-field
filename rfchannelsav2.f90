subroutine rfchannel(Bfield,bx,wL,nblock,ndimension,Thdress,Pdress,evectorssav,indextest,isav)
!the inputs Bfield, bx, wL are in atomic units
!this subroutine defines each floquet block as the nearly energy degenerate
!states. 
implicit double precision (a-h,o-z)

double precision::alpha1(100,6),alphasym2b(-4:4,100,6),alphasym2btemp(100,6)
integer NumHypfStates1
double precision i1,s1
double precision Pi,Gauss,Bfield
double precision Delta,Lambda,xi1
double precision gE1,muE1,gN1,muN1
double precision MHz
integer :: NumSym2bStates(-4:4)

integer::mfab

double precision Thresholds(-4:4,50),Thresholdstemp(50)
double precision VQcal(-4:4,0:1,50,50),VQcaltemp(0:1,50,50)
double precision evectorssav(400,400)

double precision::Pdress(0:1,400,400),Thdress(400),hcouple(50,50),hcouple1b(8,8)!,Pdress2(0:1,400,400)
integer:: Index(-4:4,20),Indextemp(20),indextest(60)
integer:: antiIndextemp(20),antiindex(-4:4,20)
double precision, allocatable:: energies2b(:)
double precision, allocatable:: work2(:)
double precision, allocatable:: Hcomplete(:,:), energiesrf(:),evectors(:,:),Hfloquetblk(:,:),Hinterfloquet(:,:)

double precision, allocatable:: Energies1(:), U1(:,:)
double precision, allocatable:: QcalTemp(:,:,:)!,hcpack(:,:,:)
double precision:: hcpack(-4:3,8,8)

double precision clebsh,Sym2b,Qcal,QcalSym2b
external clebsh,Sym2b,Qcal,QcalSym2b
!save hcouple,antiindex,VQcal,Thresholds, n2bstates, energies2b
save thresholds,VQcal,numsym2bstates,ndfloquet,hcpack

Gauss = (4.254382547308656e-10)
GHz = 1.519829846006321d-7
MHz =1.519829846006321d-10
pi=dacos(-1d0)
if (isav==1) go to 999
!nphoton=1
isav=1

call Parameters_87Rb(i1,s1,xi1,gE1,muE1,gN1,muN1,NumHypfStates1)

allocate(Energies1(NumHypfStates1))
allocate(U1(NumHypfStates1,NumHypfStates1))
call CalcHypfStates(Bfield,i1,s1,xi1,muN1,gN1,muE1,gE1,Energies1,U1)
!write(*,*) "NumHypfStates1=",NumHypfStates1
call  CalcAlpha(i1,s1,0,xi1,muN1,gN1,muE1,gE1,NumSym2bStates,alpha1,alphaSym2btemp)
!go to 309

do mfi=-4,4
call Hyperfine2b(mfi,Bfield,NumSym2bStates(mfi),Thresholdstemp,VQcaltemp,Indextemp,alpha1,alphaSym2btemp)
index(mfi,:)=Indextemp
Thresholds(mfi,:)=Thresholdstemp
VQcal(mfi,:,:,:)=VQcaltemp
alphaSym2b(mfi,:,:)=alphaSym2btemp
end do

!do mfi=-4,4
!do i=1,8
!write(793,'(100e22.12)') (VQcal(mfi,0,i,j), j=1,8)
!write(794,'(100e22.12)') (VQcal(mfi,1,i,j), j=1,8)
!end do
!end do 

!write(*,*) "index0(1)=",index0(1)
!write(*,*) "indexm1(1)=",indexm1(1)
!write(*,*) "|",alpha1(7,1:2), ">,|",alpha1(8,1:2),">,|", alpha2(1,1:2),">"

!the array y=index(x), x is 2benergy sorted state label, y is unsorted state label

!antiindex is the inverse function of index


do mfi=-4,4
do i=1,NumSym2bStates(mfi)
antiIndextemp(Index(mfi,i))=i
enddo
antiindex(mfi,:)=antiindextemp
end do

n2bstates=0
do mfi=-4,4
n2bstates=n2bstates+NumSym2bStates(mfi)
end do
write(*,*) "n2bstates=",n2bstates

allocate(Energies2b(n2bstates))
!allocate(hcpack(-4:3,NumSym2bStates(0),NumSym2bStates(0)))

jstep=0
do mfi=-4,4
do j=1,NumSym2bStates(mfi)
jstep=jstep+1
energies2b(jstep)=Energies1(nint(alphaSym2b(mfi,j,1)))+Energies1(nint(alphaSym2b(mfi,j,2)))
enddo
enddo
309 continue
!write(*,*) energies2b(1:n2bstates)

!================================================================================================
!Construct the off-diagonal part of the asymptotic Hamiltonian
!go to 1001


!initialize the rfcouple function
rfdummy=rfcouple(bx,-1,0,alpha1,U1)

!write(*,*) 'a <->b', rfcouple(bx,1,2,alpha1,U1)/MHz
!do i=1,7
!rf2=rfcouple(bx,i,i+1,alpha1,U1)
!write(*,*) 'rf2=',rf2
!end do 
!write(*,*)'bx=',bx/gauss,'U1=',(U1(i,1),i=1,8),'alpha1=',(alpha1(i,2),i=1,8)
!write(*,*)'rfcouple=',rfcouple(bx,7,8,alpha1,U1),rfcouple(bx,2,3,alpha1,U1)
!rf2=Sym2b(7,1,1,7)* &
!Sym2b(8,1,1,8)*&
!rfcouple(bx,7,8,alpha1,U1)
!write(*,*) 'rf2=',rf2
! function rfcouple is the rf field coupling two hyperfine states of one atom,
! I use it to construct two body coupling term and symemtrize them

do i=1,8
do j=1,8
hcouple1b(i,j)=rfcouple(bx,i,j,alpha1,U1)
end do
end do 
!write(*,*) 'Rabi frequency='
do mfi=-4,3

do i=1,NumSym2bstates(mfi+1)
do j=1,NumSym2bstates(mfi)

ialpha_a=nint(alphaSym2b(mfi+1,i,1))
ialpha_b=nint(alphaSym2b(mfi+1,i,2))
ialpha_ap=nint(alphaSym2b(mfi,j,1))
ialpha_bp=nint(alphaSym2b(mfi,j,2))


xSum = 0.d0
do iat_a = 1,NumHypfStates1
mfa = nint(alpha1(iat_a,2))
do iat_b = 1,NumHypfStates1
mfb = nint(alpha1(iat_b,2))
do iat_ap = 1,NumHypfStates1
mfap = nint(alpha1(iat_ap,2))
do iat_bp = 1,NumHypfStates1
mfbp = nint(alpha1(iat_bp,2))
if (mfa+mfb.eq.Mfi+1) then
if (mfap+mfbp.eq.Mfi) then

if (iat_b == iat_bp) then
xSum = xSum+Sym2b(iat_a,iat_b,ialpha_a,ialpha_b)* &
Sym2b(iat_ap,iat_bp,ialpha_ap,ialpha_bp)* &
hcouple1b(iat_a,iat_ap)

else if (iat_a == iat_ap) then
xSum = xSum+Sym2b(iat_a,iat_b,ialpha_a,ialpha_b)* &
Sym2b(iat_ap,iat_bp,ialpha_ap,ialpha_bp)* &
hcouple1b(iat_b,iat_bp)
end if



endif
endif
enddo
enddo
enddo
enddo

hcpack(mfi,i,j)=xsum

end do
end do


end do
write(*,*)'hcpack(0,2,7)(MHz)=',hcpack(0,2,7)/MHz,'hcpack(0,2,8)(MHz)=',hcpack(0,2,8)/MHz,'hcpack(-1,5,4)(MHz)=',hcpack(-1,5,4)/(MHz)
!*****************************************************
!do mfi=-4,3
!do i=1,8
!write(245,'(8e20.10)')(hcpack(mfi,i,j)/MHz,j=1,8)
!end do 
!end do 
!================================================================================================
!Constuct the asymptotic Hamiltonian in bare state basis
ndfloquet=0
!nblock=3
do mfi=-4,4
ndfloquet=ndfloquet+NumSym2bStates(mfi)
end do
999 continue
ndimension=ndfloquet*nblock

write(*,*) "ndimension=",ndimension


allocate(Hcomplete(ndimension,ndimension),energiesrf(ndimension),&
evectors(ndimension,ndimension),hfloquetblk(ndfloquet,ndfloquet),Hinterfloquet(ndfloquet,ndfloquet))

Hcomplete=0d0
energiesrf=0d0
evectors=0d0
Hfloquetblk=0d0
Hinterfloquet=0d0


jcount=0
do mfi=-4,4
do i=1,numsym2bstates(mfi)
jcount=jcount+1
if (mod(mfi,2)==0) then
Hfloquetblk(jcount,jcount)=thresholds(mfi,i)+0d0*wL
else
Hfloquetblk(jcount,jcount)=thresholds(mfi,i)+1d0*wL
end if
end do
end do

isum1=0
isum2=0
do mfi=-4,3
isum1=isum1+numsym2bstates(mfi)
do i=1,numsym2bstates(mfi+1)
do j=1,numsym2bstates(mfi)

ihf=isum1+i
jhf=isum2+j

Hfloquetblk(ihf,jhf)=hcpack(mfi,index(mfi+1,i),index(mfi,j))
Hfloquetblk(jhf,ihf)=Hfloquetblk(ihf,jhf)

if (mod(mfi,2)==0) then
Hinterfloquet(ihf,jhf)=hcpack(mfi,index(mfi+1,i),index(mfi,j))
else
Hinterfloquet(jhf,ihf)=hcpack(mfi,index(mfi+1,i),index(mfi,j))
end if

end do
end do
isum2=isum2+numsym2bstates(mfi)
end do

do i=1,nblock
Hcomplete(1+(i-1)*ndfloquet:i*ndfloquet,1+(i-1)*ndfloquet:i*ndfloquet)=Hfloquetblk
end do

do i=1,nblock-1
Hcomplete(1+(i-1)*ndfloquet:i*(ndfloquet),1+(i)*ndfloquet:(i+1)*ndfloquet)=Hinterfloquet
Hcomplete(1+(i)*ndfloquet:(i+1)*ndfloquet,1+(i-1)*ndfloquet:i*(ndfloquet))=transpose(Hinterfloquet)
end do

do i=1,nblock
do j=1,ndfloquet
Hcomplete(j+ndfloquet*(i-1),j+ndfloquet*(i-1))=Hcomplete(j+ndfloquet*(i-1),j+ndfloquet*(i-1))+2d0*wL*(i)
enddo
enddo

!1003 continue
!================================================================================================
!do i=1,ndimension
!write(211,'(400e20.10)') (Hcomplete(i,j),j=1,ndimension)
!end do 
!stop
!================================================================================================
!Diagonalize the asymptotic Hamiltonian
lwork=ndimension*5

evectors=Hcomplete

deallocate(Hcomplete,hfloquetblk,hinterfloquet)

allocate(work2(lwork))
call dsyev('V','U',ndimension,evectors,ndimension,Energiesrf,work2,lwork,info)
if (abs(info)>1d-8) then
write(*,*)"info=/=0",info
stop
end if

!if (abs(Bfield/gauss-9d0)<1d-6) then
!open(456,file='eveccomp.dat')
!open(455,file='hcouplecomp.dat')
!do i=1,ndimension
!write(456,'(400e20.10)')(evectors(i,j),j=1,ndimension)
!end do 
!do i=1,n2bstates
!write(455,'(50e20.10)')(hcouple(i,j),j=1,n2bstates)
!end do 
!end if
!write(457,'(400e20.10)')Bfield/gauss,energiesrf(1:ndimension)
!stop

!ndimension =n2bstates
!write(*,*)"n2bstates=",n2bstates,"ndimension=",ndimension
Thdress(1:ndimension)=Energiesrf
evectorssav(1:ndimension,1:ndimension)=evectors

!Thdress(1:20)=Energiesrf(21:40)

!write(*,*)"ndimension=",ndimension
!================================================================================================


!================================================================================================
!Construct the projection matrix
Pdress=0d0
do idress=1,ndimension
do jdress=1,ndimension

icount=0
jcount=0

do iblk=1,nblock
do mfi=-4,4
do istate=1,numsym2bstates(mfi)
icount=icount+1
jcount=0

do jblk=1,nblock
do mfj=-4,4
do jstate=1,numsym2bstates(mfj)
jcount=jcount+1

if (iblk==jblk .and. mfi==mfj) then
Pdress(0,idress,jdress)=Pdress(0,idress,jdress)+evectors(icount,idress)*evectors(jcount,jdress)*VQcal(mfi,0,istate,jstate)
Pdress(1,idress,jdress)=Pdress(1,idress,jdress)+evectors(icount,idress)*evectors(jcount,jdress)*VQcal(mfi,1,istate,jstate)
end if

end do
end do
end do

end do
end do
end do

end do
end do
!write(*,*) "good"

!write(10,*) Energiesrf

!do i=1,ndimension
!write(12,103) (evectors(i,j),j=1,ndimension)
!end do

!write(*,*) "Press is correct"
103 format(105e20.12)

do i=1,ndimension
do j=i+1,ndimension
if(abs(Pdress(0,i,j)+Pdress(1,i,j))>1d-8) write(11,*) i,j,Pdress(0,i,j)+Pdress(1,i,j)
end do
end do

!do j=1,10
!write(*,*) j,Pdress(0,j,j)+Pdress(1,j,j)
!end do
deallocate(evectors,energiesrf)
!================================================================================================

end subroutine


!================================================================================================
double precision function rfcouple(bx,ialpha_a,ialpha_ap,alpha1,U1)
implicit double precision (a-h,o-z)

double precision::alpha1(100,6),U1(8,8)
double precision::i1,s1,gE1,muE1,gN1,muN1
save i1,s1,gE1,muE1,gN1,muN1


if (ialpha_a==-1) then
call Parameters_87Rb(i1,s1,xi1,gE1,muE1,gN1,muN1,NumHypfStates1)
rfcouple=0d0
return
end if


if (abs(abs(alpha1(ialpha_a,2)-alpha1(ialpha_ap,2))-1d0)>1d-10) then
rfcouple=0d0
return
end if
!write(*,*) (U1(i,1),i=1,8)
!write(*,*) i1,s1,ge1,mue1,gn1,mun1

xSum = 0.d0
mis_a = 0
do mi_a = -nint(2.d0*i1),nint(2.d0*i1),2
do ms_a = -nint(2.d0*s1),nint(2.d0*s1),2
mis_a = mis_a+1
mis_ap = 0

do mi_ap = -nint(2.d0*i1),nint(2.d0*i1),2
do ms_ap = -nint(2.d0*s1),nint(2.d0*s1),2
mis_ap = mis_ap+1


if  (mi_a.eq.mi_ap+2 .and. ms_a.eq.ms_ap ) then

xSum = xSum+U1(mis_a,ialpha_a)*&
U1(mis_ap,ialpha_ap)* &
(-gN1*muN1*Bx/2d0*sqrt((i1-mi_ap/2d0)*(i1+mi_ap/2d0+1))/4d0)

else if  (mi_a.eq.mi_ap-2 .and. ms_a.eq.ms_ap ) then

xSum = xSum+U1(mis_a,ialpha_a)*&
U1(mis_ap,ialpha_ap)* &
(-gN1*muN1*Bx/2d0*sqrt((i1+mi_ap/2d0)*(i1-mi_ap/2d0+1))/4d0)

else if  (mi_a.eq.mi_ap .and. ms_a.eq.ms_ap+2 ) then

xSum = xSum+U1(mis_a,ialpha_a)*&
U1(mis_ap,ialpha_ap)* &
(-ge1*mue1*Bx/2d0*sqrt((s1-ms_ap/2d0)*(s1+ms_ap/2d0+1))/4d0)

else if  (mi_a.eq.mi_ap .and. ms_a.eq.ms_ap-2 )then

xSum = xSum+U1(mis_a,ialpha_a)*&
U1(mis_ap,ialpha_ap)* &
(-ge1*mue1*Bx/2d0*sqrt((s1+ms_ap/2d0)*(s1-ms_ap/2d0+1))/4d0)

endif


enddo
enddo
enddo
enddo

rfcouple=xSum

end function rfcouple
