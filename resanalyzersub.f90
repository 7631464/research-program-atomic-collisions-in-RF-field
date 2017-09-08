subroutine resanalyzer(bdeparray,numpoints)

implicit double precision (a-h,o-z)
double precision::bdeparray(numpoints,2)
integer::iblarge(10), ibsmall(10),ibzero(10)


open(63,file='resinfo2.dat')

iblarge=0
ibsmall=0
ibzero=0

istep1=0
istep2=0
do i=2,numpoints-1
if (bdeparray(i,2)>=bdeparray(i-1,2) .and. bdeparray(i,2)>=bdeparray(i+1,2)) then
istep1=istep1+1
iblarge(istep1)=i
endif
if (bdeparray(i,2)<=bdeparray(i-1,2) .and. bdeparray(i,2)<=bdeparray(i+1,2)) then
istep2=istep2+1
ibsmall(istep2)=i
endif
end do

do i=1,istep2
j=0
do while(bdeparray(ibsmall(i)+j,2)<0d0 .and. ibsmall(i)+j<=numpoints)
j=j+1
end do
if (j>0) then
ibzero(i)=ibsmall(i)+j
else
continue
end if
enddo


!write(*,*) istep1,istep2
if (istep2/=istep1) write(*,*) 'feshbach resonace not resolved'


do i=1,istep2
bres=(bdeparray(ibsmall(i),1)+bdeparray(iblarge(i),1))/2d0

if (ibzero(i)>0) then
b0=(bdeparray(ibzero(i),1)+bdeparray(ibzero(i)-1,1))/2d0
deltab=abs(bres-b0)
tunability=abs(bdeparray(ibsmall(i),2)-bdeparray(iblarge(i),2))
gamab=abs(bdeparray(ibsmall(i),1)-bdeparray(iblarge(i),1))
write(63,*) 'res',i, '=',bres,deltab,tunability,gamab
else
tunability=abs(bdeparray(ibsmall(i),2)-bdeparray(iblarge(i),2))
gamab=abs(bdeparray(ibsmall(i),1)-bdeparray(iblarge(i),1))
write(63,*) 'res',i, '=',bres,tunability,gamab
end if
end do 

end subroutine

subroutine stateinterpreter(istate,nblock,istateinfo)
!the floquet block is constructed as (mf=odd,n=1) and (mf=even,n=0), adjacent floquet blocks differ by 2*hv(two photon energy)
!Another way of constructing the floquet block is mf=n
!each floquet block includes all mf=-4 to 4 hyperfine states, arranged with the increment of mf. States with the same mf are in energy accending order.
implicit double precision (a-h,o-z)
integer:: f1,mf1,f2,mf2,n,istate,nblock
integer:: nsarray(-4:4)
integer:: istateinfo(5)

ifb=istate/36+1
istate1=mod(istate,36)
if (istate1==0) then
ifb=ifb-1
istate1=36
end if

nsarray=(/1,2,5,6,8,6,5,2,1/)

nsum=0
do mfi=-4,4
nsum=nsum+nsarray(mfi)
if (istate1<=nsum) then
mftarget=mfi
istate2=istate1+nsarray(mfi)-nsum
exit
end if
enddo

!Construct floquet block as mf=odd, n=1, mf=even, n=0
if (mod(mftarget,2)==0) then
n=0+(ifb-(nblock+1)/2)*2
else
n=1+(ifb-(nblock+1)/2)*2
end if


!Construct floquet block as n=mf
!n=mftarget+(ifb-(nblock+1)/2)*2



select case (mftarget)

case (-4)
istateinfo=(/2,-2,2,-2,n/)

case (-3)
select case (istate2)
case (1)
istateinfo=(/1,-1,2,-2,n/)
case (2)
istateinfo=(/2,-1,2,-2,n/)
end select

case (-2)
select case (istate2)
case (1)
istateinfo=(/1,-1,1,-1,n/)
case (2)
istateinfo=(/1,0,2,-2,n/)
case (3)
istateinfo=(/1,-1,2,-1,n/)
case (4)
istateinfo=(/2,-2,2,0,n/)
case (5)
istateinfo=(/2,-1,2,-1,n/)
end select


case (-1)
select case (istate2)
case (1)
istateinfo=(/1,0,1,-1,n/)
case (2)
istateinfo=(/1,1,2,-2,n/)
case (3)
istateinfo=(/1,0,2,-1,n/)
case (4)
istateinfo=(/1,-1,2,0,n/)
case (5)
istateinfo=(/2,-2,2,1,n/)
case (6)
istateinfo=(/2,-1,2,0,n/)
end select


case (0)
select case (istate2)
case (1)
istateinfo=(/1,0,1,0,n/)
case (2)
istateinfo=(/1,1,1,-1,n/)
case (3)
istateinfo=(/1,1,2,-1,n/)
case (4)
istateinfo=(/1,0,2,0,n/)
case (5)
istateinfo=(/1,-1,2,1,n/)
case (6)
istateinfo=(/2,-2,2,2,n/)
case (7)
istateinfo=(/2,-1,2,1,n/)
case (8)
istateinfo=(/2,0,2,0,n/)
end select

case (1)
select case (istate2)
case (1)
istateinfo=(/1,1,1,0,n/)
case (2)
istateinfo=(/1,1,2,0,n/)
case (3)
istateinfo=(/1,0,2,1,n/)
case (4)
istateinfo=(/1,-1,2,2,n/)
case (5)
istateinfo=(/2,-1,2,2,n/)
case (6)
istateinfo=(/2,0,2,1,n/)
end select

case (2)
select case (istate2)
case (1)
istateinfo=(/1,1,1,1,n/)
case (2)
istateinfo=(/1,1,2,1,n/)
case (3)
istateinfo=(/1,0,2,2,n/)
case (4)
istateinfo=(/2,0,2,2,n/)
case (5)
istateinfo=(/2,1,2,1,n/)
end select

case (3)
select case (istate2)
case (1)
istateinfo=(/1,1,2,2,n/)
case (2)
istateinfo=(/2,1,2,2,n/)
end select

case (4)
istateinfo=(/2,2,2,2,n/)

end select

end subroutine


