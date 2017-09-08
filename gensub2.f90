!************************************************************************
!      gensub.f  -- standard subroutines, to calculate wigner coeffs,   
!			interpolation, root-finding, numerical integration                  
!			Coulomb function programs also.                                     
!**********************************************************C            
!                                                                      
! WIGNER 3J, 6J, AND 9J FUNCTION SUBPROGRAMS.                           
! NOTE:  ALL ANGULAR MOMENTUM QUANTUM NUMBERS SUPPLIED TO THESE         
!        FUNCTIONS ARE INTEGERS WHICH ARE TWICE THE VALUE OF THE        
!        ACTUAL ANGULAR MOMENTUM DESIRED.  (THIS ALLOWS FOR HALF-       
!        INTEGER VALUES CONVENIENTLY.)  ALSO YOU MUST CALL SETUP        
!        ONE TIME BEFORE CALLING ANY OF THESE SO THAT THE RELE-         
!        VANT FACTORIALS CAN BE CALCULATED ONCE AND FOR ALL AND         
!        STORED, AS IN THE ABOVE EXAMPLE.                               
!                                                                       
!**********************************************************C            
!                                                                       
function coeff(l1,l2,l1p,l2p,l,k)
	implicit real*8(a-h,o-z)                                               
	front=(2*l1+1)*(2*l2+1)*(2*l1p+1)*(2*l2p+1)                            
	front=dsqrt(front)*(-1)**(l1+l1p+l)                                    
	l1d=2*l1                                                               
	l2d=2*l2                                                               
	l1pd=2*l1p                                                             
	l2pd=2*l2p                                                             
	ld=2*l                                                                 
	kd=2*k                                                                 
	iz=0                                                                   
	t1=thrj(l1d,kd,l1pd,iz,iz,iz)                                          
	t2=thrj(l2d,kd,l2pd,iz,iz,iz)                                          
	s1=sixj(l1d,l2d,ld,l2pd,l1pd,kd)                                       
	coeff=front*t1*t2*s1                                                   
	return                                                                 
	end        function
!                                                                      
      FUNCTION XNINEJ(J11,J12,J13,J21,J22,J23,J31,J32,J33)              
      IMPLICIT REAL*8(A-H,O-Z)                                          
      KMIN1 = IABS(J11-J33)                                             
      KMIN2 = IABS(J32-J21)                                             
      KMIN3 = IABS(J23-J12)                                             
      KMAX1 = J11+J33                                                   
      KMAX2 = J32+J21                                                   
      KMAX3 = J23+J12                                                   
!                                                                       
      IF(KMIN2.GT.KMIN1) KMIN1=KMIN2                                    
       IF(KMIN3.GT.KMIN1) KMIN1=KMIN3                                   
      IF(KMAX2.LT.KMAX1) KMAX1=KMAX2                                    
      IF(KMAX3.LT.KMAX1) KMAX1=KMAX3                                    
!                                                                       
      KMIN1 = KMIN1 + 1                                                 
      KMAX1 = KMAX1 + 1                                                 
       XNINEJ = 0.D0                                                    
      IF(KMIN1.GT.KMAX1) GO TO 1000                                     
      DO 100 K1 = KMIN1,KMAX1,2                                         
      K = K1 - 1                                                        
      S1 = SIXJ(J11,J21,J31,J32,J33,K)                                  
      S2 = SIXJ(J12,J22,J32,J21,K,J23)                                  
      S3 = SIXJ(J13,J23,J33,K,J11,J12)                                  
      P = (K+1)*((-1)**K)                                               
      XNINEJ = XNINEJ + P*S1*S2*S3                                      
  100 CONTINUE                                                          
 1000 CONTINUE                                                          
      RETURN                                                            
      END           function
!                                                                       
      FUNCTION THRJ(J1D,J2D,J3D,M1D,M2D,M3D)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
      X1 = J1D/2.D0                                                     
      X2 = J2D/2.D0                                                     
      X3 = J3D/2.D0                                                     
      Y1 = M1D/2.D0                                                     
      Y2 = M2D/2.D0                                                     
      Y3 = M3D/2.D0                                                     
!                                                                       
! -- NEXT COME THE TRIANGULARITY CHECKS:                                
!                                                                       
      IF(J1D+J2D-J3D.LT.0) GO TO 9998                                   
      IF(J2D+J3D-J1D.LT.0) GO TO 9998                                   
      IF(J3D+J1D-J2D.LT.0) GO TO 9998                                   
      IF(J3D.LT.IABS(J1D-J2D)) GO TO 9998                               
      IF(M1D+M2D+M3D.NE.0) GO TO 9998                                   
      LLL = J1D+J2D+J3D                                                 
      IF(2*(LLL/2) .NE. LLL) GO TO 9998                                 
!                                                                       
      KMIN = (J3D-J1D-M2D)/2                                            
      KMIN1 = KMIN                                                      
      KMIN2 = (J3D-J2D+M1D)/2                                           
      IF(KMIN2.LT.KMIN) KMIN=KMIN2                                      
      KMIN = (-1)*KMIN                                                  
      KMAX = X1+X2-X3 +0.1D0                                            
      KMAX1 = KMAX                                                      
      KMAX2 = X1-Y1                                                     
      KMAX3 = X2+Y2                                                     
      IF(KMAX2.LT.KMAX) KMAX=KMAX2                                      
      IF(KMAX3.LT.KMAX) KMAX=KMAX3                                      
      IF(KMIN.LT.0) KMIN = 0                                            
      IF(KMIN.GT.KMAX) GO TO 9998                                       
!                                                                       
      JMIN = KMIN+1                                                     
      JMAX = KMAX+1                                                     
      TERM1 = FRONTL(X1,X2,X3,Y1,Y2,Y3)                                 
	iphase=iabs((j1d-j2d-m3d)/2)                                           
	msign=(-1)**iphase                                                     
!g     MSIGN = (-1)**((J1D-J2D-M3D)/2)                                  
      SUM = 0.D0                                                        
      DO 10 I1 = JMIN,JMAX                                              
      I = I1 - 1                                                        
      TERM2 = FL(I1)+FL(KMIN1+I1)+FL(KMIN2+I1)                          
      TERM2 = TERM2+FL(KMAX1-I+1)+FL(KMAX2-I+1)+FL(KMAX3-I+1)           
      TERM= DEXP(TERM1-TERM2)                                           
      TERM = TERM*MSIGN*((-1)**I)                                       
      SUM = SUM + TERM                                                  
  10  CONTINUE                                                          
      THRJ = SUM                                                        
      GO TO 9999                                                        
 9998 THRJ = 0.D0                                                       
 9999 CONTINUE                                                          
      RETURN                                                            
      END    function
!                                                                       
      FUNCTION FL(I)                                                    
       IMPLICIT REAL*8(A-H,O-Z)                                         
!C!   DIMENSION FACL(60)                                                
      COMMON/FACTOR/FACL(200)                                           
      FL = FACL(I)                                                      
      RETURN                                                            
      END function
!                                                                       
!* ** **                                                               
!- THIS SUBROUTINE INITIALIZES BY FINDING THE LOGARITHM                
!--OF THE FIRST 199 FACTORIALS AND STORING THEM.                       
!                                                                       
      SUBROUTINE SETUP                                                  
      IMPLICIT REAL*8(A-H,O-Z)                                          
      COMMON/FACTOR/FACL(200)                                           
      N = 199                                                           
      FACL(1) = 0.D0                                                    
      DO 100 I = 2,N                                                    
      I1 = I - 1                                                        
      FACL(I) = FACL(I1) + DLOG(I1*1.D0)                                
  100 CONTINUE                                                          
      RETURN                                                            
      END      subroutine
!* ** **                                                               
!                                                                       
      FUNCTION FRONTL(X1,X2,X3,Y1,Y2,Y3)                                
      IMPLICIT REAL*8(A-H,O-Z)                                          
      L1 = X1+X2-X3 +1.1D0                                              
      L2 = X2+X3-X1 +1.1D0                                              
      L3 = X3+X1-X2 +1.1D0                                              
      L4 = X1+X2+X3+1+1.1D0                                             
      L5 = X1+Y1+1.1D0                                                  
      L6 = X1-Y1+1.1D0                                                  
      L7 = X2+Y2+1.1D0                                                  
      L8 = X2-Y2+1.1D0                                                  
      L9 = X3+Y3+1.1D0                                                  
      L10 = X3-Y3+1.1D0                                                 
      FRONTL = FL(L1)+FL(L2)+FL(L3)-FL(L4)+FL(L5)+FL(L6)                
      FRONTL = FRONTL +FL(L7)+FL(L8)+FL(L9)+FL(L10)                     
      FRONTL = FRONTL/2.D0                                              
      RETURN                                                            
      END           function
!                                                                       
      FUNCTION SIXJ(J1D,J2D,J3D,J4D,J5D,J6D)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
!                                                                       
! -- CHECK THAT TRIANGULARITY CONDITIONS ARE SATISFIED.                 
!                                                                       
      IF(J3D.LT.IABS(J1D-J2D) .OR. J3D.GT.J1D+J2D) GO TO 9998           
      IF(J6D.LT.IABS(J4D-J2D) .OR. J6D.GT.J4D+J2D) GO TO 9998           
      IF(J6D.LT.IABS(J1D-J5D) .OR. J6D.GT.J1D+J5D) GO TO 9998           
      IF(J3D.LT.IABS(J4D-J5D) .OR. J3D.GT.J4D+J5D) GO TO 9998           
      K1=J1D+J2D+J3D                                                    
      K2=J4D+J2D+J6D                                                    
      K3=J6D+J1D+J5D                                                    
      K4=J3D+J4D+J5D                                                    
      IF(2*(K1/2).NE.K1 .OR. 2*(K2/2).NE.K2) GO TO 9998                 
      IF(2*(K3/2).NE.K3 .OR. 2*(K4/2).NE.K4) GO TO 9998                 
!                                                                       
! -- NOW GO AHEAD AND CALCULATE THE SIXJ.                               
!                                                                       
      JM1 = (J1D+J2D+J3D)/2                                             
      JM2 = (J1D+J5D+J6D)/2                                             
      JM3 = (J4D+J2D+J6D)/2                                             
      JM4 = (J4D+J5D+J3D)/2                                             
      JX1 = (J1D+J2D+J4D+J5D)/2                                         
      JX2 = (J2D+J3D+J5D+J6D)/2                                         
      JX3 = (J3D+J1D+J6D+J4D)/2                                         
!                                                                       
      JM = JM1                                                          
      IF(JM2.GT.JM) JM = JM2                                            
      IF(JM3.GT.JM) JM = JM3                                            
      IF(JM4.GT.JM) JM = JM4                                            
      JX = JX1                                                          
      IF(JX2.LT.JX) JX = JX2                                            
      IF(JX3.LT.JX) JX = JX3                                            
      KM = JM+1                                                         
      KX = JX+1                                                         
      IF(KM.GT.KX) GO TO 9998                                           
      TERM1 = FRTSJL(J1D,J2D,J3D,J4D,J5D,J6D)                           
      SIXJ = 0.D0                                                       
      DO 10 I1 = KM,KX                                                  
      I = I1 - 1                                                        
      TERM2 = FL(I+2)-FL(I+1-JM1)-FL(I+1-JM2)-FL(I+1-JM3)               
      TERM2 = TERM2-FL(I+1-JM4)-FL(JX1-I+1)-FL(JX2-I+1)                 
      TERM2 = TERM2 -FL(JX3-I+1)                                        
      TERM = DEXP(TERM1+TERM2) * ((-1)**I)                              
      SIXJ = SIXJ + TERM                                                
   10 CONTINUE                                                          
      GO TO 9999                                                        
 9998 CONTINUE                                                          
      SIXJ = 0.D0                                                       
 9999 CONTINUE                                                          
      RETURN                                                            
      END       function
!                                                                       
      FUNCTION FRTSJL(J1D,J2D,J3D,J4D,J5D,J6D)                          
      IMPLICIT REAL*8(A-H,O-Z)                                          
      FRTSJL = DL(J1D,J2D,J3D) + DL(J1D,J5D,J6D)                        
      FRTSJL = FRTSJL + DL(J4D,J2D,J6D) + DL(J4D,J5D,J3D)               
      RETURN                                                            
      END    function
!                                                                       
      FUNCTION DL(J1D,J2D,J3D)                                          
      IMPLICIT REAL*8(A-H,O-Z)                                          
      L1 = (J1D+J2D-J3D)/2                                              
      L2 = (J2D+J3D-J1D)/2                                              
      L3 = (J3D+J1D-J2D)/2                                              
      L4 = (J1D+J2D+J3D)/2 + 1                                          
      DL = FL(L1+1)+FL(L2+1)+FL(L3+1)-FL(L4+1)                          
      DL = DL/2.D0                                                      
      RETURN                                                            
      END  function
!                                                                       
      FUNCTION CLEBSH(J1D,J2D,JD,M1D,M2D,MD)                            
      IMPLICIT REAL*8(A-H,O-Z)                                          
      CLEBSH = 0.D0                                                     
      IF(M1D+M2D.NE.MD) GO TO 100                                       
      MMD = -MD                                                         
      Q = THRJ(J1D,J2D,JD,M1D,M2D,MMD)                                  
      PHASE = ((-1)**(20 + (J2D-J1D-MD)/2))                             
      CLEBSH = Q*PHASE*DSQRT(JD+1.D0)                                   
  100 CONTINUE                                                          
      RETURN                                                            
      END    function
! ***********************************************************************

! *********************************************************************  
	function rint(f,na,nb,nq,h)                                            
!                                                                       
! -- function rint performs a numerical integration over the function f,
! -----  assuming an equally spaced x-axis grid with step size h and a  
! -----  total of nb mesh points.  Use na=1, nq=10, and nb>20 always.   
!                                                                       
	implicit real*8(a-h,o-z)                                               
	dimension c(55),d(10),f(nb)                                            
data c/1.d0,   2.d0,1.d0,   23.d0, 28.d0, 9.d0,   25.d0, 20.d0, 31.d0, 8.d0,   1413.d0, 1586.d0, 1104.d0, 1902.d0, 475.d0,&
1456.d0, 1333.d0, 1746.d0, 944.d0,1982.d0,459.d0,  119585.d0, 130936.d0, 89437.d0, 177984.d0, 54851.d0, 176648.d0,36799.d0,&
122175.d0, 111080.d0, 156451.d0,46912.d0,220509.d0,  29336.d0, 185153.d0, 35584.d0,     7200319.d0, 7783754.d0, 5095890.d0,&
12489922.d0,-1020160.d0,   16263486.d0, 261166.d0, 11532470.d0, 2082753.d0,  7305728.d0, 6767167.d0, 9516362.d0,1053138.d0,&
18554050.d0,       -7084288.d0, 20306238.d0, -1471442.d0, 11965622.d0,   2034625.d0/                                                     
data d/  2.d0,    2.d0, 24.d0, 24.d0,1440.d0,1440.d0,120960.d0,120960.d0, 7257600.d0, 7257600.d0/
	a=0.d0                                                                 
	l=na                                                                   
	m=nb                                                                   
	i=nq*(nq+1)/2                                                          
	do 10 j=1,nq                                                           
	a=a+c(i)*(f(l)+f(m))                                                   
	l=l+1                                                                  
	m=m-1                                                                  
   10   i=i-1                                                           
	a=a/d(nq)                                                              
	do 20 n=l,m                                                            
   20   a=a+f(n)                                                        
	rint=a*h                                                               
	return                                                                 
	end       function
!**********************************************************************
!**********************************************************************
        double precision function delt(i,j)                             
        implicit double precision (a-h,o-z)                             
        delt=0.d0                                                       
        if(i.eq.j) delt=1.d0                                            
        return                                                          
        end    function
!**********************************************************************
        subroutine root(xtry,hh,iroute,ep,xest,yest,yfcn)               
        implicit real*8(a-h,o-z)                                        
                                                                        
                                                                        
!**  subroutine root locates a root of yfcn(x,y=0),                    
!**  returning finally (xest,yest)                                     
!**  xtry is initial guess,hh is step size for the iteration.          
!**  iroute=0 means hh is kept constant in the initial (coarse)        
!**  root location loop.                                               
!**  iroute=1 means hh can change sign but not magnitude.              
!**  iroute=2 means both magnitude & sign of hh can vary.              
!**  after the initial loop, a 3-point method is used to               
!**  speed convergence.                                                
!**                                                                    
        iq = 3                                                          
        u1 = xtry                                                       
	f1=yfcn(u1)                                                            
        mtest = 0                                                       
        do 420 ki = 1,30                                                
        if(ki.eq.2 .and. iroute.eq.2) go to 421                         
        if(mtest.ne.0) go to 920                                        
        u2 = u1 + hh                                                    
  421   f2=yfcn(u2)                                                     
        if(f1*f2 .le. 0.d0) mtest = 1000                                
        if(mtest.ne.0)  go to 3055                                      
        q1 = dabs(f1)                                                   
        q2 = dabs(f2)                                                   
        if(ki.eq.1 .and. iroute.eq.2) go to 3047                        
        if(iroute-1) 3033,3044,3044                                     
 3044   if(q1.ge.q2 .and. mtest.eq.0) u1=u2                             
        if(q1.ge.q2 .and. mtest.eq.0) f1=f2                             
        if(q2.gt.q1 .and. mtest.eq.0) hh = (-1)*hh                      
        go to 3055                                                      
 3033   if(mtest.eq.0) u1 = u2                                          
        if(mtest.eq.0) f1 = f2                                          
        go to 3055                                                      
 3047   slp = (f2-f1)/(u2-u1)                                           
        uest = u1 - f1/slp                                              
        if(q2.lt.q1) f1 = f2                                            
        if(q2.lt.q1) u1 = u2                                            
        u2 = uest                                                       
 3055   continue                                                        
  420   continue                                                        
        if(mtest.ne.1000)  go to 8500                                   
!                                                                       
  920   continue                                                        
! cccccccccccccccccccccccccccccccc      ep = 1.d-06                      
        mchk = 0                                                        
        do 440 i = 1,50                                                 
        if(mchk.ne.0) go to 8600                                        
        if(i.gt.iq) go to 429                                           
        slp = (f2-f1)/(u2-u1)                                           
        uest = u1 - f1/slp                                              
        if(i.lt.iq) uest = (u1+u2)/2.d0                                 
        u3 = uest                                                       
        q3 = u3                                                         
        q2 = u2                                                         
        q1 = u1                                                         
  429   fest=yfcn(uest)                                                 
        f3 = fest                                                       
        u3 = uest                                                       
        if(i.ge.iq) go to 439                                           
        if(fest*f1 .lt. 0.d0) f2=fest                                   
        if(fest*f1 .lt. 0.d0) u2=uest                                   
        if(fest*f1 .ge. 0.d0) f1=fest                                   
        if(fest*f1 .ge. 0.d0) u1=uest                                   
        go to 440                                                       
  439   continue                                                        
!                                                                       
        aa = (u2-u1)*f3+(u1-u3)*f2+(u3-u2)*f1                           
        denom = (u3-u2)*(u2-u1)*(u2-u1)                                 
        aa = aa/denom                                                   
        bb = (u2-u1)*(2.d0*u3-u2-u1)*f3                                 
        bb = bb- (u3-u1)*(u3-u1)*f2+(u3-u2)*(u3-u2)*f1                  
        bb = bb/denom                                                   
        cc = (u3-u1)*f3/(u2-u1)                                         
        q0q = dabs(bb*bb - 4.d0*aa*cc)                                  
        q0q = dsqrt(q0q)                                                
        den1 = bb+q0q                                                   
        den2 = bb-q0q                                                   
        qden1 = dabs(den1)                                              
        qden2 = dabs(den2)                                              
        if(qden2.ge.qden1) den1 = den2                                  
        u1 = u2                                                         
        u2 = u3                                                         
        u3 = u2 - 2.d0*cc/den1                                          
        f1 = f2                                                         
        f2 = f3                                                         
        tst = dabs(u3-u2)                                               
        if(tst.lt.ep) mchk = 1000                                       
        uest = u3                                                       
  440   continue                                                        
        if(mchk.ne.1000) go to 8500                                     
        go to 8600                                                      
 8500   write(6,8501)                                                   
 8501   format(3x,'no convergence was achieved')                        
 8600   continue                                                        
        xest = uest                                                     
        yest = fest                                                     
        return                                                          
        end      subroutine
!********************************************************************** 
! ********                                                              
!	C. Greene, 6-8-87   -- modified 11-9-87                               
!  The following is a main program which shows how to use subroutine sea
!                                                                       
!       The form of the call is:                                        
!                                                                       
!           call seaton(l,eryd,r,zion,f,fp,g,gp)                             
!                                                                       
!      Here l=orbital angular momentum                                       
!           eryd = electron energy in rydbergs                               
!            r=radial distance                                          
! 	     zion=nuclear charge                                             
!      (CAUTION:  I have only tested my transformation program for zion=1,   
! 		   and I suspect it must be modified for other cases.      ) 
!           (f,g)=(regular,irregular) energy normalized Coulomb functions.   
!           (fp,gp)= (d/dr)(f,g).                                            
!           Check:  W(f,g)=f*gp-fp*g  should be 2/pi if the program works.   
!                                                                       
! *********                                                             
!                                                                       
	subroutine dummy                                                       
	implicit real*8(a-h,o-z)                                               
	data pi/3.14159265358979d0/                                            
	l=0                                                                    
	eryd=-1.d0/(1.5d0**2) + 1.d-10                                         
	zion=1.d0                                                              
	write(6,*) 'e=',eryd                                                   
	do 100 ir=1,10                                                         
		r=ir*0.5d0                       
		call seaton(l,eryd,r,zion,f,fp,g,gp)    
		wfgtst=f*gp-fp*g - 2.d0/pi       
		write(6,*) r,g,wfgtst                
  100	continue                                                          
	epos=-eryd                                                             
	write(6,*) 'e=',epos                                                   
	do 200 ir=1,10                                                         
		r=ir*0.5d0           
		call seaton(l,epos,r,zion,f,fp,g,gp)    
		wfgtst=f*gp-fp*g - 2.d0/pi         
		write(6,*) r,g,wfgtst      
  200 	continue                                                         
	stop                                                                   
	end                                                                    
!                                                                       
      subroutine seaton1(l,eryd,r,zion,f,fp,g,gp)                       
      implicit real*8(a-h,o-z)                                          
      acc=1.d-10                                                        
      rl=l                                                              
      rr=r*zion                                                         
      eps=eryd/(zion**2)                                                
      call coulfg(l,eps,rr,acc,f0,f0p,g0,g0p,k,ierr,actacc)             
      if(.not.(eryd.lt.0))goto 23000                                    
         ea=dabs(eryd)                                                  
         call ganda(a,gg,l,ea,zion,999)                                 
         goto 23001                                                     
!     else                                                              
23000    continue                                                       
         gam=1.d0/dsqrt(eps)                                            
         call gcmplx(a,gg,rl,gam)                                       
23001 continue                                                          
      a5=dsqrt(dabs(a))                                                 
      f=a5*f0                                                           
      fp=a5*f0p                                                         
      g=(g0+gg*f0)/a5                                                   
      gp=(g0p+gg*f0p)/a5                                                
!                                                                       
! ** the next five lines changed on 1-22-88 by c.greene thanks to h. gao
!                                                                       
	factor = dsqrt(zion)                                                   
	f=f/factor                                                             
	g=g/factor                                                             
	fp=fp*factor                                                           
	gp=gp*factor                                                           
!                                                                       
      return                                                            
      end        subroutine
!                                                                       
!                                                                       
      SUBROUTINE GANDA(A,G,L,E,ZION,NOPT)                               
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DATA PI/3.1415926535897932D0/                                     
      DPI = PI*2.D0                                                     
!** THIS PROGRAM RETURNS THE QDT TRANSFORMATION PARAMETERS A&G.        
      IF(ZION.EQ.0) WRITE(6,90)                                         
   90 FORMAT(1X,'***** ZION = 0')                                       
      IF(ZION.EQ.0) RETURN                                              
      E = DABS(E)                                                       
      XNUE = ZION/DSQRT(E)                                              
!** EVALUATE A(K,L) FIRST.                                             
      A = 1.D0                                                          
      IF(L.EQ.0) GO TO 109                                              
      DO 100 I = 1,L                                                    
      A = A* (1.D0 -I*I*E/(ZION**2) )                                   
  100 CONTINUE                                                          
  109 continue                                                          
!** GIVE WARNINGS IN CASE A < OR = 0 .                                 
      IF(A.LE.0) WRITE(6,5555)                                          
 5555 FORMAT(1X,'****** A < OR = 0')                                    
      IF(NOPT.EQ.1) RETURN                                              
!** CHECK WHETHER XNUE = INTEGER.                                      
      N = XNUE + 1.D-02                                                 
      Z = XNUE - N                                                      
      IF(Z.EQ.0) G = 0.D0                                               
      IF(Z.EQ.0) RETURN                                                 
!** G(K,L) IS NOW EVALUATED USING THE DIGAMMA PROGRAM.                 
      G = A*(DIGAM(L+1+XNUE) + DIGAM(-L+XNUE)    - 2.D0 * DLOG(XNUE) )/DPI
      RETURN                                                            
      END                                                               
      SUBROUTINE GCMPLX(B,G,RL,GAM)                                     
      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)                            
!CC!  COMPLEX*16 ZQ1,ZQ2,ZQ3,ZSUM,ZG                                    
      COMPLEX*16 CDLOG                                                  
      DIMENSION GG(2)                                                   
      EQUIVALENCE (GG(1),ZG)                                            
      DATA PI/3.1415926535897932D0/                                     
!                                                                       
      ZQ1 = DCMPLX(RL+1.D0,GAM)                                         
      ZQ2 = DCMPLX(-RL,GAM)                                             
      ZQ3 = DCMPLX(0.D0,GAM)                                            
      ZSUM = ZDIGAM(ZQ1) + ZDIGAM(ZQ2) + CDLOG(ZQ3)*(-2.D0)             
      PROD = 1.D0                                                       
      L = RL                                                            
      DO 100 I = 1,L                                                    
      IF(L.EQ.0) GO TO 100                                              
      PROD = PROD*( 1.D0 + I*I/(GAM*GAM) )                              
  100 CONTINUE                                                          
      A1 = PROD                                                         
      ZG = ZSUM*A1/(2.D0*PI)                                            
      G = GG(1)                                                         
      QQ = 1.D0 - DEXP(-2.D0*PI*GAM)                                    
      B = A1 / QQ                                                       
      RETURN                                                            
      END        subroutine
!                                                                       
!      function cdabs(z)                                                
!      implicit real*8(a-h,o-y),complex*16(z)                           
!      cdabs = cabs(z)                                                  
!      return                                                           
!      end                                                              
!cc!                                                                     
!      function cdlog(z)                                                     
!      implicit complex*16(a-h,o-z)                                          
!      cdlog=clog(z)                                                         
!      return                                                                
!      end                                                                   
!                                                                       
!      function cdexp(z)                                                
!      implicit complex*16(a-h,o-z)                                     
!      cdexp=cexp(z)                                                    
!      return                                                           
!      end                                                              
!                                                                       
      FUNCTION DIGAM(ARGG)                                              
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIGAM = 0.D0                                                      
      ARG5 =-0.57721566490153286D0                                      
      EN = 1.D0                                                         
      ARG  = ARGG                                                       
      ARG2 = ARGG                                                       
    1 IF(ARG2-40.D0) 2,3,3                                              
    2 DIGAM = DIGAM - 1.D0/ARG                                          
      ARG = 1.D0+ ARG                                                   
      ARG2 = 1.D0 + ARG2                                                
      GO TO 1                                                           
    3 PSI = DLOG(ARG) - 1.D0/(2.D0*ARG) - 1.D0/(12.D0*ARG**2)           
      PSI = PSI + 1.D0/(120.D0*ARG**4) - 1.D0/(252.D0*ARG**6)           
      DIGAM = DIGAM + PSI                                               
      RETURN                                                            
      END            function
!     COMPLEX FUNCTION ZDIGAM*16(ARG)                                   
	function zdigam(arg)                                                   
      IMPLICIT REAL*8(A-H,O-Y),COMPLEX*16(Z)                            
      COMPLEX*16 ARG,ZDIGAM,cdlog                                       
      COMPLEX*8 ARGG                                                    
      REAL*8 INC                                                        
      ZDIGAM = (0.D0,0.D0)                                              
      ARG5 =-0.57721566490153286D0                                      
      PI = 3.1415926535897932D0                                         
      EN = 1.D0                                                         
      ARGG = ARG                                                        
      ARG2 = REAL(ARGG)                                                 
      ARG3 = AIMAG(ARGG)                                                
      IF(ARG3) 4,1,4                                                    
    1 IF(ARG2-40.D0) 2,3,3                                              
    2 ZDIGAM = ZDIGAM - 1.D0/ARG                                        
      ARG = 1.D0+ ARG                                                   
      ARG2 = 1.D0 + ARG2                                                
      GO TO 1                                                           
    3 PSI = CDLOG(ARG)-1.D0/(2.D0*ARG)-1.D0/(12.D0*ARG**2)              
      PSI=PSI +1.D0/(120.D0*ARG**4)-1.D0/(252.D0*ARG**6)                
      ZDIGAM = ZDIGAM + PSI                                             
      GO TO 12                                                          
    4 IF(ARG2) 5,7,6                                                    
    5 ZDIGAM = ZDIGAM - 1.D0/ARG                                        
      ARG = ARG + 1.D0                                                  
      ARG2 = ARG2 + 1.D0                                                
      GO TO 4                                                           
    6 ARG = ARG - 1.D0                                                  
      ARG2 = ARG2 - 1.D0                                                
      ZDIGAM = ZDIGAM + 1.D0/ARG                                        
      GO TO 4                                                           
    7 Y = CDABS(ARG)                                                    
      ARG7 = PI*Y                                                       
      ARG4 = 0.5D0/Y + (PI/2.D0)/DTANH(ARG7)                            
      IF(Y-20.D0) 8,10,10                                               
    8 INC = Y*Y/(EN*(EN*EN+Y*Y))                                        
      ARG5 = ARG5 + INC                                                 
      IF(INC - 1.D-12) 11,11,9                                          
    9 EN = EN + 1.D0                                                    
      GO TO 8                                                           
   10 ARG5 = 1.D0/(12.D0*Y**2) + 1.D0/(120.D0*Y**4)                     
      ARG5 = ARG5 + 1.D0/(252.D0*Y**6)+ DLOG(Y)                         
   11 ZDIGAM = DCMPLX(ARG5,ARG4) + ZDIGAM                               
!     XQ1 = REAL(ZDIGAM)                                                
!     XQ2 = AIMAG(ZDIGAM)                                               
   12 continue                                                          
	return                                                                 
      END       function
!**********************************************************************
!                                                                       
!                                                                       
      subroutine coulfg(ll,eps,rho,acc,f,fp,g,gp,k,ierr,actacc)         
!                                                                       
!  calculates coulomb functions f and g and their derivatives           
!                                                                       
!  input -                                                              
!        ll=angular momentum quantum number                             
!        eps=z-scaled energy in rydbergs                                
!        rho=z-scaled radial variable in atomi! units                   
!        acc=accuracy required                                          
!                                                                       
!  output -                                                             
!        f=regular function                                             
!        fp=derivative of f                                             
!        g=irregular function                                           
!        gp=derivative of g                                             
!        k=number of terms needed in expansion                          
!        ierr=error code                                                
!        actacc=accuracy actually achieved                              
!                                                                       
!  convergence criterion -                                              
!        value of wronskian converged to accuracy of 0.5*acc            
!                                                                       
!  error codes -                                                        
!        ierr=0, converged with actacc.lt.ac!                           
!        ierr=1, converged with actacc.gt.acc                           
!        ierr=2, not converged with 101 terms in main summation         
!                                                                       
!  initialization                                                       
!                                                                       
!                                                                       
      implicit real*8 (a-h,o-z)                                         
!d    delete previous card for single precision                         
!                                                                       
      data r2pi,ps0/.159154943,-.154431330/                             
      ierr=0                                                            
      lp1=ll+1                                                          
      l2=2*ll                                                           
      l2p1=l2+1                                                         
      fl=ll                                                             
      flp1=lp1                                                          
      fl2p1=l2p1                                                        
      e2=.5*eps                                                         
      r2=2.*rho                                                         
      acc2=2.*acc                                                       
!                                                                       
!     initialize fa=factorial(2*ll+1)                                   
!     and ps=psi(2*ll+2)+psi(1)                                         
!                                                                       
      fa=1.                                                             
      ps=ps0                                                            
!                                                                       
!                                                                       
!  calculate alpha(n) and beta(n) and initialize s and sp               
!  continue calculation of fa and ps                                    
!                                                                       
!     s and sp for n=0                                                  
      x3=-l2                                                            
      x2=l2p1                                                           
      x1=-2.*r2**(-lp1)                                                 
      sp=x3*x1                                                          
      x1=r2*x1                                                          
      s=x1                                                              
!                                                                       
!     initialize for coefficients in recursion formulae                 
      p1=fl*e2                                                          
      p2=p1                                                             
      q1=-e2                                                            
!                                                                       
!     initialize alpha and beta                                         
      alp1=1.                                                           
      alp2=1.+p2                                                        
      bet1=0.                                                           
      bet2=q1                                                           
!                                                                       
      if(ll.eq.0)goto 20                                                
!                                                                       
!     s and sp for n=1                                                  
      x3=x3+2.                                                          
      x2=x2-1.                                                          
      x1=x1/x2                                                          
      sp=sp+x3*x1                                                       
      x1=r2*x1                                                          
      s=s+x1                                                            
!                                                                       
!     loop for n=2 to 2*ll                                              
      do 10 n=2,l2                                                      
!                                                                       
!     continue calculation of fa and psi                                
      fn=n                                                              
      fa=fn*fa                                                          
      ps=ps+1./fn                                                       
!                                                                       
!     continue calculation of s and sp                                  
      x3=x3+2.                                                          
      x2=x2-1.                                                          
      x1=x1/(x2*fn)                                                     
      sp=sp+x3*x1*alp2                                                  
      x1=r2*x1                                                          
      s=s+x1*alp2                                                       
!                                                                       
!     compute coefficients in recursion formulae                        
      p1=p1-e2                                                          
      p2=p2+p1                                                          
      q1=q1-e2                                                          
!     now have p2=-n*(n-2*ll-1)*eps/4                                   
!     and q1=-n*eps/2                                                   
!                                                                       
!     new alpha and beta                                                
      alp0=alp1                                                         
      alp1=alp2                                                         
      alp2=alp1+p2*alp0                                                 
      bet0=bet1                                                         
      bet1=bet2                                                         
   10 bet2=bet1+p2*bet0+q1*alp0                                         
!                                                                       
!     normalize s and sp, complete calculation of fa and ps             
      s=s*fa                                                            
      sp=sp*fa                                                          
      fa=fl2p1*fa                                                       
      ps=ps+1./fl2p1                                                    
!                                                                       
!     complete calculation of alpha and beta                            
      p1=p1-e2                                                          
      p2=p2+p1                                                          
      q1=q1-e2                                                          
      alp0=alp1                                                         
      alp1=alp2                                                         
      bet0=bet1                                                         
      bet1=bet2                                                         
      bet2=bet1+p2*bet0+q1*alp0                                         
!                                                                       
   20 continue                                                          
!     now have alp1=alpha(2*ll+1)                                       
!     and bet1=beta(2*ll+1), bet2=beta(2*ll+2)                          
!                                                                       
!     value of a=a(eps,ll)                                              
      a=alp1                                                            
      a4=4.*a                                                           
      cl=2.*a*dlog(dabs(r2))                                            
!d    for single precision replace dlog by alog and dabs by abs         
      clp=2.*a/rho                                                      
!                                                                       
!  calculate a(n) and d(n), f and fp and                                
!  complete calculation of s and sp                                     
!                                                                       
!     calculate a0,a1,d0,d1                                             
      a0=(2.**lp1)/fa                                                   
      a1=-a0/flp1                                                       
      ps=2.*ps*a                                                        
      d0=(bet1-ps)*a0                                                   
      d1=(bet2-ps-(2.+1./flp1)*a)*a1                                    
!                                                                       
!     initialize f,fp, continue calculation of s,sp                     
!     - values for n=0                                                  
      fnplp1=flp1                                                       
      c1=rho**ll                                                        
      c1p=fnplp1*c1                                                     
      fp=c1p*a0                                                         
      sp=sp+c1p*d0                                                      
      c1=c1*rho                                                         
      f=c1*a0                                                           
      s=s+c1*d0                                                         
      w1=f*(clp*f+sp)-fp*s                                              
!                                                                       
!     - values for n=1                                                  
      fnplp1=fnplp1+1.                                                  
      c1p=fnplp1*c1                                                     
      fp=fp+c1p*a1                                                      
      sp=sp+c1p*d1                                                      
      c1=c1*rho                                                         
      f=f+c1*a1                                                         
      s=s+c1*d1                                                         
      w2=f*(clp*f+sp)-fp*s                                              
      dw2=dabs(w2-w1)                                                   
!d    for single precision replace dabs by abs                          
!                                                                       
!     initialize for coefficients in recursion formulae                 
      p1=-2.*flp1                                                       
      p2=p1                                                             
      q1=a4+2.*a*fl2p1                                                  
!                                                                       
!     loop for n=2 to 100                                               
      do 40 n=2,100                                                     
!                                                                       
!     compute coefficients in recursion formulae                        
      p1=p1-2.                                                          
      p2=p2+p1                                                          
      q1=q1+a4                                                          
!     now have p2=-n*(n+2*ll+1)                                         
!     and q1=2*a*(2*n+2*ll+1)                                           
!                                                                       
!     compute a2=a(n) and d2=d(n)                                       
      a2=(2.*a1+eps*a0)/p2                                              
      d2=(2.*d1+eps*d0+q1*a2)/p2                                        
!                                                                       
!     increment fp and sp                                               
      fnplp1=fnplp1+1.                                                  
      fp=fp+a2*fnplp1*rho**(ll+2)                                       
      sp=sp+d2*fnplp1*rho**(ll+2)                                       
!                                                                       
!     increment f and s                                                 
      f=f+a2*rho**(ll+3)                                                
      s=s+d2*rho**(ll+3)                                                
!                                                                       
!     calculate wronskian                                               
      w1=w2                                                             
      dw1=dw2                                                           
      w2=f*(clp*f+sp)-fp*s                                              
      dw2=dabs(w2-w1)                                                   
!d    for single precision replace dabs by abs                          
!                                                                       
!     convergence test                                                  
      k=n+1                                                             
      if(dw1.gt.acc2)goto 30                                            
      if(dw2.gt.acc2)goto 30                                            
      goto 50                                                           
!                                                                       
!     new a0,a1,do,d1                                                   
   30 a0=a1*rho                                                         
      a1=a2*rho                                                         
      d0=d1*rho                                                         
      d1=d2*rho                                                         
!                                                                       
   40 continue                                                          
!                                                                       
!  not converged                                                        
!                                                                       
      ierr=2                                                            
      actacc=dabs(0.25*w2-1.)                                           
!d    for single precision replace dabs by abs                          
      goto 60                                                           
!                                                                       
!  converged                                                            
!                                                                       
   50 actacc=dabs(0.25*w2-1.)         
!d    for single precision replace dabs by abs                          
      if(actacc.gt.acc)ierr=1                                           
!                                                                       
!  complete calculation of g and gp                                     
!                                                                       
   60 g=(s+cl*f)*r2pi                                                   
      gp=(sp+cl*fp+clp*f)*r2pi                                          
!                                                                       
      return                                                            
      end       subroutine
	subroutine seaton(l,eryd,r,zion,f,fp,g,gp)                             
	implicit real*8(a-h,o-z)                                               
! -- assume for now that zion=1                                         
	rho=zion*r                                                             
	epsr=eryd/(zion**2)                                                    
	acc=1.d-11                                                             
	call fogo(f,fp,g,gp,epsr,l,rho,ww,acc,actacc)                          
	s2=dsqrt(2.d0)                                                         
	f=f*s2                                                                 
	fp=fp*s2                                                               
	g=g*s2                                                                 
	gp=gp*s2                                                               
!hg -- add the following lines on 1-8-91:                               
                                                                        
	factor = dsqrt(zion)                                                   
	f=f/factor                                                             
	g=g/factor                                                             
	fp=fp*factor                                                           
	gp=gp*factor                                                           
                                                                        
!hg -- end additions                                                    
                                                                        
	return                                                                 
	end                       subroutine
!ccc                                                                    
      SUBROUTINE FOGO(F0,FP0,G0,GP0,EPSR,LL,RO,WW,ACC,ACTACC)           
!                                                                       
!   INPUT EPSR (E RYDBERG) LL (L) ACC=1.E-11                            
!    OUTPUT ENERGY NORMALISED FUNCTIONS AND DERIVATIVES (RYDBERG)       
!                                                                       
!                                                                       
        IMPLICIT REAL*8 (A-H,O-Z)                                       
!hg                                                                     
	common/agqdt/aqdt,gqdt,fr0,fpr0,gr0,gpr0                               
!hg                                                                     
      DATA R2PI/.159154943091895336D0/                                
       DATA PI/ 3.1415926535897932D0/                                
!     PRINT 1013                                                        
1013  FORMAT(//////,121(1H*) ,/)                                        
                                        
       A=1.D0                                                           
       IF(LL.EQ.0) GO TO 50                                             
       TESEN=1.D-10                                                     
       IF(DABS(EPSR).LT.TESEN)  GO TO 50                                
!                                                                       
!                                                                       
! CALCUL DE A= PRODUIT( 1. + EPS * P2 )                                 
       DO 20 LP=1,LL                                                    
       B=1.D0 + LP*LP*EPSR                                              
       A=A*B                                                            
!     PRINT *,' A B ',A,B                                               
20     CONTINUE                                                         
50     CONTINUE                                                         
!     PRINT 1010 ,EPSR,LL,RO,A                                          
1010  FORMAT(/////,' ENERGIE  ',D25.16,' RYD  L= ',I3,' RHO = ',D25.16,/,' A = ',D25.16 ,//)
!                                                                       
!                                                                       
!                                                                       
! CALCUL DE G ET DE B                                                   
! *******************                                                   
      IF(DABS(EPSR).LT.TESEN )  GO TO 8000                              
      IF(EPSR.GE.TESEN) GO TO 5000                                      
!                                                                       
! CAS ENERGIE NEGATIVE                                                  
      PNU=DSQRT(-1.D0/EPSR)                                             
      PL1=DFLOAT(LL)                                                    
!hg                                                                     
	pl1=0.d0                                                               
!hg                                                                     
      IF(PNU.GT.PL1)  GO TO 61                                          
      PRINT 1062,LL,EPSR                                                
 1062 FORMAT(//////,131(1H*),/,'****** CAS ETUDIE IMPOSSIBLE  .  L = ',I4,' ENERGIE =  ',D25.13,' RYDBERG ',/, 131(1H*) )                
      STOP                                                              
 61   CONTINUE                                                          
      X=PNU+LL+1.D0                                                     
      PSI=D1IGAM(X)                                                     
!     PRINT *,' X PSI ',X,PSI                                           
      X=PNU-LL                                                          
      PSI=PSI +D1IGAM(X) -2.D0*DLOG(PNU)                                
      X2=DLOG(PNU)                                                      
      X1=D1IGAM(X)                                                      
!     PRINT *,' X PSI LOG ',X,X1,X2                                     
      G=A*R2PI*PSI                                                      
!     PRINT 1020,G,PNU,PSI                                              
1020  FORMAT('  G = ',D25.16,' KAPA ',D25.16,' PSI',D25.16)             
      B=A                                                               
!     PRINT *, ' B = ',B                                                
!hg                                                                     
	aqdt=b                                                                 
	gqdt=g                                                                 
!hg                                                                     
      GO TO 4000                                                        
5000    CONTINUE                                                        
!                                                                       
! CAS ENERGIE POSITIVE                                                  
      PGA=DSQRT(1.D0/EPSR)                                              
      SUM=0.D0                                                          
! CALCUL DE LA DEPENDANCE EN L                                          
      IF(LL.EQ.0) GO TO 90                                              
      DO 80 LP=1,LL                                                     
      S1=LP*EPSR                                                        
      SUM=SUM + S1/(1.D0 + LP*S1)                                       
!     PRINT *,'  SUM S1  ',SUM,S1                                       
 80   CONTINUE                                                          
90    CONTINUE                                                          
      IF(EPSR>0.05D0 ) GO TO 92                                   
! ENERGIE < 0.05 RYD . FORMULE 2                                        
      X1=GAMI(PGA,2)                                                    
      GO TO 99                                                          
 92   IF(EPSR>0.8D0 ) GO TO 96                                    
! 0.05 RYD < ENERGIE < 0.8 RYD . FORMULE 1                              
      X1=GAMI(PGA,1)                                                    
      GO TO 99                                                          
! O.8 RYD < ENERGIE . FORMULE 3                                         
96    X1=GAMI(PGA,3)                                                    
99    CONTINUE                                                          
      G=(SUM + X1)* A * 2.D0 * R2PI                                     
!     PRINT 1040,G,PGA,SUM,X1                                           
 1040 FORMAT('  G = ',D25.16,' GAMA = ',D25.16,' SUM = ',D25.16,' X = ',D25.16)                                            
      B=A/(1.D0 - DEXP( -PGA/R2PI))                                     
!     PRINT *, '  B = ',B                                               
!hg                                                                     
	aqdt=b                                                                 
	gqdt=g                                                                 
!                                                                      
      GO TO 4000                                                        
8000   CONTINUE                                                         
!                                                                       
! CAS ENERGIE=0                                                         
      G=0.D0                                                            
      B=1.D0                                                            
!     PRINT 1050,G                                                      
1050  FORMAT('  G = ',D25.16)                                           
!     PRINT *,'  B = ',B                                                
4000  CONTINUE                                                          
!                                                                       
! CALCUL DE F(R0) ET DE G(R0) ANALYTIQUES EN ENERGIE                    
! **************************************************                    
!********WRONSKIEN = 2 / PI                                            
!         K TERMES DANS LE CALCUL                                       
!         ACTAC! PRECISION RELATIVE DU WRONSKIEN                        
!                                                                       
      CALL COULFG(LL,EPSR,RO,ACC,FR0,FPR0,GR0,GPR0,KK,IERR,ACTACC)     
!     PRINT *                                                           
!     PRINT *                                                           
!     PRINT *,'  LL  EPSR   RO ',LL,EPSR,RO                             
!     PRINT *,' VALEURS DE F  FP  G  GP '                               
!     PRINT *,FR0,FPR0,GR0,GPR0                                         
!     PRINT *,' K  IER WRONSKIEN ',KK,IERR,ACTAC!                       
!                                                                       
! CALCUL DE F0 ET G0 NORMEES EN RYDBERG                                 
! *************************************                                 
!                                                                       
!                                                                      
	b=dabs(b)                                                              
!                                                                      
      HR0 =GR0 + G * FR0                                                
      HPR0 = GPR0 + G * FPR0                                            
      B = DSQRT ( B / 2.D0  )                                           
      F0 = FR0 * B                                                      
      FP0 = FPR0 * B                                                    
      G0 = HR0 / ( 2.D0 * B )                                           
       GP0 = HPR0 /( 2.D0 * B )                                         
      WW=(F0*GP0-FP0*G0)*PI-1.D0                                        
!       PRINT*,' FONCTION NORMEES EN RYDBERG'                           
!     PRINT *,F0,FP0,G0,GP0                                             
!     PRINT*,' 2EME WRONSKIEN ',WW                                      
!      PRINT *                                                          
      RETURN                                                            
      END            subroutine
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
      FUNCTION GAMI(X,N)                                                
      IMPLICIT REAL*8(A-H,O-Z)                                          
      DIMENSION B(7)                                                    
      DIMENSION C(20)                                                   
      DIMENSION CCOR(4)                                                 
      DATA CCOR/ .20200740065967752D0,.03692775269295320D0,.00834927738176102D0, .00200839282608212D0 /
      DATA B /8.333333333333333D-2,+8.333333333333333D-3, 3.968253968253968D-3,+4.166666666666667D-3,&
7.575757575757576D-3,+2.109279609279609D-2, 8.333333333333333D-2/
      DATA C /2.020569031595943D-1, 3.692775514336993D-2,8.349277381922827D-3, 2.008392826082214D-3,&
4.941886041194646D-4, 1.227133475784891D-4, 3.058823630702049D-5, 7.63719763789976D-6,&
1.90821271655394D-6, 4.7693298678781D-7,1.1921992596531D-7, 2.980350351465D-8,7.45071178984D-9,&
1.86265972351D-9,4.6566290650D-10, 1.1641550173D-10,2.910385044D-11,7.27595984D-12,1.81898965D-12,4.5474738D-13 /
      DATA EUL/.5772156649015329D0/                                
      DATA TES/1.D-10/                                                  
! CALCUL DE REEL(PSI(I*GAM)) - LOG(GAM)                                 
!        N=1 SOMME 1/N(N2+GAM2)                                         
!        N=2 DEVELOPPEMENT ASYMPTOTIQUE                                 
!        N=3 CAS GAM INF A 2                                            
      IF(N.GE.4) GO TO 8000                                             
      A=DABS(X)                                                         
      AA=A*A                                                            
      IF(N.EQ.2) GO TO 50                                               
      IF(N.EQ.3) GO TO 100                                              
      S1=1.D0/(1.D0 + AA)                                               
!     PRINT *, '  S1  ',S1                                              
      DO 10 IK=2,100                                                    
      DS=AA + IK*IK                                                     
      DS=1.D0/(DS*IK)                                                   
      S1=S1+DS                                                          
10    CONTINUE                                                          
!     PRINT *,' X S1 DS TES ',X,S1,DS,TES                               
!     PRINT *,' ***** DEVELOPPEMENT NON CONVERGE'                       
      GO TO 30                                                          
20    S1=S1+DS                                                          
30    CONTINUE                                                          
      RAP=DS/S1                                                         
      S2=DLOG(A)                                                        
!     PRINT *,' S1,LOG',S1,S2                                           
      GAMI=S1*AA - EUL -S2                                              
       RAP2=DABS(GAMI/(EUL+S2))                                         
!     PRINT *,' CONVERG ANULATION ',RAP,RAP2                            
!     PRINT 1020,N,X,GAMI                                               
1020  FORMAT(' OPTION N= ',I3,' X = ',D25.16,' GAMI = ',D25.16,/)      
      C1=0.D0                                                           
      X1=-X*X                                                           
      X2=-1.D0                                                          
      DO 5010 IC=1,4                                                    
      X2=X2*X1                                                          
      C2=X2*(C(IC)  - CCOR(IC))                                         
      C1=C1 + C2                                                        
      GAMI= GAMI + C2                                                   
!     PRINT *,' C1 C2 GAMI ',C1,C2,GAMI                                 
5010  CONTINUE                                                          
!     PRINT 1020,N,X,GAMI                                               
      RETURN                                                            
50    AA=1.D0/AA                                                        
      A1=1.D0                                                           
      S1=0.D0                                                           
      DO 60 IK=1,7                                                      
      A1=A1*AA                                                          
      DS=B(IK)*A1                                                       
      S1= S1 + DS                                                       
      RAP=DS/S1                                                         
!     PRINT *,' S1 DS RAP ',S1,DS,RAP                                   
60    CONTINUE                                                          
      GAMI=S1                                                           
!     PRINT 1020,N,X,GAMI                                               
      RETURN                                                            
100   A1=-1.D0                                                          
      S1=0.D0                                                           
      DO 110 IK=1,20                                                    
      A1=-A1*AA                                                         
      DS=C(IK)*A1                                                       
      S1=S1 + DS                                                        
      RAP=DS/S1                                                         
!     PRINT *,' S1 DS RAP ', S1,DS,RAP                                  
110   CONTINUE                                                          
      S2=DLOG(A)                                                        
      S3=1.D0/(1.D0 + AA)                                               
      GAMI = 1.D0 -EUL -S2+S1 -S3                                       
!     PRINT *,' S1 LOG S3 ',S1,S2,S3                                    
      RAP=DABS(S1/GAMI)                                                 
      RAP=1.D0/RAP                                                      
!     PRINT *,'  ANULATION ',RAP                                        
!     PRINT 1020,N,X,GAMI                                               
      RETURN                                                            
8000   PRINT 8010,N                                                     
      GAMI=0.D0                                                         
8010      FORMAT('  FORMULE NON PROGRAMMEE N',' ***********************',///)
      RETURN                                                            
       END          function
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
      FUNCTION D1IGAM(X)                                                
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)                              
      DIMENSION B(6)                                                    
      DATA PI/ 3.141592653589793D0/                                  
      DATA B / +8.333333333333333D-2, -8.333333333333333D-3,   +3.968253968253968D-3,&
-4.166666666666667D-3,    +7.575757575757576D-3, -2.109279609279609D-2/
! CALCUL DE PSI(X)  X REEL                                              
! PASSAGE X POSITIF PLUS GRAND QUE 15                                   
! FORMULE ASYMPTOTIQUE A 7 TERMES                                       
      A=DABS(X)                                                         
!     IF(DINT(X) + A)  ,4,                                              
      XZ=DINT(X)+A                                                      
      IF(XZ.EQ.0.D0) GO TO 4                                            
      V=A                                                               
      H=0.D0                                                            
      IF(A.GE.(15.D0)) GO TO 3                                          
      N=14 - INT(A)                                                     
      H=1.D0/V                                                          
!     IF(N)  ,2,                                                        
      IF(N.EQ.0) GO TO 2                                                
      DO 1 I=1,N                                                        
      V= V + 1.D0                                                       
1     H= H + 1.D0/V                                                     
2     V= V + 1.D0                                                       
3     R=1.D0/V**2                                                       
      D1IGAM=DLOG(V) -0.5D0/V -R*(B(1)+R*(B(2)+R*(B(3)+R* (B(4)+R*(B(5)+R*(B(6)+R*B(1))))))) - H                           
      IF(X>=(0.000D0) )RETURN    
      H=PI*A                                                            
      D1IGAM=D1IGAM + 1.D0/A + PI*DCOS(H)/DSIN(H)                       
      RETURN                                                            
! SORTIE ERREUR  :  X ENTIER NON POSITIF                                
4     PRINT 100,X                                                       
100   FORMAT(/////,131(1H*)//, '  *** DDIGAM   ARGUMENT ENTIER NON NEGATIF =',D16.9,   ' ***** ' )                                                     
      D1IGAM=0.D0                                                       
      RETURN                                                            
      END        function
                                                                        
! ********************************************

      SUBROUTINE Cubic_spline(x,y,n,yp1,ypn,y2)

      implicit none
      INTEGER::n,i,k,NMAX=5000
      real*8::yp1,ypn,x(n),y(n),y2(n),ONE=1.0d0,TWO=2.0d0,SIX=6.0d0
      real*8::XB,XF,YB,YF,P,qn,sig,un
      real*8,allocatable :: u(:)

      allocate (u(1:n))

      XF=X(3)-X(2)
      XB=X(2)-X(1)
      YF=Y(3)-Y(2)
      YB=Y(2)-Y(1)
      IF (YP1.GT..99E30) THEN
        P=TWO*XF+XB
        Y2(2)=(XB-XF)/P
        U(2)=SIX*(YF-YB*XF/XB)/((XF+XB)*P)
      ELSE
        Y2(1)=-0.5
        U(1)=(3./XB)*(YB/XB-YP1)
        P=TWO*(XF+XB)+XB*Y2(1)
        Y2(2)=-XF/P
        U(2)=(SIX*(YF/XF-YB/XB)-XB*U(1))/P
      ENDIF
      DO 11 I=3,N-1
        XF=X(I+1)-X(I)
        XB=X(I)-X(I-1)
        YF=Y(I+1)-Y(I)
        YB=Y(I)-Y(I-1)
        P=TWO*(XF+XB)+XB*Y2(I-1)
        Y2(I)=-XF/P
        U(I)=(SIX*(YF/XF-YB/XB)-XB*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        P=TWO*XF+XB
        Y2(N)=(SIX*(YF-YB*XF/XB)/(XF+XB)-P*U(N-1))  /(XF-XB+P*Y2(N-1))
      ELSE
        Y2(N)=(SIX*(YPN-YF/XF)/XF-U(N-1))/(TWO+Y2(N-1))
      ENDIF
      DO 12 K=N-1,2,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      IF (YP1.GT..99E30) THEN
        Y2(1)=Y2(2)-(X(2)-X(1))*(Y2(3)-Y2(2))/(X(3)-X(2))
        YP1=(Y(2)-Y(1))/(X(2)-X(1))  -(X(2)-X(1))*(Y2(1)/3D0+Y2(2)/6D0)
      ELSE
        Y2(1)=Y2(1)*Y2(2)+U(1)
      END IF
      IF (YPN.GT..99E30) THEN
        YPN=YF/XF+XF*(Y2(N-1)/6D0+Y2(N)/3D0)
      END IF

      deallocate (u)

      return
      END SUBROUTINE Cubic_spline

! ********************************************
!
!      CUBIC SPLINE FITTING PROGRAM
!
!  GIVEN N TABULATED FUNCTIONS XA, YA AND 2ND DERIVATIVES Y2A,
!  THIS SUBROUTINE RETURNS AT X
!              FUNCTION VALUE FOR KIND=0
!              1ST DERIVATIVE FOR KIND=1
!              2nd DERIVATIVE FOR KIND=2
!
! ********************************************

      SUBROUTINE Cubic_splint(xa,ya,y2a,n,x,y,kind)

      implicit none
      INTEGER n,k,khi,klo,kind
      real*8::x,y,xa(n),y2a(n),ya(n),a,b,h

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
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h

      if (KIND.EQ.0) THEN
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      ELSE IF (KIND.EQ.1) THEN
         Y=(YA(KHI)-YA(KLO))/H+((1.-3.*A*A)*Y2A(KLO)+(3.*B*B-1.)*Y2A(KHI))*H/6.
      ELSE
         Y= A*Y2A(KLO)+B*Y2A(KHI)
      END IF

      return
      END SUBROUTINE Cubic_splint

! ********************************************

!  there follow some (elaborate) routines from Numerical Recipes (TM)
!  for computing spherical Bessel functions

!   ********************************************

      SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)

      implicit none
      INTEGER n
      REAL*8 sj,sjp,sy,syp,x,factor,order,rj,rjp,ry,ryp,RTPIO2
!U    USES bessjy

      RTPIO2=sqrt(asin(1.d0))
      if(n < 0.d0.or.x<=0.d0)pause 'bad arguments in sphbes'
      order=n+0.5d0
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.d0*x)
      syp=factor*ryp-sy/(2.d0*x)
      return
      END SUBROUTINE sphbes

!   ********************************************

      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)

      implicit none
      INTEGER::MAXIT=10000,i,isign,l,nl
      REAL*8::rj,rjp,ry,ryp,x,xnu,PI=3.141592653589793d0
      REAL*8::EPS=1.d-16,FPMIN=1.d-30,XMIN=2.d0,q,r,sum1
!U    USES beschb
      REAL*8 a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e, &
&f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,   &
&rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum, &
&temp,w,x2,xi,xi2,xmu,xmu2,p,pimu,pimu2
      if(x<=0.d0.or.xnu < 0.d0) pause 'bad arguments in bessjy'
      if(x < XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      dr=xnu*xi
      if(dr < FPMIN)dr=FPMIN
      b=xi2*xnu
      d=0.d0
      c=dr
      do i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d) < FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c) < FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        dr=del*dr
        if(d < 0.d0)isign=-isign
        if(abs(del-1.d0) < EPS)goto 1
      end do
      print*,x
      pause 'x too large in bessjy; try asymptotic expansion'
1     continue
      rjl=isign*FPMIN
      rjpl=dr*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
      end do
      if(rjl==0.d0)rjl=EPS
      f=rjpl/rjl
      if(x < XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu) < EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e) < EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2) < EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del) < (1.d0+abs(sum))*EPS)goto 2
        end do
        pause 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di) < FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci) < FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli) < EPS)goto 3
        end do
        pause 'cf2 failed in bessjy'
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
      end do
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
      END SUBROUTINE bessjy

!   ********************************************

      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      INTEGER::NUSE1=7,NUSE2=8
      REAL*8 gam1,gam2,gammi,gampl,x,xx,c1(7),c2(8),chebev
!U    USES chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END SUBROUTINE beschb

!   ********************************************

      Real*8 FUNCTION chebev(a,b,c,m,x)
      INTEGER m,j
      REAL*8 a,b,x,c(m),d,dd,sv,y,y2
      if ((x-a)*(x-b) > 0.d0) pause 'x not in range in chebev'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
      end do
      chebev=y*d-dd+0.5d0*c(1)
      return
      END FUNCTION chebev
