      subroutine Hyperfine2b(MFab,Bfield,NumSym2bStates,Thresholds,VQcal,indexsav,alphasav,alphasym2bsav)
      implicit double precision(a-h,o-z)

      integer j,k,l,ni,nj,MFab
      integer NumHypfStates,NumSym2bStates

      double precision i,s,Bfield,Gauss
      double precision xi,gE,muE,gN,muN
      double precision alpha(100,6),alphaSym2b(100,6)
      double precision ::alphasav(100,6),alphaSym2bsav(100,6)

      double precision Thresholds(50)
      double precision VQcal(0:1,50,50)

      integer :: indexsav(20)
      integer, allocatable:: Index(:)

      double precision, allocatable:: Energies2b(:),EnergiesTemp(:)
      double precision, allocatable:: Energies(:),U(:,:)
      double precision, allocatable:: QcalTemp(:,:,:)
      double precision clebsh,Sym2b,Qcal,QcalSym2b
      external clebsh,Sym2b,Qcal,QcalSym2b

      Gauss = (4.254382547308656e-10)
      if (Bfield.eq.0.d0) Bfield = 1.d-4*Gauss

      call Parameters_87Rb(i,s,xi,gE,muE,gN,muN,NumHypfStates)
!call cpu_time(time1)
      call CalcAlpha(i,s,MFab,xi,muN,gN,muE,gE,NumSym2bStates,alpha,alphaSym2b)
!call cpu_time(time2)
!write(*,*) "time cost in calcalpha=",time2-time1
alphasav=alpha
alphaSym2bsav=alphaSym2b


!     Allocating Energies and U
      allocate(Energies2b(NumSym2bStates),EnergiesTemp(NumSym2bStates))
      allocate(Index(NumSym2bStates))
      allocate(Energies(NumHypfStates))
      allocate(U(NumHypfStates,NumHypfStates))
!call cpu_time(time1)
      call CalcHypfStates(Bfield,i,s,xi,muN,gN,muE,gE,Energies,U)
!call cpu_time(time2)
!write(*,*) "time cost in calcHypfstates=",time2-time1
!      write(200,11)Bfield/Gauss,(Energies(j),j = 1,NumHypfStates)

!     Hyperfine energies
      do j = 1,NumSym2bStates
      Energies2b(j) = Energies(nint(alphaSym2b(j,1)))+Energies(nint(alphaSym2b(j,2)))
      enddo

      EnergiesTemp = Energies2b
      do j = 1,NumSym2bStates
         Index(j) = MinLoc(EnergiesTemp,dim=1,mask=EnergiesTemp.ne.0.d0)
         indexsav(j)=index(j)
         EnergiesTemp(Index(j)) = 0.d0
      enddo
!      write(*,*) 'index',index(1:10)

!     write(100,11)Bfield/Gauss,(Energies2b(Index(j)),j=1,NumSym2bStates)
!
!     label1=(mfab+5)*100+2
!     label2=(mfab+5)*100+1
!     do j = 1,NumSym2bStates
!     write(label1,111)j,nint(alphaSym2b(Index(j),3)),nint(alphaSym2b(Index(j),4)),&
!                    nint(alphaSym2b(Index(j),5)),nint(alphaSym2b(Index(j),6))
!     enddo
!     write(label1,*)'-------------'

      do j = 1,NumSym2bStates
         Thresholds(j) = Energies2b(Index(j))
      enddo
!      write(label2,11)Bfield/Gauss,(Thresholds(j),j = 1,NumSym2bStates)
!call cpu_time(time1)      
!     QcalSym2b
      allocate(QcalTemp(0:1,NumSym2bStates,NumSym2bStates))
!     write(1000,*)Bfield/Gauss
!     write(2000,*)Bfield/Gauss
      do ni = 1,NumSym2bStates
      do nj = ni,NumSym2bStates   
      QcalTemp(0,ni,nj) = QcalSym2b(0,MFab,i,s,Index(ni),Index(nj),NumHypfStates,alpha,NumSym2bStates,alphaSym2b,U)
      QcalTemp(0,nj,ni) = QcalTemp(0,ni,nj)
      QcalTemp(1,ni,nj) = QcalSym2b(1,MFab,i,s,Index(ni),Index(nj),NumHypfStates,alpha,NumSym2bStates,alphaSym2b,U)
      QcalTemp(1,nj,ni) = QcalTemp(1,ni,nj)
      enddo
!     write(1000,11)(QcalTemp(0,ni,nj),nj=1,NumSym2bStates)
!     write(2000,11)(QcalTemp(1,ni,nj),nj=1,NumSym2bStates)
      enddo
!call cpu_time(time2)
!write(*,*)"time in Qcalsym2b=",time2-time1

      do ni = 1,NumSym2bStates
      do nj = 1,NumSym2bStates   
         VQcal(0,ni,nj) = QcalTemp(0,ni,nj)
         VQcal(1,ni,nj) = QcalTemp(1,ni,nj)
      enddo
      enddo

        do ni=1,NumSym2bStates
        do nj=ni+1,NumSym2bStates
        if (abs(VQcal(1,ni,nj)+VQcal(0,ni,nj))>1d-10) write(*,*) "VQCAL(0,0)+VQCAL(0,1)=/=0","i=",ni,"j=",nj,VQcal(1,ni,nj)+VQcal(0,ni,nj)
        end do
        end do
        do ni=1,NumSym2bStates
        if (abs(VQcal(1,ni,ni)+VQcal(0,ni,ni)-1d0)>1d-10) write(*,*) "VQCAL(0,0)+VQCAL(0,1)=/=0","i=",ni,VQcal(1,ni,ni)+VQcal(0,ni,ni)
        end do

      deallocate(QcalTemp)
      deallocate(Energies2b,EnergiesTemp)
      deallocate(Index,Energies,U)


   11 format(200e20.12)
  111 format(10X,I2,': |',I1,','I2,'>+|',I1,',',I2,'>')

      return
      end

!     ------------------------------------
!     Determining Atomic Hyperfine States
!     and Energies
!     ------------------------------------

!     using |i mi s ms > base
      subroutine CalcHypfStates(B,i,s,xi,muN,gN,muE,gE,Energies,U)
      implicit none
      integer NumHypfStates,icount,jcount,j,k,B0count
      integer mi,ms,mip,msp,mf,ni,nj
      double precision B,xi,muN,gN,muE,gE,BB
      double precision i,s,f,HypfZaux,xNorm,xSum
      double precision Energies((nint(2.d0*i)+1)*(nint(2.d0*s)+1))
      double precision U((nint(2.d0*i)+1)*(nint(2.d0*s)+1),(nint(2.d0*i)+1)*(nint(2.d0*s)+1))
      double precision EnergiesB((nint(2.d0*i)+1)*(nint(2.d0*s)+1))
      double precision UB((nint(2.d0*i)+1)*(nint(2.d0*s)+1),(nint(2.d0*i)+1)*(nint(2.d0*s)+1))
      double precision, allocatable ::HypfZ(:,:),Identity(:,:)
      double precision, allocatable ::evalues(:),evectors(:,:)

      integer lwork,info
      double precision work(3*(nint(2.d0*i)+1)*(nint(2.d0*s)+1))

      integer NonZeroTerms
      integer NonZeroLoc(100)

      double precision clebsh
      external clebsh

      B0count = 0
      U = 0.d0
      Energies = 0

!     Matrix Element for the Hamiltonian Matrix
      NumHypfStates = (nint(2.d0*i)+1)*(nint(2.d0*s)+1)
      allocate(HypfZ(NumHypfStates,NumHypfStates))
      HypfZ = 0.d0

      icount = 0
      do mi = -nint(2.d0*i),nint(2.d0*i),2
      do ms = -nint(2.d0*s),nint(2.d0*s),2
         icount = icount+1
         jcount = 0
         do mip = -nint(2.d0*i),nint(2.d0*i),2
         do msp = -nint(2.d0*s),nint(2.d0*s),2
            jcount = jcount+1
            HypfZaux = 0.d0
            f = dabs(i-s)
            do mf = -nint(2.d0*f),nint(2.d0*f),2
               HypfZaux = HypfZaux+clebsh(nint(2.d0*i),nint(2.d0*s),nint(2.d0*f),mi,ms,mf)*&
                      xi/2.d0*(f*(f+1.d0))*clebsh(nint(2.d0*i),nint(2.d0*s),nint(2.d0*f),mip,msp,mf)
            enddo
            f = i+s
            do mf = -nint(2.d0*f),nint(2.d0*f),2
               HypfZaux = HypfZaux+clebsh(nint(2.d0*i),nint(2.d0*s),nint(2.d0*f),mi,ms,mf)*&
                      xi/2.d0*(f*(f+1.d0))*clebsh(nint(2.d0*i),nint(2.d0*s),nint(2.d0*f),mip,msp,mf)
            enddo
            if (mi.eq.mip.and.ms.eq.msp) then
               HypfZaux = HypfZaux-(xi/2.d0)*(i*(i+1)+s*(s+1))
               HypfZaux = HypfZaux-B*(muN*gN*dble(mi)/2.d0+muE*gE*dble(ms)/2.d0)
            endif
         HypfZ(icount,jcount) = HypfZaux   
         enddo
         enddo
      enddo   
      enddo   

!     Identity Matrix
      allocate(Identity(NumHypfStates,NumHypfStates))
      Identity = 0.d0
      do j = 1,NumHypfStates
      do k = 1,NumHypfStates
         if (j.eq.k) then
         Identity(j,k) = 1.d0
         endif
      enddo
      enddo

!     Diagonalizing Hamiltonian
      allocate(evalues(NumHypfStates),evectors(NumHypfStates,NumHypfStates))
      
      lwork = 3*NumHypfStates
      call dsyev('V','U',NumHypfStates,HypfZ,NumHypfStates,evalues,work,lwork,info)
      evectors = HypfZ

!     Normalization
      do j = 1,NumHypfStates
      xnorm = 0.d0
         do k = 1,NumHypfStates
         xnorm = xnorm+evectors(k,j)**2
         enddo
         do k = 1,NumHypfStates
         evectors(k,j) = evectors(k,j)/dsqrt(xnorm)
         enddo
      enddo   
      do j = 1,NumHypfStates
      Energies(j) = evalues(j)
         do k = 1,NumHypfStates
         U(k,j) = evectors(k,j)
         enddo
      enddo

      deallocate(evalues,evectors,HypfZ,Identity)

!!    Imposing a Phase for the eigenvectors
!     do nj = 1,NumHypfStates
!     icount = 0
!        do ni = 1,NumHypfStates
!        if (dabs(U(ni,nj)).ge.1.d-12.and.icount.eq.0) then
!        icount = icount+1   
!        if(dsign(1.d0,U(ni,nj)).ne.1.d0) then
!          do k = ni,NumHypfStates
!             U(k,nj) = -U(k,nj)
!          enddo   
!        endif
!        endif
!        enddo
!     enddo
         
      return
      end

!     -------------------------------
!     Determining Qcal coefficients
!     -------------------------------

      double precision function Qcal(Sab,i,s,ialpha_a,ialpha_b,ialpha_ap,ialpha_bp,U)
      implicit none
      integer Sab,MSab
      integer Iab,MIab
      integer ialpha_a,ialpha_b,ialpha_ap,ialpha_bp
      integer mi_a,ms_a,mi_b,ms_b,mi_ap,ms_ap,mi_bp,ms_bp
      integer mis_a,mis_b,mis_ap,mis_bp
      double precision i,s,xSum
      double precision U((nint(2.d0*i)+1)*(nint(2.d0*s)+1),(nint(2.d0*i)+1)*(nint(2.d0*s)+1))
      double precision clebsh
      external clebsh

      xSum = 0.d0
      mis_a = 0
      do mi_a = -nint(2.d0*i),nint(2.d0*i),2
      do ms_a = -nint(2.d0*s),nint(2.d0*s),2
      mis_a = mis_a+1
      mis_b = 0
      do mi_b = -nint(2.d0*i),nint(2.d0*i),2
      do ms_b = -nint(2.d0*s),nint(2.d0*s),2
      mis_b = mis_b+1
         mis_ap = 0
         do mi_ap = -nint(2.d0*i),nint(2.d0*i),2
         do ms_ap = -nint(2.d0*s),nint(2.d0*s),2
         mis_ap = mis_ap+1
         mis_bp = 0
         do mi_bp = -nint(2.d0*i),nint(2.d0*i),2
         do ms_bp = -nint(2.d0*s),nint(2.d0*s),2
         mis_bp = mis_bp+1   
!        Option 1
            do MSab = -Sab,Sab
            if (mi_a.eq.mi_ap) then
            if (mi_b.eq.mi_bp) then
            if (U(mis_a,ialpha_a).ne.0.d0.and.U(mis_b,ialpha_b).ne.0.d0.and. &
               U(mis_ap,ialpha_ap).ne.0.d0.and.U(mis_bp,ialpha_bp).ne.0.d0) then
                xSum = xSum+U(mis_a,ialpha_a)*U(mis_b,ialpha_b)* &
                      U(mis_ap,ialpha_ap)*U(mis_bp,ialpha_bp)* &
                      clebsh(nint(2.d0*s),nint(2.d0*s),nint(2.d0*Sab),ms_a,ms_b,2*MSab)* &
                      clebsh(nint(2.d0*s),nint(2.d0*s),nint(2.d0*Sab),ms_ap,ms_bp,2*MSab)
            endif
            endif
            endif
            enddo
         go to 71
!        Option 2
            do MSab = -Sab,Sab   
            do Iab = 0,nint(2*i)
            do MIab = -Iab,Iab
            if ((ms_a+ms_b.eq.2*MSab).and.(ms_ap+ms_bp.eq.2*MSab)) then               
            if ((mi_a+mi_b.eq.2*MIab).and.(mi_ap+mi_bp.eq.2*MIab)) then
            if ((-1.d0)**(nint(2.d0*s)-Sab+nint(2.d0*i)-Iab+0).eq.1.d0) then
            if (MSab+MIab.eq.-4) then               
            if (U(mis_a,ialpha_a).ne.0.d0.and.U(mis_b,ialpha_b).ne.0.d0.and. &
               U(mis_ap,ialpha_ap).ne.0.d0.and.U(mis_bp,ialpha_bp).ne.0.d0) then
!                xSum = xSum+U(mis_a,ialpha_a)*U(mis_b,ialpha_b)*
!     .                 U(mis_ap,ialpha_ap)*U(mis_bp,ialpha_bp)*
!     .                 clebsh(nint(2.d0*s),nint(2.d0*s),nint(2.d0*Sab),ms_a,ms_b,2*MSab)*
!     .                 clebsh(nint(2.d0*s),nint(2.d0*s),nint(2.d0*Sab),ms_ap,ms_bp,2*MSab)*
!     .                 clebsh(nint(2.d0*i),nint(2.d0*i),nint(2.d0*Iab),mi_a,mi_b,2*MIab)*
!     .                 clebsh(nint(2.d0*i),nint(2.d0*i),nint(2.d0*Iab),mi_ap,mi_bp,2*MIab)
            endif 
            endif
            endif
            endif
            endif
            enddo
            enddo
            enddo
        71 continue
         enddo
         enddo
         enddo
         enddo
      enddo
      enddo
      enddo
      enddo
        
      Qcal = xSum

      return
      end

!     -----------------------------------
!     Determining QcalSym2b coefficients 
!     -----------------------------------

      double precision function QcalSym2b(Sab,MFab,i,s,ialpha2b,ialpha2bp,NumHypfStates,alpha,NumSymStates,alphaSym2b,U)
      implicit none
      integer Sab,MFab,NumHypfStates,NumSymStates,icount
      integer ialpha_a,ialpha_b,ialpha_ap,ialpha_bp
      integer iat_a,iat_b,iat_ap,iat_bp
      integer mi_a,ms_a,mi_b,ms_b,mi_ap,ms_ap,mi_bp,ms_bp
      integer mis_a,mis_b,mis_ap,mis_bp
      integer ialpha2b,ialpha2bp
      integer mfa,mfb,mfap,mfbp
      double precision alpha(100,6)
      double precision alphaSym2b(100,6)
      double precision i,s,xSum
      double precision U(NumHypfStates,NumHypfStates)
      double precision QcalPreCalc(0:1,NumHypfStates,NumHypfStates,NumHypfStates,NumHypfStates)
      double precision Qcal,Sym2b
      external Qcal,Sym2b

      ialpha_a = alphaSym2b(ialpha2b,1)
      ialpha_b = alphaSym2b(ialpha2b,2)
      ialpha_ap = alphaSym2b(ialpha2bp,1)
      ialpha_bp = alphaSym2b(ialpha2bp,2)

      xSum = 0.d0
      do iat_a = 1,NumHypfStates
      mfa = nint(alpha(iat_a,2))
      do iat_b = 1,NumHypfStates
      mfb = nint(alpha(iat_b,2))
         do iat_ap = 1,NumHypfStates
         mfap = nint(alpha(iat_ap,2))
         do iat_bp = 1,NumHypfStates
         mfbp = nint(alpha(iat_bp,2))
         if (mfa+mfb.eq.MFab) then
         if (mfap+mfbp.eq.MFab) then            
         if (Sym2b(iat_a,iat_b,ialpha_a,ialpha_b).ne.0.d0) then
         if (Sym2b(iat_ap,iat_bp,ialpha_ap,ialpha_bp).ne.0.d0) then
         xSum = xSum+Sym2b(iat_a,iat_b,ialpha_a,ialpha_b)* &
                    Sym2b(iat_ap,iat_bp,ialpha_ap,ialpha_bp)* &
               Qcal(Sab,i,s,iat_a,iat_b,iat_ap,iat_bp,U)

         endif
         endif


         endif
         endif
         enddo
         enddo
      enddo
      enddo

      QcalSym2b = xSum

      return
      end

!     -----------------------------------------
!     Calculating Symmetryzed 2b states
!     -----------------------------------------

      subroutine CalcAlpha(i,s,MFab,xi,muN,gN,muE,gE,NumSym2bStates,alpha,alphaSym2b)
      implicit none

      integer j,mFab,NumHypfStates,NumSym2bStates
      integer mf,mi,ms
      integer mip,msp
      integer icount,jcount
      integer ialpha_a,ialpha_b,ialpha_c,ialpha_ap,ialpha_bp,ialpha_cp

      double precision alpha(100,6)
      double precision alphaSym2b(100,6)
      double precision f,i,s
      double precision mfa,mfb,mfc
      double precision xi,gE,muE,gN,muN,xSum
      double precision E0(1000)
      double precision, allocatable:: Energies(:),U(:,:)
      double precision Sym2b
      external Sym2b

      NumHypfStates = (nint(2.d0*i)+1)*(nint(2.d0*s)+1)

      allocate(Energies(NumHypfStates))
      allocate(U(NumHypfStates,NumHypfStates))

!     Defining |f,mf> labels 
      call CalcHypfStates(1.d-12,i,s,xi,muN,gN,muE,gE,Energies,U)
      do j = 1,NumHypfStates
      E0(j) = Energies(j)
      enddo
      do j = 1,NumHypfStates
      icount = 0
      do mi = -nint(2.d0*i),nint(2.d0*i),2
      do ms = -nint(2.d0*s),nint(2.d0*s),2
         icount = icount+1            
         if (dabs(U(icount,j)).gt.1.d-12) then
         f = -1.d0/2.d0+dsqrt(8.d0*E0(j)+xi+4.d0*i*(i+1.d0)*xi+4.d0*s*(s+1.d0)*xi)/2.d0/dsqrt(xi)   
         alpha(j,1) = nint(f)
         alpha(j,2) = dble(mi+ms)/2.d0
         alpha(j,3) = dble(mi)/2.d0
         alpha(j,4) = dble(ms)/2.d0
         endif   
      enddo
      enddo
      enddo


      icount = 0
      do ialpha_a = 1,NumHypfStates
      do ialpha_b = ialpha_a,NumHypfStates
         xSum = 0.d0
         do ialpha_ap = 1,NumHypfStates
         do ialpha_bp = 1,NumHypfStates
         if (nint(alpha(ialpha_a,2)+alpha(ialpha_b,2)).eq.MFab) then
         if (nint(alpha(ialpha_ap,2)+alpha(ialpha_bp,2)).eq.MFab) then         
         xSum = xSum+Sym2b(ialpha_a,ialpha_b,ialpha_ap,ialpha_bp)**2
         endif
         endif
         enddo
         enddo
         if (xSum.ne.0) then
            icount = icount+1
            alphaSym2b(icount,1) = dble(ialpha_a) 
            alphaSym2b(icount,2) = dble(ialpha_b)
            alphaSym2b(icount,3) = dble(alpha(ialpha_a,1)) !fa
            alphaSym2b(icount,4) = dble(alpha(ialpha_a,2)) !mfa
            alphaSym2b(icount,5) = dble(alpha(ialpha_b,1)) !fb
            alphaSym2b(icount,6) = dble(alpha(ialpha_b,2)) !mfb
         endif   
      enddo
      enddo

      NumSym2bStates = icount

      deallocate(Energies,U)

      return
      end


!     -------------------------------
!     Diagonalization Subroutine
!     -------------------------------

      subroutine Mydggev(N,G,LDG,L,LDL,eval,evec)
      implicit none
      integer LDG,N,LDL,info
      double precision G(LDG,N),L(LDL,N),eval(N),evec(N,N)
      double precision, allocatable :: alphar(:),alphai(:),beta(:),work(:),VL(:,:)
      integer lwork,i,im,in
  
      allocate(alphar(N),alphai(N),beta(N))
            
      info = 0
      lwork = -1
      allocate(work(1))
      call dggev('N','V',N,G,LDG,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info)
!     do im = 1,N
!        alphar(im)=0.0d0
!        alphai(im)=0.0d0
!        beta(im)=0.0d0
!        eval(im)=0.0d0
!        do in = 1,N
!           evec(im,in)=0.0d0
!        enddo
!     enddo
      
      lwork=work(1)
      deallocate(work)
      allocate(work(lwork))
      call dggev('N','V',N,G,LDG,L,LDL,alphar,alphai,beta,VL,1,evec,N,work,lwork,info)

      do i = 1, N
         if (abs(alphai(i)).ge.1e-14) then
            print*, '#eigenvalue may be complex! alphai(',i,')=',alphai(i)
         endif
         if(beta(i).ne.0.0d0) then
            eval(i) = -alphar(i)/beta(i)
         endif
!         if(abs(alphar(i)).ge.1e-14) then
!         print*, alphar(i), alphai(i), beta(i)
!         endif
      enddo
      call deigsrt(eval,evec,N,N)
      do i=1,N
         eval(i)=-eval(i)
      end do
      deallocate(alphar,alphai,beta)
      deallocate(work)

      return
      end 

      subroutine deigsrt(d,v,n,np)
      implicit none
      integer n,np
      real*8 d(np),v(np,np)
      integer i,j,k
      real*8 p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue

      return
      end

      double precision function deltafun(a,b)
      implicit none
      integer a,b
      if (a.eq.b) then
      deltafun = 1.d0
      else
      deltafun = 0.d0
      endif
      return
      end

      double precision function Sym2b(a,b,ap,bp)
      implicit none
      integer a,b,ap,bp
      double precision deltafun
      external deltafun
      Sym2b = dsqrt((2.d0-deltafun(a,b))/4.d0)* &
              (deltafun(a,ap)*deltafun(b,bp)+deltafun(a,bp)*deltafun(b,ap))
      return
      end

!     -------------------------------
!     Hyperfine Parameters for 85Rb
!     -------------------------------

      subroutine Parameters_85Rb(i,s,xi,gE,muE,gN,muN,NumHypfStates)
      implicit none

      integer NumHypfStates
      double precision i,s
      double precision Hz2Au,Pi,Gauss
      double precision Delta,Lambda,xi
      double precision gE,muE,gN,muN
      
      Hz2Au = 1.519829846006d-16
      Pi = dacos(-1.d0)

!     Parameters : Rb85
      i = 5.d0/2.d0
      s = 1.d0/2.d0
      gE = 2.d0*1.001159652193d0
      muE = -1.d0/(2.d0*(1.d0)) 
      gN = 1.35302d0/i
      muN = 1.d0/(2.d0*1836.2d0)  

      if (gN.gt.0.d0) Lambda = 1.d0
      if (gN.lt.0.d0) Lambda =-1.d0
      Delta = 1.011910813d9*Hz2Au ! c
      xi = Delta !Lambda*(1.d0*Delta/(2.d0*i+1.d0))

!     Determining the number of hyperfine atomic states
      NumHypfStates = (nint(2.d0*i)+1)*(nint(2.d0*s)+1)

      return
      end

        subroutine Parameters_87Rb(i,s,xi,gE,muE,gN,muN,NumHypfStates)
        integer NumHypfStates
        double precision i,s
        double precision Hz2Au,Pi,Gauss
        double precision Delta,Lambda,xi
        double precision gE,muE,gN,muN

        Hz2Au = 1.519829846006d-16
        Pi = dacos(-1.d0)

        !     Parameters : Rb87
        i = 3.d0/2.d0
        s = 1.d0/2.d0
        gE = 2.d0*1.001159652193d0
        muE = -1.d0/(2.d0*(1.d0))
        gN = 0.0009951414d0*1836.2d0
        muN = 1.d0/(2.d0*1836.2d0)

        if (gN.gt.0.d0) Lambda = 1.d0
        if (gN.lt.0.d0) Lambda =-1.d0
        Delta = 6.834682d9/2d0*Hz2Au ! c
        xi = Delta !Lambda*(1.d0*Delta/(2.d0*i+1.d0))

        !     Determining the number of hyperfine atomic states
        NumHypfStates = (nint(2.d0*i)+1)*(nint(2.d0*s)+1)
        return
        end


