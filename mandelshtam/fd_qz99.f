        program filter_diagonalization
!
! Vladimir Mandelshtam, last modified Aug. 1998
! Please don't distribute without permission.
! The author is not responsible for any damage caused by this code.
! Report any bugs found (mandelsh@uci.edu).
! Any suggestions are welcome.
! Please let me know how it works.
!
!
!
!  This program for a signal c_n=C(t0+tau*n) n=0,1,...,Nsig
!  obtains the complex frequencies w_k and amplitudes d_k
!  by fitting to the form
!                          c_n= sum_k d_k exp(-inw_k)         (1)
!
!  The coefficients are phase corrected by d_k=d_k*exp(i*theta)
!
!  The spectrum is computed using the formula:
!
!            F(w)=Im {sum_k d_k/(w-w_k+i*Gamm)}       (2)
!
!     Cheating:    replace w_k in (2) by Re{w_k}+i*cheat*Im{w_k}
!
!     If Im{w_k}>0 it is replaced by -Im{w_k} in (2).
!
!
! Remark: For solving the generalized eigenvalue problem we use QZ-algorithm
!
! Remark: First try to use rho=1, error=1000, 
!         Gcut=10000, Gamm=0.000001, ispec=1, cheat=1
!
! Remark: Try to use such frequency windows for which Nb 
!         (the size of the matrices) becomes >3 and <300
!
! Remark: If there is no noise, i.e. (1) is exact, the choice for Nb is crucial
!         In such a case the best convergence will be achieved when Nb is 
!         just slightly greater than the actual rank of the U0 matrix.
!         However, this never happens for experimental signals.
!
      implicit none

      integer Nb,Nb2,i,j,n,k,M,Nsig,ispec,Npower,idat,
     &   Nbmax,Nsigmax,Npowermax,Nsig1max,Nskip,Dim

      parameter (Nbmax=500,Nsigmax=500001,Npowermax=10000)

      complex*16 U(Nbmax,Nbmax,0:2),uk(Nbmax),
     &     wk(Nbmax),dk(Nbmax),
     &     AL(Nbmax,Nbmax),AR(Nbmax,Nbmax),uu,phase_corr,
     &     f(2*Nbmax,2,0:2),xi,z(2*Nbmax),coef(0:Nsigmax)

      real*8 err(Nbmax),Spec(0:Npowermax),theta,cheat,
     &     pi,Omega,Gamm,t0,delt,error,wmin,wmax,Gcut,
     &     wmin_old,wmax_old,c1,c2,tt,rho,shift,delt1

      character*30 par,signal,AbsSp,ReSp

      logical choose_nb,kcosine,cheatmore

!
!       This arrays are needed for qz only    
!
      real*8 Zr(Nbmax,Nbmax,2),wr(Nbmax,3)
!
      read(5,*) signal, idat ! the input data file with the signal given as a complex column vector.
                             ! idat>0 read real signal 
                             ! idat<0 read complex signal 
                             ! |idat|=1 read c_n
                             ! |idat|=2 read t_n,c_n
                             ! idat=-3 read t_n,Re(c_n),Im(c_n)
                             ! idat=-4 read Re(c_n),Im(c_n)


c      read(5,*) delt         ! the time step should be given in the units of 1/[frequency units]
                              ! here tau is included in the file with the signal

      read(5,*) t0, theta     ! t0 --- the time delay (we assume that the first time increment is t0
                              ! theta ---- the overall phase correction   exp(i*theta)
      read(5,*) ispec         ! ispec=0  ---- Do spectral analysis with FD and FD/DFT hybrid method
                              ! ispec=1  ---- No DFT
                              ! ispec=-1 ---- DFT only

      read(5,*) par     ! file with spectral parameters

      read(5,*) ReSp    ! Absorption spectrum constructed out of the spectral parameters
                              ! theta=pi/2 will replace absorption by dispersion

      read(5,*) wmin_old,wmax_old  ! the frequency window. 
                                   ! The units are such that the Nyquist range = 2*pi/tau

      read(5,*) Nsig         ! numner of signal data points to be used in the analysis    

      read(5,*) Nskip        ! number of the signal points to skip before the analysis    

      read(5,*) rho,kcosine     ! The density of the window basis functions 
                                ! in the units of the "optimal" basis (rho=1 --- optimal)
                                ! if rho<0 Nb=-rho
                                ! kcosine=.true. - implement the cosine window
                                ! kcosine=.false. - implement the rectangular window


      read(5,*) error         !	error is used to check the errors for the U1 eigenvalues by
                              ! calculating ||(U2 B_k,B_k) - u_minus_k ||
                              ! use error=100000 to avoid missing eigenvalues

      read(5,*) Npower,Gamm,Gcut ! Npower --- Number of frequency points to plot the spectrum 
                                ! Gamm ----  the smoothing width to construct the spectrum
                              ! Gcut ---- maximum width for a pole to be used to construct the spectrum
                          
      read(5,*)cheat,cheatmore ! cheat ---- multiply all widths by cheat (be very careful or use cheat=1 !!!)
                               ! cheatmore=.true. the spectrum is computed with Im d_k
      if(Nsig.gt.Nsigmax) stop 'increase Nsigmax'
      if(Npower.gt.Npowermax) stop 'increase Npowermax'
      do j=0,Npower
         Spec(j)=0d0
      enddo
      open(15,file=signal)
      write(6,*) 'read the signal from disk'
      read(15,*) Dim,Nsig1max,delt1
      delt=delt1
      pi=dacos(-1d0)
      xi=(0d0,1d0)
      phase_corr=cdexp(dcmplx(0d0,theta))
      wmin=wmin_old*delt
      if(wmin.lt.-pi) wmin=-pi
      wmax=wmax_old*delt
      if(wmax.gt.pi) wmax=pi
      wmin_old=wmin/delt
      wmax_old=wmax/delt
      if(Nsig.gt.Nsig1max-1) Nsig=Nsig1max-1
      do n=1,Nskip
         read(15,*)
      enddo
      c2=0d0
      do n=0,Nsig
         if(idat.eq.-3) read(15,*,end=100) tt,c1,c2
         if(idat.eq.-4) read(15,*,end=100) c1,c2
         if(idat.eq.1)  read(15,*,end=100) c1
         if(idat.eq.2)  read(15,*,end=100) tt,c1
         if(idat.eq.-1) read(15,*,end=100) coef(n)
         if(idat.eq.-2) read(15,*,end=100) tt,coef(n)
         if(idat.ne.-1.and.idat.ne.-2) coef(n)=dcmplx(c1,c2)
         coef(n)=coef(n)*phase_corr
      enddo
      goto 101
 100  Nsig=n-1
 101  M=(Nsig-2)/2
      Nsig=2*M+2
      write(6,*) 'M=',M,' which means that', 2*M+3, 'c_n are used'
      if(rho.gt.0d0) then
         choose_nb=.true.
      else
         Nb=(-rho)*1.00000001
      endif
      if(ispec.ge.0) then
         If (choose_nb) then
            write(6,*) 'Use the grid with rho=',rho
            Nb=((wmax-wmin)/(2*pi))*(M+1)*rho
         endif
         shift=pi/(2*(M+1))
         if(Nb.lt.3) Nb=3
         if(Nb.gt.Nbmax) stop 'Nb>Nbmax'
         do i=1,Nb
            z(i)=wmin+(i-0.5d0)*(wmax-wmin)/Nb
            if(kcosine) then 
               z(i+Nb)=cdexp(xi*(z(i)+shift))
               z(i)=z(i)-shift
            endif
            z(i)=cdexp(xi*z(i)) ! this is actually 1/z_j as defined in the paper
         enddo
         write(6,*) 'Nb=',Nb
      endif
      write(6,*) ' wmin,wmax=',wmin_old,wmax_old
!
!_______________________________________________________________
!
      Nb2=Nb
      if(kcosine) then
         Nb2=2*Nb
         write(6,*) 'Use cosine window'
      endif
      call FD(Nsig,coef,delt,wmin,wmax,Nb,Nb2,error,
     &     z,t0,ispec,Npower,cheat,cheatmore,Gamm,Gcut,
     &     wk,dk,err,U,uk,f,AL,AR,Zr,wr,Spec)
      open(14,file=ReSp)
      do i=0,Npower
         write(14,986) (wmin+i*((wmax-wmin)/Npower))/delt,Spec(i)
      enddo
 986  format(E14.7,E12.4)
      if(ispec.lt.0) stop
      open(7,file=par)   	
      write(7,79) 
 79   format(8x,'Re w',18x,'Im w',15x,'Re d',14x,'Im d',9x,'error')
 13   format(E22.16,3E18.10,E9.2)
      call cpiksrt(Nb,wk,dk,err)
      do k=1,Nb
         if(err(k).lt.error) 
     &        write(7,13) dreal(wk(k)),dimag(wk(k)),
     &        dreal(dk(k)),dimag(dk(k)),err(k)
      enddo
      stop
      end
!
!_______________________________________________________________
!
      subroutine FD(Nsig,coef,delt,wmin,wmax,Nb,Nb2,error,
     &     z,t0,ispec,Npower,cheat,cheatmore,Gamm,Gcut,
     &     wk,dk,err,U,uk,f,AL,Ar,Zr,wr,Spec)
!
c      implicit real*8(a-h,o-z)
      implicit none
      integer Nb,Nb2,i,j,n,k,k1,M,Nsig,ispec,ierr,N0,Npower,i1,j1,l
      complex*16 U(Nb,Nb,0:2),uk(Nb),coef(0:Nsig),
     &     wk(Nb),dk(Nb),f(Nb2,2,0:2),z(Nb2),
     &     Z1,Z2,Power,xi,ss,uu,cc,f1,f2,c1,s1
      real*8 err(Nb),Gcut,cheat
      logical cheatmore
!
      real*8 AL(Nb,Nb,2),AR(Nb,Nb,2),Zr(Nb,Nb,2),wr(Nb,3),
     &     Spec(0:Npower)
!
      real*8 rho,pi,Omega,Gamm,t0,delt,yiter,error,wmin,wmax
      M=(Nsig-2)/2
      pi=dacos(-1d0)
      xi=(0d0,1d0)
      if(ispec.lt.0) goto 399
      write(6,*) 'construction of small U0,U1 and U2 matrices'
      do l=0,2
         do i=1,Nb
            do j=1,Nb
               U(j,i,l)=(0,0d0)
            enddo
         enddo
      enddo
      do i=1,Nb2
         f(i,1,1)=(0d0,0)
         f(i,2,1)=(0d0,0)
         cc=(0,0d0)
         ss=coef(1)
         s1=(0,0d0)
         c1=(0,0d0)
         f1=(0,0d0)
         f2=(0,0d0)
         do k=2,M+1
            if(mod(k,100).eq.2) then
               Z1=z(i)**(k-1)
            else
               Z1=Z1*z(i)
            endif
            f1=f1+(coef(k)*Z1)
            s1=s1+k*(coef(k)*Z1)
            f2=f2+(coef(k+M)*Z1)
            c1=c1+k*(coef(k+M)*Z1)
            if(mod(k,100).eq.2.or.k.eq.M+1) then
               f(i,1,1)=f(i,1,1)+f1
               f(i,2,1)=f(i,2,1)+f2
               ss=ss+s1
               cc=cc+c1
               s1=(0,0d0)
               c1=(0,0d0)
               f1=(0,0d0)
               f2=(0,0d0)
            endif
         enddo
         Z2=Z1**2
         f(i,2,1)=f(i,2,1)*Z1
c         ss=ss-cc*Z1+M*f(i,2,1)    !!!!!!! MISTAKE corrected 01/30/99
         ss=ss-cc*Z1+(M+2)*f(i,2,1)
         f(i,1,2)=f(i,1,1)/z(i)-coef(2)+coef(M+2)*Z1
         f(i,2,2)=f(i,2,1)/z(i)-coef(M+2)*Z1+coef(2*M+2)*Z2
         f(i,1,0)=(f(i,1,1)+coef(1)-coef(M+1)*Z1)*z(i)
         f(i,2,0)=(f(i,2,1)-coef(2*M+1)*Z2+coef(M+1)*Z1)*z(i)
         j=i
         if(i.gt.Nb) j=i-Nb
         U(j,j,1)=U(j,j,1)+ss
         U(j,j,2)=U(j,j,2)+(ss-coef(1)-f(i,1,1)+f(i,2,1))/z(i)
     &        +coef(2*M+2)*Z2
         U(j,j,0)=U(j,j,0)+(ss+coef(1)-2*coef(M+1)*Z1+f(i,1,1)
     &        -f(i,2,1))*z(i)+coef(0)
      enddo
      do i=1,Nb2
         i1=i
         if(i.gt.Nb) i1=i-Nb
         do j=1,Nb2
            if(j.ne.i) then
               j1=j
               if(j.gt.Nb) j1=j-Nb
               cc=z(i)/z(j)
               do l=0,2
                  U(i1,j1,l)=U(i1,j1,l)+coef(l)
     &                 +(f(j,1,l)-cc**(M+1)*f(j,2,l)
     &                 -cc*f(i,1,l)+cc**(-M)*f(i,2,l))/(1-cc)
               enddo
            endif
         enddo
      enddo
c
      write(6,*)
      write(6,*) 'Solving: U1 Bk=uk U0 Bk'
      write(6,*)
c
      call QZ(Nb,U(1,1,1),U(1,1,0),U(1,1,1),AL,AR,Zr,uk,wr)
      do k=1,Nb
         wk(k)=(0,1d0)*cdlog(uk(k))/delt    ! get the frequency in the correct units
      enddo
c_________________________________________________________________
c   Estimate errors using U2 
      do k=1,Nb
         z1=(0,0d0)
         do j=1,Nb
            cc=(0,0d0)
            do i=1,Nb
               cc=cc+U(i,j,2)*dcmplx(Zr(i,k,1),Zr(i,k,2))
            enddo
            z1=z1+cc*dcmplx(Zr(j,k,1),Zr(j,k,2))
         enddo
         err(k)=cdabs(wk(k)-(0,1d0)*cdlog(z1)/delt)
      enddo
      write(6,*) 'compute the coefficients'
      call coeff(Nb,Nb2,M,coef,dk,uk,err,error,Zr,z,f(1,1,1))
      do k=1,Nb
         dk(k)=dk(k)**2
     &        *cdexp((0,1d0)*wk(k)*(t0+delt)) !to adjust the t=0
      enddo
c
      if(ispec.eq.0) write(6,*) 'Split the signal'
      do k=1,Nb
         if(err(k).lt.error) then
            if(dk(k).ne.(0,0d0).and.ispec.eq.0       !subtruct the signal part
     &              .and.dimag(wk(k)).gt.-Gcut
     &              .and.dreal(wk(k)).lt.wmax/delt
     &              .and.dreal(wk(k)).gt.wmin/delt) then 
               if(dimag(wk(k)).gt.0d0) then
                  uu=cdexp(-(0,1d0)*dconjg(wk(k))*delt)
               else
                  uu=cdexp(-(0,1d0)*wk(k)*delt)
               endif
               Z1=dk(k)
     &              *cdexp(-(0,1d0)*wk(k)*(t0-delt)) !to unadjust the t=0
               do n=0,2*M+2
                  Z1=Z1*uu
                  coef(n)=coef(n)-Z1
               enddo
            endif
         endif
      enddo
c     
 399  write(6,*) 'Compute the spectrum'
      if(ispec.lt.0) write(6,*) 'Only DFT is used'
      N0=M/5
      if(ispec.lt.1) then       !before DFT use a window
         do n=0,2*M+2
            if (n.gt.N0) then
               yiter = dfloat(n-N0)/dfloat(2*M+3-N0)
               coef(n)=coef(n)*dexp(yiter)*(1D0-yiter)
            endif
         enddo
      endif
      do i=0,Npower
         Omega=wmin+i*((wmax-wmin)/Npower)
         Power=0d0
         if(ispec.lt.1) then    !do DFT
            uu=cdexp((0,1d0)*Omega)
            Z1=-delt*(0,1d0)*cdexp(dcmplx(0d0,t0*Omega/delt))
            Power=coef(0)*Z1/2 
            do n=1,2*M+2
               Z1=Z1*uu
               Power=Power+coef(n)*Z1
            enddo
         endif
         if(ispec.ge.0) then
            Omega=Omega/delt
            do k=1,Nb 
               cc=dcmplx(dreal(wk(k)),
     &              -cheat*dabs(dimag(wk(k))))
               if(-dimag(cc).lt.gamm) cc=dcmplx(dreal(cc),-gamm)
               if(err(k).lt.error.and.ispec.ge.0.and.
     &              dimag(wk(k)).gt.-Gcut ! ) then
     &              .and.dreal(wk(k)).lt.wmax/delt.and.
     &              dreal(wk(k)).gt.wmin/delt) then 
                  if(cheatmore) then
                     Spec(i)=Spec(i)-dimag(dk(k))*dimag(1/(Omega-cc))
                  else
                     Spec(i)=Spec(i)+dreal(dk(k)/(Omega-cc))
                  endif
               endif
            enddo
         endif
         Spec(i)=Spec(i)+dreal(Power)
      enddo
c
      return
      end
c
      subroutine coeff(Nb,Nb2,M,coef,dk,uk,err,error,Zr,z,f)
      implicit none
      integer Nb,Nb2,k,n,M,j,j1
      real*8 rho,err(Nb),Zr(Nb,Nb,2),error
      complex*16 uu,cc,Z1,Z2,
     &     z(Nb),uk(Nb),f(Nb2,2),dk(Nb),coef(0:100)
      do 88 k=1,Nb
         dk(k)=(0,0d0)
         if(err(k).gt.error) goto 88
         rho=cdabs(uk(k))
         if(rho.gt.0.999999999d0) then
            uu=(1d0,0)/uk(k)
         else
            uu=rho/uk(k)
         endif
         Z1 = coef(M+1)*uu
         Z2 = coef(2*M+1)*uu
         do n=M,2,-1
            Z1 = (Z1+coef(n))*uu
            Z2 = (Z2+coef(M+n))*uu
         enddo
         do j=1,Nb2
            j1=j
            if(j.gt.Nb) j1=j-Nb
            cc=uu/z(j)
            dk(k)=dk(k)+dcmplx(Zr(j1,k,1),Zr(j1,k,2))*
     &           (coef(1)+(f(j,1)-cc**(M+1)*f(j,2)-
     &           cc*Z1+z(j)**M*Z2)/(1-cc))
         enddo
         if(rho.gt.0.999999999d0) then
            dk(k)=dk(k)/(M+1)
         else
            dk(k)=dk(k)/
     &           ((1d0,0)-rho**(M+1))*((1d0,0)-rho)
         endif
 88      continue
      return
      end
c
      subroutine cpiksrt(Nb,wr,wi,err)
c the frequencies are increasing.
      implicit real*8(a-h,o-z)
      complex*16 wr(Nb),wi(Nb),war,wai
      real*8 err(Nb),wae
      do j=2,Nb
         war=wr(j)
         wai=wi(j)
         wae=err(j)
         do i=j-1,1,-1
            if(dreal(wr(i)).le.dreal(war))go to 10
            wr(i+1)=wr(i)
            wi(i+1)=wi(i)
            err(i+1)=err(i)
         enddo
         i=0
 10      continue
         wr(i+1)=war
         wi(i+1)=wai
         err(i+1)=wae
      enddo
      return
      end
C
C     ------------------------------------------------------------------
      subroutine cpiksrt1(Nb,wr,wi,err)
c the frequencies are decreasing.
      implicit real*8(a-h,o-z)
      complex*16 wr(Nb),wi(Nb),war,wai
      real*8 err(Nb),wae
      do j=2,Nb
         war=wr(j)
         wai=wi(j)
         wae=err(j)
         do i=j-1,1,-1
            if(dreal(wr(i)).ge.dreal(war))go to 10
            wr(i+1)=wr(i)
            wi(i+1)=wi(i)
            err(i+1)=err(i)
         enddo
         i=0
 10      continue
         wr(i+1)=war
         wi(i+1)=wai
         err(i+1)=wae
      enddo
      return
      end
c
      subroutine QZ(Nb,UL,UR,U0,AL,AR,Zr,uk,wr)
      implicit none
      integer Nb,i,j,k,n,ierr
      complex*16 UL(Nb,Nb),UR(Nb,Nb),U0(Nb,Nb),uk(Nb),
     &     xi,cc,uu
      real*8 AR(Nb,Nb,2),AL(Nb,Nb,2),Zr(Nb,Nb,2),wr(Nb,3)
      do i=1,Nb
         do j=1,Nb
            AL(i,j,1)=dreal(UL(i,j))
            AL(i,j,2)=dimag(UL(i,j))
         enddo
      enddo
      do i=1,Nb
         do j=1,Nb
            AR(i,j,1)=dreal(UR(i,j))
            AR(i,j,2)=dimag(UR(i,j))
         enddo
      enddo
      write(6,*) 'CQZHES'
      call CQZHES(Nb,Nb,AL(1,1,1),AL(1,1,2),
     &     Ar(1,1,1),Ar(1,1,2),.true.,Zr(1,1,1),Zr(1,1,2))
      write(6,*) 'CQZVAL'
      call CQZVAL(Nb,Nb,AL(1,1,1),AL(1,1,2),
     &     Ar(1,1,1),Ar(1,1,2),0d0,wr(1,1),wr(1,2),wr(1,3),
     &     .true.,Zr(1,1,1),Zr(1,1,2),IERR)
      if(IERR.ne.0) then
         write(6,*) '#########################'
         write(6,*) '###  IERR=',IERR
         write(6,*) '#########################'
      endif
      write(6,*) 'CQZVEC'
      call CQZVEC(Nb,Nb,AL(1,1,1),AL(1,1,2),
     &     Ar(1,1,1),Ar(1,1,2),wr(1,1),wr(1,2),wr(1,3),
     &     Zr(1,1,1),Zr(1,1,2))
c_________________________________________________________________
      do k=1,Nb
         uk(k)=dcmplx(wr(k,1),wr(k,2))/wr(k,3)
         cc=(0,0d0)
         do j=1,Nb
            uu=0.0d0
            do n=1,Nb
               uu=uu+U0(j,n)*dcmplx(Zr(n,k,1),Zr(n,k,2))
            enddo
            cc=cc+dcmplx(Zr(j,k,1),Zr(j,k,2))*uu
         enddo
         cc=(1d0,0)/cdsqrt(cc)
         do j=1,Nb
            uu=cc*dcmplx(Zr(j,k,1),Zr(j,k,2))
            Zr(j,k,1)=dreal(uu)
            Zr(j,k,2)=dimag(uu)
         enddo
      enddo
      return
      end
C
C     ------------------------------------------------------------------
C
      SUBROUTINE CQZHES(NM,N,AR,AI,BR,BI,MATZ,ZR,ZI)
C
      INTEGER I,J,K,L,N,K1,LB,L1,NM,NK1,NM1
      REAL*8 AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ZR(NM,N),ZI(NM,N)
      REAL*8 R,S,T,TI,U1,U2,XI,XR,YI,YR,RHO,U1I
      REAL*8 DSQRT,CDABS,DABS
      LOGICAL MATZ
      COMPLEX*16 DCMPLX
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-
C     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR
C     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY
C     CQZVAL  AND POSSIBLY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE
C          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE),
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
C          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
C          OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** INITIALIZE Z **********
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 I = 1, N
C
         DO 2 J = 1, N
            ZR(I,J) = 0.0
            ZI(I,J) = 0.0
    2    CONTINUE
C
         ZR(I,I) = 1.0
    3 CONTINUE
C     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH
C                TEMPORARILY REAL DIAGONAL ELEMENTS **********
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
C
      DO 100 L = 1, NM1
         L1 = L + 1
         S = 0.0
C
         DO 20 I = L, N
            S = S + DABS(BR(I,L)) + DABS(BI(I,L))
   20    CONTINUE
C
         IF (S .EQ. 0.0) GO TO 100
         RHO = 0.0
C
         DO 25 I = L, N
            BR(I,L) = BR(I,L) / S
            BI(I,L) = BI(I,L) / S
            RHO = RHO + BR(I,L)**2 + BI(I,L)**2
   25    CONTINUE
C
         R = DSQRT(RHO)
         XR = CDABS(DCMPLX(BR(L,L),BI(L,L)))
         IF (XR .EQ. 0.0) GO TO 27
         RHO = RHO + XR * R
         U1 = -BR(L,L) / XR
         U1I = -BI(L,L) / XR
         YR = R / XR + 1.0
         BR(L,L) = YR * BR(L,L)
         BI(L,L) = YR * BI(L,L)
         GO TO 28
C
   27    BR(L,L) = R
         U1 = -1.0
         U1I = 0.0
C
   28    DO 50 J = L1, N
            T = 0.0
            TI = 0.0
C
            DO 30 I = L, N
               T = T + BR(I,L) * BR(I,J) + BI(I,L) * BI(I,J)
               TI = TI + BR(I,L) * BI(I,J) - BI(I,L) * BR(I,J)
   30       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 40 I = L, N
               BR(I,J) = BR(I,J) - T * BR(I,L) + TI * BI(I,L)
               BI(I,J) = BI(I,J) - T * BI(I,L) - TI * BR(I,L)
   40       CONTINUE
C
            XI = U1 * BI(L,J) - U1I * BR(L,J)
            BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
            BI(L,J) = XI
   50    CONTINUE
C
         DO 80 J = 1, N
            T = 0.0
            TI = 0.0
C
            DO 60 I = L, N
               T = T + BR(I,L) * AR(I,J) + BI(I,L) * AI(I,J)
               TI = TI + BR(I,L) * AI(I,J) - BI(I,L) * AR(I,J)
   60       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 70 I = L, N
               AR(I,J) = AR(I,J) - T * BR(I,L) + TI * BI(I,L)
               AI(I,J) = AI(I,J) - T * BI(I,L) - TI * BR(I,L)
   70       CONTINUE
C
            XI = U1 * AI(L,J) - U1I * AR(L,J)
            AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
            AI(L,J) = XI
   80    CONTINUE
C
         BR(L,L) = R * S
         BI(L,L) = 0.0
C
         DO 90 I = L1, N
            BR(I,L) = 0.0
            BI(I,L) = 0.0
   90    CONTINUE
C
  100 CONTINUE
C     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
C                ELEMENTS, WHILE KEEPING B TRIANGULAR **********
      DO 160 K = 1, NM1
         K1 = K + 1
C     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL **********
         IF (AI(N,K) .EQ. 0.0) GO TO 105
         R = CDABS(DCMPLX(AR(N,K),AI(N,K)))
         U1 = AR(N,K) / R
         U1I = AI(N,K) / R
         AR(N,K) = R
         AI(N,K) = 0.0
C
         DO 103 J = K1, N
            XI = U1 * AI(N,J) - U1I * AR(N,J)
            AR(N,J) = U1 * AR(N,J) + U1I * AI(N,J)
            AI(N,J) = XI
  103    CONTINUE
C
         XI = U1 * BI(N,N) - U1I * BR(N,N)
         BR(N,N) = U1 * BR(N,N) + U1I * BI(N,N)
         BI(N,N) = XI
  105    IF (K .EQ. NM1) GO TO 170
         NK1 = NM1 - K
C     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     ********** ZERO A(L+1,K) **********
            S = DABS(AR(L,K)) + DABS(AI(L,K)) + AR(L1,K)
            IF (S .EQ. 0.0) GO TO 150
            U1 = AR(L,K) / S
            U1I = AI(L,K) / S
            U2 = AR(L1,K) / S
            R = DSQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            AR(L,K) = R * S
            AI(L,K) = 0.0
            AR(L1,K) = 0.0
C
            DO 110 J = K1, N
               XR = AR(L,J)
               XI = AI(L,J)
               YR = AR(L1,J)
               YI = AI(L1,J)
               AR(L,J) = U1 * XR + U1I * XI + U2 * YR
               AI(L,J) = U1 * XI - U1I * XR + U2 * YI
               AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110       CONTINUE
C
            XR = BR(L,L)
            BR(L,L) = U1 * XR
            BI(L,L) = -U1I * XR
            BR(L1,L) = -U2 * XR
C
            DO 120 J = L1, N
               XR = BR(L,J)
               XI = BI(L,J)
               YR = BR(L1,J)
               YI = BI(L1,J)
               BR(L,J) = U1 * XR + U1I * XI + U2 * YR
               BI(L,J) = U1 * XI - U1I * XR + U2 * YI
               BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  120       CONTINUE
C     ********** ZERO B(L+1,L) **********
            S = DABS(BR(L1,L1)) + DABS(BI(L1,L1)) + DABS(BR(L1,L))
            IF (S .EQ. 0.0) GO TO 150
            U1 = BR(L1,L1) / S
            U1I = BI(L1,L1) / S
            U2 = BR(L1,L) / S
            R = DSQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            BR(L1,L1) = R * S
            BI(L1,L1) = 0.0
            BR(L1,L) = 0.0
C
            DO 130 I = 1, L
               XR = BR(I,L1)
               XI = BI(I,L1)
               YR = BR(I,L)
               YI = BI(I,L)
               BR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               BI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               BR(I,L) = U1 * YR - U1I * YI - U2 * XR
               BI(I,L) = U1 * YI + U1I * YR - U2 * XI
  130       CONTINUE
C
            DO 140 I = 1, N
               XR = AR(I,L1)
               XI = AI(I,L1)
               YR = AR(I,L)
               YI = AI(I,L)
               AR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               AI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               AR(I,L) = U1 * YR - U1I * YI - U2 * XR
               AI(I,L) = U1 * YI + U1I * YR - U2 * XI
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
C
            DO 145 I = 1, N
               XR = ZR(I,L1)
               XI = ZI(I,L1)
               YR = ZR(I,L)
               YI = ZI(I,L)
               ZR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               ZI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               ZR(I,L) = U1 * YR - U1I * YI - U2 * XR
               ZI(I,L) = U1 * YI + U1I * YR - U2 * XI
  145       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
  170 RETURN
C     ********** LAST CARD OF CQZHES **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE CQZVAL(NM,N,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,
     X                                       MATZ,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,N,EN,K1,K2,LL,L1,NA,NM,ITS,KM1,LM1,
     X        ENM2,IERR,LOR1,ENORN
      REAL*8 AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),
     X       BETA(N),ZR(NM,N),ZI(NM,N)
      REAL*8 R,S,A1,A2,EP,SH,U1,U2,XI,XR,YI,YR,ANI,A1I,A33,A34,A43,A44,
     X       BNI,B11,B33,B44,SHI,U1I,A33I,A34I,A43I,A44I,B33I,B44I,
     X       EPSA,EPSB,EPS1,ANORM,BNORM,B3344,B3344I
      REAL*8 DSQRT,CDABS,DABS
      INTEGER MAX0
      LOGICAL MATZ
      COMPLEX*16 Z3
      COMPLEX*16 CDSQRT,DCMPLX
      REAL*8 DREAL,DIMAG
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,
C     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
C     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING
C     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS
C     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS
C     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY
C     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX
C          WITH REAL SUBDIAGONAL ELEMENTS,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE,
C
C        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET
C          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO
C          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,
C
C        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE
C          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN
C          THE RATIOS ((ALFR+I*ALFI)/BETA),
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF AR(J,J-1) HAS NOT BECOME
C                     ZERO AFTER 50 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     ********** COMPUTE EPSA,EPSB **********
      ANORM = 0.0
      BNORM = 0.0
C
      DO 30 I = 1, N
         ANI = 0.0
         IF (I .NE. 1) ANI = DABS(AR(I,I-1))
         BNI = 0.0
C
         DO 20 J = I, N
            ANI = ANI + DABS(AR(I,J)) + DABS(AI(I,J))
            BNI = BNI + DABS(BR(I,J)) + DABS(BI(I,J))
   20    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
C
      IF (ANORM .EQ. 0.0) ANORM = 1.0
      IF (BNORM .EQ. 0.0) BNORM = 1.0
      EP = EPS1
      IF (EP .GT. 0.0) GO TO 50
C     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
      EP = 1.0
   40 EP = EP / 2.0
      IF (1.0 + EP .GT. 1.0) GO TO 40
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
C     ********** REDUCE A TO TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR **********
      LOR1 = 1
      ENORN = N
      EN = N
C     ********** BEGIN QZ STEP **********
   60 IF (EN .EQ. 0) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- **********
   70 DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (DABS(AR(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 AR(L,LM1) = 0.0
C     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********
   95 B11 = CDABS(DCMPLX(BR(L,L),BI(L,L)))
      IF (B11     .EQ. 0.0) GO TO 98
      U1 = BR(L,L) / B11
      U1I = BI(L,L) / B11
C
      DO 97 J = L, ENORN
         XI = U1 * AI(L,J) - U1I * AR(L,J)
         AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
         AI(L,J) = XI
         XI = U1 * BI(L,J) - U1I * BR(L,J)
         BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
         BI(L,J) = XI
   97 CONTINUE
C
      BI(L,L) = 0.0
   98 IF (L .NE. EN) GO TO 100
C     ********** 1-BY-1 BLOCK ISOLATED **********
      ALFR(EN) = AR(EN,EN)
      ALFI(EN) = AI(EN,EN)
      BETA(EN) = B11
      EN = NA
      GO TO 60
C     ********** CHECK FOR SMALL TOP OF B **********
  100 L1 = L + 1
      IF (B11 .GT. EPSB) GO TO 120
      BR(L,L) = 0.0
      S = DABS(AR(L,L)) + DABS(AI(L,L)) + DABS(AR(L1,L))
      U1 = AR(L,L) / S
      U1I = AI(L,L) / S
      U2 = AR(L1,L) / S
      R = DSQRT(U1*U1+U1I*U1I+U2*U2)
      U1 = U1 / R
      U1I = U1I / R
      U2 = U2 / R
      AR(L,L) = R * S
      AI(L,L) = 0.0
C
      DO 110 J = L1, ENORN
         XR = AR(L,J)
         XI = AI(L,J)
         YR = AR(L1,J)
         YI = AI(L1,J)
         AR(L,J) = U1 * XR + U1I * XI + U2 * YR
         AI(L,J) = U1 * XI - U1I * XR + U2 * YI
         AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
         XR = BR(L,J)
         XI = BI(L,J)
         YR = BR(L1,J)
         YI = BI(L1,J)
         BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         BR(L,J) = U1 * XR + U1I * XI + U2 * YR
         BI(L,J) = U1 * XI - U1I * XR + U2 * YI
         BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110 CONTINUE
C
      LM1 = L
      L = L1
      GO TO 90
C     ********** ITERATION STRATEGY **********
  120 IF (ITS .EQ. 50) GO TO 1000
      IF (ITS .EQ. 10) GO TO 135
C     ********** DETERMINE SHIFT **********
      B33 = BR(NA,NA)
      B33I = BI(NA,NA)
      IF (CDABS(DCMPLX(B33,B33I)) .GE. EPSB) GO TO 122
      B33 = EPSB
      B33I = 0.0
  122 B44 = BR(EN,EN)
      B44I = BI(EN,EN)
      IF (CDABS(DCMPLX(B44,B44I)) .GE. EPSB) GO TO 124
      B44 = EPSB
      B44I = 0.0
  124 B3344 = B33 * B44 - B33I * B44I
      B3344I = B33 * B44I + B33I * B44
      A33 = AR(NA,NA) * B44 - AI(NA,NA) * B44I
      A33I = AR(NA,NA) * B44I + AI(NA,NA) * B44
      A34 = AR(NA,EN) * B33 - AI(NA,EN) * B33I
     X    - AR(NA,NA) * BR(NA,EN) + AI(NA,NA) * BI(NA,EN)
      A34I = AR(NA,EN) * B33I + AI(NA,EN) * B33
     X     - AR(NA,NA) * BI(NA,EN) - AI(NA,NA) * BR(NA,EN)
      A43 = AR(EN,NA) * B44
      A43I = AR(EN,NA) * B44I
      A44 = AR(EN,EN) * B33 - AI(EN,EN) * B33I - AR(EN,NA) * BR(NA,EN)
      A44I = AR(EN,EN) * B33I + AI(EN,EN) * B33 - AR(EN,NA) * BI(NA,EN)
      SH = A44
      SHI = A44I
      XR = A34 * A43 - A34I * A43I
      XI = A34 * A43I + A34I * A43
      IF (XR .EQ. 0.0 .AND. XI .EQ. 0.0) GO TO 140
      YR = (A33 - SH) / 2.0
      YI = (A33I - SHI) / 2.0
      Z3 = CDSQRT(DCMPLX(YR**2-YI**2+XR,2.0*YR*YI+XI))
      U1 = DREAL(Z3)
      U1I = DIMAG(Z3)
      IF (YR * U1 + YI * U1I .GE. 0.0) GO TO 125
      U1 = -U1
      U1I = -U1I
  125 Z3 = (DCMPLX(SH,SHI) - DCMPLX(XR,XI) / DCMPLX(YR+U1,YI+U1I))
     X   / DCMPLX(B3344,B3344I)
      SH = DREAL(Z3)
      SHI = DIMAG(Z3)
      GO TO 140
C     ********** AD HOC SHIFT **********
  135 SH = AR(EN,NA) + AR(NA,ENM2)
      SHI = 0.0
C     ********** DETERMINE ZEROTH COLUMN OF A **********
  140 A1 = AR(L,L) / B11 - SH
      A1I = AI(L,L) / B11 - SHI
      A2 = AR(L1,L) / B11
      ITS = ITS + 1
      IF (.NOT. MATZ) LOR1 = L
C     ********** MAIN LOOP **********
      DO 260 K = L, NA
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX0(K-1,L)
C     ********** ZERO A(K+1,K-1) **********
         IF (K .EQ. L) GO TO 170
         A1 = AR(K,KM1)
         A1I = AI(K,KM1)
         A2 = AR(K1,KM1)
  170    S = DABS(A1) + DABS(A1I) + DABS(A2)
         U1 = A1 / S
         U1I = A1I / S
         U2 = A2 / S
         R = DSQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
C
         DO 180 J = KM1, ENORN
            XR = AR(K,J)
            XI = AI(K,J)
            YR = AR(K1,J)
            YI = AI(K1,J)
            AR(K,J) = U1 * XR + U1I * XI + U2 * YR
            AI(K,J) = U1 * XI - U1I * XR + U2 * YI
            AR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            AI(K1,J) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(K,J)
            XI = BI(K,J)
            YR = BR(K1,J)
            YI = BI(K1,J)
            BR(K,J) = U1 * XR + U1I * XI + U2 * YR
            BI(K,J) = U1 * XI - U1I * XR + U2 * YI
            BR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            BI(K1,J) = U1 * YI + U1I * YR - U2 * XI
  180    CONTINUE
C
         IF (K .EQ. L) GO TO 240
         AI(K,KM1) = 0.0
         AR(K1,KM1) = 0.0
         AI(K1,KM1) = 0.0
C     ********** ZERO B(K+1,K) **********
  240    S = DABS(BR(K1,K1)) + DABS(BI(K1,K1)) + DABS(BR(K1,K))
         U1 = BR(K1,K1) / S
         U1I = BI(K1,K1) / S
         U2 = BR(K1,K) / S
         R = DSQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
         IF (K .EQ. NA) GO TO 245
         XR = AR(K2,K1)
         AR(K2,K1) = U1 * XR
         AI(K2,K1) = -U1I * XR
         AR(K2,K) = -U2 * XR
C
  245    DO 250 I = LOR1, K1
            XR = AR(I,K1)
            XI = AI(I,K1)
            YR = AR(I,K)
            YI = AI(I,K)
            AR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            AI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            AR(I,K) = U1 * YR - U1I * YI - U2 * XR
            AI(I,K) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(I,K1)
            XI = BI(I,K1)
            YR = BR(I,K)
            YI = BI(I,K)
            BR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            BI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            BR(I,K) = U1 * YR - U1I * YI - U2 * XR
            BI(I,K) = U1 * YI + U1I * YR - U2 * XI
  250    CONTINUE
C
         BI(K1,K1) = 0.0
         BR(K1,K) = 0.0
         BI(K1,K) = 0.0
         IF (.NOT. MATZ) GO TO 260
C
         DO 255 I = 1, N
            XR = ZR(I,K1)
            XI = ZI(I,K1)
            YR = ZR(I,K)
            YI = ZI(I,K)
            ZR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            ZI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            ZR(I,K) = U1 * YR - U1I * YI - U2 * XR
            ZI(I,K) = U1 * YI + U1I * YR - U2 * XI
  255    CONTINUE
C
  260 CONTINUE
C     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP **********
      IF (AI(EN,NA) .EQ. 0.0) GO TO 70
      R = CDABS(DCMPLX(AR(EN,NA),AI(EN,NA)))
      U1 = AR(EN,NA) / R
      U1I = AI(EN,NA) / R
      AR(EN,NA) = R
      AI(EN,NA) = 0.0
C
      DO 270 J = EN, ENORN
         XI = U1 * AI(EN,J) - U1I * AR(EN,J)
         AR(EN,J) = U1 * AR(EN,J) + U1I * AI(EN,J)
         AI(EN,J) = XI
         XI = U1 * BI(EN,J) - U1I * BR(EN,J)
         BR(EN,J) = U1 * BR(EN,J) + U1I * BI(EN,J)
         BI(EN,J) = XI
  270 CONTINUE
C
      GO TO 70
C     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT
C                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
 1000 IERR = EN
C     ********** SAVE EPSB FOR USE BY CQZVEC **********
 1001 IF (N .GT. 1) BR(N,1) = EPSB
      RETURN
C     ********** LAST CARD OF CQZVAL **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE CQZVEC(NM,N,AR,AI,BR,BI,ALFR,ALFI,BETA,ZR,ZI)
C
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN
      REAL*8 AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),
     X       BETA(N),ZR(NM,N),ZI(NM,N)
      REAL*8 R,T,RI,TI,XI,ALMI,ALMR,BETM,EPSB
      REAL*8 CDABS
      COMPLEX*16 Z3
      COMPLEX*16 DCMPLX
      REAL*8 DREAL,DIMAG
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER
C     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
C     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
C     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL
C          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS
C          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL,
C
C        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,
C
C        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT-
C
C        A IS UNALTERED,
C
C        B HAS BEEN DESTROYED,
C
C        ALFR, ALFI, AND BETA ARE UNALTERED,
C
C        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED
C          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IF (N .LE. 1) GO TO 1001
      EPSB = BR(N,1)
C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
      DO 800 NN = 2, N
         EN = N + 2 - NN
         NA = EN - 1
         ALMR = ALFR(EN)
         ALMI = ALFI(EN)
         BETM = BETA(EN)
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 700 II = 1, NA
            I = EN - II
            R = 0.0
            RI = 0.0
            M = I + 1
C
            DO 610 J = M, EN
               T = BETM * AR(I,J) - ALMR * BR(I,J) + ALMI * BI(I,J)
               TI = BETM * AI(I,J) - ALMR * BI(I,J) - ALMI * BR(I,J)
               IF (J .EQ. EN) GO TO 605
               XI = T * BI(J,EN) + TI * BR(J,EN)
               T = T * BR(J,EN) - TI * BI(J,EN)
               TI = XI
  605          R = R + T
               RI = RI + TI
  610       CONTINUE
C
            T = ALMR * BETA(I) - BETM * ALFR(I)
            TI = ALMI * BETA(I) - BETM * ALFI(I)
            IF (T .EQ. 0.0 .AND. TI .EQ. 0.0) T = EPSB
            Z3 = DCMPLX(R,RI) / DCMPLX(T,TI)
            BR(I,EN) = DREAL(Z3)
            BI(I,EN) = DIMAG(Z3)
  700    CONTINUE
C
  800 CONTINUE
C     ********** END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 2 DO -- **********
      DO 880 JJ = 2, N
         J = N + 2 - JJ
         M = J - 1
C
         DO 880 I = 1, N
C
            DO 860 K = 1, M
               ZR(I,J) = ZR(I,J) + ZR(I,K) * BR(K,J) - ZI(I,K) * BI(K,J)
               ZI(I,J) = ZI(I,J) + ZR(I,K) * BI(K,J) + ZI(I,K) * BR(K,J)
  860       CONTINUE
C
  880 CONTINUE
C     ********** NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1 **********
      DO 950 J = 1, N
         T = 0.0
C
         DO 930 I = 1, N
            R = CDABS(DCMPLX(ZR(I,J),ZI(I,J)))
            IF (R .GT. T) T = R
  930    CONTINUE
C
         DO 940 I = 1, N
            ZR(I,J) = ZR(I,J) / T
            ZI(I,J) = ZI(I,J) / T
  940    CONTINUE
C
  950 CONTINUE
C
 1001 RETURN
C     ********** LAST CARD OF CQZVEC **********
      END
