c     Updated May 19-23, 2003 by Geoff Armstrong and Delaglio.
c-------------------------------------------------------------------------
c      program fdm1d
c-------------------------------------------------------------------------
c     Purpose: implement the 1D Filter Diagonalization Method (FDM). Given
c     a 1D signal data array, c(n) = C(tau*n), n=0, 1, ..., Nsig, compute 
c     the spectral parameter list for a specified local spectral window by
c     fitting the signal to:
c     
c            c(n) = sum_k d_k exp(-i n tau w_k)        (1)
c     
c     where w_k is the complex frequency (line position & width) and d_k is
c     the complex amplitude (phase & intregral). The nonlinear fitting 
c     problem is solved by converting it into a pure linear algebraic 
c     problem of diagonalizing a data matrix in a frequency domain subspace.
c
c     In current implementation, QZ-algorithm is applied to solve the 
c     generalized eigenvalue  problem.
c
c     Authors: Jianhan Chen (jianhanc@scripps.edu)
c              V. A. Mandelshtam (mandelsh@uci.edu)
c
c     References: 
c      1. H. Hu et. al., J. Magn. Reson. 134, 76-87 (1998).
c      2. J. Chen et. al., J. Chem. Phys. 112, 4429-4437 (2000).
c      3. V. A. Mandelshtam, Prog. NMR Spect. 38, 159-196 (2001).
c      4. QZ: http://www.netlib.org/toms/535
c-------------------------------------------------------------------------
c     Inputs and outputs: 
c       
c     The definitions of arguments of subroutine fdm1d() are given in the
c     the variable declaration part of the subroutine. The input paramters
c     will not be changed on return. The output paramters and local work
c     arrys are modified on return.
c
c     Additional information is given as following:
c     
c     1. coef, delt, Nsig: 
c             coef(n) = C(t) = C(n*delt), n=0,...,Nsig. This 1D arry defines
c             the 1D FID to be processed. 
c             
c             NOTE: The array must contain even number of data points, i.e.,
c                   Nsig = 2*M +1.
c             NOTE: The dwell time in 1/sec is actually given by delt/(2*pi).
c                   In another word, the value of delt feeded to the 1D FDM 
c                   module is the actual dwell time * (2 Pi).
c             
c     2. Nb0, Nbc, Nb: 
c             defines the number of fine, coarse, and total number of basis 
c             functions to be used. (Nb=Nb0+Nbc)
c
c     3. wminl, wmaxl:
c             defines the spectral range of the local window for FDM 
c             calculation: [wminl,wmaxl]
c
c             NOTE: number of fine grid window basis functions (Nb0) is 
c             related to the window range by,
c
c                    Nb0 = (wmaxl-wminl)/SW * (Nsig+1)/2 * rho, 
c                with,
c                    spectral width SW = 2*Pi/delt     
c                    basis density rho >= 1.0 (default should be 1.1).
c
c              Typical value is about 20-500, depending on the computational
c              power and the signal length. Large Nb0 means bigger spectral 
c              range per window. The result might be more stable but the will
c              be slower (t_comp ~ Nb^3). Except for very long signal (say, 
c              Nsig > 10,000), Nb0 = 50 ~ 150 will be a good choice.  
c
c             NOTE: number of coarse grid window basis functions (Nbc) is
c              typically in the range of 10 to 50.  Adding coarse basis 
c              (Nbc>0) will usually give more stable results and more 
c              accurate linelist, especially when the signal is noisy or 
c              contains non-localized features such as broad lines, strong
c              backgrounds. A good default coud be Nbc=20. Nbc=0 corresponds
c              to the original FDM with a single-scale basis. 
c
c     4. prflag: if .true., print some message on the console. (Removed: FD)
c     5. msflag: set to .false. 
c             
c              NOTE: if set to .true., use multi-scale FT basis instead of a
c               2-scale FT basis. It is usually not necessary unless ...
c
c     6. wk, dk: {wk(i), dk(i)|i=1,...,Nb} gives the paramter list.
c             outputs of the subroutine. complex frequencies and amplitudes.
c     
c     7. U, f, g, diag, z, coef1, zz, U0, beta: local work arrays.
c
c     8. ierror: returned, 0 on success.
c
c     FINAL NOTE: further information/example on how to set up the input 
c     parameters (Nb0,Nbc,Nb,wminl,wmaxl), see the testing code (test1d.f)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fdm1d(
c     input parameters
     &     coef, Nsig, delt, Nb0, Nbc, Nb, wminl, wmaxl, msflag,
c     output parameters
     &     wk, dk,
c     local work arrays
     &     U, f, g, diag, z, coef1, zz, U0, beta, 
c     error status returned
     &     ierror)
      implicit none
      integer Nsig, Nb0, Nbc, Nb, j, k
      real*8 wminl, wmaxl, delt, beta(Nb)
      complex*16 U(Nb,Nb,0:1),coef(0:Nsig),f(Nb),g(Nb),diag(Nb),z(Nb),
     &     coef1(Nsig),wk(Nb),dk(Nb),ZZ(Nb,Nb),U0(Nb,Nb)
      logical msflag
      integer ierror
c     
      ierror = 0
 
      call umatrix1d(coef, Nsig, 0, delt, Nb0, Nbc, Nb,
     &     wminl, wmaxl, f, U, msflag,
     &     g, diag, z, wk, dk, coef1,
     &     ierror)
c      
      call CQZ(Nb,U(1,1,1),U,ZZ,wk, U0, g, beta, ierror ) 
c
      do k=1,Nb                 ! get the frequency in the correct units
         wk(k)=(0,1d0)*cdlog(wk(k))/delt    
      enddo
      do k=1,Nb                 ! compute the amplitudes
         dk(k)=(0,0d0)
         do j=1,Nb
            dk(k)=dk(k) + ZZ(j,k)*f(j)
         enddo
         dk(k) = dk(k)**2
      enddo
c
      call cpiksrt(Nb,wk,dk)    ! sort according to fequency
      return
      end

c------------------------------------------------------------------------
c     Construction of U matrices in given basis. The multi-scale FT basis
c     is implemented in a more efficient way then the original algorithm
c     described in reference 2: instead of using various Mj to multi-scale 
c     basis, complex frequency grids are used, Im(zj) <-> Mj.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine umatrix1d(coef, Nsig, Nskip, delt, Nb0, Nbc, Nb,
     &     wminl, wmaxl, f, U, MULTI_SCALE,
c     local work arrays
     &     g, diag, z, ca1, ca2, coef1, 
     &     ierror )
      implicit none
      integer Nb,i,j,n,k,M,Nsig,Nskip,l,Nb0,Nbc,Mc
      complex*16 U(Nb,Nb,0:1),coef(0:Nsig),f(Nb),g(Nb),diag(Nb),z(Nb),
     &    ca1(Nb),ca2(Nb),coef1(Nsig),Z1,Z2,cc,xi
      real*8 delt,wminl,wmaxl,wminl_new,wmaxl_new,pi,z_re,z_im,dw,
     &     r1,rhoc,dw0,ros,Omega
      logical MULTI_SCALE
      integer ierror

c
c      MULTI_SCALE=.false.
c
c     NOTE:  set MULTI_SCALE = .false. if only 2-scale needed 
c     NOTE:  for multi-scale case, diffierent rhoc(j) distribution
c            can be chosen by inserting new formula for rhoc(j) at line 180
c     NOTE:  for typical application, set MULTI_SCALE=.false. as 2-scal
c            basis is more efficient and sufficient in most cases.
c
      xi=(0,1d0)
      pi=dacos(-1d0)
      M=(Nsig-1)/2
      wminl_new=wminl*delt
      wmaxl_new=wmaxl*delt
      dw0=(wmaxl_new-wminl_new)/Nb0
      z_im=0d0
      do i=1,Nb0                ! fine grid window basis
         z_re=wminl_new+(i-0.5d0)*dw0
         z(i+Nbc)=dcmplx(z_re,z_im)   
      enddo
      if(.not.MULTI_SCALE.and.Nbc.gt.0) then ! coarse grid window basis
         r1=dfloat(Nbc)/dfloat(M+1)
         call Im_Mj(r1,z_im,ierror) 
         z_im=z_im*delt
         dw=(2d0*pi-(wmaxl_new-wminl_new))/Nbc
         if(dw.le.dw0) then
c WARNING: coarse grid is too dense, reduce Nbc!
           ierror=1
         endif
         do i=1, Nbc
            z(i)=dcmplx(wminl_new-(i-0.5)*dw,z_im)
         enddo
      else if(Nbc.gt.0) then    ! multi-scale coarse basis
         do n=1, Nbc/2+mod(Nbc,2)
            r1=1-dfloat(2*n-1)/dfloat(Nbc+1)
 180        rhoc=(dsin(pi*r1/2.0))**2
c      rhoc=dfloat(Nbc)/dfloat(Ml-Nb0)        !alternative basis distributions
c      rhoc=dsin(pi*r1/2.0)**2*(3.0/Mmin)
c      rhoc=dsin(pi*r1/2.0)
            Mc=M*rhoc
            if(Mc.lt.4) Mc=4
            rhoc=dfloat(Mc)/dfloat(M)
            call Im_Mj(rhoc,z_im,ierror)
            z_im=z_im*delt
            if(n.eq.1) then
               z_re=dreal(z(Nbc+1))-dw0/rhoc
            else
               z_re=dreal(z(n-1))-dw0/rhoc
            endif
            z(n)=dcmplx(z_re,z_im)
            if(n.eq.1) then
               z_re=dreal(z(Nb))+dw0/rhoc
            else
               z_re=dreal(z(Nbc+2-n))+dw0/rhoc
            endif
            z(Nbc+1-n)=dcmplx(z_re,z_im)
         enddo
      endif
      do n=1, Nb
         z(n)=cdexp(z(n)*xi)
      enddo
c
c      calculating the 1D FT series
      call arrays(Nsig,M,coef(Nskip),Nb,z,f,g,diag, coef1)    
c
      do l=0,1
         do i=1,Nb
            Z1=z(i)**M
            Z2=Z1**2
            if(l.eq.0) then
               U(i,i,0)=diag(i)
               ca1(i)=f(i)
               ca2(i)=g(i)
            else 
               U(i,i,l)=( U(i,i,l-1)-coef(Nskip+l-1)
     &              - ca1(i) + ca2(i)*Z1 ) / z(i)
     &              + coef(Nskip+2*M+l) * Z2
               ca1(i)=ca1(i)/z(i)-coef(Nskip+l)+coef(Nskip+M+l)*Z1
               ca2(i)=ca2(i)/z(i)-coef(Nskip+M+l)+coef(Nskip+2*M+l)*Z1
            endif
         enddo
         do i=1,Nb
            Z1=z(i)**M
            do j=1,i-1
               Z2=z(j)**M
               cc=z(i)/z(j)
               U(i,j,l)=coef(Nskip+l)
     &              +(ca1(j)-cc*Z1*ca2(j)
     &              -cc*ca1(i)+Z2*ca2(i))/(1-cc)
               U(j,i,l)=U(i,j,l)
            enddo
         enddo
      enddo
      do i=1, Nb                ! actually gsave_Nskip 
         f(i)=f(i)+coef(Nskip)
      enddo
      do n=Nskip-1,0,-1         ! back-rotate gsave_Nskip to get f_0 
         do i=1, Nb
            f(i)=f(i)*z(i)-z(i)**(M+1)*coef(M+n+1)+coef(n)
         enddo
      enddo
      return
      end
c__________________________________________________________________________
c     Required by Umatrix1D_dMS1(): computing FT series for U matrices.
c
c        f(i) = sum_i=1^M z(i)**i * c(i)
c        g(i) = sum_i=1^M z(i)**i * c(i+M)
c
      subroutine arrays(Nsig,M,coef,Nz,z,f,g,diag, coef1)
      implicit none
      integer Nz,i,n,M,Nsig
      complex*16 coef(0:Nsig),f(Nz),g(Nz),diag(Nz),z(Nz),coef1(M),Z1
      do n=1,M
         coef1(n)=(n+1)*coef(n)
      enddo
      call Slow_FT(M,coef1(1),Nz,z,diag)
      do n=1,M
         coef1(n)=(n+1)*coef(n+M)
      enddo
      call Slow_FT(M,coef1(1),Nz,z,f)
      call Slow_FT(M,coef(M+1),Nz,z,g)
      do i=1,Nz
         Z1=z(i)**M
         diag(i)=diag(i)+coef(0)+((M+2)*g(i)-f(i))*Z1
      enddo
      call Slow_FT(M,coef(1),Nz,z,f)
      return
      end
c----------------------------------------------------------------------------
c     Required by arrays(): a more accurate slow FT 
c
      subroutine Slow_FT(M,coef,Nz,z,f)
      implicit none
      integer M,n,Nz,i
      complex*16 coef(M),z(Nz),f(Nz),fff
      do i=1,Nz
         f(i) = (0,0d0)
         fff = coef(M)*z(i)
         do n=M-1,1,-1
            fff = (fff+coef(n))*z(i)
            if(mod(n,100).eq.0.or.n.eq.1) then
               f(i)=f(i)+fff*z(i)**(n-1)
               fff=(0,0d0)
            endif
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------------
c     Given t=(Mj+1)/(M+1), calculate z_im which satisfies: 
c        sum_{n=0}^M exp(-z_im*n)= Mj+1
c
      subroutine Im_Mj(t, Im, ierror) 
      implicit none
      real*8 t, x, Im
      integer ierror

      if(t.le.0.or.t.gt.1.0) then
c WARNING t must be between 0 and 1.0')
         ierror=2
         t=1d-2
      endif
      if(t.le.0.685) then
         x=0d0
         if(t.gt.1d-2) x=-dexp(-1d0/t)/t
         Im=1d0/t+x-x**2+1.5d0*x**3-8d0*x**4/3d0+125d0*x**5/24d0
      else
         Im=2.82d0*(1d0-t)
      endif
      return
      end
c------------------------------------------------------------------------
c     The frequencies are sorted in increasing order.
c
      subroutine cpiksrt(Nb,wr,wi)
      implicit real*8(a-h,o-z)
      complex*16 wr(Nb),wi(Nb),war,wai
      do j=2,Nb
         war=wr(j)
         wai=wi(j)
         do i=j-1,1,-1
            if(dreal(wr(i)).le.dreal(war)) go to 10
            wr(i+1)=wr(i)
            wi(i+1)=wi(i)
         enddo
         i=0
 10      continue
         wr(i+1)=war
         wi(i+1)=wai
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Name   : CQZ (Nb, UL, UR, ZZ, uk)
c     Purpuse: driver for solving the GEP, UL Bk = uk UR Bk
c     Inputs : Nb define the dimension of UR and UL
c              on return UL is destroyed, UR is not changed.
c     Ouputs : ZZ contains the normalized eigenvector
c              uk contains the eigenvalues
c     Author : Jianhan Chen (jianhanc@uci.edu)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine CQZ(Nb, UL, UR, ZZ, uk,
c     local arrays
     &     U0, alf, beta,
     &     ierror )
      implicit none
      integer Nb, i, j, k, n, ierr
      complex*16 UL(Nb,Nb), UR(Nb,Nb), U0(Nb,Nb), ZZ(Nb,Nb),
     &     uk(Nb), alf(Nb), cc, uu
      real*8 beta(Nb)
      integer ierror
c
      do i=1,Nb
         do j=1,Nb
            U0(i,j)=UR(i,j)
         enddo
      enddo
      call CQZHES(Nb,Nb,UL,U0,.true.,ZZ)
      call CQZVAL(Nb,Nb,UL,U0,0d0,alf,beta,.true.,ZZ,IERR)
      if(IERR.ne.0) then
         ierror=3
      endif
      call CQZVEC(Nb,Nb,UL,U0,alf,beta,ZZ)
c
      do k=1,Nb
         if(beta(k).le.1d-8) then
            uk(k)=cdexp(dcmplx(-10,-10))
            do j=1, Nb
               ZZ(j,k)=(0,0d0)
            enddo
         else
            uk(k)=alf(k)/beta(k)
            cc=(0,0d0)
            do j=1,Nb
               uu=0.0d0
               do n=1,Nb
                  uu=uu+UR(j,n)*ZZ(n,k)
               enddo
               cc=cc+ZZ(j,k)*uu
            enddo
            cc=(1d0,0)/cdsqrt(cc)
            do j=1,Nb
               ZZ(j,k)=cc*ZZ(j,k)
            enddo
         endif
      enddo
      return
      end
c---------------------------------------------------------------------------
c     END OF PROGRAM (fdm1d.f)
c---------------------------------------------------------------------------
