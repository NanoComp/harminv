c     modified for testing the 1D FDM module: prepare the FDM input
c     parameters and work arrays, call fdm1d() subroutine, and generate
c     a 1D spectrum.     Jianhan Chen, 01/19/2003
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        program fdm1d_multiwindow
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1. NOTE:    USE OR DISTRIBUTE WITH PERMISSION ONLY.
c
c     2. NAME:    fd_rrt1d.f
c
c     3. AUTHORS: Jianhan Chen (jianhanc@uci.edu)
c                 V. A. Mandelshtam (mandelsh@uci.edu)
c
c     4. DATE:    10-10-2001
c
c     5. PURPOSE: 1D FDM/RRT/DFT
c     
C              FDM = Filter Diagonalization Method
c              RRT = Regularized Resolvent Transform 
c              DFT = Discrete Fourier Transform
C
c     6. GENERAL DESCRIPTION: 
c
c      The program is based on the FDM and RRT algorithms for high resolution 
c      spectral analysis of time domain signals. See references for details.
c      Given a signal c(n) = C(t0+tau*n) n=0, 1, ..., Nsig, calculate the 
c      line list {w_k, d_k} and/or the complex spectrum F(w) by fitting c(n) 
c      to the Auto-Regressive (AR) model,
c
c                 c(n) = sum_k d_k exp(-i n w_k)        (1)
c                 
c      where, w_k is the complex frequency (line position & width)
c            d_k is the complex amplitude 30482(phase & intregral)
c     
c      1) in FDM, the complex spectrum is computed using the formula,
c     
c                 F(w) = sum_k d_k/(w-w_k+i*Gamm)       (2)
c
c         in RRT, it is estimated as,
c
c                 F(w) ~ Cv 1/R(w) Cv'                    (3)
c     
c      2). QZ-algorithm is applied to solve the generalized eigenvalue 
c      problem. ZGESV routine from LAPACK is used to invert matrix R(w).
c
c     7. INPUTS:
C
C      IT IS IMPORTANT THAT YOU UNDERSTAND SOME BASIC IDEAS OF FDM/RRT
C      IN ORDER TO UNDERSTAND THE MEANING OF SOME INPUTS SUCH AS BASIS 
C      SIZE, BASIS DENSITY. IN SUMMARY, FDM/RRT PROCESS THE TIME SIGNAL
C      BY DIVIDING THE WHOLE NYQUIST RANGE INTO ONE OR SEVERAL LOCAL
C      WINDOWS IN THE FREQUENCY DOMAIN. 
c
c      signal: input data file with appropriate first line. The data 
c              file is assumed to contain an ASCII formatted list of
c              c(n\tau). The first line must be as following:
c
c               dim { NSigmax tau }_[1,..,dim] {idat}_[1,...,dim]
c              
c              where dim    : the dimension of the signal (dim=1, 2, ...)
c                    Nsigmax: the maximum signal length
c                    tau    : time step (acquistion inteval)
c                    idat   : data format along each dimension
c               e.g/ 1D example: 1 1024 0.001 1 
c                    2D example: 2 2024 0.001 128 0.0001 1 1
c              
c              Note: 'idat' implies the data tormat. In 1D, it is defined
c              as following:
c                     idat>0 read real signal, idat<0 read complex signal 
c                     |idat| = 1 read c_n
c                     |idat| = 2 read t_n,c_n
c
c      t0    : 1st order (i.e., linear) phase correction (in seconds) 
c              NOTE: t0 is not used in this code.
c      theta : zero order (overall) phase correction
c      method: method used to compute the spectra: FDM, RRT, or DFT
c      Nsig  : number of data points to be used in computing the spectra.
c              It will be automatically reduced to Nsigmax (maximum length)
c              if specified Nsig is greater than Nsigmax.
c      wmino : (see wmaxo)
c      wmaxo : specify the frequency range of the spectra to be computed. 
c              It might be adjusted to fit harmonically with the basis size
c              specified (see Nb0, rho) and to fit in the Nyquist Range.
c
c      threshhold: see par below.
c      par   : file where the linelist will be written (FDM only).
c              Only those with |d(k)| > threshhold will be output.
c      ReSp  : file where Re[F(w)] will be written 
c      ImSp  : file where Im[F(w)] will be written
c      AbsSp : file where Abs[F(w)] will be written
c              NOTE: if file name is specified as 'none' or 'None', no
c                    corresponding output will be written.

c      rho   : basis density (default = 1.1). rho > 1.0 means putting 
c              more basis functions than default. For example, rho=1.1
c              means putting 10% more basis functions then default. 
c              rho should not be less than 1.0. rho = 1.5 should be
c              considered as the maximum. For typical data, rho = 1.1-1.2
c              will usually work the best. if rho < 0, default is used.
c      Nb0   : number of narrow band Fourier basis functions used per 
c              window (see references). Typical value is about 50-300.
c              Large Nb0 means bigger spectral range per window. The result
c              might be more stable but the calculation will be slower for
c              larger Nb0. Except for very long signal (say, Nsig > 10,000),
c              Nb0 = 150 will be a good choice. 
c      Nbc   : number of broad band Fourier basis functions (coarse basis).
c              See second reference for more details. Nbc <= 0 means no
c              coarse basis. Adding coarse basis (Nbc>0) will usually give
c              stabler results and more accurate linelist, especially when
c              the signal is noisy or contains non-localized features such
c              as broad lines, backgrounds. Typical Nbc = 10 ~ 50. A good
c              choice could be Nbc = 20.
c
c      Nsp   : number of points used in plotting the spectra. Defines the
c              digital resolution of output spectra.
c      Gamm  : smoothing parameter. Defines the smallest possible linewidth.
c              Default Gamm = 0.2 * 1/(N*tau). FDM/RRT might give some very
c              narrow lines (spikes) due to noise. A non-zero Gamm will 
c              improve the looking of resulted spectra. if Gamm < 0, default
c              will be used.
c      cheat : multiply all widths by cheat (FDM only). 
c              NOTE: be very careful to use cheat < 1.0. It might lead to 
c                    very misleading results. Use cheat=1.0 only.
c      cheatmore: if is .true. F(w) is computed with Im d_k (FDM only).
c              NOTE: be very careful again. Default: cheatmore = .false.
c
c      ros   : regularization parameter (RRT only). If ros < = 0, no 
c              regularization in RRT. See 3rd reference for more details.
c              In general, bigger ros means more 'regularization', and will
c              lead to more stable result but will lower resolution.
c              Qualitively, ros is related to the amount of noise present
c              in the signal. Noisy data requires bigger regularization.
c              For typical NMR signals with fine S/N (say, 500M magnet,
c              1 mM sample, normal probe, 4-16 scans), ros = 1d-6 should
c              a good guess. If you see artifacts (very easy to tell if 
c              you know what a typical NMR peak will look like), increase
c              ros (say, double the value) and run again. If the spectrum
c              look very clean but the resolution is very poor, it may 
c              indicate too big regularization. Decrease ros and run again.
c              NOTE: Choosing an optimal value for ros is a try-and-error
c                    game. Typically, there will be a flat regime where
c                    the spectra remains similar even changing ros by order 
c                    of a magnitude. Try following:
c
c                              too many artifacts
c                  ros=1d-6---------------  ros = ros*10 = 1d-5 --> cont.
c                            | clean, but low res. 
c                            -------------  ros = ros/10 = 1d-7 --> cont.
c                            | GREAT! 
c                            --------- done!
c
c      *** SAMPLE INPUT FILE: TO PROCESS A SIGNAL STORED IN 'signal.txt'
c
c       'signal.txt'        	                /signal
c       1.57                     		/theta
c       1                               	/ispec
c       1000	                            	/Nsig
c       -10 10 	                  		/wmin wmax
c       'par'	 	                      	/parameters file
c       'ft','none','none' 	          	/ReSp,ImSp,AbsSp files
c       1., 200, -20                          	/rho, Nb0, Nbc
c       3000, 1d-2				/Nsp, Gamm
c       1 F                             	/cheat, cheatmore
c       5d-7                             	/ros
c
c     8. OUTPUTS: can be found in the spectra and parameter files specified
c        in the input file. When multi-scale is used, the position of basis
c        functions can be found in file 'fort.11'.
c
c     9. REFERENCES
c     
c      FDM/RRT:
c       H. Hu et. al., J. Magn. Reson. 134, 76-87 (1998).
c       J. Chen and V. A. Mandelshtam, J. Chem. Phys. 112, 4429-4437 (2000). 
c       J. Chen, et.al., J. Magn. Reson. 147, 129-137 (2000). 
c       V. A. Mandelshtam, Prog. Nuc. Magn. Reson. Spect. 38, 159-196 (2001). 
c      QZ: http://www.netlib.org/toms/535
c      ZGESV/LAPACK: http://www.netlib.org/lapack/index.html
c
c     Questions and comments should be directed to V. A. Mandelshtam.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fdm1d_multiwindow(
c     input parameters
     &     coef, Nsig, delt, method, rho, Nb0, Nbc, Nb, wmino, wmaxo, 
     &       msflag, Nsp, Gamm, cheat, cheatmore, ros
c     output 
     &     Spec(0:Nsp),
c     local work arrays
     &     wk, dk, U, f, g, diag, z, coef1, zz, U0, beta, 
     &       Spec_win,
c     error status returned
     &     ierror)
      implicit none
      integer Nsig, Nb0, Nbc, Nb, j, k, M, NWin, ispec, Nsp,
     &     iminl, imaxl, nw, N0, i, n
      real*8 wminl, wmaxl, delt, beta(Nb), cheat, pi, Omega, Gamm, 
     &     wmin, wmax,
     &     wmino, wmaxo, rho, dOmega, dW, wmeanl,
     &     tt, r1, ros
      complex*16 U(Nb,Nb,0:1), coef(0:Nsig), f(Nb), g(Nb), diag(Nb), 
     &     z(Nb),
     &     coef1(Nsig), wk(Nb), dk(Nb), ZZ(Nb,Nb), U0(Nb,Nb),
     &     Spec(0:Nsp), cc, dd, uu, Z1, Z2, Z3, xi, Spec_win(0:Nsp)
      logical msflag,cheatmore
      integer ierror
      character*30 method
c
c      write(*,001)
 001  format(/,'  1D FDM/RRT/DFT with Multi-Scale option, Beta version',
     &     /,  '       J. Chen and V. A. Mandelshtam, UC-Irvine',/)
c
      xi=(0d0,1d0)
      pi=dacos(-1d0)
      r1=0d0                    ! r1 = sum_n |c(n)|^2
      do n=0,Nsig
         r1=r1+cdabs(coef(n))**2
      enddo
      M=(Nsig-1)/2
c      write(6,9991) 2*M+1
      ros=ros*r1*(M+1d0)**2     ! actual regularization parameter
c
      wmin=wmino*delt
      wmax=wmaxo*delt
      if(wmin.ge.wmax.or.dabs(wmin).gt.pi.or.dabs(wmax).gt.pi)
     &     stop 'Error: invalid spectral range, abort!'
c      write(6,9992) wmino,wmaxo
 9992 format(' Spectral range to be processed: [',F10.3,', ',F10.3,']')
c      if(Gamm.lt.0) write(6,*) 'Warning: Gamm < 0,',
c     &     ' which means line narrowing instead of smoothing!'
      if(Gamm.le.0) Gamm = 0.2*pi/(M*delt)
c      write(*,*) 'Smoothing factor Gamm = ', Gamm
      dOmega=(wmaxo-wmino)/Nsp
      do j=0,Nsp
         Spec(j)=(0d0,0)
      enddo
c
c      write(6,*)
      if(method.eq.'DFT') then
         ispec = -1
c         write(6,*) 'Discrete Fourier Transform'
         goto 246               ! DFT
      else if(method.eq.'RRT') then
         ispec = 0
c         write(6,*) 'Regularized Resolvent Transform (RRT)'
         stop 'Error: RRT not tested in current code'
      else 
         ispec = 1
c         write(6,*) 'Filter Diagonalization Method (FDM)'
      endif
c
CC      dW = 2*(wmaxo-wmino)/(NWin+1)   ! window size
CC      Nb0 = dW/2*delt/pi * (M+1)*rho ! adjusted size of window basis
c
c      if(ispec.eq.0) then
c         write(6,*) 'Regularization paramter q^2 = ', ros         
c      else if(ispec.gt.0) then
c         if(cheat.ne.1d0) write(6,*) 'using cheat = ', cheat
c         if(cheatmore) write(6,*) 'cheatmore = .true.,',
c     &        ' |dk| is used to construct F(w)'
c      endif
c_______________________________________________________________
c
c      if(checkfile(par).and.ispec.gt.0) then
c         threshhold = max(0d0, threshhold)
c         open(7,file=par)
c         write(7,8880) threshhold
 8880    format('# Parameter files generated by fd_rrt1d.f',/,
     &        '# Only lines with  |d(k)| > ',F10.7, ' are recorded',/
     &        '# w: complex frequencies, d: complex amplitudes',/,/,
     &        '#',7x,'Re w',18x,'Im w',15x,'Re d',14x,'Im d')
c         write(6,8881) threshhold, par
 8881    format(' {wk,dk} for lines with |d(k)| > ',F8.6, 
     &        ' will be writen to : ',A14)
c      endif
c      write(6,*)
      wmeanl=wmino
      do nw=1, NWin           ! multi-window implementation
         wminl=wmeanl
         wmeanl=wminl+dW
         wmaxl=wmeanl+dW
c         write(6,9995) nw, wminl, wmaxl
 9995    format('  window No.',I4,': [wmin,wmax]=[',F10.3,',',F10.3,']')
         iminl=(wminl-wmino)/dOmega
         imeanl=(wmeanl-wmino)/dOmega
         imaxl=(wmaxl-wmino)/dOmega
         do i=iminl, imaxl
            Spec_win(i) = dcmplx(0,0)
         enddo
c
c     call fdm1d() module 
c
         call fdm1d(
     &     coef, 2*M+1, delt, Nb0, Nbc, Nb, wminl, wmaxl, .true.,.true.,
     &     wk, dk,
     &     U, f, g, diag, z, coef1, zz, U0, beta)
c
c     computing the 1D spectrum
c
         do k=1,Nb              ! F(w) ~ sum_k {...}
            dd=dk(k)
            cc=wk(k)
            if(dabs(dimag(cc)).lt.Gamm) then
               cc=dcmplx(dreal(cc),-Gamm)
            else
               cc=dcmplx(dreal(cc),cheat*dimag(cc))
            endif                  
            do i=iminl, imaxl
                  Omega=wmino+i*dOmega
                  if(cheatmore) then
                     Spec_win(i)=Spec_win(i)-xi*delt*dd*
     &                    dreal(1/(1-cdexp(xi*delt*(Omega-cc))))
                  else
                     Spec_win(i)=Spec_win(i)-xi*delt*dd/
     &                    (1-cdexp(xi*delt*(Omega-cc)))
                  endif
               enddo
c               if(checkfile(par)) then
c                  if(dreal(wk(k)).ge.wminl.and.dreal(wk(k)).le.wmaxl
c     &                 .and.cdabs(dk(k)).ge.threshhold) 
c     &                 write(7,13) dreal(wk(k)),dimag(wk(k)),
c     &                 dreal(dk(k)),dimag(dk(k))
c               endif
 13            format(E22.16,3E18.10)
            enddo
c
c     add up the indivual contributions from each window
c     
            do i=iminl, imaxl
               Omega=wmino+i*dOmega
               if(iminl==0.and.Omega<0.5*(wmaxl+wminl).or.          ! the left edge 
     &             imaxl.ge.Nsp-1.and.Omega>0.5*(wmaxl+wminl) then  ! the right edge
                  Spec(i)=Spec(i)+Spec_win(i)
               else
                  Spec(i)=Spec(i)+Spec_win(i)*
     &                 dsin(pi*(Omega-wminl)/(wmaxl-wminl))**2
               endif
            enddo
         endif
      enddo
c      
 246  if(ispec.lt.0) then       ! DFT
         N0=Nsig/10
         do n=0,Nsig
            if (n.gt.N0) then
               r1 = dfloat(n-N0)/dfloat(Nsig+1-N0)
               coef(n)=coef(n)*dexp(r1)*(1D0-r1)
            endif
         enddo
         do i=0,Nsp
            Omega=wmino+i*dOmega
            uu=cdexp((0,1d0)*Omega*delt)
            Z1=-delt*xi         !*cdexp(dcmplx(0d0,t0*Omega))
            Spec(i)=Spec(i)+coef(0)*Z1/2
            do n=1,Nsig
               Z1=Z1*uu
               Spec(i)=Spec(i)+coef(n)*Z1
            enddo
         enddo
      endif
c
      if(ispec.ne.-1) then      ! first point correction
         cc=-coef(0)*delt*xi/2.0
         do i=0, Nsp
            Spec(i)=Spec(i)-cc
         enddo
      endif
c
c      write(6,*)                ! output spectra
c      if(checkfile(ReSp)) then
c         write(6,*) 'write Re[I(w)] to file: ', ReSp
c         open(14,file=ReSp)
c         do i=0,Nsp
c            write(14,986) wmino+i*dOmega,dreal(Spec(i))
c         enddo
c         close(14)
c      endif
c      if(checkfile(ImSp)) then
c         write(6,*) 'write Im[I(w)] to file: ', ImSp
c         open(14,file=ImSp)
c         do i=0,Nsp
c            write(14,986) wmino+i*dOmega,dimag(Spec(i))
c         enddo
c         close(14)
c      endif
c      if(checkfile(AbsSp)) then
c         write(6,*) 'write Abs[I(w)] to file: ', AbsSp
c         open(14,file=AbsSp)
c         do i=0,Nsp
c            write(14,986) wmino+i*dOmega,cdabs(Spec(i))
c         enddo
c         close(14)
c      endif
c 986  format(E14.7,E12.4)
c      write(6,*)
c      write(6,*) 'End of program: successful excution'
      return
      end
c      
c      logical function checkfile(file)
c      character*30 file
c      checkfile=.true.
c      if(file.eq.'none'.or.file.eq.'NONE') checkfile=.false.
c      return
c      end

      subroutine fdm_parameters(Nsig,delt,wmino,wmaxo,rho,
     &     Nb0,Nbc,Nb,NWin,dW)
!
!         Input:
!
!   [wmino;wmaxo] - big frequency window
!   rho - basis density
!   Nb0 - small window basis size  
!   Nbc - course basis size
!   delt - time step
!   Nsig - signal array is coef(0:Nsig)
!
!         Output:
!
!   [wmino;wmaxo] (truncated to the Nyquist range, if too big)
!   rho, Nb0, Nbc, Nb=Nb0+Nbc (adjusted to be consistent with the multiwindow calculation)
!   Nwin - number of small overlapping windows
!   dW - size of a small window (coinsides with wmaxo-wmino if Nwin=1)
!   Nsig=2*M+1  (truncated by one if even)
!   
!
! To force a single window calculation in [wmino;wmaxo] set Nb0=100000
!
!
      implicit none
      parameter(Nbmax=1000)
      real*8 delt,wmino,wmaxo,rho,wmin,wmax,dW
      integer Nb0,Nbc,Nb,NWin,M,Nsig,Nz
      if(rho<0d0) then  ! Use default
         rho=1.1
         Nb0=100
         Nbc=0
      endif
      pi=dacos(-1d0)
      wmin=wmino*delt
      wmax=wmaxo*delt
      if(wmax.le.wmin) stop 'error: wrong window parameters'
      if(wmax>pi) then
         wmax=pi
         wmaxo=wmax/delt
      endif
      if(wmin<pi) then
         wmin=-pi
         wmino=wmin/delt
      endif
      M=(Nsig-1)/2
      Nsig=2*M+1
      Nz=(wmaxo-wmino)/((2*pi)/delt)*(M+1)*rho
      if(Nb0<2) Nb0=2
      if(Nb0>M+1) Nb0=M+1
      if(Nb0>Nz) Nb0=Nz
      NWin = (2*Nz)/Nb0-1         ! number of small windows
      dW = (wmaxo-wmino)/NWin     ! small window size
      Nb0 = dW/(2*pi/delt) * (M+1)*rho ! adjust 
      rho=Nb0*(2*pi/delt)/dW/(M+1d0)+1d-6  ! adjust
      if(Nbc+Nb0>M+1) Nbc=M+1-Nb0
      Nb=Nbc+Nb0
c      write(6,9993) NWin, Nb0
 9993 format(' Divide the spectral range into ',I4,' windows with Nb0 ='
     &     , I4)
c      write(6,9994) Nb
 9994 format(' The actual basis size per window is Nb = ',I4)
      return
      end
