      program test_fdm1d_multiwindow
      implicit none
      integer Nsig, Nb0, Nbc, Nb, j, k, M, NWin, Nz, Nsp,
     &     i, n, Nsigmax, Nspmax, Nbmax
      parameters(Nsigmax=4896,Nspmax=10000,Nbmax=400)
      real*8 delt, beta(Nbmax), cheat, Gamm, 
     &     wmino, wmaxo, rho, ros,dW, tt,r1,r2
      complex*16 U(Nbmax,Nbmax,0:1), coef(0:Nsigmax), f(Nbmax), 
     &     z(Nbmax), g(Nbmax), diag(Nbmax), U0(Nbmax,Nbmax),
     &     coef1(Nsigmax), wk(Nbmax), dk(Nbmax), ZZ(Nbmax,Nbmax), 
     &     Spec(0:Nspmax), Spec_win(0:Nspmax)
      logical msflag,cheatmore
      integer ierror
      character*30 method

! Input:

      delt=0.628319E-02
      Nsig=2000
      Nsp=10000
      method='FDM'
      wmino=-100
      wmaxo=100
      rho=1.2
      Nb=50
      Nbc=0
      Gamm=0.0001
      
      open(1,file='../demo/1d/modelsignal')
      do n=0,Nsig
         read(1,*) tt,r1,r2
         coef(n)=dcmplx(r1,r2)
      enddo
      close(1)


      call fdm_parameters(Nsig,delt,wmino,wmaxo,rho,
     &     Nb0,Nbc,Nb,NWin,dW)

      write(6,*) 'Nsig=',Nsig
      write(6,*) 'delt=',delt
      write(6,*) 'wmino,wmaxo=',wmino,wmaxo
      write(6,*) 'Nwin=',Nwin, ' dW=', dW
      write(6,*) 'rho=',rho
      write(6,*) 'Nb0=',Nb0
      write(6,*) 'Nbc=',Nbc
      write(6,*) 'Nb=',Nb
      write(6,*) 'Gamm=',Gamm


       call fdm1d_multiwindow(
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


      open(2,file='../demo/1d/spectrum')
      do n=0,Nsp
         write(2,*) wmino+n*((wmaxo-wmino)/Nsp),dreal(Spec(n))
      enddo
      close(2)


      stop
      end
