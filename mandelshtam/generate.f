      real*8 pi
      parameter (K=50,tau=1E-3,N=65536,imax=16384)
      complex*16 w(3*K),d(3*K),coef(0:N),Z1,Z
      integer*4 ir,ii
      pi=dacos(-1d0)
      delt=tau*(2*pi) !/2
      w(K)=(1,-0.001d0)/tau/2
      w(2*K)=(0.995,-0.001d0)/tau/2
      w(3*K)=(0.990,-0.001d0)/tau/2
      d(K)=(1d0,0)
      d(2*K)=(2d0,0)
      d(3*K)=(1d0,0)
      do j=K-1,1,-1
         w(j)=w(j+1)*0.9
         d(j)=d(j+1)   !*0.9
      enddo
      do j=2*K-1,K+1,-1
         w(j)=w(j+1)*0.9
         d(j)=d(j+1)   !*0.9
      enddo
      do j=3*K-1,2*K+1,-1
         w(j)=w(j+1)*0.9
         d(j)=d(j+1)   !*0.9
      enddo
      open(2,file='param')
      open(1,file='signal')
      do j=1,K
         write(2,20) real(w(2*K+j)),imag(w(2*K+j)),
     &        real(d(2*K+j)),imag(d(2*K+j))
         write(2,20) real(w(K+j)),imag(w(K+j)),
     &        real(d(K+j)),imag(d(K+j))
         write(2,20) real(w(j)),imag(w(j)),real(d(j)),imag(d(j))
      enddo
 20   format(4E15.7)
      idim=1 
      write(1,*) idim,N,delt
      do n1=0,N
         coef(n1)=(0,0d0)
      enddo
      do j=1,3*K
         Z=cdexp(-(0,1d0)*w(j)*delt)
         Z1=d(j)/Z
         do n1=0,N
            Z1=Z1*Z
            coef(n1)=coef(n1)+Z1
         enddo
      enddo
      Z1=imax/coef(0)
      do n1=0,N
         cr=real(coef(n1)*Z1)
         ir=abs(cr)
         if(cr.lt.0E0) ir=-ir
         ci=imag(coef(n1)*Z1)
         ii=abs(ci)
         if(ci.lt.0E0) ii=-ii
         write(1,10) n1*delt,ir,ii
      enddo
 10   format(E14.6,2I8)
      end







