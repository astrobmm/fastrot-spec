      program normspec

C     This program normalizes the spectrum created by spec3d

      implicit none

      integer i,j,np,nc,nsub
      integer idot
      
      real*4 w(50000),f(50000)
      real*4 wsub(5000),fsub(5000)
      real*4 c1(25),c2(25)
      real*4 wm,fm,wmax(25),fmax(25)
      real*4 fm1,fmn,fm2(25)
      real*4 wint,fint
      
      character*60 ifile,ofile,dummystr

      open(unit=1,file='input.spec3d',form='formatted',status='old')
      read(1,5) dummystr
      read(1,5) ifile
      close(unit=1)
 5    format(a60)
      
      i=1
      open(unit=1,file=ifile,form='formatted',status='old')
 1    read(1,*,end=1000) w(i),f(i)
      i=i+1
      goto 1

 1000 np=i-1
      close(unit=1)

      j=1
      open(unit=1,file='clips.normspec',form='formatted',status='old')
 2    read(1,*,end=1001) c1(j),c2(j)
      j=j+1
      goto 2
 1001 nc=j-1
      close(unit=1)

      do j=1,nc
         nsub=0
         do i=1,np
         if(w(i).gt.c2(nc)) goto 10
         if(w(i).gt.c1(j).and.w(i).le.c2(j)) then
               nsub=nsub+1
               wsub(nsub)=w(i)
               fsub(nsub)=f(i)
         end if
         end do
               call maximum(nsub,wsub,fsub,wm,fm)
               wmax(j)=wm
               fmax(j)=fm
      end do

 10   call dota(ifile,idot)
         
      ofile=ifile(1:idot-1)//'.n1'
      open(unit=2,file=ofile,form='formatted',status='unknown')
      
      fm1=(fmax(2)-fmax(1))/(wmax(2)-wmax(1))
      fmn=(fmax(nc)-fmax(nc-1))/(wmax(nc)-wmax(nc-1))

      call spline(wmax,fmax,nc,fm1,fmn,fm2)

      do i=1,np
         wint=w(i)
         call splint(wmax,fmax,fm2,nc,wint,fint)
         write(2,100) wint,f(i)/fint
 100  format(f11.5,2x,1pe13.6)         
      end do

      close(unit=2)

      write(*,6) ofile
 6    format(/,'Normalized spectrum: ',a60,/)
      
      stop
      end
      
      subroutine maximum(nsub,wsub,fsub,wm,fm)

      implicit none

      integer i,nsub

      real*4 wsub(nsub),fsub(nsub),wm,fm

      fm=-1.0e+6
      do i=1,nsub
         if(fsub(i).gt.fm) then
            fm=fsub(i)
            wm=wsub(i)
         end if
      end do
     
      return
      end

      subroutine dota(ifile,idot)

      implicit none

      integer i,idot

      character*60 ifile
      
      do i=1,60
         if(ifile(i:i+1).eq.'.a') then
            idot=i
            return
         end if
      end do

      return
      end
      
      subroutine spline(x,y,n,yp1,ypn,y2)

      implicit real*4 (a-h,o-z)
      parameter(nmax=30000)
      dimension x(n),y(n),y2(n),u(nmax)
      if (yp1.gt.0.99e30) then
          y2(1)=0.
          u(1)=0.
      else
          y2(1)=-0.5
          u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      end if

      do 11 i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*y2(i-1)+2.
          y2(i)=(sig-1.)/p
          u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  11  continue
      
      if (ypn.gt.0.99e30) then
          qn=0.
          un=0.
      else
          qn=0.5
          un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      end if

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      
      do 12 k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
  12  continue
      
      return
      end

      subroutine splint(xa,ya,y2a,n,x,y)
      
      implicit real*4 (a-h,o-z)                 
      dimension xa(n),ya(n),y2a(n)
      common/bracket/klo,khi

      klo=1
      khi=n

  1    if (khi-klo.gt.1) then
           k=(khi+klo)/2
           if (xa(k).gt.x) then
               khi=k
           else
               klo=k
           end if
       goto 1
       end if
      
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'Bad xa input!'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end

     
