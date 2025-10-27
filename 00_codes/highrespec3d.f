      program highrespect3d

      implicit none
      
      integer i,j,k,ib
      integer nlines
      integer nreg
      integer ilat,ilats,nlat
      integer npts
      integer ivsini
      
      real*4 pi,c,gamma
      real*4 xr,yr,ith,ithold,areapn,vp
      real*4 c1,c2
      real*4 mu,xrmax,yrmax
      real*4 limbdark,a1,a2,a3,a4
      real*4 ld1(5),ld2(5),ld3(5),ld4(5)
      real*4 ld1w,ld2w,ld3w,ld4w
      real*4 lat(200),t(200),lg(200)
      real*4 mint,maxt,ming,maxg
      real*4 wlat(200,45000),flat(200,45000)
      real*4 kt,kg
      real*4 temp,logg
      real*4 w4(45000),f4(45000)
      real*4 w(45000),f(45000),wmin,wden,wratio
      real*4 wvp(45000),fareapn(45000)
      real*4 fa1,fan,f2n(45000),wint,fint
      real*4 fmapped(45000),ftot(45000)

      real*4 dum1
      
      character*1  sfilter
      character*3  kcode
      character*4  metal
      character*10 m(4)
      character*45 stin,stdum,sttpole,stlgpole,strpole,stom    
      character*60 outspec,koutspec
      character*65 s(1000),string(1000)
      
      parameter(c=299792.458)

      pi=acos(-1.0)
     
      open(unit=3,file='input.spec3d',form='formatted',status='old')
      read(3,101) metal
 101  format(a4)
      read(3,102) outspec
 102  format(a60)
      read(3,103) nlat
 103  format(i3)
      close(unit=3)

      if(metal.eq.'-2.5') kcode='m25'
      if(metal.eq.'-2.0') kcode='m20'
      if(metal.eq.'-1.5') kcode='m15'
      if(metal.eq.'-1.0') kcode='m10'
      if(metal.eq.'-0.5') kcode='m05'
      if(metal.eq.'+0.0') kcode='p00'
      if(metal.eq.'+0.2') kcode='p02'
      if(metal.eq.'+0.5') kcode='p05'                        
      
      write(*,99)
 99   format(' ')

C     spec3d.lat:  Here we read the four HIGH-RESOLUTION models 
C
C                  (T_low,g_low)   (T_high,g_low)
C                  (T_low,g_high)  (T_high,g_high)
C      
C     corresponding to each latitude, that bracket the values of (T,g),
C     and the constants c1, c2:
C
C     (T,g_low)  = c1*(T_low,g_low)  + (1-c1)*(T_high,g_low)
C     (T,g_high) = c1*(T_low,g_high) + (1-c1)*(T_high,g_high)      
C
C     (T,g) = c2*(T,g_low) + (1-c2)*(T,g_high)

      do i=1,nlat
        lat(i)=0.0
        t(i)=0.0
        lg(i)=0.0
      end do

      mint=1.0e+6
      maxt=-1.0e+6
      ming=6.0
      maxg=-6.0
      
      open(unit=1,file='spec3d.lat',form='formatted',status='old')      
      do ilat=1,nlat/2
      read(1,11,end=1100) lat(ilat),t(ilat),lg(ilat),(m(j),j=1,4),c1,c2
 11   format(1x,f6.2,1x,f7.1,1x,f5.3,4(1x,a10),1x,2(1x,f6.4))
      
      write(*,12) abs(lat(ilat)),int(t(ilat)),lg(ilat)
 12   format('Building Kurucz HIGH-RESOLUTION models for th=+-',
     +        f5.2,' T=',i5,' log g=',f4.2)
      
      ilats=nlat-(ilat-1)
      lat(ilats)=-lat(ilat)
      t(ilats)=t(ilat)
      lg(ilats)=lg(ilat)

      mint=min(mint,t(ilat))
      maxt=max(maxt,t(ilat))
      ming=min(ming,lg(ilat))
      maxg=max(maxg,lg(ilat))
      
C     HIGH-RESOLUTION models
      
         call mod4lat(kcode,m,c1,c2,npts,w4,f4)
         do k=1,npts
            wlat(ilat,k)=w4(k)
            flat(ilat,k)=f4(k)
            wlat(ilats,k)=w4(k)         
            flat(ilats,k)=f4(k)      
         end do
      end do
      
 1100 close(unit=1)
      
      write(*,99)

C     Writing a subset of the data for the limb-darkening coefficients
C     using only the specific max/min temperatures and gravities
C     of the star.

      call subsetld(kcode,mint,maxt,ming,maxg,nlines,s)
      do i=1,nlines
         string(i)=s(i)
      end do
      
C     ------------------------------------------------
C     HIGH RESOLUTION SPECTRUM BETWEEN 3600 AND 5500 A
C     ------------------------------------------------
      
C     ---- HI-RES spectrum: initializing total flux
      
      do i=1,npts
         ftot(i)=0.0
      end do

      ithold=-999.0

C     tiles.a:  the file containing information on the VISIBLE tiles
C     is read out. Relevant parameters are the latitude of each 
C     surface element, the projected area and the velocity shift to be
C     applied to each individual spectrum.
C
C     [dum1 is a dummy variable, not used]
  
      open(unit=3,file='tiles.a',form='formatted',status='old')
 3    read(3,30,end=300) xr,yr,ith,dum1,areapn,vp,mu
 30   format(7(1pe14.6,2x))

C     The HI-RES fluxes w(ilat),f(ilat,*) are assigned to the tile      
C     according to its latitude

      if(ith.eq.ithold) goto 5
      write(*,31) ith
 31   format('Computing models for individual tiles at latitude=',
     +        f6.2,' deg')

      do ilat=1,nlat
         if(ith.eq.lat(ilat)) then
            temp=t(ilat)
            logg=lg(ilat)
            do i=1,npts
              w(i)=wlat(ilat,i)
              f(i)=flat(ilat,i)
            end do
C
C     Computing the limb darkening coefficients for the (T,log g)
C     corresponding to each latitude
C
           sfilter='U'
           call ld(nlines,string,temp,logg,a1,a2,a3,a4,sfilter)
           ld1(1)=a1
           ld2(1)=a2
           ld3(1)=a3
           ld4(1)=a4
           sfilter='B'
           call ld(nlines,string,temp,logg,a1,a2,a3,a4,sfilter)
           ld1(2)=a1
           ld2(2)=a2
           ld3(2)=a3
           ld4(2)=a4      
           sfilter='V'
           call ld(nlines,string,temp,logg,a1,a2,a3,a4,sfilter)      
           ld1(3)=a1
           ld2(3)=a2
           ld3(3)=a3
           ld4(3)=a4
           sfilter='R'
           call ld(nlines,string,temp,logg,a1,a2,a3,a4,sfilter)      
           ld1(4)=a1
           ld2(4)=a2
           ld3(4)=a3
           ld4(4)=a4
           sfilter='I'
           call ld(nlines,string,temp,logg,a1,a2,a3,a4,sfilter)      
           ld1(5)=a1
           ld2(5)=a2
           ld3(5)=a3
           ld4(5)=a4
           goto 5
        end if
      end do

C     
C     Limb darkening and shift in radial velocity corresponding
C     to the specific tile. 

 5    continue
      
      gamma=(1.0d0+vp/c)
      
      do i=1,npts
         if(w(i).le.4400.0) then
            ib=1
            wmin=3600.0
            wden=800.0
         else if(w(i).gt.4400.and.w(i).le.5500) then
            ib=2
            wmin=4400.0
            wden=1100.0
         else if(w(i).gt.5500.and.w(i).le.6900) then
            ib=3
            wmin=5500.0
            wden=1400.0
         else if(w(i).gt.6900.and.w(i).le.9500) then
            ib=4
            wmin=6900.0
            wden=2600.0
         end if
         
         wratio=(w(i)-wmin)/wden
         ld1w=ld1(ib)+wratio*(ld1(ib+1)-ld1(ib))
         ld2w=ld2(ib)+wratio*(ld2(ib+1)-ld2(ib))
         ld3w=ld3(ib)+wratio*(ld3(ib+1)-ld3(ib))
         ld4w=ld4(ib)+wratio*(ld4(ib+1)-ld4(ib))
         limbdark=1.0-(ld1w*(1.0-mu**0.5)+ld2w*(1.0-mu)
     +        +ld3w*(1.0-mu**1.5)+ld4w*(1.0-mu**2))
 
         wvp(i)=gamma*w(i)
         fareapn(i)=areapn*limbdark*f(i)
      end do

C     HI-RES spectrum: the wavelength mesh of the individual spectra
C     for each tile, is put into a common mesh to add all the spectra

      fa1=(fareapn(2)-fareapn(1))/(wvp(2)-wvp(1))
      fan=(fareapn(npts)-fareapn(npts-1))/(wvp(npts)-wvp(npts-1))
      
      call spline(wvp,fareapn,npts,fa1,fan,f2n)
      
      do i=1,npts
        wint=w(i)
        call splint(wvp,fareapn,f2n,npts,wint,fint)
        fmapped(i)=fint
      end do
      
C     HI-RES spectrum: the individual spectra of each tile
C     are added up to create the final spectrum

      do i=1,npts
         ftot(i)=ftot(i)+fmapped(i)
      end do
      
C     -----
      
      ithold=ith
      goto 3
      
 300  close(unit=3)

C     STORING HI-RES spectrum
      
      open(unit=13,file=outspec,form='formatted',status='unknown')      
      do i=1,npts
         write(13,130) w(i),ftot(i)
      end do
      close(unit=13)
 130  format(1pe14.8,2x,1pe13.6)
      
      write(*,131) outspec
 131  format(/,'Final HI-RES  spectrum is ',a60)

      end

      subroutine mod4lat(kcode,m,c1,c2,nptm4,w4,f4)

      implicit none

      integer i,nptm4

      real*4 wmin,wmax
      real*4 c1,c2
      real*4 wave,f(4)
      real*4 w4(45000),f4(45000)

      character*3  kcode
      character*10 m(4)
      character*16 suff
      character*13 path

      character*86 modext(4)

      path='../00_models/'
      suff=kcode//'_360_550_v0.m'
      
C
C     |  m(i)  |
C     ----------
C     t09900g450p00_360_550_v0.m
C
         wmin=3600.0
         wmax=5500.0

      do i=1,4
        modext(i)=path//m(i)//suff
      end do

      i=0
      open(unit=11,file=modext(1),form='formatted',status='old')
      open(unit=12,file=modext(2),form='formatted',status='old')
      open(unit=13,file=modext(3),form='formatted',status='old')
      open(unit=14,file=modext(4),form='formatted',status='old')
 1    read(11,*,end=100) wave,f(1)
      read(12,*,end=100) wave,f(2)
      read(13,*,end=100) wave,f(3)
      read(14,*,end=100) wave,f(4)
      if(wave.lt.wmin.or.wave.gt.wmax) goto 1
      i=i+1
      w4(i)=wave
      f4(i)=c2*(c1*f(1)+(1.0-c1)*f(2))+(1.0-c2)*(c1*f(3)+(1.0-c1)*f(4))     
      goto 1

 100  close(unit=11)
      close(unit=12)
      close(unit=13)
      close(unit=14)
      nptm4=i-1

      return
      end

C     LIMB DARKENING ROUTINES 
      
      subroutine subsetld(kcode,mint,maxt,ming,maxg,nlines,s)

C     This routine creates a subset of data from the files
C     00_ld/ld.coeffs. .txt containing only the relevant lines
C     to cover the range of temperatures and gravities of the
C     specific star that is modelled.
      
      implicit none

      integer nlines

      real*4 mint,maxt,ming,maxg
      real*4 t,g
      
      character*3  kcode
      character*4  sab
      character*4  sg
      character*5  st
      character*26 ifile
      character*65 string,s(1000)
      
C        10        20        30        40        50        60
C12345678901234567890123456789012345678901234567890123456789012345      
C 5.00 40000.  1.0  2.0  0.9591 -1.5398  1.3559 -0.4532 J  L ATLAS
C 3.50 30000. -5.0  2.0  1.1409 -1.2284  0.9583 -0.3052 Kp L A
C
C logg   T    [M/H} vt  a1      a2      a3      a4      F  Met Mod
C
      ifile='../00_ld/ld.coeffs.'//kcode//'.txt'
      
      nlines=0
      mint=mint-1100.0
      maxt=maxt+1100.0
      ming=ming-0.60
      maxg=maxg+0.60
      open(unit=1,file=ifile,form='formatted',status='old')
 1    read(1,15,end=100) string
 15   format(a65)
      st=string(7:11)
      sg=string(2:5)
      read(st,*) t
      read(sg,*) g

      if (t.lt.mint.or.t.gt.maxt.or.
     +    g.lt.ming.or.g.gt.maxg) goto 1
      if (string(20:22).ne.'2.0')   goto 1
      if (string(56:57).eq.'C ')    goto 1
      if (string(56:56).eq.'J')     goto 1
      if (string(56:56).eq.'H')     goto 1
      if (string(56:56).eq.'K')     goto 1
      if (string(56:56).eq.'u')     goto 1
      if (string(56:56).eq.'b')     goto 1
      if (string(56:56).eq.'v')     goto 1
      if (string(56:56).eq.'y')     goto 1
      if (string(56:56).eq.'g')     goto 1
      if (string(56:56).eq.'r')     goto 1
      if (string(56:56).eq.'i')     goto 1
      if (string(56:56).eq.'z')     goto 1
      if (string(56:57).eq.'Kp')    goto 1
      if (string(56:56).eq.'S')     goto 1
      if (string(61:65).ne.'ATLAS') goto 1
      nlines=nlines+1
      s(nlines)=string
      goto 1
 100  close(unit=1)

      return
      end

      subroutine ld(nlines,string,tint,gint,a1,a2,a3,a4,sf)

      implicit none

      integer i,j
      integer nlines
      integer i11,i12,i21,i22

      real*4 tint,gint
      real*4 t(nlines),g(nlines)
      real*4 t1,t2,g1,g2
      real*4 a1,a2,a3,a4

      character*1  sf
      character*18 ifile
      character*65 string(nlines)

      character*6  st1,st2
      character*4  sg1,sg2
      character*11 cstring11,cstring12,cstring21,cstring22  

      do i=1,nlines
        read(string(i)(2:5),*)   g(i)
        read(string(i)(7:12),*)  t(i)
      end do

C     Bracketing the temperature: the subset of coefficients contains
C     information for 4 filters: UBVR, with the same number of
C     lines log g Teff for each filter, so to bracket the temperature
C     the first nlines/8 are explored.
      
      do j=1,nlines/5-1
         if(tint.ge.t(j).and.tint.lt.t(j+1)) goto 1001
      end do
      
 1001 t1=t(j)
      t2=t(j+1)

C     Bracketing the gravity

      do j=1,nlines/5-1
         if(gint.ge.g(j).and.gint.lt.g(j+1)) goto 1002
      end do

 1002 g1=g(j)
      g2=g(j+1)
      
C     Now the lines corresponding to:
C
C     g1  t1  ...
C     g1  t2  ...
C     g2  t1  ...
C     g2  t2  ...
C
C     are identified:

      write(st1,1003) t1
      write(st2,1003) t2
 1003 format(f6.0)
      write(sg1,1004) g1
      write(sg2,1004) g2
 1004 format(f4.2)

      cstring11=' '//sg1//' '//st1
      cstring12=' '//sg1//' '//st2
      cstring21=' '//sg2//' '//st1
      cstring22=' '//sg2//' '//st2

      i11=0
      i12=0
      i21=0
      i22=0
      do j=1,nlines/5
         if(cstring11.eq.string(j)(1:11)) i11=j
         if(cstring12.eq.string(j)(1:11)) i12=j
         if(cstring21.eq.string(j)(1:11)) i21=j
         if(cstring22.eq.string(j)(1:11)) i22=j
      end do

C     Since only the lines 1 to nlines/5 -which correspond to filter U-
C     have been explored, the relevant lines corresponding to other
C     filters are easily found as i11+ , i12+ ,i21+ , i22+(nlines/5)*(n-1)
C     where n=2,...,5 correspond to BVRI

C     The following routine ldinterp computes the coefficients
     
      call ldinterp(tint,gint,t1,t2,g1,g2,nlines,
     +     string,sf,i11,i12,i21,i22,a1,a2,a3,a4)

      return
      end

      subroutine ldinterp(tint,gint,t1,t2,g1,g2,nlines,
     +           string,sf,i11,i12,i21,i22,a1,a2,a3,a4)

      implicit none

      integer i,j
      integer nlines
      integer i11,i12,i21,i22
      integer nf

      real*4 tint,gint
      real*4 t1,t2,g1,g2
      real*4 a(4,4)

      real*4 dg2g1,dgg1,dt2t1,dtt1
      real*4 a1g1,a2g1,a3g1,a4g1,a1g2,a2g2,a3g2,a4g2
      real*4 a1,a2,a3,a4
      
      character*1  sf
      character*65 string(nlines)

      if(sf.eq.'U') nf=1
      if(sf.eq.'B') nf=2
      if(sf.eq.'V') nf=3
      if(sf.eq.'R') nf=4
      if(sf.eq.'I') nf=5

      i11=i11+(nlines/5)*(nf-1)
      i12=i12+(nlines/5)*(nf-1)
      i21=i21+(nlines/5)*(nf-1)
      i22=i22+(nlines/5)*(nf-1)

      read(string(i11)(24:30),*) a(1,1)
      read(string(i11)(32:38),*) a(2,1)
      read(string(i11)(40:46),*) a(3,1)
      read(string(i11)(48:54),*) a(4,1)

      read(string(i12)(24:30),*) a(1,2)
      read(string(i12)(32:38),*) a(2,2)
      read(string(i12)(40:46),*) a(3,2)
      read(string(i12)(48:54),*) a(4,2)
     
      read(string(i21)(24:30),*) a(1,3)
      read(string(i21)(32:38),*) a(2,3)
      read(string(i21)(40:46),*) a(3,3)
      read(string(i21)(48:54),*) a(4,3)

      read(string(i22)(24:30),*) a(1,4)
      read(string(i22)(32:38),*) a(2,4)
      read(string(i22)(40:46),*) a(3,4)
      read(string(i22)(48:54),*) a(4,4)

C     Notation: a1g1 = coeff a1 for the lower gravity
C               a3g2 = coeff a3 for the higher gravity

      dt2t1=t2-t1
      dtt1 =tint-t1
       
      a1g1=a(1,1)+(a(1,2)-a(1,1))/dt2t1*dtt1
      a2g1=a(2,1)+(a(2,2)-a(2,1))/dt2t1*dtt1
      a3g1=a(3,1)+(a(3,2)-a(3,1))/dt2t1*dtt1
      a4g1=a(4,1)+(a(4,2)-a(4,1))/dt2t1*dtt1
                                     
      a1g2=a(1,3)+(a(1,4)-a(1,3))/dt2t1*dtt1
      a2g2=a(2,3)+(a(2,4)-a(2,3))/dt2t1*dtt1
      a3g2=a(3,3)+(a(3,4)-a(3,3))/dt2t1*dtt1
      a4g2=a(4,3)+(a(4,4)-a(4,3))/dt2t1*dtt1      

C     Final values

      dg2g1=g2-g1
      dgg1 =gint-g1
      
      a1=a1g1+(a1g2-a1g1)/dg2g1*dgg1
      a2=a2g1+(a2g2-a2g1)/dg2g1*dgg1
      a3=a3g1+(a3g2-a3g1)/dg2g1*dgg1
      a4=a4g1+(a4g2-a4g1)/dg2g1*dgg1
      
      return
      end


C     INTERPOLATION ROUTINES
C
C     spline
C     splint
C
C     from Numerical Recipes in FORTRAN

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
 11   continue
      
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
 12   continue
      
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
