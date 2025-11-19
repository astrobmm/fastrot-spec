      program spec3d

C     This program computes the synthetic spectrum of an oblate
C     star with the input parameters specified in file input.star3d
C     and a grid on the stellar surface with tiles of sizes given
C     in file input.grid3d

      implicit none
      
      integer i,j,k,ib
      integer itile,ncells,ncstep,ncmod
      integer ilat,ilats,nlat,npts
      integer ivsini
      
      real*4 pi,c,gamma
      real*4 xr,yr,ith,ithold,areapn,vp
      real*4 c1,c2
      real*4 mu,xrmax,yrmax
      real*4 limbdark,a1,a2,a3,a4
      real*4 ld1(5),ld2(5),ld3(5),ld4(5)
      real*4 ld1w,ld2w,ld3w,ld4w
      real*4 lat(200),t(200),lg(200),wlat(200,45000),flat(200,45000)
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
      character*60 outspec
      
      parameter(c=299792.458)
      
      pi=acos(-1.0)
C     
      open(unit=3,file='input.spec3d',form='formatted',status='old')
      read(3,101) metal
 101  format(a4)
      read(3,102) outspec
 102  format(a60)
      read(3,103) nlat
 103  format(i3)
      read(3,104) ncells
 104  format(i6)
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

C     spec3d.lat:  Here we read the four models 
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
            
      open(unit=1,file='spec3d.lat',form='formatted',status='old')
      do ilat=1,nlat/2
      read(1,10)  lat(ilat),t(ilat),lg(ilat),(m(j),j=1,4),c1,c2
      write(*,11) abs(lat(ilat)),int(t(ilat)),lg(ilat)      
 10   format(1x,f6.2,1x,f7.1,1x,f5.3,4(1x,a10),1x,2(1x,f6.4))
 11   format('Building models for th=+-',f5.2,' T=',i5,' log g=',f4.2)
      
      ilats=nlat-(ilat-1)
      lat(ilats)=-lat(ilat)
      t(ilats)=t(ilat)
      lg(ilats)=lg(ilat)
      
      call mod4lat(kcode,m,c1,c2,npts,w4,f4)
      
        do j=1,npts
           wlat(ilat,j)=w4(j)
           flat(ilat,j)=f4(j)
           wlat(ilats,j)=w4(j)         
           flat(ilats,j)=f4(j)      
        end do     
      end do
      
      close(unit=1)

      write(*,99)
      
C     tiles.a:  Here we read the file containing information on the
C     visible tiles. Relevant parameters are the latitude of each 
C     surface element, the projected area and the velocity shift to be
C     applied to each individual spectrum.
C
C     [dum1 is a dummy variable, not used]
      
      do i=1,npts
         ftot(i)=0.0
      end do

      ithold=-999.0
      itile=0
      write(*,25) ncells
 25   format(/,'Number of visible cells: ',i6,/)
      ncstep=ncells/20
      
      open(unit=3,file='tiles.a',form='formatted',status='old')
 3    read(3,30,end=300) xr,yr,ith,dum1,areapn,vp,mu
 30   format(7(1pe14.6,2x))
      itile=itile+1
      ncmod=mod(itile,ncstep)
      if(ncmod.eq.0) write(*,31) itile,ncells
 31   format('Spectra computed for ',i6,' out of',i6,' cells') 
C     Here we assign the fluxes w(ilat),f(ilat,*) to the tile      
C     according to its latitude

      if(ith.eq.ithold) goto 5

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
           call ld(temp,logg,a1,a2,a3,a4,kcode,sfilter)
           ld1(1)=a1
           ld2(1)=a2
           ld3(1)=a3
           ld4(1)=a4
           sfilter='B'
           call ld(temp,logg,a1,a2,a3,a4,kcode,sfilter)
           ld1(2)=a1
           ld2(2)=a2
           ld3(2)=a3
           ld4(2)=a4      
           sfilter='V'
           call ld(temp,logg,a1,a2,a3,a4,kcode,sfilter)      
           ld1(3)=a1
           ld2(3)=a2
           ld3(3)=a3
           ld4(3)=a4
           sfilter='R'
           call ld(temp,logg,a1,a2,a3,a4,kcode,sfilter)
           ld1(4)=a1
           ld2(4)=a2
           ld3(4)=a3
           ld4(4)=a4
           sfilter='I'
           call ld(temp,logg,a1,a2,a3,a4,kcode,sfilter)
           ld1(5)=a1
           ld2(5)=a2
           ld3(5)=a3
           ld4(5)=a4           
           goto 5
         end if
      end do
C     
C     Limb darkening and shift in radial velocity corresponding
C     to the specific tile

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
     +               +ld3w*(1.0-mu**1.5)+ld4w*(1.0-mu**2))
         wvp(i)=gamma*w(i)
         fareapn(i)=areapn*limbdark*f(i)
      end do

C     Here we put the wavelength mesh of the individual spectra
C     for each tile, into a common mesh to add all the spectra
      
      fa1=(fareapn(2)-fareapn(1))/(wvp(2)-wvp(1))
      fan=(fareapn(npts)-fareapn(npts-1))/(wvp(npts)-wvp(npts-1))
      
      call spline(wvp,fareapn,npts,fa1,fan,f2n)

      do i=1,npts
        wint=w(i)
        call splint(wvp,fareapn,f2n,npts,wint,fint)
        fmapped(i)=fint
      end do

C     Here we start accummulating the individual spectrum of each
C     tile to create the final spectrum

      do i=1,npts
         ftot(i)=ftot(i)+fmapped(i)
      end do

      ithold=ith
      goto 3
 300  close(unit=3)

      open(unit=13,file=outspec,form='formatted',status='unknown')      
      do i=1,npts
         write(13,130) w(i),ftot(i)
      end do
      close(unit=13)
 130  format(f11.5,2x,1pe13.6)
      
      write(*,131) outspec
 131  format(/,'Final spectrum is ',a60,/)

      stop
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

      subroutine ld(tint,gint,a1,a2,a3,a4,kcode,sf)

      implicit none

      integer i,j
      integer ntgrid

      real*4 g(11)
      real*4 trest(61),tp00(76),tp05(59)
      real*4 t(80)
      real*4 gint,tint
      real*4 gold,told,gnew,tnew
      real*4 g1,g2,t1,t2
      real*4 dg2g1,dgg1,dt2t1,dtt1
      real*4 a(4,4)
      real*4 a1g1,a2g1,a3g1,a4g1,a1g2,a2g2,a3g2,a4g2
      real*4 a1,a2,a3,a4

      character*1  sf
      character*3  kcode
      character*4  sg1,sg2
      character*5  st1,st2
      character*11 string(4)
      character*26 ifile
      character*65 longstring,cstring(4)

      data g/0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0/

C     In the tables of coefficients for limb darkening the grid
C     of temperatures for m25,m20,m15,m10,m05,p02 has 61 elements,
C     they are under the common variable trest

      data trest/3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250,
     +           5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250,
     +           7500, 7750, 8000, 8250, 8500, 8750, 9000, 9250,
     +           9500, 9750,10000,10500,11000,11500,12000,12500,
     +          13000,14000,15000,16000,17000,18000,19000,20000,
     +          21000,22000,23000,24000,25000,26000,27000,28000,
     +          29000,30000,31000,32000,33000,34000,35000,37500,
     +          40000,42500,45000,47500,50000/
      
C    For p00 there are 76 temperatures
      
      data tp00/3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250,
     +          5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250,
     +          7500, 7750, 8000, 8250, 8500, 8750, 9000, 9250,
     +          9500, 9750,10000,10250,10500,10750,11000,11250,
     +         11500,11750,12000,12250,12500,12750,13000,14000,
     +         15000,16000,17000,18000,19000,20000,21000,22000,
     +         23000,24000,25000,26000,27000,28000,29000,30000,
     +         31000,32000,33000,34000,35000,36000,37000,38000,
     +         39000,40000,41000,42000,43000,44000,45000,46000,
     +         47000,48000,49000,50000/
      
C    For p05 there are 59 temperatures
      
      data tp05/3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250,
     +          5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250,
     +          7500, 7750, 8000, 8250, 8500, 8750, 9000, 9250,
     +          9500, 9750,10000,10500,11000,11500,12000,12500,
     +         13000,14000,15000,16000,17000,18000,19000,20000,
     +         21000,22000,23000,24000,25000,26000,27000,28000,
     +         29000,30000,31000,32000,33000,34000,35000,37500,
     +         40000,42500,45000/
      
      if(kcode.eq.'m25'.or.kcode.eq.'m20'.or.kcode.eq.'m15'.or.
     +     kcode.eq.'m05'.or.kcode.eq.'p02') then
           ntgrid=61
           do i=1,ntgrid
              t(i)=trest(i)
           end do
        else if(kcode.eq.'p00') then
           ntgrid=76
           do i=1,ntgrid
              t(i)=tp00(i)
           end do
        else if(kcode.eq.'p05') then
           ntgrid=59
           do i=1,ntgrid
              t(i)=tp05(i)
           end do
        end if

      gold=g(1)
      told=t(1)
  
      do i=1,11
        gnew=g(i)
        if(gint.ge.gold.and.gint.lt.gnew) then
           g1=gold
           g2=gnew
           goto 10
        else
           gold=gnew
        end if
      end do
      
 10   do i=1,ntgrid
         tnew=t(i)
         if(tint.ge.told.and.tint.lt.tnew) then
            t1=told
            t2=tnew
            goto 11
          else
            told=tnew
         end if
      end do

 11   write(sg1,100) g1
      write(sg2,100) g2
 100  format(f4.2)
      write(st1,101) int(t1)
      write(st2,101) int(t2)
 101  format(i5)

C     Building the four strings to be extracted from the
C     file ld.coeffs.txt
C             1         2         3         4         5         6 
C    12345678901234567890123456789012345678901234567890123456789012345
C     0.00  3500.  0.0  2.0  0.5317  0.8077 -1.0003  0.5773 U  L ATLAS
C
      string(1)=' '//sg1//' '//st1
      string(2)=' '//sg1//' '//st2
      string(3)=' '//sg2//' '//st1
      string(4)=' '//sg2//' '//st2
      
      ifile='../00_ld/ld.coeffs.'//kcode//'.txt'
      j=0
      open(unit=1,file=ifile,form='formatted',status='old')
 1    read(1,104,end=1000) longstring
 104  format(a65)
       if(longstring(56:56).ne.sf) goto 1
       if(longstring(1:11).eq.string(1).or.
     +    longstring(1:11).eq.string(2).or.
     +    longstring(1:11).eq.string(3).or.
     +    longstring(1:11).eq.string(4)) then
          j=j+1
          cstring(j)=longstring
          read(cstring(j)(24:30),105) a(1,j)
          read(cstring(j)(32:38),105) a(2,j)      
          read(cstring(j)(40:46),105) a(3,j)
          read(cstring(j)(48:54),105) a(4,j)
  105     format(f7.4) 
       end if
       goto 1
 1000  close(unit=1)

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
      
C      write(*,106) a1,a2,a3,a4
C 106  format(4(f7.4,2x))

      return
      end
     
      subroutine spline(x,y,n,yp1,ypn,y2)

      implicit real*4 (a-h,o-z)
      parameter(nmax=45000)
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

  
