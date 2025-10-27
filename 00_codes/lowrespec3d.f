      program lowrespec3d

      implicit none
      
      integer i,j,k,ib
      integer nlines
      integer nreg
      integer ilat,ilats,nlat
      integer nkur
      integer ivsini
      
      real*4 pi,c,rsuncm,gamma
      real*4 xr,yr,ith,ithold,areapn,vp
      real*4 c1,c2
      real*4 mu,xrmax,yrmax
      real*4 limbdark,a1,a2,a3,a4
      real*4 ld1(8),ld2(8),ld3(8),ld4(8)
      real*4 ld1w,ld2w,ld3w,ld4w
      real*4 lat(200),t(200),lg(200),aringsi(200)
      real*4 mint,maxt,ming,maxg
      real*4 kw4(1221),kf4(1221)
      real*4 integral,totalflux(200)
      real*4 kflat(200,1221)
      real*4 kt,kg
      real*4 temp,logg
      real*4 lum,lsunsi
      real*4 wmin,wden,wratio

      real*4 kw(1221),kf(1221)
      real*4 kwvp(1221),kfareapn(1221)
      real*4 kfa1,kfan,kf2n(1221),kwint,kfint
      real*4 kfmapped(1221),kftot(1221)      
      real*4 dum1
      
      character*1  sfilter
      character*3  kcode
      character*4  metal
      character*10 m(4)
      character*45 stin,stdum,sttpole,stlgpole,strpole,stom    
      character*60 outspec,koutspec
      character*65 s(1000),string(1000)
      
      parameter(c=299792.458)
      parameter(nkur=1221)
      parameter(lsunsi=3.827e+26)
      pi=acos(-1.0)
C     
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
      read(1,11,end=1100) lat(ilat),t(ilat),lg(ilat),
     +                    (m(j),j=1,4),c1,c2,aringsi(ilat)
 11   format(1x,f6.2,1x,f7.1,1x,f5.3,4(1x,a10),1x,
     +       2(1x,f6.4),2x,1pe16.8)
      
      write(*,12) abs(lat(ilat)),int(t(ilat)),lg(ilat)
 12   format('Building Kurucz LOW-RES models for th=+-',
     +        f5.2,' T=',i5,' log g=',f4.2)
      
      ilats=nlat-(ilat-1)
      lat(ilats)=-lat(ilat)
      t(ilats)=t(ilat)
      lg(ilats)=lg(ilat)
      mint=min(mint,t(ilat))
      maxt=max(maxt,t(ilat))
      ming=min(ming,lg(ilat))
      maxg=max(maxg,lg(ilat))
      
C     LOW-RESOLUTION models

         kt=t(ilats)
         kg=lg(ilats)
         call kur4lat(kcode,kt,kg,kw4,kf4)
         do k=1,nkur
            kw(k)=kw4(k)
            kflat(ilat,k)=kf4(k)
            kflat(ilats,k)=kf4(k)
         end do

C     Computing the integral of the models. The fluxes in the
C     Kurucz models are given in erg cm-2 s-1 A-1, therefore
C     the integral would be in units erg cm-2 s-1. Since the
C     parameter aringsi has been computed in SI units, the integral
C     is divided by 1.0e+3 to convert it into watts m-2

         call integrate_trap(kw4,kf4,nkur,integral)
         totalflux(ilat)=integral/1000.0
         
      end do
      
 1100 close(unit=1)

C     Multiplying the integral the model at each latitude 
C     in the grid by the area of the corona (aringsi) f
C     for one hemisphere

      lum=0.0
      do i=1,nlat/2
         lum=lum+aringsi(i)*totalflux(i)
      end do
 
C     Since the star is simmetric with respect to the
C     equator, we multiply by 2.0 the result for one
C     hemisphere.

      lum=2.0*lum/lsunsi
      
      write(*,99)

C     Writing a subset of the data for the limb-darkening coefficients
C     using only the specific max/min temperatures and gravities
C     of the star.

      call subsetld(kcode,mint,maxt,ming,maxg,nlines,s)
      do i=1,nlines
         string(i)=s(i)
      end do
      
C     ----- LOW-RES spectrum: initializing total flux
      do i=1,nkur
         kftot(i)=0.0
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
            do i=1,nkur
               kf(i)=kflat(ilat,i)
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
           sfilter='J'
           call ld(nlines,string,temp,logg,a1,a2,a3,a4,sfilter)      
           ld1(6)=a1
           ld2(6)=a2
           ld3(6)=a3
           ld4(6)=a4
           sfilter='H'
           call ld(nlines,string,temp,logg,a1,a2,a3,a4,sfilter)      
           ld1(7)=a1
           ld2(7)=a2
           ld3(7)=a3
           ld4(7)=a4
           sfilter='K'
           call ld(nlines,string,temp,logg,a1,a2,a3,a4,sfilter)      
           ld1(8)=a1
           ld2(8)=a2
           ld3(8)=a3
           ld4(8)=a4
           goto 5
        end if
      end do

C     
C     Limb darkening corresponding to the specific tile. 

 5    continue

C     ----- LOW-RES spectrum: Wavelengths and fluxes

      do i=1,nkur
         if(kw(i).le.1200.or.kw(i).gt.22000) then
            limbdark=1.0
            goto 6
         end if
         if(kw(i).gt.1200.and.kw(i).le.3600) then
            ib=1
            wmin=1200.0
            wden=2400.0
         else if(kw(i).gt.3600.and.kw(i).le.4400.0) then
            ib=1
            wmin=3600.0
            wden=800.0
         else if(kw(i).gt.4400.and.kw(i).le.5500) then
            ib=2
            wmin=4400.0
            wden=1100.0
         else if(kw(i).gt.5500.and.kw(i).le.6820) then
            ib=3
            wmin=5500.0
            wden=1320.0
         else if(kw(i).gt.6820.and.kw(i).le.8660) then
            ib=4
            wmin=6820.0
            wden=1840.0
         else if(kw(i).gt.8660.and.kw(i).le.12320) then
            ib=5
            wmin=8660.0
            wden=3660.0
         else if(kw(i).gt.12320.and.kw(i).le.16400) then
            ib=6
            wmin=12320.0
            wden=4080.0
         else if(kw(i).gt.16400.and.kw(i).le.22000) then
            ib=7
            wmin=16400.0
            wden=5600.0
         end if

         wratio=(kw(i)-wmin)/wden
         ld1w=ld1(ib)+wratio*(ld1(ib+1)-ld1(ib))
         ld2w=ld2(ib)+wratio*(ld2(ib+1)-ld2(ib))
         ld3w=ld3(ib)+wratio*(ld3(ib+1)-ld3(ib))
         ld4w=ld4(ib)+wratio*(ld4(ib+1)-ld4(ib))
         limbdark=1.0-(ld1w*(1.0-mu**0.5)+ld2w*(1.0-mu)
     +               +ld3w*(1.0-mu**1.5)+ld4w*(1.0-mu**2))
 6       continue
         kfareapn(i)=areapn*limbdark*kf(i)
      end do
      
C     LOW-RES spectrum: the individual spectra of each tile
C     are added up to create the final spectrum

      do i=1,nkur
         kftot(i)=kftot(i)+kfareapn(i)
      end do
      
C     -----
      
      ithold=ith
      goto 3
      
 300  close(unit=3)
      
C     ----- STORING LOW-RES spectrum
      
      koutspec='low-res'//outspec
      open(unit=13,file=koutspec,form='formatted',status='unknown')      
      do i=1,nkur
         write(13,130) kw(i),kftot(i)
      end do
      close(unit=13)
 130  format(1pe14.8,2x,1pe13.6)
      
      write(*,131) koutspec,lum
 131  format(/,86('-'),/,
     +         ' Final LOW-RES spectrum is ',a60,/,
     +         ' Stellar luminosity L_star=',f8.3,' L_sun',/,
     +         86('-'),/) 
      end

C     ----------------------------------------------------------

C     KURUCZ LOW-RESOLUTION MODELS
      
      subroutine kur4lat(kcode,tint,gint,kw4,kf4)

      implicit none

      integer i,npt

      real*4 tint,gint    
      real*4 tmod(76)
      real*4 oldt,newt,t1,t2
      real*4 oldg,newg,g1,g2
      real*4 mh
      real*4 ct,cg

      real*4 wave(1221),flux(1221)
      real*4 w(1221),f(4,1221)
      real*4 ft1t2g1(1221),ft1t2g2(1221)
      real*4 kw4(1221),kf4(1221)

      character*3 kcode
      character*7 st1,st2
      character*4 sg1,sg2

      data tmod/ 3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250,
     +     5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250, 7500,
     +     7750, 8000, 8250, 8500, 8750, 9000, 9250, 9500, 9750,
     +    10000,10250,10500,10750,11000,11250,11500,11750,12000,
     +    12250,12500,12750,13000,14000,15000,16000,17000,18000,
     +    19000,20000,21000,22000,23000,24000,25000,26000,27000,
     +    28000,29000,30000,31000,32000,33000,34000,35000,36000,
     +    37000,38000,39000,40000,41000,42000,43000,44000,45000,
     +    46000,47000,48000,49000,50000/
      
C     Clipping the temperature
      
      oldt=3500.0
      do i=1,76
      newt=tmod(i)
      if(tint.ge.oldt.and.tint.lt.newt) then
         t1=oldt
         t2=newt
         goto 1000
      end if
      oldt=newt
      end do

 1000 continue

C     Clipping the gravity
      
      do i=1,10
         oldg=0.0+0.5*float(i)
         newg=0.0+0.5*float(i+1)
      if(gint.ge.oldg.and.gint.lt.newg) then
         g1=oldg
         g2=newg
         goto 1001
      end if
      end do

 1001 continue

C     Building some strings to look into the big files containing
C     containing all the models, and for naming the output model

C     Temperature
      
      write(st1,12) t1
      write(st2,12) t2
 12   format(f7.1)
      
C     log g * 100

      write(sg1,13) g1
      write(sg2,13) g2
 13   format(f4.2)
      
C     Interpolation constants

      ct=(tint-t2)/(t1-t2)
      cg=(gint-g2)/(g1-g2)

C     Extracting the models
      
      call a9read(st1,sg1,kcode,wave,flux)
      
      do i=1,1221
         w(i)=wave(i)
         f(1,i)=flux(i)
      end do

      call a9read(st1,sg2,kcode,wave,flux)
      do i=1,1221
         f(2,i)=flux(i)
      end do     

      call a9read(st2,sg1,kcode,wave,flux)
      do i=1,1221
         f(3,i)=flux(i)
      end do

      call a9read(st2,sg2,kcode,wave,flux)
      do i=1,1221
         f(4,i)=flux(i)
      end do
      
C     Building the interpolated model
      
C     1: interpolation in T

      do i=1,1221
         ft1t2g1(i)=ct*f(1,i)+(1.0-ct)*f(3,i)
         ft1t2g2(i)=ct*f(2,i)+(1.0-ct)*f(4,i)
      end do

C     2: interpolation in log g and storing the results

      do i=1,1221
         kw4(i)=w(i)
         kf4(i)=cg*ft1t2g1(i)+(1.0-cg)*ft1t2g2(i)
      end do
      
      return
      end

      subroutine integrate_trap(x, y, n, integral)
c------------------------------------------------------------
c  integrates y(x) using the trapezoidal rule for irregular x
c
c  input:
c     x(n)   - array of abscissae (not equally spaced)
c     y(n)   - array of ordinates
c     n      - number of data points
c
c  output:
c     integral - approximate integral of y dx from x(1) to x(n)
c
c  note: x must be in ascending order
c------------------------------------------------------------
      integer n, i
      real*4 x(n), y(n), integral, dx

      integral = 0.0
      do 10 i = 1, n-1
         dx = x(i+1) - x(i)
         integral = integral + 0.5 * dx * (y(i+1) + y(i))
10    continue

      return
      end

      subroutine a9read(st,sg,smh,wave,flambda)

      implicit none

      integer i

      real*4 c,pi
      real*4 wave,hnu,hnucont,flambda
      
      character*3  smh
      character*4  sg
      character*7  st
      character*31 ifile
      character*80 line
      
      dimension wave(1221),flambda(1221)
      dimension hnu(1221),hnucont(1221)

C     The speed of light in cm s-1

      parameter(c=29979245800.0,pi=3.1415926534)

      ifile='../00_pckfiles/f'//smh//'k2odfnew.pck'
      open(unit=10,file=ifile,form='formatted')

C     Here we read the wavelength in nm

      read(10,1) wave
    1 format(8f10.2)
      
C     Here we look for the header identifying the effective temperature
C     and log gravity for which we want to extract the fluxes

C  TEFF   3500.  GRAVITY 0.00000   [0.0] VTURB=2  L/H=1.25 NOVER NEW ODF
C  1234567890123456789012345678901234567890
      
 10   read(10,100,end=1000) line
      if(line(7:11).eq.st(1:5).and.line(23:26).eq.sg) goto 1010
      goto 10
      
 100  format(a80)

C     Here we read hnu      (erg cm-2 s-1 Hz-1 ster-1)
C                  hnucont  (erg cm-2 s-1 Hz-1 ster-1)

 1010 read(10,110) hnu
      read(10,110) hnucont
 110  format(8e10.4)

 1000 close(unit=10)

      do i=1,1221
         flambda(i)=4.0*pi*hnu(i)*(c*1.0e8)/wave(i)/wave(i)/100.0
         wave(i)=wave(i)*10
      end do

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

C     Bracketing the temperature: the file temp.ld.coeffs.txt contains
C     information for 8 filters: UBVRIJHK, with the same number of
C     lines log g Teff for each filter, so to bracket the temperature
C     the first nlines/8 are explored.
      
      do j=1,nlines/8-1
         if(tint.ge.t(j).and.tint.lt.t(j+1)) goto 1001
      end do
      
 1001 t1=t(j)
      t2=t(j+1)

C     Bracketing the gravity

      do j=1,nlines/8-1
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
      do j=1,nlines/8
         if(cstring11.eq.string(j)(1:11)) i11=j
         if(cstring12.eq.string(j)(1:11)) i12=j
         if(cstring21.eq.string(j)(1:11)) i21=j
         if(cstring22.eq.string(j)(1:11)) i22=j          
      end do

C     Since only the lines 1 to nlines/8 -which correspond to filter U-
C     have been explored, the relevant lines corresponding to other
C     filters are easily found as i11+ , i12+ ,i21+ , i22+(nlines/8)*(n-1)
C     where n=2,...,8 correspond to BVRIJHK

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
      if(sf.eq.'J') nf=6
      if(sf.eq.'H') nf=7
      if(sf.eq.'K') nf=8  

      i11=i11+(nlines/8)*(nf-1)
      i12=i12+(nlines/8)*(nf-1)
      i21=i21+(nlines/8)*(nf-1)
      i22=i22+(nlines/8)*(nf-1)

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
