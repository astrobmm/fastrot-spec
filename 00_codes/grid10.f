      program grid10

C     Creates a file with a grid of parallels and meridians
C     separated by 10 degrees. Points are calculated with steps
C     of 0.5 deg both in latitude and in longitude.

      implicit none

      real*8 pi,deg2rad,rad2deg

      real*8 om
      real*8 in,al,sal,cal
      real*8 ith,iph,th,ph,abscolatr,vth
      real*8 sth,cth,sph,cph
      real*8 sithpar,siphmer,sithmid,siphmid
      real*8 ithmod,iphmod
      real*8 x,y,z,xr,yr,zr
      real*8 r

      character*10 stin,stom
      character*10 sdummy
      
      real*8 radius,angle
      external radius,angle

      common/comega/om
      common/cradius/r
      
      pi=acos(-1.0d0)
      deg2rad=pi/180.0d0
      rad2deg=180.0d0/pi

      open(unit=1,file='input.star3d',form='formatted',status='old')
 91   read(1,101) stin
      if(stin(1:1).eq.'*') goto 91
      read(1,101) sdummy
      read(1,101) sdummy
      read(1,101) sdummy
      read(1,101) stom
      close(unit=1)
 101  format(a10)
C     
      read(stin,*) in
      read(stom,*) om
C     
      siphmer=0.5
      sithpar=0.5
      
      al=(90.0d0-in)*deg2rad
      sal=sin(al)
      cal=cos(al)

      open(unit=1,file='parallels.a',form='formatted',status='unknown')
      open(unit=11,file='par10.a',form='formatted',status='unknown')
      open(unit=12,file='par10i90.a',form='formatted',status='unknown')
      
      do 1 ith=-90.0d0,90.0d0,sithpar
         ithmod=mod(ith,10.0d0)
         r=radius(ith)
         th=deg2rad*ith
         sth=sin(th)
         cth=cos(th)
         do 2 iph=-180.0d0,180.0d0,siphmer
            ph=deg2rad*iph
            sph=sin(ph)
            cph=cos(ph)
            x=r*cth*sph
            y=r*sth
            z=r*cth*cph
            xr=x
            yr=y*cal-z*sal
            zr=y*sal+z*cal
C         
            if(zr.ge.0.0d0) write(1,1001) xr,yr
            if(zr.ge.0.0d0.and.ithmod.eq.0.0d0) write(11,1001) xr,yr
            if(z.ge.0.0d0.and.ithmod.eq.0.0d0) write(12,1001) x,y
 1001       format(2(1pe16.8,3x))
 2       continue
 1    continue

      close(unit=1)
      close(unit=11)
      close(unit=12)
      
C     Meridians
C     Step in degrees in longitude (meridians) siphmer

      open(unit=2,file='meridians.a',form='formatted',status='unknown')
      open(unit=22,file='mer10.a',form='formatted',status='unknown')
      open(unit=23,file='mer10i90.a',form='formatted',status='unknown')
      
      do 3 iph=-180.0d0,180.0d0,siphmer
         iphmod=mod(iph,10.0d0)
         ph=deg2rad*iph
         sph=sin(ph)
         cph=cos(ph)
         do 4 ith=-90.0d0,90.0d0,sithpar
            r=radius(ith)
            th=deg2rad*ith
            sth=sin(th)
            cth=cos(th)
            x=r*cth*sph
            y=r*sth
            z=r*cth*cph
            xr=x
            yr=y*cal-z*sal
            zr=y*sal+z*cal
            if(zr.ge.0.0d0) write(2,1001) xr,yr
            if(zr.ge.0.0d0.and.iphmod.eq.0.0d0) write(22,1001) xr,yr
            if(z.ge.0.0d0.and.ithmod.eq.0.0d0) write(23,1001) x,y
 4       continue
 3    continue

      close(unit=2)
      close(unit=22)
      close(unit=23)
      
      stop
      end
      
      function radius(latd)

      implicit none

      real*8 radius
      real*8 pi,deg2rad
      real*8 om
      real*8 latd,thd,thr
      real*8 olds,news
      real*8 r,stepr

      real*8 frad
      external frad

      common/comega/om
      common/ccolatitude/thr
      
      pi=dacos(-1.0d0)

      deg2rad=pi/180.0d0

      thd=90.0-abs(latd)
      thr=thd*deg2rad  
      r=0.01d0
      stepr=0.05d0

      olds=int(abs(frad(r))/frad(r))
      r=r+stepr
      news=int(abs(frad(r))/frad(r))
      
 1    if(news.eq.olds) then
         olds=news
         r=r+stepr
         news=int(abs(frad(r))/frad(r))
      else
         olds=news
         stepr=-stepr/2.0d0
         r=r+stepr
         news=int(abs(frad(r))/frad(r))
      end if

      if(abs(stepr).lt.1.0d-6) goto 100 

      goto 1

 100  radius=r
      
      return
      end

      function angle(latd)

      implicit none

      real*8 angle
      real*8 pi,deg2rad
      real*8 om
      real*8 r
      real*8 latd,thd,thr
      real*8 olds,news
      real*8 vth,stepvth

      real*8 fang
      external fang

      common/comega/om
      common/ccolatitude/thr
      
      pi=dacos(-1.0d0)

      deg2rad=pi/180.0d0

      thd=90.0-abs(latd)
      thr=thd*deg2rad  
      
C     Computation of vth

      if(thr.eq.0.0d0) then
         vth=0.0d0
         angle=vth
         return
      end if
      
      if(thr.eq.pi/2.0d0) then
         r=1.0d0
         angle=thr-pi/2.0d0/360.0d0
         return
      end if

      vth=pi/360.0d0
      stepvth=0.01

      olds=int(abs(fang(vth))/fang(vth))
      vth=vth+stepvth
      news=int(abs(fang(vth))/fang(vth))
      
 1    if(news.eq.olds) then
         olds=news
         vth=vth+stepvth
         news=int(abs(fang(vth))/fang(vth))
      else
         olds=news
         stepvth=-stepvth/2.0d0
         vth=vth+stepvth
         news=int(abs(fang(vth))/fang(vth))
      end if

      if(abs(stepvth).lt.1.0d-6) goto 100

      goto 1

 100  angle=vth
      
      return
      end
      
      function frad(r)
      
      implicit none
      
      real*8 frad,r
      real*8 om,thr
      
      common/comega/om
      common/ccolatitude/thr
      
      frad=1.0/om/om/r+0.5d0*(r*sin(thr))**2-(1.0d0/om/om+0.5d0)
      
      return
      end
      
      function fang(vth)

      implicit none
      real*8 fang,vth,r,pi
      real*8 om,thr
      
      common/comega/om
      common/cradius/r
      common/ccolatitude/thr

      pi=acos(-1.0)
      
      fang=cos(vth)+log(tan(vth/2.0d0))
     *     -(1.0d0/3.0d0)*om**2*(r*cos(thr))**3-cos(thr)
     *     -log(tan(thr/2.0d0))
      
      return
      end

      subroutine temp(om,r,ith,vth,t)

      implicit none

      real*8 oe,pi,deg2rad
      real*8 om,r,ith,colat,vth
      real*8 t,ratan
      
      common/cratan/ratan
      
      pi=acos(-1.0d0)
      deg2rad=pi/180.0d0
      
      colat=(90.0d0-abs(ith))*deg2rad
      
      oe=1.0d0/8.0d0

      t=(1.0d0/r**4+om**4*(r*sin(colat))**2
     *  -2.0d0*(om*sin(colat))**2/r)**oe
     *  *sqrt(tan(vth)/tan(colat))
      
      ratan=sqrt(tan(vth)/tan(colat))
      
      return
      end
      
      subroutine bounds(t,lg)

      integer tl,th,lgl,lgh

      real*8 t,lg
      real*8 c1,c2

      character*4 slgl,slgh,slglfinal,slghfinal
      character*6 stl,sth,stlfinal,sthfinal
      character*10 model1,model2,model3,model4

      common/cmodels/model1,model2,model3,model4,c1,c2
      
      tl=int(t/100.0d0)*100
      th=tl+100

      write(stl,1) tl
      write(sth,1) th
 1    format(i5)
      stlfinal='t0'//stl(2:5)
      sthfinal='t0'//sth(2:5)
      if(int(tl).ge.10000) stlfinal='t'//stl(1:5)
      if(int(th).ge.10000) sthfinal='t'//sth(1:5)
      
      if(lg.eq.3.5d0) then
         lgl=350
         lgh=350
      else if(lg.eq.4.0d0) then
         lgl=400
         lgh=400
      else if(lg.eq.4.5d0) then
         lgl=450
         lgh=450
      else if(lg.lt.3.5d0) then
         lgl=350
         lgh=400
      else if(lg.gt.3.5d0.and.lg.lt.4.0d0) then
         lgl=350
         lgh=400
      else if(lg.gt.4.0d0.and.lg.lt.4.5d0) then
         lgl=400
         lgh=450
      else if(lg.gt.4.5d0) then
         lgl=400
         lgh=450
      end if

      write(slgl,2) lgl
      write(slgh,2) lgh
 2    format(i3)
      slglfinal='g'//slgl
      slghfinal='g'//slgh
      
      model1=stlfinal//slglfinal
      model2=sthfinal//slglfinal
      model3=stlfinal//slghfinal
      model4=sthfinal//slghfinal
      
      c1=(float(th)-t)/(float(th)-float(tl))
      c2=(float(lgh)-lg*100.0)/(float(lgh)-float(lgl))

      return
      end

      subroutine last(string,ilast)

      implicit none

      integer i,ilast
      character*10 string 

      do i=1,10
         if(string(i:i).eq.' ') goto 10
      end do

 10   ilast=i-1

      return
      end
      

