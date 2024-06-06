      program star3d

      implicit none

      integer ilat,k,ktot,ivisible
      integer ii,im,it,ir,io,irsi
      integer ivsini
      
      real*8 pi,deg2rad,rad2deg
      
      real*8 om,omega,ocrit,v,vp,veq
      real*8 in,al,sal,cal
      real*8 ith,iph,th,ph,abscolatr,vth
      real*8 latc1,vthc1,tc1
      real*8 sth,cth,sph,cph
      real*8 sithpar,siphmer,sithmid,siphmid
      real*8 ithmod,iphmod
      real*8 x,y,z,xr,yr,zr
      real*8 ith0,ith1,ith2,r1,r2,ar1r2
      real*8 x1,y1,x2,y2,dr1r2
      real*8 xi,cthxi,sthxi,thplusxi
      real*8 area,areap,areapn
      real*8 sumarea,sumareapn,sumtareapn,sumlgareapn
      real*8 ld,mu
      real*8 lgmean,lg,lgareapn,tareapn

      real*8 gnewton,sigma,msun,lsun,rsun
      real*8 r,r_si,rsi2,rsimean,teffmean
      real*8 lum,sigmat4
      real*8 req,req_si
      real*8 tpole,teq,tpolenew
      real*8 t,tcons,t_k,t_k_vz,ratiot,ratior
      real*8 ming,maxg
      real*8 mstar,m_si
      real*8 lstar,l_si
      real*8 g,gr,gt,gpole
      real*8 ratan

      real*8 c1,c2
      
      character*1 tmodify
      character*3 stvsini,stfinalvsini
      character*4 metal
      character*10 m1,m2,m3,m4,model(2000)
      character*10 stin,stsiphmer,stsithpar
      character*10 stmstar,sttpole,streq
      character*10 stom,smetal
      character*60 outspec
      character*70 outpars

      real*8 radius,angle
      external radius,angle
      
      common/comega/om
      common/cradius/r
      common/temperature/t
      common/cmodels/m1,m2,m3,m4,c1,c2
      common/cratan/ratan

C     G, Msun and Rsun in SI units
     
      parameter(gnewton=6.67259d-11)
      parameter(sigma=5.670374d-8)
      parameter(msun=1.9899d+30)
      parameter(lsun=3.827d+26)
      parameter(rsun=6.9599d+8)
     
      pi=acos(-1.0d0)
      deg2rad=pi/180.0d0
      rad2deg=180.0d0/pi

C     Input parameters

C     File input.star3d
      
C     String      Numerical  Description of parameter            
C     ----------  ---------  ----------------------------------  
C     stin:       in         inclination  [deg, i=90 equator on] 
C     stmstar:    mstar      stellar mass       [solar units]    
C     sttpole:    tpole      polar temperature  [K]              
C     streq:      req        equatorial radius  [solar units]    
C     stom:       om         omega              [Omega/Omega_K]  
C     smetal                 [M/H]              [-2.5 ... +0.5]        
C      
      open(unit=1,file='input.star3d',form='formatted',status='old')
 91   read(1,101) stin
      if(stin(1:1).eq.'*') goto 91
      read(1,101) stmstar
      read(1,101) sttpole
      read(1,101) streq
      read(1,101) stom
      read(1,101) smetal  
      close(unit=1)

C     File input.grid3d

C     stsiphmer:  siphmer    deltaphi     [deg, tile separation] 
C     stsithpar:  sithpar    deltatheta   [deg, tile separation] 

      open(unit=1,file='input.grid3d',form='formatted',status='old')
 92   read(1,101) stsiphmer
      if(stsiphmer(1:1).eq.'*') goto 92
      read(1,101) stsithpar
      close(unit=1)
 101  format(a10)
      
C     
      read(stin,*)      in
      read(stmstar,*)   mstar
      read(sttpole,*)   tpole 
      read(streq,*)     req
      read(stom,*)      om
      read(stsiphmer,*) siphmer
      read(stsithpar,*) sithpar
      metal=smetal(1:4)

C     Dimensional quantities

      m_si=mstar*msun
      req_si=req*rsun

C     -----
C     Computing the constant that replaces the first parenthesis
C     in the second expression of eq. (31) of the paper
C     ER11: (L/4\pi\sigma R_e^2)^{1/4) which is basically T_eq
C
C     First we compute ratiot,the ratio T_eq/T_pole (eq. 32, ER11)
      
      ratiot=sqrt(2.0d0/(2.0d0+om*om))*(1.0d0-om*om)**(1.0d0/12.0d0)*
     *     exp(-4.0d0*om*om/3.0d0/(2.0d0+om*om)**3)

C     Since T_pole is an input parameter we compute T_eq:
      
      teq=tpole*ratiot
      
C     Now we take the first cell from the pole c1, compute the corresponding
C     radius, value of the adimensional t, and the angle vth (vartheta)
      
      latc1=-90.0d0+sithpar/2.0d0
      r=radius(latc1)
      vthc1=angle(latc1)
      call temp(om,r,latc1,vthc1,tc1)

C     The constant in front of eq. (31) substituting that parenthesis is
C     such as: T_eff(theta)=tcons*t(theta), therefore, for c1:
C     T_eff(c1)=tcons*t(c1)
 
      tcons=tpole/tc1

C     ...in that way, when we start the computation of t for all
C     colatitudes:
C
C     T_eff=tpole/t(c1)*[t: the rest of eq. 31, coming from routine temp)
C
C     For c1, T_eff=T_pole, for cn, the cell close to the equator
C     T_eff=tpole/t(c1)*t(eq)=tpole/
      
      ocrit=sqrt(gnewton*m_si/req_si**3)
      omega=om*ocrit
C
      al=(90.0d0-in)*deg2rad
      sal=sin(al)
      cal=cos(al)

C     Midpoints of cells

      sithmid=sithpar
      siphmid=siphmer
      sumarea=0.0d0
      sumareapn=0.0d0
      sumtareapn=0.0d0
      sumlgareapn=0.0d0
      lum=0.0d0
      ming=+99999.0d0
      maxg=-99999.0d0

      open(unit=1,file='tiles.a',form='formatted',status='unknown')
      open(unit=2,file='xyparams.a',form='formatted',status='unknown')
      open(unit=3,file='colat-r-t.a',form='formatted',status='unknown')
      open(unit=10,file='spec3d.lat',form='formatted',status='unknown')
C      open(unit=11,file='spec3d.vz',form='formatted',status='unknown')

      ilat=0
      irsi=0
      ivisible=0
      do 5 ith=-90.0d0+sithmid/2.0d0,90.0d0-sithmid/2.0d0,sithmid
         ilat=ilat+1
         
C     Computation of xi, the angle between the normal to the
C     differential surface element and the z direction before
C     the star is tilted by an angle i
         
         ith1=ith+sithpar/2.0d0
         r1=radius(ith1)
         x1=r1*cos(ith1*deg2rad)
         y1=r1*sin(ith1*deg2rad)
         
         ith2=ith-sithpar/2.0d0
         r2=radius(ith2)
         x2=r2*cos(ith2*deg2rad)
         y2=r2*sin(ith2*deg2rad)         

C     Length of the segment joining the points in the surface
C     with radii r1 and r2
         
         dr1r2=sqrt((y2-y1)**2+(x2-x1)**2)
         
         ar1r2=abs((r2-r1)/((r1+r2)/2.0d0))
         xi=atan(ar1r2/(sithpar*deg2rad))*rad2deg
C
         r=radius(ith)
         vth=angle(ith)
         
         if(ilat.eq.1) ratior=1.0d0/r

         if(ith.eq.0.0d0) then
            ith0=(pi/2.0d0-vth)*rad2deg
            call temp(om,r,ith0,vth,t)
         else
            call temp(om,r,ith,vth,t)
         end if
            
C     Dimensional quantities

         abscolatr=(90.0d0-abs(ith))*deg2rad 
         r_si=r*req_si
         v=omega*r_si*sin(abscolatr)/1000.0d0      
         gr=-gnewton*m_si/(r_si**2)+r_si*(omega*sin(abscolatr))**2
         gt=r_si*(omega**2)*sin(abscolatr)*cos(abscolatr)         
         g=sqrt(gr*gr+gt*gt)
         lg=log10(g*100.0d0)
         t_k=tcons*t
         rsi2=r_si**2
         rsimean=rsimean+r_si
         irsi=irsi+1
         sigmat4=sigma*t_k**4
         ming=min(ming,lg)
         maxg=max(maxg,lg)

C     Keeping the values of g at the pole to compute the temperature
C     structure using the Van Zeipel approximation (T ~ g**0.25)

         if(ilat.eq.1) then
            gpole=g
         end if

         t_k_vz=tpole*(g/gpole)**0.25
         
C
         th=ith*deg2rad
         thplusxi=(ith+xi)*deg2rad
         sth=sin(th)
         cth=cos(th)
         sthxi=sin(thplusxi)
         cthxi=cos(thplusxi)

C     The area of all cells at a given latitude (radius r) is computed
C     by multiplying the length of the segment joining the points
C     on the surface with radii r1, r2 (r1 < r < r2) by the length
C     of the cell measured on the corresponding parallel
         
         area=dr1r2*(r*siphmid*deg2rad*cth)      
         
         do 6 iph=-180.0d0+siphmid/2.0d0,180.0d0-siphmid/2.0d0,siphmid
           ph=deg2rad*iph
           sph=sin(ph) 
           cph=cos(ph)
           x=r*cth*sph
           y=r*sth
           z=r*cth*cph
C     
C     The variable  lum  after beeing accummulated at the end of the
C     do loops, will be the actual luminosity of the star.
C
           lum=lum+area*req_si*req_si*sigmat4
           xr=x
           yr=y*cal-z*sal
           zr=y*sal+z*cal
C           
C     The variable  area  is the area of each of the tiles on the
C     surface of the star. To find the area of the rotated element
C     projected onto the plane  x y  (areapn) we take the z component
C     of the rotated vector normal to the surface element.

C     areapn:  projected area of each tile on the sky (plane x y)
C     vp:      radial velocity of each tile to be applied to each
C              individual spectrum
C     sumarea: total area of the star (both zr>=0 and zr<0)          

           areap=area*(sth*sal+cth*cph*cal)
           areapn=area*(sthxi*sal+cthxi*cph*cal)
           sumarea=sumarea+area
           
C     ld:    angle between the normal to each tile and the line
C            of sight, needed to compute the limb darkening 

           mu=sthxi*sal+cthxi*cph*cal
           ld=acos(mu)*rad2deg

C     Converting variables into dimensional quantities and
C     storing into files for further computing

           vp=v*cal*sph         
           
           if(zr.ge.0.0d0) then
              ivisible=ivisible+1
              write(1,1001) xr,yr,th/deg2rad,ph/deg2rad,
     +             areapn,vp,mu
              write(2,1002) xr,yr,t_k,lg,ld,vp,areapn,xi
 1001         format(7(1pe14.6,2x))
 1002         format(8(1pe13.6,2x))
              
C     Computing an average value of T and log g (step 1)
              
              tareapn=t_k*areapn
              lgareapn=lg*areapn
              sumtareapn=sumtareapn+tareapn
              sumlgareapn=sumlgareapn+lgareapn
              sumareapn=sumareapn+areapn   
           end if
           
 6      continue
      
        call bounds(t_k,lg)        
        k=4*(ilat-1)+1        
        write(3,1003) (90.0d0-ith)*deg2rad,r,t_k
        write(10,1004) ith,t_k,lg,m1,m2,m3,m4,c1,c2
        model(k)=m1
        model(k+1)=m2
        model(k+2)=m3
        model(k+3)=m4
        ktot=4*(ilat-1)+4
C        call bounds(t_k_vz,lg)
C        write(11,1004) ith,t_k_vz,lg,m1,m2,m3,m4,c1,c2
 1003   format(f9.6,1x,f8.6,1x,f10.3)
 1004   format(1x,f6.2,1x,f7.1,1x,f5.3,4(1x,a10),1x,2(1x,f6.4))
        
 5    continue
      
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=10)
C      close(unit=11)

C     Computing an average value of T and log g (step 2)

      tareapn=sumtareapn/sumareapn
      lgareapn=sumlgareapn/sumareapn

C     Checking the value of v(equator)

      veq=omega*req_si/1000.0d0

C     Luminosity in solar units

      lum=lum/lsun

C     Mean value of the variable r_si and computation of
C     of the corresponding "mean" effective temperature from
C     the luminosity and the average r_si

      rsimean=rsimean/float(irsi)
      teffmean=((lum*lsun)/4.0d0/pi/sigma/rsimean**2)**0.25d0

C     Value of log g using the mean value of the variable r_si

      lgmean=4.44d0+log10(mstar)-2.0d0*log10(rsimean/rsun)
      
C     Building the name of the output spectrum that will 
C     be computed with spec3d.f, using the values of the
C     file input.star3d
C
      call last(stin,ii)
      call last(stmstar,im)
      if(tpole.lt.10000.0d0) write(sttpole,1111) int(tpole)
      if(tpole.ge.10000.0d0) write(sttpole,1112) int(tpole)
 1111 format(i4)
 1112 format(i5)
      call last(sttpole,it)
      call last(streq,ir)
      call last(stom,io)
      
      ivsini=int(veq*cal)
      if(veq*cal-ivsini.ge.0.5d0) ivsini=ivsini+1
      write(stvsini,1113) ivsini
 1113 format(i3)
      if(ivsini.le.9) stfinalvsini='00'//stvsini(3:3)
      if(ivsini.ge.10.and.ivsini.le.99) stfinalvsini='0'//stvsini(2:3)
      if(ivsini.ge.100) stfinalvsini=stvsini(1:3)
      
      outspec='in'//stin(1:ii)//'_m'//stmstar(1:im)//'_tpole'
     *    //sttpole(1:it) //'_req'//streq(1:ir)//'_mh'//metal
     *     //'_om'//stom(1:io)//'_v'//stfinalvsini//'.a'

      outpars='in'//stin(1:ii)//'_m'//stmstar(1:im)//'_tpole'
     *    //sttpole(1:it) //'_req'//streq(1:ir)//'_mh'//metal
     *     //'_om'//stom(1:io)//'_v'//stfinalvsini//'.params'

C     Module to print input and output parameters into file outpars

      open(unit=12,file=outpars,form='formatted',status='unknown')
      
      write(12,1100)
 1100 format(/,52('-'),/,'Input parameters:',/)
      write(12,1102) stsiphmer
      write(12,1103) stsithpar     
      write(12,1101) stin
      write(12,1104) stmstar 
      write(12,1105) sttpole 
      write(12,1106) streq   
      write(12,1107) stom    
      write(12,1108) smetal
 1101 format(' Inclination:    ',a10,' deg')
 1102 format(' Delta phi:      ',a10,' deg')
 1103 format(' Delta theta:    ',a10,' deg',/)
 1104 format(' M_star:         ',a10,' M_sun')
 1105 format(' T_pole:         ',a10,' K')
 1106 format(' R_equator:      ',a10,' R_sun')
 1107 format(' omega/omega_k:  ',a10)
 1108 format(' Metallicity:   ',a10,' ([M/H])')

      write(12,1110) ocrit,int(tpole*ratiot),int(tpole),ratiot,ratior,
     +     ming,maxg,lum,veq,veq*cal,int(tareapn),lgareapn,
     +     int(teffmean),lgmean,sumarea
 1110 format(/,52('-'),/,
     +     ' Surface  parameters:',//,
     +     ' omega_c:               ',1pe10.4,' s-1',/,0p,
     +     ' T_equator:            ',i5,' K',/,
     +     ' T_pole:               ',i5,' K',/,
     +     ' T_equator/T_pole:     ',f7.4,/,
     +     ' R_equator/R_pole:     ',f6.3,/,
     +     ' log g_equator:        ',f6.3,/,
     +     ' log g_pole:           ',f6.3,/,
     +     ' L_star:               ',f8.2,' L_sun',/,
     +     ' v_equator:            ',f5.1,' km/s',/,
     +     ' v_equator*sin i:      ',f5.1,' km/s',//,
     +     ' Average T and log g, weighting the values by areapn:',/,
     +     ' T=',i5,' K,   log g=',f5.3,//,
     +     ' T and log g from L_star, M_star and the average of R:',/,
     +     ' T=',i5,' K    log g=',f5.3,//,
     +     ' Surface of the star (r normalized to R(eq))=',f6.3,/,
     +     ' [A=4*pi=12.566 for a spherical star of r=1.0]')

      close(unit=12)
      
C     Writing a small input file for spec3d.f
      
      open(unit=13,file='input.spec3d',
     +             form='formatted',status='unknown')
      write(13,1114) metal
 1114 format(a4)
      write(13,1115) outspec
 1115 format(a60)
      write(13,1116) ilat
 1116 format(i3)
      write(13,1117) ivisible
 1117 format(i6)
      close(unit=13)

      write(*,1120) outpars
 1120 format(/,40('-'),/,
     +     'Input and surface parameters are in file: ',/,40('-'),/,
     +     a70,//,
     +     26('-'),/,'Additional files generated:',/,26('-')/,
     +     'input.spec3d: a few inputs required by program spec3d',/,
     +     'spec3d.lat:   information to be used by spec3d to',/,
     +     '              compute spectra at each latitudinal strip',/,
     +     'tiles.a:      parameters for each individual cell',/,
     +     'xyparams.a:   several parameters of the star to plot',/,
     +     '              3D graphs in 2D',/)
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
      

