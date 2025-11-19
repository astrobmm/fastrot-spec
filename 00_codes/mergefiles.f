      program mergefiles

      implicit none

      character*120 st1
      character*45  st2
      character*13  stxr,styr,sttm,stlg,stvr,stld
      
      open(unit=1,file='xyparams.a',form='formatted',status='old')
      open(unit=2,file='limbdark.a',form='formatted',status='old')
      open(unit=3,file='4plots.a',form='formatted',status='unknown')

C         10        20        30        40        50        60        70        80        90       100       110     118
C 1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678      
C -4.141236E-01  -9.015491E-01   8.502148E+03   3.700267E+00   8.438375E+01  -1.762929E+01   1.164587E-04   5.552018E+00

 1    read(1,100,end=1000) st1
      stxr=st1(1:13)
      styr=st1(16:28)
      sttm=st1(32:43)
      stlg=st1(47:58)
      stvr=st1(76:88)
 100  format(a120)

C         10        20        30        40 43 
C 1234567890123456789012345678901234567890123    
C -4.141236E-01  -9.015491E-01   4.265046E-01

 1000 read(2,200,end=2000) st2
      stld=st2(32:43)
 200  format(a45)

      write(3,300) stxr,styr,sttm,stlg,stvr,stld
 300  format(6(a13,1x))

      goto 1

 2000 close(unit=1)
      close(unit=2)
      close(unit=3)

      end
