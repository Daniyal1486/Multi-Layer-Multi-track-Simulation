      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1 TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME
      REAL thermcon,alpha,vel,R0,R,pi,dist
     & x,y,z,power,lamda,abs,AI,M,Eff,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,u
      power=300.0
      vel=0.5
C ********************
      Eff=2. ! 11 this added t othis code
      pi=3.141592654
      abs=0.41
      R0=0.00100
      u=0.01
      t1=u/vel
      t2=t1+u/vel
      t3=t2+u/vel
      t4=t3+u/vel
      t5=t4+u/vel
      t6=t5+u/vel
      t7=t6+u/vel
      t8=t7+u/vel
      t9=t8+u/vel
      t10=t9+u/vel
      AI=(Eff*power)/(pi*(R0*R0)) ! 11
      IF (Time(1) .GT. 0 .AND. TIME(1) .LE. t1) THEN          ! first pass
      dist=vel*Time(1)
      x=COORDS(1)-dist
      y=COORDS(3)
      z=COORDS(2)-0.0006
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t1 .AND. TIME(1) .LE. t2) THEN    ! second pass
      dist=vel*t1+vel*(Time(1)-t1)
      x=COORDS(1)+0.01-dist
      y=COORDS(3)
      z=COORDS(2)-0.0012
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t2 .AND. TIME(1) .LE. t3) THEN   ! third pass
      dist=vel*t1+vel*(t2-t1)+vel*(Time(1)-t2)
      x=COORDS(1)+0.02-dist
      y=COORDS(3)
      z=COORDS(2)-0.0018
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t3 .AND. TIME(1) .LE. t4) THEN   ! fourth pass           rotate
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(Time(1)-t3)
      x=COORDS(1)-0.004
      y=COORDS(3)
      z=COORDS(2)+0.03-dist
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t4 .AND. TIME(1) .LE. t5) THEN   ! fifth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(Time(1)-t4)
      x=COORDS(1)-0.0046
      y=COORDS(3)
      z=COORDS(2)+0.04-dist
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t5 .AND. TIME(1) .LE. t6) THEN   ! sixth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(Time(1)-t5)
      x=COORDS(1)-0.0052
      y=COORDS(3)
      z=COORDS(2)+0.05-dist
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t6 .AND. TIME(1) .LE. t7) THEN   ! seventh pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(Time(1)-t6)
      x=COORDS(1)-0.0058
      y=COORDS(3)
      z=COORDS(2)+0.06-dist
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t7 .AND. TIME(1) .LE. t8) THEN   ! eighth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(t7-t6)+vel*(Time(1)-t7)
      x=COORDS(1)-0.0064
      y=COORDS(3)
      z=COORDS(2)+0.07-dist
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t8 .AND. TIME(1) .LE. t9) THEN   ! ninth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(t7-t6)+vel*(t8-t7)+vel*(Time(1)-t8)
      x=COORDS(1)-0.0070
      y=COORDS(3)
      z=COORDS(2)+0.08-dist
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t9 .AND. TIME(1) .LE. t10) THEN   ! tenth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(t7-t6)+vel*(t8-t7)+vel*(t9-t8)+vel*(Time(1)-t9)
      x=COORDS(1)-0.0076
      y=COORDS(3)
      z=COORDS(2)+0.09-dist
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      FLUX(2)=0.
      JLTYP=0. 
      END IF
      IF (Time(1) .GT. t10) THEN
      FLUX(1)=0
      FLUX(2)=0.
      JLTYP=0.
      END IF
      RETURN
      END