      SUBROUTINE BETACA(NVAR,NPAR,NCOR,CORREL,NWORK,ITYPE,SPAR,
     -   BWORK,X,ALFA,U,PSENTIV,P,PSENTIP,STONAM,PARNAM,MAXIT,
     -   BETA,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,Ela,Ruo,Sigc,Sigt,
     -   Area,Exc,WYYS,RuoL,II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag
     -   ,NOFG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/GLSTOC/JNAME,IBNAME

      CHARACTER*8 STONAM(NVAR),PARNAM(NPAR),JNAME
      DIMENSION CORREL(NCOR,NCOR),ITYPE(NVAR),SPAR(NVAR,4),
     -    BWORK(NWORK),X(NVAR),ALFA(NVAR),U(NVAR),
     -    PSENTIV(NVAR,4),P(NPAR),PSENTIP(NPAR)
c   结构分析需用的变量
      REAL*8 XX,YY,Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,WY,YL
      DIMENSION XX(nJD+1),YY(nJD+1),Ela(nDY+1),Ruo(nDY+1),Sigc(nDY+1),
     -          Sigt(nDY+1),Area(nDY+1),Exc(nExc+1),WYYS(nWYYS+1),
     -          RuoL(nDY+1),WY(nJD*2+1),YL(nDY+1)
	INTEGER*4 II,JJ,IP,JP,IExc,IWYYS
      DIMENSION II(nDY+1),JJ(nDY+1),IP(nJD+1),JP(nJD+1),
     -	      IExc(nExc+1),IWYYS(nWYYS+1)
      INCLUDE 'COMMY.INC'
C
      EPS = 1.D-14

      IWRITE=2  !可改，为0表不写入MRE2中和屏幕输出，为1表写入，为2表写入更多
      
	IF(NOFG.EQ.0)THEN
      DO 46 I=1,NVAR
   46 U(I)=0.0D0    !初值，可改
c	u(1)=-0.25	
c	u(2)=-0.25	
c	u(3)=1.0	
c	u(4)=0.25 !!XULIN4_3见文献31
	ENDIF
C
      WRITE(MRE1,*)
      WRITE(MRE1,*)
C
	CALL BETAE0 (NVAR,NPAR,U,BETA,ALFA,X,SPAR,ITYPE,NCOR,
     -    CORREL,BWORK,NWORK,EPS,IWRITE,PSENTIV,IER,P,STONAM,
     -    PARNAM,PSENTIP,MAXIT,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -    Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,II,JJ,IP,JP,
     -    IExc,IWYYS,WY,YL,IFlag,NOFG)
C
C    .WRITE RESULT OF RELIABILITY ANALYSIS TO FILE .RE1.
C
       print *, '!!!!!!!!!!!!!!!!!BETACA!!!!!!!!!!!'
      DO  I = 1,NVAR
      print *, ' I', I
      print *, 'PSENTIV(I,1,2)',PSENTIV(I,1),PSENTIV(I,2)
      
      END DO
      print *, '!!!!!!!!!!!!!!!!!BETACA!!!!!!!!!!!'
      IF (IER.NE.0) THEN
	WRITE(MRE1,4000) IER
      ELSE
	WRITE(MRE1,4010) BETA,PHI(-BETA)
	WRITE(MRE1,1000)
	DO 40 I=1,NVAR
	  WRITE(MRE1,905) I,U(I),X(I)
   40   CONTINUE
	WRITE(MRE1,4020) 
	DO 62 I = 1,NVAR
   62	WRITE(MRE1,4022) STONAM(I),ALFA(I),ALFA(I)**2,
     -                PSENTIV(I,1),PSENTIV(I,2)
	WRITE(MRE1,4040)
	DO 63 I = 1,NPAR
   63	WRITE(MRE1,4030) PARNAM(I),PSENTIP(I)
C
c      IF(IWRITE.GE.2)THEN
c	WRITE(*,4010) BETA,PHI(-BETA)
c	WRITE(*,4020) 
c	DO 64 I = 1,NVAR
c   64   WRITE(*,4022) STONAM(I),ALFA(I),ALFA(I)**2,
c     -                PSENTIV(I,1),PSENTIV(I,2)
c	WRITE(*,4040)
c	DO 65 I = 1,NPAR
c   65   WRITE(*,4030) PARNAM(I),PSENTIP(I)
c      ENDIF

      ENDIF

      CALL QDATE(MRE1)
C
  905 FORMAT(I5,1F9.4,1PE12.4)
 1000 FORMAT('    I   ',7X,'U',8X,'X')
 4000 FORMAT(' Error in beta-algorithm : IER = ',I4)
 4010 FORMAT(' Reliability index: ',F8.3,
     -       '        Probability of failure:',
     -        E12.4)
 4020 FORMAT(/,'              Alfa-coeff.  Alfa**2 ',
     -         '   d beta/d mu    d beta/d sigma')
 4022 FORMAT(' ',A8,'  ',F12.3,' ',F12.3,'  ',E14.6,' ',E14.6)
 4030 FORMAT(' ',A8,'  ',E14.6,'  ',E14.6,'  ',F8.3)
 4040 FORMAT(/,'             d beta/d p')
      RETURN
      END
C
C
      SUBROUTINE BETAE0(NVAR,NPAR,U,BETA,ALFA,X,SPAR,
     -    ITYPE,NCOR,CORREL,BWORK,NWORK,EPS,IWRITE,PSENTI,IER,P,
     -    STONAM,PARNAM,PSENTIP,MAXIT,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -    Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,II,JJ,IP,JP,
     -    IExc,IWYYS,WY,YL,IFlag,NOFG)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*8 STONAM(NVAR),PARNAM(NPAR)
      COMMON/CMACHE/EPSS
    
      INCLUDE 'COMMY.INC'
      DIMENSION X(NVAR),U(NVAR),SPAR(NVAR,4),CORREL(NCOR,NCOR),
     -          ITYPE(NVAR),ALFA(NVAR),
     -          PSENTI(NVAR,4),BWORK(NWORK),P(NPAR),PSENTIP(NPAR)
      DIMENSION AX(200)
c   结构分析需用的变量
      REAL*8 XX,YY,Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,WY,YL
      DIMENSION XX(nJD+1),YY(nJD+1),Ela(nDY+1),Ruo(nDY+1),Sigc(nDY+1),
     -          Sigt(nDY+1),Area(nDY+1),Exc(nExc+1),WYYS(nWYYS+1),
     -          RuoL(nDY+1),WY(nJD*2+1),YL(nDY+1)
	INTEGER*4 II,JJ,IP,JP,IExc,IWYYS
      DIMENSION II(nDY+1),JJ(nDY+1),IP(nJD+1),JP(nJD+1),
     -	      IExc(nExc+1),IWYYS(nWYYS+1)

      IER=0
      EPSS=EPS
      NDP=NPAR+NVAR
      M1=1
      M2=M1+NVAR
      M3=M2+NVAR*NVAR

      IF(M3.GT.NWORK)THEN
	WRITE(*,*) ' NWORK TOO SMALL 1',M3,NWORK
	STOP
      ENDIF
C
      MA0=1
      MA1=MA0+NDP
      IF(MA1.GT.200)THEN
	WRITE(*,*) ' AX TOO SMALL'
	STOP
      ENDIF
C
      CALL BETAE1(NVAR,NPAR,X,U,SPAR,BWORK(M1),ITYPE,CORREL,NCOR,
     -  BETA,G0,SQ2,IWRITE,MRE2,ITER,IER,P,STONAM,PARNAM,AX(MA0),
     -  NDP,BWORK(M2),MAXIT,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,Ela,Ruo,
     -  Sigc,Sigt,Area,Exc,WYYS,RuoL,II,JJ,IP,JP,IExc,IWYYS,WY,YL,
     -  IFlag,NOFG)
c
      CALL SENANA(NVAR,NPAR,X,U,SPAR,ITYPE,CORREL,NCOR,
     -  BETA,SQ2,PSENTI,IER,PSENTIP,AX(MA0),NDP,BWORK(M1))
      IF(IWRITE.GE.1)THEN
	PF=PHI(-BETA)
	WRITE(MRE2,1060) ITER,BETA,PF
	WRITE(MRE2,1000)
	DO 40 I=1,NVAR
	ALFA(I)=U(I)/BETA
	IF(ABS(ALFA(I)).GT.0.001)THEN
	  WRITE(MRE2,905) I,ALFA(I),U(I),X(I)
	ENDIF
   40   CONTINUE
	  WRITE(MRE2,1070)
	  DO 45 I=1,NVAR
   45     WRITE(MRE2,1080) I,(PSENTI(I,J),J=1,4)
	  WRITE(MRE2,1075)
	  DO 55 I=1,NPAR
   55     WRITE(MRE2,1080) I,PSENTIP(I)
	ENDIF
C
  905 FORMAT(I5,2F9.4,1PE12.4)
 1000 FORMAT('    I   ALFA',7X,'U',8X,'X')
 1060 FORMAT(/,' ITERATIONS : ',I5,' BETA : ',F10.5,' PF :',1PE12.4)
 1070 FORMAT(/,
     -' SENSITIVITY COEFFICIENTS FOR STATISTICAL PARAMETERS:',
     -' (D BETA / D PAR,I) ',/,
     -'        PAR,1       PAR,2       PAR,3       PAR,4')
 1075 FORMAT(/,
     -' SENSITIVITY COEFFICIENTS FOR CONSTANTS:(D BETA / D PAR,I) ',/,
     -'        PAR,1       ')
 1080 FORMAT(I4,4(1PE12.4))
      RETURN
      END
C
C
      SUBROUTINE BETAE1(N,NP,X,U,SPAR,DuG,IT,RL,NCOR,BETA,G0,
     -  SQ2,IWRITE,MRE2,ITER,IER,P,STONAM,PARNAM,DP,NDP,XJU,
     -  MAXIT,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,Ela,Ruo,Sigc,Sigt,
     -  Area,Exc,WYYS,RuoL,II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
C********************************************************************
C
C     .CALCULATION OF RELIABILITY INDEX.
C
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 BETAK
      CHARACTER*8 STONAM(N),PARNAM(NP)
      DIMENSION X(N),U(N),SPAR(N,4),DuG(N),RL(NCOR,NCOR)
      DIMENSION IT(N),P(NP),DP(NDP), XJU(N,N)
      COMMON/CMACHE/EPS
	DIMENSION ALFA(N), gu(N)

c   结构分析需用的变量
      REAL*8 XX,YY,Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,WY,YL
      DIMENSION XX(nJD+1),YY(nJD+1),Ela(nDY+1),Ruo(nDY+1),Sigc(nDY+1),
     -          Sigt(nDY+1),Area(nDY+1),Exc(nExc+1),WYYS(nWYYS+1),
     -          RuoL(nDY+1),WY(nJD*2+1),YL(nDY+1)
	INTEGER*4 II,JJ,IP,JP,IExc,IWYYS
      DIMENSION II(nDY+1),JJ(nDY+1),IP(nJD+1),JP(nJD+1),
     -	      IExc(nExc+1),IWYYS(nWYYS+1)
C
      IER=0
      ZE=0.D0
      ON=1.D0
      
      IF(IWRITE.GE.2)THEN
	WRITE(MRE2,1000)
	DO 6 I=1,N
	 WRITE(MRE2,905) I,0.,U(I)
    6 CONTINUE
	 WRITE(MRE2,*)
	ENDIF
C
      DO 50 ITER=1,MAXIT
	IJK = 0
	PRINT*,'BETAE1 START'
      CALL GFUNCU(N,NP,U,DuG,SS2,G0,NCOR,RL,IT,SPAR,X,P,STONAM,PARNAM,
     -    DP,NDP,XJU,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,Ela,Ruo,Sigc,Sigt,
     -    Area,Exc,WYYS,RuoL,II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
      sq2=dsqrt(ss2)
      PRINT*,'sq2',sq2
	dw1=0.d0
	S=ZE
	DO 15 I=1,N
	 dw1 = dw1 + DuG(i) * u(i)
   15 CONTINUE
      BETAK= (g0 - dw1)/sq2
      PRINT*,'BETAK',BETAK
      PRINT*,'dw1',dw1
       PRINT*,'BETAE1 BEFORE-U(I)',(U(I),I=1,N)
       do 116 i=1,n
       u(i) = (-BETAK*DuG(i)) /sq2  
  116 CONTINUE
        PRINT*,'-----------------------------------'
       PRINT*,'BETAE1 BETAK-U(I)',(U(I),I=1,N)
       PRINT*,'-----------------------------------'
      do 16 i=1,n
       u(i) = ((dw1- g0)*DuG(i)) / ss2
   16 CONTINUE
100	dw2 = 0.D0
	do 17 i = 1,n 
17	dw2 = dw2 + u(i)*u(i)
      PRINT*,'-----------------------------------'
      PRINT*,'BETAE1 U(I)',(U(I),I=1,N)
      PRINT*,'-----------------------------------'
      PRINT*,'dw2',dw2
      CALL ZUTOX(N,NCOR,U,RL,SPAR,IT,X)
      PRINT*,'BETAE1 X(I)',(X(I),I=1,N)
      CALL GFUNCX(N,X,NP,P,G,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -            Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
     -            II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)

1111	 temp = dsqrt(dw2)
       PRINT*,'BETAE1 temp',temp
       PRINT*,'BETAE1 beta',beta
	if( dw1 > 0.0 ) temp = -temp
	PRINT*,'BETAE1 abs(beta/temp-1.0)',abs(beta/temp-1.0)
      if(abs(beta/temp-1.0)<1.d-4) then !收敛准则，可改
       
	   beta = temp
	   PRINT*,'BETAE1 A beta',beta
		goto 60
	 end if
!	 DO 15 I=1,N
!       ALFA(I) = u(i) - g0*DuG(i)/ss2
!	 dw1 = dw1 + DuG(i) * u(i)
!   15 CONTINUE
!      do 16 i=1,n
!	 gu(i) = dw1*DuG(i) / ss2 - u(i)
!!      u(i) = dw1*DuG(i) / ss2 - g0*DuG(i)/ss2 !瞎改
!   16 CONTINUE
!       PRINT*,'BETAE1 gu(i)',(gu(i),I=1,N)
!100	dw2 = 0.D0
!	do 17 i = 1,n
!	  u(i) = alfa(I) + 0.5D0**IJK* gu(i) 
!17	dw2 = dw2 + u(i)*u(i)
!      PRINT*,'-----------------------------------'
!       PRINT*,'BETAE1 u(i)',(u(I),I=1,N)
!       PRINT*,'-----------------------------------'
!       PRINT*,'BETAE1 dw2',dw2
!      CALL ZUTOX(N,NCOR,U,RL,SPAR,IT,X)
!      PRINT*,'BETAE1 -J'
!      PRINT*,'BETAE1 X(I)',(X(I),I=1,N)
!      CALL GFUNCX(N,X,NP,P,G,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
!     -            Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
!     -            II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
!			print *, 'alfa,gu,u',alfa,gu,u
!			print *, 'IJK, g1,g ',IJK,g,g0
!      if( abs(g0) > 1.0d-10.and.abs(g) > 0.5*abs(g0) ) then  !这里若灰掉表示没有用许林提出的改进
!	  write(mre2,920) IJK,u
!c      改:加一避免死循环的语句
!	  if(IJK>10)goto 1111
!c      改:加一避免死循环的语句
!	  IJK = IJK + 1
!	  goto 100
!	end if
!1111	 temp = dsqrt(dw2)
!      PRINT*,'BETAE1 temp',temp
!	if( dw1 > 0.0 ) temp = -temp
!      if(abs(beta/temp-1.0)<1.d-4) then !收敛准则，可改
!      PRINT*,'BETAE1 abs(beta/temp-1.0)',abs(beta/temp-1.0)
!	   beta = temp	   
!		goto 60
!	 end if
	BETA=temp
	
      IF(IWRITE.GE.2)THEN
	 PF=PHI(-BETA)
	 PRINT*,'BETAE1 PF',PF
	 WRITE(MRE2,1060) ITER,BETA,PF
	 WRITE(*,1060) ITER,BETA,PF
	 WRITE(MRE2,1000)
	 DO 35 I=1,N
	 WRITE(MRE2,905) I,DuG(I)/sq2,U(I),X(I)
   35  CONTINUE
	 WRITE(MRE2,*)
      ENDIF
   50 CONTINUE
   60 CONTINUE
      
      CALL ZUTOX(N,NCOR,U,RL,SPAR,IT,X)
      PRINT*,'------------------'
      WRITE(*,1060) ITER,BETA,PF
      PRINT*,'FINISH BETAE1'
      WRITE(MRE2,*) '*******  STOP IN BETAE1!  *********'
  905 FORMAT(I5,2F9.4,1PE12.4)
  920 FORMAT(' Inner ITER(IJK) = ',I3,3X,' U = ',100F9.5)
 1000 FORMAT('    I    ALFA',4X,'U',11X,'X')
 1060 FORMAT('Outter ITERATION= ',I5,' BETA : ',F10.5,' PF :',1PE12.4)
      RETURN
      END
C
C
	SUBROUTINE SENANA(N,NP,X,U,SPAR,IT,CORREL,NCOR,
     -  BETA,SQ2,PSENTI,IER,PSENTIP,DP,NDP,DuG)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),U(N),SPAR(N,4),CORREL(NCOR,NCOR),IT(N),DuG(N),
     -          PSENTI(N,4),DP(NDP),PSENTIP(NP),STAT(4),UN(N),STAP(4)
      IER=0
      DEPS=1.D-3
       print *,'------------------------------'
       print *,'SENANA'
       print *,'------------------------------'
      CALL ZXTOU(N,NCOR,X,CORREL,SPAR,IT,U)
      print *,'SENANA -A'
      print *,'X(I)=',(X(I),I=1,N)
      print *,'U(I)=',(U(I),I=1,N)
      DO 100 IVAR=1,N
      DO 3 I=1,4
    3 STAP(I)=SPAR(IVAR,I)
!      print *,' STAP(I)=',( STAP(I),I=1,4)
      IF(IT(IVAR).GT.0)THEN
	DO 90 J=1,4
	IF(J.GT.1 .AND. SPAR(IVAR,J).EQ.0.D0)THEN
	  PSENTI(IVAR,J)=0.D0
	  ELSE
	  CALL PARCAI(IT(IVAR),STAT(1),STAT(2),STAT(3),STAT(4),
     -            SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4))
	  CALL PARCAL(MRE2,IT(IVAR),STAT(1),STAT(2),STAT(3),
     -      STAT(4),SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),
     -      SPAR(IVAR,4),2)
	  CALL ZXTOU(N,NCOR,X,CORREL,SPAR,IT,U)
	  print *,'SENANA -B'
      print *,'U(I)=',(U(I),I=1,N)
      print *,' STAT(I)=',( STAT(I),I=1,4)
	  DSTAT=STAT(J)
!	   print *,' J,STAT(J)=',J,STAT(J)
	  IF(DSTAT.NE.0.D0)THEN
	    DS=DEPS*DABS(DSTAT)
	  ELSE
	    DS=DEPS
	  ENDIF
	  STAT(J)=DSTAT+DS
	  print *,'SENANA -C'
	  print *,' STAT(I)=',( STAT(I),I=1,4)
	  CALL PARCAL(MRE2,IT(IVAR),STAT(1),STAT(2),STAT(3),
     -      STAT(4),SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),
     -      SPAR(IVAR,4),2)
	  CALL ZXTOU(N,NCOR,X,CORREL,SPAR,IT,UN)
	  print *,'SENANA -D'
      print *,'UN(I)=',(UN(I),I=1,N)
 
	  DB=0.D0
	  DO 10 I=1,N
   10     DB=DB+U(I)*(UN(I)-U(I))/DS
	  PSENTI(IVAR,J)=DB/BETA
	   print *,'IVAR,J,PSENTI(IVAR,J)=',IVAR,J,PSENTI(IVAR,J)
	  STAT(J)=DSTAT
	  CALL PARCAL(MRE2,IT(IVAR),STAT(1),STAT(2),STAT(3),
     -           STAT(4),
     -           SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4),2)
	ENDIF
   90   CONTINUE
      ELSE
       print *,'SENANA -E'
	PSENTI(IVAR,2) = 0.D0
	PSENTI(IVAR,3) = 0.D0
	PSENTI(IVAR,4) = 0.D0
	PSENTI(IVAR,1) = DP(IVAR)/SQ2
      ENDIF
      DO 93 I=1,4
   93 SPAR(IVAR,I)=STAP(I)
  100 CONTINUE
      CALL ZUTOX(N,NCOR,U,CORREL,SPAR,IT,X)
      DO 110 IPAR=1,NP
  110 PSENTIP(IPAR) = DP(N+IPAR)/SQ2
      CALL ZUTOX(N,NCOR,U,CORREL,SPAR,IT,X)
      RETURN
      END