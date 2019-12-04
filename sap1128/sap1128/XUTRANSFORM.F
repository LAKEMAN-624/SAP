      SUBROUTINE PARCAL(MRE2,IT,P1,P2,P3,P4,S1,S2,S3,S4,ISW)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRUNCN/XMU,XVA,A,B,MUU,VAR
      COMMON/CMACHE/EPSS
      DIMENSION X(2)
      REAL*8 MUU
      ZE=0.D0
      ON=1.D0
      TW=2.D0
      PI=3.1415926535898D0
      S1=ZE
      S2=ZE
      S3=ZE
      S4=ZE
      ITT=IABS(IT)
      IF(ITT.EQ.0)GOTO 5
      IF(ITT.EQ.1)GOTO 10
      IF(ITT.EQ.2)GOTO 20
      IF(ITT.EQ.3)GOTO 30
      IF(ITT.EQ.4)GOTO 40
      IF(ITT.EQ.5)GOTO 50
      IF(ITT.EQ.6)GOTO 60
      IF(ITT.EQ.7)GOTO 70
      IF(ITT.EQ.8)GOTO 80
      IF(ITT.EQ.9)GOTO 90
      IF(ITT.EQ.10)GOTO 100
      IF(ITT.EQ.11)GOTO 110
C
      WRITE(*,*) ' ERROR IN DISTRIBUTION TYPE'
      STOP
    5 IF(ISW.EQ.1)THEN
	WRITE(MRE2,905) J,P1
      ENDIF
  905 FORMAT(I5,' DETERMISTIC ',' P   = ',1PE12.5)
      S1=P1
      S2=ZE
      GOTO 7777
   10 IF(IT.GT.0)THEN
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,910) J,P1,P2
	ENDIF
  910   FORMAT(I5,' UNIFORM     ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
	S1=P1-SQRT(3.D0)*P2
	S2=P1+SQRT(3.D0)*P2
      ELSE
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,911) J,P1,P2
	ENDIF
  911   FORMAT(I5,' UNIFORM     ',' A   = ',1PE12.5,' B    =',1PE12.5)
	S1=P1
	S2=P2
      ENDIF
      GOTO 7777
   20 IF(ISW.EQ.1)THEN
	WRITE(MRE2,920) J,P1,P2
      ENDIF
  920 FORMAT(I5,' NORMAL      ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
      S1=P1
      S2=P2
      GOTO 7777
   30 IF(IT.GT.0)THEN
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,930) J,P1,P2
	ENDIF
  930   FORMAT(I5,' LOGNORMAL   ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
	S2=SQRT(DLOG((P2/P1)**2+1.D0))
	S1=DLOG(P1)-.5D0*S2**2
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,931) S1,S2
	ENDIF
  931   FORMAT(5X,'             ',' BETA= ',1PE12.5,' ALFA=',1PE12.5)
      ELSE
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,932) J,P1,P2
	ENDIF
  932   FORMAT(I5,' LOGNORMAL   ',' BETA= ',1PE12.5,' ALFA=',1PE12.5)
	S1=P1
	S2=P2
      ENDIF
      GOTO 7777
   40 CONTINUE
      IF(IT.GT.0)THEN
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,940) J,P1,P2
	ENDIF
  940   FORMAT(I5,' EXPONENTIAL ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
       S1=ON/P2
	S2=P1-P2
      ELSE
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,941) J,P1,P2
	ENDIF
  941   FORMAT(I5,' EXPONENTIAL ',' LAMB= ',1PE12.5,' EPS = ',1PE12.5)
	S1=P1
	S2=P2
      ENDIF
      GOTO 7777
   50 CONTINUE
      IF(IT.GT.0)THEN
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,950) J,P1,P2
	ENDIF
  950   FORMAT(I5,' GAMMA       ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
	S1=P1/P2**2
	S2=P1**2/P2**2
      ELSE
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,951) J,P1,P2
	ENDIF
  951   FORMAT(I5,' GAMMA       ',' LAMB= ',1PE12.5,' K   = ',1PE12.5)
	S1=P1
	S2=P2
      ENDIF
      GOTO 7777
   60 IF(IT.GT.0)THEN
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,960) J,P1,P2
	ENDIF
  960   FORMAT(I5,' GUMBEL-MAX  ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
	S1=PI/(SQRT(6.D0)*P2)
	S2=P1-.57722D0/S1
      ELSE
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,961) J,P1,P2
	ENDIF
  961   FORMAT(I5,' GUMBEL-MAX   ',' ALFA= ',1PE12.5,' U    =',1PE12.5)
	S1=P1
	S2=P2
      ENDIF
      GOTO 7777
   70 IF(IT.GT.0)THEN
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,970) J,P1,P2
	ENDIF
  970   FORMAT(I5,' GUMBEL-MIN  ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
	S1=PI/(SQRT(6.D0)*P2)
	S2=P1+.57722D0/S1
      ELSE
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,971) J,P1,P2
	ENDIF
  971   FORMAT(I5,' GUMBEL-MIN   ',' ALFA= ',1PE12.5,' U    =',1PE12.5)
	S1=P1
	S2=P2
      ENDIF
      GOTO 7777
   80 IF(IT.GT.0)THEN
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,980) J,P1,P2
	ENDIF
  980 FORMAT(I5,' FRECHET     ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
	VY=P2/(P1)
	B1=VY**(-1.09)
	B2=B1*1.1
	F11=GAMMA(ON-TW/B1)
	F12=GAMMA(ON-ON/B1)
	F1=F11/F12**2-ON-VY**2
	F21=GAMMA(ON-TW/B2)
	F22=GAMMA(ON-ON/B2)
	F2=F21/F22**2-ON-VY**2
	F20=F2
C
	DO 82 I=1,1000
	AL=(F2-F1)/(B2-B1)
	B3=B2-F2/AL
	B1=B2
	F1=F2
	B2=B3
	F21=GAMMA(ON-TW/B2)
	F22=GAMMA(ON-ON/B2)
	F2=F21/F22**2-ON-VY**2
	IF(ABS(F2/F20).LT.1.D-6)GOTO 83
   82   CONTINUE
   83   S1=B2
	S2=(P1)/F22
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,981) S1,S2
	ENDIF
  981 FORMAT('                  ',' K   = ',1PE12.5,' U   = ',1PE12.5)
C
      ELSE
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,982) J,P1,P2
	ENDIF
  982 FORMAT(I5,' WEIBULL     ',' K   = ',1PE12.5,' U   = ',1PE12.5)
	S1=P1
	S2=P2
	Q1=GAMMA(ON-ON/S1)
	Q2=GAMMA(ON-TW/S1)
	PP1=(S2)*Q1
	PP2=(S2)*SQRT(Q2-Q1**2)
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,984) PP1,PP2
	ENDIF
      ENDIF
  984 FORMAT('                  ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
      GOTO 7777
   90 IF(IT.GT.0)THEN
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,990) J,P1,P2
	ENDIF
  990 FORMAT(I5,' WEIBULL     ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
	VY=P2/P1
	B1=VY**(-1.09D0)
	B2=B1*1.1D0
	F11=GAMMA(ON+TW/B1)
	F12=GAMMA(ON+ON/B1)
	F1=F11/F12**2-ON-VY**2
	F21=GAMMA(ON+TW/B2)
	F22=GAMMA(ON+ON/B2)
	F2=F21/F22**2-ON-VY**2
	F20=F2
C
	DO 92 I=1,1000
	AL=(F2-F1)/(B2-B1)
	B3=B2-F2/AL
	B1=B2
	F1=F2
	B2=B3
	F21=GAMMA(ON+TW/B2)
	F22=GAMMA(ON+ON/B2)
	F2=F21/F22**2-ON-VY**2
	IF(ABS(F2/F20).LT.1.D-6)GOTO 93
   92   CONTINUE
	WRITE(*,*) ' ERROR IN ROOT FINDING IN WEIBULL DIST. IN PARCAL'
	STOP
   93   S1=B2
	S2=P1/F22
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,991) S1,S2
	ENDIF
  991 FORMAT('                  ',' K   = ',1PE12.5,' U   = ',1PE12.5)
C
      ELSE
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,992) J,P1,P2
	ENDIF
  992 FORMAT(I5,' WEIBULL     ',' K   = ',1PE12.5,' U   = ',1PE12.5)
	S1=P1
	S2=P2
	Q1=GAMMA(ON+ON/S1)
	Q2=GAMMA(ON+TW/S1)
	PP1=(S2)*Q1
	PP2=(S2)*SQRT(Q2-Q1**2)
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,993) PP1,PP2
	ENDIF
      ENDIF
  993 FORMAT('                  ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
      GOTO 7777
  100 IF(IT.GT.0)THEN
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,1000) J,P1,P2
	ENDIF
 1000   FORMAT(I5,' RAYLEIGH    ',' MY  = ',1PE12.5,' SIG = ',1PE12.5)
	S1=P2/SQRT(TW-PI/TW)
	S2=P1-P2*SQRT(PI/(4.D0-PI))
      ELSE
	IF(ISW.EQ.1)THEN
	  WRITE(MRE2,1001) J,P1,P2
	ENDIF
 1001   FORMAT(I5,' RAYLEIGH     ',' ALFA= ',1PE12.5,' EPS  =',1PE12.5)
	S1=P1
	S2=P2
      ENDIF
      GOTO 7777
CXX
  110 IF(P2.EQ.0.D0)THEN
	  WRITE(MRE2,1111) P2
	  WRITE(*,1111) P2
 1111   FORMAT(//,' ERROR IN TRUNCATED NORMAL VARIABEL. ',
     -            ' STD.DEV HAS TO BE > 0',
     -            ' (IS = ',E12.4,')')
	  STOP
      ENDIF
      IF(P3.GE.P4)THEN
	  WRITE(MRE2,1112) P3,P4
	  WRITE(*,1112) P3,P4
 1112   FORMAT(//,' ERROR IN TRUNCATED NORMAL VARIABEL. ',
     -            ' B HAS TO BE > A',
     -            ' (A = ',E12.4,' B = ',E12.4,')')
	  STOP
      ENDIF
      XMU=P1/P2
      XVA=1.D0
      A=P3/P2
      B=P4/P2
      IF(IT.GT.0)THEN
	  IF(ISW.EQ.1)THEN
	   WRITE(MRE2,1020) J,P1,P2,P3,P4
        ENDIF
 1020   FORMAT(I5,' TRUNC. NORMAL',' MYX = ',1PE12.5,' SIGX= ',1PE12.5,
     -   ' A   = ',1PE12.5,' B    = ',1PE12.5)
      EPSS=1.D-14
      X(1) =XMU
      X(2) =1.3D0
      CALL DFPMIN(X,2,EPSS,F,ITER,IER)
      IF(ABS(MUU-XMU).GT.0.001 .OR.
     -  ABS(VAR-1.D0).GT.0.001 .OR.
     -  IER.NE.0)THEN
 	 WRITE(*,*) ' MU,SIG',XMU,SQRT(XVA)
	 WRITE(*,*) ' MU,SIG',MUU,SQRT(VAR)
	 WRITE(*,*)  
     -  ' ERROR IN CALCULATING PARAMETERS IN TRUNCATED NORMAL',
     -  ' DISTRIBUTION'
	 STOP
      ENDIF                          
	  S1=X(1)*P2
	  S2=X(2)*P2
	  S3=P3
	  S4=P4
      ELSE
	  IF(ISW.EQ.1)THEN
	   WRITE(MRE2,1021) J,P1,P2,P3,P4
	  ENDIF
 1021   FORMAT(I5,' TRUNC. NORMAL',' MY  = ',1PE12.5,' SIG = ',1PE12.5,
     -   ' A   = ',1PE12.5,' B    = ',1PE12.5)
	  S1=P1
	  S2=P2
	  S3=P3
	  S4=P4
      ENDIF
 7777 J=0
      RETURN
      END
C
C
      SUBROUTINE TRFUNC(N,X,F)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TRUNCN/XMU,XVA,A,B,MUU,VAR
      DIMENSION X(N)
      REAL*8 MU,MUU,K
      PI=3.1415926535898D0
C
      MU=X(1)
      SIG=X(2)
      AL=((A-MU)/SIG)**2/2.D0
      BL=((B-MU)/SIG)**2/2.D0
      UA=(A-MU)/SIG
      UB=(B-MU)/SIG
      K=PHI(UB)-PHI(UA)
      IF(K.LE.1.D-15)K=1.D-15
      QA=GAMMIN(1.5D0,AL)*GAMMA(1.5D0)
      QB=GAMMIN(1.5D0,BL)*GAMMA(1.5D0)
      WA=GAMMIN(1.D0,AL)
      WB=GAMMIN(1.D0,BL)
      IF(A.GE.MU)THEN
	 SA=1.D0
      ELSE
	 SA=-1.D0
      ENDIF
      IF(B.GE.MU)THEN
	 SB=1.D0
      ELSE
	 SB=-1.D0
      ENDIF
      MUU=MU + SIG*(SB*WB+SA*WA)/K/SQRT(2.D0*PI)
      CX=SIG**2*(SB*QB-SA*QA)/K/SQRT(PI)
      VAR=CX+2.D0*MU*MUU-MU**2-MUU**2
C      WRITE(*,*) ' VAR,MUU',VAR,MUU
      F=((MUU-XMU)/XVA)**2 + (VAR/XVA/XVA-1.D0)**2
      F=F*1.D6
      RETURN
      END
C
C
      SUBROUTINE TRGRAD(N,X,DF)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DF(N),X(N)
      ON=1.D0
      EPS=1.D-5
      CALL TRFUNC(N,X,F)
      DO1 I=1,N
      XEPS=EPS*(DMAX1(ON,ABS(X(I))))
      XEPSI=ON/XEPS
      X(I)=X(I) + XEPS
      CALL TRFUNC(N,X,FEPS)
      DF(I)=(FEPS - F)*XEPSI
    1 X(I)=X(I) - XEPS
C     .WRITE GRADIENTS - NOT ACTIVE !!
      IF(I.LT.0)THEN
	 WRITE(*,*)
	 WRITE(*,900) 0,F,(DF(I),I=1,N)
      ENDIF
  900 FORMAT(I3,4(1PE12.4))
      RETURN
      END
C
C
      SUBROUTINE PARCAI(IT,P1,P2,P3,P4,S1,S2,S3,S4)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MU,MUU,K
      PI=3.1415926535898D0
      ZE=0.D0
      ON=1.D0
      TW=2.D0
 
      P1=ZE
      P2=ZE
      P3=ZE
      P4=ZE
      IF(IABS(IT).EQ.0)GOTO 5
      IF(IABS(IT).EQ.1)GOTO 10
      IF(IABS(IT).EQ.2)GOTO 20
      IF(IABS(IT).EQ.3)GOTO 30
      IF(IABS(IT).EQ.4)GOTO 40
      IF(IABS(IT).EQ.5)GOTO 50
      IF(IABS(IT).EQ.6)GOTO 60
      IF(IABS(IT).EQ.7)GOTO 70
      IF(IABS(IT).EQ.8)GOTO 80
      IF(IABS(IT).EQ.9)GOTO 90
      IF(IABS(IT).EQ.10)GOTO 100
      IF(IABS(IT).EQ.11)GOTO 110
      WRITE(*,*) ' ERROR IN DISTRIBUTION TYPE'
      STOP
    5 P1=S1
      P2=ZE
      GOTO 7777
   10 IF(IT.GT.0)THEN
	 P1=(S1+S2)/TW
	 P2=(S2-S1)/SQRT(12.D0)
      ELSE
	 P1=S1
	 P2=S2
      ENDIF
      GOTO 7777
   20 P1=S1
      P2=S2
      GOTO 7777
   30 IF(IT.GT.0)THEN
	 P1=EXP(S1+.5D0*S2**2)
	 P2=SQRT(P1**2*(EXP(S2**2)-1.D0))
      ELSE
	 P1=S1
	 P2=S2
      ENDIF
      GOTO 7777
   40 IF(IT.GT.0)THEN
	 P2=ON/S1
	 P1=S2+P2
      ELSE
	 P1=S1
	 P2=S2
      ENDIF
      GOTO 7777
   50 IF(IT.GT.0)THEN
	 P1=S2/S1
	 P2=SQRT(S2)/S1
      ELSE
	 P1=S1
	 P2=S2
      ENDIF
      GOTO 7777
   60 IF(IT.GT.0)THEN
	 P1=S2+.57722D0/S1
	 P2=PI/(SQRT(6.D0)*S1)
      ELSE
	 P1=S1
	 P2=S2
      ENDIF
      GOTO 7777
   70 IF(IT.GT.0)THEN
	 P1=S2-.57722D0/S1
	 P2=PI/(SQRT(6.D0)*S1)
      ELSE
	 P1=S1
	 P2=S2
      ENDIF
      GOTO 7777
   80 IF(IT.GT.0)THEN
	 Q1=GAMMA(ON-ON/S1)
	 Q2=GAMMA(ON-TW/S1)
	 P1=S2*Q1
	 P2=S2*SQRT(Q2-Q1**2)
      ELSE
	 P1=S1
	 P2=S2
      ENDIF
      GOTO 7777
   90 IF(IT.GT.0)THEN
	 Q1=GAMMA(ON+ON/S1)
	 Q2=GAMMA(ON+TW/S1)
	 P1=S2*Q1
	 P2=S2*SQRT(Q2-Q1**2)
      ELSE
	 P1=S1
	 P2=S2
      ENDIF
      GOTO 7777
  100 IF(IT.GT.0)THEN
	 P1=S2+S1*SQRT(PI/TW)
	 P2=S1*SQRT(TW-PI/TW)
      ELSE
	 P1=S1
	 P2=S2
      ENDIF
      GOTO 7777
CXX
  110 IF(IT.GT.0)THEN
	 MU=S1
	 SIG=S2
	 A=S3
	 B=S4
	 AL=((A-MU)/SIG)**2/2.D0
	 BL=((B-MU)/SIG)**2/2.D0
	 UA=(A-MU)/SIG
	 UB=(B-MU)/SIG
	 K=PHI(UB)-PHI(UA)
	 QA=GAMMIN(1.5D0,AL)*GAMMA(1.5D0)
	 QB=GAMMIN(1.5D0,BL)*GAMMA(1.5D0)
	 WA=GAMMIN(1.D0,AL)
	 WB=GAMMIN(1.D0,BL)
	 IF(A.GE.MU)THEN
	  SA=1.D0
	 ELSE
	  SA=-1.D0
	 ENDIF
	 IF(B.GE.MU)THEN
	  SB=1.D0
	 ELSE
	  SB=-1.D0
	 ENDIF
	 MUU=MU + SIG*(SB*WB+SA*WA)/K/SQRT(2.D0*PI)
	 CX=SIG**2*(SB*QB-SA*QA)/K/SQRT(PI)
	 VAR=CX+2.D0*MU*MUU-MU**2-MUU**2
 	 P1=MUU
	 P2=SQRT(VAR)
	 P3=S3
	 P4=S4
      ELSE
	 P1=S1
	 P2=S2
	 P3=S3
	 P4=S4
      ENDIF
 7777 RETURN
      END
C
C
      SUBROUTINE ZUTOX(N,NCOR,U,RL,SPAR,ITYPE,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),U(N),RL(NCOR,NCOR),SPAR(N,4),ITYPE(N)
      DIMENSION PAR(4)
      ZE=0.D0
      DO 20 I=1,N
      Y=U(I)
      IF(I.GT.NCOR)GOTO 11
      S=ZE
      DO 10 J=1,NCOR
      IF(J.GT.I)GOTO 10
      S=S+RL(I,J)*U(J)
   10 CONTINUE
      Y=S
   11 DO 12 J=1,4
   12 PAR(J)=SPAR(I,J)
      IF(ITYPE(I).EQ. 2) THEN
	   X(I)=PAR(2)*Y+PAR(1)
	   GOTO 20
      ENDIF
      IF(ITYPE(I).EQ. 3) THEN
	   X(I)=EXP(PAR(2)*Y+PAR(1))
	   GOTO 20
      ENDIF
      FU=PHI(Y)
      CALL IDISTR(FU,ITYPE(I),PAR,X(I))
   20 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE ZXTOU(N,NCOR,X,RL,SPAR,ITYPE,U)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),U(N),RL(NCOR,NCOR),SPAR(N,4),ITYPE(N)
      DIMENSION PAR(4),XS(200)
C
      DO 40 I=1,N
        Y=X(I)
        XS(I)=X(I)
        DO 22 J=1,4
   22   PAR(J)=SPAR(I,J)
        IF(ITYPE(I).EQ.0) THEN
	   X(I)=0.D0
	   GOTO 40
        ENDIF
        IF(ITYPE(I).EQ. 2) THEN
	   X(I)=(Y-PAR(1))/PAR(2)
	   FX=PHI(U(I))
	   GOTO 40
        ENDIF
        IF(ITYPE(I).EQ. 3) THEN
	   X(I)=(DLOG(Y)-PAR(1))/PAR(2)
	   FX=PHI(U(I))
	   GOTO 40
        ENDIF
        CALL CDISTR(Y,ITYPE(I),PAR,FX)
        X(I)=PHIINV(FX)
   40 CONTINUE
      DO 20 I=1,N
      U(I)=X(I)
      IF(I.GT.NCOR)GOTO 20
      S=X(I)
      DO 10 J=1,NCOR
      IF(J.GE.I)GOTO 10
      S=S-RL(I,J)*U(J)
   10 CONTINUE
      U(I)=S/RL(I,I)
   20 CONTINUE
      DO 104 I=1,N
      X(I)=XS(I)
 104  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE IDISTR(FX,IT,PAR,X)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PAR(4)
*      PI=3.1415926535898D0
      ON=1.D0
      TW=2.D0
      IF(IABS(IT).EQ.0)GOTO 5
      IF(IABS(IT).EQ.1)GOTO 10
      IF(IABS(IT).EQ.2)GOTO 20
      IF(IABS(IT).EQ.3)GOTO 30
      IF(IABS(IT).EQ.4)GOTO 40
      IF(IABS(IT).EQ.5)GOTO 50
      IF(IABS(IT).EQ.6)GOTO 60
      IF(IABS(IT).EQ.7)GOTO 70
      IF(IABS(IT).EQ.8)GOTO 80
      IF(IABS(IT).EQ.9)GOTO 90
      IF(IABS(IT).EQ.10)GOTO 100
      IF(IABS(IT).EQ.11)GOTO 110
      WRITE(*,*) ' ERROR IN DISTRIBUTION TYPE',IT
      STOP 
    5 X=PAR(1)
      GOTO 500
   10 X=PAR(1)+(PAR(2)-PAR(1))*FX
      GOTO 500
   20 X=PAR(1)+PAR(2)*PHIINV(FX)
      GOTO 500
   30 X=EXP(PAR(1)+PAR(2)*PHIINV(FX))
      GOTO 500
   40 X=-DLOG(ON-FX)/PAR(1)+PAR(2)
      GOTO 500
   50 X1=0.01D0
      X2=10.D0*PAR(2)
      XACC=1.D-4
      FACC=1.D-4
      X=RTSAFE(X1,X2,XACC,FACC,IABS(IT),PAR,FX,IER)
      IF(IER.NE.0)THEN
	   STOP ' ERROR IN TRANSFORMATION OF GAMMA VARIABLE'
      ENDIF
      GOTO 500
   60 X=PAR(2)-DLOG(-DLOG(FX))/PAR(1)
      GOTO 500
   70 X=PAR(2)+DLOG(-DLOG(ON-FX))/PAR(1)
      GOTO 500
   80 X=PAR(2) / (-DLOG( FX ))**(ON/PAR(1)) 
      GOTO 500
   90 X=PAR(2) * (-DLOG( ON-FX ))**(ON/PAR(1)) 
      GOTO 500
  100 X=PAR(2)+PAR(1)*SQRT(-TW*DLOG(ON-FX))
      GOTO 500
CXX
  110 SIG=PAR(2)
      XMU=PAR(1)
      A=(PAR(3)-XMU)/SIG
      B=(PAR(4)-XMU)/SIG
      FA=PHI(A)
      FB=PHI(B)
      QQ=FX*(FB-FA)+FA
      X=PHIINV(QQ)*SIG+XMU
C      F1=(X-XMU)/SIG
C      DF1=EXP(-F1*F1/2.D0)/SQRT(2.D0*PI)/SIG/(FB-FA)
      GOTO 500
  500 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CDISTR(X,IT,PAR,FX)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PARAMM/ PAR1,PAR2,PAR3,PAR4
      DIMENSION PAR(4)
      ZE=0.D0
      ON=1.D0
      IF(IABS(IT).EQ.0)GOTO 5
      IF(IABS(IT).EQ.1)GOTO 10
      IF(IABS(IT).EQ.2)GOTO 20
      IF(IABS(IT).EQ.3)GOTO 30
      IF(IABS(IT).EQ.4)GOTO 40
      IF(IABS(IT).EQ.5)GOTO 50
      IF(IABS(IT).EQ.6)GOTO 60
      IF(IABS(IT).EQ.7)GOTO 70
      IF(IABS(IT).EQ.8)GOTO 80
      IF(IABS(IT).EQ.9)GOTO 90
      IF(IABS(IT).EQ.10)GOTO 100
      IF(IABS(IT).EQ.11)GOTO 110
      WRITE(*,*) ' ERROR IN DISTRIBUTION TYPE'
      STOP
    5 FX=0.5D0
      GOTO 500
   10 IF(X.GE.PAR(2))THEN
	   FX=ON
	   GOTO 500
      ENDIF
      IF(X.LE.PAR(1))THEN
	   FX=ZE
	   GOTO 500
      ENDIF
      FX=(X-PAR(1))/(PAR(2)-PAR(1))
      GOTO 500
   20 U=(X-PAR(1))/PAR(2)
      FX=PHI(U)
      GOTO 500
   30 U=(DLOG(X)-PAR(1))/PAR(2)
      FX=PHI(U)
      GOTO 500
   40 FX=ON-EXP(-(X-PAR(2))*PAR(1))
      GOTO 500
   50 P1=PAR(2)
      P2=X*PAR(1)
      FX=GAMMIN(P1,P2)
      GOTO 500
   60 FE=EXP(-PAR(1)*(X-PAR(2)))
      FX=EXP(-FE)
      GOTO 500
   70 FE=EXP(PAR(1)*(X-PAR(2)))
      FX=ON-EXP(-FE)
      GOTO 500
   80 FX=EXP( - ((PAR(2)/X))**PAR(1) )
      GOTO 500
   90 FX=ON-EXP( - (X/PAR(2))**PAR(1) )
      GOTO 500
  100 FX=ON-EXP( - 0.5D0 * ((X-PAR(2))/PAR(1))**2 )
      GOTO 500
CXX
  110 SIG=PAR(2)
      XMU=PAR(1)
      A=(PAR(3)-XMU)/SIG
      B=(PAR(4)-XMU)/SIG
      XX=(X-XMU)/SIG
      FA=PHI(A)
      FB=PHI(B)
      FX=PHI(XX)
      IF(X.LE.PAR(3))THEN
	   FX=0.D0
      ELSE IF(X.GE.PAR(4))THEN
	   FX=1.D0
      ELSE
	   FX=(FX-FA)/(FB-FA) 
      ENDIF
  500 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CHOLES(N,NN,R,ISW)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(NN,NN)
      ZE=0.D0
      ON=1.D0
      ISW=0
      DO 10 I=1,N
      DO 10 J=I,N
   10 R(J,I)=ZE
      DO 40 K=1,N
      print *, 'K',K
      K1=K-1
      print *, 'K1',K1
      DO 40 I=K,N
      print *, 'I',I
      C=R(K,I)
      print *, 'C',C
      IF(I.EQ.K)C=ON
      IF(K1.EQ.0) GOTO 35
      DO 30 J=1,K1
   30 C=C-R(I,J)*R(K,J)
      
   35 IF(I.EQ.K)GOTO 38
      R(I,K)=C/BKK
      GOTO 40
   38 IF(C.GT.ZE)GOTO 39
      ISW=1
      RETURN
   39 BKK=SQRT(C)
      R(K,K)=BKK
   40 CONTINUE
       print *, ' R(2,2)', R(2,2)
       print *,'R(J,I)',((R(J,I),I=1,N),J=1,N)
      RETURN
      END
C
C
      REAL*8 FUNCTION GAMMIN(A,X)
      IMPLICIT REAL*8(A-H,O-Z)
      IF(X.LT.0..OR.A.LE.0.)THEN
	   WRITE(*,*) ' ERROR IN INCOMPLETE GAMMA FUNCTION'
	   WRITE(*,*) 'A = ',A,' X = ',X
	   STOP
      ENDIF
      IF(X.LT.A+1.)THEN
	   CALL GSER(GAMMIN,A,X,GLN)
      ELSE
	   CALL GCF(GAMMCF,A,X,GLN)
	   GAMMIN=1.-GAMMCF
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ITMAX=200,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
	   GAMSER=0.D0
	   RETURN
      ENDIF
      AP=A
      SUM=1.D0/A
      DEL=SUM
      DO 11 N=1,ITMAX
	   AP=AP+1.D0
	   DEL=DEL*X/AP
	   SUM=SUM+DEL
	   IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      WRITE(*,*) 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*DLOG(X)-GLN)
      RETURN
      END
C
C
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ITMAX=300,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.D0
      A0=1.D0
      A1=X
      B0=0.D0
      B1=1.D0
      FAC=1.D0
      DO 11 N=1,ITMAX
	   AN=FLOAT(N)
	   ANA=AN-A
	   A0=(A1+A0*ANA)*FAC
	   B0=(B1+B0*ANA)*FAC
	   ANF=AN*FAC
	   A1=X*A0+ANF*A1
	   B1=X*B0+ANF*B1
	   IF(A1.NE.0.D0)THEN
	    FAC=1.D0/A1
	    G=B1*FAC
	    IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
	    GOLD=G
	   ENDIF
11    CONTINUE
      WRITE(*,*) ' A,X',A,X
      WRITE(*,*) 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*DLOG(X)-GLN)*G
      RETURN
      END
C
C
      REAL*8 FUNCTION GAMMA(XX)
      IMPLICIT REAL*8(A-H,O-Z)
      GAMMA=EXP(GAMMLN(XX))
      RETURN
      END
C
C
      REAL*8 FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,XX,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
	   X=X+ONE
	   SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C
C
      REAL*8 FUNCTION RTSAFE(X1,X2,XACC,FACC,ITYPE,PAR,FU,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PAR(4)
      PARAMETER (MAXIT=100)
      IER=0
      CALL DIROOT(X1,FL,DF,ITYPE,PAR,FU,XACC)
      CALL DIROOT(X2,FH,DF,ITYPE,PAR,FU,XACC)
      IF(FL*FH.GE.0.) THEN
	   WRITE(*,*) 'root must be bracketed'
	   IER=1
	   RETURN
      ENDIF
      IF(FL.LT.0.)THEN
	   XL=X1
	   XH=X2
      ELSE
	   XH=X1
	   XL=X2
	   SWAP=FL
	   FL=FH
	   FH=SWAP
      ENDIF
      RTSAFE=.5*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL DIROOT(RTSAFE,F,DF,ITYPE,PAR,FU,XACC)
      DO 11 J=1,MAXIT
	   IF(((RTSAFE-XH)*DF-F)*((RTSAFE-XL)*DF-F).GE.0.
     *      .OR. ABS(2.*F).GT.ABS(DXOLD*DF) ) THEN
	    DXOLD=DX
	    DX=0.5*(XH-XL)
	    RTSAFE=XL+DX
	    IF(XL.EQ.RTSAFE)RETURN
	   ELSE
	    DXOLD=DX
	    DX=F/DF
	    TEMP=RTSAFE
	    RTSAFE=RTSAFE-DX
	    IF(TEMP.EQ.RTSAFE)RETURN
	   ENDIF
	   IF(ABS(DX).LT.XACC .AND. F.LT.FACC) RETURN
	   CALL DIROOT(RTSAFE,F,DF,ITYPE,PAR,FU,XACC)
	   IF(F.LT.0.) THEN
	    XL=RTSAFE
	    FL=F
	   ELSE
	    XH=RTSAFE
	    FH=F
	   ENDIF
11    CONTINUE
      WRITE(*,*) 'RTSAFE exceeding maximum iterations'
      IER=2
      RETURN
      END
C
C
      SUBROUTINE DIROOT(X,F,DF,ITYPE,PAR,FU,XACC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION PAR(4)
      CALL CDISTR(X,ITYPE,PAR,FX)
      F=FX-FU
      ACC=XACC*1.D-3
      XX=X+ACC
      CALL CDISTR(XX,ITYPE,PAR,FX)
      FX=FX-FU
      DF=(FX-F)/ACC
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION PHI (X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST/ EPSMIN,EPSMAX,EMAX,SQ2PI
      DIMENSION B(6)
      B(1)= 0.319381530D0
      B(2)=-0.356563782D0
      B(3)= 1.781477937D0
      B(4)=-1.821255978D0
      B(5)= 1.330274429D0
      B(6)= 0.2316419D0
      IF (X.LT.0.D0) THEN
	 ISW = -1
      ELSE
	 ISW = 1
      ENDIF
      XX=ABS(X)
      IF (XX.GT.EMAX) THEN
	 IF (X.GT.0.D0) THEN
	  PHI = 1.D0
	 ELSE
	  PHI = 0.D0
	 ENDIF
	 GOTO 999
      ELSE
	 IF (XX.LT.EPSMIN) XX = 0.D0
      ENDIF
      IF (XX.EQ.0.D0) THEN
	 PHI = 0.5D0
	 GOTO 999
      ENDIF
      IF (XX.LE.4.D0) THEN
	 T = 1.D0/(1.D0+B(6)*XX)
	 PHI = EXP(-XX*XX/2.D0)*SQ2PI*
     -     ((((T*B(5) + B(4))*T + B(3))*T + B(2))*T + B(1))*T
      ELSE
	 PHI = (SQRT(4.D0 + XX*XX) - XX)/2.D0*SQ2PI*
     -     EXP(-XX*XX/2.D0)
      ENDIF
      IF (ISW.EQ.1) PHI = 1.D0 - PHI
  999 RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION PHIINV (P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST/ EPSMIN,EPSMAX,EMAX,SQ2PI
      DIMENSION C(3),D(3)
      C(1)= 2.515517D0
      C(2)= 0.802853D0
      C(3)= 0.010328D0
      D(1)= 1.432788D0
      D(2)= 0.189269D0
      D(3)= 0.001308D0
      IF (P.EQ.0.D0) THEN
	 PHIINV = -EMAX
	 GOTO 999
      ENDIF
      IF (P.EQ.1.D0) THEN
	 PHIINV = EMAX
	 GOTO 999
      ENDIF
      IF (P.GT.0.5D0) THEN
	 ISW = 1
	 PP = 1.D0 - P
      ELSE
	 ISW = -1
	 PP = P
      ENDIF
      IF (PP.LT.SQRT(EPSMIN)) PP = SQRT(EPSMIN)
      IF (PP.GT.1.D-4) THEN
	 T = SQRT(DLOG(1.D0/(PP*PP)))
	 PHIINV = T - ((T*C(3) + C(2))*T + C(1))/
     -      (((T*D(3) + D(2))*T + D(1))*T + 1.D0)
      ELSE
	 X = 3.D0
	 DO 10 J = 1,4
	  X = SQRT(2.D0*(DLOG((SQRT(4.D0 + X*X) - X)/2.D0*
     -      SQ2PI) - DLOG(PP)))
   10  CONTINUE
	 PHIINV = X
      ENDIF
      PHIINV = ISW*PHIINV
  999 RETURN
      END
C
C
      SUBROUTINE INPREA(MINP,CARR,NC,IARR,NI,RARR,NR,QARR,IER)
c   MINP为读文件名, CARR大概是字符串数组,NC为字符串个数,IARR为整数数组,NI为整数个数,
c   RARR为实数数组,NR为实数个数,QARR不知道,IER指示读入数据过程中是否有错.
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 CHAR
      CHARACTER*8 CARR
      CHARACTER*60 QARR
      CHARACTER*20 CHAR2
      DIMENSION IARR(10),RARR(5),carr(20)
      NR=0
      NI=0
      NC=0
      IER=0
  900 FORMAT(A80)
   20 READ(MINP,900,END=8888) CHAR
      print *, 'CHAR',CHAR
      IF(CHAR(1:1).EQ.'*')THEN
	 GOTO 20
c     ELSE IF(CHAR(1:1).EQ.'$')THEN
c	 NC=1
c	 CARR(1)=CHAR(2:9)
      ELSE
	 I0=1
	 DO 40 I=1,78
	 IF(CHAR(I:I).EQ.',' .OR. I.EQ.78)THEN
	  print *, 'I',I
	   print *, 'I0',I0
	  ISW=0
	  DO 25 K=I0,(I-1)
   25   IF(CHAR(K:K).EQ.'.')ISW=1
	  CHAR2=CHAR(I0:(I-1))
	   print *, 'CHAR2',CHAR2
	    print *, 'I0',I0
	  IF(ISW.EQ.1)THEN  
	   READ(CHAR2,800,END=30,ERR=30) REA
  800    FORMAT(E20.10)
         print *, 'REA',REA
	   NR=NR+1
	   print *, 'NR',NR
	   RARR(NR)=REA
	   GOTO 35
	  ELSE 
	   READ(CHAR2,810,END=30,ERR=30) INT
  810    FORMAT(I20)
         print *, 'INT',INT
	   NI=NI+1
	   print *, 'NI',NI
	   IARR(NI)=INT
	   print *, 'IARR(NI)',IARR(NI)
	   GOTO 35
	  ENDIF
   30   IF(I.LT.78 .OR. I0.EQ.1)THEN
	   NC=NC+1
	   CARR(NC)=CHAR(I0:(I-1))//'       '
	    print *, 'CARR(NC)',CARR(NC)
	  ELSE
	   QARR=CHAR(I0:78)
	   print *, 'QARR',QARR
	  ENDIF
   35   I0=I+1
	 ENDIF
   40  CONTINUE
      ENDIF
      GOTO 7777
 8888 IER=7777
 7777 RETURN
      END
C
C
      SUBROUTINE UPTRAN (NAME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 NAME
      DO 20 I = 1,8
      JCHAR = ICHAR(NAME(I:I))
       print *, 'I',I
        print *, 'NAME(I:I)',NAME(I:I)
       print *, 'JCHAR',JCHAR
      IF ((JCHAR.EQ.32).OR.(JCHAR.EQ.254)) THEN
	 DO 10 J = I,8
   10   NAME(J:J) = ' '
	 GOTO 30
      ELSE
	 IF (JCHAR.GE.97) NAME(I:I) = CHAR(JCHAR-32)
      ENDIF
      print *, ' CHAR(JCHAR-32)', CHAR(JCHAR-32)
   20 CONTINUE
   30 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE UPDATX (N,NP,STONAM,PARNAM,X,P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CHANAM/ NN,NNP,CVARNA(200),XVARNA(200)
      CHARACTER*8 STONAM(N),PARNAM(NP),CVARNA
      DIMENSION X(N),P(NP)
      DO 30 I = 1,N
      XVARNA(I) = X(I)
   30 CVARNA(I) = STONAM(I)
      IF(NP.GT.0)THEN
	 DO 40 I = 1,NP
	  XVARNA(N+I) = P(I)
   40   CVARNA(N+I) = PARNAM(I)
      ENDIF
      RETURN
      END
C
C
      SUBROUTINE QDATE(MOUT)
C     .     VAX
c      CHARACTER*9 DATO
c      CHARACTER*8 TID
c      CALL DATE(DATO)
c      CALL TIME(TID)
c      WRITE(MOUT,880) DATO,TID
c  880 FORMAT(/,' DATE : ',A9,' TIME : ',A8)
C
C     .     POWER FORTRAN
      INTEGER*2 IYR,IMON,IDAY,IHR,IMIN,ISEC,I100TH
      CALL GETDAT(IYR,IMON,IDAY)
      CALL GETTIM(IHR,IMIN,ISEC,I100TH)
      WRITE(MOUT,900) IDAY,IMON,IYR,IHR,IMIN,ISEC
  900 FORMAT(/,' DATE : ',I2,'-',I2.2,'-',I4.4,
     -         ' TIME : ',I2,':',I2.2,':',I2.2)
C
      RETURN
      END
C
C
      SUBROUTINE JOBNAM(IBNAME,JNAME)
      CHARACTER*60 CHAMOD
      CHARACTER*15 JNAME
      JNAME='        '
C     .     LAHEY F77L:
*      CALL GETCL(CHAMOD)
*      IBNAME=INDEX(CHAMOD(2:60),' ')
*      JNAME=CHAMOD(1:ibname)
c
C     .     POWER FORTRAN:
!      CALL GETARG(1,CHAMOD,ier)
      CALL GET_COMMAND(CHAMOD) 
!      CALL GETARG(0,CHAMOD,ier)
      IBNAME=INDEX(CHAMOD(1:60),' ')-1
      JNAME=CHAMOD(1:ibname)
      RETURN
      END

C
C
      SUBROUTINE dfpmin(p,n,gtol,FP,ITER,ier)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(N)
      PARAMETER (NMAX=50,ITMAX=200,STPMX=100.,EPS=3.e-8,TOLX=4.*EPS)
      LOGICAL check
      DIMENSION dg(NMAX),g(NMAX),hdg(NMAX),hessin(NMAX,NMAX),
     -          pnew(NMAX),xi(NMAX)
      ier=0
      CALL trfunc(N,p,FP)
      call TRGRAD(N,p,g)
      sum=0.D0
      do 12 i=1,n
	   do 11 j=1,n
	    hessin(i,j)=0.D0
11     continue
	   hessin(i,i)=1.D0
	   xi(i)=-g(i)
	   sum=sum+p(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum),float(n))
      do 27 its=1,ITMAX
	   iter=its
	   call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check)
	   fp=fret
	   do 13 i=1,n
	    xi(i)=pnew(i)-p(i)
	    p(i)=pnew(i)
13     continue
	   test=0.D0
	   do 14 i=1,n
	    temp=abs(xi(i))/max(abs(p(i)),1.D0)
	    if(temp.gt.test)test=temp
14     continue
	   if(test.lt.TOLX)return
	   do 15 i=1,n
	    dg(i)=g(i)
15     continue
	   call TRGRAD(N,p,g)
	   test=0.d0
	   den=max(fret,1.D0)
	   do 16 i=1,n
	    temp=abs(g(i))*max(abs(p(i)),1.D0)/den
	    if(temp.gt.test)test=temp
16     continue
	   if(test.lt.gtol)return
	   do 17 i=1,n
	    dg(i)=g(i)-dg(i)
17     continue
	   do 19 i=1,n
	    hdg(i)=0.D0
	    do 18 j=1,n
	     hdg(i)=hdg(i)+hessin(i,j)*dg(j)
18      continue
19     continue
	   fac=0.D0
	   fae=0.D0
	   sumdg=0.D0
	   sumxi=0.D0
	   do 21 i=1,n
	    fac=fac+dg(i)*xi(i)
	    fae=fae+dg(i)*hdg(i)
	    sumdg=sumdg+dg(i)**2
	    sumxi=sumxi+xi(i)**2
21     continue
	   if(fac**2.gt.EPS*sumdg*sumxi)then
	    fac=1.D0/fac
	    fad=1.D0/fae
	    do 22 i=1,n
 	     dg(i)=fac*xi(i)-fad*hdg(i)
22      continue
	    do 24 i=1,n
	     do 23 j=1,n
	      hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+
     *      fae*dg(i)*dg(j)
23       continue
24      continue
	   endif
	   do 26 i=1,n
	    xi(i)=0.D0
	    do 25 j=1,n
	     xi(i)=xi(i)-hessin(i,j)*g(j)
25      continue
26     continue
27    continue
      ier=1
      return
      END
C
C
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL check
      DIMENSION g(n),p(n),x(n),xold(n)
      PARAMETER (ALF=1.e-4,TOLX=1.e-7)
      check=.false.
      sum=0.D0
      do 11 i=1,n
	   sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
	   do 12 i=1,n
	    p(i)=p(i)*stpmax/sum
12     continue
      endif
      slope=0.D0
      do 13 i=1,n
	   slope=slope+g(i)*p(i)
13    continue
      test=0.D0
      do 14 i=1,n
	   temp=abs(p(i))/max(abs(xold(i)),1.D0)
	   if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.D0
1     continue
	  do 15 i=1,n
	   x(i)=xold(i)+alam*p(i)
15    continue
	  CALL TRfunc(N,x,F)
	  if(alam.lt.alamin)then
	   do 16 i=1,n
	    x(i)=xold(i)
16     continue
	   check=.true.
	   return
	  else if(f.le.fold+ALF*alam*slope)then
	   return
	  else
	   if(alam.eq.1.D0)then
	    tmplam=-slope/(2.D0*(f-fold-slope))
	  else
	    rhs1=f-fold-alam*slope
	    rhs2=f2-fold2-alam2*slope
	    a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
	    b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
	    if(a.eq.0.D0)then
	      tmplam=-slope/(2.D0*b)
	    else
	      disc=b*b-3.D0*a*slope
	      tmplam=(-b+sqrt(disc))/(3.D0*a)
	    endif
	    if(tmplam.gt..5D0*alam)tmplam=.5D0*alam
	   endif
	  endif
	  alam2=alam
	  f2=f
	  fold2=fold
	  alam=max(tmplam,.1D0*alam)
      goto 1
      END
