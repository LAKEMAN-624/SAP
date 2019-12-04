      SUBROUTINE INPUT01(nZYD, nJD, nDY, nExc, nWYYS, IER)
C     INPUT OF STRUCTURES ANALYSIS
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 CARR(20)
      CHARACTER*60 QARR
      DIMENSION IARR(10),RARR(5)
      INCLUDE 'COMMY.INC'

      IER=0
C
	WRITE(MRE3,2888) 
      WRITE(MRE3,901)
  901 FORMAT(' Run with structures analysis of plain truss, ')
	WRITE(MRE3,2888) 
      WRITE(MRE3,902)
  902 FORMAT(/,' Structures general information:')
      CALL INPREA(MFEM,CARR,NC,IARR,NI,RARR,NRR,QARR,IER)
      IF(NC.NE.0 .OR. NI.NE.5 .OR. NRR.NE.0)THEN
	WRITE(*,*) ' ERROR IN INPUT OF Structures information'
	IER=111
	RETURN
      ENDIF
      nJD=IARR(1)
      nZYD=IARR(2)
      nDY=IARR(3)
      nExc=IARR(4)
      nWYYS=IARR(5)

      WRITE(MRE3,903)nJD,nZYD,nDY,nExc,nWYYS
  903 FORMAT(' nJD=',I4,'  nZYD=',I4,'  nDY=',I4,'  nExc=',I4,
     &'  nWYYS=',I4/)
      IF (nJD.GT.80.OR.nDY.GT.100.) THEN
	WRITE(*,*) 'Too many nodes and elements!'
	IER=111
      ENDIF

 2888 FORMAT(' ***************************************************',
     - '*******************')
      RETURN
      END
c
c
      SUBROUTINE INPUT1(NVAR,NPAR,NCOR,MAXIT,IER)
C**********************************************************************
C
C     .INPUT OF INFORMATION ON STOCHASTIC MODEL.
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 CARR(20)
      CHARACTER*60 QARR
      DIMENSION IARR(10),RARR(5)
      INCLUDE 'COMMY.INC'
      COMMON /CONST/ EPSMIN,EPSMAX,EMAX,SQ2PI
      COMMON /YCMACH/ CEMIN,CEMAX,COGEN,COMPI,ZERO,ONE
      EPSMIN = 0.538D-292
      EPSMAX = 0.692D+293
      EMAX  = SQRT(-2.D0*LOG(EPSMIN))
      PI    = 3.14159265358979323846D0
      SQ2PI = 1/SQRT(2.D0*PI)
      CEMIN = .538D-292
      CEMAX = .692D+293
      COGEN = 1.D-15
      ZERO  = 0.D0
      ONE   = 1.D0
      COMPI = DACOS(-ONE)
C
      IER=0

	WRITE(MRE2,2888) 
      WRITE(MRE2,901)
  901 FORMAT(' run with RELiability part of PRADSS. J.D. Sorensen, ',
     - 'Aalborg University')
	WRITE(MRE2,2888) 
      WRITE(MRE2,902)
  902 FORMAT(/,' General information:')
      CALL INPREA(MREL,CARR,NC,IARR,NI,RARR,NR,QARR,IER)
      IF(NC.NE.0 .OR. NI.NE.4 .OR. NR.NE.0)THEN
	WRITE(*,*) ' ERROR IN INPUT : CONTROL DATA'
	STOP
      ENDIF
      NVAR=IARR(1)
      NPAR=IARR(2)
      NCOR=IARR(3)
      MAXIT=IARR(4)

      WRITE(MRE2,903)NVAR,NPAR,NCOR,MAXIT
  903 FORMAT(' NVAR=',I4,'  NPAR=',I4,'  NCOR=',I4,'  MAXIT=',I5/)
  !2018/10/09 MAXIT=',I4改成MAXIT=',I5
      IF (NVAR+NPAR.GT.200) THEN
	WRITE(*,*) 'Too many variables and parameters!'
	IER=111
      ENDIF

 2888 FORMAT(' ***************************************************',
     - '*******************')
      RETURN
      END
C
C
      SUBROUTINE INPUT2(NVAR,NPAR,NCOR,CORREL,ITYPE,
     -                SPAR,P,PL,PU,STONAM,PARNAM,IER)
C**********************************************************************
C
C     .MAIN SUBROUTINE FOR INPUT OF STATISTICAL PARAMETERS
C
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 STONAM(NVAR),PARNAM(NPAR),NAME
      CHARACTER*8 CARR(20)
      CHARACTER*60 QARR
      CHARACTER*78 TEXT
      DIMENSION IARR(10),RARR(5)
      DIMENSION CORREL(NCOR,NCOR),ITYPE(NVAR),SPAR(NVAR,4)
      DIMENSION P(NPAR),PL(NPAR),PU(NPAR)
      INCLUDE 'COMMY.INC'
C
      IER=0
C
C     .READ STOCHASTIC VARIABLES
C
      WRITE(MRE2,900)
  900 FORMAT(' Stochastic variables:')
      DO 10 I = 1,NVAR
      CALL INPREA(MREL,CARR,NC,IARR,NI,RARR,NR,QARR,IER)
       print *, 'NC,NI,NR',NC,NI,NR
      IF(NC.NE.1 .OR. NI.NE.2 .OR. NR.NE.4)THEN
      
	WRITE(*,*) ' ERROR IN INPUT2 OF STOCHASTIC VARIABLES'
	IER=111
	RETURN
      ENDIF
      IVAR=IARR(1)
!      print *, 'IVAR',IVAR
      IF(IVAR.NE.I)THEN
	WRITE(*,*) ' WRONG ORDER OF STOCHASTIC VARIABLES IN INPUT'
	IER=111
	RETURN
      ENDIF
      NAME=CARR(1)
!      print *, 'NAME',NAME
      CALL UPTRAN (NAME)
      STONAM(I) = NAME
!      print *, ' I', I
!      print *, ' STONAM(I)', STONAM(I)
      JTYPE=IARR(2)
      ITYPE(I) = JTYPE
      WRITE(MRE2,910)I,CARR(1)
  910 FORMAT(' NO: ',I4,'   ',A8)
C
      SPAR1=RARR(1)
      SPAR2=RARR(2)
      SPAR3=RARR(3)
      SPAR4=RARR(4)
      print *, ' ------------------------'
      print *, ' input2'
      print *, ' ------------------------'
      CALL PARCAL(MRE2,ITYPE(I),SPAR1,SPAR2,SPAR3,SPAR4,
     -             SPAR(I,1),SPAR(I,2),SPAR(I,3),SPAR(I,4),1)
   10 CONTINUE
C
C     .READ CORRELATION COEFFICIENTS
C
      IF(NCOR.GT.0)THEN 
      DO 25 I=1,NCOR
      DO 25 J=1,NCOR
      IF(I.NE.J)THEN
	CORREL(I,J)=0.D0
      ELSE
	CORREL(I,I)=1.D0
      ENDIF
   25 CONTINUE
      CALL INPREA(MREL,CARR,NC,IARR,NI,RARR,NRR,QARR,IER)
      IF(NC.NE.0 .OR. NI.NE.1 .OR. NRR.NE.0)THEN
	WRITE(*,*) ' ERROR IN INPUT2 : CORREL'
	IER=111
	RETURN
      ENDIF
      NCORI=IARR(1)
	DO 30 I=1,NCORI
	CALL INPREA(MREL,CARR,NC,IARR,NI,RARR,NR,QARR,IER)
	IF(NC.NE.0 .OR. NI.NE.2 .OR. NR.NE.1)THEN
	  WRITE(*,*) ' ERROR IN INPUT2 : CORREL'
	  IER=111
	  RETURN
	ENDIF
	I1=IARR(1)
	 print *, ' I1', I1
	I2=IARR(2)
	 print *, ' I2', I2
	IF(I1.GT.NCOR .OR. I2.GT.NCOR)THEN
	  WRITE(*,*) ' ERROR IN INPUT2 OF CORR.COEFF. MAX NUMBER OF',
     -               ' CORRELATED VARIABLES:',NCOR
	  WRITE(*,*) ' VARIABEL 1, VARIABEL 2:',I1,I2
	  STOP
	ENDIF
	CORREL(I1,I2)= RARR(1)
	CORREL(I2,I1)= RARR(1)
	 print *, ' CORREL(I1,I2)', CORREL(I1,I2)
   30   CONTINUE
C
	WRITE(MRE2,1032)
	DO 35 I=1,NCOR
   35   WRITE(MRE2,922) (CORREL(I,J),J=1,NCOR)
 1032   FORMAT(/,' Correlation coefficients:')
  922   FORMAT(10F7.3)
	CALL CHOLES(NCOR,NCOR,CORREL,IER)
C加了输出CHOLES处理后的相关系数矩阵
	WRITE(MRE2,1033)
	DO 36 I=1,NCOR
   36   WRITE(MRE2,923) (CORREL(I,J),J=1,NCOR)
 1033   FORMAT(/,' Correlation coefficients1111:')
  923   FORMAT(10F7.3)
	IF(IER.NE.0)THEN
	  WRITE(*,*) ' CORR.COEFF. MATRIX NOT POS. DEFINITE'
	  STOP
	ENDIF
      ENDIF
C
C     .READ PARAMETERS
C
      WRITE(MRE2,930)
  930 FORMAT(//,' Deterministic parameters:')
      DO 20 I = 1,NPAR
      CALL INPREA(MREL,CARR,NC,IARR,NI,RARR,NR,QARR,IER)
      IF(NC.NE.1 .OR. NI.NE.1 .OR. NR.NE.3)THEN
	  WRITE(*,*) ' ERROR IN INPUT2 OF Deterministic PARAMETS'
	  IER=111
	  RETURN
      ENDIF
      IVAR=IARR(1)
      IF(IVAR.NE.I)THEN
	WRITE(*,*) ' WRONG ORDER OF STOCHASTIC VARIABLES IN INPUT'
	IER=111
	RETURN
      ENDIF
      NAME=CARR(1)
      CALL UPTRAN (NAME)
      PARNAM(I)=NAME
      P(I)=RARR(1)
      PL(I)=RARR(2)
      PU(I)=RARR(3)
      WRITE(MRE2,940)I,CARR(1),P(I),PL(I),PU(I)
  940 FORMAT(' NO: ',I4,'   ',A8,'(P, PL, PU): ',1PE12.4,'  ',
     -     1PE12.4,'  ',1PE12.4)
   20 CONTINUE
C
	WRITE(MRE1,2999) 
       CALL QDATE(MRE1)
	WRITE(MRE1,2888) 
      WRITE(*,2999)
      WRITE(*,2888)
C
      ISQ=1
      REWIND (MREL)
      DO WHILE (ISQ.EQ.1)
	 READ(MREL,888,END=55) TEXT
	 WRITE(MRE1,888) TEXT
      ENDDO
  888 FORMAT(A78)
   55 CONTINUE
C
 2999 FORMAT(' run with RELiability part of PRADSS. J.D. Sorensen, ',
     - 'Aalborg University')
 2888 FORMAT(' ***************************************************',
     - '*******************',/)
      WRITE(MRE2,1053)
 1053 FORMAT(//,' Reliability index iterations:',/)
      RETURN
      END
C
C
      SUBROUTINE INPUT02(nZYD, nJD, nDY, nExc, nWYYS, X, Y, Ela,Ruo,  
     -        Sigc, Sigt,Area,Exc,WYYS, RuoL,II,JJ,IP,JP,IExc,IWYYS,IER)
C**********************************************************************
C
C     .MAIN SUBROUTINE FOR INPUT OF STRUCTURES ANALYSIS
C
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 CARR(20)
      CHARACTER*60 QARR
      DIMENSION IARR(10),RARR(5)
      INCLUDE 'COMMY.INC'
      DIMENSION X(nJD),Y(nJD),Ela(nDY),Ruo(nDY),Sigc(nDY),Sigt(nDY)
      DIMENSION Area(nDY),RuoL(nDY),Exc(nExc),WYYS(nWYYS)
      DIMENSION II(nDY),JJ(nDY),IP(nJD),JP(nJD),IExc(nExc),IWYYS(nWYYS)
C
C     .READ Nodal Information
C
      WRITE(MRE3,910)
  910 FORMAT(' Nodal Information:')
      DO 10 I = 1,nJD
      CALL INPREA(MFEM,CARR,NC,IARR,NI,RARR,NRR,QARR,IER)
      IF(NC.NE.0 .OR. NI.NE.3 .OR. NRR.NE.2)THEN
	WRITE(*,*) ' ERROR IN INPUT : PARAMET'
	IER=111
	RETURN
      ENDIF
      IVAR=IARR(1)
      IF(IVAR.NE.I)THEN
	WRITE(*,*) ' WRONG ORDER OF STOCHASTIC VARIABLES IN INPUT'
	IER=111
	RETURN
      ENDIF
	X(I)=RARR(1)
	Y(I)=RARR(2)
	IP(I)=IARR(2)
	JP(I)=IARR(3)
      WRITE(MRE3,920)I,X(I),Y(I),IP(I),JP(I)
  920 FORMAT(' NO:',I3,'  ',1PE12.4,'  ',1PE12.4,'  ',I4,'  ',I4)
   10 CONTINUE
C
C     .READ Element Information
C
      WRITE(MRE3,930)
  930 FORMAT(/,' Element Information:')
      DO 20 I = 1,nDY
      CALL INPREA(MFEM,CARR,NC,IARR,NI,RARR,NRR,QARR,IER)
      IF(NC.NE.0 .OR. NI.NE.3 .OR. NRR.NE.4)THEN
	WRITE(*,*) ' ERROR IN INPUT : PARAMET'
	IER=111
	RETURN
      ENDIF
      IVAR=IARR(1)
      IF(IVAR.NE.I)THEN
	WRITE(*,*) ' WRONG ORDER OF STOCHASTIC VARIABLES IN INPUT'
	IER=111
	RETURN
      ENDIF
	II(I)=IARR(2)
	JJ(I)=IARR(3)
	Ela(I)=RARR(1)
	Ruo(I)=RARR(2)
	Sigc(I)=RARR(3)
	Sigt(I)=RARR(4)
      WRITE(MRE3,940)I,II(I),JJ(I),Ela(I),Ruo(I),Sigc(I),Sigt(I)
  940 FORMAT(' NO:',I3,'  ',I4,'  ',I4,'  ',1PE12.4,'  ',
     -        1PE12.4,'  ',1PE12.4,'  ',1PE12.4)
   20 CONTINUE
C
C     .READ Nodal Excitations Information
C
      WRITE(MRE3,950)
  950 FORMAT(/,' Nodal Excitations Information:')
      DO 30 I = 1,nExc
      CALL INPREA(MFEM,CARR,NC,IARR,NI,RARR,NRR,QARR,IER)
      IF(NC.NE.0 .OR. NI.NE.2 .OR. NRR.NE.1)THEN
	WRITE(*,*) ' ERROR IN INPUT : PARAMET'
	IER=111
	RETURN
      ENDIF
      IVAR=IARR(1)
      IF(IVAR.NE.I)THEN
	WRITE(*,*) ' WRONG ORDER OF STOCHASTIC VARIABLES IN INPUT'
	IER=111
	RETURN
      ENDIF
	Exc(I)=RARR(1)
	IExc(I)=IARR(2)
      WRITE(MRE3,960)I,Exc(I),IExc(I)
  960 FORMAT(' NO:'I3,'  ',1PE12.4,'  ',I4)
   30 CONTINUE
C
C     .READ Section Areas
C
      WRITE(MRE3,970)
  970 FORMAT(/,' Section Areas:')
      DO 40 I = 1,nDY
      CALL INPREA(MFEM,CARR,NC,IARR,NI,RARR,NRR,QARR,IER)
      IF(NC.NE.0 .OR. NI.NE.1 .OR. NRR.NE.1)THEN
	WRITE(*,*) ' ERROR IN INPUT : PARAMET'
	IER=111
	RETURN
      ENDIF
      IVAR=IARR(1)
      IF(IVAR.NE.I)THEN
	WRITE(*,*) ' WRONG ORDER OF STOCHASTIC VARIABLES IN INPUT'
	IER=111
	RETURN
      ENDIF
	Area(I)=RARR(1)
      WRITE(MRE3,980)I,Area(I)
  980 FORMAT(' NO:',I3,'  ',1PE12.4)
   40 CONTINUE
C
C     .READ Displacement Constrains
C
      WRITE(MRE3,990)
  990 FORMAT(/,' Displacement Constrains:')
      DO 50 I = 1,nWYYS
      CALL INPREA(MFEM,CARR,NC,IARR,NI,RARR,NRR,QARR,IER)
      IF(NC.NE.0 .OR. NI.NE.2 .OR. NRR.NE.1)THEN
	WRITE(*,*) ' ERROR IN INPUT : PARAMET'
	IER=111
	RETURN
      ENDIF
      IVAR=IARR(1)
      IF(IVAR.NE.I)THEN
	WRITE(*,*) ' WRONG ORDER OF STOCHASTIC VARIABLES IN INPUT'
	IER=111
	RETURN
      ENDIF
	WYYS(I)=RARR(1)
      IWYYS(I)=IARR(2)
      WRITE(MRE3,1000)I,WYYS(I),IWYYS(I)
 1000 FORMAT(' NO:',I3,'  ',1PE12.4,'  ',I4)
   50 CONTINUE
C
      DO 60 I = 1,nDY
      i1=II(I)
      i2=JJ(I)
      x12=X(i2)-X(i1)
      y12=Y(i2)-Y(i1)
      RuoL(I)=sqrt(x12*x12+y12*y12)*Ruo(I)
   60 CONTINUE

      WRITE(MRE3,1053)
 1053 FORMAT(//,' Structural analysis results:')
      RETURN
      END
