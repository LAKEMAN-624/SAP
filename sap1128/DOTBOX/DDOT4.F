
      SUBROUTINE DOT203(NDV,NCON,NGMAX,NGT,G,IC,CT,CTMIN,NAC,NVC,IGMAX)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(*),IC(*)
      NAC=0
      NVC=0
      NGT=0
      DO 60 I=1,NCON
         G1=G(I)
         IF (G1.LT.CT) GO TO 60
         IF (G1.GE.CT.AND.G1.LE.CTMIN) NAC=NAC+1
         IF (G1.GT.CTMIN) NVC=NVC+1
         NGT=NGT+1
         IC(NGT)=I
60    CONTINUE
      IF(IGMAX.EQ.0) RETURN
      NG1=MIN(NGMAX,NDV+5)
      N1=NG1-NGT
      IF(NCON.LT.NG1) N1=NCON-NGT
      IF(N1.LT.1) RETURN
      DO 90 I=1,N1
         GMAX=-1.0E+15
         KK=0
         DO 80 J=1,NCON
            DO 70 K=1,NGT
               IF(IC(K).EQ.J) GO TO 80
70          CONTINUE
            GG=G(J)
            IF(GG.GT.GMAX) THEN
               GMAX=GG
               KK=J
            ENDIF
80       CONTINUE
         IF(KK.EQ.0) GO TO 100
         NGT=NGT+1
         IC(NGT)=KK
90    CONTINUE
100   CONTINUE
      RETURN
      END
      SUBROUTINE DOT204 (NDB,B,NRB)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(NRB,*)
      ZRO=1.0E-10
      ZRO10=10.*ZRO
      DO 60 I=1,NDB
         BMX=B(I,I)
         IF (ABS(BMX).GT.ZRO10) GO TO 20
         DO 10 J=1,NDB
            B(I,J)=0.
            B(J,I)=0.
10       CONTINUE
         GO TO 60
20       BMX=1./BMX
         DO 30 J=1,NDB
            B(I,J)=BMX*B(I,J)
30       CONTINUE
         DO 50 II=1,NDB
         IF (II.EQ.I) GO TO 50
         B1=B(II,I)
         DO 40 J=1,NDB
            B(II,J)=B(II,J)-B1*B(I,J)
40       CONTINUE
         B(II,I)=-BMX*B1
50        CONTINUE
         B(I,I)=BMX
60    CONTINUE
      RETURN
      END
      SUBROUTINE DOT205 (NDV,NAVC,NAVCR,A,B,XTMP,GTMP,IC,DG)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NDV,*),B(*),XTMP(*),GTMP(*),IC(*),DG(*)
      CALL DOT214(NAVCR,NAVCR,1,B,GTMP)
      DO 40 I=1,NDV
         XTMP(I)=0.
40    CONTINUE
      JJ=0
      KK=0
      DO 50 J=1,NAVC
         JJ=JJ+1
         K=IC(J)
         IF(DG(K).LT.-0.001) GO TO 50
         KK=KK+1
         GG=GTMP(KK)
         DO 60 I=1,NDV
            XTMP(I)=XTMP(I)-A(I,JJ)*GG
60       CONTINUE
50    CONTINUE
      RETURN
      END
      SUBROUTINE DOT206 (N1,N2,A,NDV,GGL)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NDV,*),GGL(*)
      ZRO=1.0E-10
      DO 20 J=1,N2
         AMXI=ZRO
         DO 10 I=1,N1
            AI=ABS(A(I,J))
            IF (AI.GT.AMXI) AMXI=AI
10       CONTINUE
         GGL(J)=1./AMXI
20    CONTINUE
      RETURN
      END
      SUBROUTINE DOT207 (NDV,NAVC,A,B)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NDV,*),B(*)
      IJ=0
      DO 20 J=1,NAVC
         DO 20 I=1,J
            IJ=IJ+1
            B(IJ)=0.
            DO 10 K=1,NDV
               B(IJ)=B(IJ)+A(K,I)*A(K,J)
10          CONTINUE
20    CONTINUE
      RETURN
      END
      SUBROUTINE DOT208 (OBJ,OBJ1,DEL,DAB,IDEL,IDAB)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      P1M10=1.0E-10
      DENOM=MAX(ABS(OBJ1),P1M10)
      DEL1=1.-ABS(OBJ)/DENOM
      IDEL=IDEL+1
      IF (ABS(DEL1).GT.DEL) IDEL=0
      DAB1=ABS(OBJ1)-ABS(OBJ)
      IDAB=IDAB+1
      IF (ABS(DAB1).GT.DAB) IDAB=0
      RETURN
      END
      SUBROUTINE DOT209 (NDV,NAC,NVC,NACS,CT,CTMIN,PHI,GGL,G,
     1DF,A,IC,ISIDC,ITMP,B,NRB,SCALF,S,THETA,C,BETA,RSDMX,THTZ)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(*),DF(*),A(NDV,*),IC(*),S(*),C(*),ISIDC(*),
     1THETA(*),ITMP(*),B(NRB,*),GGL(*)
      ZRO=1.0E-10
      NAVC=NAC+NVC
      NACT=NAVC+NACS
      IF (NAVC.GT.0) CALL DOT211 (NAVC,CT,CTMIN,G,IC,THETA,THTZ)
      NACT1=NACT+1
      DO 10 I=1,NACT1
         C(I)=0.
10    CONTINUE
      DO 60 J=1,NDV
         D=DF(J)*SCALF
         C(NACT1)=C(NACT1)+D*D
         IF (NAVC.EQ.0) GO TO 30
         DO 20 I=1,NAVC
            C(I)=C(I)+D*A(J,I)*GGL(I)
20       CONTINUE
30       IF (NAVC.GE.NACT) GO TO 60
         I=NAVC
         DO 50 II=1,NDV
            ID=ISIDC(II)
            IF (ID.EQ.0) GO TO 50
            IF (ID.EQ.1) GO TO 40
            I=I+1
            IF (II.EQ.J) C(I)=C(I)-D
            IF (ID.LT.0) GO TO 50
40          CONTINUE
            I=I+1
            IF (II.EQ.J) C(I)=C(I)+D
50       CONTINUE
60    CONTINUE
      IF (NAVC.EQ.0) GO TO 130
      DO 80 I=1,NAVC
         GA=GGL(I)
         DO 80 J=I,NAVC
            B(J,I)=0.
            GG=GA*GGL(J)
            DO 70 K=1,NDV
               B(J,I)=B(J,I)-A(K,I)*A(K,J)*GG
70          CONTINUE
80    CONTINUE
      IF (NACS.EQ.0) GO TO 130
      J=NAVC
      DO 120 II=1,NDV
         ID=ISIDC(II)
         IF (ID.EQ.0) GO TO 120
         IF (ID.EQ.1) GO TO 100
         J=J+1
         DO 90 I=1,NAVC
            B(J,I)=A(II,I)*GGL(I)
90       CONTINUE
         IF (ID.LT.0) GO TO 120
100      CONTINUE
         J=J+1
         DO 110 I=1,NAVC
            B(J,I)=-A(II,I)*GGL(I)
110   CONTINUE
120   CONTINUE
130   CONTINUE
      IF (NACS.EQ.0) GO TO 170
      J=NAVC
      DO 160 II=1,NDV
         ID=ISIDC(II)
         IF (ID.EQ.0) GO TO 160
         JD=IABS(ID)
         DO 150 K=1,JD
            J=J+1
            DO 140 I=J,NACT
               B(I,J)=0.
               IF (ID.EQ.2) B(I,J)=1.
140         CONTINUE
         B(J,J)=-1.
150      CONTINUE
160   CONTINUE
170   CONTINUE
      IF (NAVC.EQ.0) GO TO 190
      DO 180 I=1,NAVC
         DO 180 J=1,I
            B(I,J)=B(I,J)-THETA(I)*THETA(J)
180   CONTINUE
190   CONTINUE
      IF (NVC.GT.0) GO TO 210
      NDB=NACT+1
      SAVC=C(NDB)
      DO 200 J=1,NDB
         B(NDB,J)=-C(J)
         C(J)=0.
         IF (J.GT.NAVC) GO TO 200
         B(NDB,J)=B(NDB,J)-THETA(J)
         C(J)=-THETA(J)
200   CONTINUE
      B(NDB,NDB)=-SAVC-1.
      C(NDB)=-1.
      GO TO 240
210   CONTINUE
      NDB=NACT
      IF (NAVC.EQ.0) GO TO 230
      DO 220 I=1,NAVC
         C(I)=C(I)-PHI*THETA(I)
220   CONTINUE
230   CONTINUE
240   CONTINUE
      DO 250 I=1,NDB
         DO 250 J=I,NDB
            B(I,J)=B(J,I)
250   CONTINUE
      CALL DOT210 (NDB,B,NRB,C,ITMP)
      AMULT=C(NDB)
      IF (NVC.GT.0) AMULT=1.
      AMULT=AMULT*SCALF
      DO 270 I=1,NDV
         S(I)=-AMULT*DF(I)
         IF (NAVC.EQ.0) GO TO 270
         DO 260 J=1,NAVC
            S(I)=S(I)-A(I,J)*C(J)*GGL(J)
260      CONTINUE
270   CONTINUE
      IF (NACS.EQ.0) GO TO 300
      J=NAVC
      DO 290 II=1,NDV
         ID=ISIDC(II)
         IF (ID.EQ.0) GO TO 290
         IF (ID.EQ.1) GO TO 280
         J=J+1
         S(II)=S(II)+C(J)
         IF (ID.LT.0) GO TO 290
280      CONTINUE
         J=J+1
         S(II)=S(II)-C(J)
290   CONTINUE
300   CONTINUE
      BETA=1.
      IF (NVC.GT.0) GO TO 320
      BETA=1.-C(NDB)
      IF (NAVC.EQ.0) GO TO 320
      DO 310 I=1,NAVC
         BETA=BETA-THETA(I)*C(I)
310   CONTINUE
320   CONTINUE
      CALL DOT403 (NDV,S,RSDMX)
      DO 330 I=1,NDV
         IF (ABS(S(I)).LE.ZRO) S(I)=0.
330   CONTINUE
      IF (NACS.EQ.0) GO TO 350
      DO 340 I=1,NDV
         I1=ISIDC(I)
         IF (I1.EQ.2) S(I)=0.
         IF (I1.GT.0.AND.S(I).GT.0.) S(I)=0.
         IF (I1.LT.0.AND.S(I).LT.0.) S(I)=0.
340   CONTINUE
350   CONTINUE
      IF (RSDMX.LE.ZRO) BETA=0.
      IF (NAVC.EQ.0) GO TO 370
      DO 360 I=1,NAVC
         C(I)=C(I)*GGL(I)/SCALF
360   CONTINUE
370   IF (NACS.EQ.0) GO TO 390
      NAVC1=NAVC+1
      DO 380 I=NAVC1,NACT
         C(I)=C(I)/SCALF
380   CONTINUE
390   CONTINUE
      RETURN
      END
      SUBROUTINE DOT210 (NDB,B,NRB,C,ITMP)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(NRB,*),C(*),ITMP(*)
      NMAX=5*NDB
      ZP9=0.9
      EPS=-1.0E-4
      CBMIN=0.
      ONEPM4=1.0E-4
      ONEPM6=1.0E-6
      ONEPMT=1.0E-10
      DO 10 I=1,NDB
         ITMP(I)=0
         CBMAX=0.
         BI=B(I,I)
         IF(BI.LT.-ONEPM6) CBMAX=C(I)/BI
         IF(BI.GT.EPS) EPS=BI
         IF(CBMAX.GT.CBMIN) CBMIN=CBMAX
10    CONTINUE
      EPS=ONEPM4*EPS
      IF(EPS.GT.-ONEPMT) EPS=-ONEPMT
      CBMIN=CBMIN*ONEPM6
      IF(CBMIN.LT.ONEPMT) CBMIN=ONEPMT
      ITER=0
20    ITER=ITER+1
      CBMAX=ZP9*CBMIN
      ICHK=0
      DO 30 I=1,NDB
         C1=C(I)
         BI=B(I,I)
         IF (BI.LE.EPS.AND.C1.LE.0.) THEN
            CB=C1/BI
            IF (CB.GT.CBMAX) THEN
               ICHK=I
               CBMAX=CB
            ENDIF
         ENDIF
30    CONTINUE
      IF (CBMAX.GT.CBMIN.AND.ITER.LE.NMAX) THEN
         IF (ICHK.GT.0) THEN
            JJ=ITMP(ICHK)
            IF (JJ.EQ.0) THEN
               ITMP(ICHK)=ICHK
            ELSE
               ITMP(ICHK)=0
            ENDIF
         BB=1./B(ICHK,ICHK)
         DO 40 J=1,NDB
            B(J,ICHK)=BB*B(J,ICHK)
40       CONTINUE
         C(ICHK)=CBMAX
         B(ICHK,ICHK)=BB
         DO 60 I=1,NDB
            IF (I.NE.ICHK) THEN
               BB1=B(ICHK,I)
               B(ICHK,I)=0.
               DO 50 J=1,NDB
                  B(J,I)=B(J,I)-BB1*B(J,ICHK)
50             CONTINUE
               C(I)=C(I)-BB1*CBMAX
            ENDIF
60    CONTINUE
      GO TO 20
         ENDIF
      ENDIF
      DO 80 I=1,NDB
         B(I,1)=C(I)
80    CONTINUE
      DO 90 I=1,NDB
         C(I)=0.
         J=ITMP(I)
         IF (J.GT.0) THEN
            C(I)=B(J,1)
         ELSE
            C(I)=0.
         ENDIF
90    CONTINUE
      RETURN
      END
      SUBROUTINE DOT211 (NAVC,CT,CTMIN,G,IC,THETA,THTZ)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION G(*),IC(*),THETA(*)
      THTMAX=1.
      DO 20 I=1,NAVC
         NCI=IC(I)
         GI=G(NCI)
         THT1=1.-GI/CT
         THT=THTZ
         IF(GI.GT.CTMIN) THEN
            THT=(1.0+THTZ)*THT1
         ELSE
            THT=THTZ*(THT1**2)
         ENDIF
         THETA(I)=THT
         THTMAX=MAX(THT,THTMAX)
20    CONTINUE
      IF(THTMAX.GT.1.0E-6) THEN
         DO 30 I=1,NAVC
            THETA(I)=THETA(I)/THTMAX
30       CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE DOT212 (NDV,NACS,X,XL,XU,ISIDC,SCALX)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),ISIDC(*),SCALX(*)
      NACS=0
      P9999=0.9999
      EPS=1.0E-10
      DO 30 I=1,NDV
         II=0
         JJ=0
         XI=X(I)
         SCLI=SCALX(I)
         XLI=XL(I)*SCLI
         IF(SCLI.GT.P9999) THEN
            DX=(XI-XLI)/(ABS(XLI)+EPS)
         ELSE
            DX=(XI-XLI)/(ABS(XLI)+.001)
         ENDIF
         IF (DX.LE.0.001) THEN
            NACS=NACS+1
            II=1
            JJ=-1
         ENDIF
         XUI=XU(I)*SCLI
         IF(SCLI.GT.P9999) THEN
            DX=(XUI-XI)/(ABS(XUI)+EPS)
         ELSE
            DX=(XUI-XI)/(ABS(XUI)+.001)
         ENDIF
         IF (DX.LE.0.001) THEN
            NACS=NACS+1
            JJ=1
            IF (II.GT.0) JJ=2
         ENDIF
      ISIDC(I)=JJ
30    CONTINUE
      RETURN
      END
      SUBROUTINE DOT213 (N,AK,NER)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AK(*)
      AMAX=1.0E-10
      AMIN=1.0E+10
      I1=0
      DO 5 I=1,N
      I1=I1+1
      AMIN=MIN(AK(I1),AMIN)
      AMAX=MAX(AK(I1),AMAX)
5     CONTINUE
      EPS=1.0E-4
      EPS=AMIN*EPS
      IF(EPS.LT.1.0E-10) EPS=1.0E-10
      AMAX=AMAX*1.0E+6
      IF(AMAX.GT.1.0E+10) AMAX=1.0E+10
      JJ=0
      NER=0
      IF(AK(1).LT.EPS) THEN
         NER=1
         AK(1)=AMAX
      ENDIF
      IF (N.LE.1) RETURN
      N1=N-1
      II=0
      DO 30 I=1,N1
      II=II+I
      JJ=II
      IF (AK(II).LT.EPS) THEN
         NER=NER+1
         AK(II)=AMAX
      ENDIF
      AA=AK(II)
      AA=1./AA
      I1=I+1
      DO 20 J=I1,N
      JJ=JJ+J
      IJ=JJ+I-J
      BB=AA*AK(IJ)
      IF (ABS(BB).LT.EPS) GO TO 20
      KK=JJ-J
      DO 10 K=J,N
      IK=I+KK
      JK=J+KK
      AK(JK)=AK(JK)-AK(IK)*BB
      KK=KK+K
10    CONTINUE
20    CONTINUE
30    CONTINUE
      IF (AK(JJ).LT.EPS) THEN
         NER=NER+1
         AK(JJ)=AMAX
      ENDIF
      RETURN
      END
      SUBROUTINE DOT214 (N,NMAX,NLC,AK,PL)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AK(*),PL(*)
      N1=N-1
      NP1=N+1
      EPS=1.0E-10
      M=1-NMAX
      AA=1./AK(1)
      DO 10 J=1,NLC
      M=M+NMAX
10    PL(M)=AA*PL(M)
      IF (N.LE.1) RETURN
      II=1
      DO 60 I=2,N
      JJ=II
      II=II+I
      I1=I-1
      KK=JJ+1
      AA=1./AK(II)
      IL=I-NMAX
      KK=KK+I-II
      DO 30 J=KK,I1
      JJ=JJ+1
      BB=AK(JJ)
      IF (ABS(BB).LT.EPS) GO TO 30
      IK=IL
      JK=J-NMAX
      DO 20 K=1,NLC
      IK=IK+NMAX
      JK=JK+NMAX
20    PL(IK)=PL(IK)-BB*PL(JK)
30    CONTINUE
40    CONTINUE
      IK=IL
      DO 50 K=1,NLC
      IK=IK+NMAX
50    PL(IK)=AA*PL(IK)
60    CONTINUE
      II=0
      DO 70 I=1,N1
      II=II+I
      AA=AK(II)
      M=I-NMAX
      DO 70 K=1,NLC
      M=M+NMAX
70    PL(M)=AA*PL(M)
      NN=1+N-NMAX
      DO 100 I=2,N
      IL=NN-I
      IJ=NP1-I
      II=IJ*(IJ+1)/2
      AA=1./AK(II)
      I1=I-1
      DO 90 J=1,I1
      JJ=NP1-J
      JI=JJ*(JJ+1)/2+IJ-JJ
      JK=JJ-1
      JK=JK*(JK+1)/2
      IF (JI.LE.JK) GO TO 90
      BB=AK(JI)
      IF (ABS(BB).LT.EPS) GO TO 90
      IK=IL
      JK=NN-J
      DO 80 K=1,NLC
      IK=IK+NMAX
      JK=JK+NMAX
80    PL(IK)=PL(IK)-BB*PL(JK)
90    CONTINUE
      IK=IL
      DO 100 K=1,NLC
      IK=IK+NMAX
100   PL(IK)=AA*PL(IK)
      RETURN
      END
      SUBROUTINE DOT217 (IGOTO,INFO,NDV,NCON,OBJ,X,DF,S,XL,XU,
     2IPRINT,WK,IWK,H,DX,DFG,ISIDC,G,RPRM,IPRM,IPOINT,
     3MINMAX,SCALX,NSCAL,ALPHA,ALPMAX,ALPSAV,OBJ1,SLOPF,IDAB,IDEL,
     4ITER,JSCAL,LGOTO,METHOD,ICOUNT,DF2)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),S(*),DF(*),XL(*),XU(*),WK(*),
     1IWK(*),DX(*),DFG(*),H(*),ISIDC(*),G(*),RPRM(*),IPRM(*),
     2IPOINT(*),SCALX(*)
      ZRO=1.0E-10
      DABOBJ=RPRM(3)
      DELOBJ=RPRM(4)
      DOBJ1=RPRM(5)
      DOBJ2=RPRM(6)
      DX1=RPRM(7)
      DX2=RPRM(8)
      ICNDIR=NDV
      ISCAL=IPRM(2)
      ITMAX=IPRM(3)
      ITRM=IPRM(4)
      IWRITE=IPRM(5)
      IPRNT2=IPRM(12)
      JWRITE=IPRM(13)
      IF(ISCAL.LT.0) ISCAL=ICNDIR
      IGOTO1=IGOTO+1
      GO TO (10,20,40,170), IGOTO1
10    CONTINUE
      IF (IPRINT.GT.1.AND.METHOD.LE.1) WRITE(IWRITE,300)
      IF (IPRINT.GT.1.AND.METHOD.GT.1) WRITE(IWRITE,301)
      ITER=0
      ICOUNT=0
      IDEL=0
      IDAB=0
      ALPHA=0.
      ALPSAV=0.
20    CONTINUE
      IF(JWRITE.GT.0) CALL DOT506 (ITER,NDV,NCON,
     *JWRITE,OBJ,X,SCALX,G)
30    CONTINUE
      ITER=ITER+1
      IF(ICOUNT.GE.ICNDIR.OR.NSCAL.EQ.1) ICOUNT=0
      ICOUNT=ICOUNT+1
      IF (IPRINT.GT.2) WRITE (IWRITE,340) ITER
      IGOTO=2
      INFO=2
      IF (ITER.GT.1.AND.ALPHA.LE.ZRO) GO TO 40
      RETURN
40    CONTINUE
      OBJ1=OBJ
      CALL DOT212 (NDV,NACS,X,XL,XU,ISIDC,SCALX)
      IF(IPRINT.GE.4) THEN
         ICALC=3
         JSCAL=1
         CALL DOT502(ICALC,DF,NDV,IWRITE,JSCAL,SCALX)
      ENDIF
50    CONTINUE
      IF(ICOUNT.GT.1.AND.METHOD.LE.1) THEN
         DO 110 I=1,NDV
            DFG(I)=DF(I)-DFG(I)
            DX(I)=X(I)-DX(I)
110      CONTINUE
      ENDIF
      CALL DOT218 (NDV,ICOUNT,DF,DX,DFG,H,S,SCALS,ISIDC,METHOD,DF2)
      JSCAL=0
      ICALC=2
      IF (IPRINT.GT.4) CALL DOT502 (ICALC,S,NDV,IWRITE,JSCAL,SCALX)
      CALL DOT404 (NDV,DF,S,SLOPF)
      ZROM=-ZRO
      IF (SLOPF.LT.ZROM.OR.ICOUNT.EQ.1) GO TO 150
      ICOUNT=0
      DO 140 I=1,NDV
         IF (ISIDC(I).EQ.0) GO TO 50
140   CONTINUE
      GO TO 200
150   CONTINUE
      IF(ITER.EQ.1) THEN
         CALL DOT219(OBJ,SLOPF,DOBJ1,DOBJ2,NDV,DX1,DX2,X,S,ALPHA)
         ALPSAV=ALPHA
      ENDIF
      ALPHA=.5*(ALPHA+ALPSAV)
      IF (ALPHA.LT..001) ALPHA=.001
      CALL DOT112 (NDV,X,S,XL,XU,ALPMAX,SCALX)
      IF (ALPMAX.LE.ZRO) GO TO 200
      IF (ALPHA.GT.ALPMAX) ALPHA=ALPMAX
      CALL DOT401 (NDV,X,DX)
      CALL DOT401 (NDV,DF,DFG)
      LGOTO=0
      IF (IPRNT2.GT.0) WRITE (IWRITE,220) ALPHA
170   CONTINUE
      CALL DOT100 (LGOTO,NDV,NCON,OBJ,X,XL,XU,G,WK,
     1IWK,ALPHA,ALPMAX,RPRM,IPRM,IPOINT,IPRNT2,SLOPF,MINMAX)
      IF (LGOTO.EQ.0) GO TO 180
      IGOTO=3
      INFO=1
      RETURN
180   CONTINUE
      ALPSAV=ALPHA
      IF (ALPHA.LE.ZRO) ICOUNT=0
      OBJT=OBJ
      IF(MINMAX.GT.0) OBJT=-OBJT
      IF (IPRINT.EQ.2) WRITE (IWRITE,350) ITER,OBJT
      IF (IPRINT.GE.3) WRITE (IWRITE,330) OBJT
      IF (IPRNT2.GT.0) WRITE (IWRITE,360) ALPHA,OBJT
      IF (IPRINT.GE.2) THEN
         ICALC=1
         JSCAL=1
         CALL DOT502 (ICALC,X,NDV,IWRITE,JSCAL,SCALX)
      ENDIF
      IF(JWRITE.GT.0) CALL DOT506 (ITER,NDV,NCON,
     *JWRITE,OBJ,X,SCALX,G)
      CALL DOT208 (OBJ,OBJ1,DELOBJ,DABOBJ,IDEL,IDAB)
      IF (ITER.GE.ITMAX) GO TO 200
      IF(IDEL.EQ.1.OR.IDAB.EQ.1) ICOUNT=0
      IF (IDEL.LT.ITRM.AND.IDAB.LT.ITRM) GO TO 30
200   CONTINUE
      IF (IPRINT.EQ.0) GO TO 210
      WRITE (IWRITE,370)
      WRITE (IWRITE,230)
      IF (ITER.GE.ITMAX) WRITE (IWRITE,240)
      IF (IDEL.GE.ITRM) WRITE (IWRITE,250) ITRM
      IF (IDAB.GE.ITRM) WRITE (IWRITE,260) ITRM
      IF (ALPMAX.LE.ZRO) WRITE (IWRITE,270)
210   CONTINUE
      IGOTO=0
      RETURN
220   FORMAT (/,2X,'PROPOSED ALPHA =',1PE12.5)
230   FORMAT (/,2X,'TERMINATION CRITERIA')
240   FORMAT (/,2X,'NUMBER OF ITERATIONS = ITMAX')
250   FORMAT (/,2X,'RELATIVE CONVERGENCE CRITERION WAS MET FOR',I3,
     1' CONSECUTIVE ITERATIONS')
260   FORMAT (/,2X,'ABSOLUTE CONVERGENCE CRITERION WAS MET FOR',I3,
     1' CONSECUTIVE ITERATIONS')
270   FORMAT (/,2X,'KUHN-TUCKER CONDITIONS ARE SATISFIED')
300   FORMAT (///,2X,
     1'--- BEGIN UNCONSTRAINED OPTIMIZATION: BFGS METHOD')
301   FORMAT (///,2X,
     1'--- BEGIN UNCONSTRAINED OPTIMIZATION: FLETCHER-REEVES METHOD')
330   FORMAT (/,2X,'OBJECTIVE =',1PE12.5)
340   FORMAT (//,2X,'-- BEGIN ITERATION',I3)
350   FORMAT (//,2X,'-- ITERATION',I5,3X,'OBJ =',1PE12.5)
360   FORMAT (/,2X,'CALCULATED ALPHA =',1PE12.5,//,2X,'OBJECTIVE =',
     *1PE12.5)
370   FORMAT (////,2X,'-- OPTIMIZATION IS COMPLETE')
      END
      SUBROUTINE DOT218 (NDV,ICOUNT,DF,DX,DFG,H,S,SCALS,ISIDC,METHOD,
     1DF2)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DF(*),DX(*),DFG(*),H(*),S(*),ISIDC(*)
      ZRO=1.0E-10
      IF (ICOUNT.GT.1) GO TO 40
10    CONTINUE
      IF(METHOD.LE.1) THEN
         CALL DOT406(NDV,H)
      ELSE
         DF2=0.
         DO 20 I=1,NDV
            DFD=DF(I)**2
            IF(ISIDC(I).GE.1.AND.DF(I).GE.0.) DFD=0.
            IF(ISIDC(I).LE.-1.AND.DF(I).LE.0.) DFD=0.
            DF2=DF2+DFD
20       CONTINUE
      ENDIF
      CALL DOT401 (NDV,DF,S)
      ICOUNT=1
      GO TO 90
40    CONTINUE
      IF(METHOD.LE.1) THEN
         CALL DOT405 (NDV,H,DFG,S)
         CALL DOT404 (NDV,DX,DFG,SIGMA)
         IF(ABS(SIGMA).LT.ZRO) GO TO 10
         CALL DOT404 (NDV,DFG,S,TAU)
         SIGMB=1./SIGMA
         BB=(1.+TAU*SIGMB)*SIGMB
         NA=0
         DO 50 J=1,NDV
            DO 50 I=1,J
               NA=NA+1
            H(NA)=H(NA)+BB*DX(I)*DX(J)-SIGMB*(S(I)*DX(J)+S(J)*DX(I))
50       CONTINUE
         CALL DOT405 (NDV,H,DF,S)
      ELSE
         DFDF=0.
         DO 70 I=1,NDV
            DFD=DF(I)**2
            IF(ISIDC(I).GE.1.AND.DF(I).GE.0.) DFD=0.
            IF(ISIDC(I).LE.-1.AND.DF(I).LE.0.) DFD=0.
            DFDF=DFDF+DFD
70       CONTINUE
         BETA=SCALS*DFDF/DF2
         DF2=DFDF
         DO 60 I=1,NDV
            S(I)=DF(I)-BETA*S(I)
60       CONTINUE
      ENDIF
90    CONTINUE
      DO 100 I=1,NDV
         S(I)=-S(I)
         IF (ISIDC(I).EQ.-1.AND.S(I).LT.0.) S(I)=0.
         IF (ISIDC(I).EQ.1.AND.S(I).GT.0.) S(I)=0.
         IF (ISIDC(I).EQ.2) S(I)=0.
100   CONTINUE
      CALL DOT403(NDV,S,SCALS)
      RETURN
      END
      SUBROUTINE DOT219 (OBJ,SLOPF,DOBJ1,DOBJ2,NDV,DX1,DX2,X,S,ALPHA)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),S(*)
      ZRO=1.0E-10
      DENOM=SLOPF
      IF (ABS(SLOPF).LT.ZRO) DENOM=ZRO
      ANUM=(DOBJ1*OBJ)
      IF (ANUM.LT.0.01) ANUM=0.01
      ALP1=ANUM/DENOM
      ALPHA=ABS(ALP1)
      ALP1=DOBJ2/DENOM
      ALP1=ABS(ALP1)
      ALPHA=MIN(ALPHA,ALP1)
      DO 10 I=1,NDV
      SI=ABS(S(I))
      IF (SI.LT.ZRO) SI=ZRO
      DX=X(I)
      ALP=DX*DX1/SI
      ALP=ABS(ALP)
      ALP2=DX2/SI
      ALP2=ABS(ALP2)
      IF (ALP2.GT.ALP) ALP=ALP2
      ALPHA=MIN(ALPHA,ALP)
10    CONTINUE
      RETURN
      END
      SUBROUTINE DOT220 (NDV,A,NGT,IC,NAC,NVC,G,IWRITE,IPRINT,CTMIN,GGL,
     *IER)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NDV,*),IC(*),G(*),GGL(*)
      IER=0
      NAVC=NAC+NVC
      NAC1=NAC
      NVC1=NVC
      ICHK=0
      DO 5 J=1,NAVC
         SCALJ=GGL(J)
         DO 6 I=1,NDV
            IF(ABS(A(I,J)*SCALJ).GT.1.0E-10) GO TO 5
6        CONTINUE
         ICJ=IC(J)
         IF(G(ICJ).LE.CTMIN) NAC=NAC-1
         IF(G(ICJ).GT.CTMIN) NVC=NVC-1
         IF(G(ICJ).GT.CTMIN) IER=1
         IC(J)=0
         ICHK=1
         IF(IPRINT.GT.3.OR.G(ICJ).GT.CTMIN) WRITE(IWRITE,110)ICJ
         IF(G(ICJ).GT.CTMIN) WRITE(IWRITE,120)
5     CONTINUE
      IF(NAVC.LE.1) GO TO 80
      EPS=1.0E-7
      DO 30 J=2,NAVC
         JG=IC(J)
         IF(JG.EQ.0) GO TO 30
         GJ=G(JG)
         JM1=J-1
         SCALJ=GGL(J)
         DO 20 K=1,JM1
            KG=IC(K)
            IF(KG.EQ.0) GO TO 20
            SCALK=GGL(K)
            DO 10 I=1,NDV
               AIJ=A(I,J)*SCALJ
               AIK=A(I,K)*SCALK
               IF(ABS(AIJ).LT.EPS.AND.ABS(AIK).LT.EPS) GO TO 10
               IF(AIJ*AIK.LT.0.) GO TO 20
               IF(ABS(AIJ-AIK).GT.0.01) GO TO 20
10          CONTINUE
            ICHK=ICHK+1
            GK=G(KG)
            IF(GK.LT.GJ.OR.JG.EQ.0) THEN
               IF(IPRINT.GT.3) WRITE(IWRITE,100)IC(J),IC(K),IC(K)
               IC(K)=0
               IF(GK.LE.CTMIN) NAC=NAC-1
               IF(GK.GT.CTMIN) NVC=NVC-1
            ELSE
               IF(IPRINT.GT.3) WRITE(IWRITE,100)IC(J),IC(K),IC(J)
               IC(J)=0
               IF(GJ.LE.CTMIN) NAC=NAC-1
               IF(GJ.GT.CTMIN) NVC=NVC-1
               GO TO 30
            ENDIF
20       CONTINUE
30    CONTINUE
      IF(ICHK.EQ.0) GO TO 80
      ICHK=0
      LCHK=0
      DO 50 J=1,NGT
         IF(IC(J).EQ.0) THEN
            LCHK=LCHK-1
            IC(J)=LCHK
         ENDIF
         JCHK=IC(J)
         IF(JCHK.LT.0) THEN
            ICHK=JCHK
         ELSE
            IF(ICHK.EQ.0) GO TO 50
            JM=J+ICHK
            IF(JM.LE.0) GO TO 50
            DO 40 I=1,NDV
               A(I,JM)=A(I,J)
40          CONTINUE
            GGL(JM)=GGL(J)
         ENDIF
50    CONTINUE
      ICHK=0
      DO 70 J=1,NGT
         JCHK=IC(J)
         IF(JCHK.LT.0) THEN
            ICHK=JCHK
         ELSE
            IF(ICHK.EQ.0) GO TO 70
            JM=J+ICHK
            IC(JM)=IC(J)
         ENDIF
70    CONTINUE
80    CONTINUE
      NGT=NGT+NAC+NVC-NAC1-NVC1
      RETURN
100   FORMAT(/,2X,'THE GRADIENT OF CONSTRAINT NUMBERS',I10,
     *' AND',I10,' ARE DEPENDENT'/,2X,
     *'CONSTRAINT NUMBER',I10,' IS REMOVED FROM THE ACTIVE SET')
110   FORMAT(/,2X,'THE GRADIENT OF CONSTRAINT NUMBER',I10,
     *' IS ZERO'/,2X,'CONSTRAINT IS REMOVED FROM THE ACTIVE SET')
120   FORMAT(2X,'** THIS IS A VIOLATED CONSTRAINT',/,2X,
     *'** A FATAL ERROR WILL BE GENERATED')
      END
      SUBROUTINE DOT221 (NGT,IC,G,NGMAX,NAC,NVC,CTMIN)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IC(*),G(*)
      IF(NGT.LE.NGMAX) RETURN
      NN=NGT-NGMAX
      DO 20 I=1,NN
         NGTM1=NGT-1
         GMIN=1.0E+30
         IJ=0
         DO 10 J=1,NGT
            ICJ=IC(J)
            IF(ICJ.GT.0.AND.G(ICJ).LT.GMIN) THEN
               GMIN=G(ICJ)
               IJ=J
            ENDIF
10       CONTINUE
         IF(GMIN.GT.CTMIN) THEN
            NVC=NVC-1
         ELSE
            NAC=NAC-1
         ENDIF
         IF(IJ.LT.NGT) THEN
            DO 15 J=IJ,NGTM1
               IC(J)=IC(J+1)
15          CONTINUE
            IC(NGT-I+1)=0
         ELSE
            IC(NGT)=0
         ENDIF
20    CONTINUE
      NGT=NGMAX
      RETURN
      END
      SUBROUTINE DOT222 (NDV,NAVC,A,B,IC,DG,NAVCR)
C
C                 ***** (C) COPYRIGHT *****
C     ***** VANDERPLAATS RESEARCH & DEVELOPMENT, INC. *****
C         ***** ALL RIGHTS RESERVED, WORLDWIDE *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NDV,*),B(*),IC(*),DG(*)
      IJ=0
      JJ=0
      NAVCR=NAVC
      DO 30 J=1,NAVC
         K=IC(J)
         JJ=JJ+1
         IF(DG(K).LT.-0.001) THEN
            NAVCR=NAVCR-1
            GO TO 30
         ENDIF
         II=0
         DO 20 I=1,J
            K=IC(I)
            II=II+1
            IF(DG(K).LT.-0.001) GO TO 20
            IJ=IJ+1
            B(IJ)=0.
            DO 10 K=1,NDV
               B(IJ)=B(IJ)+A(K,II)*A(K,JJ)
10          CONTINUE
20       CONTINUE
30    CONTINUE
      RETURN
      END
