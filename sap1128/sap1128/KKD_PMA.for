      SUBROUTINE BETAE0_PMA(NVAR,NPAR,NCOR,CORREL,NWORK,ITYPE,   !�ó���Ϊ�����ʹ��ܶ�����������ȼ�������
     -    SPAR,BWORK,X,ALFA,U,PSENTI,P,PSENTIP,STONAM,PARNAM,
     -    MAXIT,BETA,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,           !MAXITΪ�ڲ�ɿ���Լ�����������������Ϊ1����ʾSAP
     -    Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,II,JJ,IP,JP,
     -    IExc,IWYYS,WY,YL,IFlag,NOFG,NOBETA,G)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*8 STONAM(NVAR),PARNAM(NPAR)
      COMMON/CMACHE/EPSS
      INCLUDE 'COMMY.INC'
      DIMENSION X(NVAR),U(NVAR),SPAR(NVAR,4),CORREL(NCOR,NCOR),
     -          ITYPE(NVAR),ALFA(NVAR),
     -          PSENTI(NVAR,4),BWORK(NWORK),P(NPAR),PSENTIP(NPAR)
c   �ṹ�������õı���
      REAL*8 XX,YY,Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,WY,YL
      DIMENSION XX(nJD+1),YY(nJD+1),Ela(nDY+1),Ruo(nDY+1),Sigc(nDY+1),
     -          Sigt(nDY+1),Area(nDY+1),Exc(nExc+1),WYYS(nWYYS+1),
     -          RuoL(nDY+1),WY(nJD*2+1),YL(nDY+1)
	INTEGER*4 II,JJ,IP,JP,IExc,IWYYS
      DIMENSION II(nDY+1),JJ(nDY+1),IP(nJD+1),JP(nJD+1),
     -	      IExc(nExc+1),IWYYS(nWYYS+1)
      EPS = 1.D-14

      IWRITE=2  !�ɸģ�Ϊ0��д��MRE2�к���Ļ�����Ϊ1��д�룬Ϊ2��д����࣬��Ϊ������Ӧ��0

	 IF(NOFG.EQ.0)THEN
       DO 46 I=1,NVAR
   46   U(I)=0.0D0  !��ֵ���ɸ�,-1.1D0, 0.1D0����0.0Ч������õ�
	 ENDIF
C
      EPSS=EPS
      NDP=NPAR+NVAR
      M1=1
c       DuG
      M2=M1+NVAR
c       XJU
      M3=M2+NVAR*NVAR
c       DP
	M4=M3+NDP
      IF(M4.GT.NWORK)THEN
	WRITE(*,*) ' NWORK TOO SMALL 1',M4,NWORK
	STOP
      ENDIF
C
	NOBETA=NOBETA+1
      CALL BETAE1_PMA(NVAR,NPAR,X,U,SPAR,BWORK(M1),ITYPE,CORREL,NCOR,
     -  BETA,G,SQ2,IWRITE,MRE2,ITER,P,STONAM,PARNAM,BWORK(M3),
     -  NDP,BWORK(M2),MAXIT,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,Ela,Ruo,
     -  Sigc,Sigt,Area,Exc,WYYS,RuoL,II,JJ,IP,JP,IExc,IWYYS,WY,YL,
     -  IFlag,NOFG,NOBETA,BETAs)
C
      CALL SENANA_PMA_X(NVAR,NPAR,X,U,SPAR,ITYPE,CORREL,NCOR, !�ɸ�ΪSENANA_PMA_U
     -  BETAs,SQ2,PSENTI,PSENTIP,BWORK(M3),NDP,BWORK(M1))
Cд��MRE2�к���Ļ���
      IF(IWRITE.GE.1)THEN
      WRITE(MRE2,*) '*******  STOP IN BETAE0!  *********'
	WRITE(MRE2,1060) ITER,NOFG,G,BETAs
	WRITE(MRE2,1000)
	DO 40 I=1,NVAR
	  ALFA(I)=U(I)/BETA
	  WRITE(MRE2,905) I,ALFA(I),U(I),X(I)
   40   CONTINUE
	WRITE(MRE2,1070)
	DO 45 I=1,NVAR
   45   WRITE(MRE2,1080) I,(PSENTI(I,J),J=1,4)
	WRITE(MRE2,1075)
	DO 55 I=1,NPAR
   55   WRITE(MRE2,1080) I,PSENTIP(I)
C      CALL QDATE(MRE2)
c
	WRITE(*,1060) ITER,NOFG,G,BETAs
	WRITE(*,4020) 
	DO 56 I = 1,NVAR
   56   WRITE(*,4022) ALFA(I),U(I),X(I),
     -                PSENTI(I,1),PSENTI(I,2)
	WRITE(*,4040)
	DO 57 I = 1,NPAR
   57   WRITE(*,4030) PARNAM(I),PSENTIP(I)
	ENDIF
Cд��MRE1��  
c	WRITE(MRE1,1059) (P(I),I=1,NPAR)!ΪYang05�������˸����
	WRITE(MRE1,1060) ITER,NOFG,G,BETAs
	WRITE(MRE1,4020) 
	DO 62 I = 1,NVAR
   62	WRITE(MRE1,4022) ALFA(I),U(I),X(I),
     -                PSENTI(I,1),PSENTI(I,2)
	WRITE(MRE1,4040)
	DO 63 I = 1,NPAR
   63	WRITE(MRE1,4030) PARNAM(I),PSENTIP(I)
      CALL QDATE(MRE1)
C
  905 FORMAT(I5,2F9.4,1PE12.4)
 1000 FORMAT('    I   ALFA',7X,'U',8X,'X')
 1059 FORMAT(/,'design : ',4F8.4)
 1060 FORMAT(/,'ITER: ',I5,' NOFG: ',I5,' PerformanceM: ',1PE12.4,
     -' Betas :',1PE12.4)
 1070 FORMAT(/,
     -' SENSITIVITY COEFFICIENTS FOR STATISTICAL PARAMETERS:',
     -' (D Gp / D PAR,I) ',/,                                          !�ı�
     -'        PAR,1       PAR,2       PAR,3       PAR,4')
 1075 FORMAT(/,
     -' SENSITIVITY COEFFICIENTS FOR CONSTANTS:(d Gp / D PAR,I) ',/,     !�ı�
     -'        PAR,1       ')
 1080 FORMAT(I4,4(1PE12.4))
 4020 FORMAT(/,'   Alfa-coeff.          U           X      ',
     -         '   d Gp/d mu    d Gp/d sigma')                           !�ı�
 4022 FORMAT('',E14.5,'  ',E14.5,' ',E14.5,'  ',E14.5,'  ',E14.5)
 4040 FORMAT(/,'             d Gp/d p')                                 !�ı�
 4030 FORMAT(' ',A8,'  ',E14.5,'  ',E14.5,'  ',F8.3)
      RETURN
      END
C
C
      SUBROUTINE BETAE1_PMA(N,NP,X,U,SPAR,DuG,IT,RL,NCOR,BETA,G,
     -  SQ2,IWRITE,MRE2,ITER,P,STONAM,PARNAM,DP,NDP,XJU,MAXIT,
     -  nZYD,nJD,nDY,nExc,nWYYS,XX,YY,Ela,Ruo,Sigc,Sigt,Area,Exc,
     -  WYYS,RuoL,II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG,NOBETA,BETAs)
C********************************************************************
C
C     .CALCULATION OF RELIABILITY INDEX.
C
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*8 STONAM(N),PARNAM(NP)
      DIMENSION X(N),U(N),SPAR(N,4),DuG(N),RL(NCOR,NCOR),NORM1(N)
      DIMENSION IT(N),P(NP),DP(NDP), XJU(N,N),XTT(N),UTT(N)
      COMMON/CMACHE/EPS
      COMMON/TIMES/NOB1,NOB2,NOB3
	DIMENSION ALFA(N), gu(N),DIRK(N),DIRK_1(N),DIRK_2(N),DDD(N)
c   �ṹ�������õı���
      REAL*8 XX,YY,Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,WY,YL
      DIMENSION XX(nJD+1),YY(nJD+1),Ela(nDY+1),Ruo(nDY+1),Sigc(nDY+1),
     -          Sigt(nDY+1),Area(nDY+1),Exc(nExc+1),WYYS(nWYYS+1),
     -          RuoL(nDY+1),WY(nJD*2+1),YL(nDY+1)
	INTEGER*4 II,JJ,IP,JP,IExc,IWYYS,IHUHAO,C,lambda
      DIMENSION II(nDY+1),JJ(nDY+1),IP(nJD+1),JP(nJD+1),
     -	      IExc(nExc+1),IWYYS(nWYYS+1)

!      Ϊ��C������޸�      
      REAL*8 CCC,TTT,namd,NORM1,NORM,ALPHA1,ALPHA2,C1,B0
      DIMENSION CCC(N,N),TTT(N)  !ע������ֻ����2�������������(N=2)
      DO 250 I=1,N
      DO 250 J=1,N
      IF(I.NE.J)THEN
	CCC(I,J)=0.D0
      ELSE
	CCC(I,I)=1.D0
      ENDIF
  250 CONTINUE
      CCC(1,1)=1.D0
      CCC(2,2)=1.D0
	ZETA=0.D0
	THETA1=0.0
	THETA2=0.0
	nChaos=0 !�ɸģ�0-AMV, 1-HMV��2-ԭ����cc��3-MCC,4-HCC,5-ACC  ����������PPT������1����Gp��ⷽ������ȷ��
      namd=0.1D0  !�ɸģ�namd=1.0��CΪ��λ�󣬱�ʾ�����ƣ�����Ϊ����,
      C1=2.5
      lambda=10
      NORM=0.D0
      ALPHA1=0.D0
      ALPHA2=0.D0
C
      IF(IWRITE.GE.2)THEN
	WRITE(MRE2,1000)
	DO 6 I=1,N
	 WRITE(MRE2,905) I,0.,U(I)
    6 CONTINUE
	 WRITE(MRE2,*)
	ENDIF
C
      DO 50 ITER=1,MAXIT
      CALL GFUNCU(N,NP,U,DuG,SS2,G0,NCOR,RL,IT,SPAR,X,P,STONAM,PARNAM,
     -    DP,NDP,XJU,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,Ela,Ruo,Sigc,Sigt,
     -    Area,Exc,WYYS,RuoL,II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
      SQ2=dsqrt(SS2)
C
c	IF(ITER.EQ.1)THEN  !�жϿɿ�����������������ָ������Ҫ��������,PMA��⹫ʽ��һ����
      IHUHAO=1
c	DO 7 I=1,N
c7	  UTT(I)=0.0D0
c      CALL ZUTOX(N,NCOR,UTT,RL,SPAR,IT,XTT)
c      CALL GFUNCX (N,XTT,NP,P,GTT,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
c     -             Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
c     -             II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
c      if(GTT<0.0)IHUHAO=-1 
c	ENDIF
      IF(nChaos==6) THEN
      NORM=0
	DO  I=1,N
	 NORM1(I)=U(I)-lambda*DuG(I)
      END DO    
      DO  I=1,N
	 NORM=NORM+NORM1(I)*NORM1(I)
      END DO 
       NORM=SQRT( NORM)
       IF(ITER.EQ.1)THEN      
	DO  I=1,N
	 DIRK(I)=BETA*NORM1(I)/NORM
	U(I)=DIRK(I)
      END DO    
	ELSE IF(ITER.EQ.2)THEN
	DO  I=1,N
	 DIRK_1(I)=DIRK(I)
	 DIRK(I)= BETA*NORM1(I)/NORM
	 U(I)=DIRK(I)
      END DO    
	ELSE IF(ITER.GE.3)THEN	
	DO  I=1,N
	 DIRK_2(I)=DIRK_1(I)
	 DIRK_1(I)=DIRK(I)
	 DIRK(I)= BETA*NORM1(I)/NORM
       U(I)=DIRK(I)
      END DO  
      END IF	
	IF(ITER.GE.3)THEN
	DO  I=1,N
	ALPHA1=ALPHA1+(DIRK_1(I)-DIRK(I))**2
      ALPHA2=ALPHA2+(DIRK_2(I)-DIRK_1(I))**2
	END DO
	ALPHA1=ALPHA1/(BETA**2)
	ALPHA2=ALPHA2/(BETA**2)
	IF(ALPHA1.GE.ALPHA2)THEN
	lambda=lambda/C1
	END IF
	END IF
      IF(lambda.LE.0.01)THEN
      lambda=0.01
      END IF
      GOTO 3333
      END IF
	IF(ITER.EQ.1)THEN  !�ɸģ�Ϊ��CMV��HMV���ĳɸþ�, 2014�Ķ�Ϊ��һֱ�ô˾䣬AMV��HMV������
!	IF(ITER.GE.1)THEN  !�����ô˾䣬����AMV
	
	DO 13 I=1,N
	 DIRK(I)= - DuG(I)/SQ2
	 DDD(I)=DIRK(I)
   13 CONTINUE    
	GOTO 1111
C
	ELSE IF(ITER.EQ.2)THEN
	
	DO 16 I=1,N
	 DIRK_1(I)=DIRK(I)
	 DIRK(I)= - DuG(I)/SQ2
	 DDD(I)=DIRK(I)
   16 CONTINUE
	IF(nChaos.NE.5) GOTO 1111
	THETA1=0.0             !14��ΪACC���ӣ������������������� ��֪Ϊʲô������Щ�����ܵ�������������ڲ���ACCʱҪ��ִ����Щ���
	DO 916 I=1,N   
  916	 THETA1=THETA1+DIRK_1(I)*DIRK(I)
	THETA1=DACOS (THETA1)   !14��ΪACC���ӣ�������������������
      GOTO 1111
C
	ELSE IF(ITER.GE.3)THEN  !�����Ӧ��ENDIF@@@@@@@@@@@
	
	DO 17 I=1,N
	 DIRK_2(I)=DIRK_1(I)
	 DIRK_1(I)=DIRK(I)
	 DIRK(I)= - DuG(I)/SQ2
   	 DDD(I)=DIRK(I)   
   17 CONTINUE
	IF(nChaos.NE.5) GOTO 918
	THETA2=THETA1          !14��ΪACC���ӣ�������������������
	THETA1=0.0             
	DO 917 I=1,N   
  917	 THETA1=THETA1+DIRK_1(I)*DIRK(I)
	THETA1=DACOS (THETA1)   !14��ΪACC���ӣ�������������������
  918	ZETA=0.D0
	DO 18 I=1,N
	 ZETA=ZETA+(DIRK(I)-DIRK_1(I))*(DIRK_1(I)-DIRK_2(I))
   18 CONTINUE
      IF(nChaos.eq.1)THEN  !2014Ϊ�������ģ� ���������HMV
      IF(ZETA<0.0)THEN  !CONCAVE,��CMV��������AMV
	DMODE=0.D0 
	DO 20 I=1,N
   20	 DMODE= DMODE+(DIRK(I)+DIRK_1(I)+DIRK_2(I))**2
      DMODE=DSQRT(DMODE)
	DO 21 I=1,N
   21	 DDD(I)= (DIRK(I)+DIRK_1(I)+DIRK_2(I))/DMODE 
	ENDIF              !CONCAVE,��CMV��������AMV
	DO 138 I=1,N
  138	 U(I)=IHUHAO*BETA* DDD(I)   
	GOTO 3333   
	ENDIF                  !2014Ϊ�������ģ� ���������HMV
      GOTO 1111
	ENDIF                   !�����Ӧ��ENDIF@@@@@@@@@@@

1111  continue
***
	IF (nChaos.EQ.2.OR.nChaos.EQ.3) THEN !Chaos  MCC
	DO 136 I=1,N
  136	 U(I)= U(I)+namd*(IHUHAO*BETA* DDD(I)-U(I))  !�������,CΪ��λ����namdΪ1��������ơ���IHUHAO=1��������
	IF(nChaos.EQ.3) THEN    ! 2014��ģ�MCC   
	DMODE=0.D0 
	DO 1420 I=1,N
 1420	 DMODE= DMODE+U(I)**2
      DMODE=DSQRT(DMODE)
	DO 1421 I=1,N
 1421	 U(I)= BETA*U(I)/DMODE      
	ENDIF                   ! 2014��ģ�MCC
	GOTO 3333
***
	ELSEIF(ITER<3.OR.nChaos.EQ.0)THEN !AMV
	DO 137 I=1,N
  137	 U(I)=IHUHAO*BETA* DDD(I) 
	GOTO 3333
***
	ELSEIF (nChaos.EQ.4.OR. nChaos.EQ.5) THEN !HCC,ACC 2014�ģ����׵�ENDIF%%%%%%%%%%%%%
      IF(ZETA.GE.0.0)THEN       !CONVEX,��AMV������MCC
	DO 139 I=1,N
  139	 U(I)=IHUHAO*BETA* DDD(I)   
      ELSE                      !CONVEX,��AMV������MCC     
	IF( nChaos.EQ.5)  THEN  !ACC����Ҫ�Ĳ���namd   ��������������������
	   IF((0.2*THETA1)>THETA2)  THEN
	       namd=0.2*namd
	   ELSEIF(THETA1>THETA2) THEN
	       namd=namd*THETA2/THETA1
	   ENDIF
	ENDIF                   !ACC����Ҫ�Ĳ���namd   ��������������������
	DO 140 I=1,N
	 U(I)= U(I)+namd*(IHUHAO*BETA* DDD(I)-U(I))  !�������,CΪ��λ����namdΪ1��������ơ���IHUHAO=1��������
  140 CONTINUE
c 2014��ģ������������������   
	DMODE=0.D0 
	DO 1422 I=1,N
 1422	 DMODE= DMODE+U(I)**2
      DMODE=DSQRT(DMODE)
	DO 1423 I=1,N
 1423	 U(I)= BETA*U(I)/DMODE  
c 2014��ģ������������������   
	ENDIF                    !CONVEX,��AMV������MCC
	goto 3333   
	ENDIF                                  !���׵�ENDIF%%%%%%%%%%%%%
c
!    �ĳ�C����Ϊ�����������
	DO 14 J=1,N
      TTT(J)=0.0
	DO 15 I=1,N
	 TTT(J)= TTT(J)+CCC(J,I)*(IHUHAO*BETA* DDD(I)-U(I))
   15 CONTINUE    
   14 CONTINUE    
	DO 11 I=1,N
	   U(I)= U(I)+namd*TTT(I) 
   11 CONTINUE    
!    �ĳ�C����Ϊ�����������

3333	CALL ZUTOX(N,NCOR,U,RL,SPAR,IT,X)
      CALL UPDATX (N,NP,STONAM,PARNAM,X,P)
      IF(nChaos==6) goto 3131
      CALL GFUNCX (N,X,NP,P,G,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -             Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
     -             II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
 3131 continue 
  	BETAs = 0.D0
	DO 22 I=1,N
!22	BETAs = BETAs + U(I)*DIRK(I)   !���ʲô��˼��������������Ҳ���ԣ����������������
22	BETAs = BETAs + U(I)*U(I)      !���ʲô��˼��������������Ҳ���ԣ����������������
	BETAs = DSQRT(BETAs)           !���ʲô��˼��������������Ҳ���ԣ����������������
!      IF(MAX(MAX(ABS(BETAs-BETA),ABS(G0-G)),ABS(1.0-G0/G))<1.d-3) !�Ķ�Ϊ�¾䣬��ΪBETAs����Ϊ��
!      IF(MAX(MAX(ABS(ABS(BETAs)-BETA),ABS(G0-G)),ABS(1.0-G0/G))<1.d-3) 
!     -	GOTO 60           !����׼��Ӧ�ÿɸģ�����������130    14��ģ�ΪʲôBETAs�ǵ�ǰ�Ŀɿ�ָ�ꣿ
C
       B0=0
       DO  I=1,N
       B0=B0+(DIRK_1(I)-U(I))*(DIRK_1(I)-U(I))
      END DO  
      B0=SQRT(B0)
      IF(B0<1.d-3)  GOTO 60
      IF(IWRITE.GE.2)THEN
	 WRITE(MRE2,1000)
	 DO 35 I=1,N
	 WRITE(MRE2,905) I,DuG(I)/SQ2,U(I),X(I)
   35  CONTINUE
	 WRITE(MRE2,1060) ITER,BETAs,G
	 WRITE(*,1060) ITER,BETAs,G
	 WRITE(MRE2,*)
      ENDIF
   50 CONTINUE
   60 CONTINUE
!20181216ΪSLA����
      CALL GFUNCX (N,X,NP,P,G,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -             Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
     -             II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
      WRITE(MRE2,*) '*******  STOP IN BETAE1!  *********'
  905 FORMAT(I5,2F9.4,1PE12.4)
 1000 FORMAT('    I    ALFA',4X,'U',11X,'X')
 1060 FORMAT('Outter ITERATION= ',I5,' BETAs : ',1PE12.4,' G :',1PE12.4)
      RETURN
      END
C
C
	SUBROUTINE SENANA_PMA_X(N,NP,X,U,SPAR,IT,CORREL,NCOR,
     -        BETA,SQ2,PSENTI,PSENTIP,DP,NDP,DuG)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),U(N),SPAR(N,4),CORREL(NCOR,NCOR),IT(N),
     -          PSENTI(N,4),DP(NDP),PSENTIP(NP),STAT(4),UN(N),
     -          STAP(4),DuG(N),XN(N)
      INCLUDE 'COMMY.INC'   !2016�¼ӵ�
      DEPS=1.D-6
      
      DO 100 IVAR=1,N
      IF(IT(IVAR).GT.0)THEN
c	DO 90 J=1,1  !ΪYang�����Ķ�Ϊ�þ䣬������������
	DO 90 J=1,4   
	IF(J.GT.1 .AND. SPAR(IVAR,J).EQ.0.D0)THEN
	  PSENTI(IVAR,J)=0.D0
	ELSE
	  CALL PARCAI(IT(IVAR),STAT(1),STAT(2),STAT(3),STAT(4),
     -            SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4))
	  CALL PARCAL(MRE2,IT(IVAR),STAT(1),STAT(2),STAT(3),STAT(4),
     -            SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4),2)
      CALL ZUTOX(N,NCOR,U,CORREL,SPAR,IT,X)
      print *,'SENANA -B'
           print *,'X(I)=',(X(I),I=1,N)
           print *,'U(I)=',(U(I),I=1,N)
           print *,' STAT(I)=',( STAT(I),I=1,4)
	  DSTAT=STAT(J)
	  IF(DSTAT.NE.0.D0)THEN
	    DS=DEPS*DABS(DSTAT)
	  ELSE
	    DS=DEPS
	  ENDIF
	  STAT(J)=DSTAT+DS
	   print *,'SENANA -C'
	    print *,' STAT(I)=',( STAT(I),I=1,4)
	  CALL PARCAL(MRE2,IT(IVAR),STAT(1),STAT(2),STAT(3),STAT(4),
     -        SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4),2)
      CALL ZUTOX(N,NCOR,U,CORREL,SPAR,IT,XN)
      print *,'SENANA -D'
           print *,'XN(I)=',(XN(I),I=1,N)
           print *,'U(I)=',(U(I),I=1,N)
	  DB=0.D0
	  print *,'DP(I)=',(DP(I),I=1,N)
	  DO 10 I=1,N
   10     DB=DB+DP(I)*(XN(I)-X(I))/DS
       print *,'DB=',DB
	  PSENTI(IVAR,J)=DB
	  STAT(J)=DSTAT
	  CALL PARCAL(MRE2,IT(IVAR),STAT(1),STAT(2),STAT(3),STAT(4),
     -         SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4),2)
	ENDIF
   90   CONTINUE
c	PSENTI(IVAR,2) = 0.D0 !ΪYang�����ӣ�������DO 90 J=1,4��1ʱ�����
c	PSENTI(IVAR,3) = 0.D0 !ΪYang�����ӣ�������DO 90 J=1,4��1ʱ�����
c	PSENTI(IVAR,4) = 0.D0 !ΪYang�����ӣ�������DO 90 J=1,4��1ʱ�����
      ELSE
	PSENTI(IVAR,1) = DP(IVAR)
	PSENTI(IVAR,2) = 0.D0
	PSENTI(IVAR,3) = 0.D0
	PSENTI(IVAR,4) = 0.D0
      ENDIF
  100 CONTINUE
      DO 110 IPAR=1,NP
  110 PSENTIP(IPAR) = DP(N+IPAR)
	CALL ZXTOU(N,NCOR,X,CORREL,SPAR,IT,U)
      RETURN
      END
C
C
C��SENANA_PMA_U��SENANA_PMA_X�к�����
	SUBROUTINE SENANA_PMA_U(N,NP,X,U,SPAR,IT,CORREL,NCOR,
     -        BETA,SQ2,PSENTI,PSENTIP,DP,NDP,DuG)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),U(N),SPAR(N,4),CORREL(NCOR,NCOR),IT(N),
     -          PSENTI(N,4),DP(NDP),PSENTIP(NP),STAT(4),UN(N),
     -          STAP(4),DuG(N)
      DEPS=1.D-6
      
      DO 100 IVAR=1,N
      IF(IT(IVAR).GT.0)THEN
	DO 90 J=1,4
	IF(J.GT.1 .AND. SPAR(IVAR,J).EQ.0.D0)THEN
	  PSENTI(IVAR,J)=0.D0
	ELSE
	  CALL PARCAI(IT(IVAR),STAT(1),STAT(2),STAT(3),STAT(4),
     -            SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4))
	  CALL PARCAL(MRE2,IT(IVAR),STAT(1),STAT(2),STAT(3),STAT(4),
     -            SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4),2)
	  CALL ZXTOU(N,NCOR,X,CORREL,SPAR,IT,U)
	  DSTAT=STAT(J)
	  IF(DSTAT.NE.0.D0)THEN
	    DS=DEPS*DABS(DSTAT)
	  ELSE
	    DS=DEPS
	  ENDIF
	  STAT(J)=DSTAT+DS
	  CALL PARCAL(MRE2,IT(IVAR),STAT(1),STAT(2),STAT(3),STAT(4),
     -        SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4),2)
	  CALL ZXTOU(N,NCOR,X,CORREL,SPAR,IT,UN)
	  DB=0.D0
	  DO 10 I=1,N
   10     DB=DB+DuG(I)*(UN(I)-U(I))/DS
	  PSENTI(IVAR,J)=DB
	  STAT(J)=DSTAT
	  CALL PARCAL(MRE2,IT(IVAR),STAT(1),STAT(2),STAT(3),STAT(4),
     -         SPAR(IVAR,1),SPAR(IVAR,2),SPAR(IVAR,3),SPAR(IVAR,4),2)
	ENDIF
   90   CONTINUE
      ELSE
	PSENTI(IVAR,1) = DP(IVAR)
	PSENTI(IVAR,2) = 0.D0
	PSENTI(IVAR,3) = 0.D0
	PSENTI(IVAR,4) = 0.D0
      ENDIF
  100 CONTINUE
      DO 110 IPAR=1,NP
  110 PSENTIP(IPAR) = DP(N+IPAR)
      CALL ZUTOX(N,NCOR,U,CORREL,SPAR,IT,X)
      RETURN
      END

