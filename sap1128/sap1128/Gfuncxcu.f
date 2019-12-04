      SUBROUTINE GFUNCX(N,X,NP,P,G,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -                  Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
     -                  II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
C     .CALCULATION OF FAILURE FUNCTION  G(X) : MAIN int *IWYYS) SUBROUTINE.
C     .INPUT :
C     .  N     : NUMBER OF STOCHASTIC VARIABLES
C     .  X(N)  : STOCHASTIC VARIABLES IN BASIC SPACE
C     .  NP   : NUMBER OF CONSTANT PARAMETERS
C     .  P(NP): DETERMINISTIC PARAMETERS
C     .OUTPUT :
C     .  G     : VALUE OF FAILURE FUNCTION
!      INTERFACE TO SUBROUTINE JGFX[C,ALIAS:'_JGFX'] (nZYD,nJD,
!     &	nDY,nExc,X,Y,Ela,Area,Exc,II,JJ,IP,JP,IExc,WY, YL)
!      INTEGER*4 nZYD   [VALUE]
!      INTEGER*4 nJD    [VALUE]
!      INTEGER*4 nDY    [VALUE]
!      INTEGER*4 nExc    [VALUE]
!      REAL*8    X      [REFERENCE]
!      REAL*8    Y      [REFERENCE]
!      REAL*8    Ela    [REFERENCE]
!      REAL*8    Area   [REFERENCE]
!      REAL*8    Exc    [REFERENCE]
!      INTEGER*4 II     [REFERENCE]
!      INTEGER*4 JJ     [REFERENCE]
!      INTEGER*4 IP     [REFERENCE]
!      INTEGER*4 JP     [REFERENCE]
!      INTEGER*4 IExc   [REFERENCE]
!      REAL*8    WY     [REFERENCE]
!      REAL*8    YL     [REFERENCE]
!      END

      REAL*8 X,P,G,MAXT,MAXC,Pu,Pv,A12,sig,P2,LL,Iz1,Iz2,Iz,sig1,tuo1,PI !可改，增加了很多临时变量
      DIMENSION X(N),P(NP)
c   结构分析需用的变量
      REAL*8 XX,YY,Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,WY,YL
      DIMENSION XX(nJD+1),YY(nJD+1),Ela(nDY+1),Ruo(nDY+1),Sigc(nDY+1),
     -          Sigt(nDY+1),Area(nDY+1),Exc(nExc+1),WYYS(nWYYS+1),
     -          RuoL(nDY+1),WY(nJD*2+1),YL(nDY+1)
	INTEGER*4 II,JJ,IP,JP,IExc,IWYYS
      DIMENSION II(nDY+1),JJ(nDY+1),IP(nJD+1),JP(nJD+1),
     -	      IExc(nExc+1),IWYYS(nWYYS+1)
c   临时变量

      REAL*8 XX1,YY1,Ela1,Ruo1,Sigc1,Sigt1,Area1,Exc1,WYYS1,RuoL1
      DIMENSION XX1(nJD+1),YY1(nJD+1),Ela1(nDY+1),Ruo1(nDY+1),
     -          Sigc1(nDY+1),Sigt1(nDY+1),Area1(nDY+1),Exc1(nExc+1),
     -          WYYS1(nWYYS+1),RuoL1(nDY+1)
C
      GOTO (101,102,103,104,105,106,107,108,109,110),IFlag   
C


!C  许林论文算例4-1
!  101 G=1.0-(4*X(2))/(P(1)*P(2)**2*X(3))-(X(1)**2)/((P(1)*P(2)*X(3))**2)
  101 G=1.0-(4*X(3))/(P(1)*P(2)**2*X(2))-(X(1)**2)/((P(1)*P(2)*X(2))**2)
!  101 G=0.2*P(1)*P(2)*((X(2))**2)-X(1)
!  101 G=(X(1))**3+ ((X(1))**2)*X(2)+(X(2))**3-18
!  101 G=-X(1)*sin(4*X(1))-1.1*X(2)*sin(2*X(2))   
!  101 G=-EXP(X(1)-7)-X(2)+10
!  101 G=0.3*X(1)*X(1)*X(2)-X(2)+0.8*X(1)+1
!  101 G=X(3)-(600D0*X(2)/(P(1)*P(2)*P(2))+600D0*X(1)/(P(1)*P(1)*P(2)))
!  101 G=X(1)/P(2)
!  101 G=X(1)-X(2)*X(3)*X(3)
!  101 G=X(3)-(600.*X(2)/2/4/4+600.*X(1)/2/2/4)
!  101 G=X(1)*X(2)-78.12*X(3)
!  101 G=X(1)*X(1)+X(2)*X(2)-18.
!  101 G=X(1)*X(1)*X(1)+X(1)*X(1)*X(2)+X(2)*X(2)*X(2)-18.
!  101 G=X(1)*X(1)*X(1)*X(1)+2*X(2)*X(2)*X(2)*X(2)-20.
!  101 G=X(1)*X(1)*X(1)+X(1)*X(1)*X(2)+X(2)*X(2)*X(2)-18.
!  101 G=48*X(1)*X(2)-100*X(3)*36
!  101 G=X(1)*X(1)*X(2)/20.-1.
!  101 G=0.5*P(1)*P(2)*X(2)*X(2)-X(1)
      GOTO 999
C
C  下面两算例为Yang2004 算例2
C  102 G=X(3)-(600.*X(2)/P(1)/P(2)/P(2)+600.*X(1)/P(1)/P(1)/P(2))
!  102 G=1.0-(600.*X(2)/P(1)/P(2)/P(2)+600.*X(1)/P(1)/P(1)/P(2))/X(3)
!  102 G=2.5-(4.*1000000./X(4)/P(1)/P(2))*(SQRT((X(2)/P(2)/P(2))**2
!     -   +(X(1)/P(1)/P(1))**2))
  102 G=((X(1)+X(2)-5.)**2)/30+((X(1)-X(2)-12.)**2)/120-1
!  102 G=X(1)+X(2)-3
!  102 G=1.-(0.9063*X(1)+0.4226*X(2)-6)**2-(0.9063*X(1)+0.4226*X(2)-6)**3
!     -   +0.6*(0.9063*X(1)+0.4226*X(2)-6)**4-0.4226*X(1)+0.9063*X(2)
      GOTO 999
C
C  许林论文算例4-2
  103 G=80.0/(X(1)**2+8.0*X(2)+5.0)-1.0 
!  103 PI= 3.14159265358979323846D0
!      S=7500.0D0
!	F=X(2)+X(3)+X(4)
!	As=2.0*X(5)*X(6)
!	Us=X(5)*X(6)*X(7)
!	Ui=0.5*X(5)*X(6)*(X(7)**2)
!	Eb=PI**2*X(9)*Ui/(S**2)
!      G=X(1)-F*(1.0/As+X(8)*Eb/Us/(Eb-F))
      GOTO 999

c以下三个G（X）函数为K.K.Choi2004算例1, 就是博士论文3-1
  104 G=X(1)*X(1)*X(2)/20.0-1.0
      GOTO 999
  105 G=(X(1)+X(2)-5.0)**2/30.0+(X(1)-X(2)-12.0)**2/120.0-1.0
      GOTO 999
  106 G=80.0/(X(1)**2+8.0*X(2)+5.0)-1.0
      GOTO 999

c以下G（X）函数为许林论文算例3.5
  107 G=X(1)*X(1)*X(1)+X(2)*X(2)*X(2)-18.0
c  107 G=X(1)*X(1)+X(2)*X(2)-18.0
C  107 G=X(1)*X(1)*X(1)*X(1)+X(2)*X(2)*X(2)*X(2)-18
c  107 G=X(1)*X(1)*X(1)*X(1)+2*X(2)*X(2)*X(2)*X(2)-20.0
c  107 G=X(1)*X(1)*X(1)+X(1)*X(1)*X(2)+X(2)*X(2)*X(2)-18.0
c  107 sig1=sqrt(1+P(2)*P(2))*(X(2)/P(1)+X(1)/P(1)/P(2))/200
c	sig2=sqrt(1+P(2)*P(2))*(X(2)/P(1)-X(1)/P(1)/P(2))/200
c	G=1-sig1
c      if(sig1<sig2)  G=1-sig2
      GOTO 999

c以下G（X）函数为许林论文算例4-3
  108 G=(1.0-X(4)*P(3)*X(1)/P(1)/P(2)/X(2))*P(3)*P(2)*X(1)-X(3)
      GOTO 999


c以下一个G（X）函数为K.K.Choi2005 Adaptive probability analysis 算例2
  109 G=(4-(X(1)+0.25)**2+(X(1)+0.25)**3+(X(1)+0.25)**4-X(2))
      GOTO 999

c以下G（X）函数为许林论文算例4.4
  110 G=1.0-4*X(1)/(P(1)*P(2)*P(2)*X(4))-4*X(2)/(P(1)*P(1)*P(2)*X(4))
     -     -X(3)*X(3)/((P(1)*P(2)*X(4))**2)
      GOTO 999

  999 NOFG=NOFG+1
      RETURN
      END

C上面的GFUNCX子程序经常需要根据所求例题的功能函数形式进行修改，
c下面的子程序则无需改动
c
C
      DOUBLE PRECISION FUNCTION phiLIT (X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST/ EPSMIN,EPSMAX,EMAX,SQ2PI
	 phiLIT=SQ2PI*EXP(-X*X/2.D0) 
      RETURN
      END

C
C
      SUBROUTINE GFUNCU(N,NP,U,DuG,SS2,G,NCOR,RL,ITYPE,SPAR,X,P,
     -    STONAM,PARNAM,DP,NDP,XJU,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -    Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,II,JJ,IP,JP,
     -    IExc,IWYYS,WY,YL,IFlag,NOFG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 STONAM(N),PARNAM(NP)
      DIMENSION U(N),DuG(N),RL(NCOR,NCOR),ITYPE(N),SPAR(N,4),X(N),P(NP),
     -    DP(NDP)
      DIMENSION XJU(N,N),XN(N)
         INTEGER*4 INIT1
      COMMON/XXOO/ INIT1
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
      
      
      PRINT*,'GFUNCU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      PRINT*,'U(I)',(U(I),I=1,N)
      PRINT*,'X(I)',(X(I),I=1,N)
      PRINT*,'RL(I,J)',((RL(I,J),J=1,NCOR),I=1,NCOR)
      CALL ZUTOX(N,NCOR,U,RL,SPAR,ITYPE,X)
!      PRINT*,'GFUNCU -D'
      PRINT*,'U(I)',(U(I),I=1,N)
      PRINT*,'X(I)',(X(I),I=1,N)
!      PRINT*,'GFUNCU -P'
!      PRINT*,'P(I)',(P(I),I=1,NP)
      CALL UPDATX (N,NP,STONAM,PARNAM,X,P)
      CALL GFUNCX (N,X,NP,P,G,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -             Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
     -             II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
!      PRINT*,'GFUNCU -P2'
!      PRINT*,'P(I)',(P(I),I=1,NP)
!      DEPS=1.D-6  !可改，可能为1.D-4  20181130
      DO 30 J=1,N
      DEPS = dmax1( 1.D-6, 1.D-6*dabs(u(j)) ) 
!      PRINT*,'J,DEPS',J,DEPS
      
      U(J)=U(J)+DEPS
      CALL ZUTOX(N,NCOR,U,RL,SPAR,ITYPE,XN)
!      PRINT*,'GFUNCU -E'
!      PRINT*,'XN(I)',(XN(I),I=1,N)
      DO 28 I=1,N
   28 XJU(I,J)=(XN(I)-X(I))/DEPS
      
   30 U(J)=U(J)-DEPS
      PRINT*,'GFUNCU -F'
      PRINT*,'XJU(I,J)',((XJU(I,J),J=1,N),I=1,N)
      DO 34 I=1,NDP
   34 DP(I)=0.D0
      DO 35 I=1,N
   35 DuG(I)=0.D0
      DO 40 I=1,N
      XG=X(I)
!      PRINT*,'GFUNCU -G'
      
!      PRINT*,'I,X(I)',I,X(I)
      IF(ABS(XG).GT.1.D-10)THEN
!	 DXEPS=DEPS*DABS(XG)
	 DXEPS=1.D-4*DABS(XG) !2018.12.03
      ELSE
!	 DXEPS=1.D-6
	 DXEPS=1.D-4
      ENDIF
      X(I)=X(I)+DXEPS
      CALL UPDATX (N,NP,STONAM,PARNAM,X,P)
      CALL GFUNCX (N,X,NP,P,GG,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -             Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
     -             II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
      PRINT*,'G',G
      PRINT*,'GG',GG
      DO 36 J=1,N
   36 DuG(J)=DuG(J) + XJU(I,J) * (GG-G)/DXEPS
   
	   DP(I)=(GG-G)/DXEPS
   40 X(I)=XG
      PRINT*,'DP',(DP(J),J=1,N)
       PRINT*,'DuG',(DuG(J),J=1,N)
!       
      SS2=0.D0
      DO 41 I=1,N
   41 SS2=SS2+DuG(I)**2
C
C
       PRINT*,'------------------------'
      PRINT*,'GFUNCU DP(N+I)'
      PRINT*,'INIT1',INIT1
!------------------------------------------------------
!      IF(INIT1==1) GOTO 70  !与aform PMA中666处对应
!--------------------------------------------------------
      IF(NP.NE.N)THEN    !当设计变量d和随机变量x个数同，但d又不是x的统计参数时要隐去，否则出错
      DO 60 I=1,NP
      PG=P(I)
      PRINT*,'I,P(I)',I,P(I)
      IF(ABS(PG).GT.1.D-10)THEN
!	 DPEPS=DEPS*DABS(PG)
	 DPEPS=1.D-4*DABS(PG)
      ELSE
!	 DPEPS=1.D-6
	 DPEPS=1.D-4
      ENDIF
      P(I)=P(I)+DPEPS
      PRINT*,'I,P(I)+DPEPS',I,P(I)
      CALL ZUTOX(N,NCOR,U,RL,SPAR,ITYPE,X)
      CALL UPDATX (N,NP,STONAM,PARNAM,X,P)
      PRINT*,'P',P
       PRINT*,'X',X
      CALL GFUNCX (N,X,NP,P,GG,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
     -             Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
     -             II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
      PRINT*,'DPEPS',DPEPS
      PRINT*,'G',G
      PRINT*,'GG',GG
      DP(N+I)=(GG-G)/DPEPS
      PRINT*,'I,DP(N+I)',I,DP(N+I)
      P(I)=PG
   
   60 CONTINUE
   70   CONTINUE
      PRINT*,'GFUNCU'
      PRINT*,'DP(N+I)',(DP((N+I)),I=1,NP)
      ENDIF    !当设计变量d和随机变量x个数同，但d又不是x的统计参数时要隐去，否则出错

      CALL ZUTOX(N,NCOR,U,RL,SPAR,ITYPE,X)
       PRINT*,'X(I)',(X(I),I=1,N)
       PRINT*,'GFUNCU -FINISH'
      
      RETURN
      END
