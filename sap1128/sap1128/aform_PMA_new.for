c2010.6.6改，试图算三个桁架算例（上海特邀报告）  
c2010.8.4改为RIA和PMA混用

c原形为d:\Ph.D\YP_AFORM\MIX\AFORM.for.原形定稿时间:2004.1.5
c改动的地方:在计算可靠度的基础上加优化功能
c           定稿时间:2004.1.5
C
C    void JGFX(int nZYD,int nJD,int nDY,int nExc,double *X,double *Y,
C              double *Ela,double *Area,double *Exc, int *II,int *JJ,
C              int *IP,int *JP,int *IExc,double WY[],double YL[])

!      IMPLICIT REAL*8(A-H,O-Z)
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

C      PROGRAM OPTANA
      INTEGER*4 OPTANA
C.... Master routine
      INCLUDE 'COMMY.INC'
C
C       结构分析需用的常数和变量
	PARAMETER (N_JD=80,N_DY=100,N_Exc=10,N_WYYS=50,N_BDK=30)
      INTEGER*4 IRET,IJW1,IJW2,IJW3,IJW4,IJW5,IJW6,IJW7
      INTEGER*4 IFW1,IFW2,IFW3,IFW4,IFW5,IFW6,IFW7,IFW8,IFW9,IFW10
	INTEGER*4 IFW11,IFW12,IFW13
      INTEGER*4  nZYD, nJD, nDY, nExc, nWYYS
      INTEGER*4 JWORK(N_JD*2+N_DY*2+N_Exc+N_WYYS+1)
      REAL*8 FWORK(N_JD*4+N_DY*7+N_Exc+N_WYYS+1)
C
C       可靠度计算需用的常数和变量
      PARAMETER (NRWORK=140000,NIWORK=200,NCWORK=400,NUMPF=200)
c 随机变量最多200个，设计变量最多200个，概率约束最多200个
	REAL*8 RWORK, XVARNA, BETA,  GU
	REAL(kind=8) BETAT
      COMMON /CHANAM/ NVAR,NPAR,CVARNA(200),XVARNA(200)
      COMMON/QQQQ/NNEW
!!!!!!!!!!!     
      INTEGER*4 INIT
      COMMON/XXOO/ INIT
      DIMENSION IWORK(NIWORK),RWORK(NRWORK)
      CHARACTER*100 CWORK(NCWORK),CVARNA,JNAME
      CHARACTER*100 INAME
C  计算G函数计算的次数
	INTEGER*4 NOFG,NoBETA,IRIAPMA,IFEM,IJKM
      DIMENSION NOFG(NUMPF) 
      common/nfea/nctg,nctf,nctdg,nctdf,ngfun,nffun 
C
C       调用DOT进行优化需用的常数和变量
      PARAMETER (NPARE=NCWORK-NIWORK)!,NUMG=N_DY+N_WYYS)  !10年7-25改
      PARAMETER (NRWK=8000,NRIWK=2000)  !可改，原为800和200
	REAL*8 RPRM, P, PL, PU, OBJ, G, WK    !RPRM, P, PL, PU, OBJ, G, WK要作为DDOT的参变量使用,必须为双精度
c	REAL*4 RPRM, P, PL, PU, OBJ, G, WK   !为Yang05算例改为DOT,成而为单精度，可改
      DIMENSION RPRM(20),IPRM(20),P(NPARE),PL(NPARE),PU(NPARE),
     -          G(NUMPF),WK(NRWK),IWK(NRIWK)
	REAL*8 x12, y12, RuoLA,SPAR1,SPAR2,SPAR3,SPAR4  !SPAR1-SPAR4要作为PARCAI, PARCAL的参变量使用,必须为双精度
	REAL*8 AA, BB, DS, DEPS, DSAT
	DIMENSION AA(NPARE,NUMPF), BB(NPARE,NUMPF),BETAT(NUMPF)

	REAL*8   PI

      PI=3.1415926535898D0
      nJD=0
      nZYD=0
      nExc=0
      nWYYS=0
      nDY=0
      IFW1=2

C  
C  
      DO 9 I=1,NUMPF
9        NOFG(I)=0
      INIT=1
      IER=0
	NoDOT=0  !调用DOT次数
	nffun=0  !总的函数计算次数，如结构分析次数
	NoBETA=0 !？？？？？？？？？？？？？？？？
C     .FILE UNITS
      MFEM = 40
      MREL = 41
      MRE1 = 42
      MRE2 = 43
      MRE3 = 44
      MRE4 = 45   !为三种桁架增加,批处理时使用
      MRE5 = 46   !为批处理增加，一个方便直接抄答案到文档

C     .GET JOBNAME
      CALL JOBNAM(IBNAME,JNAME)
      IF(JNAME.EQ.'      ')THEN
	WRITE(*,*) ' NO JOBNAME SPECIFIED'
       STOP
      ENDIF
!---------------------------------------------------------------------------------------------
c      注意三个系数的设置有先后
	IFEM=0    !该系数为 1 基于桁架结构的分析， 0  显式功能函数
	IJKM=3  !该系数为 1 仅做桁架结构分析， 2 仅做可靠度分析   3   确定性桁架结构优化或RBDO 
	IRIAPMA=1 !该系数为 0 确定性桁架结构优化 1 RIA方法   2 PMA方法  
!---------------------------------------------------------------------------------------------
C     .         OPEN FILES  2016：要弄清楚每个文件里写的是什么
      INAME=JNAME(1:IBNAME)//'.REL'
      OPEN(MREL,FILE=INAME,STATUS='OLD',ERR=888)
      INAME=JNAME(1:IBNAME)//'.RE1'
      OPEN(MRE1,FILE=INAME)
      INAME=JNAME(1:IBNAME)//'.RE2'
      OPEN(MRE2,FILE=INAME)
      INAME=JNAME(1:IBNAME)//'.RE3'
      OPEN(MRE3,FILE=INAME)
      OPEN(MRE4,FILE='result.re1',ACCESs='APPEND')  !为yang05算例改,该语句应是批处理时使用,可改
      OPEN(MRE5,FILE='result.re2',ACCESs='APPEND')  !为yang05算例改,该语句应是批处理时使用,可改

	
      IF(IFEM.EQ.0)THEN
      
      IJW1=2;IJW2=2;IJW3=2;IJW4=2;IJW5=2;IJW6=2;IJW7=2;
       IFW1=2;IFW2=2;IFW3=2;IFW4=2;IFW5=2;IFW6=2;IFW7=2;IFW8=2;IFW9=2;
       IFW10=2;IFW11=2;IFW12=2;IFW13=2;
      END IF
      GOTO 2999  !表没有桁架，显式功能函数
      
C**************************************
C              结构分析
C**************************************
      INAME=JNAME(1:IBNAME)//'.FEM' 
      OPEN(MFEM,FILE=INAME,STATUS='OLD',ERR=888)

      CALL QDATE(MRE3)
      CALL INPUT01(nZYD, nJD, nDY, nExc, nWYYS,IER)  
      IF(IER.NE.0)then 
	WRITE(*,*) 'Wrong numbers in .FEM file!'
	GOTO 999
	endif
       IJW1=2
!      IJW1= 2   !2014年改，灰掉换下句  !所有JWORK应该要改的，-1变+1，2016年：不改
!     IJW1=0 
!c        I nodal number
      IJW2=IJW1+nDY  
!c        J nodal number
      IJW3=IJW2+nDY  
!c        No of X displacement
      IJW4=IJW3+nJD  
c        No of Y displacement
      IJW5=IJW4+nJD  
c        No of excitation
      IJW6=IJW5+nExc  
c        No of displacement constrain
      IJW7=IJW6+nWYYS  

      IFW1=2   !2014年改，灰掉换下句  !所有FWORK应该要改的，-1变+1，2016年：不改
c      IFW1=0
c        X coordinate
      IFW2=IFW1+nJD  
c        Y coordinate
      IFW3=IFW2+nJD 
c        Elastic
      IFW4=IFW3+nDY 
c        Density or gravity
      IFW5=IFW4+nDY  
c        Compression
      IFW6=IFW5+nDY 
c        Tension
      IFW7=IFW6+nDY 
c        Section area
      IFW8=IFW7+nDY 
c        excitation
      IFW9=IFW8+nExc  
c        displacement constrain
      IFW10=IFW9+nWYYS  
c        Density * length
      IFW11=IFW10+nDY  
c        Displacement 
      IFW12=IFW11+nJD*2  
c        Stress
      IFW13=IFW12+nDY  

      CALL INPUT02(nZYD, nJD, nDY, nExc, nWYYS, FWORK(IFW1),FWORK(IFW2),   !这里应该要改的，+2，2014年，2016: 还是IJW1=2,不用改了
     -  FWORK(IFW3),FWORK(IFW4),FWORK(IFW5),FWORK(IFW6),FWORK(IFW7),
     -  FWORK(IFW8),FWORK(IFW9),FWORK(IFW10),JWORK(IJW1),JWORK(IJW2),
     -  JWORK(IJW3),JWORK(IJW4),JWORK(IJW5),JWORK(IJW6),IER)
      IF(IER.NE.0)then
	WRITE(*,*) 'Wrong values in .FEM file!'
	GOTO 999
	endif
!      CALL JGFX(nZYD, nJD, nDY, nExc, FWORK(IFW1-1),FWORK(IFW2-1),   !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
!     -   FWORK(IFW3-1),FWORK(IFW7-1),FWORK(IFW8-1),JWORK(IJW1-1),
!     -   JWORK(IJW2-1),JWORK(IJW3-1),JWORK(IJW4-1),JWORK(IJW5-1),
!     -   FWORK(IFW11-1), FWORK(IFW12-1))

      WRITE(MRE3,910)
  910 FORMAT(/,' Nodal displacement of X and Y direction:')
      DO 10 I = 1,nJD
      WRITE(MRE3,920)I,FWORK(IFW11-1+2*I-1), FWORK(IFW11-1+2*I)  !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
  920 FORMAT(' NO:',I3,'  ',1PE12.4,'  ',1PE12.4)
   10 CONTINUE
      WRITE(MRE3,930)
  930 FORMAT(/,' Stress of element:')
      DO 11 I = 1,nDY
      WRITE(MRE3,940)I,FWORK(IFW12-1+I) !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
  940 FORMAT(' NO:',I3,'  ',E12.4)
   11 CONTINUE
      CALL QDATE(MRE3)
	WRITE(*,*) 'Finish structural analysis!'
      IF(IJKM.EQ.1) GOTO 999     !仅进行结构分析   2016：好像不合适，若制作确定性桁架优化，则不应进行可靠度文件读入

2999	CONTINUE
C**************************************
C                可靠度计算
C**************************************
      CALL QDATE(MRE2)
      CALL INPUT1(NVAR,NPAR,NCOR,MAXIT,IER)
      IF(IER.NE.0)then
	WRITE(*,*) 'Wrong numbers in .REL file!'
	GOTO 999
	endif
C     .            DEFINE WORKING ARRAYS
      NWORK=7500+2*NVAR
C     .            real arrays
      IWR1  = 1      
C     .       P
      IWR2  = IWR1  + NPAR
c     .       PL
      IWR3 =  IWR2 + NPAR
c     .       PU
      IWR4 = IWR3 + NPAR
C     .       CORREL
      IWR5  = IWR4  + NCOR*NCOR 
C     .       SPAR
      IWR6  = IWR5  + NVAR*4
C     .       WORKING ARRAY IN BETAE0
      IWR7  = IWR6  + NWORK
c下面为一个极限状态函数需要用的空间
C     .       X
      IWR8  = IWR7  + NVAR
C     .       ALFA
      IWR9  = IWR8  + NVAR
C     .       U
      IWR10 = IWR9  + NVAR
C     .       PSENTI
      IWR11 = IWR10  + NVAR*4
C     .       PSENTIP
      IWR12 = IWR11 + NPAR

      IWR14=NVAR*7+NPAR
      IWR13=IWR12+(NUMPF-1)*IWR14 !最多算NUMPF个概率约束

      IF (IWR13.GT.NRWORK) THEN
	WRITE(*,*) ' Real Working array RWORK too small. ',
     -      'It should at least be:',IWR13
	GOTO 999
      ENDIF
C    .             integer arrays
      IWI1=1
C     .       ITYPE
      IWI2  = IWI1  + NVAR 
c
      IF (IWI2.GT.NIWORK) THEN
	WRITE(*,*) ' Integer Working array IWORK too small. ',
     -      'It should at least be:',IWI2
	GOTO 999
      ENDIF
C     .            character arrays
      IWC1  = 1
C     .       STONAM
      IWC2  = IWC1 + NVAR
C      .       PARNAM
      IWC3  = IWC2 + NPAR
C
      IF (IWC3.GT.NCWORK) THEN
	WRITE(*,*) ' Character Working array CWORK too small. ',
     -      'It should at least be:', IWC3
	GOTO 999
      ENDIF
C
      CALL INPUT2(NVAR,NPAR,NCOR,RWORK(IWR4),IWORK(IWI1),
     -  RWORK(IWR5),RWORK(IWR1),RWORK(IWR2),RWORK(IWR3),
     -  CWORK(IWC1),CWORK(IWC2),IER)
      IF(IER.NE.0)then
	WRITE(*,*) 'Wrong values in .REL file!'
	GOTO 999
	endif
C
	WRITE(MRE1,1011)
 1011   FORMAT(/,70('*'),/,'    OUTPUT:')

      IF(IJKM.EQ.3) GOTO 3999     !RBDO，所以先不用进行下面的可靠性分析
c	DO 12345 III=1, NPAR !为三种桁架而加  2016：否则，应灰掉？
	IF (IRIAPMA.EQ.1)THEN
	
	CALL BETACA(NVAR,NPAR,NCOR,RWORK(IWR4),NWORK,IWORK(IWI1),
     -  RWORK(IWR5),RWORK(IWR6),RWORK(IWR7),RWORK(IWR8),
     -  RWORK(IWR9),RWORK(IWR10),RWORK(IWR1),RWORK(IWR11),
     -  CWORK(IWC1),CWORK(IWC2),MAXIT,BETA,nZYD,nJD,nDY,nExc,nWYYS,
     -  FWORK(IFW1-1),FWORK(IFW2-1),FWORK(IFW3-1),FWORK(IFW4-1), !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
     -  FWORK(IFW5-1),FWORK(IFW6-1),FWORK(IFW7-1),FWORK(IFW8-1),
     -  FWORK(IFW9-1),FWORK(IFW10-1),JWORK(IJW1-1),JWORK(IJW2-1),
     -  JWORK(IJW3-1),JWORK(IJW4-1),JWORK(IJW5-1), JWORK(IJW6-1),
     -  FWORK(IFW11-1),FWORK(IFW12-1),1,NOFG(1)) !记得修改功能函数序号
	ELSEIF (IRIAPMA.EQ.2)THEN
	CALL BETAE0_PMA(NVAR,NPAR,NCOR,RWORK(IWR4),NWORK,IWORK(IWI1),
     -  RWORK(IWR5),RWORK(IWR6),RWORK(IWR7),RWORK(IWR8),RWORK(IWR9),
     -  RWORK(IWR10),RWORK(IWR1),RWORK(IWR11),CWORK(IWC1),CWORK(IWC2),
     -  MAXIT,6.00D0,nZYD,nJD,nDY,nExc,nWYYS,FWORK(IFW1+1),   !记得修改目标可靠指标
     -  FWORK(IFW2+1),FWORK(IFW3+1),FWORK(IFW4+1),FWORK(IFW5+1),
     -  FWORK(IFW6+1),FWORK(IFW7+1),FWORK(IFW8+1),FWORK(IFW9+1), !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
     -  FWORK(IFW10+1),JWORK(IJW1+1),JWORK(IJW2+1),JWORK(IJW3+1),  
     -  JWORK(IJW4+1),JWORK(IJW5+1), JWORK(IJW6+1),FWORK(IFW11+1),
     -  FWORK(IFW12+1),1,NOFG(1),NoBETA,GU)  !记得修改功能函数序号,曾为III
	ELSE
	  write(*,*) 'please check IRIAPMA'
      ENDIF
      CALL QDATE(MRE2)
c12345 CONTINUE !为三种桁架而加

	WRITE(*,*) 'Finish reliability analysis!'

      IF(IJKM.EQ.2) GOTO 999     !仅进行可靠度求解

3999	CONTINUE
C**************************************
C          调用DOT进行优化
C**************************************
      NDV=NPAR
      IPRINT=3
      MINMAX=-1 !-1表min， +1表max
      METHOD=3  !3 SQP，2 SLP,  一般选3就可以，但计算144个算例时选2,这个数太重要，算三种桁架时就卡在这
C
C   经常需要修改的参数:  NCON以及目标可靠指标BETAT.
C
      NCON=1 ! 3 约束个数经常要改 NDV
      DO 12 I=1,NCON
12       BETAT(I)=2.5D0 !目标可靠指标  经常要改 3.0  3.71   3.012
      print*,'!!!!!!!!!!!!!!!!!!'
      print*,'BETAT',(BETAT(I),I=1,NCON)
      DO 13 I=1,20
         RPRM(I)=0.0
13       IPRM(I)=0
	IF(IJKM.EQ.3.and. IRIAPMA.NE.0) IPRM(1)=1 !此时为RBDO，改为1表自己提供灵敏度，否则可缺省为0，表DOT有限差分法求
c	IPRM(8)=200 !JTMAX，缺省值为50，表SLP或SQP会迭代运行几个sub-optimization-program
c         RPRM(12)=0.01   !SLP，SQP里目标函数的绝对收敛标准，缺省值为0.0001
c         RPRM(13)=0.01   !SLP，SQP里目标函数的相对收敛标准，缺省值为0.001
c         RPRM(3)=0.01   !目标函数的绝对收敛标准，缺省值为0.0001
c	   RPRM(4)=0.01   !目标函数的相对收敛标准，缺省值为0.001
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!      G(NCON+1)=0.5-P(1)/P(2)
!      G(NCON+2)=P(1)/P(2)-2.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      DO 14 I=1,NDV
         P(I)=RWORK(IWR1+I-1)
         PL(I)=RWORK(IWR2+I-1)
14       PU(I)=RWORK(IWR3+I-1)
      print *, 'P(I)=',(P(I),I=1,NDV)
      INFO=0
     
777    print *, '-A'
      print *, '2,INFO=',INFO
      CALL DOT (INFO,METHOD,IPRINT,NDV,NCON,P,PL,PU,
     *          OBJ,MINMAX,G,RPRM,IPRM,WK,NRWK,IWK,NRIWK)
	NoDOT=NoDOT+1
	 print *, 'iprm(19)',iprm(19)
	  print *, 'NoDOT',NoDOT
	 print *, '1,INFO=',INFO
	 print *, 'P(I)=',(P(I),I=1,NDV)
      IF(INFO.EQ.0)GOTO 1999 !是否优化结束
	IF(INFO.EQ.2)GOTO 666 !是否自己提供灵敏度
C****************************************************
C                以下为经常需要修改的地方， 改目标函数和功能函数
C****************************************************
C      目标函数
      nffun=nffun+1   !若为yang05例，则该句要去掉
c	OBJ=10.0*P(2)-P(1)
	OBJ=P(1)*P(2)                    !xulin4_4, 还有yang_2,xulin4_1
!	OBJ=P(1)+P(2)
!      OBJ=(P(1))**2+(P(2))**2
!      OBJ=(P(1)-3.7)**2+(P(2)-4.0)**2
!      OBJ=-((P(1)+P(2)-10)**2)/30-((P(1)-P(2)+10)**2)/120
	print *, 'obj=',OBJ
	
c	OBJ=(-(P(1)+P(2)-10)**2/30-(P(1)-P(2)+10)**2/120)   !Youn2005算例1
c	OBJ=(P(3)+2)*P(2)*P(1)**2        !Thanedar1992算例2
c	OBJ=2.828*P(1)+P(2)              !Thanedar1992算例3
c	OBJ=P(1)*sqrt(1+P(2)*P(2))       !论文算例5-5
c	OBJ=P(1)+P(2)                    !论文算例3-5，也即kkc04_1, 是易平大论文算例3.1
c	goto 19
!      DO 15 I=1,nDY
!   15  FWORK(IFW7-1+I)=P(I)     !原始设计变量：单元截面积  !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
!!   15  FWORK(IFW7-1+I)=1/P(I)  !倒数设计变量
!!	OBJ=0.D0
!	DO 16 I=1,nDY
!16	OBJ=OBJ+FWORK(IFW10-1+I)*P(I)    !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
!16	OBJ=OBJ+FWORK(IFW10-1+I)/P(I)  !倒数设计变量
c
C      确定性约束函数---无
c
C       近似可靠度约束函数  
19    CONTINUE
!------------------------------------------------------------------
      IF(NVAR.NE.NPAR)GOTO 17  !2016:特别注意，经常要改，设计变量完全独立于随机变量时，则跳过下面语句，goto 17
!	 print *, 'PMA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)'
	DO 20 I=1,NPAR
!      print *, 'I',I
!      print *, 'RWORK(IWR5-1+I)',RWORK(IWR5-1+I)
!      print *, 'RWORK(IWR5+NVAR-1+I)',RWORK(IWR5+NVAR-1+I)
!      print *, 'RWORK(IWR5+2*NVAR-1+I)',RWORK(IWR5+2*NVAR-1+I)
!      print *, 'RWORK(IWR5+3*NVAR-1+I)',RWORK(IWR5+3*NVAR-1+I)
      RWORK(IWR1+I-1)=P(I)
!       print *, 'I,P(I) RWORK(IWR1+I-1)',I,P(I),RWORK(IWR1+I-1)
 	CALL PARCAI(IWORK(IWI1-1+I),SPAR1,SPAR2,SPAR3,SPAR4,   !仅对均值为设计变量时这之后到20 continue之间的语句才有用，否则全灰掉
     -       RWORK(IWR5-1+I),RWORK(IWR5+NVAR-1+I),
     -       RWORK(IWR5+2*NVAR-1+I),RWORK(IWR5+3*NVAR-1+I))
!       print *, 'I,PARCAI',I
!      print *, 'SPAR1',SPAR1
!      print *, 'SPAR2',SPAR2
!      print *,  'SPAR3',SPAR3
!      print *,  'SPAR4',SPAR4
C	CALL PARCAI(IWORK(IWI1-1+I+4),SPAR1,SPAR2,SPAR3,SPAR4, !仅对xulin4_2，该题的第5、6、7个随机变量的均值作为设计变量
c     -       RWORK(IWR5-1+I+4),RWORK(IWR5+NVAR-1+I+4),
c     -       RWORK(IWR5+2*NVAR-1+I+4),RWORK(IWR5+3*NVAR-1+I+4))
	SPAR1=P(I)
!	print *, 'I,P(I)SPAR1',I,SPAR1
c	SPAR2=P(I)*0.1/1.732   !仅对均值为设计变量，且标准差随均值变，而变异系数不变,为KKC99_1而改,就是论文例5-5
	CALL PARCAL(MRE2,IWORK(IWI1-1+I),SPAR1,SPAR2,SPAR3,SPAR4,
     -       RWORK(IWR5-1+I),RWORK(IWR5+NVAR-1+I),
     -       RWORK(IWR5+2*NVAR-1+I),RWORK(IWR5+3*NVAR-1+I),2)
!      print *, 'I,PARCAL',I
!      print *, 'RWORK(IWR5-1+I)',RWORK(IWR5-1+I)
!      print *, 'RWORK(IWR5+NVAR-1+I)',RWORK(IWR5+NVAR-1+I)
!      print *, 'RWORK(IWR5+2*NVAR-1+I)',RWORK(IWR5+2*NVAR-1+I)
!      print *, 'RWORK(IWR5+3*NVAR-1+I)',RWORK(IWR5+3*NVAR-1+I)
c	CALL PARCAL(MRE2,IWORK(IWI1-1+I+4),SPAR1,SPAR2,SPAR3,SPAR4, !仅对xulin4_2
c     -       RWORK(IWR5-1+I+4),RWORK(IWR5+NVAR-1+I+4),
c     -       RWORK(IWR5+2*NVAR-1+I+4),RWORK(IWR5+3*NVAR-1+I+4),2)
20    CONTINUE
!------------------------------------------------------------------
c      G(1)=-(P(1)*P(1)*P(2)/20.0-1.0)   !以下三个为做KKC04_1的确定性优化而做  2016：也是08年Liang的SLA的算例1，里面有确定性优化结果
c      G(2)=-((P(1)+P(2)-5.0)**2/30.0+(P(1)-P(2)-12.0)**2/120.0-1.0) !若不是做确定性优化，这4行全隐藏
c      G(3)=-(80.0/(P(1)**2+8.0*P(2)+5.0)-1.0)
c      GO TO 777

      
      
17      DO I=1,NCON !1   !!Thanedar1992算例3,否则改为NCON

c   ****************为5、10、15杆桁架的确定性优化而加****************
c   **********此时，IPRM(1)要改为缺省值0，DOT自动求灵敏度************
!
!	IF (IRIAPMA.EQ.0)THEN !确定性结构优化
!      CALL JGFX(nZYD, nJD, nDY, nExc, FWORK(IFW1-1),FWORK(IFW2-1),   !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
!     -   FWORK(IFW3-1),FWORK(IFW7-1),FWORK(IFW8-1),JWORK(IJW1-1),  
!     -   JWORK(IJW2-1),JWORK(IJW3-1),JWORK(IJW4-1),JWORK(IJW5-1),
!     -   FWORK(IFW11-1), FWORK(IFW12-1))
!c   	IF(FWORK(IFW12-1+I)>0) THEN !可改，灰掉表不考虑压杆稳定 为赵秀丽论文而灰掉
!	MAXT=FWORK(IFW6-1+I)  !许可拉应力
!c	ELSE 
!c	MAXT=PI*FWORK(IFW3-1+I)*FWORK(IFW7-1+I)/4.0/
!c     -                     ((FWORK(IFW10-1+I))**2)  !压杆稳定临界应力
!c	ENDIF
!	G(I)=(ABS(FWORK(IFW12-1+I))-MAXT )/1.0E10
!	ENDIF
c   **********此时，IPRM(1)要改为缺省值0，DOT自动求灵敏度************
c   ****************为5、10、15杆桁架的确定性优化而加****************

      IMM=(I-1)*IWR14
      print *, '1,IMM',IMM
	IF (IRIAPMA.EQ.1)THEN
	PRINT*,'INIT',INIT              !
	CALL BETACA(NVAR,NPAR,NCOR,RWORK(IWR4),NWORK,IWORK(IWI1),
     -  RWORK(IWR5),RWORK(IWR6),RWORK(IWR7+IMM),RWORK(IWR8+IMM),
     -  RWORK(IWR9+IMM),RWORK(IWR10+IMM),P,RWORK(IWR11+IMM),
     -  CWORK(IWC1),CWORK(IWC2),MAXIT,BETA,nZYD,nJD,nDY,nExc,nWYYS,
     -  FWORK(IFW1-1),FWORK(IFW2-1),FWORK(IFW3-1),FWORK(IFW4-1),   !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
     -  FWORK(IFW5-1),FWORK(IFW6-1),FWORK(IFW7-1),FWORK(IFW8-1),
     -  FWORK(IFW9-1),FWORK(IFW10-1),JWORK(IJW1-1),JWORK(IJW2-1),
     -  JWORK(IJW3-1),JWORK(IJW4-1),JWORK(IJW5-1), JWORK(IJW6-1),
     -  FWORK(IFW11-1),FWORK(IFW12-1),I,NOFG(I))                 !记得修改功能函数序号
      
      G(I)=1.0-BETA/BETAT(I)
!       G(I)=1.0-abs(BETA)/BETAT(I)
      print *, 'PMA -FINISH BETACA'
      print *, 'I',I
      print *, 'BETA',BETA
      print *, 'BETAT',BETAT(I)
      print *, 'G(I)',G(I)
	ELSEIF (IRIAPMA.EQ.2)THEN
	CALL BETAE0_PMA(NVAR,NPAR,NCOR,RWORK(IWR4),NWORK,IWORK(IWI1),
     -  RWORK(IWR5),RWORK(IWR6),RWORK(IWR7+IMM),RWORK(IWR8+IMM),
     -  RWORK(IWR9+IMM),RWORK(IWR10+IMM),P,RWORK(IWR11+IMM),
     -  CWORK(IWC1),CWORK(IWC2),MAXIT,BETAT(I),nZYD,nJD,nDY,nExc,nWYYS, 
     -  FWORK(IFW1+1),FWORK(IFW2+1),FWORK(IFW3+1),FWORK(IFW4+1),   !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
     -  FWORK(IFW5+1),FWORK(IFW6+1),FWORK(IFW7+1),FWORK(IFW8+1),
     -  FWORK(IFW9+1),FWORK(IFW10+1),JWORK(IJW1+1),JWORK(IJW2+1),
     -  JWORK(IJW3+1),JWORK(IJW4+1),JWORK(IJW5+1), JWORK(IJW6+1),
     -  FWORK(IFW11+1),FWORK(IFW12+1),I,NOFG(I),NoBETA,GU)         !记得修改功能函数序号   这里NOFG错误？应该用NOFG(I)，14年写
      G(I)=-1.0*GU           !后四个变量含义分别是指示功能函数序号，该功能函数计算次数，NoBETA，概率功能度量
	ENDIF
	PRINT*,'-----------------------------------'
	print *, 'PMA -B'
	print *, 'iprm(19) -2'
      print *, 'iprm(19)',iprm(19)
      ENDDO
      INIT=0

      print *, 'iprm(19) -3'
      print *, 'iprm(19)',iprm(19)
   
      GO TO 777

666	CONTINUE !自己提供灵敏度
      PRINT*,'-----------------------------------'
	print *, 'PMA 666'
	print *, 'P(I)',(P(I),I=1,NPAR)
	print *, 'G(I)',(G(I),I=1,NCON)
	print *, 'OBJ',OBJ
	WRITE(MRE1,*)'*** P ***'
	WRITE(MRE1,905) (P(I),I=1,NPAR)
c	WRITE(MRE1,905) (1/P(I),I=1,NPAR)  !!倒数设计变量
	WRITE(MRE1,*)'*** CONSTRAIN ***'
	WRITE(MRE1,905) (G(I),I=1,NCON)
	WRITE(MRE1,*)'*** OBJECT FUNCTION= ', OBJ
  905 FORMAT(E12.4)!输出控制，要根据算例随时调整
!与Gfuncxcu中INIT==1相对应
!-------------------------------------------------       
      
!       DO I=1,NCON
!        IMM=(I-1)*IWR14
!      print *, '11,IMM',IMM
!      IF (IRIAPMA.EQ.1)THEN
!      CALL BETACA(NVAR,NPAR,NCOR,RWORK(IWR4),NWORK,IWORK(IWI1),
!     -  RWORK(IWR5),RWORK(IWR6),RWORK(IWR7+IMM),RWORK(IWR8+IMM),
!     -  RWORK(IWR9+IMM),RWORK(IWR10+IMM),P,RWORK(IWR11+IMM),
!     -  CWORK(IWC1),CWORK(IWC2),MAXIT,BETA,nZYD,nJD,nDY,nExc,nWYYS,
!     -  FWORK(IFW1-1),FWORK(IFW2-1),FWORK(IFW3-1),FWORK(IFW4-1),   !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
!     -  FWORK(IFW5-1),FWORK(IFW6-1),FWORK(IFW7-1),FWORK(IFW8-1),
!     -  FWORK(IFW9-1),FWORK(IFW10-1),JWORK(IJW1-1),JWORK(IJW2-1),
!     -  JWORK(IJW3-1),JWORK(IJW4-1),JWORK(IJW5-1), JWORK(IJW6-1),
!     -  FWORK(IFW11-1),FWORK(IFW12-1),I,NOFG(I))                 !记得修改功能函数序号
!      print *, '666 BETACA'
!       print *, 'BETA',BETA
!      print *, 'BETAT',BETAT(I)
!      print *, 'G(I)',G(I)
!      ELSEIF (IRIAPMA.EQ.2)THEN
!	CALL BETAE0_PMA(NVAR,NPAR,NCOR,RWORK(IWR4),NWORK,IWORK(IWI1),
!     -  RWORK(IWR5),RWORK(IWR6),RWORK(IWR7+IMM),RWORK(IWR8+IMM),
!     -  RWORK(IWR9+IMM),RWORK(IWR10+IMM),P,RWORK(IWR11+IMM),
!     -  CWORK(IWC1),CWORK(IWC2),MAXIT,BETAT(I),nZYD,nJD,nDY,nExc,nWYYS, 
!     -  FWORK(IFW1+1),FWORK(IFW2+1),FWORK(IFW3+1),FWORK(IFW4+1),   !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
!     -  FWORK(IFW5+1),FWORK(IFW6+1),FWORK(IFW7+1),FWORK(IFW8+1),
!     -  FWORK(IFW9+1),FWORK(IFW10+1),JWORK(IJW1+1),JWORK(IJW2+1),
!     -  JWORK(IJW3+1),JWORK(IJW4+1),JWORK(IJW5+1), JWORK(IJW6+1),
!     -  FWORK(IFW11+1),FWORK(IFW12+1),I,NOFG(I),NoBETA,GU)         !记得修改功能函数序号   这里NOFG错误？应该用NOFG(I)，14年写
!      G(I)=-1.0*GU           !后四个变量含义分别是指示功能函数序号，该功能函数计算次数，NoBETA，概率功能度量
!	ENDIF
!       ENDDO
!------------------------------------------------- 
c     目标函数对设计变量的灵敏度
!	DO 21 I=1,NDV         
!21	WK(I)=FWORK(IFW10-1+I)  !XULIN4_5 对一般桁架算例     !这里应该要改的，-1变+1，2014年，2016: 还是IJW1=2,不用改了
c21	WK(I)=-FWORK(IFW10-1+I)/P(I)/P(I)  !!倒数设计变量
	WK(1)=P(2)            !xulin4_4, 还有yang_2,xulin4_1
	WK(2)=P(1)
!	WK(1)=2*P(1)            
!	WK(2)=2*P(2)
!       WK(1)=-(P(1)+P(2)-10)/15-(P(1)-P(2)+10)/60
!       WK(2)=-(P(1)+P(2)-10)/15+(P(1)-P(2)+10)/60
c	WK(1)=sqrt(1+P(2)*P(2)) !论文算例5.5 两杆桁架算例
C	WK(2)=-P(1)*P(2)/sqrt(1+P(2)*P(2))
c	WK(1)=2*(P(3)+2)*P(2)*P(1) !Thanedar1992算例2
C	WK(2)=(P(3)+2)*P(1)**2
C	WK(3)=P(2)*P(1)**2
c	WK(1)=2.828                !Thanedar1992算例3
c	WK(2)=1.0
! 	WK(1)=1.0                  !论文算例3-5，也即kkc04_1，是易平大论文算例3.1
!	WK(2)=1.0
c	WK(1)=-(P(1)+P(2)-10)/15-(P(1)-P(2)+10)/60  !Youn2005算例1
c	WK(2)=-(P(1)+P(2)-10)/15+(P(1)-P(2)+10)/60
	print *, 'WK(I)',(WK(I),I=1,NDV)
c下面灵敏度为XULIN4_3
cc
C      确定性约束函数对设计变量的灵敏度---无
c
C       近似可靠度约束函数对设计变量的灵敏度
	NGT=IPRM(20)
	print *, 'PMA NGT',NGT
	IF(NGT==0)goto 777
      DO J=1,NCON !1   !!Thanedar1992算例3,否则改为NCON
      IMM=(J-1)*IWR14
      print *, '2,IMM',IMM
	DO 34 I=1, NDV
	BB(I,J)=-RWORK(IWR11+I-1+IMM) !设计变量与随机变量无关
!	 
!	BB(I,J)=-RWORK(IWR10+I-1+IMM) !随机变量均值为设计变量
!	BB(I,J)=-RWORK(IWR10+I-1+4+IMM)
      print *, 'BB',I,J,BB(I,J)
	IF (IRIAPMA.EQ.1) BB(I,J)=BB(I,J)/BETAT(J)
	print *, 'BB/BETAT(J)',I,J,BB(I,J)
34	CONTINUE
c      BB(11,J)=-RWORK(IWR11+11-1+IMM)!专为Lee3而改!!严重注意这里可能会出错，可改，
c      BB(12,J)=-RWORK(IWR11+12-1+IMM)
      ENDDO

C      STORE APPROPRIATE GRADIENTS IN ARRAY AA
      print *, 'NGT=', NGT
      DO 40 K=1,NGT
      DO 41 I=1,NDV
	J=IWK(K)
	print *, 'K=', K
      print *, 'J=', J
41	AA(I,K)=BB(I,J)
40    CONTINUE
C       PUT THE GRADIENTS IN THE WK ARRAY
	N1=NDV
	DO 60 K=1,NGT
	DO 50 I=1,NDV
50    WK(I+N1)=AA(I,K)
	N1=N1+NDV
60     CONTINUE
       print *,' ------------------'
        print *,'NGT',NGT
      print *,' wk(i+ndv)=',(wk(i+ndv),i=1,ndv*ngt)
      print *,' wk(i)',( wk(i),i=1,ndv*(ngt+1))
      print *,' iprm',iprm
      print *,' IWK(K)',(IWK(K),k=1,NGT)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      DO 160 I=1,NP
!      PG=P(I)
!      IF(ABS(PG).GT.1.D-10)THEN
!	 DPEPS=1.D-4*DABS(PG)
!      ELSE
!!	 DPEPS=1.D-6
!	 DPEPS=1.D-4
!      ENDIF
!      P(I)=P(I)+DPEPS
!      CALL ZUTOX(N,NCOR,U,RL,SPAR,ITYPE,X)
!      CALL UPDATX (N,NP,STONAM,PARNAM,X,P)
!      CALL GFUNCX (N,X,NP,P,GG,nZYD,nJD,nDY,nExc,nWYYS,XX,YY,
!     -             Ela,Ruo,Sigc,Sigt,Area,Exc,WYYS,RuoL,
!     -             II,JJ,IP,JP,IExc,IWYYS,WY,YL,IFlag,NOFG)
!      DP(N+I)=(GG-G)/DPEPS
!     
!      P(I)=PG
!   
!  160 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	GOTO 777
	
      GOTO 1999
  888 WRITE(*,*) 'File "',INAME,'" does not exist'

 1999 CONTINUE !优化结束
	WRITE(*,*) 'Finish optimization!'
      NumFG=0
      DO 8 I=1,NUMPF
8        NumFG=NumFG+NOFG(I)
      WRITE(*,*) 'NoDOT= ',NoDOT
      WRITE(*,*) 'nffun= ',nffun, ' NumFG= ',NumFG

      WRITE(MRE4,*) '优化计算结果 P(I)'  !输出控制，要根据算例随时修改
	WRITE(MRE4,906) (P(I),I=1,NPAR)
      WRITE(MRE4,*) '约束 G(I)'  
	WRITE(MRE4,906) (G(I),I=1,NCON)
      WRITE(MRE4,*) '目标函数  DOT Functioncalls   函数计算次数'  
	WRITE(MRE4,907)  OBJ, nffun,NumFG
      WRITE(MRE5,*) '优化计算结果 P(I)'  !输出控制，要根据算例随时修改
	WRITE(MRE5,908) (P(I),I=1,NPAR)
      WRITE(MRE5,*) '约束 G(I)'  
	WRITE(MRE5,908) (G(I),I=1,NCON)
      WRITE(MRE5,*) '目标函数  DOT Function calls   函数计算次数'  
	WRITE(MRE5,907)  OBJ, nffun,NumFG
  906 FORMAT(E12.4,',',E12.4,',',E12.4,',',E12.4,',',E12.4,',')  
  907 FORMAT(E12.4,     I5,  I11)  
  908 FORMAT(E12.4)  
C  906 FORMAT(E12.4, 2F8.3,F8.5,I5, I11)  
c	OBJ=-(4.558+1.9645-10)**2/30-(4.558-1.9645+10)**2/120   !Youn2005算例1
c      WRITE(*,*) 'OBJ= ',OBJ
  999 CONTINUE
      STOP
      END
      
 
     