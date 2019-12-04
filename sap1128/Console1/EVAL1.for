      PROGRAM EX12
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(2),XL(2),XU(2),G(2),
     *WK(800),IWK(200),RPRM(20),IPRM(20)
      NRWK=800
      NRIWK=200
      DO 10 I=1,20
          RPRM(I)=0.0
10      IPRM(I)=0
      METHOD=1
      NDV=2
      NCON=2
      DO 20 I=1,NDV
          X(I)=1.0
          XL(I)=0.1
 20     XU(I)=100.
      IPRINT=3
      MINMAX=-1
      INFO=0
100   CALL DOT (INFO, METHOD, IPRINT,  NDV, NCON, X,  XL,  XU,  
     * OBJ,  MINMAX,  G, RPRM, IPRM,  WK, NRWK,  IWK,  NRIWK)
      IF(INFO.EQ.0) then
      STOP
      END IF
      CALL EVAL(OBJ,X,G)
      GO TO 100
      END
      
       SUBROUTINE EVAL(OBJ,X,G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),G(*)
      
      OBJ=2.0*SQRT(2.)*X(1)+X(2)
      G(1)=(2.*X(1)+SQRT(2.)*X(2))/(2.*X(1)*(X(1)+
     *SQRT(2.0)*X(2)))-1.
      G(2)=1./(2.*(X(1)+SQRT(2.)*X(2)))-1.
      RETURN
      END SUBROUTINE