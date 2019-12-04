      SUBROUTINE EVAL(OBJ,X,G)
      !DEC$ ATTRIBUTES DLLEXPORT::EVAL
       !DEC$ ATTRIBUTES STDCALL,ALIAS:'EVAL'::EVAL
      !DEC$ ATTRIBUTES REFERENCE::OBJ,X,G
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),G(*)
      
      OBJ=2.0*SQRT(2.)*X(1)+X(2)
      G(1)=(2.*X(1)+SQRT(2.)*X(2))/(2.*X(1)*(X(1)+
     *SQRT(2.0)*X(2)))-1.
      G(2)=1./(2.*(X(1)+SQRT(2.)*X(2)))-1.
      RETURN
      END