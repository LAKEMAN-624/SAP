      SUBROUTINE EVAL(OBJ,X,G)
      !DEC$ ATTRIBUTES DLLEXPORT::EVAL
       !DEC$ ATTRIBUTES STDCALL,ALIAS:'EVAL'::EVAL
      !DEC$ ATTRIBUTES REFERENCE::OBJ,X,G
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(*),G(*)
    
      OBJ=2.0*X(2)*X(1)+2.0*X(3)*X(1)+4.0*X(2)*X(3)
      G(1)=1.0-0.5*X(1)*X(2)*X(3)
      RETURN
      END