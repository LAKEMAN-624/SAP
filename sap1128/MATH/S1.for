      SUBROUTINE  DEALARRAY(ARRAY, LENGTH, OARRTY)
      !DEC$ ATTRIBUTES DLLEXPORT::DEALARRAY
      !DEC$ ATTRIBUTES STDCALL,ALIAS:'DEALARRAY'::DEALARRAY
      INTEGER :: LENGTH
      DOUBLE PRECISION ::ARRAY(LENGTH)
      DOUBLE PRECISION ::OARRTY(LENGTH)
      
      CALL TOARRAY(ARRAY, LENGTH, OARRTY)
      OARRTY(1)=0;
      RETURN
      END