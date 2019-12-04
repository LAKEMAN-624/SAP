      PROGRAM EX1
       INTEGER,PARAMETER :: LENS=5
      DOUBLE PRECISION ::A(LENS)
      DOUBLE PRECISION ::O(LENS)
      
      DATA A /1,2,3,4,5/
      CALL DEALARRAY(A, LENS, O)
      O=A+O
      
      END