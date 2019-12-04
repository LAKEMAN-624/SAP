        DOUBLE PRECISION FUNCTION ADD(A,B)
        !DEC$ ATTRIBUTES DLLEXPORT::ADD
        !DEC$ ATTRIBUTES STDCALL,ALIAS:'Add'::ADD
            DOUBLE PRECISION:: A,B
            ADD=A+B
        END

        FUNCTION SORTANDFINDMAX(ARRAY,LENGTH)
        !DEC$ ATTRIBUTES DLLEXPORT::SORTANDFINDMAX
        !DEC$ ATTRIBUTES STDCALL,ALIAS:'Sortandfindmax'::SORTANDFINDMAX
        DOUBLE PRECISION ::ARRAY(LENGTH)
        INTEGER::I,J
        DOUBLE PRECISION::SORTANDFINDMAX,TEMP
        SORTANDFINDMAX=ARRAY(1)
        DO I=1,LENGTH-1
            DO J=I+1,LENGTH
                IF(ARRAY(I).GT.ARRAY(J)) THEN
                    TEMP=ARRAY(I)
                    ARRAY(I)=ARRAY(J)
                    ARRAY(J)=TEMP
                    SORTANDFINDMAX=ARRAY(J)
                    END IF
            END DO
            END DO
        END