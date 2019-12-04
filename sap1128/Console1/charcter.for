        SUBROUTINE S(ka, m)

        !DEC$ ATTRIBUTESDLLEXPORT :: S

        CHARACTER(m)::str

        DIMENSION ka(m)

        INTEGER i

        INTEGER x, y

        DO i=1, m

        str(i : i) = char(ka(i))

        END DO

        END