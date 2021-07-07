!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Subroutine generate the EOS table that used for !
! solving the initial star according to the user input !
! which assume a ideal degenerate fermi gas EOS        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EOSTABLE
USE DEFINITION
IMPLICIT NONE

! The dummy variable for density and pressure !
REAL(DP) :: den, pree1, pree2

! Integer Parameter !
INTEGER :: i, j, start, end

! The width of the table for each order of magnitue of density !
INTEGER :: width

! Initialize the width !
width = 3000

! We initialize the starting and ending order of magnitude for density !
start = -22
end = -7

! We open the EOS table for input !
OPEN (UNIT = 99, FILE = 'EOS_Table1.eos', STATUS = 'REPLACE')
OPEN (UNIT = 100, FILE = 'EOS_Table2.eos', STATUS = 'REPLACE')

! First assign the fermi constant !
CALL FINDCONST

! We do the loop the get the pressure of Fermi Gas !
! For each density input !
DO i = start, end
	DO j = 0, width - 1

		! We initialize the density input !
		den = (10.0E0_DP ** (DBLE(i+1)) - 10.0E0_DP ** (DBLE(i))) * DBLE(j) / DBLE(width) + 10.0E0_DP ** (DBLE(i))

		! We get the pressure !
		CALL GETRHO_EOSRTOP (pree1, den, gs1, mb1, me1, ye1, 1)
		CALL GETRHO_EOSRTOP (pree2, den, gs2, mb2, me2, ye2_old, 2)

		! We print the pressure and density to EOS table !
		WRITE (99, *) pree1, den
		WRITE (100, *) pree2, den

	END DO
END DO

! We assign the eosline number !
eoslineno = width * (end - start + 1)

! We close the file !
close(99)
close(100)

END SUBROUTINE