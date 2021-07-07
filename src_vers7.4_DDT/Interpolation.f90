!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do 2-point interpolation for a given point of existing !
!datum and output the quantity that the user wished to interpolate 	 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LINEAR(x0, x1, y0, y1, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: x0, x1, y0, y1, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: nu0, nu1
REAL (DP) :: de0, de1
REAL (DP) :: l0, l1

! Assign numerator !
nu0 = (x_in - x1)
nu1 = (x_in - x0)

! Assign denominator !
de0 = (x0 - x1)
de1 = (x1 - x0)

! Assign polynominal coefficient !
l0 = nu0/de0
l1 = nu1/de1

! Compute the output !
y_out = l0*y0 + l1*y1

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do 2-point interpolation for a given point of existing !
!datum and output the quantity that the user wished to interpolate 	 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BILINEAR(x0, x1, y0, y1, f00, f10, f01, f11, x_in, y_in, f_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: x_in, y_in
REAL (DP), INTENT(IN) :: x0, x1, y0, y1
REAL (DP), INTENT(IN) :: f00, f01, f10, f11

! Input parameter !
REAL (DP), INTENT(OUT) :: f_out

! Intermediate polynominal !
REAL (DP) :: a0, a1
REAL (DP) :: b0, b1

! Assign coefficient !
a0 = (x1 - x_in)/(x1 - x0)
a1 = (x_in - x0)/(x1 - x0)
b0 = (y1 - y_in)/(y1 - y0)
b1 = (y_in - y0)/(y1 - y0)

! Compute the output !
f_out = b0*(a0*f00 + a1*f10) + b1*(a0*f01 + a1*f11)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine do 4-point interpolation for a given point of existing !
!datum and output the quantity that the user wished to interpolate 	 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CUBIC(x0, x1, x2, x3, y0, y1, y2, y3, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: x0, x1, x2, x3, y0, y1, y2, y3, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: nu0, nu1, nu2, nu3
REAL (DP) :: de0, de1, de2, de3
REAL (DP) :: l0, l1, l2, l3

! Assign numerator !
nu0 = (x_in - x1)*(x_in - x2)*(x_in - x3)
nu1 = (x_in - x0)*(x_in - x2)*(x_in - x3)
nu2 = (x_in - x0)*(x_in - x1)*(x_in - x3)
nu3  = (x_in - x0)*(x_in - x1)*(x_in - x2)

! Assign denominator !
de0 = (x0 - x1)*(x0 - x2)*(x0 - x3)
de1 = (x1 - x0)*(x1 - x2)*(x1 - x3)
de2 = (x2 - x0)*(x2 - x1)*(x2 - x3)
de3 = (x3 - x0)*(x3 - x1)*(x3 - x2)

! Assign polynominal coefficient !
l0 = nu0/de0
l1 = nu1/de1
l2 = nu2/de2
l3 = nu3/de3

! Compute the output !
y_out = l0*y0 + l1*y1 + l2*y2 + l3*y3

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Akima spline interpolation. See Hiroshi Akima 1970 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AKIMA(xm2, xm1, x0, xp1, xp2, xp3, ym2, ym1, y0, yp1, yp2, yp3, x_in, y_out)
USE DEFINITION, ONLY : DP
IMPLICIT NONE

! Input parameter !
REAL (DP), INTENT(IN) :: xm2, xm1, x0, xp1, xp2, xp3, ym2, ym1, y0, yp1, yp2, yp3, x_in

! Input parameter !
REAL (DP), INTENT(OUT) :: y_out

! Intermediate polynominal !
REAL (DP) :: dm2, dm1, d0, dp1, dp2
REAL (DP) :: s0, s1

! Weights !
REAL (DP) :: w1, w2, w3, w4

! Coefficient of polynominal !
REAL (DP) :: p0, p1, p2, p3

! Temporal arrays !
REAL (DP) :: diff, temp

! Assign slopes !
dm2 = (ym1 - ym2)/(xm1 - xm2)
dm1 = (y0 - ym1)/(x0 - xm1)
d0 = (yp1 - y0)/(xp1 - x0)
dp1 = (yp2 - yp1)/(xp2 - xp1)
dp2 = (yp3 - yp2)/(xp3 - xp2)

! Assign weights !
w1 = abs(dp1 - d0)
w2 = abs(dm1 - dm2)
w3 = abs(dp2 - dp1)
w4 = abs(d0 - dm1)

! assign slopes !
IF(w1 == 0.0D0 .AND. w2 == 0.0D0) THEN
	s0 = 0.5D0*(dm1 + d0)
ELSE
	s0 = (w1*dm1 + w2*d0)/(w1 + w2)
END IF
IF(w3 == 0.0D0 .AND. w4 == 0.0D0) THEN
	s1 = 0.5D0*(d0 + dm1)
ELSE
	s1 = (w3*d0 + w4*dp1)/(w3 + w4)
END IF

! Assign temp !
diff = xp1 - x0
temp = x_in - x0

! assign coefficients !
p0 = y0
p1 = s0
p2 = (3.0D0*d0 - 2.0D0*s0 - s1)/diff
p3 = (s0 + s1 - 2.0D0*d0)/diff**2

! Output the interpolation !
y_out = p0 + p1*temp + p2*temp**2 + p3*temp**3

END SUBROUTINE