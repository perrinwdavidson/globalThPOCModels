!-----------------------------------------------------------------------
!                                  llf
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to calculate the log likelihood of
! the Box-Cox parameter. This function is intended to be MINIMIZED, meaning
! that we will negate the value at the end.
!
! Derivation:
!
!
! Inputs:
!
!
! Dummy variables:
!
!
! Outputs:
!
!
! References:
! [1]
!-----------------------------------------------------------------------
! Define function ------------------------------------------------------
FUNCTION LLF(lambda, y, n)

    ! Define variables -------------------------------------------------
    ! Dimension:
    INTEGER :: n

    ! Data:
    DOUBLE PRECISION :: lambda
    DOUBLE PRECISION, DIMENSION(n) :: y

    ! Dummy:
    DOUBLE PRECISION, DIMENSION(n) :: ylambda
    DOUBLE PRECISION :: ymean
    DOUBLE PRECISION, DIMENSION(n) :: yresid2
    DOUBLE PRECISION :: sigma2
    DOUBLE PRECISION :: jacobian
    DOUBLE PRECISION :: likelihood

    ! Output ::
    DOUBLE PRECISION :: LLF

    ! Calculate log likelihood -----------------------------------------
    ! Transform data:
    IF (lambda .EQ. 0) THEN

        ylambda = LOG(y)

    ELSE

        ylambda = ((y ** lambda) - 1) / lambda

    ENDIF

    ! Calculate mean:
    ymean = SUM(ylambda) / n

    ! Calculate squared residuals:
    yresid2 = (ylambda - ymean) ** 2

    ! Estimate variance:
    sigma2 = SUM(yresid2) / n

    ! Calculate Jacobian:
    jacobian = (lambda - 1) * SUM(LOG(y))

    ! Calculate likelihood:
    likelihood = jacobian - ((n / 2) * LOG(sigma2))

    ! Assign:
    LLF = -likelihood

    ! Return value:
    RETURN

ENDFUNCTION LLF

!-----------------------------------------------------------------------
