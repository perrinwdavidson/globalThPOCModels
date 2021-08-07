!-----------------------------------------------------------------------
!                             normalize0mean
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to normalize and subtract the mean from
! a dataset in preparation for kriging.
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
! Define subroutine ----------------------------------------------------
SUBROUTINE NORMALIZE0MEAN(lambdahat, y, n, ylambda, ymean)

    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Define variables -------------------------------------------------
    ! Dimension:
    INTEGER, INTENT(IN) :: n

    ! Arrays:
    DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: y
    DOUBLE PRECISION, DIMENSION(n), INTENT(OUT) :: ylambda

    ! Scalars:
    DOUBLE PRECISION, INTENT(IN) :: lambdahat
    DOUBLE PRECISION, INTENT(OUT) :: ymean

    ! Define variavles for f2py ----------------------------------------
!f2py INTENT(IN) y
!f2py INTENT(IN) lambdahat
!f2py INTENT(OUT) ylambda
!f2py INTENT(OUT) ymean

    ! Normalize and subtract mean --------------------------------------
    ! Normalize data:
    ylambda = ((y ** lambdahat) - 1) / lambdahat

    ! Calculate mean:
    ymean = SUM(ylambda) / n

    ! Subtract mean:
    ylambda = ylambda - ymean

END SUBROUTINE NORMALIZE0MEAN
!-----------------------------------------------------------------------
