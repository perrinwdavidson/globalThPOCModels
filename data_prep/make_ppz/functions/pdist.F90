!-----------------------------------------------------------------------
!                      calculate pairwise distances
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to calculate the pairwise distances
! between coordinates using the great circle distance.
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
SUBROUTINE PDIST(x, y, z, n, dist, vals, m, maxdist, numbins)

    ! Define variables -------------------------------------------------
    ! Dimension:
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: m

    ! Distance calculation:
    DOUBLE PRECISION :: GREATCIRCLE
    EXTERNAL :: GREATCIRCLE

    ! Arrays:
    DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: x
    DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: y
    DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: z
    DOUBLE PRECISION, DIMENSION(m), INTENT(OUT) :: dist
    DOUBLE PRECISION, DIMENSION(m, 2), INTENT(OUT) :: vals

    ! Scalars:
    DOUBLE PRECISION, INTENT(OUT) :: maxdist
    DOUBLE PRECISION, INTENT(OUT) :: numbins

    ! Indices:
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k

    ! Define variables for f2py ----------------------------------------
!f2py INTENT(IN) n
!f2py INTENT(IN) m
!f2py INTENT(IN) x
!f2py INTENT(IN) y
!f2py INTENT(IN) z
!f2py INTENT(OUT) dist
!f2py INTENT(OUT) vals
!f2py INTENT(OUT) maxdist
!f2py INTENT(OUT) numbins

    ! Calculate distance -----------------------------------------------
    ! Set distance counter:
    k = 1

    ! Loop through all points:
    DO i = 1, n, 1

        DO j = i + 1, n, 1

            ! Calculate distance:
            dist(k) = GREATCIRCLE(x(i), y(i), x(j), y(j))

            ! Get coorespinding values:
            vals(k, 1) = z(i)
            vals(k, 2) = z(j)

            ! Iterate:
            k = k + 1

        ENDDO

    ENDDO

    ! Calculate parameters ---------------------------------------------
    ! Calculate diameter:
    maxdist = MAXVAL(dist)

    ! Calculate number of bins with Sturges Rule:
    numbins = INT(CEILING(2 * (LOG(REAL(m)) / LOG(2.0)) + 1))

ENDSUBROUTINE PDIST

!-----------------------------------------------------------------------
!                    great circle distance - scalar
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to calculate the great circle distance
! between two points.
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
FUNCTION GREATCIRCLE(lon1, lat1, lon2, lat2)

    ! Define variables -------------------------------------------------
    ! coordinates:
    DOUBLE PRECISION :: lon1
    DOUBLE PRECISION :: lat1
    DOUBLE PRECISION :: lon2
    DOUBLE PRECISION :: lat2

    ! Distance:
    DOUBLE PRECISION :: GREATCIRCLE

    ! Scalars:
    DOUBLE PRECISION, PARAMETER :: r = 6371.0088

    ! Calculate distance -----------------------------------------------
    GREATCIRCLE = (2 * r) * ASIN(SQRT((SIN((lat1 - lat2) / 2) ** 2) &
                                 + (COS(lat1) * COS(lat2)           &
                                 * (SIN((lon1 - lon2) / 2) ** 2))))

    ! Return values:
    RETURN

ENDFUNCTION GREATCIRCLE
!-----------------------------------------------------------------------
