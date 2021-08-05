!-----------------------------------------------------------------------
!                             interp2d
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 04/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to interpolate regular, gridded data
! to either a secondary grid with different coordinates or an irregular
! set of points. The method by which we are doing this is bilinear
! interpolation.
!
! Derivation:
! We start with the 1D case, remembering that this is simply the 1D case
! implemented twice. The 1D case is derived by noting that for a set of
! known, gridded coordinates (x1, y1) and (x2, y2) and the known position
! xk for which we are trying to estimate the corresponding yk that the
! following relation holds:
!       (yk - y1)/(xk - x1) = (y2 - yk)/(x2 - xk).                   (1)
! The reason for this is that we assume that (xq, yq) falls on the line
! between (x1, y1) and (x2, y2). This line is the called the interpolant,
! or specifically the 'linear interpolant'. Therefore, the slopes are
! the same!
!   Now, consider this for two x lines at two different y values
! as is the case in bilinear interpolation. Then, you can solve this
! set of equations, first in the x direction and then the y (x or y
! first doesn't matter), and get the equation that we use in this
! subroutine:
!       zk = (1 - rx)*(1 - ry)*zij + rx*(1 - ry)*zi1j
!          + (1 - rx)*ry*zij1 + rx*ry*zi1j1,
! where
!       rx = (xk - xi)/(xi1 - xi), and
!       ry = (yk - yj)/(yj1 - yj)
! Also note that this is an unbiased method, as we can view the rx and ry
! as weights that sum to 1.
!
! Inputs:
!   x - x sample coordinates vect0r of dimension (nx)
!   y - y sample coordinates vect0r of dimension (ny)
!   z - sample value array of dimensions (nx, ny)
!   xq - x query coordinate vector of dimension (nq)
!   yq - y query coordinate vector of dimension (nq)
!   nx - dimension of x sample vector
!   ny - dimension of y sample vector
!   nq - dimension of query vectors
!   fillvalue - the fill value given to extrapolated data
!
! Dummy variables:
!   minx - the minimum value of the x sample coordinate vector
!   miny - the minimum value of the y sample coordinate vector
!   maxx - the maximum value of the x sample coordinate vector
!   maxy - the maximum value of the y sample coordinate vector
!   xk - instance of query x coordinate
!   yk - instance of query y coordinate
!   xi - instance of sample x coordinate at index (i)
!   xi1 - instance of sample x coordinate at index (i + 1)
!   yj - instance of sample y coordinate at index (j)
!   yj1 - instance of sample y coordinate at index (j + 1)
!   zij - instance of sample z value at index (i, j)
!   zi1j - instance of sample z value at index (i + 1, j)
!   zij1 - instance of sample z value at index (i, j + 1)
!   zi1j1 - instance of sample z value at index (i + 1, j + 1)
!   rx - interpolation weights in x direction
!   ry - interpolation weights in y direction
!   i - x coordinate index
!   j - y coordinate index
!   k - query coordinate index
!   sx - success indicator of finding sample x coordinates
!   sy - success indicator of finding sample y coordinates
!
! Outputs:
!   zq - estimated z values at the query coordinates of for each pair (xk, yk) with dimension (nq)
!
! References:
! [1] Glover, Jenkins, and Doney (2011): 'Modeling Methods for Marine Science'
!-----------------------------------------------------------------------
! Define subroutine ----------------------------------------------------
SUBROUTINE interp2d(x, y, z, xq, yq, zq, nx, ny, nq, fillvalue)

    ! Ensure all variables declared:
    IMPLICIT NONE

    ! Define variables -------------------------------------------------
    ! Dimensions:
    INTEGER, INTENT(IN) :: nx
    INTEGER, INTENT(IN) :: ny
    INTEGER , INTENT(IN):: nq

    ! Arrays:
    DOUBLE PRECISION, DIMENSION(nx), INTENT(IN) :: x
    DOUBLE PRECISION, DIMENSION(ny), INTENT(IN) :: y
    DOUBLE PRECISION, DIMENSION(nx, ny), INTENT(IN) :: z
    DOUBLE PRECISION, DIMENSION(nq), INTENT(IN) :: xq
    DOUBLE PRECISION, DIMENSION(nq), INTENT(IN) :: yq
    DOUBLE PRECISION, DIMENSION(nq), INTENT(OUT) :: zq

    ! Fill value:
    INTEGER , INTENT(IN):: fillvalue

    ! Dummy:
    DOUBLE PRECISION :: minx
    DOUBLE PRECISION :: miny
    DOUBLE PRECISION :: maxx
    DOUBLE PRECISION :: maxy
    DOUBLE PRECISION :: xk
    DOUBLE PRECISION :: yk
    DOUBLE PRECISION :: xi
    DOUBLE PRECISION :: xi1
    DOUBLE PRECISION :: yj
    DOUBLE PRECISION :: yj1
    DOUBLE PRECISION :: zij
    DOUBLE PRECISION :: zi1j
    DOUBLE PRECISION :: zij1
    DOUBLE PRECISION :: zi1j1
    DOUBLE PRECISION :: rx
    DOUBLE PRECISION :: ry

    ! Indices:
    INTEGER :: i
    INTEGER :: j
    INTEGER :: k
    INTEGER :: sx
    INTEGER :: sy

    ! Define variables for f2py ----------------------------------------
    ! Inputs:
!f2py intent(in) nx
!f2py intent(in) ny
!f2py intent(in) nq
!f2py intent(in) x
!f2py intent(in) y
!f2py intent(in) z
!f2py intent(in) xq
!f2py intent(in) yq
!f2py intent(in) fillvalue

    ! Outputs:
!f2py intent(out) zq

    ! Interpolate ------------------------------------------------------
    ! Preliminarily calculate mins and maxes of arrays:
    minx = MINVAL(x)
    miny = MINVAL(y)
    maxx = MAXVAL(x)
    maxy = MAXVAL(y)

    ! Initialize values:
    sx = 0
    sy = 0
    xi = 0
    xi1 = 0
    yj = 0
    yj1 = 0

    ! Loop through all query locations:
    DO k = 1, nq, 1

        ! Find query points -------------------------------
        ! Get data:
        xk = xq(k)
        yk = yq(k)

        ! Skip and fill if extrapolation in x:
        IF (xk .LT. minx .OR. xk .GT. maxx) THEN

            zq(k) = fillvalue
            CYCLE

        ENDIF

        ! Skip and fill if extrapolation in y:
        IF (yk .LT. miny .OR. yk .GT. maxy) THEN

            zq(k) = fillvalue
            CYCLE

        ENDIF

        ! Find corresponding sample x values --------------
        ! Note that the index that results from this is that for
        ! x(i + 1).
        ! Loop through all possible x values:
        i = 0
        DO WHILE (i .LE. (nx - 1))

            ! Calculate index:
            i = i + 1

            ! Find x values:
            xi = x(i)
            xi1 = x(i + 1)

            ! If bounds query value, exit loop:
            IF (xi .LE. xk .AND. xi1 .GE. xk) THEN

                ! Let me know it was a success:
                sx = 1

                ! Move on:
                GOTO 10

            ENDIF

        ENDDO

        ! Fill and cycle if unsuccessful:
        IF (sx .NE. 1) THEN

            zq(k) = fillvalue
            CYCLE

        ELSE

            sx = 0

        ENDIF

        ! Success - let's move on:
10      CONTINUE

        ! Find corresponding sample y values --------------
        ! Note that the index that results from this is that for
        ! y(j + 1).
        ! Loop through all possible x values:
        j = 0
        DO WHILE (j .LE. (ny - 1))

            ! Calculate index:
            j = j + 1

            ! Find x values:
            yj = y(j)
            yj1 = y(j + 1)

            ! If bounds query value, exit loop:
            IF (yj .LE. yk .AND. yj1 .GE. yk) THEN

                ! Let me know it was a success:
                sy = 1

                ! Move on:
                GOTO 20

            ENDIF

        ENDDO

        ! Fill and cycle if unsuccessful:
        IF (sy .NE. 1) THEN

            zq(k) = fillvalue
            CYCLE

        ELSE

            sy = 0

        ENDIF

        ! Success - let's move on:
20      CONTINUE

        ! Perform calculations ----------------------------
        ! Find correspinding sample z values:
        zij = z(i, j)
        zi1j = z(i + 1, j)
        zij1 = z(i, j + 1)
        zi1j1 = z(i + 1, j + 1)

        ! Calculate unbaised, weighted coefficients:
        rx = (xk - xi) / (xi1 - xi)
        ry = (yk - yj) / (yj1 - yj)

        ! Interpolate values:
        zq(k) = ((1 - rx) * (1 - ry) * zij) &
                + (rx * (1 - ry) * zi1j)    &
                + ((1 - rx) * ry * zij1)    &
                + (rx * ry * zi1j1)

    ENDDO

END SUBROUTINE interp2d
!-----------------------------------------------------------------------
