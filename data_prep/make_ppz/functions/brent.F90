!-----------------------------------------------------------------------
!                               brent
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to minimize a log likelihood function.
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
FUNCTION BRENT(data, n, ax, bx, cx, f, tol)

    ! Define variables -------------------------------------------------
    ! Numerical routine parameters:
    INTEGER, PARAMETER :: ITMAX = 100
    DOUBLE PRECISION, PARAMETER :: CGOLD = 0.3819660
    DOUBLE PRECISION, PARAMETER :: ZEPS = 1.0E-10

    ! Bracket:
    DOUBLE PRECISION :: ax
    DOUBLE PRECISION :: bx
    DOUBLE PRECISION :: cx

    ! Data input:
    INTEGER :: n
    DOUBLE PRECISION, DIMENSION(n) :: data

    ! Precision of estimate:
    DOUBLE PRECISION :: tol

    ! Function to minimize:
    DOUBLE PRECISION :: f
    EXTERNAL :: f

    ! Output:
    DOUBLE PRECISION :: xmin
    DOUBLE PRECISION :: BRENT

    ! Indices:
    INTEGER :: iter

    ! Dummy:
    DOUBLE PRECISION :: a
    DOUBLE PRECISION :: b
    DOUBLE PRECISION :: d
    DOUBLE PRECISION :: e
    DOUBLE PRECISION :: etemp
    DOUBLE PRECISION :: fu
    DOUBLE PRECISION :: fv
    DOUBLE PRECISION :: fw
    DOUBLE PRECISION :: fx
    DOUBLE PRECISION :: p
    DOUBLE PRECISION :: q
    DOUBLE PRECISION :: r
    DOUBLE PRECISION :: tol1
    DOUBLE PRECISION :: tol2
    DOUBLE PRECISION :: u
    DOUBLE PRECISION :: v
    DOUBLE PRECISION :: w
    DOUBLE PRECISION :: x
    DOUBLE PRECISION :: xm

    ! Minimize function ------------------------------------------------
    ! Initialize variables and functions:
    a = MIN(ax, cx)
    b = MAX(ax, cx)
    v = bx
    w = v
    x = v
    e = 0.
    fx = f(x, data, n)
    fv = fx
    fw = fx

    ! Start main loop:
    DO iter = 1, ITMAX, 1

        xm = 0.5 * (a + b)
        tol1 = (tol * ABS(x)) + ZEPS
        tol2 = 2. * tol1

        IF (ABS(x - xm) .LE. (tol2 - (.5 * (b - a)))) THEN

             GOTO 3

        ENDIF

        IF (ABS(e) .GT. tol1) THEN

            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * (q - ((x - w) * r))
            q = 2. * (q - r)

            IF (q .GT. 0.) THEN

                p = -p

            ENDIF

            q = ABS(q)
            etemp = e
            e = d

            IF (ABS(p) .GE. ABS(.5 * q * etemp) &
               .OR. p .LE. (q * (a - x))        &
               .OR. p .GE. (q * (b - x))) THEN

               GOTO 1

            ENDIF

            d = p / q
            u = x + d

            IF ((u - a) .LT. tol2 .OR. (b - u) .LT. tol2) THEN

                d = SIGN(tol1, (xm - x))

            ENDIF

            GOTO 2

        ENDIF

1       IF (x .GE. xm) THEN

            e = a - x

        ELSE

            e = b - x

        ENDIF

        d = CGOLD * e

2       IF (ABS(d) .GE. tol1) THEN

            u = x + d

        ELSE

            u = x + SIGN(tol1, d)

        ENDIF

        fu = f(u, data, n)

        IF (fu .LE. fx) THEN

            IF (u .GE. x) THEN

                a = x

            ELSE

                b = x

            ENDIF

            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu

        ELSE

            IF (u .LT. x) THEN

                a = u

            ELSE

                b = u

            ENDIF

            IF (fu .LE. fw .OR. w .EQ. x) THEN

                v = w
                fv = fw
                w = u
                fw = fu

            ELSEIF (fu .LE. fv .OR. v .EQ. x .OR. v .EQ. w) THEN

                v = u
                fv = fu

            ENDIF

        ENDIF

    ENDDO

    PRINT*, '*-- Brent minimization exceeded maximum iterations --*'

3   xmin = x
    BRENT = xmin  ! could be equal to fx

    RETURN

ENDFUNCTION BRENT

!-----------------------------------------------------------------------
