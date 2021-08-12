!-----------------------------------------------------------------------
!                             calcvariogram
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! This program calculates the experimental and fitted variogram for a
! given set of geographic data.
!-----------------------------------------------------------------------
! Define program -------------------------------------------------------
PROGRAM CALCVARIOGRAM

    ! Ensure all variables defined:
    IMPLICIT NONE

    ! Set includes:
    INCLUDE '/usr/local/include/mpif.h'

    ! Define variables -------------------------------------------------
    ! Dimensions:
    INTEGER :: rows
    INTEGER :: columns
    INTEGER :: numdist
    DOUBLE PRECISION :: numbins
    INTEGER :: numsend

    ! MPI:
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER, PARAMETER :: root_process = 0
    INTEGER, PARAMETER :: data_tag0 = 1
    INTEGER, PARAMETER :: data_tag1 = 2
    INTEGER, PARAMETER :: data_tag2 = 3
    INTEGER, PARAMETER :: data_tag3 = 4
    INTEGER, PARAMETER :: data_tag4 = 5
    INTEGER, PARAMETER :: data_tag5 = 6
    INTEGER, PARAMETER :: data_tag6 = 7
    INTEGER :: ierr
    INTEGER :: proc_id
    INTEGER :: num_procs
    INTEGER :: iid
    INTEGER :: sender

    ! Arrays:
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: data_array
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mld
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mldnorm
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lon
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lonrad
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lat
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: latrad
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dist
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: vals
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: binedges
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bins
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bin1
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bin2
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: semivariance
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: semivariancevar

    ! Filename:
    CHARACTER(15) :: filename

    ! Functions:
    DOUBLE PRECISION :: LLF
    DOUBLE PRECISION :: BRENT
    DOUBLE PRECISION :: GREATCIRCLE

    ! Set externals:
    EXTERNAL :: LLF
    EXTERNAL :: BRENT
    EXTERNAL :: GREATCIRCLE

    ! Scalars:
    DOUBLE PRECISION :: ax
    DOUBLE PRECISION :: bx
    DOUBLE PRECISION :: cx
    DOUBLE PRECISION :: tol
    DOUBLE PRECISION :: lambdahat
    DOUBLE PRECISION :: mldmean
    DOUBLE PRECISION :: maxdist
    DOUBLE PRECISION :: lon1
    DOUBLE PRECISION :: lon2
    DOUBLE PRECISION :: lat1
    DOUBLE PRECISION :: lat2
    DOUBLE PRECISION :: semivar
    DOUBLE PRECISION :: semivar_root
    DOUBLE PRECISION :: semivar_child
    DOUBLE PRECISION :: semivarvar
    DOUBLE PRECISION :: semivarvar_root
    DOUBLE PRECISION :: semivarvar_child

    ! Conversion:
    DOUBLE PRECISION, PARAMETER :: deg2rad = 3.14159265359 / 180

    ! Indices:
    INTEGER :: i
    INTEGER :: iter
    INTEGER :: iter_received
    INTEGER :: countbins

    ! MPI initialization -----------------------------------------------
    ! Now replicate this process to create parallel processes. From this
    ! point on, every process executes a separate copy of this program:
    CALL MPI_INIT(ierr)

    ! Find out MY process ID (that core reading this part of the code
    ! after the split above):
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, ierr)

    ! Find out how many processes were started:
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

    ! Ensure only root processes this data:
    IF (proc_id .EQ. root_process) THEN

        ! Print out how many processes started:
        PRINT *, 'There were ', num_procs, ' processes started.'

        ! Open and Read data -------------------------------------------
        ! Set filename:
        filename = 'inputs/mld1.txt'

        ! Get data dimensions:
        CALL GETDIMS(filename, rows, columns)

        ! Allocate data array:
        ALLOCATE(data_array(rows, columns))

        ! Open file:
        OPEN(UNIT=1, FILE=TRIM(filename), STATUS='old', ACTION='read')

        ! Read data:
        DO i = 1, rows, 1

            READ(1, *) data_array(i, :)

        ENDDO

        ! Close file:
        CLOSE(1)

        ! Assign data --------------------------------------------------
        ! Allocate data:
        ALLOCATE(mld(rows))
        ALLOCATE(lon(rows))
        ALLOCATE(lat(rows))

        ! Assign:
        mld = data_array(:, 1)
        lon = data_array(:, 2)
        lat = data_array(:, 3)

        ! Normalize data -----------------------------------------------
        ! Allocate:
        ALLOCATE(mldnorm(rows))

        ! Set Brent minimization parameters:
        ax = -3.0
        bx = 0.0
        cx = 3.0
        tol = 1.48E-8

        ! Determine optimal Box-Cox parameter:
        lambdahat = BRENT(mld, rows, ax, bx, cx, LLF, tol)

        ! Normalize and subtract mean:
        CALL NORMALIZE0MEAN(lambdahat, mld, rows, mldnorm, mldmean)

        ! Make bins ----------------------------------------------------
        ! Calculate dimensions:
        numdist = (rows * (rows - 1)) / 2

        ! Allocate:
        ALLOCATE(lonrad(rows))
        ALLOCATE(latrad(rows))
        ALLOCATE(dist(numdist))
        ALLOCATE(vals(numdist, 2))

        ! Calculate coordinates in radians:
        lonrad = lon * deg2rad
        latrad = lat * deg2rad

        ! Calculate pairwise distances:
        CALL PDIST(lonrad, latrad, mldnorm, rows, dist, vals, numdist, &
                   maxdist, numbins)

        ! Redefine max distance:
        lon1 = 0.0
        lon2 = 0.0
        lat1 = 0.0
        lat2 = 30.0
        maxdist = GREATCIRCLE(lon1, lat1, lon2, lat2)

        ! Allocate more space:
        ALLOCATE(binedges(INT(numbins) + 1))
        ALLOCATE(bins(INT(numbins)))

        ! Make the bins:
        CALL MAKEBINS(INT(numbins), INT(numbins + 1), maxdist, &
                      binedges, bins)

    ! End only the root process doing work:
    ENDIF

    ! Calculate experimental variogram ---------------------------------
    ! Initialize bin counter:
    countbins = 0

    ! Count bins:
    DO

        ! Split up root and children processes:
        IF (proc_id .EQ. root_process) THEN

            ! Initialization:
            IF (countbins .EQ. 0) THEN

                ! Allocate space:
                ALLOCATE(semivariance(INT(6))) ! numbins
                ALLOCATE(semivariancevar(INT(6))) ! numbins

                ! Set counter:
                iter = 1

            ENDIF

            ! Send data -----------------------------------
            DO iid = 1, (num_procs - 1), 1

                ! Bin data:
                bin1 = PACK(vals(:, 1), dist .GE. binedges(iter) &
                            .AND. dist .LE. binedges(iter + 1))
                bin2 = PACK(vals(:, 2), dist .GE. binedges(iter) &
                            .AND. dist .LE. binedges(iter + 1))
                numsend = size(bin1)

                ! Send to processes:
                CALL MPI_SEND(iter, 1, MPI_INTEGER, iid, &
                              data_tag5, MPI_COMM_WORLD, ierr)
                CALL MPI_SEND(numsend, 1, MPI_INTEGER, iid, &
                              data_tag0, MPI_COMM_WORLD, ierr)
                CALL MPI_SEND(bin1, numsend, MPI_DOUBLE_PRECISION, &
                              iid, data_tag1, MPI_COMM_WORLD, ierr)
                CALL MPI_SEND(bin2, numsend, MPI_DOUBLE_PRECISION, &
                              iid, data_tag2, MPI_COMM_WORLD, ierr)

                ! Iterate through bins:
                iter = iter + 1

            ENDDO

            ! Calculate own data --------------------------
            ! Bin data:
            bin1 = PACK(vals(:, 1), dist .GE. binedges(iter) &
                        .AND. dist .LE. binedges(iter + 1))
            bin2 = PACK(vals(:, 2), dist .GE. binedges(iter) &
                        .AND. dist .LE. binedges(iter + 1))
            numsend = size(bin1)

            ! Calculate semivariance:
            CALL CALCEXPVARIO(bin1, bin2, numsend, semivar_root, &
                              semivarvar_root)

            ! Store data:
            semivariance(iter) = semivar_root
            semivariancevar(iter) = semivarvar_root

            ! Print out -----------------------------------
            PRINT *, 'Process: ', 0, ' done with bin ', iter, &
                     '| Results: ', semivar_root, semivarvar_root

            ! Iterate through bins ------------------------
            iter = iter + 1

            ! Recieve calculations ------------------------
            DO iid = 1, (num_procs - 1), 1

                ! Recieve number_input:
                CALL MPI_RECV(iter_received, 1, MPI_INTEGER, &
                              MPI_ANY_SOURCE, data_tag6, &
                              MPI_COMM_WORLD, status, ierr)
                CALL MPI_RECV(semivar, 1, MPI_DOUBLE_PRECISION, &
                              MPI_ANY_SOURCE, data_tag3, &
                              MPI_COMM_WORLD, status, ierr)
                CALL MPI_RECV(semivarvar, 1, MPI_DOUBLE_PRECISION, &
                              MPI_ANY_SOURCE, data_tag4, &
                              MPI_COMM_WORLD, status, ierr)

                ! Get ID:
                sender = status(MPI_SOURCE)

                ! Store data:
                semivariance(iter_received) = semivar
                semivariancevar(iter_received) = semivarvar

                ! Print out:
                PRINT *, 'Process: ', sender, ' done with bin ', &
                         iter_received, '| Results: ', semivar, &
                         semivarvar

            ENDDO

        ELSE

            ! Receive data --------------------------------
            ! Get dimension:
            CALL MPI_RECV(iter, 1, MPI_INTEGER, root_process, &
                          data_tag5, MPI_COMM_WORLD, status, ierr)
            CALL MPI_RECV(numsend, 1, MPI_INTEGER, root_process, &
                          data_tag0, MPI_COMM_WORLD, status, ierr)

            ! Allocate space:
            ALLOCATE(bin1(numsend))
            ALLOCATE(bin2(numsend))

            ! Get arrays:
            CALL MPI_RECV(bin1, numsend, MPI_DOUBLE_PRECISION, &
                          root_process, data_tag1, MPI_COMM_WORLD, &
                          status, ierr)
            CALL MPI_RECV(bin2, numsend, MPI_DOUBLE_PRECISION, &
                          root_process, data_tag2, MPI_COMM_WORLD, &
                          status, ierr)

            ! Calculate own data --------------------------
            ! Calculate:
            CALL CALCEXPVARIO(bin1, bin2, numsend, semivar, semivarvar)

            ! Get values:
            semivar_child = semivar
            semivarvar_child = semivarvar

            ! Send data -----------------------------------
            CALL MPI_SEND(iter, 1, MPI_INTEGER, root_process, &
                          data_tag6, MPI_COMM_WORLD, ierr)
            CALL MPI_SEND(semivar_child, 1, MPI_DOUBLE_PRECISION, &
                          root_process, data_tag3, MPI_COMM_WORLD, ierr)
            CALL MPI_SEND(semivarvar_child, 1, MPI_DOUBLE_PRECISION, &
                          root_process, data_tag4, MPI_COMM_WORLD, ierr)

            ! Deallocate:
            DEALLOCATE(bin1)
            DEALLOCATE(bin2)

        ENDIF

        ! Wait for all to catch up ------------------------
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        ! Iterate through:
        countbins = countbins + num_procs

        ! Move on if done:
        IF (countbins .GE. 6) THEN ! numbins

            IF (proc_id .EQ. 0) THEN

                PRINT*, 'Done with experimental variogram calculations.'

            ENDIF

            GO TO 10

        ENDIF

    ENDDO

    ! Continue with the program:
10  CONTINUE

    ! Stop program -----------------------------------------------------
    ! Stop all processes doing these calculations:
    CALL MPI_FINALIZE(ierr)

    ! Stop program:
    STOP

ENDPROGRAM

!-----------------------------------------------------------------------
!                                getdims
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! This program calculates the dimenions of a file to read.
!-----------------------------------------------------------------------
! Define program -------------------------------------------------------
SUBROUTINE GETDIMS(filename, rows, columns)

    ! Variables -------------------------------------------
    ! Inputs:
    CHARACTER(15), INTENT(IN) :: filename

    ! Dummy:
    INTEGER :: i
    INTEGER :: io
    CHARACTER(128) :: buffer
    INTEGER :: string_length

    ! Outputs:
    INTEGER, INTENT(OUT) :: rows
    INTEGER, INTENT(OUT) :: columns

    ! Open text file --------------------------------------
    OPEN(UNIT=1, FILE=TRIM(filename), STATUS='old', ACTION='read')

    ! Count the number of columns -------------------------
    ! Read first line with spaces included:
    READ(1, '(a)') buffer

    ! Go back to the file beginning:
    REWIND(1)

    ! Find the length of the string without ending spaces:
    string_length = len(buffer)
    DO WHILE (buffer(string_length:string_length) == ' ')

        string_length = string_length - 1

    ENDDO

    ! Count the number of columns in the first line:
    columns = 0
    DO i = 0, string_length, 1

        IF (buffer(i:i) == ' ') THEN

            columns = columns + 1

        ENDIF

    ENDDO
    columns = columns + 1  ! there is one more column than space, as removed the last one above

    ! Count the number of rows ----------------------------
    ! Count the number of lines in the file:
    rows = 0
    DO

        ! Read data:
        READ(1, *, IOSTAT=io)

        ! Exit if no more rows:
        IF (io .NE. 0) THEN

            EXIT

        ENDIF

        ! Count rows:
        rows = rows + 1

    ENDDO

    ! Close file ------------------------------------------
    CLOSE(1)

ENDSUBROUTINE GETDIMS

!-----------------------------------------------------------------------
!                               brent
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to minimize a log likelihood function.
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
!                             normalize0mean
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to normalize and subtract the mean from
! a dataset in preparation for kriging.
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

    ! Normalize and subtract mean --------------------------------------
    ! Normalize data:
    ylambda = ((y ** lambdahat) - 1) / lambdahat

    ! Calculate mean:
    ymean = SUM(ylambda) / n

    ! Subtract mean:
    ylambda = ylambda - ymean

END SUBROUTINE NORMALIZE0MEAN

!-----------------------------------------------------------------------
!                                pdist
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to calculate the pairwise distances
! between coordinates using the great circle distance.
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
!                                makebins
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to make bins for an experimental
! variogram estimation.
!-----------------------------------------------------------------------
! Define subroutine ----------------------------------------------------
SUBROUTINE MAKEBINS(numbins, numbins1, maxdist, binedges, bins)

    ! Ensure variables defined:
    IMPLICIT NONE

    ! Define variables -------------------------------------------------
    ! Dimension:
    INTEGER, INTENT(IN):: numbins
    INTEGER, INTENT(IN):: numbins1

    ! Arrays:
    DOUBLE PRECISION, DIMENSION(numbins1), INTENT(OUT) :: binedges
    DOUBLE PRECISION, DIMENSION(numbins), INTENT(OUT) :: bins

    ! Scalars:
    DOUBLE PRECISION, INTENT(IN) :: maxdist

    ! Indices:
    INTEGER :: i

    ! Make bins --------------------------------------------------------
    ! Make numbins with equal spacing:
    DO i = 1, numbins1, 1

        binedges(i) = 0 + (maxdist * (i - 1) / (numbins1 - 1))

    ENDDO

    ! Find center of bins:
    bins = ((binedges(2:numbins1) - binedges(1:numbins)) / 2) &
           + binedges(1:numbins)

ENDSUBROUTINE MAKEBINS

!-----------------------------------------------------------------------
!                              calcexpvario
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to estimate the experimental
! variogram of a given dataset.
!-----------------------------------------------------------------------
! Define subroutine ----------------------------------------------------
SUBROUTINE CALCEXPVARIO(bin1, bin2, k, semivar, semivarvar)

    ! Ensure variables defined:
    IMPLICIT NONE

    ! Define variables -------------------------------------------------
    ! Dimensions:
    INTEGER, INTENT(IN) :: k

    ! Inputs:
    DOUBLE PRECISION, DIMENSION(k), INTENT(IN) :: bin1
    DOUBLE PRECISION, DIMENSION(k), INTENT(IN) :: bin2

    ! Variance calculation:
    DOUBLE PRECISION :: CALCVAR
    EXTERNAL :: CALCVAR

    ! Scalars:
    DOUBLE PRECISION :: lagvalssum
    DOUBLE PRECISION :: normcorrect

    ! Outputs:
    DOUBLE PRECISION, INTENT(OUT) :: semivar
    DOUBLE PRECISION, INTENT(OUT) :: semivarvar

    ! Calculate semivariance ------------------------------
    ! Calculate square-rooted distances:
    lagvalssum = SUM(SQRT(ABS(bin1 - bin2)))

    ! Calculate Cressie estimator denominator:
    normcorrect = 0.457 + (0.494 / k) + (0.045 / (k ** 2))

    ! Calculate semivariance for bin:
    semivar = (((lagvalssum / k) ** 4) / 2) / normcorrect

    ! Estimate variance of semivariance -------------------
    semivarvar = CALCVAR(bin1, bin2, semivar, k)

ENDSUBROUTINE CALCEXPVARIO

!-----------------------------------------------------------------------
!                                calcvar
!-----------------------------------------------------------------------
! Author: Perrin W. Davidson
! Contact: perrinwdavidson@gmail.com
! Last updated: 05/08/2021
!-----------------------------------------------------------------------
! About:
! The purpose of this function is to calculate the variance of the
! semivariance estimate.
!-----------------------------------------------------------------------
! Define function ------------------------------------------------------
FUNCTION CALCVAR(bin1, bin2, semivar, n)

    ! Define variables -------------------------------------------------
    ! Dimension:
    INTEGER :: n

    ! Data:
    DOUBLE PRECISION, DIMENSION(n) :: bin1
    DOUBLE PRECISION, DIMENSION(n) :: bin2
    DOUBLE PRECISION :: semivar

    ! Dummy:
    DOUBLE PRECISION :: C
    DOUBLE PRECISION :: zi
    DOUBLE PRECISION :: zi1
    DOUBLE PRECISION :: zj
    DOUBLE PRECISION :: zj1

    ! Indices:
    INTEGER :: i
    INTEGER :: j

    ! Output:
    DOUBLE PRECISION :: CALCVAR

    ! Calculate variance ----------------------------------
    ! Initialize:
    C = 0

    ! Loop:
    DO i = 1, n, 1

        ! Get values:
        zi = bin1(i)
        zi1 = bin2(i)

        DO j = 1, n, 1

            ! Get values:
            zj = bin1(j)
            zj1 = bin2(j)

            ! Calculate covariance:
            C = C + ((((zi - zi1) ** 2) * ((zj - zj1) ** 2)) &
                - (semivar ** 2))

        ENDDO

    ENDDO

    ! Calculate variance:
    CALCVAR = C / (2 * (n ** 2))

    ! Return value:
    RETURN

ENDFUNCTION CALCVAR

!-----------------------------------------------------------------------
