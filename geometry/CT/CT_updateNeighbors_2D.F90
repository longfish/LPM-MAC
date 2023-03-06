INCLUDE 'mkl_pardiso.f90'
PROGRAM MAIN
USE mkl_pardiso
USE OMP_lib
USE IFPORT
IMPLICIT NONE

CHARACTER(LEN=100):: line, temp
INTEGER:: ndim, i, t, nparticle, nbond, nneighbors, ngrain
INTEGER, ALLOCATABLE, DIMENSION(:,:):: bonds, layer, neighbors
INTEGER, ALLOCATABLE, DIMENSION(:):: grainsign
DOUBLE PRECISION, PARAMETER :: PI = 3.1415926
DOUBLE PRECISION:: Box(2,2), H_matrix(9), Radius, start, finish, offset, threshold
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: Positions

!Dimensionality of the problem, 3D
ndim = 2
!The dimension of the physical domain using the format:
!3D: [Xmin, Ymin, Xmax, Ymax]
Box = RESHAPE((/ 0.D0, 0.D0, 40D0, 40D0 /), (/ 2, 2 /))

particles_per_row = 31

!Number of neighbors
nneighbors = 36

!record the CPU time
start = OMP_GET_WTIME()
!----------------------------------------------------------------------------------------------------------
!initialize the particle position and allocate memory for the global matrices
!----------------------------------------------------------------------------------------------------------
CALL  Initialization

!search neighbors and write bond information
nbond = 0
CALL  Searchneighbor

CALL  Output

CALL MKL_FREE_BUFFERS

!record the CPU time
finish = OMP_GET_WTIME()
WRITE(*,*) 'Total computation time: ', finish-start, 'seconds.'

!SUBROUTINES
CONTAINS

!==================================================SUBROUTINE Initialization============================================

SUBROUTINE Initialization
IMPLICIT NONE

DOUBLE PRECISION:: delta, x, y
INTEGER:: i, j, n

delta = ( Box(1,2) - Box(1,1) ) / particles_per_row
Radius = delta / 2.0
rows = FLOOR ( ( Box(2,2) - Box(2,1) ) / delta )
nparticle = particles_per_row*rows
ALLOCATE(Positions( nparticle,ndim))
Positions = 0.

n = 0
DO j = 1, rows
    y = Box(2,1) + delta * ( j - 1 ) + Radius
    DO i = 1, particles_per_row
        x = Box(1,1) + delta * ( i - 1 ) + Radius
        n = n + 1
        IF  (n .LE. nparticle) THEN
            Positions(n,1) = x
            Positions(n,2) = y
        ENDIF
    ENDDO
ENDDO

CLOSE(UNIT=13)
WRITE(*,*) 'Number of particles:', nparticle


END SUBROUTINE Initialization

!================================================END Initialization======================================================

!===============================================SUBROUTINE Searchneighbor============================================

SUBROUTINE  Searchneighbor
IMPLICIT NONE

DOUBLE PRECISION:: dis
INTEGER:: i, j, index1

!finding the first and second neighbors for each particle
!$OMP PARALLEL DO PRIVATE(index1, dis) &
!$OMP& REDUCTION(+:nbond) 
DO i = 1, nparticle
    index1 = 0
    DO j = 1, nparticle
        dis = SQRT((Positions(j,1) - Positions(i,1))**2 + (Positions(j,2) - Positions(i,2))**2 + (Positions(j,3) - Positions(i,3))**2)
        IF(dis  .LT. 2.01*Radius .AND. j .NE. i) THEN !the first nearest neighbors
            index1 = index1 + 1
            neighbors(i,index1) = j
            layer(i,index1) = 1
        ELSEIF(dis .GT. 2.01*Radius .AND. dis .LT. 2.01*SQRT(2.)*Radius)THEN !the second nearest neighbors            
            index1 = index1 + 1
            neighbors(i,index1) = j
            layer(i,index1) = 2
        ENDIF
    ENDDO
    nbond = nbond + index1
ENDDO
!$OMP END PARALLEL DO

WRITE(*,*) 'Number of bonds:', nbond

ALLOCATE(bonds(nbond, 3))
bonds = 0

END SUBROUTINE Searchneighbor
!=============================================END Searchneighbor=====================================================

!==============================================SUBROUTINE Output=====================================================

SUBROUTINE Output
    IMPLICIT NONE
    
    INTEGER:: t
    
    OPEN(UNIT = 13, FILE = 'polycrystal.dump', ACTION = 'WRITE', STATUS = 'REPLACE', POSITION = 'APPEND')
    WRITE(13, '(A)') 'ITEM: TIMESTEP'
    WRITE(13, '(I0)') 0
    WRITE(13, '(A)') 'ITEM: NUMBER OF ATOMS'
    WRITE(13, '(I0)') nparticle
    WRITE(13, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
    WRITE(13, '(2F8.1)') Box(1, 1), Box(1, 2) ! unit is mm
    WRITE(13, '(2F8.1)') Box(2, 1), Box(2, 2)
    WRITE(13, '(2F8.1)') 0, 0
    WRITE(13, '(A)') 'ITEM: ATOMS id type x y z '
    DO t = 1, nparticle
        WRITE(13, '(I0, I6, 3F8.4)') t, grainsign(t), Positions(t, 1), Positions(t, 2), 0
    ENDDO
    CLOSE(UNIT = 13)

    
    END SUBROUTINE Output
    
    !==============================================END Output===========================================================
    
END PROGRAM Main
