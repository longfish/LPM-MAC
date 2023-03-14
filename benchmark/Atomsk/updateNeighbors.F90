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
DOUBLE PRECISION:: Box(3,2), H_matrix(9), Radius, start, finish, offset, threshold
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: Positions

!Dimensionality of the problem, 3D
ndim = 3
!The dimension of the physical domain using the format:
!3D: [Xmin, Ymin, Zmin, Xmax, Ymax, Zmax]
Box = RESHAPE((/ 0.D0, 0.D0, 0.D0, 1D0, 1D0, 1D0 /), (/ 3, 2 /))

Radius = 1.43

!model configuration parameters
offset = 1.5*Radius
threshold = 1.3*Radius

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

INTEGER:: i, j
DOUBLE PRECISION:: s1, s2, s3

!Read in the Voronoi Cell seeds data
OPEN(UNIT=13, FILE="final.cfg", STATUS='old', ACTION='read')
READ(13, '(A)') line
READ(line,*) temp, temp, temp, temp, nparticle

ALLOCATE(Positions( nparticle,ndim))
ALLOCATE(neighbors(nparticle,nneighbors))
neighbors = 0
ALLOCATE(layer(nparticle, nneighbors))
layer = 0
ALLOCATE(grainsign(nparticle))
grainsign = 0

DO i = 1, 2
    READ(13, '(A)') line
ENDDO

DO i = 1, 9
    READ(13, '(A)') line
    do j = 1, len_trim( line )
        if ( line(j:j) == "(" .or. line(j:j) == ")" .or. line(j:j) == "," ) line(j:j) = " "
    enddo

    READ(line,*) temp, temp, temp, temp, H_matrix(i)
ENDDO

Box(1,2) = H_matrix(1)
Box(2,2) = H_matrix(5)
Box(3,2) = H_matrix(9)

DO i = 1, 5
    READ(13, '(A)') line
ENDDO

DO i = 1, nparticle
    READ(13,*) s1, s2, s3, grainsign(i)
    Positions(i,1) = (s1*H_matrix(1)+s2*H_matrix(4)+s3*H_matrix(7))
    Positions(i,2) = (s1*H_matrix(2)+s2*H_matrix(5)+s3*H_matrix(8))
    Positions(i,3) = (s1*H_matrix(3)+s2*H_matrix(6)+s3*H_matrix(9))
    ! Positions(i,1) = s1
    ! Positions(i,2) = s2
    ! Positions(i,3) = s3
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
    WRITE(13, '(2F8.1)') Box(3, 1), Box(3, 2)
    WRITE(13, '(A)') 'ITEM: ATOMS id type x y z '
    DO t = 1, nparticle
        WRITE(13, '(I0, I6, 3F8.4)') t, grainsign(t), Positions(t, 1), Positions(t, 2), Positions(t, 3)
    ENDDO
    CLOSE(UNIT = 13)

    
    END SUBROUTINE Output
    
    !==============================================END Output===========================================================
    
END PROGRAM Main
