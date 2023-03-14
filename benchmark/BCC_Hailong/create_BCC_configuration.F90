INCLUDE 'mkl_pardiso.f90'
PROGRAM MAIN
USE mkl_pardiso
USE OMP_lib
USE IFPORT
IMPLICIT NONE

INTEGER:: ndim, i, t, nparticle, nneighbors, nneighbors1, nneighbors2, nneighbors_afem, nneighbors_afem1, nneighbors_afem2, ngrain
INTEGER, ALLOCATABLE, DIMENSION(:,:):: K_pointer, nsign, neighbors, neighbors1, neighbors2, Conn, neighbors_afem, neighbors_afem1, neighbors_afem2, Sneighbor
INTEGER, ALLOCATABLE, DIMENSION(:):: IK, JK, Front, Back, Left, Right, Top, Bottom, grainsign
DOUBLE PRECISION:: Box(3,2), C(3), Mapping(3,3), KnTv(3), U_stepx, U_stepy, U_stepz, Kn1, Kn2, Tv, start, finish, offset, threshold, SRadius
DOUBLE PRECISION, PARAMETER :: PI = 3.1415926D0, Radius = 1.5D-4
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:):: distance, origindistance, dL, csx, csy, csz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: Positions, VCseeds, Ori
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: K_global, P

!Dimensionality of the problem, 3D
ndim = 3
!The dimension of the physical domain using the format:
!3D: [Xmin, Ymin, Zmin, Xmax, Ymax, Zmax]
Box = RESHAPE((/ 0.D0, 0.D0, 0.D0, 0.01D0, 0.01D0, 0.01D0 /), (/ 3, 2 /))
!model configuration parameters
offset = 1.3*Radius
threshold = 1.3*Radius

!Number of neighbors
nneighbors = 16 !14
nneighbors1 = 10 !8
nneighbors2 = 10 !6
!Number of AFEM neighbors
nneighbors_afem = 46 !40
nneighbors_afem1 = 44 !34
nneighbors_afem2 = 40 !24
!Read in the Voronoi Cell seeds data
OPEN(UNIT=13, FILE="../Vonoroi/0.1_5seeds.dat", STATUS='old', ACTION='read')
READ(13,*) ngrain
ALLOCATE(VCseeds(ngrain, ndim))
DO i = 1, ngrain
    READ(13,*) VCseeds(i,1), VCseeds(i,2), VCseeds(i,3)
ENDDO
CLOSE(UNIT=13)
!Rescale the seeds data, original data is based on 1
VCseeds = VCseeds*0.01D0
!the neighboring seeds search radius
SRadius = 2.*MAX(Box(1,2)-Box(1,1), Box(2,2)-Box(2,1), Box(3,2)-Box(3,1))/INT(ngrain**(1./3.))
ALLOCATE(Sneighbor(ngrain, 100))
Sneighbor = 0
!The lattice rotation angles, thess rotation angles are counter clockwise, thus should have minus sign in rotation matrix.
! Rotation sequence: around x, y, z, respectively
ALLOCATE(Ori(ngrain + 1, ndim))
OPEN(UNIT=13, FILE="../Vonoroi/0.1_5ori.dat", STATUS='old', ACTION='read')
DO i = 1, ngrain + 1
    READ(13,*) Ori(i,1), Ori(i,2), Ori(i,3)
ENDDO
CLOSE(UNIT=13)

!the unrotated stiffness matrix, has following structure: 1=11, 2=12, 3=44
C = [230.D9, 135.D9, 117.D9] !Fe, bcc

!the mapping matrix from C to Kn and Tv
Mapping = RESHAPE((/            0.,                             0.,                   4.*SQRT(3.), &
                                          4.*SQRT(3.)/3.,        -4.*SQRT(3.)/3.,                0.,&
                                                      0.,                    2.*SQRT(3.)/7.,     -2.*SQRT(3.)/7./),(/3,3/));
!the model parameters
KnTv = MATMUL(TRANSPOSE(Mapping), C)/4. !using the whole bond length change, thus needs a dividion by 4

!record the CPU time
start = OMP_GET_WTIME()
!----------------------------------------------------------------------------------------------------------
!initialize the particle position and allocate memory for the global matrices
!----------------------------------------------------------------------------------------------------------
CALL  Initialization
!----------------------------------------------------------------------------------------------------------
!Save the initial positions to file
!----------------------------------------------------------------------------------------------------------

OPEN(UNIT = 13, FILE = 'BCC_5grains.dump', ACTION = 'WRITE', STATUS = 'REPLACE', POSITION = 'APPEND')
WRITE(13, '(A)') 'ITEM: TIMESTEP'
WRITE(13, '(I0)') 0
WRITE(13, '(A)') 'ITEM: NUMBER OF ATOMS'
WRITE(13, '(I0)') nparticle
WRITE(13, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
WRITE(13, '(2F8.1)') 1000*Box(1, 1), 1000*Box(1, 2) ! unit is mm
WRITE(13, '(2F8.1)') 1000*Box(2, 1), 1000*Box(2, 2)
WRITE(13, '(2F8.1)') 1000*Box(3, 1), 1000*Box(3, 2)
WRITE(13, '(A)') 'ITEM: ATOMS id type x y z '
DO t = 1, nparticle
    WRITE(13, '(I0, I6, 3F12.8)') t-1, grainsign(t)-1, 1000*Positions(t, 1), 1000*Positions(t, 2), 1000*Positions(t, 3)
ENDDO
CLOSE(UNIT = 13)

CALL MKL_FREE_BUFFERS

!record the CPU time
finish = OMP_GET_WTIME()
WRITE(*,*) 'Total computation time: ', finish-start, 'seconds.'

!SUBROUTINES
CONTAINS

!==================================================SUBROUTINE Initialization============================================

SUBROUTINE Initialization
IMPLICIT NONE

INTEGER:: i, j, k, m, mm, n, nn
DOUBLE PRECISION:: dis

!The N Near Neighbor search for each VC seeds
DO i = 1, ngrain
    k = 0
    DO j = 1, ngrain
        dis = SQRT((VCseeds(i,1) - VCseeds(j,1))**2 + (VCseeds(i,2) - VCseeds(j,2))**2 + (VCseeds(i,3) - VCseeds(j,3))**2)
        IF(dis .LE. SRadius .AND. j .NE. i)THEN
            k = k + 1
            Sneighbor(i,k) = j
        ENDIF
    ENDDO
ENDDO

nparticle = 0
DO i = 1, ngrain
    CALL Mbuilder(i)
ENDDO

WRITE(*,*) 'Number of particles:', nparticle

!allocate memory for the global matrices
!----------------------------------------------------------------------------
ALLOCATE(distance(nparticle,nneighbors1,2))
distance = 0.D0
ALLOCATE(origindistance(nparticle,nneighbors1,2))
origindistance = 0.D0
ALLOCATE(csx(nparticle,nneighbors1,2))
csx = 0.D0
ALLOCATE(csy(nparticle,nneighbors1,2))
csy = 0.D0
ALLOCATE(csz(nparticle,nneighbors1,2))
csz = 0.D0
ALLOCATE(dL(nparticle,nneighbors1,2))
dL = 0.D0
ALLOCATE(neighbors(nparticle,nneighbors))
neighbors = 0
ALLOCATE(neighbors1(nparticle,nneighbors1))
neighbors1 = 0
ALLOCATE(neighbors2(nparticle,nneighbors2))
neighbors2 = 0
ALLOCATE(neighbors_afem(nparticle, nneighbors_afem))
neighbors_afem = 0
ALLOCATE(neighbors_afem1(nparticle, nneighbors_afem1))
neighbors_afem1 = 0
ALLOCATE(neighbors_afem2(nparticle, nneighbors_afem2))
neighbors_afem2 = 0
ALLOCATE(K_pointer(nparticle+1, 2))
K_pointer = 0
ALLOCATE(Conn(nparticle, nneighbors_afem +1))
Conn = 0
ALLOCATE(nsign(nparticle, nneighbors))
nsign = 0
ALLOCATE(P(ndim*nparticle))
P = 0.D0

END SUBROUTINE Initialization

!================================================END Initialization======================================================

!================================================SUBROUTINE Mbuilder=================================================

SUBROUTINE Mbuilder(grain_index)
IMPLICIT NONE

DOUBLE PRECISION:: delta, x, y, z, a, Box_t(3,2), RMatrix(3,3), dis1, dis2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: Positions_t, Positions_temp
INTEGER:: i, j, k, m, n, nb, np, grain_index, particles_first_row, rows, layers, nparticle_t, nparticle_last
INTEGER, ALLOCATABLE, DIMENSION(:):: in_out, grainsign_temp

!record total numbers of particles for previous grains
IF(grain_index .EQ. 1)THEN
    nparticle_last = 0
ELSE
    nparticle_last = nparticle
ENDIF

!defined the size of the grains to be rotated, centered at the vonoroi cell
Box_t = RESHAPE((/ -SRadius, -SRadius, -SRadius, SRadius, SRadius, SRadius /), (/ 3, 2 /))

delta = 4.D0/SQRT(3.)*Radius
particles_first_row = 1 + FLOOR(( Box_t(1,2) - Box_t(1,1) ) / delta)
rows = 1 + FLOOR ( ( Box_t(2,2) - Box_t(2,1) ) / (delta/2.D0))
layers = 1 + FLOOR ( ( Box_t(3,2) - Box_t(3,1) ) / (delta/2.D0))
nparticle_t = FLOOR((layers+1)/2.0)*(particles_first_row * rows) + FLOOR((layers)/2.0)*((particles_first_row -1 ) * (rows-1))

!calculate the positions of each particle
ALLOCATE(Positions_t( nparticle_t,ndim), in_out(nparticle_t))
Positions_t = 0.
in_out = 0

n = 0
DO k = 1, layers
    z = Box_t(3,1) + delta/2. * ( k - 1 )
    IF(MOD ( k, 2 ) .EQ. 1)THEN
        DO j = 1, rows
            y = Box_t(2,1) + delta * ( j - 1 )
                DO i = 1, particles_first_row
                    x = Box_t(1,1) + delta * ( i - 1 )
                    n = n + 1
                    IF  (n .LE. nparticle_t) THEN
                        Positions_t(n,1) = x
                        Positions_t(n,2) = y
                        Positions_t(n,3) = z
                    ENDIF
                ENDDO
        ENDDO
    ELSE
        DO j = 1, rows - 1
            y = Box_t(2,1) + delta * ( j - 1 ) + delta/2.
                DO i = 1, particles_first_row - 1
                    x = Box_t(1,1) + delta * ( i - 1 ) + delta/2.
                    n = n + 1
                    IF  (n .LE. nparticle_t) THEN
                        Positions_t(n,1) = x
                        Positions_t(n,2) = y
                        Positions_t(n,3) = z
                    ENDIF
                ENDDO
        ENDDO
    ENDIF
ENDDO

!rotate the particle system by angle Ori(i,1), Ori(i,2), Ori(i,3)
!The transformatiom matrix
Rmatrix = RESHAPE((/COS(Ori(grain_index, 2))*COS(Ori(grain_index, 3)),        COS(Ori(grain_index, 1))*SIN(Ori(grain_index, 3)) + SIN(Ori(grain_index, 1))*SIN(Ori(grain_index, 2))*COS(Ori(grain_index, 3)),          SIN(Ori(grain_index, 1))*SIN(Ori(grain_index, 3)) - COS(Ori(grain_index, 1))*SIN(Ori(grain_index, 2))*COS(Ori(grain_index, 3)),&
                                         -COS(Ori(grain_index, 2))*SIN(Ori(grain_index, 3)),         COS(Ori(grain_index, 1))*COS(Ori(grain_index, 3)) - SIN(Ori(grain_index, 1))*SIN(Ori(grain_index, 2))*SIN(Ori(grain_index, 3)),           SIN(Ori(grain_index, 1))*COS(Ori(grain_index, 3)) + COS(Ori(grain_index, 1))*SIN(Ori(grain_index, 2))*SIN(Ori(grain_index, 3)),&
                                                        SIN(Ori(grain_index, 2)),                                                                                                       -SIN(Ori(grain_index, 1))*COS(Ori(grain_index, 2)),                                                                                                                                            COS(Ori(grain_index, 1))*COS(Ori(grain_index, 2))/),(/3,3/));
np = 0
DO i = 1, nparticle_t
    n = 0
    x = Positions_t(i,1)*Rmatrix(1,1) + Positions_t(i,2)*Rmatrix(1,2) + Positions_t(i,3)*Rmatrix(1,3)
    y = Positions_t(i,1)*Rmatrix(2,1) + Positions_t(i,2)*Rmatrix(2,2) + Positions_t(i,3)*Rmatrix(2,3)
    z = Positions_t(i,1)*Rmatrix(3,1) + Positions_t(i,2)*Rmatrix(3,2) + Positions_t(i,3)*Rmatrix(3,3)
    !check whether this particle falls within the vonoroi cell
    dis1 = SQRT(x**2 + y**2 + z**2)
    nb = COUNT(Sneighbor(grain_index,:) .NE. 0)
    DO j = 1, nb
        dis2 = SQRT((x + VCseeds(grain_index,1) - VCseeds(Sneighbor(grain_index,j),1))**2 + (y + VCseeds(grain_index,2) - VCseeds(Sneighbor(grain_index,j),2))**2 + (z + VCseeds(grain_index,3) - VCseeds(Sneighbor(grain_index,j),3))**2)
        IF(dis2 .GT. dis1 - offset)THEN
            n = n + 1
        ENDIF
    ENDDO
    IF (n .EQ. nb) THEN
        np = np + 1
        in_out(i) = 1
    ENDIF
ENDDO

!save the particle positions fall in the vonoroi cell temporarily
ALLOCATE(Positions_temp(np,ndim))
k = 0
DO i = 1, nparticle_t
    IF (in_out(i) .EQ. 1) THEN
        
        x = Positions_t(i,1)*Rmatrix(1,1) + Positions_t(i,2)*Rmatrix(1,2) + Positions_t(i,3)*Rmatrix(1,3)
        y = Positions_t(i,1)*Rmatrix(2,1) + Positions_t(i,2)*Rmatrix(2,2) + Positions_t(i,3)*Rmatrix(2,3)
        z = Positions_t(i,1)*Rmatrix(3,1) + Positions_t(i,2)*Rmatrix(3,2) + Positions_t(i,3)*Rmatrix(3,3)
        
        k = k + 1
        Positions_temp(k,1) = x + VCseeds(grain_index, 1)
        Positions_temp(k,2) = y + VCseeds(grain_index, 2)
        Positions_temp(k,3) = z + VCseeds(grain_index, 3)
    ENDIF
ENDDO

!save particles within the system boundary
k = 0
DO i = 1, np
    IF(Positions_temp(i,1) .LE. Box(1,2) .AND. Positions_temp(i,1) .GE. Box(1,1) .AND. Positions_temp(i,2) .LE. Box(2,2) .AND. Positions_temp(i,2) .GE. Box(2,1) .AND. Positions_temp(i,3) .LE. Box(3,2) .AND. Positions_temp(i,3) .GE. Box(3,1))THEN
        k = k + 1
    ENDIF
ENDDO

DEALLOCATE(Positions_t, in_out)
ALLOCATE(Positions_t(k,ndim), in_out(k))
in_out = 1

k = 0
DO i = 1, np
    IF(Positions_temp(i,1) .LE. Box(1,2) .AND. Positions_temp(i,1) .GE. Box(1,1) .AND. Positions_temp(i,2) .LE. Box(2,2) .AND. Positions_temp(i,2) .GE. Box(2,1) .AND. Positions_temp(i,3) .LE. Box(3,2) .AND. Positions_temp(i,3) .GE. Box(3,1))THEN
        k = k + 1
        Positions_t(k,:) = Positions_temp(i,:)
    ENDIF
ENDDO

np = k
DEALLOCATE(Positions_temp)

!check the overlapping of the temporary particles with other particle
!$OMP PARALLEL DO
DO i = 1, np
    DO j = 1, nparticle_last
        IF(SQRT((Positions(j,1) - Positions_t(i,1))**2 + (Positions(j,2) - Positions_t(i,2))**2 + (Positions(j,3) - Positions_t(i,3))**2) .LT. threshold)THEN
            !Positions(j,:) = (Positions(j,:) + Positions_t(i,:))/2.D0
            grainsign(j) = ngrain + 1
            in_out(i) = 0
        ENDIF
    ENDDO
ENDDO
!$OMP END PARALLEL DO

nparticle = nparticle + COUNT(in_out .EQ. 1)
IF(grain_index .EQ. 1)THEN
    ALLOCATE(Positions(nparticle, ndim), grainsign(nparticle))
ELSE
    ALLOCATE(Positions_temp(nparticle_last,ndim), grainsign_temp(nparticle_last))
    Positions_temp = Positions
    grainsign_temp = grainsign
    DEALLOCATE(Positions, grainsign)
    ALLOCATE(Positions(nparticle, ndim), grainsign(nparticle))
    Positions(1:nparticle_last,:) = Positions_temp
    grainsign(1:nparticle_last) = grainsign_temp
ENDIF

nparticle = nparticle - COUNT(in_out .EQ. 1)
DO i = 1, np
    !save the particle within the vonoroi cell
    IF (in_out(i) .EQ. 1) THEN
        nparticle = nparticle + 1
        Positions(nparticle,:) = Positions_t(i,:)
        grainsign(nparticle) = grain_index
    ENDIF
ENDDO

RETURN
END SUBROUTINE Mbuilder

!================================================END Mbuilder=========================================================


END PROGRAM Main