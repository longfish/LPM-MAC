!INCLUDE 'mkl_pardiso.f90'
PROGRAM MAIN
USE mkl_pardiso
USE OMP_lib
USE IFPORT
IMPLICIT NONE

INTEGER:: ndim, i, t, nparticle, nneighbors, nneighbors1, nneighbors2, nneighbors_afem, nneighbors_afem1, nneighbors_afem2, ngrain
INTEGER, ALLOCATABLE, DIMENSION(:,:):: K_pointer, nsign, neighbors, neighbors1, neighbors2, Conn, neighbors_afem, neighbors_afem1, neighbors_afem2, Sneighbor
INTEGER, ALLOCATABLE, DIMENSION(:):: IK, JK, Front, Back, Left, Right, Top, Bottom, grainsign
DOUBLE PRECISION:: Box(3,2), C(3), Mapping(3,3), KnTv(3), U_stepx, U_stepy, U_stepz, Kn1, Kn2, Tv12, start, finish, offset, threshold, SRadius
DOUBLE PRECISION, PARAMETER :: PI = 3.1415926, Radius = 1.0D-4
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:):: distance, origindistance, dL, csx, csy, csz, Kn, Tv
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: Positions, TdL_total, Tcs_sumx, Tcs_sumy, Tcs_sumz, VCseeds, Ori
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: K_global, P

!Dimensionality of the problem, 3D
ndim = 3
!The dimension of the physical domain using the format:
!3D: [Xmin, Ymin, Zmin, Xmax, Ymax, Zmax]
Box = RESHAPE((/ 0.D0, 0.D0, 0.D0, 0.01D0, 0.01D0, 0.01D0 /), (/ 3, 2 /))
!model configuration parameters

offset = 0.5D0*Radius
threshold = 1.95D0*Radius
!Number of neighbor
nneighbors = 18
nneighbors1 = 12
nneighbors2 = 6
!Number of AFEM neighbor
nneighbors_afem = 60
nneighbors_afem1 = 54
nneighbors_afem2 = 24
!Read in the Voronoi Cell seeds data
OPEN(UNIT=13, FILE="D:\Results\Polycrystals\EBSD\3Dseeds.txt", STATUS='old', ACTION='read')
READ(13,*) ngrain
ALLOCATE(VCseeds(ngrain, ndim))
DO i = 1, ngrain
    READ(13,*) VCseeds(i,1), VCseeds(i,2), VCseeds(i,3)
ENDDO
CLOSE(UNIT=13)
!Rescale the seeds data, original data is based on 1
VCseeds = VCseeds/1.D2
!the neighboring seeds search radius
SRadius = 2.*MAX(Box(1,2)-Box(1,1), Box(2,2)-Box(2,1), Box(3,2)-Box(3,1))/INT(ngrain**(1./3.))
ALLOCATE(Sneighbor(ngrain, 100))
Sneighbor = 0
!The lattice rotation angles, thess rotation angles are counter clockwise, thus should have minus sign in rotation matrix.
! Rotation sequence: around x, y, z, respectively
ALLOCATE(Ori(ngrain + 1, ndim))
OPEN(UNIT=13, FILE="D:\Results\Polycrystals\EBSD\ori_3D.txt", STATUS='old', ACTION='read')
DO i = 1, ngrain + 1
    READ(13,*) Ori(i,1), Ori(i,2), Ori(i,3)
ENDDO
CLOSE(UNIT=13)
!the unrotated stiffness matrix, has following structure: 1=11, 2=12, 3=44
C = [108.D9, 61.D9, 29.D9] !Al, fcc
!C = [169.D9, 122.D9, 75.3D9] !Cu, fcc

!the mapping matrix from C to Kn and Tv
Mapping = RESHAPE((/            0.,                             0.,                   4.*SQRT(2.), &
                                                SQRT(2.),               -SQRT(2.),              -SQRT(2.),&
                                                      0.,                       SQRT(2.)/6.,          -SQRT(2.)/6./),(/3,3/));
!the model parameters
KnTv = MATMUL(TRANSPOSE(Mapping), C)/4. !using the whole bond length change, thus needs a dividion by 4

!record the CPU time
start = OMP_GET_WTIME()
!----------------------------------------------------------------------------------------------------
!initialize the particle position and allocate memory for the global matrices
!----------------------------------------------------------------------------------------------------
CALL  Initialization
!----------------------------------------------------------------------------------------------------
Kn1 = Radius*KnTv(1)
Kn2 = Radius*KnTv(2)
Tv12 = Radius*KnTv(3)
!----------------------------------------------------------------------------------------------------------
!Save the initial positions to file
!--------------------------------------------------------------------------------------------------------------------------------------
OPEN(UNIT=13, FILE="D:\initialP.dat", ACTION="write", STATUS="replace", POSITION='append')
DO i = 1, nparticle
    WRITE(13,"(I0, I8, 3F8.4)") i, grainsign(i), Positions(i,1)*100, Positions(i,2)*100, Positions(i,3)*100
ENDDO
CLOSE(UNIT=13)
!---------------------------------------------------------------------------------------------------------------------------------------
!Search the neighbors for each particle at the initial configuration
!----------------------------------------------------------------------------------------
CALL  Searchneighbor
!----------------------------------------------------------------------------------------
!applying the boundary condition using different steps
!----------------------------------------------------------------------------------------
U_stepx = 1.D-5
U_stepy = 1.D-5
U_stepz = 1.D-5

CALL Get_K_P

CALL PARDISOsolver

CALL Output

CALL MKL_FREE_BUFFERS

!record the CPU time
finish = OMP_GET_WTIME()
WRITE(*,*) 'Total computation time: ', finish - start, 'seconds.'

!SUBROUTINES
CONTAINS

!==================================================SUBROUTINE Initialization============================================

SUBROUTINE  Initialization
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

!specify the boundary particles with applied boundary conditions
!---------------------------------------------------------------------------------------------------------
j = 0
k = 0
m = 0
mm = 0
n = 0
nn = 0
DO i = 1, nparticle
    IF(Positions(i,3) .LE. Box(3,1))THEN
        j = j + 1
    ENDIF
    IF(Positions(i,3) .GE. Box(3,2))THEN
        k = k + 1
    ENDIF
    IF(Positions(i,1) .LE. Box(1,1))THEN
        m = m + 1
    ENDIF
    IF(Positions(i,1) .GE. Box(1,2))THEN
        mm = mm + 1
    ENDIF
    IF(Positions(i,2) .LE. Box(2,1))THEN
        n = n + 1
    ENDIF
    IF(Positions(i,2) .GE. Box(2,2))THEN
        nn = nn + 1
    ENDIF
ENDDO

ALLOCATE(Bottom(j), Top(k), Left(m), Right(mm), Front(n), Back(nn))

j = 0
k = 0
m = 0
mm = 0
n = 0
nn = 0
DO i = 1, nparticle
    IF(Positions(i,3) .LE. Box(3,1))THEN
        j = j + 1
        Bottom(j) = i
    ENDIF
    IF(Positions(i,3) .GE. Box(3,2))THEN
        k = k + 1
        Top(k) = i
    ENDIF
    IF(Positions(i,1) .LE. Box(1,1))THEN
        m = m + 1
        Left(m) = i
    ENDIF
    IF(Positions(i,1) .GE. Box(1,2))THEN
        mm = mm + 1
        Right(mm) = i
    ENDIF
    IF(Positions(i,2) .LE. Box(2,1))THEN
        n = n + 1
        Front(n) = i
    ENDIF
    IF(Positions(i,2) .GE. Box(2,2))THEN
        nn = nn + 1
        Back(nn) = i
    ENDIF
ENDDO

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

delta = 2.*SQRT(2.)*Radius
particles_first_row = 1 + FLOOR(( Box_t(1,2) - Box_t(1,1) ) / delta)
rows = 1 + FLOOR ( ( Box_t(2,2) - Box_t(2,1) ) / (delta/2.0))
layers = 1 + FLOOR ( ( Box_t(3,2) - Box_t(3,1) ) / (delta/2.0))

nparticle_t = FLOOR((layers+1)/2.0)*(particles_first_row * FLOOR ( ( rows + 1) / 2.0) + ( particles_first_row - 1 ) * FLOOR ( rows / 2.0 )) &
                + FLOOR(layers/2.0)*((particles_first_row -1 ) * FLOOR ( ( rows + 1) / 2.0) + particles_first_row * FLOOR ( rows / 2.0 ))

!calculate the positions of each particle
ALLOCATE(Positions_t( nparticle_t,ndim), in_out(nparticle_t))
Positions_t = 0.
in_out = 0

n = 0
DO k = 1, layers
    z = Box_t(3,1) + delta/2. * ( k - 1 )
    IF(MOD ( k, 2 ) .EQ. 1)THEN
        DO j = 1, rows
            y = Box_t(2,1) + delta/2. * ( j - 1 )
            IF(MOD(j,2) .EQ. 1)THEN
                DO i = 1, particles_first_row
                    x = Box_t(1,1) + delta * ( i - 1 )
                    n = n + 1
                    IF  (n .LE. nparticle_t) THEN
                        Positions_t(n,1) = x
                        Positions_t(n,2) = y
                        Positions_t(n,3) = z
                    ENDIF
                ENDDO
            ELSE
                DO i = 1, particles_first_row - 1
                    x = Box_t(1,1) + delta * ( i - 1 ) + delta/2.
                    n = n + 1
                    IF  (n .LE. nparticle_t) THEN
                        Positions_t(n,1) = x
                        Positions_t(n,2) = y
                        Positions_t(n,3) = z
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
    ELSE
        DO j = 1, rows
            y = Box_t(2,1) + delta/2. * ( j - 1 )
            IF(MOD(j,2) .EQ. 1)THEN
                DO i = 1, particles_first_row - 1
                    x = Box_t(1,1) + delta * ( i - 1 ) + delta/2.
                    n = n + 1
                    IF  (n .LE. nparticle_t) THEN
                        Positions_t(n,1) = x
                        Positions_t(n,2) = y
                        Positions_t(n,3) = z
                    ENDIF
                ENDDO
            ELSE
                DO i = 1, particles_first_row
                    x = Box_t(1,1) + delta * ( i - 1 )
                    n = n + 1
                    IF  (n .LE. nparticle_t) THEN
                        Positions_t(n,1) = x
                        Positions_t(n,2) = y
                        Positions_t(n,3) = z
                    ENDIF
                ENDDO
            ENDIF
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
    IF(Positions_temp(i,1) .LE. Box(1,2) + 2.0*Radius .AND. Positions_temp(i,1) .GE. Box(1,1) - 2.0*Radius .AND. Positions_temp(i,2) .LE. Box(2,2) + 2.0*Radius .AND. Positions_temp(i,2) .GE. Box(2,1) - 2.0*Radius .AND. Positions_temp(i,3) .LE. Box(3,2) + 2.0*Radius .AND. Positions_temp(i,3) .GE. Box(3,1) - 2.0*Radius)THEN
        k = k + 1
    ENDIF
ENDDO

DEALLOCATE(Positions_t, in_out)
ALLOCATE(Positions_t(k,ndim), in_out(k))
in_out = 1

k = 0
DO i = 1, np
    IF(Positions_temp(i,1) .LE. Box(1,2) + 2.0*Radius .AND. Positions_temp(i,1) .GE. Box(1,1) - 2.0*Radius .AND. Positions_temp(i,2) .LE. Box(2,2) + 2.0*Radius .AND. Positions_temp(i,2) .GE. Box(2,1) - 2.0*Radius .AND. Positions_temp(i,3) .LE. Box(3,2) + 2.0*Radius .AND. Positions_temp(i,3) .GE. Box(3,1) - 2.0*Radius)THEN
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

!===============================================SUBROUTINE Searchneighbor============================================

SUBROUTINE  Searchneighbor
IMPLICIT NONE

DOUBLE PRECISION:: dis
INTEGER:: i, j, k, m, nb1, nb2, index1, index2, index3, index4, temp(nneighbors_afem+1)
INTEGER, ALLOCATABLE, DIMENSION(:,:):: collection
ALLOCATE(collection(nparticle, nneighbors*nneighbors))
collection = 0

K_pointer(1, 2) = 1

!finding the first and second neighbors for each particle
!$OMP PARALLEL DO PRIVATE(index1, index2, index3, dis)
DO i = 1, nparticle
    index1 = 0
    index2 = 0
    index3 = 0
    DO j = 1, nparticle
        dis = SQRT((Positions(j,1) - Positions(i,1))**2 + (Positions(j,2) - Positions(i,2))**2 + (Positions(j,3) - Positions(i,3))**2)
        IF(dis  .LT. 2.01*Radius .AND. j .NE. i) THEN !the first nearest neighbors
            index1 = index1 + 1
            neighbors1(i,index1) = j
            origindistance(i,index1,1) = dis
            
            !assigning the spring stiffness
            IF(grainsign(i) .EQ. grainsign(j))THEN
                Kn(i,index1,1) = Kn1
                Tv(i,index1,1) = Tv12
            ELSE
                Kn(i,index1,1) =  (Kn1 + Kn2)/2.! the arithmetic average
                Tv(i,index1,1) =  (Tv12 + Tv12)/2.
            ENDIF
            
            index3 = index3 + 1
            neighbors(i,index3) = j
            nsign(i,index3) = 1
        ELSEIF(dis .GT. 2.01*Radius .AND. dis .LT. 2.01*SQRT(2.)*Radius)THEN !the second nearest neighbors
            index2 = index2 + 1
            neighbors2(i,index2) = j
            origindistance(i,index2,2) = dis
            
            !assigning the spring stiffness
            IF(grainsign(i) .EQ. grainsign(j))THEN
                Kn(i,index2,2) = Kn2
                Tv(i,index2,2) = Tv12
            ELSE
                Kn(i,index2,2) = (Kn1 + Kn2)/2. ! the arithmetic mean
                Tv(i,index2,2) = (Tv12 + Tv12)/2.
            ENDIF
            
            index3 = index3 + 1
            neighbors(i,index3) = j
            nsign(i,index3) = 2
        ENDIF
    ENDDO
ENDDO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------------------------neighbors_afem1------------------------------------------------------------------------
!collecting the first afem neighbors for each particle
!$OMP PARALLEL DO PRIVATE(nb1, nb2, index1)
DO i = 1, nparticle
    index1 = 0
    nb1 = COUNT(neighbors1(i,:) .NE. 0)
    DO j = 1, nb1
        nb2 = COUNT(neighbors1(neighbors1(i,j),:) .NE. 0)
        index1 = index1 + 1
        collection(i,index1) = neighbors1(i,j)
        DO k = 1, nb2
            IF(neighbors1(neighbors1(i,j),k) .NE. i)THEN
                index1 = index1 + 1
                collection(i,index1) = neighbors1(neighbors1(i,j),k)
            ENDIF
        ENDDO
    ENDDO
ENDDO
!$OMP END PARALLEL DO

!remove the duplicates in the collection1 matrix, obtain the neighbors_afem1 matrix
!$OMP PARALLEL DO PRIVATE(nb1, nb2)
DO i = 1, nparticle
    nb1 = 1
    neighbors_afem1(i,1) = collection(i,1)
    nb2 = COUNT(collection(i,:) .NE. 0)
outer1: DO j = 2, nb2
        DO k =1,nb1
            IF(neighbors_afem1(i,k) .EQ. collection(i,j))THEN
                CYCLE outer1
            ENDIF
        ENDDO
        nb1 = nb1 + 1
        neighbors_afem1(i,nb1) = collection(i,j)
    ENDDO outer1
ENDDO
!$OMP END PARALLEL DO

collection = 0
!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------neighbors_afem2------------------------------------------------------------------------
!collecting the second afem neighbors for each particle
!$OMP PARALLEL DO PRIVATE(nb1, nb2, index1)
DO i = 1, nparticle
    index1 = 0
    nb1 = COUNT(neighbors2(i,:) .NE. 0)
    DO j = 1, nb1
        nb2 = COUNT(neighbors2(neighbors2(i,j),:) .NE. 0)
        index1 = index1 + 1
        collection(i,index1) = neighbors2(i,j)
        DO k = 1, nb2
            IF(neighbors2(neighbors2(i,j),k) .NE. i)THEN
                index1 = index1 + 1
                collection(i,index1) = neighbors2(neighbors2(i,j),k)
            ENDIF
        ENDDO
    ENDDO
ENDDO
!$OMP END PARALLEL DO

!remove the duplicates in the collection2 matrix, obtain the neighbors_afem2 matrix
!$OMP PARALLEL DO PRIVATE(nb1, nb2)
DO i = 1, nparticle
    nb1 = 1
    neighbors_afem2(i,1) = collection(i,1)
    nb2 = COUNT(collection(i,:) .NE. 0)
outer2: DO j = 2, nb2
        DO k =1,nb1
            IF(neighbors_afem2(i,k) .EQ. collection(i,j))THEN
                CYCLE outer2
            ENDIF
        ENDDO
        nb1 = nb1 + 1
        neighbors_afem2(i,nb1) = collection(i,j)
    ENDDO outer2
ENDDO
!$OMP END PARALLEL DO

collection = 0
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------neighbors_afem-------------------------------------------------------------------------------
!collecting the afem neighbors for each particle
!$OMP PARALLEL DO PRIVATE(nb1, nb2, index1, index3, index4)
DO i = 1, nparticle
    index1 = 0
    index3 = 0
    index4 = 0
    nb1 = COUNT(neighbors(i,:) .NE. 0)
    DO j = 1, nb1
        IF(nsign(i,j) .EQ. 1)THEN
            index3 = index3 + 1
            nb2 = COUNT(neighbors1(neighbors1(i,index3),:) .NE. 0)
            index1 = index1 + 1
            collection(i,index1) = neighbors1(i,index3)
            DO k = 1, nb2
                IF(neighbors1(neighbors1(i,index3),k) .NE. i)THEN
                    index1 = index1 + 1
                    collection(i,index1) = neighbors1(neighbors1(i,index3),k)
                ENDIF
            ENDDO
        ELSEIF(nsign(i,j) .EQ. 2)THEN
            index4 = index4 + 1
            nb2 = COUNT(neighbors2(neighbors2(i,index4),:) .NE. 0)
            index1 = index1 + 1
            collection(i,index1) = neighbors2(i,index4)
            DO m = 1, nb2
                IF(neighbors2(neighbors2(i,index4),m) .NE. i)THEN
                    index1 = index1 + 1
                    collection(i,index1) = neighbors2(neighbors2(i,index4),m)
                ENDIF
            ENDDO
        ENDIF
    ENDDO
ENDDO
!$OMP END PARALLEL DO

!remove the duplicates in the collection matrix, obtain the neighbors_afem matrix
!$OMP PARALLEL DO PRIVATE(nb1, nb2)
DO i = 1, nparticle
    nb1 = 1
    neighbors_afem(i,1) = collection(i,1)
    nb2 = COUNT(collection(i,:) .NE. 0)
outer3: DO j = 2, nb2
        DO k =1,nb1
            IF(neighbors_afem(i,k) .EQ. collection(i,j))THEN
                CYCLE outer3
            ENDIF
        ENDDO
        nb1 = nb1 + 1
        neighbors_afem(i,nb1) = collection(i,j)
    ENDDO outer3
ENDDO
!$OMP END PARALLEL DO

collection = 0
!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------Conn--------------------------------------------------------------------------------------------
!collecting the afem neighbors for each particle
!$OMP PARALLEL DO PRIVATE(nb1, nb2, index2, index3, index4)
DO i = 1, nparticle
    index2 = 0
    index3 = 0
    index4 = 0
    nb1 = COUNT(neighbors(i,:) .NE. 0)
    DO j = 1, nb1
        IF(nsign(i,j) .EQ. 1)THEN
            index3 = index3 + 1
            nb2 = COUNT(neighbors1(neighbors1(i,index3),:) .NE. 0)
            index2 = index2 + 1
            collection(i,index2) = neighbors1(i,index3)
            DO k = 1, nb2
                index2 = index2 + 1
                collection(i,index2) = neighbors1(neighbors1(i,index3),k)
            ENDDO
        ELSEIF(nsign(i,j) .EQ. 2)THEN
            index4 = index4 + 1
            nb2 = COUNT(neighbors2(neighbors2(i,index4),:) .NE. 0)
            index2 = index2 + 1
            collection(i,index2) = neighbors2(i,index4)
            DO m = 1, nb2
                index2 = index2 + 1
                collection(i,index2) = neighbors2(neighbors2(i,index4),m)
            ENDDO
        ENDIF
    ENDDO
ENDDO
!$OMP END PARALLEL DO

!remove the duplicates in the collection matrix, obtain the Conn matrix
!$OMP PARALLEL DO PRIVATE(nb1, nb2)
DO i = 1, nparticle
    nb1 = 1
    Conn(i,1) = collection(i,1)
    nb2 = COUNT(collection(i,:) .NE. 0)
outer4: DO j = 2, nb2
        DO k =1,nb1
            IF(Conn(i,k) .EQ. collection(i,j))THEN
                CYCLE outer4
            ENDIF
        ENDDO
        nb1 = nb1 + 1
        Conn(i,nb1) = collection(i,j)
    ENDDO outer4
ENDDO
!$OMP END PARALLEL DO

DEALLOCATE(collection)

!sorting the Conn matrix
DO i = 1, nparticle
    temp = nparticle + 1
    nb1 = COUNT(Conn(i,:) .NE. 0)
    temp(1:nb1) = Conn(i,1:nb1)
    CALL SORTQQ (LOC(temp), nneighbors_afem + 1, SRT$INTEGER4)
    Conn(i,1:nb1) = temp(1:nb1)
ENDDO
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

DO i = 1, nparticle
    nb1 = COUNT(Conn(i,:) .NE. 0)
    index1 = 0
    DO j = 1, nb1
        IF(Conn(i,j) .GE. i)THEN
            index1 = index1 + 1
        ENDIF    
    ENDDO
    K_pointer(i,1) = index1  !numbers of conn larger than its own index for each particle
    K_pointer(i + 1, 2) = K_pointer(i, 2) + ndim*ndim*index1 - 3 !start index for each particle in the global stiffness matrix
ENDDO

!initialize the sparse stiffness matrix, using CSR storing format for the sparse matrix

ALLOCATE(JK(K_pointer(nparticle+1,2)-1))
ALLOCATE(IK(ndim*nparticle+1))
ALLOCATE(K_global(K_pointer(nparticle+1,2)-1))

END SUBROUTINE Searchneighbor

!=============================================END Searchneighbor=====================================================

!====================================================SUBROUTINE Get_K_P============================================

SUBROUTINE Get_K_P
IMPLICIT NONE

INTEGER:: i, j, k, m, n, nb, num1, num2, num3
DOUBLE PRECISION:: K_local(ndim*ndim)
K_local = 0.D0

!the stiffness matrix is a SYMMETRIC SPARSE matrix which will be stored using the CSR (Compressed Sparse Row) format
JK = 0
IK = 0
K_global = 0.D0
P = 0.D0
!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
distance = 0.D0
csx = 0.D0
csy = 0.D0
csz = 0.D0
Tcs_sumx = 0.D0
Tcs_sumy = 0.D0
Tcs_sumz = 0.D0
dL = 0.D0
TdL_total = 0.D0

!$OMP PARALLEL DO PRIVATE(nb, num1, num2)
DO i = 1, nparticle
    num1 = 0
    num2 = 0
    nb = COUNT(neighbors(i,:) .NE. 0)
    DO j = 1,nb
        IF(nsign(i,j) .EQ. 1)THEN
            num1 = num1 + 1
            distance(i,num1,1)  = SQRT((Positions(i,1) - Positions(neighbors1(i,num1),1))**2 + (Positions(i,2) - Positions(neighbors1(i,num1),2))**2 + (Positions(i,3) - Positions(neighbors1(i,num1),3))**2)
            csx(i,num1,1) = (Positions(i,1) - Positions(neighbors1(i,num1),1))/distance(i,num1,1)
            csy(i,num1,1) = (Positions(i,2) - Positions(neighbors1(i,num1),2))/distance(i,num1,1)
            csz(i,num1,1) = (Positions(i,3) - Positions(neighbors1(i,num1),3))/distance(i,num1,1)
            dL(i,num1,1) =  distance(i,num1,1) - origindistance(i,num1,1)
            Tcs_sumx(i,1) = Tcs_sumx(i,1) + Tv(i,num1,1)*csx(i,num1,1)
            Tcs_sumy(i,1) = Tcs_sumy(i,1) + Tv(i,num1,1)*csy(i,num1,1)
            Tcs_sumz(i,1) = Tcs_sumz(i,1) + Tv(i,num1,1)*csz(i,num1,1)
            TdL_total(i,1) = TdL_total(i,1) + Tv(i,num1,1)*dL(i,num1,1)
        ELSEIF(nsign(i,j) .EQ. 2)THEN
            num2 = num2 + 1
            distance(i,num2,2)  = SQRT((Positions(i,1) - Positions(neighbors2(i,num2),1))**2 + (Positions(i,2) - Positions(neighbors2(i,num2),2))**2 + (Positions(i,3) - Positions(neighbors2(i,num2),3))**2)
            csx(i,num2,2) = (Positions(i,1) - Positions(neighbors2(i,num2),1))/distance(i,num2,2)
            csy(i,num2,2) = (Positions(i,2) - Positions(neighbors2(i,num2),2))/distance(i,num2,2)
            csz(i,num2,2) = (Positions(i,3) - Positions(neighbors2(i,num2),3))/distance(i,num2,2)
            dL(i,num2,2) =  distance(i,num2,2) - origindistance(i,num2,2)
            Tcs_sumx(i,2) = Tcs_sumx(i,2) + Tv(i,num2,2)*csx(i,num2,2)
            Tcs_sumy(i,2) = Tcs_sumy(i,2) + Tv(i,num2,2)*csy(i,num2,2)
            Tcs_sumz(i,2) = Tcs_sumz(i,2) + Tv(i,num2,2)*csz(i,num2,2)
            TdL_total(i,2) = TdL_total(i,2) + Tv(i,num2,2)*dL(i,num2,2)
        ENDIF
    ENDDO
ENDDO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(nb, m, n, num1, num2, K_local)
DO i = 1, nparticle
    num1 = 0
    !Computing and assembling the stiffness matrix
    nb = COUNT(conn(i,:) .NE. 0) !number of nonzeros in conn(i,:)
    DO j = 1,nb
        IF(Conn(i,j) .EQ. i )THEN
            K_local = 0.5D0*du2dxyz(i,j)
        ELSEIF(ANY(neighbors_afem1(i,:) .EQ. Conn(i,j)) .AND. ANY(neighbors_afem2(i,:) .EQ. Conn(i,j)))THEN
            K_local = 0.5D0*du2dxyz1(i,j) + 0.5D0*du2dxyz2(i,j)
        ELSEIF(ANY(neighbors_afem1(i,:) .EQ. Conn(i,j)))THEN
            K_local = 0.5D0*du2dxyz1(i,j)
        ELSEIF(ANY(neighbors_afem2(i,:) .EQ. Conn(i,j)))THEN
            K_local = 0.5D0*du2dxyz2(i,j)
        ENDIF
        !---------------------------------------------------------------------------------------------------------row assembling------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        IF(conn(i,j) .EQ. i)THEN
            K_global(K_pointer(i,2):K_pointer(i,2)+2) = K_global(K_pointer(i,2):K_pointer(i,2)+2) + 2.*K_local(1:3)
            K_global(K_pointer(i,2)+ndim*K_pointer(i,1):K_pointer(i,2)+ndim*K_pointer(i,1)+1) = K_global(K_pointer(i,2) + ndim*K_pointer(i,1) : K_pointer(i,2) + ndim*K_pointer(i,1) + 1) + 2.*K_local(5:6)
            K_global(K_pointer(i,2)+2*ndim*K_pointer(i,1)-1) = K_global(K_pointer(i,2)+2*ndim*K_pointer(i,1) -1) + 2.*K_local(9)
        ELSEIF(conn(i,j) .GT. i)THEN
            num1 = num1 + 1
            K_global(K_pointer(i,2) + ndim*num1 : K_pointer(i,2) + ndim*num1 + 2) = K_global(K_pointer(i,2) + ndim*num1 : K_pointer(i,2) + ndim*num1 + 2) + K_local(1:3)
            K_global(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 - 1 : K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 + 1) = K_global(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 - 1 : K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 + 1) + K_local(4:6)
            K_global(K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1-3:K_pointer(i,2) + 2*ndim*K_pointer(i,1)+ndim*num1 - 1) = K_global(K_pointer(i,2) + 2*ndim*K_pointer(i,1)+ndim*num1 - 3 : K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1 - 1) + K_local(7:9)
        ELSE
            !-------------------------------------------------------------------------------------------------------column assembling---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            m = COUNT(conn(conn(i,j),:) .NE. 0)
            n = COUNT(conn(conn(i,j),:) .GT. conn(i,j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(conn(conn(i,j),k) .EQ. i)THEN
                    K_global(K_pointer(conn(i,j),2) + ndim*num2 : K_pointer(conn(i,j),2) + ndim*num2 + 2) = K_global(K_pointer(conn(i,j),2) + ndim*num2 : K_pointer(conn(i,j),2) + ndim*num2 + 2) + K_local(1:7:3)
                    K_global(K_pointer(conn(i,j),2) + ndim*K_pointer(conn(i,j),1) + ndim*num2 - 1 : K_pointer(conn(i,j),2) + ndim*K_pointer(conn(i,j),1) + ndim*num2 + 1) = K_global(K_pointer(conn(i,j),2) + ndim*K_pointer(conn(i,j),1) + ndim*num2 - 1 : K_pointer(conn(i,j),2) + ndim*K_pointer(conn(i,j),1) + ndim*num2 +1) + K_local(2:8:3)
                    K_global(K_pointer(conn(i,j),2) + 2*ndim*K_pointer(conn(i,j),1) + ndim*num2 - 3 : K_pointer(conn(i,j),2) + 2*ndim*K_pointer(conn(i,j),1) + ndim*num2 - 1) = K_global(K_pointer(conn(i,j),2) + 2*ndim*K_pointer(conn(i,j),1) + ndim*num2 - 3 : K_pointer(conn(i,j),2) + 2*ndim*K_pointer(conn(i,j),1) + ndim*num2-1) + K_local(3:9:3)
                ENDIF
            ENDDO
        ENDIF
        !the JK array
        !------------------------------------------------------------------------------------------
        IF(conn(i,j) .EQ. i)THEN
            JK(K_pointer(i,2)) =  ndim*conn(i,j) - 2
            JK(K_pointer(i,2) + 1) =  ndim*conn(i,j) - 1
            JK(K_pointer(i,2) + 2) =  ndim*conn(i,j)
            JK(K_pointer(i,2) + ndim*K_pointer(i,1)) = ndim*conn(i,j) - 1
            JK(K_pointer(i,2) + ndim*K_pointer(i,1)+1) = ndim*conn(i,j)
            JK(K_pointer(i,2) + 2*ndim*K_pointer(i,1)-1) = ndim*conn(i,j)
        ELSEIF(conn(i,j) .GT. i)THEN
            JK(K_pointer(i,2) + ndim*num1) =  ndim*conn(i,j) - 2
            JK(K_pointer(i,2) + ndim*num1 + 1) =  ndim*conn(i,j) - 1
            JK(K_pointer(i,2) + ndim*num1 + 2) =  ndim*conn(i,j)
            JK(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1-1) = ndim*conn(i,j) - 2
            JK(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1) = ndim*conn(i,j) - 1
            JK(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 + 1) = ndim*conn(i,j)
            JK(K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1 - 3) = ndim*conn(i,j) - 2
            JK(K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1 - 2) = ndim*conn(i,j) - 1
            JK(K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1 - 1) = ndim*conn(i,j)
        ENDIF
        !-----------------------------------------------------------------------------------------
    ENDDO
    !the IK array
    !------------------------------------------------------------------------
    IK(ndim*i-2) = K_pointer(i,2)
    IK(ndim*i-1) = K_pointer(i,2) + ndim*K_pointer(i,1)
    IK(ndim*i) = K_pointer(i,2) + 2*ndim*K_pointer(i,1) - 1
    !------------------------------------------------------------------------
    !the force vector
    P(ndim*i - 2 : ndim*i) = - dudxyz(i)
ENDDO
!$OMP END PARALLEL DO
IK(ndim*nparticle + 1) = K_pointer(nparticle + 1,2)


!-------------------------------------------------------------------------------------------------------------------
!----------------------------------ESSENTIAL BOUDARY CONDITION------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!enforcing the applied displacement boundary conditions: Left Negative x
num1 = SIZE(Left)
DO i = 1, num1
    nb = COUNT(Conn(Left(i),:) .NE. 0)
    num3 = 0
    DO j = 1, nb
        IF(Conn(Left(i),j) .EQ. Left(i))THEN
            !modify the force vector
            P(ndim*Left(i) - 2) =  - U_stepx
            P(ndim*Left(i) - 1) = P(ndim*Left(i) - 1) + U_stepx*K_global(K_pointer(Left(i),2) + 1)
            P(ndim*Left(i)) = P(ndim*Left(i)) + U_stepx*K_global(K_pointer(Left(i),2) + 2)
            !zero out the columns corresponding to the x component
            K_global(K_pointer(Left(i),2)) = 1.D0
            K_global(K_pointer(Left(i),2) + 1) = 0.D0
            K_global(K_pointer(Left(i),2) + 2) = 0.D0
        ELSEIF(Conn(Left(i),j) .LT. Left(i))THEN
            m = COUNT(Conn(Conn(Left(i),j),:) .NE. 0)
            n = COUNT(Conn(Conn(Left(i),j),:) .GT. Conn(Left(i),j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(Conn(Conn(Left(i),j),k) .EQ. Left(i))THEN
                    !modify the force vector
                    P(ndim*Conn(Left(i),j) - 2) = P(ndim*Conn(Left(i),j) - 2) + U_stepx*K_global(K_pointer(Conn(Left(i),j),2) + ndim*num2)
                    P(ndim*Conn(Left(i),j) - 1) = P(ndim*Conn(Left(i),j) - 1) + U_stepx*K_global(K_pointer(Conn(Left(i),j),2) + ndim*K_pointer(Conn(Left(i),j),1) + ndim*num2 - 1)
                    P(ndim*Conn(Left(i),j)) = P(ndim*Conn(Left(i),j)) + U_stepx*K_global(K_pointer(Conn(Left(i),j),2) + 2*ndim*K_pointer(Conn(Left(i),j),1) + ndim*num2 - 3)
                    !zero out the columns corresponding to the x component
                    K_global(K_pointer(Conn(Left(i),j),2) + ndim*num2) = 0.D0
                    K_global(K_pointer(Conn(Left(i),j),2) + ndim*K_pointer(Conn(Left(i),j),1) + ndim*num2 - 1) = 0.D0
                    K_global(K_pointer(Conn(Left(i),j),2) + 2*ndim*K_pointer(Conn(Left(i),j),1) + ndim*num2 - 3) = 0.D0
                ENDIF
            ENDDO
        ELSEIF(Conn(Left(i),j) .GT. Left(i))THEN
            !modify the force vector
            P(ndim*Conn(Left(i),j) - 2) = P(ndim*Conn(Left(i),j) - 2) + U_stepx*K_global(K_pointer(Left(i),2) + ndim*num3 + 3)
            P(ndim*Conn(Left(i),j) - 1) = P(ndim*Conn(Left(i),j) - 1) + U_stepx*K_global(K_pointer(Left(i),2)  + ndim*num3 + 4)
            P(ndim*Conn(Left(i),j)) = P(ndim*Conn(Left(i),j)) + U_stepx*K_global(K_pointer(Left(i),2) + ndim*num3 + 5)
            ! zero out the row components correspoding to the x component
            K_global(K_pointer(Left(i),2) + ndim*num3 + 3) = 0.D0
            K_global(K_pointer(Left(i),2) + ndim*num3 + 4) = 0.D0
            K_global(K_pointer(Left(i),2) + ndim*num3 + 5) = 0.D0
            
            num3 = num3 + 1
        ENDIF
    ENDDO
ENDDO

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!enforcing the applied displacement boundary conditions: Right Positive x
num1 = SIZE(Right)
DO i = 1, num1
    nb = COUNT(Conn(Right(i),:) .NE. 0)
    num3 = 0
    DO j = 1, nb
        IF(Conn(Right(i),j) .EQ. Right(i))THEN
            !modify the force vector
            P(ndim*Right(i) - 2) = U_stepx
            P(ndim*Right(i) - 1) = P(ndim*Right(i) - 1) - U_stepx*K_global(K_pointer(Right(i),2) + 1)
            P(ndim*Right(i)) = P(ndim*Right(i)) - U_stepx*K_global(K_pointer(Right(i),2) + 2)
            !zero out the columns corresponding to the x component
            K_global(K_pointer(Right(i),2)) = 1.D0
            K_global(K_pointer(Right(i),2) + 1) = 0.D0
            K_global(K_pointer(Right(i),2) + 2) = 0.D0
        ELSEIF(Conn(Right(i),j) .LT. Right(i))THEN
            m = COUNT(Conn(Conn(Right(i),j),:) .NE. 0)
            n = COUNT(Conn(Conn(Right(i),j),:) .GT. Conn(Right(i),j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(Conn(Conn(Right(i),j),k) .EQ. Right(i))THEN
                    !modify the force vector
                    P(ndim*Conn(Right(i),j) - 2) = P(ndim*Conn(Right(i),j) - 2) - U_stepx*K_global(K_pointer(Conn(Right(i),j),2) + ndim*num2)
                    P(ndim*Conn(Right(i),j) - 1) = P(ndim*Conn(Right(i),j) - 1) - U_stepx*K_global(K_pointer(Conn(Right(i),j),2) + ndim*K_pointer(Conn(Right(i),j),1) + ndim*num2 - 1)
                    P(ndim*Conn(Right(i),j)) = P(ndim*Conn(Right(i),j)) - U_stepx*K_global(K_pointer(Conn(Right(i),j),2) + 2*ndim*K_pointer(Conn(Right(i),j),1) + ndim*num2 - 3)
                    !zero out the columns corresponding to the x component
                    K_global(K_pointer(Conn(Right(i),j),2) + ndim*num2) = 0.D0
                    K_global(K_pointer(Conn(Right(i),j),2) + ndim*K_pointer(Conn(Right(i),j),1) + ndim*num2 - 1) = 0.D0
                    K_global(K_pointer(Conn(Right(i),j),2) + 2*ndim*K_pointer(Conn(Right(i),j),1) + ndim*num2 - 3) = 0.D0
                ENDIF
            ENDDO
        ELSEIF(Conn(Right(i),j) .GT. Right(i))THEN
            !modify the force vector
            P(ndim*Conn(Right(i),j) - 2) = P(ndim*Conn(Right(i),j) - 2) - U_stepx*K_global(K_pointer(Right(i),2) + ndim*num3 + 3)
            P(ndim*Conn(Right(i),j) - 1) = P(ndim*Conn(Right(i),j) - 1) - U_stepx*K_global(K_pointer(Right(i),2) + ndim*num3 + 4)
            P(ndim*Conn(Right(i),j)) = P(ndim*Conn(Right(i),j)) - U_stepx*K_global(K_pointer(Right(i),2) + ndim*num3 + 5)
            ! zero out the row components correspoding to the x component
            K_global(K_pointer(Right(i),2) + ndim*num3 + 3) = 0.D0
            K_global(K_pointer(Right(i),2) + ndim*num3 + 4) = 0.D0
            K_global(K_pointer(Right(i),2) + ndim*num3 + 5) = 0.D0
            
            num3 = num3 + 1
        ENDIF
    ENDDO
ENDDO

!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!enforcing the applied displacement boundary conditions: Front Negative y
num1 = SIZE(Front)
DO i = 1, num1
    nb = COUNT(Conn(Front(i),:) .NE. 0)
    num3 = 0
    DO j = 1, nb
        IF(Conn(Front(i),j) .EQ. Front(i))THEN
            !modify the force vector
            P(ndim*Front(i) - 2) = P(ndim*Front(i) - 2) + U_stepy*K_global(K_pointer(Front(i),2) + 1)
            P(ndim*Front(i) - 1) = - U_stepy
            P(ndim*Front(i)) = P(ndim*Front(i)) + U_stepy*K_global(K_pointer(Front(i),2) + ndim*K_pointer(Front(i),1) + 1)
            !zero out the columns corresponding to the y component
            K_global(K_pointer(Front(i),2) + 1) = 0.D0
            K_global(K_pointer(Front(i),2) + ndim*K_pointer(Front(i),1)) = 1.D0
            K_global(K_pointer(Front(i),2) + ndim*K_pointer(Front(i),1) + 1) = 0.D0
        ELSEIF(Conn(Front(i),j) .LT. Front(i))THEN
            m = COUNT(Conn(Conn(Front(i),j),:) .NE. 0)
            n = COUNT(Conn(Conn(Front(i),j),:) .GT. Conn(Front(i),j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(Conn(Conn(Front(i),j),k) .EQ. Front(i))THEN
                    !modify the force vector
                    P(ndim*Conn(Front(i),j) - 2) = P(ndim*Conn(Front(i),j) - 2) + U_stepy*K_global(K_pointer(Conn(Front(i),j),2) + ndim*num2 + 1)
                    P(ndim*Conn(Front(i),j) - 1) = P(ndim*Conn(Front(i),j) - 1) + U_stepy*K_global(K_pointer(Conn(Front(i),j),2) + ndim*K_pointer(Conn(Front(i),j),1) + ndim*num2)
                    P(ndim*Conn(Front(i),j)) = P(ndim*Conn(Front(i),j)) + U_stepy*K_global(K_pointer(Conn(Front(i),j),2) + 2*ndim*K_pointer(Conn(Front(i),j),1) + ndim*num2 - 2)
                    !zero out the columns corresponding to the y component
                    K_global(K_pointer(Conn(Front(i),j),2) + ndim*num2 + 1) = 0.D0
                    K_global(K_pointer(Conn(Front(i),j),2) + ndim*K_pointer(Conn(Front(i),j),1) + ndim*num2) = 0.D0
                    K_global(K_pointer(Conn(Front(i),j),2) + 2*ndim*K_pointer(Conn(Front(i),j),1) + ndim*num2 - 2) = 0.D0
                ENDIF
            ENDDO
        ELSEIF(Conn(Front(i),j) .GT. Front(i))THEN
            !modify the force vector
            P(ndim*Conn(Front(i),j) - 2) = P(ndim*Conn(Front(i),j) - 2) + U_stepy*K_global(K_pointer(Front(i),2) + ndim*K_pointer(Front(i),1) + ndim*num3 + 2)
            P(ndim*Conn(Front(i),j) - 1) = P(ndim*Conn(Front(i),j) - 1) + U_stepy*K_global(K_pointer(Front(i),2) + ndim*K_pointer(Front(i),1) + ndim*num3 + 3)
            P(ndim*Conn(Front(i),j)) = P(ndim*Conn(Front(i),j)) + U_stepy*K_global(K_pointer(Front(i),2) + ndim*K_pointer(Front(i),1) + ndim*num3 + 4)
            ! zero out the row components correspoding to the y component
            K_global(K_pointer(Front(i),2) + ndim*K_pointer(Front(i),1) + ndim*num3 + 2) = 0.D0
            K_global(K_pointer(Front(i),2) + ndim*K_pointer(Front(i),1) + ndim*num3 + 3) = 0.D0
            K_global(K_pointer(Front(i),2) + ndim*K_pointer(Front(i),1) + ndim*num3 + 4) = 0.D0
            
            num3 = num3 + 1
        ENDIF
    ENDDO
ENDDO

!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!enforcing the applied displacement boundary conditions: Back Positive y
num1 = SIZE(Back)
DO i = 1, num1
    nb = COUNT(Conn(Back(i),:) .NE. 0)
    num3 = 0
    DO j = 1, nb
        IF(Conn(Back(i),j) .EQ. Back(i))THEN
            !modify the force vector
            P(ndim*Back(i) - 2) = P(ndim*Back(i) - 2) - U_stepy*K_global(K_pointer(Back(i),2) + 1)
            P(ndim*Back(i) - 1) = U_stepy
            P(ndim*Back(i)) = P(ndim*Back(i)) - U_stepy*K_global(K_pointer(Back(i),2) + ndim*K_pointer(Back(i),1) + 1)
            !zero out the columns corresponding to the y component
            K_global(K_pointer(Back(i),2) + 1) = 0.D0
            K_global(K_pointer(Back(i),2) + ndim*K_pointer(Back(i),1)) = 1.D0
            K_global(K_pointer(Back(i),2) + ndim*K_pointer(Back(i),1) + 1) = 0.D0
        ELSEIF(Conn(Back(i),j) .LT. Back(i))THEN
            m = COUNT(Conn(conn(Back(i),j),:) .NE. 0)
            n = COUNT(Conn(conn(Back(i),j),:) .GT. Conn(Back(i),j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(Conn(Conn(Back(i),j),k) .EQ. Back(i))THEN
                    !modify the force vector
                    P(ndim*Conn(Back(i),j) - 2) = P(ndim*Conn(Back(i),j) - 2) - U_stepy*K_global(K_pointer(Conn(Back(i),j),2) + ndim*num2 + 1)
                    P(ndim*Conn(Back(i),j) - 1) = P(ndim*Conn(Back(i),j) - 1) - U_stepy*K_global(K_pointer(Conn(Back(i),j),2) + ndim*K_pointer(Conn(Back(i),j),1) + ndim*num2)
                    P(ndim*Conn(Back(i),j)) = P(ndim*Conn(Back(i),j)) - U_stepy*K_global(K_pointer(Conn(Back(i),j),2) + 2*ndim*K_pointer(Conn(Back(i),j),1) + ndim*num2 - 2)
                    !zero out the columns corresponding to the y component
                    K_global(K_pointer(Conn(Back(i),j),2) + ndim*num2 + 1) = 0.D0
                    K_global(K_pointer(Conn(Back(i),j),2) + ndim*K_pointer(Conn(Back(i),j),1) + ndim*num2) = 0.D0
                    K_global(K_pointer(Conn(Back(i),j),2) + 2*ndim*K_pointer(Conn(Back(i),j),1) + ndim*num2 - 2) = 0.D0
                ENDIF
            ENDDO
        ELSEIF(Conn(Back(i),j) .GT. Back(i))THEN
            !modify the force vector
            P(ndim*Conn(Back(i),j) - 2) = P(ndim*Conn(Back(i),j) - 2) - U_stepy*K_global(K_pointer(Back(i),2) + ndim*K_pointer(Back(i),1) + ndim*num3 + 2)
            P(ndim*Conn(Back(i),j) - 1) = P(ndim*Conn(Back(i),j) - 1) - U_stepy*K_global(K_pointer(Back(i),2) + ndim*K_pointer(Back(i),1) + ndim*num3 + 3)
            P(ndim*Conn(Back(i),j)) = P(ndim*Conn(Back(i),j)) - U_stepy*K_global(K_pointer(Back(i),2) + ndim*K_pointer(Back(i),1) + ndim*num3 + 4)
            ! zero out the row components correspoding to the y component
            K_global(K_pointer(Back(i),2) + ndim*K_pointer(Back(i),1) + ndim*num3 + 2) = 0.D0
            K_global(K_pointer(Back(i),2) + ndim*K_pointer(Back(i),1) + ndim*num3 + 3) = 0.D0
            K_global(K_pointer(Back(i),2) + ndim*K_pointer(Back(i),1) + ndim*num3 + 4) = 0.D0
            
            num3 = num3 + 1
        ENDIF
    ENDDO
ENDDO

!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!enforcing the applied displacement boundary conditions: Top Positive z
num1 = SIZE(Top)
DO i = 1, num1
    nb = COUNT(Conn(Top(i),:) .NE. 0)
    num3 = 0
    DO j = 1, nb
        IF(Conn(Top(i),j) .EQ. Top(i))THEN
            !modify the force vector
            P(ndim*Top(i) - 2) = P(ndim*Top(i) - 2) - U_stepz*K_global(K_pointer(Top(i),2) + 2)
            P(ndim*Top(i) - 1) = P(ndim*Top(i) - 1) - U_stepz*K_global(K_pointer(Top(i),2) + ndim*K_pointer(Top(i),1) + 1)
            P(ndim*Top(i)) = U_stepz
            !zero out the columns corresponding to the z component
            K_global(K_pointer(Top(i),2) + 2) = 0.D0
            K_global(K_pointer(Top(i),2) + ndim*K_pointer(Top(i),1) + 1) = 0.D0
            K_global(K_pointer(Top(i),2) + 2*ndim*K_pointer(Top(i),1) - 1) = 1.D0
        ELSEIF(Conn(Top(i),j) .LT. Top(i))THEN
            m = COUNT(Conn(Conn(Top(i),j),:) .NE. 0)
            n = COUNT(Conn(Conn(Top(i),j),:) .GT. Conn(Top(i),j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(Conn(Conn(Top(i),j),k) .EQ. Top(i))THEN
                    !modify the force vector
                    P(ndim*Conn(Top(i),j) - 2) = P(ndim*Conn(Top(i),j) - 2) - U_stepz*K_global(K_pointer(Conn(Top(i),j),2) + ndim*num2 + 2)
                    P(ndim*Conn(Top(i),j) - 1) = P(ndim*Conn(Top(i),j) - 1) - U_stepz*K_global(K_pointer(Conn(Top(i),j),2) + ndim*K_pointer(Conn(Top(i),j),1) + ndim*num2 + 1)
                    P(ndim*Conn(Top(i),j)) = P(ndim*Conn(Top(i),j)) - U_stepz*K_global(K_pointer(Conn(Top(i),j),2) + 2*ndim*K_pointer(Conn(Top(i),j),1) + ndim*num2 - 1)
                    !zero out the columns corresponding to the z component
                    K_global(K_pointer(Conn(Top(i),j),2) + ndim*num2 + 2) = 0.D0
                    K_global(K_pointer(Conn(Top(i),j),2) + ndim*K_pointer(Conn(Top(i),j),1) + ndim*num2 + 1) = 0.D0
                    K_global(K_pointer(Conn(Top(i),j),2) + 2*ndim*K_pointer(Conn(Top(i),j),1) + ndim*num2 - 1) = 0.D0
                ENDIF
            ENDDO
        ELSEIF(Conn(Top(i),j) .GT. Top(i))THEN
            !modify the force vector
            P(ndim*Conn(Top(i),j) - 2) = P(ndim*Conn(Top(i),j) - 2) - U_stepz*K_global(K_pointer(Top(i),2) + 2*ndim*K_pointer(Top(i),1) + ndim*num3)
            P(ndim*Conn(Top(i),j) - 1) = P(ndim*Conn(Top(i),j) - 1) - U_stepz*K_global(K_pointer(Top(i),2) + 2*ndim*K_pointer(Top(i),1) + ndim*num3 + 1)
            P(ndim*Conn(Top(i),j)) = P(ndim*Conn(Top(i),j)) - U_stepz*K_global(K_pointer(Top(i),2) + 2*ndim*K_pointer(Top(i),1) + ndim*num3 + 2)
            ! zero out the row components correspoding to the z component
            K_global(K_pointer(Top(i),2) + 2*ndim*K_pointer(Top(i),1) + ndim*num3) = 0.D0
            K_global(K_pointer(Top(i),2) + 2*ndim*K_pointer(Top(i),1) + ndim*num3 + 1) = 0.D0
            K_global(K_pointer(Top(i),2) + 2*ndim*K_pointer(Top(i),1) + ndim*num3 + 2) = 0.D0
            
            num3 = num3 + 1
        ENDIF
    ENDDO
ENDDO

!-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!enforcing the applied displacement boundary conditions: Bottom negative z
num1 = SIZE(Bottom)
DO i = 1, num1
    nb = COUNT(Conn(Bottom(i),:) .NE. 0)
    num3 = 0
    DO j = 1, nb
        IF(Conn(Bottom(i),j) .EQ. Bottom(i))THEN
            !modify the force vector
            P(ndim*Bottom(i) - 2) = P(ndim*Bottom(i) - 2) + U_stepz*K_global(K_pointer(Bottom(i),2) + 2)
            P(ndim*Bottom(i) - 1) = P(ndim*Bottom(i) - 1) + U_stepz*K_global(K_pointer(Bottom(i),2) + ndim*K_pointer(Bottom(i),1) + 1)
            P(ndim*Bottom(i)) = - U_stepz
            !zero out the columns corresponding to the z component
            K_global(K_pointer(Bottom(i),2) + 2) = 0.D0
            K_global(K_pointer(Bottom(i),2) + ndim*K_pointer(Bottom(i),1) + 1) = 0.D0
            K_global(K_pointer(Bottom(i),2) + 2*ndim*K_pointer(Bottom(i),1) - 1) = 1.D0
        ELSEIF(Conn(Bottom(i),j) .LT. Bottom(i))THEN
            m = COUNT(Conn(Conn(Bottom(i),j),:) .NE. 0)
            n = COUNT(Conn(Conn(Bottom(i),j),:) .GT. Conn(Bottom(i),j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(Conn(Conn(Bottom(i),j),k) .EQ. Bottom(i))THEN
                    !modify the force vector
                    P(ndim*Conn(Bottom(i),j) - 2) = P(ndim*Conn(Bottom(i),j) - 2) + U_stepz*K_global(K_pointer(Conn(Bottom(i),j),2) + ndim*num2 + 2)
                    P(ndim*Conn(Bottom(i),j) - 1) = P(ndim*Conn(Bottom(i),j) - 1) + U_stepz*K_global(K_pointer(Conn(Bottom(i),j),2) + ndim*K_pointer(Conn(Bottom(i),j),1) + ndim*num2 + 1)
                    P(ndim*Conn(Bottom(i),j)) = P(ndim*Conn(Bottom(i),j)) + U_stepz*K_global(K_pointer(Conn(Bottom(i),j),2) + 2*ndim*K_pointer(Conn(Bottom(i),j),1) + ndim*num2 - 1)
                    !zero out the columns corresponding to the z component
                    K_global(K_pointer(Conn(Bottom(i),j),2) + ndim*num2 + 2) = 0.D0
                    K_global(K_pointer(Conn(Bottom(i),j),2) + ndim*K_pointer(Conn(Bottom(i),j),1) + ndim*num2 + 1) = 0.D0
                    K_global(K_pointer(Conn(Bottom(i),j),2) + 2*ndim*K_pointer(Conn(Bottom(i),j),1) + ndim*num2 - 1) = 0.D0
                ENDIF
            ENDDO
        ELSEIF(Conn(Bottom(i),j) .GT. Bottom(i))THEN
            !modify the force vector
            P(ndim*Conn(Bottom(i),j) - 2) = P(ndim*Conn(Bottom(i),j) - 2) + U_stepz*K_global(K_pointer(Bottom(i),2) + 2*ndim*K_pointer(Bottom(i),1) + ndim*num3)
            P(ndim*Conn(Bottom(i),j) - 1) = P(ndim*Conn(Bottom(i),j) - 1) + U_stepz*K_global(K_pointer(Bottom(i),2) + 2*ndim*K_pointer(Bottom(i),1) + ndim*num3 + 1)
            P(ndim*Conn(Bottom(i),j)) = P(ndim*Conn(Bottom(i),j)) + U_stepz*K_global(K_pointer(Bottom(i),2) + 2*ndim*K_pointer(Bottom(i),1) + ndim*num3 + 2)
            ! zero out the row components correspoding to the z component
            K_global(K_pointer(Bottom(i),2) + 2*ndim*K_pointer(Bottom(i),1) + ndim*num3) = 0.D0
            K_global(K_pointer(Bottom(i),2) + 2*ndim*K_pointer(Bottom(i),1) + ndim*num3 + 1) = 0.D0
            K_global(K_pointer(Bottom(i),2) + 2*ndim*K_pointer(Bottom(i),1) + ndim*num3 + 2) = 0.D0
            
            num3 = num3 + 1
        ENDIF
    ENDDO
ENDDO

END SUBROUTINE Get_K_P

!====================================================END Get_K_P===================================================

!====================================================FUNCTION du2dxyz================================================

FUNCTION du2dxyz(i, j)
IMPLICIT NONE

INTEGER:: i, j, k, nb
DOUBLE PRECISION:: du2dxyz(9)

du2dxyz = 0.D0

nb = COUNT(neighbors1(i,:) .NE. 0)
DO k = 1, nb
    !d2udxidxi, 1 1
    du2dxyz(1) = du2dxyz(1) + (Kn(i,k,1)*csx(i,k,1) + 1./2.*Tcs_sumx(i,1) + 1./2.*Tv(i,k,1)*SUM(csx(i,:,1)))*csx(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*(1. - csx(i,k,1)**2)/distance(i,k,1) + (Kn(i,k,1) + Tv(i,k,1))*csx(i,k,1)**2 + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*(1. - csx(i,k,1)**2)/distance(i,k,1)
    !d2udxidyi, 1 2
    du2dxyz(2) = du2dxyz(2) + (Kn(i,k,1)*csy(i,k,1) + 1./2.*Tcs_sumy(i,1) + 1./2.*Tv(i,k,1)*SUM(csy(i,:,1)))*csx(i,k,1) - (Kn(i,k,1)*dL(i,k,1) +1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*csy(i,k,1)*csx(i,k,1)/distance(i,k,1) + (Kn(i,k,1) + Tv(i,k,1))*csy(i,k,1)*csx(i,k,1) - (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)* csy(i,k,1)/distance(i,k,1)
    !d2udxidzi, 1 3
    du2dxyz(3) = du2dxyz(3) + (Kn(i,k,1)*csz(i,k,1) + 1./2.*Tcs_sumz(i,1) + 1./2.*Tv(i,k,1)*SUM(csz(i,:,1)))*csx(i,k,1) - (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*csz(i,k,1)*csx(i,k,1)/distance(i,k,1) + (Kn(i,k,1) + Tv(i,k,1))*csz(i,k,1)*csx(i,k,1) - (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)* csz(i,k,1)/distance(i,k,1)
    !d2udyidyi, 2 2
    du2dxyz(5) = du2dxyz(5) + (Kn(i,k,1)*csy(i,k,1) + 1./2.*Tcs_sumy(i,1) + 1./2.*Tv(i,k,1)*SUM(csy(i,:,1)))*csy(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*(1. - csy(i,k,1)**2)/distance(i,k,1) + (Kn(i,k,1) + Tv(i,k,1))*csy(i,k,1)**2 + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*(1. - csy(i,k,1)**2)/distance(i,k,1)
    !d2udyidzi, 2 3
    du2dxyz(6) = du2dxyz(6) + (Kn(i,k,1)*csz(i,k,1) + 1./2.*Tcs_sumz(i,1) + 1./2.*Tv(i,k,1)*SUM(csz(i,:,1)))*csy(i,k,1) - (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*csz(i,k,1)*csy(i,k,1)/distance(i,k,1) + (Kn(i,k,1) + Tv(i,k,1))*csz(i,k,1)*csy(i,k,1) - (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*csy(i,k,1)* csz(i,k,1)/distance(i,k,1)
    !d2udzidzi, 3 3
    du2dxyz(9) = du2dxyz(9) + (Kn(i,k,1)*csz(i,k,1) + 1./2.*Tcs_sumz(i,1) + 1./2.*Tv(i,k,1)*SUM(csz(i,:,1)))*csz(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*(1. - csz(i,k,1)**2)/distance(i,k,1) + (Kn(i,k,1) + Tv(i,k,1))*csz(i,k,1)**2 + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*(1. - csz(i,k,1)**2)/distance(i,k,1)
ENDDO

nb = COUNT(neighbors2(i,:) .NE. 0)
DO k = 1, nb
    !d2udxidxi, 1 1
    du2dxyz(1) = du2dxyz(1) + (Kn(i,k,2)*csx(i,k,2) + 1./2.*Tcs_sumx(i,2) + 1./2.*Tv(i,k,2)*SUM(csx(i,:,2)))*csx(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*(1. - csx(i,k,2)**2)/distance(i,k,2) + (Kn(i,k,2) + Tv(i,k,2))*csx(i,k,2)**2 + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*(1. - csx(i,k,2)**2)/distance(i,k,2)
    !d2udxidyi, 1 2
    du2dxyz(2) = du2dxyz(2) + (Kn(i,k,2)*csy(i,k,2) + 1./2.*Tcs_sumy(i,2) + 1./2.*Tv(i,k,2)*SUM(csy(i,:,2)))*csx(i,k,2) - (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*csy(i,k,2)*csx(i,k,2)/distance(i,k,2) + (Kn(i,k,2) + Tv(i,k,2))*csy(i,k,2)*csx(i,k,2) - (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)* csy(i,k,2)/distance(i,k,2)
    !d2udxidzi, 1 3
    du2dxyz(3) = du2dxyz(3) + (Kn(i,k,2)*csz(i,k,2) + 1./2.*Tcs_sumz(i,2) + 1./2.*Tv(i,k,2)*SUM(csz(i,:,2)))*csx(i,k,2) - (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*csz(i,k,2)*csx(i,k,2)/distance(i,k,2) + (Kn(i,k,2) + Tv(i,k,2))*csz(i,k,2)*csx(i,k,2) - (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)* csz(i,k,2)/distance(i,k,2)
    !d2udyidyi, 2 2
    du2dxyz(5) = du2dxyz(5) + (Kn(i,k,2)*csy(i,k,2) + 1./2.*Tcs_sumy(i,2) + 1./2.*Tv(i,k,2)*SUM(csy(i,:,2)))*csy(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*(1. - csy(i,k,2)**2)/distance(i,k,2) + (Kn(i,k,2) + Tv(i,k,2))*csy(i,k,2)**2 + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*(1. - csy(i,k,2)**2)/distance(i,k,2)
    !d2udyidzi, 2 3
    du2dxyz(6) = du2dxyz(6) + (Kn(i,k,2)*csz(i,k,2) + 1./2.*Tcs_sumz(i,2) + 1./2.*Tv(i,k,2)*SUM(csz(i,:,2)))*csy(i,k,2) - (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*csz(i,k,2)*csy(i,k,2)/distance(i,k,2) + (Kn(i,k,2) + Tv(i,k,2))*csz(i,k,2)*csy(i,k,2) - (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*csy(i,k,2)* csz(i,k,2)/distance(i,k,2)
    !d2udzidzi, 3 3
    du2dxyz(9) = du2dxyz(9) + (Kn(i,k,2)*csz(i,k,2) + 1./2.*Tcs_sumz(i,2) + 1./2.*Tv(i,k,2)*SUM(csz(i,:,2)))*csz(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*(1. - csz(i,k,2)**2)/distance(i,k,2) + (Kn(i,k,2) + Tv(i,k,2))*csz(i,k,2)**2 + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*(1. - csz(i,k,2)**2)/distance(i,k,2)
ENDDO

!d2udyidxi, 2 1
du2dxyz(4) = du2dxyz(2)
!d2udzidxi, 3 1
du2dxyz(7) = du2dxyz(3)
!d2udzidyi, 3 2
du2dxyz(8) = du2dxyz(6)

RETURN
END FUNCTION du2dxyz

!====================================================END du2dxyz=====================================================

!====================================================FUNCTION du2dxyz1================================================

FUNCTION du2dxyz1(i, j)
IMPLICIT NONE

INTEGER:: i, j, k, kk, m, n, nb1, nb2
DOUBLE PRECISION:: du2dxyz1(9)

du2dxyz1 = 0.D0

nb1 = COUNT(neighbors1(i,:) .NE. 0)
! the off diagonal part
! the first and second layers of particles in AFEM, i.e. the neighbors
DO kk = 1, nb1
    IF(neighbors1(i,kk) .EQ. Conn(i,j)) THEN
        DO k = 1, nb1 
            IF(k .EQ. kk)THEN
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxyz1(1) = du2dxyz1(1) - (Kn(i,k,1) + Tv(i,k,1))*csx(i,k,1)**2 + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*(-1. + csx(i,k,1)**2)/distance(i,k,1) - (Kn(i,k,1)*csx(i,k,1) - 1./2.*Tcs_sumx(neighbors1(i,k),1) - 1./2.*Tv(i,k,1)*SUM(csx(neighbors1(i,k),:,1)))*csx(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*(-1. + csx(i,k,1)**2)/distance(i,k,1)
                !d2udxidyj, 1 2
                du2dxyz1(2) = du2dxyz1(2) - (Kn(i,k,1) + Tv(i,k,1))*csy(i,k,1)*csx(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1) - (Kn(i,k,1)*csy(i,k,1) - 1./2.*Tcs_sumy(neighbors1(i,k),1) - 1./2.*Tv(i,k,1)*SUM(csy(neighbors1(i,k),:,1)))*csx(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1)
                !d2udxidzj, 1 3
                du2dxyz1(3) = du2dxyz1(3) - (Kn(i,k,1) + Tv(i,k,1))*csz(i,k,1)*csx(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*csx(i,k,1)*csz(i,k,1)/distance(i,k,1) - (Kn(i,k,1)*csz(i,k,1) - 1./2.*Tcs_sumz(neighbors1(i,k),1) - 1./2.*Tv(i,k,1)*SUM(csz(neighbors1(i,k),:,1)))*csx(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csz(i,k,1)/distance(i,k,1)
                !d2udyidxj, 2 1
                du2dxyz1(4) = du2dxyz1(4) - (Kn(i,k,1) + Tv(i,k,1))*csy(i,k,1)*csx(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1) - (Kn(i,k,1)*csx(i,k,1) - 1./2.*Tcs_sumx(neighbors1(i,k),1) - 1./2.*Tv(i,k,1)*SUM(csx(neighbors1(i,k),:,1)))*csy(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1)
                !d2udyidyj, 2 2
                du2dxyz1(5) = du2dxyz1(5) - (Kn(i,k,1) + Tv(i,k,1))*csy(i,k,1)**2 + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*(-1. + csy(i,k,1)**2)/distance(i,k,1) - (Kn(i,k,1)*csy(i,k,1) - 1./2.*Tcs_sumy(neighbors1(i,k),1) - 1./2.*Tv(i,k,1)*SUM(csy(neighbors1(i,k),:,1)))*csy(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*(-1. + csy(i,k,1)**2)/distance(i,k,1)
                !d2udyidzj, 2 3
                du2dxyz1(6) = du2dxyz1(6) - (Kn(i,k,1) + Tv(i,k,1))*csy(i,k,1)*csz(i,k,1)+ (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*csz(i,k,1)*csy(i,k,1)/distance(i,k,1) - (Kn(i,k,1)*csz(i,k,1) - 1./2.*Tcs_sumz(neighbors1(i,k),1) - 1./2.*Tv(i,k,1)*SUM(csz(neighbors1(i,k),:,1)))*csy(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*csz(i,k,1)*csy(i,k,1)/distance(i,k,1)
                !d2udzidxj, 3 1
                du2dxyz1(7) = du2dxyz1(7) - (Kn(i,k,1) + Tv(i,k,1))*csz(i,k,1)*csx(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*csx(i,k,1)*csz(i,k,1)/distance(i,k,1) - (Kn(i,k,1)*csx(i,k,1) - 1./2.*Tcs_sumx(neighbors1(i,k),1) - 1./2.*Tv(i,k,1)*SUM(csx(neighbors1(i,k),:,1)))*csz(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csz(i,k,1)/distance(i,k,1)
                !d2udzidyj, 3 2
                du2dxyz1(8) = du2dxyz1(8) - (Kn(i,k,1) + Tv(i,k,1))*csz(i,k,1)*csy(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*csy(i,k,1)*csz(i,k,1)/distance(i,k,1) - (Kn(i,k,1)*csy(i,k,1) - 1./2.*Tcs_sumy(neighbors1(i,k),1) - 1./2.*Tv(i,k,1)*SUM(csy(neighbors1(i,k),:,1)))*csz(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*csy(i,k,1)*csz(i,k,1)/distance(i,k,1)
                !d2udzidzj, 3 3
                du2dxyz1(9) = du2dxyz1(9) - (Kn(i,k,1) + Tv(i,k,1))*csz(i,k,1)**2 + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(i,1) + 1./2.*Tv(i,k,1)*SUM(dL(i,:,1)))*(-1. + csz(i,k,1)**2)/distance(i,k,1) - (Kn(i,k,1)*csz(i,k,1) - 1./2.*Tcs_sumz(neighbors1(i,k),1) - 1./2.*Tv(i,k,1)*SUM(csz(neighbors1(i,k),:,1)))*csz(i,k,1) + (Kn(i,k,1)*dL(i,k,1) + 1./2.*TdL_total(neighbors1(i,k),1) + 1./2.*Tv(i,k,1)*SUM(dL(neighbors1(i,k),:,1)))*(-1. + csz(i,k,1)**2)/distance(i,k,1)
                !-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            ELSE
                !----------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxyz1(1) = du2dxyz1(1) - (1./2.*Tv(i,k,1) + 1./2.*Tv(i,kk,1))*csx(i,kk,1)*csx(i,k,1)
                !d2udxidyj, 1 2
                du2dxyz1(2) = du2dxyz1(2) - (1./2.*Tv(i,k,1) + 1./2.*Tv(i,kk,1))*csy(i,kk,1)*csx(i,k,1)
                !d2udxidzj, 1 3
                du2dxyz1(3) = du2dxyz1(3) - (1./2.*Tv(i,k,1) + 1./2.*Tv(i,kk,1))*csz(i,kk,1)*csx(i,k,1)
                !d2udyidxj, 2 1
                du2dxyz1(4) = du2dxyz1(4) - (1./2.*Tv(i,k,1) + 1./2.*Tv(i,kk,1))*csx(i,kk,1)*csy(i,k,1)
                !d2udyidyj, 2 2
                du2dxyz1(5) = du2dxyz1(5) - (1./2.*Tv(i,k,1) + 1./2.*Tv(i,kk,1))*csy(i,kk,1)*csy(i,k,1)
                !d2udyidzj, 2 3
                du2dxyz1(6) = du2dxyz1(6) - (1./2.*Tv(i,k,1) + 1./2.*Tv(i,kk,1))*csz(i,kk,1)*csy(i,k,1)
                !d2udzidxj, 3 1
                du2dxyz1(7) = du2dxyz1(7) - (1./2.*Tv(i,k,1) + 1./2.*Tv(i,kk,1))*csx(i,kk,1)*csz(i,k,1)
                !d2udzidyj, 3 2
                du2dxyz1(8) = du2dxyz1(8) - (1./2.*Tv(i,k,1) + 1./2.*Tv(i,kk,1))*csy(i,kk,1)*csz(i,k,1)
                !d2udzidzj, 3 3
                du2dxyz1(9) = du2dxyz1(9) - (1./2.*Tv(i,k,1) + 1./2.*Tv(i,kk,1))*csz(i,kk,1)*csz(i,k,1)
                !----------------------------------------------------------------------------------
                nb2 = COUNT(neighbors1(neighbors1(i,k),:) .NE. 0)
                DO n = 1,nb2
                    IF(neighbors1(neighbors1(i,k),n) .EQ. Conn(i,j))THEN
                        !---------------------------------------------------------------------------------------------------
                        !d2udxidxj, 1 1
                        du2dxyz1(1) = du2dxyz1(1) - (1./2.*Tv(neighbors1(i,k),n,1) + 1./2.*Tv(i,k,1))*csx(neighbors1(i,k),n,1)*csx(i,k,1)
                        !d2udxidyj, 1 2
                        du2dxyz1(2) = du2dxyz1(2) - (1./2.*Tv(neighbors1(i,k),n,1) + 1./2.*Tv(i,k,1))*csy(neighbors1(i,k),n,1)*csx(i,k,1)
                        !d2udxidzj, 1 3
                        du2dxyz1(3) = du2dxyz1(3) - (1./2.*Tv(neighbors1(i,k),n,1) + 1./2.*Tv(i,k,1))*csz(neighbors1(i,k),n,1)*csx(i,k,1)
                        !d2udyidxj, 2 1
                        du2dxyz1(4) = du2dxyz1(4) - (1./2.*Tv(neighbors1(i,k),n,1) + 1./2.*Tv(i,k,1))*csx(neighbors1(i,k),n,1)*csy(i,k,1)
                        !d2udyidyj, 2 2
                        du2dxyz1(5) = du2dxyz1(5) - (1./2.*Tv(neighbors1(i,k),n,1) + 1./2.*Tv(i,k,1))*csy(neighbors1(i,k),n,1)*csy(i,k,1)
                        !d2udyidzj, 2 3
                        du2dxyz1(6) = du2dxyz1(6) - (1./2.*Tv(neighbors1(i,k),n,1) + 1./2.*Tv(i,k,1))*csz(neighbors1(i,k),n,1)*csy(i,k,1)
                        !d2udzidxj, 3 1
                        du2dxyz1(7) = du2dxyz1(7) - (1./2.*Tv(neighbors1(i,k),n,1) + 1./2.*Tv(i,k,1))*csx(neighbors1(i,k),n,1)*csz(i,k,1)
                        !d2udzidyj, 3 2
                        du2dxyz1(8) = du2dxyz1(8) - (1./2.*Tv(neighbors1(i,k),n,1) + 1./2.*Tv(i,k,1))*csy(neighbors1(i,k),n,1)*csz(i,k,1)
                        !d2udzidzj, 3 3
                        du2dxyz1(9) = du2dxyz1(9) - (1./2.*Tv(neighbors1(i,k),n,1) + 1./2.*Tv(i,k,1))*csz(neighbors1(i,k),n,1)*csz(i,k,1)
                        !---------------------------------------------------------------------------------------------------
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        GOTO 777
    ENDIF
ENDDO
! the third and fourth nearest neighbor in AFEM
DO kk = 1, nb1
    nb2 = COUNT(neighbors1(neighbors1(i,kk),:) .NE. 0)
    DO m = 1, nb2
        IF(neighbors1(neighbors1(i,kk),m) .EQ. Conn(i,j))THEN
            !--------------------------------------------------------------------------------------------------------
            !d2udxidxj, 1 1
            du2dxyz1(1) = du2dxyz1(1) - (1./2.*Tv(neighbors1(i,kk),m,1) + 1./2.*Tv(i,kk,1))*csx(i,kk,1)*csx(neighbors1(i,kk),m,1)
            !d2udxidyj, 1 2
            du2dxyz1(2) = du2dxyz1(2) - (1./2.*Tv(neighbors1(i,kk),m,1) + 1./2.*Tv(i,kk,1))*csx(i,kk,1)*csy(neighbors1(i,kk),m,1)
            !d2udxidzj, 1 3
            du2dxyz1(3) = du2dxyz1(3) - (1./2.*Tv(neighbors1(i,kk),m,1) + 1./2.*Tv(i,kk,1))*csx(i,kk,1)*csz(neighbors1(i,kk),m,1)
            !d2udyidxj, 2 1
            du2dxyz1(4) = du2dxyz1(4) - (1./2.*Tv(neighbors1(i,kk),m,1) + 1./2.*Tv(i,kk,1))*csy(i,kk,1)*csx(neighbors1(i,kk),m,1)
            !d2udyidyj, 2 2
            du2dxyz1(5) = du2dxyz1(5) - (1./2.*Tv(neighbors1(i,kk),m,1) + 1./2.*Tv(i,kk,1))*csy(i,kk,1)*csy(neighbors1(i,kk),m,1)
            !d2udyidzj, 2 3
            du2dxyz1(6) = du2dxyz1(6) - (1./2.*Tv(neighbors1(i,kk),m,1) + 1./2.*Tv(i,kk,1))*csy(i,kk,1)*csz(neighbors1(i,kk),m,1)
            !d2udzidxj, 3 1
            du2dxyz1(7) = du2dxyz1(7) - (1./2.*Tv(neighbors1(i,kk),m,1) + 1./2.*Tv(i,kk,1))*csz(i,kk,1)*csx(neighbors1(i,kk),m,1)
            !d2udzidyj, 3 2
            du2dxyz1(8) = du2dxyz1(8) - (1./2.*Tv(neighbors1(i,kk),m,1) + 1./2.*Tv(i,kk,1))*csz(i,kk,1)*csy(neighbors1(i,kk),m,1)
            !d2udzidzj, 3 3
            du2dxyz1(9) = du2dxyz1(9) - (1./2.*Tv(neighbors1(i,kk),m,1) + 1./2.*Tv(i,kk,1))*csz(i,kk,1)*csz(neighbors1(i,kk),m,1)
            !---------------------------------------------------------------------------------------------------------
        ENDIF
    ENDDO
ENDDO

777   RETURN
END FUNCTION du2dxyz1

!====================================================END du2dxyz1=====================================================

!====================================================FUNCTION du2dxyz2================================================

FUNCTION du2dxyz2(i, j)
IMPLICIT NONE

INTEGER:: i, j, k, kk, m, n, nb1, nb2
DOUBLE PRECISION:: du2dxyz2(9)

du2dxyz2 = 0.D0

nb1 = COUNT(neighbors2(i,:) .NE. 0)
! the off diagonal part
! the first and second layers of particles in AFEM, i.e. the neighbors
DO kk = 1, nb1
    IF(neighbors2(i,kk) .EQ. Conn(i,j)) THEN
        DO k = 1, nb1 
            IF(k .EQ. kk)THEN
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxyz2(1) = du2dxyz2(1) - (Kn(i,k,2) + Tv(i,k,2))*csx(i,k,2)**2 + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*(-1. + csx(i,k,2)**2)/distance(i,k,2) - (Kn(i,k,2)*csx(i,k,2) - 1./2.*Tcs_sumx(neighbors2(i,k),2) - 1./2.*Tv(i,k,2)*SUM(csx(neighbors2(i,k),:,2)))*csx(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*(-1. + csx(i,k,2)**2)/distance(i,k,2)
                !d2udxidyj, 1 2
                du2dxyz2(2) = du2dxyz2(2) - (Kn(i,k,2) + Tv(i,k,2))*csy(i,k,2)*csx(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2) - (Kn(i,k,2)*csy(i,k,2) - 1./2.*Tcs_sumy(neighbors2(i,k),2) - 1./2.*Tv(i,k,2)*SUM(csy(neighbors2(i,k),:,2)))*csx(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2)
                !d2udxidzj, 1 3
                du2dxyz2(3) = du2dxyz2(3) - (Kn(i,k,2) + Tv(i,k,2))*csz(i,k,2)*csx(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*csx(i,k,2)*csz(i,k,2)/distance(i,k,2) - (Kn(i,k,2)*csz(i,k,2) - 1./2.*Tcs_sumz(neighbors2(i,k),2) - 1./2.*Tv(i,k,2)*SUM(csz(neighbors2(i,k),:,2)))*csx(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csz(i,k,2)/distance(i,k,2)
                !d2udyidxj, 2 1
                du2dxyz2(4) = du2dxyz2(4) - (Kn(i,k,2) + Tv(i,k,2))*csy(i,k,2)*csx(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2) - (Kn(i,k,2)*csx(i,k,2) - 1./2.*Tcs_sumx(neighbors2(i,k),2) - 1./2.*Tv(i,k,2)*SUM(csx(neighbors2(i,k),:,2)))*csy(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2)
                !d2udyidyj, 2 2
                du2dxyz2(5) = du2dxyz2(5) - (Kn(i,k,2) + Tv(i,k,2))*csy(i,k,2)**2 + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*(-1. + csy(i,k,2)**2)/distance(i,k,2) - (Kn(i,k,2)*csy(i,k,2) - 1./2.*Tcs_sumy(neighbors2(i,k),2) - 1./2.*Tv(i,k,2)*SUM(csy(neighbors2(i,k),:,2)))*csy(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*(-1. + csy(i,k,2)**2)/distance(i,k,2)
                !d2udyidzj, 2 3
                du2dxyz2(6) = du2dxyz2(6) - (Kn(i,k,2) + Tv(i,k,2))*csy(i,k,2)*csz(i,k,2)+ (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*csz(i,k,2)*csy(i,k,2)/distance(i,k,2) - (Kn(i,k,2)*csz(i,k,2) - 1./2.*Tcs_sumz(neighbors2(i,k),2) - 1./2.*Tv(i,k,2)*SUM(csz(neighbors2(i,k),:,2)))*csy(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*csz(i,k,2)*csy(i,k,2)/distance(i,k,2)
                !d2udzidxj, 3 1
                du2dxyz2(7) = du2dxyz2(7) - (Kn(i,k,2) + Tv(i,k,2))*csz(i,k,2)*csx(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*csx(i,k,2)*csz(i,k,2)/distance(i,k,2) - (Kn(i,k,2)*csx(i,k,2) - 1./2.*Tcs_sumx(neighbors2(i,k),2) - 1./2.*Tv(i,k,2)*SUM(csx(neighbors2(i,k),:,2)))*csz(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csz(i,k,2)/distance(i,k,2)
                !d2udzidyj, 3 2
                du2dxyz2(8) = du2dxyz2(8) - (Kn(i,k,2) + Tv(i,k,2))*csz(i,k,2)*csy(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*csy(i,k,2)*csz(i,k,2)/distance(i,k,2) - (Kn(i,k,2)*csy(i,k,2) - 1./2.*Tcs_sumy(neighbors2(i,k),2) - 1./2.*Tv(i,k,2)*SUM(csy(neighbors2(i,k),:,2)))*csz(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*csy(i,k,2)*csz(i,k,2)/distance(i,k,2)
                !d2udzidzj, 3 3
                du2dxyz2(9) = du2dxyz2(9) - (Kn(i,k,2) + Tv(i,k,2))*csz(i,k,2)**2 + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(i,2) + 1./2.*Tv(i,k,2)*SUM(dL(i,:,2)))*(-1. + csz(i,k,2)**2)/distance(i,k,2) - (Kn(i,k,2)*csz(i,k,2) - 1./2.*Tcs_sumz(neighbors2(i,k),2) - 1./2.*Tv(i,k,2)*SUM(csz(neighbors2(i,k),:,2)))*csz(i,k,2) + (Kn(i,k,2)*dL(i,k,2) + 1./2.*TdL_total(neighbors2(i,k),2) + 1./2.*Tv(i,k,2)*SUM(dL(neighbors2(i,k),:,2)))*(-1. + csz(i,k,2)**2)/distance(i,k,2)
                !-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            ELSE
                !----------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxyz2(1) = du2dxyz2(1) - (1./2.*Tv(i,k,2) + 1./2.*Tv(i,kk,2))*csx(i,kk,2)*csx(i,k,2)
                !d2udxidyj, 1 2
                du2dxyz2(2) = du2dxyz2(2) - (1./2.*Tv(i,k,2) + 1./2.*Tv(i,kk,2))*csy(i,kk,2)*csx(i,k,2)
                !d2udxidzj, 1 3
                du2dxyz2(3) = du2dxyz2(3) - (1./2.*Tv(i,k,2) + 1./2.*Tv(i,kk,2))*csz(i,kk,2)*csx(i,k,2)
                !d2udyidxj, 2 1
                du2dxyz2(4) = du2dxyz2(4) - (1./2.*Tv(i,k,2) + 1./2.*Tv(i,kk,2))*csx(i,kk,2)*csy(i,k,2)
                !d2udyidyj, 2 2
                du2dxyz2(5) = du2dxyz2(5) - (1./2.*Tv(i,k,2) + 1./2.*Tv(i,kk,2))*csy(i,kk,2)*csy(i,k,2)
                !d2udyidzj, 2 3
                du2dxyz2(6) = du2dxyz2(6) - (1./2.*Tv(i,k,2) + 1./2.*Tv(i,kk,2))*csz(i,kk,2)*csy(i,k,2)
                !d2udzidxj, 3 1
                du2dxyz2(7) = du2dxyz2(7) - (1./2.*Tv(i,k,2) + 1./2.*Tv(i,kk,2))*csx(i,kk,2)*csz(i,k,2)
                !d2udzidyj, 3 2
                du2dxyz2(8) = du2dxyz2(8) - (1./2.*Tv(i,k,2) + 1./2.*Tv(i,kk,2))*csy(i,kk,2)*csz(i,k,2)
                !d2udzidzj, 3 3
                du2dxyz2(9) = du2dxyz2(9) - (1./2.*Tv(i,k,2) + 1./2.*Tv(i,kk,2))*csz(i,kk,2)*csz(i,k,2)
                !----------------------------------------------------------------------------------
                nb2 = COUNT(neighbors2(neighbors2(i,k),:) .NE. 0)
                DO n = 1,nb2
                    IF(neighbors2(neighbors2(i,k),n) .EQ. Conn(i,j))THEN
                        !---------------------------------------------------------------------------------------------------
                        !d2udxidxj, 1 1
                        du2dxyz2(1) = du2dxyz2(1) - (1./2.*Tv(neighbors2(i,k),n,2) + 1./2.*Tv(i,k,2))*csx(neighbors2(i,k),n,2)*csx(i,k,2)
                        !d2udxidyj, 1 2
                        du2dxyz2(2) = du2dxyz2(2) - (1./2.*Tv(neighbors2(i,k),n,2) + 1./2.*Tv(i,k,2))*csy(neighbors2(i,k),n,2)*csx(i,k,2)
                        !d2udxidzj, 1 3
                        du2dxyz2(3) = du2dxyz2(3) - (1./2.*Tv(neighbors2(i,k),n,2) + 1./2.*Tv(i,k,2))*csz(neighbors2(i,k),n,2)*csx(i,k,2)
                        !d2udyidxj, 2 1
                        du2dxyz2(4) = du2dxyz2(4) - (1./2.*Tv(neighbors2(i,k),n,2) + 1./2.*Tv(i,k,2))*csx(neighbors2(i,k),n,2)*csy(i,k,2)
                        !d2udyidyj, 2 2
                        du2dxyz2(5) = du2dxyz2(5) - (1./2.*Tv(neighbors2(i,k),n,2) + 1./2.*Tv(i,k,2))*csy(neighbors2(i,k),n,2)*csy(i,k,2)
                        !d2udyidzj, 2 3
                        du2dxyz2(6) = du2dxyz2(6) - (1./2.*Tv(neighbors2(i,k),n,2) + 1./2.*Tv(i,k,2))*csz(neighbors2(i,k),n,2)*csy(i,k,2)
                        !d2udzidxj, 3 1
                        du2dxyz2(7) = du2dxyz2(7) - (1./2.*Tv(neighbors2(i,k),n,2) + 1./2.*Tv(i,k,2))*csx(neighbors2(i,k),n,2)*csz(i,k,2)
                        !d2udzidyj, 3 2
                        du2dxyz2(8) = du2dxyz2(8) - (1./2.*Tv(neighbors2(i,k),n,2) + 1./2.*Tv(i,k,2))*csy(neighbors2(i,k),n,2)*csz(i,k,2)
                        !d2udzidzj, 3 3
                        du2dxyz2(9) = du2dxyz2(9) - (1./2.*Tv(neighbors2(i,k),n,2) + 1./2.*Tv(i,k,2))*csz(neighbors2(i,k),n,2)*csz(i,k,2)
                        !---------------------------------------------------------------------------------------------------
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        GOTO 888
    ENDIF
ENDDO
! the third and fourth nearest neighbor in AFEM
DO kk = 1, nb1
    nb2 = COUNT(neighbors2(neighbors2(i,kk),:) .NE. 0)
    DO m = 1, nb2
        IF(neighbors2(neighbors2(i,kk),m) .EQ. Conn(i,j))THEN
            !--------------------------------------------------------------------------------------------------------
            !d2udxidxj, 1 1
            du2dxyz2(1) = du2dxyz2(1) - (1./2.*Tv(neighbors2(i,kk),m,2) + 1./2.*Tv(i,kk,2))*csx(i,kk,2)*csx(neighbors2(i,kk),m,2)
            !d2udxidyj, 1 2
            du2dxyz2(2) = du2dxyz2(2) - (1./2.*Tv(neighbors2(i,kk),m,2) + 1./2.*Tv(i,kk,2))*csx(i,kk,2)*csy(neighbors2(i,kk),m,2)
            !d2udxidzj, 1 3
            du2dxyz2(3) = du2dxyz2(3) - (1./2.*Tv(neighbors2(i,kk),m,2) + 1./2.*Tv(i,kk,2))*csx(i,kk,2)*csz(neighbors2(i,kk),m,2)
            !d2udyidxj, 2 1
            du2dxyz2(4) = du2dxyz2(4) - (1./2.*Tv(neighbors2(i,kk),m,2) + 1./2.*Tv(i,kk,2))*csy(i,kk,2)*csx(neighbors2(i,kk),m,2)
            !d2udyidyj, 2 2
            du2dxyz2(5) = du2dxyz2(5) - (1./2.*Tv(neighbors2(i,kk),m,2) + 1./2.*Tv(i,kk,2))*csy(i,kk,2)*csy(neighbors2(i,kk),m,2)
            !d2udyidzj, 2 3
            du2dxyz2(6) = du2dxyz2(6) - (1./2.*Tv(neighbors2(i,kk),m,2) + 1./2.*Tv(i,kk,2))*csy(i,kk,2)*csz(neighbors2(i,kk),m,2)
            !d2udzidxj, 3 1
            du2dxyz2(7) = du2dxyz2(7) - (1./2.*Tv(neighbors2(i,kk),m,2) + 1./2.*Tv(i,kk,2))*csz(i,kk,2)*csx(neighbors2(i,kk),m,2)
            !d2udzidyj, 3 2
            du2dxyz2(8) = du2dxyz2(8) - (1./2.*Tv(neighbors2(i,kk),m,2) + 1./2.*Tv(i,kk,2))*csz(i,kk,2)*csy(neighbors2(i,kk),m,2)
            !d2udzidzj, 3 3
            du2dxyz2(9) = du2dxyz2(9) - (1./2.*Tv(neighbors2(i,kk),m,2) + 1./2.*Tv(i,kk,2))*csz(i,kk,2)*csz(neighbors2(i,kk),m,2)
            !---------------------------------------------------------------------------------------------------------
        ENDIF
    ENDDO
ENDDO

888   RETURN
END FUNCTION du2dxyz2

!====================================================END du2dxyz2=====================================================

!====================================================FUNCTION dudxyz=================================================

FUNCTION dudxyz(i)
IMPLICIT NONE

INTEGER, INTENT(IN):: i
INTEGER:: j, nb
DOUBLE PRECISION:: dudxyz(3)

dudxyz = 0.D0

nb = COUNT(neighbors1(i,:) .NE. 0)
DO j = 1,nb
    !dudx
    dudxyz(1) = dudxyz(1) + (2.*Kn(i,j,1)*dL(i,j,1) + 1./2.*(TdL_total(i,1) + TdL_total(neighbors1(i,j),1)) + 1./2.*Tv(i,j,1)*(SUM(dL(i,:,1)) + SUM(dL(neighbors1(i,j),:,1))))*csx(i,j,1)
    !dudy
    dudxyz(2) = dudxyz(2) + (2.*Kn(i,j,1)*dL(i,j,1) + 1./2.*(TdL_total(i,1) + TdL_total(neighbors1(i,j),1)) + 1./2.*Tv(i,j,1)*(SUM(dL(i,:,1)) + SUM(dL(neighbors1(i,j),:,1))))*csy(i,j,1)
    !dudz
    dudxyz(3) = dudxyz(3) + (2.*Kn(i,j,1)*dL(i,j,1) + 1./2.*(TdL_total(i,1) + TdL_total(neighbors1(i,j),1)) + 1./2.*Tv(i,j,1)*(SUM(dL(i,:,1)) + SUM(dL(neighbors1(i,j),:,1))))*csz(i,j,1)
ENDDO

nb = COUNT(neighbors2(i,:) .NE. 0)
DO j = 1,nb
    !dudx
    dudxyz(1) = dudxyz(1) + (2.*Kn(i,j,2)*dL(i,j,2) + 1./2.*(TdL_total(i,2) + TdL_total(neighbors2(i,j),2)) + 1./2.*Tv(i,j,2)*(SUM(dL(i,:,2)) + SUM(dL(neighbors2(i,j),:,2))))*csx(i,j,2)
    !dudy
    dudxyz(2) = dudxyz(2) + (2.*Kn(i,j,2)*dL(i,j,2) + 1./2.*(TdL_total(i,2) + TdL_total(neighbors2(i,j),2)) + 1./2.*Tv(i,j,2)*(SUM(dL(i,:,2)) + SUM(dL(neighbors2(i,j),:,2))))*csy(i,j,2)
    !dudz
    dudxyz(3) = dudxyz(3) + (2.*Kn(i,j,2)*dL(i,j,2) + 1./2.*(TdL_total(i,2) + TdL_total(neighbors2(i,j),2)) + 1./2.*Tv(i,j,2)*(SUM(dL(i,:,2)) + SUM(dL(neighbors2(i,j),:,2))))*csz(i,j,2)
ENDDO

RETURN
END FUNCTION dudxyz

!====================================================END dudxyz======================================================

!==============================================SUBROUTINE PARDISOsolver====================================================

SUBROUTINE PARDISOsolver
IMPLICIT NONE

!.. Internal solver memory pointer 
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
!.. All other variables
INTEGER maxfct, mnum, mtype, phase, error, error1, msglvl
INTEGER, ALLOCATABLE :: iparm( : )
INTEGER i, idum(1)
DOUBLE PRECISION:: ddum(1), du(ndim*nparticle)
du = 0.D0
!.. Fill all arrays containing matrix data.
maxfct = 1
mnum = 1
!.. Set up PARDISO control parameter
ALLOCATE( iparm ( 64 ) )
iparm = 0

iparm(1) = 1 ! no solver default
iparm(2) = 3 ! parallel (OpenMP) nested dissection algorithm from METIS, =0 minimum degree algorithm, =2 nested dissection algorithm from METIS
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
iparm(6) = 0 ! =0 solution on the first n compoments of x
iparm(8) = 0 ! numbers of iterative refinement steps
iparm(10) = 13 ! perturbe the pivot elements with 1E-13
iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(13) = 0 ! =0, maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
iparm(28) = 0 ! = 0, the input arrays and all internal arrays are double precision, default value. = 1, single precision
iparm(60) = 0 != 0 (default value), then In Core PARDISO is used.
!iparm(60) - version of PARDISO.
!iparm(60) controls what version of PARDISO - out-of-core (OC) version or in-core (IC) version - is used. The OC PARDISO can solve very large problems by holding the matrix factors in files on the disk. Because of that the amount of main memory required by OC PARDISO is significantly reduced.
!If iparm(60) = 0 (default value), then IC PARDISO is used.
!If iparm(60) = 1 - then IC PARDISO is used if the total memory of RAM (in megabytes) needed for storing the matrix factors is less than sum of two values of the environment variables: MKL_PARDISO_OOC_MAX_CORE_SIZE (its default value is 2000 MB) and MKL_PARDISO_OOC_MAX_SWAP_SIZE (its default value is 0 MB); otherwise OOC PARDISO is used.
!In this case amount of RAM used by OOC PARDISO can not exceed the value of MKL_PARDISO_OOC_MAX_CORE_SIZE.
!If iparm(60) = 2 - then OOC PARDISO is used.
!If iparm(60) is equal to 1 or 2, and the total peak memory needed for storing the local arrays is more than MKL_PARDISO_OOC_MAX_CORE_SIZE, the program stops with error -9. In this case, increase MKL_PARDISO_OOC_MAX_CORE_SIZE.

error  = 0 ! initialize error flag
!error        Information
!0               no error
!-1              input inconsistent
!-2              not enough memory
!-3              reordering problem
!-4              zero pivot, numerical factorization or iterative refinement problem
!-5              unclassified (internal) error
!-6              reordering failed (matrix types 11 and 13 only)
!-7              diagonal matrix is singular
!-8              32-bit integer overflow problem
!-9              not enough memory for OOC
!-10             problems with opening OOC temporary files
!-11             read/write problems with the OOC data file

!msglvl = 1 ! print statistical information
msglvl = 0 ! no print statistical information
!mtype  = -2 ! real, symmetric, indefinite
!mtype  = 1 ! structurally symmetric
mtype  = 2! real and symmetric positive definite
!mtype  = -2! real and symmetric indifinite
!mtype  = 11! real unsymmetric
ALLOCATE ( pt ( 64 ) )

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!.. Initiliaze the internal solver memory pointer. This is only necessary for the FIRST call of the PARDISO solver.
DO i = 1, 64
   pt( i )%DUMMY =  0
ENDDO

!CALL MKL_SET_NUM_THREADS(8)
!CALL OMP_SET_NUM_THREADS(8)
!.. Reordering and Symbolic Factorization, This step also allocates all memory that is necessary for the factorization
phase = 11 ! only reordering and symbolic factorization
 CALL PARDISO (pt, maxfct, mnum, mtype, phase, ndim*nparticle, K_global, IK, JK, idum, 1, iparm, msglvl, ddum, ddum, error)
WRITE(*,*) 'PARDISO: Size of factors(MB): ', MAX(iparm(15), iparm(16) + iparm(17))/1000

!phase           Solver Execution Steps
!11           Analysis
!12           Analysis, numerical factorization
!13           Analysis, numerical factorization, solve, iterative refinement
!22           Numerical factorization
!23           Numerical factorization, solve, iterative refinement
!33           Solve, iterative refinement
!331         like phase=33, but only forward substitution
!332         like phase=33, but only diagonal substitution
!333         like phase=33, but only backward substitution
!0              Release internal memory for L and U matrix number mnum
!-1             Release all internal memory for all matrices

!.. Factorization
phase = 22 ! only factorization
CALL PARDISO (pt, maxfct, mnum, mtype, phase, ndim*nparticle, K_global, IK, JK, idum, 1, iparm, msglvl, ddum, ddum, error)

!.. Back substitution and iterative refinement
phase = 33 ! only factorization
CALL PARDISO(pt, maxfct, mnum, mtype, phase, ndim*nparticle, K_global, IK, JK, idum, 1, iparm, msglvl, P, du, error)

IF (error .NE. 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
ENDIF
1000 CONTINUE
!.. Termination and release of memory
phase = -1 ! release internal memory
CALL PARDISO(pt, maxfct, mnum, mtype, phase, ndim*nparticle, ddum, idum, idum, idum, 1, iparm, msglvl, ddum, ddum, error1)

IF (error1 .NE. 0) THEN
   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
   STOP 1
ENDIF
IF ( error .NE. 0 ) STOP 1
!---------------------------------------------------------------------------------------------------------------------------------------------
!components of position
Positions(:,1) = Positions(:,1) + du(1:ndim*nparticle:3)
Positions(:,2) = Positions(:,2) + du(2:ndim*nparticle:3)
Positions(:,3) = Positions(:,3) + du(3:ndim*nparticle:3)

END SUBROUTINE PARDISOsolver

!==============================================END PARDISOsolver===========================================================

!==============================================SUBROUTINE Output=======================================================

SUBROUTINE Output
IMPLICIT NONE

INTEGER:: i, j, nb, num1, num2
DOUBLE PRECISION:: Energy

!Calculate the strain energy of the system
Energy = 0.
TdL_total = 0.

DO i = 1, nparticle
    num1 = 0
    num2 = 0
    nb = COUNT(neighbors(i,:) .NE. 0)
    DO j = 1, nb
        IF(nsign(i,j) .EQ. 1)THEN
            num1 = num1 + 1
            distance(i,num1,1)  = SQRT((Positions(i,1) - Positions(neighbors1(i,num1),1))**2 + (Positions(i,2) - Positions(neighbors1(i,num1),2))**2 + (Positions(i,3) - Positions(neighbors1(i,num1),3))**2)
            dL(i,num1,1) =  distance(i,num1,1) - origindistance(i,num1,1)
	        Energy = Energy + 2.*Kn(i,num1,1)*(dL(i,num1,1)/2.)**2
	        TdL_total(i,1) = TdL_total(i,1) + 4.*Tv(i,num1,1)*dL(i,num1,1)/2.
        ELSEIF(nsign(i,j) .EQ. 2)THEN
            num2 = num2 + 1
            distance(i,num2,2)  = SQRT((Positions(i,1) - Positions(neighbors2(i,num2),1))**2 + (Positions(i,2) - Positions(neighbors2(i,num2),2))**2 + (Positions(i,3) - Positions(neighbors2(i,num2),3))**2)
            dL(i,num2,2) =  distance(i,num2,2) - origindistance(i,num2,2)
	        Energy = Energy + 2.*Kn(i,num2,2)*(dL(i,num2,2)/2.)**2
	        TdL_total(i,2) = TdL_total(i,2) + 4.*Tv(i,num2,2)*dL(i,num2,2)/2.
        ENDIF
    ENDDO
        Energy = Energy + 1./2.*TdL_total(i,1)*SUM(dL(i,:,1)/2.) + 1./2.*TdL_total(i,2)*SUM(dL(i,:,2)/2.)
ENDDO

WRITE(*,*) "Starin energy of the system is(J): ", Energy

END SUBROUTINE Output

!==============================================END Output=======================================================


END PROGRAM main