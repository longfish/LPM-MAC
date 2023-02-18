
PROGRAM MAIN
USE OMP_lib
USE IFPORT
IMPLICIT NONE
INCLUDE 'mkl_rci.fi'


INTEGER:: ndim, particles_per_row, nparticle, rows, layers, steps, nneighbors, nneighbors1, nneighbors2, nneighbors_afem, nneighbors_afem1, nneighbors_afem2
INTEGER, ALLOCATABLE, DIMENSION(:,:):: K_pointer, neighbors, neighbors1, neighbors2, Conn, nsign, neighbors_afem,neighbors_afem1, neighbors_afem2
INTEGER, ALLOCATABLE, DIMENSION(:):: IK, JK, U_fix, U_appz
DOUBLE PRECISION:: Box(3,2), Radius, U_stepz, F_stepy, Kn1, Kn2, Tv, start, finish
DOUBLE PRECISION, PARAMETER :: E = 6.9D10, mu = 0.3
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:):: distance, origindistance, dL, csx, csy, csz
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: Positions, initialP
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: K_global, P

!Dimensionality of the problem, 3D
ndim = 3
!Number of neighbors
nneighbors = 18
nneighbors1 = 6
nneighbors2 = 12
nneighbors_afem = 60
nneighbors_afem1 = 24
nneighbors_afem2 = 54
!The dimension of the physical domain using the format:
! [Xmin, Ymin, Zmin, Xmax, Ymax, Zmax]
Box = RESHAPE((/ 0., 0., 0., 0.01, 0.01, 0.03 /), (/ 3, 2 /))
!Number of particles in the first row
particles_per_row = 21

!record the CPU time
start = OMP_GET_WTIME()
!----------------------------------------------------------------------------------------------------
!initialize the particle position and allocate memory for the global matrices
!----------------------------------------------------------------------------------------------------
CALL  Initialization
!----------------------------------------------------------------------------------------------------
Kn1 = Radius*E/(2.*(1. + mu))
Kn2 = Radius*E/(2.*(1. + mu))
Tv = Radius*E*(4.*mu - 1)/(36.*(1. + mu)*(1. - 2.*mu))
!----------------------------------------------------------------------------------------------------
!Save the initial positions to file
!----------------------------------------------
initialP = Positions
!----------------------------------------------------------------------------------------------------
!Search the neighbors for each particle at the initial configuration
!----------------------------------------------------------------------------------------------------
CALL  Searchneighbor
!----------------------------------------------------------------------------------------------------
!applying the boundary condition using different steps
!----------------------------------------------------------------------------------------------------

U_stepz = -1.D-5

CALL Get_K_P

CALL CGsolver

CALL  Output

CALL MKL_FREE_BUFFERS

!record the CPU time
finish = OMP_GET_WTIME()
WRITE(*,*) 'Total computation time: ', finish-start, 'seconds.'

!SUBROUTINES
CONTAINS

!==================================================SUBROUTINE Initialization========================================

SUBROUTINE  Initialization
IMPLICIT NONE

DOUBLE PRECISION:: delta, x, y, z
INTEGER:: i, j, k, n

delta = ( Box(1,2) - Box(1,1) ) / ( particles_per_row -  1 )
Radius = delta / 2.0
rows = 1 + FLOOR ( ( Box(2,2) - Box(2,1) ) / delta )
layers = 1 + FLOOR ( ( Box(3,2) - Box(3,1) ) / delta )
nparticle = layers*particles_per_row*rows

!calculate the positions of each particle
ALLOCATE(Positions( nparticle,ndim))
Positions = 0.
n = 0
DO k = 1, layers
    z = Box(3,1) + delta * ( k - 1 )
    DO j = 1, rows
        y = Box(2,1) + delta * ( j - 1 )
		DO i = 1, particles_per_row
			x = Box(1,1) + delta * ( i - 1 )
			n = n + 1
			Positions(n,1) = x
			Positions(n,2) = y
			Positions(n,3) = z
		ENDDO
    ENDDO
ENDDO

!specify the boundary particles
!---------------------------------------------------------------------------------------------------------
j = 0
k = 0
DO i = 1, nparticle
	IF(Positions(i,3) .GT. Positions(nparticle,3) - Radius)THEN
		j = j + 1
	ENDIF
	IF(Positions(i,3) .LT. Positions(1,3) + Radius)THEN
		k = k + 1
	ENDIF
ENDDO

ALLOCATE(U_fix(j), U_appz(k))

j = 0
k = 0
DO i = 1, nparticle
	IF(Positions(i,3) .GT. Positions(nparticle,3) - Radius)THEN
		j = j + 1
		U_fix(j) = i
	ENDIF
	IF(Positions(i,3) .LT. Positions(1,3) + Radius)THEN
		k = k + 1
		U_appz(k) = i
	ENDIF
ENDDO

!allocate memory for the global matrices
!----------------------------------------------------------------------------
ALLOCATE(initialP( nparticle,ndim))
initialP = 0.
ALLOCATE(distance(nparticle,nneighbors2,2))
distance = 0.D0
ALLOCATE(origindistance(nparticle,nneighbors2,2))
origindistance = 0.D0
ALLOCATE(csx(nparticle,nneighbors2,2))
csx = 0.D0
ALLOCATE(csy(nparticle,nneighbors2,2))
csy = 0.D0
ALLOCATE(csz(nparticle,nneighbors2,2))
csz = 0.D0
ALLOCATE(dL(nparticle,nneighbors2,2))
dL = 0.D0
ALLOCATE(neighbors(nparticle,nneighbors))
neighbors = 0
ALLOCATE(neighbors1(nparticle,nneighbors1))
neighbors1 = 0
ALLOCATE(neighbors2(nparticle,nneighbors2))
neighbors2 = 0
ALLOCATE(nsign(nparticle,nneighbors))
nsign = 0
ALLOCATE(neighbors_afem(nparticle, nneighbors_afem))
neighbors_afem = 0
ALLOCATE(neighbors_afem1(nparticle, nneighbors_afem1))
neighbors_afem1 = 0
ALLOCATE(neighbors_afem2(nparticle, nneighbors_afem2))
neighbors_afem2 = 0
ALLOCATE(K_pointer(nparticle+1, 2))
K_pointer = 0
ALLOCATE(conn(nparticle, nneighbors_afem+1))
Conn = 0
ALLOCATE(P(ndim*nparticle))
P = 0.D0

END SUBROUTINE Initialization

!================================================END Initialization=================================================

!===============================================SUBROUTINE Searchneighbor===========================================

SUBROUTINE  Searchneighbor
IMPLICIT NONE

DOUBLE PRECISION:: dis
INTEGER:: i, j, k, m, nb1, nb2, index1, index2, index3, index4, temp(nneighbors_afem+1)
INTEGER, ALLOCATABLE, DIMENSION(:,:):: collection
ALLOCATE(collection(nparticle, nneighbors*nneighbors))
collection = 0

K_pointer(1, 2) = 1

!$OMP PARALLEL DO PRIVATE(index1, index2, index3, dis)
DO i = 1, nparticle
    index1 = 0
    index2 = 0
    index3 = 0
    DO j = i - 2*particles_per_row*rows, i + 2*particles_per_row*rows
		IF(j .GE. 1 .AND. j .LE. nparticle .AND. j .NE. i)THEN
			dis = SQRT((Positions(j,1) - Positions(i,1))**2 + (Positions(j,2) - Positions(i,2))**2 + (Positions(j,3) - Positions(i,3))**2)
			IF(dis .LT. 2.01*Radius) THEN !the first nearest neighbors
				index1 = index1 + 1
				neighbors1(i,index1) = j
				origindistance(i,index1,1) = dis
            
				index3 = index3 + 1
				nsign(i,index3) = 1
				neighbors(i,index3) = j
			ELSEIF(dis .GT. 2.01*Radius .AND. dis .LT. 2.01*SQRT(2.)*Radius)THEN !the second nearest neighbors
				index2 = index2 + 1
				neighbors2(i,index2) = j
				origindistance(i,index2,2) = dis
            
				index3 = index3 + 1
				nsign(i,index3) = 2
				neighbors(i,index3) = j
			ENDIF
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

!remove the duplicates in the collection matrix, obtain the neighbors_afem1 matrix
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

!remove the duplicates in the collection matrix, obtain the neighbors_afem2 matrix
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
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------neighbors_afem--------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------Conn--------------------------------------------------------------------------
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

!=============================================END Searchneighbor====================================================

!====================================================SUBROUTINE Get_K_P=============================================

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
dL = 0.D0

!$OMP PARALLEL DO PRIVATE(nb, num1, num2)
DO i = 1, nparticle
    num1 = 0
    num2 = 0
    nb = COUNT(neighbors(i,:) .NE. 0) !number of nonzeros in neighbors(i,:)
    DO j = 1,nb
        IF(nsign(i,j) .EQ. 1)THEN
            num1 = num1 + 1
            distance(i,num1,1)  = SQRT((Positions(i,1) - Positions(neighbors1(i,num1),1))**2 + (Positions(i,2) - Positions(neighbors1(i,num1),2))**2 + (Positions(i,3) - Positions(neighbors1(i,num1),3))**2)
            csx(i,num1,1) = (Positions(i,1) - Positions(neighbors1(i,num1),1))/distance(i,num1,1)
            csy(i,num1,1) = (Positions(i,2) - Positions(neighbors1(i,num1),2))/distance(i,num1,1)
            csz(i,num1,1) = (Positions(i,3) - Positions(neighbors1(i,num1),3))/distance(i,num1,1)
            dL(i,num1,1) =  distance(i,num1,1) - origindistance(i,num1,1)
        ELSEIF(nsign(i,j) .EQ. 2)THEN
            num2 = num2 + 1
            distance(i,num2,2)  = SQRT((Positions(i,1) - Positions(neighbors2(i,num2),1))**2 + (Positions(i,2) - Positions(neighbors2(i,num2),2))**2 + (Positions(i,3) - Positions(neighbors2(i,num2),3))**2)
            csx(i,num2,2) = (Positions(i,1) - Positions(neighbors2(i,num2),1))/distance(i,num2,2)
            csy(i,num2,2) = (Positions(i,2) - Positions(neighbors2(i,num2),2))/distance(i,num2,2)
            csz(i,num2,2) = (Positions(i,3) - Positions(neighbors2(i,num2),3))/distance(i,num2,2)
            dL(i,num2,2) =  distance(i,num2,2) - origindistance(i,num2,2)
        ENDIF
    ENDDO
ENDDO
!$OMP END PARALLEL DO

!--------------------------------------------------------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(nb, m, n, num1, num2, K_local)
DO i = 1, nparticle
    num1 = 0
    !Computing and assembling the stiffness matrix
    nb = COUNT(Conn(i,:) .NE. 0) !number of nonzeros in conn(i,:)
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
        IF(Conn(i,j) .EQ. i)THEN
            K_global(K_pointer(i,2):K_pointer(i,2)+2) = K_global(K_pointer(i,2):K_pointer(i,2)+2) + 2.*K_local(1:3)
            K_global(K_pointer(i,2)+ndim*K_pointer(i,1):K_pointer(i,2)+ndim*K_pointer(i,1)+1) = K_global(K_pointer(i,2) + ndim*K_pointer(i,1) : K_pointer(i,2) + ndim*K_pointer(i,1) + 1) + 2.*K_local(5:6)
            K_global(K_pointer(i,2)+2*ndim*K_pointer(i,1)-1) = K_global(K_pointer(i,2)+2*ndim*K_pointer(i,1) -1) + 2.*K_local(9)
        ELSEIF(Conn(i,j) .GT. i)THEN
            num1 = num1 + 1
            K_global(K_pointer(i,2) + ndim*num1 : K_pointer(i,2) + ndim*num1 + 2) = K_global(K_pointer(i,2) + ndim*num1 : K_pointer(i,2) + ndim*num1 + 2) + K_local(1:3)
            K_global(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 - 1 : K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 + 1) = K_global(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 - 1 : K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 + 1) + K_local(4:6)
            K_global(K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1-3:K_pointer(i,2) + 2*ndim*K_pointer(i,1)+ndim*num1 - 1) = K_global(K_pointer(i,2) + 2*ndim*K_pointer(i,1)+ndim*num1 - 3 : K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1 - 1) + K_local(7:9)
        ELSE
            !-------------------------------------------------------------------------------------------------------column assembling---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            m = COUNT(Conn(Conn(i,j),:) .NE. 0)
            n = COUNT(Conn(Conn(i,j),:) .GT. Conn(i,j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(Conn(Conn(i,j),k) .EQ. i)THEN
                    K_global(K_pointer(Conn(i,j),2) + ndim*num2 : K_pointer(Conn(i,j),2) + ndim*num2 + 2) = K_global(K_pointer(Conn(i,j),2) + ndim*num2 : K_pointer(Conn(i,j),2) + ndim*num2 + 2) + K_local(1:7:3)
                    K_global(K_pointer(Conn(i,j),2) + ndim*K_pointer(Conn(i,j),1) + ndim*num2 - 1 : K_pointer(Conn(i,j),2) + ndim*K_pointer(Conn(i,j),1) + ndim*num2 + 1) = K_global(K_pointer(conn(i,j),2) + ndim*K_pointer(conn(i,j),1) + ndim*num2 - 1 : K_pointer(conn(i,j),2) + ndim*K_pointer(conn(i,j),1) + ndim*num2 +1) + K_local(2:8:3)
                    K_global(K_pointer(Conn(i,j),2) + 2*ndim*K_pointer(Conn(i,j),1) + ndim*num2 - 3 : K_pointer(Conn(i,j),2) + 2*ndim*K_pointer(Conn(i,j),1) + ndim*num2 - 1) = K_global(K_pointer(conn(i,j),2) + 2*ndim*K_pointer(conn(i,j),1) + ndim*num2 - 3 : K_pointer(conn(i,j),2) + 2*ndim*K_pointer(conn(i,j),1) + ndim*num2-1) + K_local(3:9:3)
                ENDIF
            ENDDO
        ENDIF
        !the JK array
        !------------------------------------------------------------------------------------------
        IF(Conn(i,j) .EQ. i)THEN
            JK(K_pointer(i,2)) =  ndim*Conn(i,j) - 2
            JK(K_pointer(i,2) + 1) =  ndim*Conn(i,j) - 1
            JK(K_pointer(i,2) + 2) =  ndim*Conn(i,j)
            JK(K_pointer(i,2) + ndim*K_pointer(i,1)) = ndim*Conn(i,j) - 1
            JK(K_pointer(i,2) + ndim*K_pointer(i,1)+1) = ndim*Conn(i,j)
            JK(K_pointer(i,2) + 2*ndim*K_pointer(i,1)-1) = ndim*Conn(i,j)
        ELSEIF(Conn(i,j) .GT. i)THEN
            JK(K_pointer(i,2) + ndim*num1) =  ndim*Conn(i,j) - 2
            JK(K_pointer(i,2) + ndim*num1 + 1) =  ndim*Conn(i,j) - 1
            JK(K_pointer(i,2) + ndim*num1 + 2) =  ndim*Conn(i,j)
            JK(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1-1) = ndim*Conn(i,j) - 2
            JK(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1) = ndim*Conn(i,j) - 1
            JK(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 + 1) = ndim*Conn(i,j)
            JK(K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1 - 3) = ndim*Conn(i,j) - 2
            JK(K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1 - 2) = ndim*Conn(i,j) - 1
            JK(K_pointer(i,2) + 2*ndim*K_pointer(i,1) + ndim*num1 - 1) = ndim*Conn(i,j)
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
!enforcing the fixed boundary conditions
num1 = SIZE(U_fix)
DO i = 1, num1
    nb = COUNT((conn(U_fix(i),:) .LE. U_fix(i)) .AND. (conn(U_fix(i),:) .NE. 0))
    DO j = 1, nb
        IF(conn(U_fix(i),j) .EQ. U_fix(i))THEN
            K_global(K_pointer(U_fix(i),2) : K_pointer(U_fix(i)+1,2) - 1) = 0.D0 ! zero out the entire row
            K_global(K_pointer(U_fix(i),2)) = 1.D0
            K_global(K_pointer(U_fix(i),2) + ndim*K_pointer(U_fix(i),1)) = 1.D0
            K_global(K_pointer(U_fix(i),2) + 2*ndim*K_pointer(U_fix(i),1) - 1) = 1.D0
        ELSE
            m = COUNT(conn(conn(U_fix(i),j),:) .NE. 0)
            n = COUNT(conn(conn(U_fix(i),j),:) .GT. conn(U_fix(i),j))
            num2 = 0
            DO k = m-n+1, m
                num2 = num2 + 1
                IF(conn(conn(U_fix(i),j),k) .EQ. U_fix(i))THEN
                    K_global(K_pointer(conn(U_fix(i),j),2) + ndim*num2:K_pointer(conn(U_fix(i),j),2) + ndim*num2+2) = 0.D0
                    K_global(K_pointer(conn(U_fix(i),j),2) + ndim*K_pointer(conn(U_fix(i),j),1) + ndim*num2 - 1 : K_pointer(conn(U_fix(i),j),2) + ndim*K_pointer(conn(U_fix(i),j),1) + ndim*num2 + 1) = 0.D0
                    K_global(K_pointer(conn(U_fix(i),j),2) + 2*ndim*K_pointer(conn(U_fix(i),j),1) + ndim*num2 - 3 : K_pointer(conn(U_fix(i),j),2) + 2*ndim*K_pointer(conn(U_fix(i),j),1) + ndim*num2 - 1) = 0.D0
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    P(ndim*U_fix(i) - 2) = 0.D0
    P(ndim*U_fix(i) - 1) = 0.D0
    P(ndim*U_fix(i) ) = 0.D0
ENDDO

!!enforcing the fixed boundary conditions: z component
!num1 = SIZE(U_fixz)
!DO i = 1, num1
!    nb = COUNT((conn(U_fixz(i),:) .LE. U_fixz(i)) .AND. (conn(U_fixz(i),:) .NE. 0))
!    DO j = 1, nb
!        IF(conn(U_fixz(i),j) .EQ. U_fixz(i))THEN
!            K_global(K_pointer(U_fixz(i),2) + 2) = 0.D0
!            K_global(K_pointer(U_fixz(i),2) + ndim*K_pointer(U_fixz(i),1) + 1) = 0.D0
!            K_global(K_pointer(U_fixz(i),2) + 2*ndim*K_pointer(U_fixz(i),1) - 1 : K_pointer(U_fixz(i)+1,2) - 1) = 0.D0 ! zero out the entire row
!            K_global(K_pointer(U_fixz(i),2) + 2*ndim*K_pointer(U_fixz(i),1) - 1) = 1.D0
!        ELSE
!            m = COUNT(conn(conn(U_fixz(i),j),:) .NE. 0)
!            n = COUNT(conn(conn(U_fixz(i),j),:) .GT. conn(U_fixz(i),j))
!            num2 = 0
!            DO k = m-n+1, m
!                num2 = num2 + 1
!                IF(conn(conn(U_fixz(i),j),k) .EQ. U_fixz(i))THEN
!                    K_global(K_pointer(conn(U_fixz(i),j),2) + ndim*num2 + 2) = 0.D0
!                    K_global(K_pointer(conn(U_fixz(i),j),2) + ndim*K_pointer(conn(U_fixz(i),j),1) + ndim*num2 + 1) = 0.D0
!                    K_global(K_pointer(conn(U_fixz(i),j),2) + 2*ndim*K_pointer(conn(U_fixz(i),j),1) + ndim*num2 - 1) = 0.D0
!                ENDIF
!            ENDDO
!        ENDIF
!    ENDDO
!    P(ndim*U_fixz(i) ) = 0.D0
!ENDDO

!enforcing the applied displacement boundary conditions
num1 = SIZE(U_appz)
DO i = 1, num1
    nb = COUNT(conn(U_appz(i),:) .NE. 0)
    num3 = 0
    DO j = 1, nb
        IF(conn(U_appz(i),j) .EQ. U_appz(i))THEN
            !modify the force vector
            P(ndim*U_appz(i) - 2) = P(ndim*U_appz(i) - 2) - U_stepz*K_global(K_pointer(U_appz(i),2) + 2)
            P(ndim*U_appz(i) - 1) = P(ndim*U_appz(i) - 1) - U_stepz*K_global(K_pointer(U_appz(i),2) + ndim*K_pointer(U_appz(i),1) + 1)
            P(ndim*U_appz(i)) = U_stepz
            !zero out the columns corresponding to the z component
            K_global(K_pointer(U_appz(i),2) + 2) = 0.D0
            K_global(K_pointer(U_appz(i),2) + ndim*K_pointer(U_appz(i),1) + 1) = 0.D0
            K_global(K_pointer(U_appz(i),2) + 2*ndim*K_pointer(U_appz(i),1) - 1) = 1.D0
        ELSEIF(conn(U_appz(i),j) .LT. U_appz(i))THEN
            m = COUNT(conn(conn(U_appz(i),j),:) .NE. 0)
            n = COUNT(conn(conn(U_appz(i),j),:) .GT. conn(U_appz(i),j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(conn(conn(U_appz(i),j),k) .EQ. U_appz(i))THEN
                    !modify the force vector
                    P(ndim*conn(U_appz(i),j) - 2) = P(ndim*conn(U_appz(i),j) - 2) - U_stepz*K_global(K_pointer(conn(U_appz(i),j),2) + ndim*num2 + 2)
                    P(ndim*conn(U_appz(i),j) - 1) = P(ndim*conn(U_appz(i),j) - 1) - U_stepz*K_global(K_pointer(conn(U_appz(i),j),2) + ndim*K_pointer(conn(U_appz(i),j),1) + ndim*num2 + 1)
                    P(ndim*conn(U_appz(i),j)) = P(ndim*conn(U_appz(i),j)) - U_stepz*K_global(K_pointer(conn(U_appz(i),j),2) + 2*ndim*K_pointer(conn(U_appz(i),j),1) + ndim*num2 - 1)
                    !zero out the columns corresponding to the z component
                    K_global(K_pointer(conn(U_appz(i),j),2) + ndim*num2 + 2) = 0.D0
                    K_global(K_pointer(conn(U_appz(i),j),2) + ndim*K_pointer(conn(U_appz(i),j),1) + ndim*num2 + 1) = 0.D0
                    K_global(K_pointer(conn(U_appz(i),j),2) + 2*ndim*K_pointer(conn(U_appz(i),j),1) + ndim*num2 - 1) = 0.D0
                ENDIF
            ENDDO
        ELSEIF(conn(U_appz(i),j) .GT. U_appz(i))THEN
            !modify the force vector
            P(ndim*conn(U_appz(i),j) - 2) = P(ndim*conn(U_appz(i),j) - 2) - U_stepz*K_global(K_pointer(U_appz(i),2) + 2*ndim*K_pointer(U_appz(i),1) + ndim*num3)
            P(ndim*conn(U_appz(i),j) - 1) = P(ndim*conn(U_appz(i),j) - 1) - U_stepz*K_global(K_pointer(U_appz(i),2) + 2*ndim*K_pointer(U_appz(i),1) + ndim*num3 + 1)
            P(ndim*conn(U_appz(i),j)) = P(ndim*conn(U_appz(i),j)) - U_stepz*K_global(K_pointer(U_appz(i),2) + 2*ndim*K_pointer(U_appz(i),1) + ndim*num3 + 2)
            ! zero out the row components correspoding to the z component
            K_global(K_pointer(U_appz(i),2) + 2*ndim*K_pointer(U_appz(i),1) + ndim*num3) = 0.D0
            K_global(K_pointer(U_appz(i),2) + 2*ndim*K_pointer(U_appz(i),1) + ndim*num3 + 1) = 0.D0
            K_global(K_pointer(U_appz(i),2) + 2*ndim*K_pointer(U_appz(i),1) + ndim*num3 + 2) = 0.D0
            
            num3 = num3 + 1
        ENDIF
    ENDDO
ENDDO

!!---------------------------------------NATURAL BOUDARY CONDITION---------------------------------
!!enforcing the applied force boundary conditions
!num1 = SIZE(F_appz)
!DO i = 1, num1
!    P(ndim*F_appy(i)  - 1) = P(ndim*F_appy(i)  - 1) + F_stepy ! y component
!ENDDO


END SUBROUTINE Get_K_P

!====================================================END Get_K_P====================================================

!====================================================FUNCTION du2dxyz===============================================

FUNCTION du2dxyz(i, j)
IMPLICIT NONE

INTEGER:: i, j, k, nb
DOUBLE PRECISION:: du2dxyz(9)

du2dxyz = 0.D0

nb = COUNT(neighbors1(i,:) .NE. 0)
DO k = 1, nb
    !d2udxidxi, 1 1
    du2dxyz(1) = du2dxyz(1) + (Kn1*csx(i,k,1) + Tv*SUM(csx(i,:,1)))*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(1. - csx(i,k,1)**2)/distance(i,k,1) + (Kn1 + Tv)*csx(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(1. - csx(i,k,1)**2)/distance(i,k,1)
    !d2udxidyi, 1 2
    du2dxyz(2) = du2dxyz(2) + (Kn1*csy(i,k,1) + Tv*SUM(csy(i,:,1)))*csx(i,k,1) - (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csy(i,k,1)*csx(i,k,1)/distance(i,k,1) + (Kn1 + Tv)*csy(i,k,1)*csx(i,k,1) - (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)* csy(i,k,1)/distance(i,k,1)
    !d2udxidzi, 1 3
    du2dxyz(3) = du2dxyz(3) + (Kn1*csz(i,k,1) + Tv*SUM(csz(i,:,1)))*csx(i,k,1) - (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csz(i,k,1)*csx(i,k,1)/distance(i,k,1) + (Kn1 + Tv)*csz(i,k,1)*csx(i,k,1) - (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)* csz(i,k,1)/distance(i,k,1)
    !d2udyidyi, 2 2
    du2dxyz(5) = du2dxyz(5) + (Kn1*csy(i,k,1) + Tv*SUM(csy(i,:,1)))*csy(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(1. - csy(i,k,1)**2)/distance(i,k,1) + (Kn1 + Tv)*csy(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(1. - csy(i,k,1)**2)/distance(i,k,1)
    !d2udyidzi, 2 3
    du2dxyz(6) = du2dxyz(6) + (Kn1*csz(i,k,1) + Tv*SUM(csz(i,:,1)))*csy(i,k,1) - (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csz(i,k,1)*csy(i,k,1)/distance(i,k,1) + (Kn1 + Tv)*csz(i,k,1)*csy(i,k,1) - (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csy(i,k,1)* csz(i,k,1)/distance(i,k,1)
    !d2udzidzi, 3 3
    du2dxyz(9) = du2dxyz(9) + (Kn1*csz(i,k,1) + Tv*SUM(csz(i,:,1)))*csz(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(1. - csz(i,k,1)**2)/distance(i,k,1) + (Kn1 + Tv)*csz(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(1. - csz(i,k,1)**2)/distance(i,k,1)
ENDDO

nb = COUNT(neighbors2(i,:) .NE. 0)
DO k = 1, nb
    !d2udxidxi, 1 1
    du2dxyz(1) = du2dxyz(1) + (Kn2*csx(i,k,2) + Tv*SUM(csx(i,:,2)))*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(1. - csx(i,k,2)**2)/distance(i,k,2) + (Kn2 + Tv)*csx(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(1. - csx(i,k,2)**2)/distance(i,k,2)
    !d2udxidyi, 1 2
    du2dxyz(2) = du2dxyz(2) + (Kn2*csy(i,k,2) + Tv*SUM(csy(i,:,2)))*csx(i,k,2) - (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csy(i,k,2)*csx(i,k,2)/distance(i,k,2) + (Kn2 + Tv)*csy(i,k,2)*csx(i,k,2) - (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)* csy(i,k,2)/distance(i,k,2)
    !d2udxidzi, 1 3
    du2dxyz(3) = du2dxyz(3) + (Kn2*csz(i,k,2) + Tv*SUM(csz(i,:,2)))*csx(i,k,2) - (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csz(i,k,2)*csx(i,k,2)/distance(i,k,2) + (Kn2 + Tv)*csz(i,k,2)*csx(i,k,2) - (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)* csz(i,k,2)/distance(i,k,2)
    !d2udyidyi, 2 2
    du2dxyz(5) = du2dxyz(5) + (Kn2*csy(i,k,2) + Tv*SUM(csy(i,:,2)))*csy(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(1. - csy(i,k,2)**2)/distance(i,k,2) + (Kn2 + Tv)*csy(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(1. - csy(i,k,2)**2)/distance(i,k,2)
    !d2udyidzi, 2 3
    du2dxyz(6) = du2dxyz(6) + (Kn2*csz(i,k,2) + Tv*SUM(csz(i,:,2)))*csy(i,k,2) - (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csz(i,k,2)*csy(i,k,2)/distance(i,k,2) + (Kn2 + Tv)*csz(i,k,2)*csy(i,k,2) - (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csy(i,k,2)* csz(i,k,2)/distance(i,k,2)
    !d2udzidzi, 3 3
    du2dxyz(9) = du2dxyz(9) + (Kn2*csz(i,k,2) + Tv*SUM(csz(i,:,2)))*csz(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(1. - csz(i,k,2)**2)/distance(i,k,2) + (Kn2 + Tv)*csz(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(1. - csz(i,k,2)**2)/distance(i,k,2)
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

!====================================================FUNCTION du2dxyz1===============================================

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
                du2dxyz1(1) = du2dxyz1(1) - (Kn1 + Tv)*csx(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(-1. + csx(i,k,1)**2)/distance(i,k,1) - (Kn1*csx(i,k,1) - Tv*SUM(csx(neighbors1(i,k),:,1)))*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(-1. + csx(i,k,1)**2)/distance(i,k,1)
                !d2udxidyj, 1 2
                du2dxyz1(2) = du2dxyz1(2) - (Kn1 + Tv)*csy(i,k,1)*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1) - (Kn1*csy(i,k,1) - Tv*SUM(csy(neighbors1(i,k),:,1)))*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1)
                !d2udxidzj, 1 3
                du2dxyz1(3) = du2dxyz1(3) - (Kn1 + Tv)*csz(i,k,1)*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csx(i,k,1)*csz(i,k,1)/distance(i,k,1) - (Kn1*csz(i,k,1) - Tv*SUM(csz(neighbors1(i,k),:,1)))*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csz(i,k,1)/distance(i,k,1)
                !d2udyidxj, 2 1
                du2dxyz1(4) = du2dxyz1(4) - (Kn1 + Tv)*csy(i,k,1)*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1) - (Kn1*csx(i,k,1) - Tv*SUM(csx(neighbors1(i,k),:,1)))*csy(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1)
                !d2udyidyj, 2 2
                du2dxyz1(5) = du2dxyz1(5) - (Kn1 + Tv)*csy(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(-1. + csy(i,k,1)**2)/distance(i,k,1) - (Kn1*csy(i,k,1) - Tv*SUM(csy(neighbors1(i,k),:,1)))*csy(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(-1. + csy(i,k,1)**2)/distance(i,k,1)
                !d2udyidzj, 2 3
                du2dxyz1(6) = du2dxyz1(6) - (Kn1 + Tv)*csy(i,k,1)*csz(i,k,1)+ (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csz(i,k,1)*csy(i,k,1)/distance(i,k,1) - (Kn1*csz(i,k,1) - Tv*SUM(csz(neighbors1(i,k),:,1)))*csy(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csz(i,k,1)*csy(i,k,1)/distance(i,k,1)
                !d2udzidxj, 3 1
                du2dxyz1(7) = du2dxyz1(7) - (Kn1 + Tv)*csz(i,k,1)*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csx(i,k,1)*csz(i,k,1)/distance(i,k,1) - (Kn1*csx(i,k,1) - Tv*SUM(csx(neighbors1(i,k),:,1)))*csz(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csz(i,k,1)/distance(i,k,1)
                !d2udzidyj, 3 2
                du2dxyz1(8) = du2dxyz1(8) - (Kn1 + Tv)*csz(i,k,1)*csy(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csy(i,k,1)*csz(i,k,1)/distance(i,k,1) - (Kn1*csy(i,k,1) - Tv*SUM(csy(neighbors1(i,k),:,1)))*csz(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csy(i,k,1)*csz(i,k,1)/distance(i,k,1)
                !d2udzidzj, 3 3
                du2dxyz1(9) = du2dxyz1(9) - (Kn1 + Tv)*csz(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(-1. + csz(i,k,1)**2)/distance(i,k,1) - (Kn1*csz(i,k,1) - Tv*SUM(csz(neighbors1(i,k),:,1)))*csz(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(-1. + csz(i,k,1)**2)/distance(i,k,1)
                !-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            ELSE
                !----------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxyz1(1) = du2dxyz1(1) - Tv*csx(i,kk,1)*csx(i,k,1)
                !d2udxidyj, 1 2
                du2dxyz1(2) = du2dxyz1(2) - Tv*csy(i,kk,1)*csx(i,k,1)
                !d2udxidzj, 1 3
                du2dxyz1(3) = du2dxyz1(3) - Tv*csz(i,kk,1)*csx(i,k,1)
                !d2udyidxj, 2 1
                du2dxyz1(4) = du2dxyz1(4) - Tv*csx(i,kk,1)*csy(i,k,1)
                !d2udyidyj, 2 2
                du2dxyz1(5) = du2dxyz1(5) - Tv*csy(i,kk,1)*csy(i,k,1)
                !d2udyidzj, 2 3
                du2dxyz1(6) = du2dxyz1(6) - Tv*csz(i,kk,1)*csy(i,k,1)
                !d2udzidxj, 3 1
                du2dxyz1(7) = du2dxyz1(7) - Tv*csx(i,kk,1)*csz(i,k,1)
                !d2udzidyj, 3 2
                du2dxyz1(8) = du2dxyz1(8) - Tv*csy(i,kk,1)*csz(i,k,1)
                !d2udzidzj, 3 3
                du2dxyz1(9) = du2dxyz1(9) - Tv*csz(i,kk,1)*csz(i,k,1)
                !----------------------------------------------------------------------------------
                nb2 = COUNT(neighbors1(neighbors1(i,k),:) .NE. 0)
                DO n = 1,nb2
                    IF(neighbors1(neighbors1(i,k),n) .EQ. Conn(i,j))THEN
                        !---------------------------------------------------------------------------------------------------
                        !d2udxidxj, 1 1
                        du2dxyz1(1) = du2dxyz1(1) - Tv*csx(neighbors1(i,k),n,1)*csx(i,k,1)
                        !d2udxidyj, 1 2
                        du2dxyz1(2) = du2dxyz1(2) - Tv*csy(neighbors1(i,k),n,1)*csx(i,k,1)
                        !d2udxidzj, 1 3
                        du2dxyz1(3) = du2dxyz1(3) - Tv*csz(neighbors1(i,k),n,1)*csx(i,k,1)
                        !d2udyidxj, 2 1
                        du2dxyz1(4) = du2dxyz1(4) - Tv*csx(neighbors1(i,k),n,1)*csy(i,k,1)
                        !d2udyidyj, 2 2
                        du2dxyz1(5) = du2dxyz1(5) - Tv*csy(neighbors1(i,k),n,1)*csy(i,k,1)
                        !d2udyidzj, 2 3
                        du2dxyz1(6) = du2dxyz1(6) - Tv*csz(neighbors1(i,k),n,1)*csy(i,k,1)
                        !d2udzidxj, 3 1
                        du2dxyz1(7) = du2dxyz1(7) - Tv*csx(neighbors1(i,k),n,1)*csz(i,k,1)
                        !d2udzidyj, 3 2
                        du2dxyz1(8) = du2dxyz1(8) - Tv*csy(neighbors1(i,k),n,1)*csz(i,k,1)
                        !d2udzidzj, 3 3
                        du2dxyz1(9) = du2dxyz1(9) - Tv*csz(neighbors1(i,k),n,1)*csz(i,k,1)
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
            du2dxyz1(1) = du2dxyz1(1) - Tv*csx(i,kk,1)*csx(neighbors1(i,kk),m,1)
            !d2udxidyj, 1 2
            du2dxyz1(2) = du2dxyz1(2) - Tv*csx(i,kk,1)*csy(neighbors1(i,kk),m,1)
            !d2udxidzj, 1 3
            du2dxyz1(3) = du2dxyz1(3) - Tv*csx(i,kk,1)*csz(neighbors1(i,kk),m,1)
            !d2udyidxj, 2 1
            du2dxyz1(4) = du2dxyz1(4) - Tv*csy(i,kk,1)*csx(neighbors1(i,kk),m,1)
            !d2udyidyj, 2 2
            du2dxyz1(5) = du2dxyz1(5) - Tv*csy(i,kk,1)*csy(neighbors1(i,kk),m,1)
            !d2udyidzj, 2 3
            du2dxyz1(6) = du2dxyz1(6) - Tv*csy(i,kk,1)*csz(neighbors1(i,kk),m,1)
            !d2udzidxj, 3 1
            du2dxyz1(7) = du2dxyz1(7) - Tv*csz(i,kk,1)*csx(neighbors1(i,kk),m,1)
            !d2udzidyj, 3 2
            du2dxyz1(8) = du2dxyz1(8) - Tv*csz(i,kk,1)*csy(neighbors1(i,kk),m,1)
            !d2udzidzj, 3 3
            du2dxyz1(9) = du2dxyz1(9) - Tv*csz(i,kk,1)*csz(neighbors1(i,kk),m,1)
            !---------------------------------------------------------------------------------------------------------
        ENDIF
    ENDDO
ENDDO

777   RETURN
END FUNCTION du2dxyz1

!====================================================END du2dxyz1====================================================

!====================================================FUNCTION du2dxyz2===============================================

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
                du2dxyz2(1) = du2dxyz2(1) - (Kn2 + Tv)*csx(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(-1. + csx(i,k,2)**2)/distance(i,k,2) - (Kn2*csx(i,k,2) - Tv*SUM(csx(neighbors2(i,k),:,2)))*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(-1. + csx(i,k,2)**2)/distance(i,k,2)
                !d2udxidyj, 1 2
                du2dxyz2(2) = du2dxyz2(2) - (Kn2 + Tv)*csy(i,k,2)*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2) - (Kn2*csy(i,k,2) - Tv*SUM(csy(neighbors2(i,k),:,2)))*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2)
                !d2udxidzj, 1 3
                du2dxyz2(3) = du2dxyz2(3) - (Kn2 + Tv)*csz(i,k,2)*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csx(i,k,2)*csz(i,k,2)/distance(i,k,2) - (Kn2*csz(i,k,2) - Tv*SUM(csz(neighbors2(i,k),:,2)))*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csz(i,k,2)/distance(i,k,2)
                !d2udyidxj, 2 1
                du2dxyz2(4) = du2dxyz2(4) - (Kn2 + Tv)*csy(i,k,2)*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2) - (Kn2*csx(i,k,2) - Tv*SUM(csx(neighbors2(i,k),:,2)))*csy(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2)
                !d2udyidyj, 2 2
                du2dxyz2(5) = du2dxyz2(5) - (Kn2 + Tv)*csy(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(-1. + csy(i,k,2)**2)/distance(i,k,2) - (Kn2*csy(i,k,2) - Tv*SUM(csy(neighbors2(i,k),:,2)))*csy(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(-1. + csy(i,k,2)**2)/distance(i,k,2)
                !d2udyidzj, 2 3
                du2dxyz2(6) = du2dxyz2(6) - (Kn2 + Tv)*csy(i,k,2)*csz(i,k,2)+ (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csz(i,k,2)*csy(i,k,2)/distance(i,k,2) - (Kn2*csz(i,k,2) - Tv*SUM(csz(neighbors2(i,k),:,2)))*csy(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csz(i,k,2)*csy(i,k,2)/distance(i,k,2)
                !d2udzidxj, 3 1
                du2dxyz2(7) = du2dxyz2(7) - (Kn2 + Tv)*csz(i,k,2)*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csx(i,k,2)*csz(i,k,2)/distance(i,k,2) - (Kn2*csx(i,k,2) - Tv*SUM(csx(neighbors2(i,k),:,2)))*csz(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csz(i,k,2)/distance(i,k,2)
                !d2udzidyj, 3 2
                du2dxyz2(8) = du2dxyz2(8) - (Kn2 + Tv)*csz(i,k,2)*csy(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csy(i,k,2)*csz(i,k,2)/distance(i,k,2) - (Kn2*csy(i,k,2) - Tv*SUM(csy(neighbors2(i,k),:,2)))*csz(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csy(i,k,2)*csz(i,k,2)/distance(i,k,2)
                !d2udzidzj, 3 3
                du2dxyz2(9) = du2dxyz2(9) - (Kn2 + Tv)*csz(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(-1. + csz(i,k,2)**2)/distance(i,k,2) - (Kn2*csz(i,k,2) - Tv*SUM(csz(neighbors2(i,k),:,2)))*csz(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(-1. + csz(i,k,2)**2)/distance(i,k,2)
                !-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            ELSE
                !----------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxyz2(1) = du2dxyz2(1) - Tv*csx(i,kk,2)*csx(i,k,2)
                !d2udxidyj, 1 2
                du2dxyz2(2) = du2dxyz2(2) - Tv*csy(i,kk,2)*csx(i,k,2)
                !d2udxidzj, 1 3
                du2dxyz2(3) = du2dxyz2(3) - Tv*csz(i,kk,2)*csx(i,k,2)
                !d2udyidxj, 2 1
                du2dxyz2(4) = du2dxyz2(4) - Tv*csx(i,kk,2)*csy(i,k,2)
                !d2udyidyj, 2 2
                du2dxyz2(5) = du2dxyz2(5) - Tv*csy(i,kk,2)*csy(i,k,2)
                !d2udyidzj, 2 3
                du2dxyz2(6) = du2dxyz2(6) - Tv*csz(i,kk,2)*csy(i,k,2)
                !d2udzidxj, 3 1
                du2dxyz2(7) = du2dxyz2(7) - Tv*csx(i,kk,2)*csz(i,k,2)
                !d2udzidyj, 3 2
                du2dxyz2(8) = du2dxyz2(8) - Tv*csy(i,kk,2)*csz(i,k,2)
                !d2udzidzj, 3 3
                du2dxyz2(9) = du2dxyz2(9) - Tv*csz(i,kk,2)*csz(i,k,2)
                !----------------------------------------------------------------------------------
                nb2 = COUNT(neighbors2(neighbors2(i,k),:) .NE. 0)
                DO n = 1,nb2
                    IF(neighbors2(neighbors2(i,k),n) .EQ. Conn(i,j))THEN
                        !---------------------------------------------------------------------------------------------------
                        !d2udxidxj, 1 1
                        du2dxyz2(1) = du2dxyz2(1) - Tv*csx(neighbors2(i,k),n,2)*csx(i,k,2)
                        !d2udxidyj, 1 2
                        du2dxyz2(2) = du2dxyz2(2) - Tv*csy(neighbors2(i,k),n,2)*csx(i,k,2)
                        !d2udxidzj, 1 3
                        du2dxyz2(3) = du2dxyz2(3) - Tv*csz(neighbors2(i,k),n,2)*csx(i,k,2)
                        !d2udyidxj, 2 1
                        du2dxyz2(4) = du2dxyz2(4) - Tv*csx(neighbors2(i,k),n,2)*csy(i,k,2)
                        !d2udyidyj, 2 2
                        du2dxyz2(5) = du2dxyz2(5) - Tv*csy(neighbors2(i,k),n,2)*csy(i,k,2)
                        !d2udyidzj, 2 3
                        du2dxyz2(6) = du2dxyz2(6) - Tv*csz(neighbors2(i,k),n,2)*csy(i,k,2)
                        !d2udzidxj, 3 1
                        du2dxyz2(7) = du2dxyz2(7) - Tv*csx(neighbors2(i,k),n,2)*csz(i,k,2)
                        !d2udzidyj, 3 2
                        du2dxyz2(8) = du2dxyz2(8) - Tv*csy(neighbors2(i,k),n,2)*csz(i,k,2)
                        !d2udzidzj, 3 3
                        du2dxyz2(9) = du2dxyz2(9) - Tv*csz(neighbors2(i,k),n,2)*csz(i,k,2)
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
            du2dxyz2(1) = du2dxyz2(1) - Tv*csx(i,kk,2)*csx(neighbors2(i,kk),m,2)
            !d2udxidyj, 1 2
            du2dxyz2(2) = du2dxyz2(2) - Tv*csx(i,kk,2)*csy(neighbors2(i,kk),m,2)
            !d2udxidzj, 1 3
            du2dxyz2(3) = du2dxyz2(3) - Tv*csx(i,kk,2)*csz(neighbors2(i,kk),m,2)
            !d2udyidxj, 2 1
            du2dxyz2(4) = du2dxyz2(4) - Tv*csy(i,kk,2)*csx(neighbors2(i,kk),m,2)
            !d2udyidyj, 2 2
            du2dxyz2(5) = du2dxyz2(5) - Tv*csy(i,kk,2)*csy(neighbors2(i,kk),m,2)
            !d2udyidzj, 2 3
            du2dxyz2(6) = du2dxyz2(6) - Tv*csy(i,kk,2)*csz(neighbors2(i,kk),m,2)
            !d2udzidxj, 3 1
            du2dxyz2(7) = du2dxyz2(7) - Tv*csz(i,kk,2)*csx(neighbors2(i,kk),m,2)
            !d2udzidyj, 3 2
            du2dxyz2(8) = du2dxyz2(8) - Tv*csz(i,kk,2)*csy(neighbors2(i,kk),m,2)
            !d2udzidzj, 3 3
            du2dxyz2(9) = du2dxyz2(9) - Tv*csz(i,kk,2)*csz(neighbors2(i,kk),m,2)
            !---------------------------------------------------------------------------------------------------------
        ENDIF
    ENDDO
ENDDO

888   RETURN
END FUNCTION du2dxyz2

!====================================================END du2dxyz2====================================================

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
    dudxyz(1) = dudxyz(1) + (2.*Kn1*dL(i,j,1) + Tv*(SUM(dL(i,:,1)) + SUM(dL(neighbors1(i,j),:,1))))*csx(i,j,1)
    !dudy
    dudxyz(2) = dudxyz(2) + (2.*Kn1*dL(i,j,1) + Tv*(SUM(dL(i,:,1)) + SUM(dL(neighbors1(i,j),:,1))))*csy(i,j,1)
    !dudz
    dudxyz(3) = dudxyz(3) + (2.*Kn1*dL(i,j,1) + Tv*(SUM(dL(i,:,1)) + SUM(dL(neighbors1(i,j),:,1))))*csz(i,j,1)
ENDDO

nb = COUNT(neighbors2(i,:) .NE. 0)
DO j = 1,nb
    !dudx
    dudxyz(1) = dudxyz(1) + (2.*Kn2*dL(i,j,2) + Tv*(SUM(dL(i,:,2)) + SUM(dL(neighbors2(i,j),:,2))))*csx(i,j,2)
    !dudy
    dudxyz(2) = dudxyz(2) + (2.*Kn2*dL(i,j,2) + Tv*(SUM(dL(i,:,2)) + SUM(dL(neighbors2(i,j),:,2))))*csy(i,j,2)
    !dudz
    dudxyz(3) = dudxyz(3) + (2.*Kn2*dL(i,j,2) + Tv*(SUM(dL(i,:,2)) + SUM(dL(neighbors2(i,j),:,2))))*csz(i,j,2)
ENDDO

RETURN
END FUNCTION dudxyz

!====================================================END dudxyz======================================================

!==============================================SUBROUTINE CGSolver===================================================

SUBROUTINE CGSolver
IMPLICIT NONE

INTEGER:: RCI_request, itercount, ipar(128)
DOUBLE PRECISION::  dpar(128), du(ndim*nparticle), TMP(ndim*nparticle,4)

! The initial guess
du = 0.D0

! Initialize the solver
CALL dcg_init(ndim*nparticle, du, P, RCI_request, ipar, dpar, TMP)
IF (RCI_request .NE. 0 ) GOTO 10002

! Modify the initialized solver parameters
ipar(5) = ndim*nparticle ! specifies the maximum number of iterations. The default value is min(150, ndim*nparticle)
ipar(9) = 1 !default value is 0, does not perform the residual stopping test; otherwise, perform the test
ipar(10) = 0 !default value is 1, perform user defined stopping test; otherwise, does not perform the test
dpar(1) = 1.D-20 !specifies the relative tolerance. The default value is 1.0D-6
dpar(2) = 1.D-20 !specifies the absolute tolerance. The default value is 0.0.

! Check the correctness and consistency of the newly set parameters
CALL dcg_check(ndim*nparticle, du, P, RCI_request, ipar, dpar, TMP)
    
IF (RCI_request .NE. 0 ) GOTO 10002
        
! Compute the solution by RCI (P)CG solver without preconditioning
! Reverse Communications starts here
!---------------------------------------------------------------------------
10000   CALL dcg(ndim*nparticle, du, P, RCI_request, ipar, dpar, TMP)
!---------------------------------------------------------------------------
! If RCI_request=0, then the solution was found with the required precision
!---------------------------------------------------------------------------
IF (RCI_request .EQ. 0) THEN
    GOTO 10001
!---------------------------------------------------------------------------
! If RCI_request=1, then compute the vector A*TMP(:,1) and put the result in vector TMP(:,2)
!---------------------------------------------------------------------------
ELSEIF (RCI_request .EQ. 1) THEN
    CALL MKL_DCSRSYMV('U', ndim*nparticle, K_global, IK, JK, TMP(:,1), TMP(:,2))
    GOTO 10000
ELSE
!---------------------------------------------------------------------------
! If RCI_request=anything else, then dcg subroutine failed to compute the solution vector: solution(N)
!---------------------------------------------------------------------------
    GOTO 10002
ENDIF
!---------------------------------------------------------------------------
! Reverse Communication ends here, get the current iteration number
!---------------------------------------------------------------------------
10001 CALL dcg_get(ndim*nparticle, du, P, RCI_request, ipar, dpar, TMP, itercount)
!---------------------------------------------------------------------------
WRITE(*, '(A,I8)') ' The system has been solved with total number of iteration: ', itercount

!components of position
Positions(:,1) = Positions(:,1) + du(1:ndim*nparticle:3)
Positions(:,2) = Positions(:,2) + du(2:ndim*nparticle:3)
Positions(:,3) = Positions(:,3) + du(3:ndim*nparticle:3)

GOTO 10003

!---------------------------------------------------------------------------
10002 WRITE( *,'(A,I2)') 'The solution trial was FAILED with RCI_request =  ', RCI_request
10003 CONTINUE
        
END SUBROUTINE CGSolver

!==============================================SUBROUTINE CGSolver===================================================

!==============================================SUBROUTINE Output=====================================================

SUBROUTINE Output
IMPLICIT NONE

INTEGER:: i

!Save the displacement to file
OPEN(UNIT = 13, FILE = 'initialP.dump', ACTION = 'WRITE', STATUS = 'REPLACE', POSITION = 'APPEND')
WRITE(13, '(A)') 'ITEM: TIMESTEP'
WRITE(13, '(I0)') 0
WRITE(13, '(A)') 'ITEM: NUMBER OF ATOMS'
WRITE(13, '(I0)') nparticle
WRITE(13, '(A)') 'ITEM: BOX BOUNDS pp pp pp'
WRITE(13, '(2F4.1)') 100*Box(1, 1), 100*Box(1, 2)
WRITE(13, '(2F4.1)') 100*Box(2, 1), 100*Box(2, 2)
WRITE(13, '(2F4.1)') 100*Box(3, 1), 100*Box(3, 2)
WRITE(13, '(A)') 'ITEM: ATOMS id x y z '
DO i = 1, nparticle
    WRITE(13, '(I0, 3F8.4)') i, 100*initialP(i, 1), 100*initialP(i, 2), 100*initialP(i, 3)
ENDDO
CLOSE(UNIT = 13)


OPEN(UNIT=13, FILE="finalP.dump", ACTION="write", STATUS="replace", POSITION='append')
WRITE(13,"(A)") "ITEM: TIMESTEP"
WRITE(13,"(I0)") 1
WRITE(13,"(A)") "ITEM: NUMBER OF ATOMS"
WRITE(13,"(I0)") nparticle
WRITE(13,"(A)") "ITEM: BOX BOUNDS pp pp pp"
WRITE(13,"(2F8.4)") 100*Box(1,1), 100*Box(1,2)
WRITE(13,"(2F8.4)") 100*Box(2,1), 100*Box(2,2)
WRITE(13,"(2F8.4)") 100*Box(3,1), 100*Box(3,2)
WRITE(13,"(A)") "ITEM: ATOMS id x y z dx dy dz"
DO i = 1, nparticle
    WRITE(13,"(I0, 3F12.4, 3E12.4E2)") i, 100*Positions(i,1), 100*Positions(i,2), 100*Positions(i,3), 100*(Positions(i,1)-initialP(i,1)), 100*(Positions(i,2)-initialP(i,2)), 100*(Positions(i,3)-initialP(i,3))
ENDDO
CLOSE(UNIT=13)


END SUBROUTINE Output

!==============================================END Output===========================================================

END PROGRAM MAIN