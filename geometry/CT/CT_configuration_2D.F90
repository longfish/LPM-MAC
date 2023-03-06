!INCLUDE 'mkl_pardiso.f90'
PROGRAM MAIN
USE mkl_pardiso
USE OMP_lib
USE IFPORT
IMPLICIT NONE

INTEGER:: t, ndim, particles_per_row, nparticle, rows, steps, nneighbors, nneighbors1, nneighbors2, nneighbors_afem, nneighbors_afem1, nneighbors_afem2
INTEGER, ALLOCATABLE, DIMENSION(:,:):: K_pointer, neighbors, neighbors1, neighbors2, Conn, nsign, neighbors_afem, neighbors_afem1, neighbors_afem2
INTEGER, ALLOCATABLE, DIMENSION(:):: IK, JK, U_fix, U_appy, F_appy
DOUBLE PRECISION:: Box(2,2), Radius, U_stepy, F_stepy, Kn1, Kn2, Tv, start, finish
DOUBLE PRECISION, PARAMETER :: E = 6.9D10, mu = 0.3
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:):: distance, origindistance, dL, csx, csy
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: Positions
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:):: K_global, P

!Dimensionality of the problem, 2D
ndim = 2
!The dimension of the physical domain using the format: [Xmin Xmax; Ymin Ymax]
Box = RESHAPE((/ 0., 0., 0.01, 0.03 /), (/ 2, 2 /))
!Number of particles in the first row
particles_per_row = 31
!Number of neighbors
nneighbors = 8
nneighbors1 = 4
nneighbors2 = 4
nneighbors_afem = 16
nneighbors_afem1 = 12
nneighbors_afem2 = 12

!-------------------------------------------------------------------------------------------
!record the CPU time
CALL CPU_TIME(start)
!initialize the particle position and allocate memory for the global matrices
!----------------------------------------------------------------------------------------------------
CALL  Initialization
!----------------------------------------------------------------------------------------------------
!plane strain
Kn1 = E/(2.*(1. + mu))
Kn2 = E/(4.*(1. + mu))
Tv = E*(4.*mu - 1.)/(24.*(1. + mu)*(1. - 2.*mu))
!!plane stress
!Kn1 = E/(2.*(1. + mu))
!Kn2 = E/(4.*(1. + mu))
!Tv = E*(3.*mu - 1.)/(24.*(1. + mu)*(1. - mu))
!---------------------------------------------------------------------------------------------------
!Save the initial positions to file
!-------------------------------------------------------------------------------------------------------------------------------------------
OPEN(UNIT=13, FILE="D:\FortranLPM\VCPM_2D_elasticity_afem_sqr_2\VCPM_2D_elasticity_afem_sqr_2\Results\initialP.dat", ACTION="write", STATUS="replace", POSITION='append')
DO t = 1, nparticle
    WRITE(13,"(I0, I6, 2F15.10)") t, 1, Positions(t,1), Positions(t,2)
ENDDO
CLOSE(UNIT=13)
!-------------------------------------------------------------------------------------------------------------------------------------------
!Search the neighbors for each particle at the initial configuration
!----------------------------------------------------------------------------------------
CALL  Searchneighbor
!----------------------------------------------------------------------------------------
!applying the boundary condition using different steps
!---------------------------------------------------------------------------------------

U_stepy = -1.D-5
!F_stepy = -2.0D3/particles_per_row

CALL PARDISOsolver

!-------------------------------------------------------------------------------------------------------------------------------------------
!Save the displacement to file
OPEN(UNIT=13, FILE="D:\FortranLPM\VCPM_2D_elasticity_afem_sqr_2\VCPM_2D_elasticity_afem_sqr_2\Results\finalP.dat", ACTION="write", STATUS="replace", POSITION='append')
DO t = 1, nparticle
    WRITE(13,"(I0, I6, 2F15.10)") t, 1, Positions(t,1), Positions(t,2)
ENDDO
CLOSE(UNIT=13)

!record the CPU time
CALL CPU_TIME(finish)
WRITE(*,*) 'Total computation time: ', finish-start, 'seconds.'
!--------------------------------------------------------------------------------------------------------------------------------------------

CALL MKL_FREE_BUFFERS

!SUBROUTINES
CONTAINS


!==================================================SUBROUTINE Initialization============================================

SUBROUTINE  Initialization
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

!specify the boundary particles
!---------------------------------------------------------------------------------------------------
ALLOCATE(U_fix(particles_per_row), U_appy(particles_per_row), F_appy(particles_per_row))
U_fix = (/ nparticle - particles_per_row + 1:nparticle /)
U_appy = (/ 1: particles_per_row /)
F_appy = U_appy

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
ALLOCATE(dL(nparticle,nneighbors1,2))
dL = 0.D0
ALLOCATE(neighbors(nparticle,nneighbors))
neighbors = 0
ALLOCATE(neighbors1(nparticle,nneighbors1))
neighbors1 = 0
ALLOCATE(neighbors2(nparticle,nneighbors2))
neighbors2 = 0
ALLOCATE(neighbors_afem(nparticle,nneighbors_afem))
neighbors_afem = 0
ALLOCATE(neighbors_afem1(nparticle,nneighbors_afem1))
neighbors_afem1 = 0
ALLOCATE(neighbors_afem2(nparticle,nneighbors_afem2))
neighbors_afem2 = 0
ALLOCATE(K_pointer(nparticle+1, 2))
K_pointer = 0
ALLOCATE(conn(nparticle, nneighbors_afem + 1))
Conn = 0
ALLOCATE(nsign(nparticle,nneighbors))
nsign = 0
ALLOCATE(P(ndim*nparticle))
P = 0.D0

END SUBROUTINE Initialization

!================================================END Initialization======================================================

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
        dis = SQRT((Positions(j,1) - Positions(i,1))**2 + (Positions(j,2) - Positions(i,2))**2)
        IF(j .NE. i .AND. dis  .LT. 2.01*Radius)THEN !the first nearest neighbors
            index1 = index1 + 1
            neighbors1(i,index1) = j
            origindistance(i,index1,1) = dis
            
            index3 = index3 + 1
            neighbors(i,index3) = j
            nsign(i,index3) = 1
        ELSEIF(dis  .LT. 2.01*SQRT(2.)*Radius .AND. dis  .GT. 2.01*Radius)THEN !the second nearest neighbors
            index2 = index2 + 1
            neighbors2(i,index2) = j
            origindistance(i,index2,2) = dis
            
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

!--------------------------------------------------------------------------------------neighbors_afem---------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------Conn------------------------------------------------------------------------
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

!sorting the conn matrix
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
    K_pointer(i + 1, 2) = K_pointer(i, 2) + ndim*ndim*index1 - 1 !start index for each particle in the global stiffness matrix
ENDDO

!initialize the sparse stiffness matrix, using SCR storing format for the sparse matrix

ALLOCATE(JK(K_pointer(nparticle + 1, 2) - 1))
ALLOCATE(IK(ndim*nparticle + 1))
ALLOCATE(K_global(K_pointer(nparticle + 1, 2) - 1))

END SUBROUTINE Searchneighbor

!=============================================END Searchneighbor=====================================================

!====================================================SUBROUTINE Get_K_P============================================

SUBROUTINE Get_K_P
IMPLICIT NONE
INTEGER:: i, j, k, m, n, nb, num1, num2, num3
DOUBLE PRECISION:: dis, K_local(ndim*ndim)
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
dL = 0.D0

!$OMP PARALLEL DO PRIVATE(nb, num1, num2)
DO i = 1, nparticle
    nb = COUNT(neighbors(i,:) .NE. 0)
    num1 = 0
    num2 = 0
    DO j = 1, nb
        IF(nsign(i,j) .EQ. 1)THEN
            num1 = num1 + 1
            distance(i,num1,1)  = SQRT((Positions(i,1) - Positions(neighbors1(i,num1),1))**2 + (Positions(i,2) - Positions(neighbors1(i,num1),2))**2)
            csx(i,num1,1) = (Positions(i,1) - Positions(neighbors1(i,num1),1))/distance(i,num1,1)
            csy(i,num1,1) = (Positions(i,2) - Positions(neighbors1(i,num1),2))/distance(i,num1,1)
            dL(i,num1,1) =  distance(i,num1,1) - origindistance(i,num1,1)
        ELSEIF(nsign(i,j) .EQ. 2)THEN
            num2 = num2 + 1
            distance(i,num2,2)  = SQRT((Positions(i,1) - Positions(neighbors2(i,num2),1))**2 + (Positions(i,2) - Positions(neighbors2(i,num2),2))**2)
            csx(i,num2,2) = (Positions(i,1) - Positions(neighbors2(i,num2),1))/distance(i,num2,2)
            csy(i,num2,2) = (Positions(i,2) - Positions(neighbors2(i,num2),2))/distance(i,num2,2)
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
    nb = COUNT(Conn(i,:) .NE. 0)
    DO j = 1,nb
        IF(Conn(i,j) .EQ. i )THEN
            K_local = 0.5D0*du2dxy(i,j)
        ELSEIF(ANY(neighbors_afem1(i,:) .EQ. Conn(i,j)) .AND. ANY(neighbors_afem2(i,:) .EQ. Conn(i,j)))THEN
            K_local = 0.5D0*du2dxy1(i,j) + 0.5D0*du2dxy2(i,j)
        ELSEIF(ANY(neighbors_afem1(i,:) .EQ. Conn(i,j)))THEN
            K_local = 0.5D0*du2dxy1(i,j)
        ELSEIF(ANY(neighbors_afem2(i,:) .EQ. Conn(i,j)))THEN
            K_local = 0.5D0*du2dxy2(i,j)
        ENDIF
        !--------------------------------------------------------------------------------------------------------------row assembling------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        IF(Conn(i,j) .EQ. i)THEN
            K_global(K_pointer(i,2) : K_pointer(i,2) + 1) = K_global(K_pointer(i,2) : K_pointer(i,2) + 1) + 2.*K_local(1:2)
            K_global(K_pointer(i,2) + ndim*K_pointer(i,1)) = K_global(K_pointer(i,2) + ndim*K_pointer(i,1)) + 2.*K_local(4)
        ELSEIF(Conn(i,j) .GT. i)THEN
            num1 = num1 + 1
            K_global(K_pointer(i,2) + ndim*num1 : K_pointer(i,2) + ndim*num1 + 1) = K_global(K_pointer(i,2) + ndim*num1 : K_pointer(i,2) + ndim*num1 + 1) + K_local(1:2)
            K_global(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 - 1 : K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1) = K_global(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 - 1 : K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1) + K_local(3:4)
        ELSEIF(Conn(i,j) .LT. i)THEN
            !-------------------------------------------------------------------------------------------------------column assembling---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            m = COUNT(Conn(Conn(i,j),:) .NE. 0)
            n = COUNT(Conn(Conn(i,j),:) .GT. Conn(i,j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(conn(conn(i,j),k) .EQ. i)THEN
                    K_global(K_pointer(Conn(i,j),2) + ndim*num2 : K_pointer(Conn(i,j),2) + ndim*num2 + 1) = K_global(K_pointer(Conn(i,j),2) + ndim*num2 : K_pointer(Conn(i,j),2) + ndim*num2 + 1) + K_local(1:3:2)
                    K_global(K_pointer(Conn(i,j),2) + ndim*K_pointer(Conn(i,j),1) + ndim*num2 - 1 : K_pointer(Conn(i,j),2) + ndim*K_pointer(Conn(i,j),1) + ndim*num2) = K_global(K_pointer(Conn(i,j),2) + ndim*K_pointer(Conn(i,j),1) + ndim*num2 - 1 : K_pointer(Conn(i,j),2) + ndim*K_pointer(Conn(i,j),1) + ndim*num2) + K_local(2:4:2)
                ENDIF
            ENDDO
        ENDIF
        !the JK array
        !------------------------------------------------------------------------------------------
        IF(Conn(i,j) .EQ. i)THEN
            JK(K_pointer(i,2)) =  ndim*Conn(i,j) - 1
            JK(K_pointer(i,2) + 1) =  ndim*Conn(i,j)
            JK(K_pointer(i,2) + ndim*K_pointer(i,1)) = ndim*Conn(i,j)
        ELSEIF(Conn(i,j) .GT. i)THEN
            JK(K_pointer(i,2) + ndim*num1) =  ndim*Conn(i,j) - 1
            JK(K_pointer(i,2) + ndim*num1 + 1) =  ndim*Conn(i,j)
            JK(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1 - 1) = ndim*Conn(i,j) - 1
            JK(K_pointer(i,2) + ndim*K_pointer(i,1) + ndim*num1) = ndim*Conn(i,j)
        ENDIF
        !-----------------------------------------------------------------------------------------
    ENDDO
    !the IK array
    !------------------------------------------------------------------------
    IK(ndim*i-1) = K_pointer(i,2)
    IK(ndim*i) = K_pointer(i,2) + ndim*K_pointer(i,1)
    !------------------------------------------------------------------------
    !the force vector
    P(ndim*i - 1 : ndim*i) = - dudxy(i)
ENDDO
!$OMP END PARALLEL DO
IK(ndim*nparticle + 1) = K_pointer(nparticle + 1,2)

do i=1,100
    write(*,*) K_global(i)
enddo

!-------------------------------------------------------------------------------------------------------------------
!---------------------------------------NATURAL BOUDARY CONDITION---------------------------------
!enforcing the applied force boundary conditions
!num1 = SIZE(F_appy)
!DO i = 1, num1
!    P(ndim*F_appy(i)  - 1) = P(ndim*F_appy(i)  - 1) + F_stepy ! y component
!ENDDO
!----------------------------------ESSENTIAL BOUDARY CONDITION------------------------------
!enforcing the fixed boundary conditions
num1 = SIZE(U_fix)
DO i = 1, num1
    nb = COUNT((Conn(U_fix(i),:) .LE. U_fix(i)) .AND. (Conn(U_fix(i),:) .NE. 0))
    DO j = 1, nb
        IF(Conn(U_fix(i),j) .EQ. U_fix(i))THEN
            K_global(K_pointer(U_fix(i),2) : K_pointer(U_fix(i)+1,2) - 1) = 0.D0 ! zero out the entire row
            K_global(K_pointer(U_fix(i),2)) = 1.D0
            K_global(K_pointer(U_fix(i),2) + ndim*K_pointer(U_fix(i),1)) = 1.D0
        ELSE
            m = COUNT(Conn(conn(U_fix(i),j),:) .NE. 0)
            n = COUNT(Conn(conn(U_fix(i),j),:) .GT. Conn(U_fix(i),j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(Conn(Conn(U_fix(i),j),k) .EQ. U_fix(i))THEN
                    K_global(K_pointer(Conn(U_fix(i),j),2) + ndim*num2 : K_pointer(Conn(U_fix(i),j),2) + ndim*num2 + 1) = 0.D0
                    K_global(K_pointer(Conn(U_fix(i),j),2) + ndim*K_pointer(Conn(U_fix(i),j),1) + ndim*num2 - 1 : K_pointer(Conn(U_fix(i),j),2) + ndim*K_pointer(Conn(U_fix(i),j),1) + ndim*num2) = 0.D0
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    P(ndim*U_fix(i) - 1) = 0.D0
    P(ndim*U_fix(i) ) = 0.D0
ENDDO

!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!!enforcing the applied displacement boundary conditions
num1 = SIZE(U_appy)
DO i = 1, num1
    nb = COUNT(conn(U_appy(i),:) .NE. 0)
    num3 = 0
    DO j = 1, nb
        IF(conn(U_appy(i),j) .EQ. U_appy(i))THEN
            !modify the force vector
            P(ndim*U_appy(i) - 1) = P(ndim*U_appy(i) - 1) - U_stepy*K_global(K_pointer(U_appy(i),2) + 1)
            P(ndim*U_appy(i)) = U_stepy
            !zero out the columns corresponding to the y component
            K_global(K_pointer(U_appy(i),2) + 1) = 0.D0
            K_global(K_pointer(U_appy(i),2) + ndim*K_pointer(U_appy(i),1)) = 1.D0
        ELSEIF(Conn(U_appy(i),j) .LT. U_appy(i))THEN
            m = COUNT(Conn(Conn(U_appy(i),j),:) .NE. 0)
            n = COUNT(Conn(Conn(U_appy(i),j),:) .GT. Conn(U_appy(i),j))
            num2 = 0
            DO k = m - n + 1, m
                num2 = num2 + 1
                IF(Conn(Conn(U_appy(i),j),k) .EQ. U_appy(i))THEN
                    !modify the force vector
                    P(ndim*Conn(U_appy(i),j) - 1) = P(ndim*Conn(U_appy(i),j) - 1) - U_stepy*K_global(K_pointer(Conn(U_appy(i),j),2) + ndim*num2 + 1)
                    P(ndim*Conn(U_appy(i),j)) = P(ndim*Conn(U_appy(i),j)) - U_stepy*K_global(K_pointer(Conn(U_appy(i),j),2) + ndim*K_pointer(conn(U_appy(i),j),1) + ndim*num2)
                    !zero out the columns corresponding to the y component
                    K_global(K_pointer(Conn(U_appy(i),j),2) + ndim*num2 + 1) = 0.D0
                    K_global(K_pointer(Conn(U_appy(i),j),2) + ndim*K_pointer(Conn(U_appy(i),j),1) + ndim*num2) = 0.D0
                ENDIF
            ENDDO
        ELSEIF(Conn(U_appy(i),j) .GT. U_appy(i))THEN
            !modify the force vector
            P(ndim*Conn(U_appy(i),j) - 1) = P(ndim*Conn(U_appy(i),j) - 1) - U_stepy*K_global(K_pointer(U_appy(i),2) + ndim*K_pointer(U_appy(i),1) + ndim*num3 + 1)
            P(ndim*Conn(U_appy(i),j)) = P(ndim*Conn(U_appy(i),j)) - U_stepy*K_global(K_pointer(U_appy(i),2) + ndim*K_pointer(U_appy(i),1) + ndim*num3 + 2)
            ! zero out the row components correspoding to the y component
            K_global(K_pointer(U_appy(i),2) + ndim*K_pointer(U_appy(i),1) + ndim*num3 + 1) = 0.D0
            K_global(K_pointer(U_appy(i),2) + ndim*K_pointer(U_appy(i),1) + ndim*num3 + 2) = 0.D0
            num3 = num3 + 1
        ENDIF
    ENDDO
ENDDO

END SUBROUTINE Get_K_P

!====================================================END Get_K_P===================================================

!====================================================FUNCTION du2dxy================================================

FUNCTION du2dxy(i, j)
IMPLICIT NONE

INTEGER:: i, j, k, nb
DOUBLE PRECISION:: du2dxy(4)
du2dxy = 0.D0

!first layer neighbors
nb = COUNT(neighbors1(i,:) .NE. 0)
DO k = 1, nb
    !d2udxidxi, 1 1
    du2dxy(1) = du2dxy(1) + (Kn1*csx(i,k,1) + Tv*SUM(csx(i,:,1)))*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(1. - csx(i,k,1)**2)/distance(i,k,1) + (Kn1 + Tv)*csx(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(1. - csx(i,k,1)**2)/distance(i,k,1)
    !d2udxidyi, 1 2
    du2dxy(2) = du2dxy(2) + (Kn1*csy(i,k,1) + Tv*SUM(csy(i,:,1)))*csx(i,k,1) - (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csy(i,k,1)*csx(i,k,1)/distance(i,k,1) + (Kn1 + Tv)*csy(i,k,1)*csx(i,k,1) - (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)* csy(i,k,1)/distance(i,k,1)
    !d2udyidyi, 2 2
    du2dxy(4) = du2dxy(4) + (Kn1*csy(i,k,1) + Tv*SUM(csy(i,:,1)))*csy(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(1. - csy(i,k,1)**2)/distance(i,k,1) + (Kn1 + Tv)*csy(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(1. - csy(i,k,1)**2)/distance(i,k,1)
ENDDO

! second layer neighbors
nb = COUNT(neighbors2(i,:) .NE. 0)
DO k = 1, nb
    !d2udxidxi, 1 1
    du2dxy(1) = du2dxy(1) + (Kn2*csx(i,k,2) + Tv*SUM(csx(i,:,2)))*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(1. - csx(i,k,2)**2)/distance(i,k,2) + (Kn2 + Tv)*csx(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(1. - csx(i,k,2)**2)/distance(i,k,2)
    !d2udxidyi, 1 2
    du2dxy(2) = du2dxy(2) + (Kn2*csy(i,k,2) + Tv*SUM(csy(i,:,2)))*csx(i,k,2) - (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csy(i,k,2)*csx(i,k,2)/distance(i,k,2) + (Kn2 + Tv)*csy(i,k,2)*csx(i,k,2) - (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)* csy(i,k,2)/distance(i,k,2)
    !d2udyidyi, 2 2
    du2dxy(4) = du2dxy(4) + (Kn2*csy(i,k,2) + Tv*SUM(csy(i,:,2)))*csy(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(1. - csy(i,k,2)**2)/distance(i,k,2) + (Kn2 + Tv)*csy(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(1. - csy(i,k,2)**2)/distance(i,k,2)
ENDDO

!d2udyidxi, 2 1
du2dxy(3) = du2dxy(2)

RETURN
END FUNCTION du2dxy

!====================================================END du2dxy=====================================================

!====================================================FUNCTION du2dxy1================================================

FUNCTION du2dxy1(i, j)
IMPLICIT NONE

INTEGER:: i, j, k, kk, m, n, nb1, nb2
DOUBLE PRECISION:: du2dxy1(4)
du2dxy1 = 0.D0

nb1 = COUNT(neighbors1(i,:) .NE. 0)
    ! the first nearest neighbor in AFEM
DO kk = 1, nb1
    IF(neighbors1(i,kk) .EQ. Conn(i,j)) THEN
        DO k = 1, nb1
            IF(k .EQ. kk)THEN
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxy1(1) = du2dxy1(1) - (Kn1 + Tv)*csx(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(-1. + csx(i,k,1)**2)/distance(i,k,1) - (Kn1*csx(i,k,1) - Tv*SUM(csx(neighbors1(i,k),:,1)))*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(-1. + csx(i,k,1)**2)/distance(i,k,1)
                !d2udxidyj, 1 2
                du2dxy1(2) = du2dxy1(2) - (Kn1 + Tv)*csy(i,k,1)*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1) - (Kn1*csy(i,k,1) - Tv*SUM(csy(neighbors1(i,k),:,1)))*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1)
                !d2udyidxj, 2 1
                du2dxy1(3) = du2dxy1(3) - (Kn1 + Tv)*csy(i,k,1)*csx(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1) - (Kn1*csx(i,k,1) - Tv*SUM(csx(neighbors1(i,k),:,1)))*csy(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*csx(i,k,1)*csy(i,k,1)/distance(i,k,1)
                !d2udyidyj, 2 2
                du2dxy1(4) = du2dxy1(4) - (Kn1 + Tv)*csy(i,k,1)**2 + (Kn1*dL(i,k,1) + Tv*SUM(dL(i,:,1)))*(-1. + csy(i,k,1)**2)/distance(i,k,1) - (Kn1*csy(i,k,1) - Tv*SUM(csy(neighbors1(i,k),:,1)))*csy(i,k,1) + (Kn1*dL(i,k,1) + Tv*SUM(dL(neighbors1(i,k),:,1)))*(-1. + csy(i,k,1)**2)/distance(i,k,1)
                !-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            ELSE
                !----------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxy1(1) = du2dxy1(1) - Tv*csx(i,kk,1)*csx(i,k,1)
                !d2udxidyj, 1 2
                du2dxy1(2) = du2dxy1(2) - Tv*csy(i,kk,1)*csx(i,k,1)
                !d2udyidxj, 2 1
                du2dxy1(3) = du2dxy1(3) - Tv*csx(i,kk,1)*csy(i,k,1)
                !d2udyidyj, 2 2
                du2dxy1(4) = du2dxy1(4) - Tv*csy(i,kk,1)*csy(i,k,1)
                !----------------------------------------------------------------------------------
                nb2 = COUNT(neighbors1(neighbors1(i,k),:) .NE. 0)
                DO n = 1,nb2
                    IF(neighbors1(neighbors1(i,k),n) .EQ. Conn(i,j))THEN
                        !---------------------------------------------------------------------------------------------------
                        !d2udxidxj, 1 1
                        du2dxy1(1) = du2dxy1(1) - Tv*csx(neighbors1(i,k),n,1)*csx(i,k,1)
                        !d2udxidyj, 1 2
                        du2dxy1(2) = du2dxy1(2) - Tv*csy(neighbors1(i,k),n,1)*csx(i,k,1)
                        !d2udyidxj, 2 1
                        du2dxy1(3) = du2dxy1(3) - Tv*csx(neighbors1(i,k),n,1)*csy(i,k,1)
                        !d2udyidyj, 2 2
                        du2dxy1(4) = du2dxy1(4) - Tv*csy(neighbors1(i,k),n,1)*csy(i,k,1)
                        !---------------------------------------------------------------------------------------------------
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        GOTO 777
    ENDIF
ENDDO
! the second nearest neighbor in AFEM
DO kk = 1, nb1
    nb2 = COUNT(neighbors1(neighbors1(i,kk),:) .NE. 0)
    DO m = 1, nb2
        IF(neighbors1(neighbors1(i,kk),m) .EQ. Conn(i,j))THEN
            !--------------------------------------------------------------------------------------------------------
            !d2udxidxj, 1 1
            du2dxy1(1) = du2dxy1(1) - Tv*csx(i,kk,1)*csx(neighbors1(i,kk),m,1)
            !d2udxidyj, 1 2
            du2dxy1(2) = du2dxy1(2) - Tv*csx(i,kk,1)*csy(neighbors1(i,kk),m,1)
            !d2udyidxj, 2 1
            du2dxy1(3) = du2dxy1(3) - Tv*csy(i,kk,1)*csx(neighbors1(i,kk),m,1)
            !d2udyidyj, 2 2
            du2dxy1(4) = du2dxy1(4) - Tv*csy(i,kk,1)*csy(neighbors1(i,kk),m,1)
            !---------------------------------------------------------------------------------------------------------
        ENDIF
    ENDDO
ENDDO

777 RETURN
END FUNCTION du2dxy1

!====================================================END du2dxy1=====================================================

!====================================================FUNCTION du2dxy2================================================

FUNCTION du2dxy2(i, j)
IMPLICIT NONE

INTEGER:: i, j, k, kk, m, n, nb1, nb2
DOUBLE PRECISION:: du2dxy2(4)
du2dxy2 = 0.D0

nb1 = COUNT(neighbors2(i,:) .NE. 0)
    ! the first nearest neighbor in AFEM
DO kk = 1, nb1
    IF(neighbors2(i,kk) .EQ. Conn(i,j)) THEN
        DO k = 1, nb1
            IF(k .EQ. kk)THEN
                !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxy2(1) = du2dxy2(1) - (Kn2 + Tv)*csx(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(-1. + csx(i,k,2)**2)/distance(i,k,2) - (Kn2*csx(i,k,2) - Tv*SUM(csx(neighbors2(i,k),:,2)))*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(-1. + csx(i,k,2)**2)/distance(i,k,2)
                !d2udxidyj, 1 2
                du2dxy2(2) = du2dxy2(2) - (Kn2 + Tv)*csy(i,k,2)*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2) - (Kn2*csy(i,k,2) - Tv*SUM(csy(neighbors2(i,k),:,2)))*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2)
                !d2udyidxj, 2 1
                du2dxy2(3) = du2dxy2(3) - (Kn2 + Tv)*csy(i,k,2)*csx(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2) - (Kn2*csx(i,k,2) - Tv*SUM(csx(neighbors2(i,k),:,2)))*csy(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*csx(i,k,2)*csy(i,k,2)/distance(i,k,2)
                !d2udyidyj, 2 2
                du2dxy2(4) = du2dxy2(4) - (Kn2 + Tv)*csy(i,k,2)**2 + (Kn2*dL(i,k,2) + Tv*SUM(dL(i,:,2)))*(-1. + csy(i,k,2)**2)/distance(i,k,2) - (Kn2*csy(i,k,2) - Tv*SUM(csy(neighbors2(i,k),:,2)))*csy(i,k,2) + (Kn2*dL(i,k,2) + Tv*SUM(dL(neighbors2(i,k),:,2)))*(-1. + csy(i,k,2)**2)/distance(i,k,2)
                !-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            ELSE
                !----------------------------------------------------------------------------------
                !d2udxidxj, 1 1
                du2dxy2(1) = du2dxy2(1) - Tv*csx(i,kk,2)*csx(i,k,2)
                !d2udxidyj, 1 2
                du2dxy2(2) = du2dxy2(2) - Tv*csy(i,kk,2)*csx(i,k,2)
                !d2udyidxj, 2 1
                du2dxy2(3) = du2dxy2(3) - Tv*csx(i,kk,2)*csy(i,k,2)
                !d2udyidyj, 2 2
                du2dxy2(4) = du2dxy2(4) - Tv*csy(i,kk,2)*csy(i,k,2)
                !----------------------------------------------------------------------------------
                nb2 = COUNT(neighbors2(neighbors2(i,k),:) .NE. 0)
                DO n = 1,nb2
                    IF(neighbors2(neighbors2(i,k),n) .EQ. Conn(i,j))THEN
                        !---------------------------------------------------------------------------------------------------
                        !d2udxidxj, 1 1
                        du2dxy2(1) = du2dxy2(1) - Tv*csx(neighbors2(i,k),n,2)*csx(i,k,2)
                        !d2udxidyj, 1 2
                        du2dxy2(2) = du2dxy2(2) - Tv*csy(neighbors2(i,k),n,2)*csx(i,k,2)
                        !d2udyidxj, 2 1
                        du2dxy2(3) = du2dxy2(3) - Tv*csx(neighbors2(i,k),n,2)*csy(i,k,2)
                        !d2udyidyj, 2 2
                        du2dxy2(4) = du2dxy2(4) - Tv*csy(neighbors2(i,k),n,2)*csy(i,k,2)
                        !---------------------------------------------------------------------------------------------------
                    ENDIF
                ENDDO
            ENDIF
        ENDDO
        GOTO 888
    ENDIF
ENDDO
! the second nearest neighbor in AFEM
DO kk = 1, nb1
    nb2 = COUNT(neighbors2(neighbors2(i,kk),:) .NE. 0)
    DO m = 1, nb2
        IF(neighbors2(neighbors2(i,kk),m) .EQ. Conn(i,j))THEN
            !--------------------------------------------------------------------------------------------------------
            !d2udxidxj, 1 1
            du2dxy2(1) = du2dxy2(1) - Tv*csx(i,kk,2)*csx(neighbors2(i,kk),m,2)
            !d2udxidyj, 1 2
            du2dxy2(2) = du2dxy2(2) - Tv*csx(i,kk,2)*csy(neighbors2(i,kk),m,2)
            !d2udyidxj, 2 1
            du2dxy2(3) = du2dxy2(3) - Tv*csy(i,kk,2)*csx(neighbors2(i,kk),m,2)
            !d2udyidyj, 2 2
            du2dxy2(4) = du2dxy2(4) - Tv*csy(i,kk,2)*csy(neighbors2(i,kk),m,2)
            !---------------------------------------------------------------------------------------------------------
        ENDIF
    ENDDO
ENDDO

888 RETURN
END FUNCTION du2dxy2

!====================================================END du2dxy2=====================================================

!====================================================FUNCTION dudxy=================================================

FUNCTION dudxy(i)
IMPLICIT NONE
 
INTEGER, INTENT(IN):: i
INTEGER:: j, nb
DOUBLE PRECISION:: dudxy(2)
dudxy = 0.D0

nb = COUNT(neighbors1(i,:) .NE. 0)
DO j = 1,nb
    !dudx
    dudxy(1) = dudxy(1) + (2.*Kn1*dL(i,j,1) + Tv*(SUM(dL(i,:,1)) + SUM(dL(neighbors1(i,j),:,1))))*csx(i,j,1)
    !dudy
    dudxy(2) = dudxy(2) + (2.*Kn1*dL(i,j,1) + Tv*(SUM(dL(i,:,1)) + SUM(dL(neighbors1(i,j),:,1))))*csy(i,j,1)
ENDDO

nb = COUNT(neighbors2(i,:) .NE. 0)
DO j = 1,nb
    !dudx
    dudxy(1) = dudxy(1) + (2.*Kn2*dL(i,j,2) + Tv*(SUM(dL(i,:,2)) + SUM(dL(neighbors2(i,j),:,2))))*csx(i,j,2)
    !dudy
    dudxy(2) = dudxy(2) + (2.*Kn2*dL(i,j,2) + Tv*(SUM(dL(i,:,2)) + SUM(dL(neighbors2(i,j),:,2))))*csy(i,j,2)
ENDDO

RETURN
END FUNCTION dudxy

!====================================================END dudxy======================================================

!==============================================SUBROUTINE PARDISOsolver=============================================

SUBROUTINE PARDISOsolver
IMPLICIT NONE

!.. Internal solver memory pointer 
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
!.. All other variables
INTEGER maxfct, mnum, mtype, phase, error, error1, msglvl
INTEGER, ALLOCATABLE :: iparm( : )
INTEGER i, idum(1)
DOUBLE PRECISION:: ddum(1), du(ndim*nparticle)
du = 0
!.. Fill all arrays containing matrix data.
maxfct = 1
mnum = 1
!.. Set up PARDISO control parameter
ALLOCATE( iparm ( 64 ) )
iparm = 0

iparm(1) = 1 ! no solver default
iparm(2) = 3 ! fill-in reordering from METIS
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
iparm(6) = 0 ! =0 solution on the first n compoments of x
iparm(8) = 0 ! numbers of iterative refinement steps
iparm(10) = 13 ! perturbe the pivot elements with 1E-13
iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(13) = 0 ! =0, maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy

error  = 0 ! initialize error flag
msglvl = 1 ! print statistical information
!msglvl = 0 ! no print statistical information
mtype  = 2 ! real and symmetric positive definite
!mtype  = 1 ! structurally symmetric
!mtype  = -2 ! symmetric, indefinite
!mtype  = 11! real unsymmetric
ALLOCATE ( pt ( 64 ) )

CALL Get_K_P
!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!.. Initiliaze the internal solver memory pointer. This is only necessary for the FIRST call of the PARDISO solver.
DO i = 1, 64
   pt( i )%DUMMY =  0 
ENDDO
!.. Reordering and Symbolic Factorization, This step also allocates all memory that is necessary for the factorization
phase = 11 ! only reordering and symbolic factorization
CALL PARDISO (pt, maxfct, mnum, mtype, phase, ndim*nparticle, K_global, IK, JK, idum, 1, iparm, msglvl, ddum, ddum, error)
Write(*,*) 'PARDISO: Size of factors(MB): ', MAX(iparm(15), iparm(16) + iparm(17))/1000

!.. Factorization.
phase = 22 ! only factorization
CALL PARDISO (pt, maxfct, mnum, mtype, phase, ndim*nparticle, K_global, IK, JK, idum, 1, iparm, msglvl, ddum, ddum, error)

!.. Back substitution and iterative refinement
iparm(8) = 9 ! max numbers of iterative refinement steps
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
Positions(:,1) = Positions(:,1) + du(1:ndim*nparticle:2)
Positions(:,2) = Positions(:,2) + du(2:ndim*nparticle:2)


END SUBROUTINE PARDISOsolver

!==============================================END PARDISOsolver====================================================

END PROGRAM MAIN