/******************************************************************************
 * Implemented the solving process of nonlocal lattice particle method(LPM) in
 * general cases(2 layers of neighbors)
 * This code can deal with general elastic and plasticity problems under careful
 * modifications
 * 
 * Supported lattice types: 
 *     2D square, 2D hexagon, 3D simple cubic, 3D BCC & FCC crystal(crystal 
 *     plasticity simulation)
 * 
 * Contact information:
 *     Changyu Meng, PhD candidate, Arizona State University, cmeng12@asu.edu
******************************************************************************/

#include "lpm.h"

/* definition of global variables */
/* int */
int nparticle;
int nbreak;

int *IK, *JK, *type, *dispBC_index, *fix_index, *lacknblist, *pl_flag;
int *nb, *nb_initial, *nb_conn;
int **neighbors, **neighbors1, **neighbors2, **neighbors_AFEM;
int **K_pointer, **conn, **nsign;

/* double precision float */
double critical_bstrain;

double *K_global, *residual, *Pin, *Pex, *Pex_temp, *disp;
double *reaction_force, *damage_visual;

double **xyz, **xyz_initial, **xyz_temp, **distance, **distance_initial, **KnTve, **F, **csx, **csy, **csz;
double **dL_total, **TdL_total, **csx_initial, **csy_initial, **csz_initial, **Ce;
double **stress_tensor;
double **strain_tensor, **ddL_total, **TddL_total, **F_temp;
double **dL, **ddL, **bond_stress, **damage_broken, **damage_w, **bond_vector;

double **Kn, **Tv;
double ***slip_vector, ***damage_D;


/************************************************************************/
/****************************** Main procedure **************************/
/************************************************************************/
int main(int argc, char *argv[])
{
    printf("\n==================================================\n");
    printf("            Nonlocal LPM Program in C             \n");
    printf("==================================================\n");

    const int nt = omp_get_max_threads(); /* maximum number of threads provided by the computer */
    const int nt_force = 2;               /* number of threads for general force calculation */

    omp_set_num_threads(nt);
    printf("OpenMP with %d/%d threads for bond force calculation\n", nt_force, nt);

    double start = omp_get_wtime(); // record the CPU time, begin

    // lattice type -> lattice number, int
    // square -> 0; hexagon -> 1; simple cubic -> 2; face-centered cubic -> 3
    // body-centered cubic with 1 types of slip systems -> 4
    // body-centered cubic with 2 types of slip systems -> 5
    // body-centered cubic with 3 types of slip systems -> 6
    int lattice = 2;
    double radius = 0.25;
    UnitCell cell(lattice, radius); /*lattice type is 2, simple cubic*/

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 45.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    double box[] = {-0.2, 10.2, -0.2, 10.2, -0.2, 10.2};
    
    createCuboid(box, cell, R_matrix);

    // move the particles coordinates
    double offset[] = {-0., -0., -0.};
    moveParticle(box, offset);

    // initialize the necessary matrices
    initMatrices(cell);
    copyDouble2D(xyz_initial, xyz, nparticle, NDIM);

    printf("\nParticle number is %d\n", nparticle);

    // search neighbor
    searchNormalNeighbor(cell);
    int len_K_global = searchAFEMNeighbor(cell);


    // assign types for particles located in different rectangular regions
    // xlo, xhi, ylo, yhi, zlo, zhi, type
    int ntype = 0;
    type = allocInt1D(nparticle, ntype++); // set particle type as 0 in default

    setTypeRect(-100.0, 100.0, -100.0, 100.0, 10 - 1.2 * radius, 100.0, ntype++); // top layer, type 1,
    setTypeRect(-100.0, 100.0, -100.0, 100.0, -100.0, 1.2 * radius, ntype++);     // lower layer, type 2

    setTypeFullNeighbor(ntype++, cell); // set the particle type as 3 if the neighbor list is full

    // material elastic parameters setting, MPa
    double C11, C12, C44;
    double E0 = 146e3, mu0 = 0.3;
    // plane strain or 3D
    C11 = E0 * (1.0 - mu0) / (1.0 + mu0) / (1.0 - 2.0 * mu0);
    C12 = E0 * mu0 / (1.0 + mu0) / (1.0 - 2.0 * mu0);
    C44 = E0 / 2.0 / (1.0 + mu0);
    // plane stress
    // C11 = E0 / (1.0 - mu0) / (1.0 + mu0);
    // C12 = E0 * mu0 / (1.0 - mu0) / (1.0 + mu0);
    // C44 = E0 / 2.0 / (1.0 + mu0);
    // 3D cases
    // C11 = 108.2e3; // in Hailong's paper
    // C12 = 61.3e3;
    // C44 = 28.5e3;
    Ce = allocDouble2D(ntype, 3, 0.);
    for (int i = 0; i < ntype; i++)
    {
        Ce[i][0] = C11; // C11
        Ce[i][1] = C12; // C12
        Ce[i][2] = C44; // C44
    }

    // elasticity or plasticity model settings
    // (6) elastic (brittle) material
    int plmode = 6;

    printf("Constitutive mode is %d\n", plmode);

    // damage settings
    nbreak = 20;               // limit the broken number of bonds in a single iteration, should be an even number
    critical_bstrain = 1.0e-2; // critical bond strain value at which bond will break

    // define the crack
    // defineCrack(ca1, 5. + w, ch);

    // output parameters
    char aflag = 'z';      // output directional axis
    int dtype = 2;         // particle type used to output displacement or force
    int outdisp_flag = 1;  // output displacement, 0 or 1
    int outforce_flag = 1; // output force, 0 or 1
    int out_step = 1;      // print initial strain and dump info, output step is set as 1
    char dumpflag = 'm';   // 's' for single step; 'm' for multiple step

    char dispFile[] = "result_disp.txt";
    char dispFile1[] = "result_disp_CMOD1.txt";
    char dispFile2[] = "result_disp_CMOD2.txt";
    char forceFile[] = "result_force.txt";
    char strainFile[] = "result_strain.txt";
    char stressFile[] = "result_stress.txt";
    char dumpFile[] = "result_position.dump";
    char actFile[] = "result_Jact.txt";
    char slipFile[] = "result_slipRSS.txt";
    char damageFile[] = "result_damage.txt";
    char bondFile[] = "result_brokenbonds.txt";
    char stretchFile[] = "result_stretch.txt";
    char dlambdaFile[] = "result_dlambda.txt";
    char neighborFile[] = "result_neighbor.txt";
    char cabFile[] = "result_Cab.txt";
    char bforceFile[] = "result_bforce.txt";

    // boundary conditions and whole simulation settings
    int n_steps = 2;           // number of loading steps
    double step_size = -2000.0; // step size for force or displacement loading
    // int n_steps = 10;        // number of loading steps
    // double step_size = -2e-3; // step size for force or displacement loading

    int nbd = 0, nbf = 0;     // total number of disp or force boundary conditions
    char cal_method[] = "cg"; // calculation method, pardiso or conjugate gradient
    struct DispBCs dBP[MAXLINE][MAXSMALL] = {0};
    struct ForceBCs fBP[MAXLINE][MAXSMALL] = {0};
    int load_indicator[MAXLINE] = {0}; // tension (1) or compressive (-1) loading condition for uniaxial loading

    // displace boundary conditions
    // for (int i = 0; i < n_steps; i++)
    // {
    //     nbd = 0;
    //     load_indicator[i] = 1;
    //     dBP[i][nbd].type = 1, dBP[i][nbd].flag = 'z', dBP[i][nbd++].step = 0.0;
    //     dBP[i][nbd].type = 2, dBP[i][nbd].flag = 'z', dBP[i][nbd++].step = step_size;
    // }

    // force boundary conditions
    for (int i = 0; i < n_steps; i++)
    {
        load_indicator[i] = 1;

        nbd = 0;
        dBP[i][nbd].type = 1, dBP[i][nbd].flag = 'z', dBP[i][nbd++].step = 0.0;

        nbf = 0;
        fBP[i][nbf].type = 2;
        fBP[i][nbf].flag1 = 'x', fBP[i][nbf].step1 = 0.0;
        fBP[i][nbf].flag2 = 'y', fBP[i][nbf].step2 = 0.0;
        fBP[i][nbf].flag3 = 'z', fBP[i][nbf++].step3 = step_size;
    }

    /************************** Simulation begins *************************/
    // compute necessary into before loading starts
    calcKnTv(ntype, cell);
    computedL();

    // output files
    if (outdisp_flag)
    {
        writeStrain(strainFile, 3, 0);
        writeDisp(dispFile, aflag, dtype, 0);
        // writeDisp(dispFile1, 'y', 5, 0);
        // writeDisp(dispFile2, 'y', 6, 0);
    }
    if (outforce_flag)
    {
        writeStress(stressFile, 3, 0);
        writeReaction(forceFile, aflag, dtype, 0);
        // writeForce(forceFile, aflag, 0.3, 0); // lower half body force
    }
    writeDump(dumpFile, 0, dumpflag, box, plmode);

    // writeCab(cabFile, 0);
    // writeCab(cabFile, 288);
    // writeCab(cabFile, 4227);
    // writeCab(cabFile, 4228);
    // writeCab(cabFile, 4229);

    // compute the elastic stiffness matrix
    // omp_set_num_threads(nt_force);
    // if (dim == 2)
    //     calcStiffness2DFiniteDifference(6);
    // else if (dim == 3)
    //     calcStiffness3DFiniteDifference(6);
    // omp_set_num_threads(nt);

    // omp_set_num_threads(nt_force);
    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    // incremental loading procedure
    for (int i = 0; i < n_steps; i++)
    {
        double startrun = omp_get_wtime();

        printf("######################################## Loading step %d ######################################\n", i + 1);
        copyDouble2D(xyz_temp, xyz, nparticle, NDIM);
        copyDouble2D(F_temp, F, nparticle, cell.nneighbors);
        copyDouble1D(Pex_temp, Pex, cell.dim * nparticle);

        omp_set_num_threads(nt_force);
        // compute the elastic stiffness matrix
        if (cell.dim == 2)
            calcStiffness2DFiniteDifference(6, cell);
        else if (cell.dim == 3)
            calcStiffness3DFiniteDifference(6, cell);

        double time_t1 = omp_get_wtime();
        printf("Stiffness matrix calculation costs %f seconds\n", time_t1 - startrun);

    label_wei:
        setDispBC(nbd, dBP[i], cell);  // update displacement BC
        setForceBC(nbf, fBP[i], cell); // update force BC

        computeBondForceGeneral(plmode, load_indicator[i], cell); // incremental updating
        omp_set_num_threads(nt);

    label_broken_bond:
        updateRR(cell); // residual force vector and reaction force (RR)

        // compute the Euclidean norm (L2 norm)
        double norm_residual = cblas_dnrm2(cell.dim * nparticle, residual, 1);
        double norm_reaction_force = cblas_dnrm2(countNEqual(dispBC_index, nparticle * cell.dim, 1), reaction_force, 1);
        double tol_multiplier = MAX(norm_residual, norm_reaction_force);
        char tempChar1[] = "residual", tempChar2[] = "reaction";
        printf("\nNorm of residual is %.5e, norm of reaction is %.5e, tolerance criterion is based on ", norm_residual, norm_reaction_force);
        if (norm_residual > norm_reaction_force)
            printf("%s force\n", tempChar1);
        else
            printf("%s force\n", tempChar2);

        // global Newton iteration starts
        int ni = 0;
        while (norm_residual > TOLITER * tol_multiplier && ni < MAXITER)
        {
            printf("Step-%d, iteration-%d: ", i + 1, ++ni);

            // time_t1 = omp_get_wtime();
            // compute the stiffness matrix, then modify it for displacement boundary condition
            if (cell.dim == 2)
            {
                // calcStiffness2DFiniteDifference(6);
                setDispBC_stiffnessUpdate2D(cell);
            }
            else if (cell.dim == 3)
            {
                // calcStiffness3DFiniteDifference(6);
                setDispBC_stiffnessUpdate3D(cell);
            }
            // time_t2 = omp_get_wtime();
            // printf("Modify stiffness costs %f seconds\n", time_t2 - time_t1);

            // solve for the incremental displacement
            if (strcmp(cal_method, "pardiso") == 0)
                solverPARDISO(cell);
            else if (strcmp(cal_method, "cg") == 0)
                solverCG(cell);
            // time_t1 = omp_get_wtime();
            // printf("Solve the linear system costs %f seconds\n", time_t1 - time_t2);

            omp_set_num_threads(nt_force);
            computeBondForceGeneral(4, load_indicator[i], cell); // update the bond force
            omp_set_num_threads(nt);

            // time_t2 = omp_get_wtime();
            // printf("Update bond force costs %f seconds\n", time_t2 - time_t1);

            // writeDlambda(dlambdaFile, 9039, 9047, i + 1, ni);

            updateRR(cell); /* update the RHS risidual force vector */
            norm_residual = cblas_dnrm2(cell.dim * nparticle, residual, 1);
            printf("Norm of residual is %.3e, residual ratio is %.3e\n", norm_residual, norm_residual / tol_multiplier);
        }
        computeStrain(cell);

        /* accumulate damage, and break bonds when damage reaches critical values */
        int broken_bond = updateDamageGeneral(bondFile, i + 1, plmode, cell);
        updateCrack(cell);

        printf("Loading step %d has finished in %d iterations\n\nData output ...\n", i + 1, ni);

        /* ----------------------- data output setting section ------------------------ */

        if ((i + 1) % out_step == 0)
        {
            if (outdisp_flag)
            {
                writeStrain(strainFile, 3, i + 1);
                writeDisp(dispFile, aflag, dtype, i + 1);
                // writeDisp(dispFile1, 'y', 5, i + 1);
                // writeDisp(dispFile2, 'y', 6, i + 1);
            }
            if (outforce_flag)
            {
                writeStress(stressFile, 3, i + 1);
                writeReaction(forceFile, aflag, dtype, i + 1);
                // writeForce(forceFile, aflag, 0.3, i + 1);
            }
            writeDump(dumpFile, i + 1, dumpflag, box, plmode); /* print particle position info */
        }

        // check if breakage happens, recompute stiffness matrix (this is same as a new loading step)
        if (broken_bond > 0)
        {
            // switchStateV(0); // copy the last converged state variable [1] into current one [0]
            // computeBondForceGeneral(plmode, load_indicator[i]); // update the bond force

            startrun = omp_get_wtime();
            omp_set_num_threads(nt_force);
            // recompute stiffness matrix
            if (cell.dim == 2)
                calcStiffness2DFiniteDifference(6, cell);
            else if (cell.dim == 3)
                calcStiffness3DFiniteDifference(6, cell);
            omp_set_num_threads(nt);
            time_t1 = omp_get_wtime();
            printf("\nStiffness matrix calculation costs %f seconds\n", time_t1 - startrun);

            goto label_broken_bond;
        }

        double finishrun = omp_get_wtime();
        printf("Time costed for step %d: %f seconds\n\n", i + 1, finishrun - startrun);
    }

    writeNeighbor(neighborFile);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);

    // frees unused memory allocated by the Intel MKL Memory Allocator
    mkl_free_buffers();

    return 0;
}

/************************************************************************/
/*************************** End main procedures ************************/
/************************************************************************/

/* compute the bond force and the state variables using constitutive relationship determined by plmode */
void computeBondForceGeneral(int plmode, int t, UnitCell cell)
{
    // double el_radius = 7.0;
    // double pc1[3] = {10.0, 13.0, 0.0}, pc2[3] = {10.0, 35.0, 0.0};

    if (plmode == 0)
    {
#pragma omp parallel for
        for (int iID = 0; iID < nparticle; iID++)
        {
            computeBondForceElastic(iID, cell);
        }
    }
    else if (plmode == 4)
    {
#pragma omp parallel for
        for (int iID = 0; iID < nparticle; iID++)
            computeBondForceIncrementalUpdating(iID, cell);
    }
    else if (plmode == 6)
    {
#pragma omp parallel for
        for (int iID = 0; iID < nparticle; iID++)
            computeBondForceElastic(iID, cell);
    }

    computeStress(cell);
}

/* update the damage variables */
int updateDamageGeneral(const char *dataName, int tstep, int plmode, UnitCell cell)
{
    int broken = 0;

    if (plmode == 6)
        broken = updateBrittleDamage(dataName, tstep, nbreak);

    return broken;
}

/* allocate memories for some global matrices */
void initMatrices(UnitCell cell)
{
    disp = allocDouble1D(cell.dim * nparticle, 0);        /* global displacement vector */
    xyz_initial = allocDouble2D(nparticle, NDIM, 0); /* store coordinate information as 3D */
    xyz_temp = allocDouble2D(nparticle, NDIM, 0);

    distance = allocDouble2D(nparticle, cell.nneighbors, 0.0);
    distance_initial = allocDouble2D(nparticle, cell.nneighbors, 0.0);
    csx = allocDouble2D(nparticle, cell.nneighbors, 0.0);
    csy = allocDouble2D(nparticle, cell.nneighbors, 0.0);
    csz = allocDouble2D(nparticle, cell.nneighbors, 0.0);
    csx_initial = allocDouble2D(nparticle, cell.nneighbors, 0.0);
    csy_initial = allocDouble2D(nparticle, cell.nneighbors, 0.0);
    csz_initial = allocDouble2D(nparticle, cell.nneighbors, 0.0);

    // bond geometric measure
    dL = allocDouble2D(nparticle, cell.nneighbors, 0.0); /* total bond stretch */
    dL_total = allocDouble2D(nparticle, 2, 0);
    TdL_total = allocDouble2D(nparticle, 2, 0);
    ddL = allocDouble2D(nparticle, cell.nneighbors, 0.0);
    ddL_total = allocDouble2D(nparticle, 2, 0);
    TddL_total = allocDouble2D(nparticle, 2, 0);
    bond_stress = allocDouble2D(nparticle, cell.nneighbors, 0); /* bond stress, projection of stress tensor */

    pl_flag = allocInt1D(nparticle, 0); /* denote whether the plastic deformtion has been calculated, to avoid repetition */

    nb = allocInt1D(nparticle, -1);         /* number of normal bonds */
    nb_initial = allocInt1D(nparticle, -1); /* initial number of normal bonds */
    nb_conn = allocInt1D(nparticle, -1);    /* number of connections, including itself */
    neighbors = allocInt2D(nparticle, cell.nneighbors, -1);
    neighbors1 = allocInt2D(nparticle, cell.nneighbors1, -1);
    neighbors2 = allocInt2D(nparticle, cell.nneighbors2, -1);
    nsign = allocInt2D(nparticle, cell.nneighbors, -1);
    conn = allocInt2D(nparticle, cell.nneighbors_AFEM + 1, -1); /* connection of AFEM particles */

    K_pointer = allocInt2D(nparticle + 1, 2, 0);

    residual = allocDouble1D(cell.dim * nparticle, 0);  /* right hand side, residual */
    Pin = allocDouble1D(NDIM * nparticle, 0);      /* total internal force, fixed 3 dimension */
    Pex = allocDouble1D(cell.dim * nparticle, 0);       /* total external force */
    Pex_temp = allocDouble1D(cell.dim * nparticle, 0);  /* temp external force */
    dispBC_index = allocInt1D(cell.dim * nparticle, 1); /* disp info for each degree of freedom, 0 as being applied disp BC */
    fix_index = allocInt1D(cell.dim * nparticle, 1);    /* fix info for each degree of freedom, 0 as being fixed */

    Kn = allocDouble2D(nparticle, cell.nneighbors, 0.0);
    Tv = allocDouble2D(nparticle, cell.nneighbors, 0.0);

    stress_tensor = allocDouble2D(nparticle, 2 * NDIM, 0.);
    strain_tensor = allocDouble2D(nparticle, 2 * NDIM, 0.);

    // damage parameters
    damage_broken = allocDouble2D(nparticle, cell.nneighbors, 1.); /* crack index parameter */
    damage_D = allocDouble3D(nparticle, cell.nneighbors, 2, 0.);   /* bond-associated damage parameter, initially be 0, state variable */
    damage_w = allocDouble2D(nparticle, cell.nneighbors, 1.);      /* intact parameter, apply on bond force */
    damage_visual = allocDouble1D(nparticle, 0.0);            /* damage indicator for visualization */

    F = allocDouble2D(nparticle, cell.nneighbors, 0.);          /* FIJ is total bond force between I and J */
    F_temp = allocDouble2D(nparticle, cell.nneighbors, 0.);     /* F_tempIJ stores the bond force at last time step */
}