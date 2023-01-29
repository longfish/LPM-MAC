#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <mkl.h>
#include <omp.h>

#include "lpm.h"
#include "utilities.h"

struct BondStrain
{
    int i_index;
    int j_index;
    double bstrain;
};

/* Incremental updating (elastic) in general and for Wei's algorithm, trial */
void computeBondForceIncrementalUpdating(int ii, UnitCell cell)
{
    // create a temporary list to store particle ii and all of its neighbors ID
    int *temporary_nb = allocInt1D(nb[ii] + 1, ii);
    int s = 1;
    for (int k = 1; k < cell.nneighbors + 1; k++)
    {
        if (damage_broken[ii][k - 1] > EPS && neighbors[ii][k - 1] != -1)
            temporary_nb[s++] = neighbors[ii][k - 1];
    }

    //#pragma omp parallel for
    // initialize the dilatational change for the above particle list
    for (int k = 0; k < nb[ii] + 1; k++)
    {
        int i = temporary_nb[k];
        ddL_total[i][0] = 0, TddL_total[i][0] = 0;
        ddL_total[i][1] = 0, TddL_total[i][1] = 0;

        for (int j = 0; j < nb_initial[i]; j++)
        {
            double dis0 = sqrt(pow(xyz_temp[i][0] - xyz_temp[neighbors[i][j]][0], 2) +
                               pow(xyz_temp[i][1] - xyz_temp[neighbors[i][j]][1], 2) +
                               pow(xyz_temp[i][2] - xyz_temp[neighbors[i][j]][2], 2));
            double dis1 = sqrt(pow(xyz[i][0] - xyz[neighbors[i][j]][0], 2) +
                               pow(xyz[i][1] - xyz[neighbors[i][j]][1], 2) +
                               pow(xyz[i][2] - xyz[neighbors[i][j]][2], 2));
            ddL[i][j] = damage_broken[i][j] * (dis1 - dis0);
            ddL_total[i][nsign[i][j]] += ddL[i][j];
            TddL_total[i][nsign[i][j]] += Tv[i][j] * ddL[i][j];
            // csx[i][j] = (xyz[i][0] - xyz[neighbors[i][j]][0]) / dis1;
            // csy[i][j] = (xyz[i][1] - xyz[neighbors[i][j]][1]) / dis1;
            // csz[i][j] = (xyz[i][2] - xyz[neighbors[i][j]][2]) / dis1;
        }
    }

    // update the bond force for the above particle list
    int i = ii;
    // memset(stress_tensor[i], 0.0, 2 * NDIM * sizeof(double));              // initialize the stress tensor
    Pin[i * NDIM] = 0.0, Pin[i * NDIM + 1] = 0.0, Pin[i * NDIM + 2] = 0.0; // initialize the internal force vector

    //#pragma omp parallel for
    for (int j = 0; j < nb_initial[i]; j++)
    {
        /* update Fij */
        F[i][j] = F_temp[i][j] + 2.0 * Kn[i][j] * ddL[i][j] + 0.5 * (TddL_total[i][nsign[i][j]] + TddL_total[neighbors[i][j]][nsign[i][j]]) + 0.5 * Tv[i][j] * (ddL_total[i][nsign[i][j]] + ddL_total[neighbors[i][j]][nsign[i][j]]);
        F[i][j] *= damage_broken[i][j];
        // F[i][j] *= damage_w[i][j];
        // F[i][j] *= (1.0 - damage_D[i][j][0]);

        /* compute internal forces */
        Pin[i * NDIM] += csx[i][j] * F[i][j];
        Pin[i * NDIM + 1] += csy[i][j] * F[i][j];
        Pin[i * NDIM + 2] += csz[i][j] * F[i][j]; // 0 for 2D
    }

    // free the array memory of particle ii
    free(temporary_nb);
}

/* Elastic material law */
void computeBondForceElastic(int ii, UnitCell cell)
{
    // create a temporary list to store particle ii and all of its neighbors ID
    int *temporary_nb = allocInt1D(nb[ii] + 1, ii);
    int s = 1;
    for (int k = 1; k < cell.nneighbors + 1; k++)
    {
        if (damage_broken[ii][k - 1] > EPS && neighbors[ii][k - 1] != -1)
            temporary_nb[s++] = neighbors[ii][k - 1];
    }

    //#pragma omp parallel for
    // compute the dilatational change for the above particle list
    for (int k = 0; k < nb[ii] + 1; k++)
    {
        int i = temporary_nb[k];
        dL_total[i][0] = 0, TdL_total[i][0] = 0;
        dL_total[i][1] = 0, TdL_total[i][1] = 0;

        for (int j = 0; j < nb_initial[i]; j++)
        {
            double dis = sqrt(pow(xyz[i][0] - xyz[neighbors[i][j]][0], 2) +
                              pow(xyz[i][1] - xyz[neighbors[i][j]][1], 2) +
                              pow(xyz[i][2] - xyz[neighbors[i][j]][2], 2));
            dL[i][j] = dis - distance_initial[i][j];
            dL[i][j] *= damage_broken[i][j];
            dL_total[i][nsign[i][j]] += dL[i][j];
            TdL_total[i][nsign[i][j]] += Tv[i][j] * dL[i][j];
            csx[i][j] = (xyz[i][0] - xyz[neighbors[i][j]][0]) / dis;
            csy[i][j] = (xyz[i][1] - xyz[neighbors[i][j]][1]) / dis;
            csz[i][j] = (xyz[i][2] - xyz[neighbors[i][j]][2]) / dis;
        }
    }

    // update the bond force
    int i = ii;
    // memset(stress_tensor[i], 0.0, 2 * NDIM * sizeof(double));              // initialize the stress tensor
    Pin[i * NDIM] = 0.0, Pin[i * NDIM + 1] = 0.0, Pin[i * NDIM + 2] = 0.0; // initialize the internal force vector
    for (int j = 0; j < nb_initial[i]; j++)
    {
        /* update Fij */
        F[i][j] = 2.0 * Kn[i][j] * dL[i][j] + 0.5 * (TdL_total[i][nsign[i][j]] + TdL_total[neighbors[i][j]][nsign[i][j]]) + 0.5 * Tv[i][j] * (dL_total[i][nsign[i][j]] + dL_total[neighbors[i][j]][nsign[i][j]]);
        F[i][j] *= damage_broken[i][j];
        // F[i][j] *= damage_w[i][j];
        // F[i][j] *= (1.0 - damage_D[i][j][0]);

        /* compute internal forces */
        Pin[i * NDIM] += csx[i][j] * F[i][j];
        Pin[i * NDIM + 1] += csy[i][j] * F[i][j];
        Pin[i * NDIM + 2] += csz[i][j] * F[i][j]; // 0 for 2D
    }

    // finish the loop for particle ii
    free(temporary_nb);
}

/* update the damage state after bonds broken */
void updateCrack(UnitCell cell)
{
    // #pragma omp parallel for
    for (int i = 0; i < nparticle; i++)
    {
        nb[i] = nb_initial[i];
        damage_visual[i] = 0.0;
        Pin[i * NDIM] = 0.0, Pin[i * NDIM + 1] = 0.0, Pin[i * NDIM + 2] = 0.0; // initialize the internal force vector
        for (int j = 0; j < nb_initial[i]; j++)
        {
            if (damage_broken[i][j] <= EPS)
                nb[i] -= 1;

            damage_visual[i] += damage_broken[i][j];

            /* update the internal forces */
            F[i][j] *= damage_w[i][j];
            // F[i][j] *= (1.0 - damage_D[i][j][0]);
            // F[i][j] *= damage_broken[i][j];
            Pin[i * NDIM] += csx[i][j] * F[i][j];
            Pin[i * NDIM + 1] += csy[i][j] * F[i][j];
            Pin[i * NDIM + 2] += csz[i][j] * F[i][j]; // 0 for 2D
        }

        if (nb[i] < 1)
        {
            // fix this particle's position, if all bonds are broken
            fix_index[cell.dim * i] = 0;
            fix_index[cell.dim * i + 1] = 0;
            if (cell.dim == 3)
                fix_index[cell.dim * i + 2] = 0;
        }

        damage_visual[i] = 1 - damage_visual[i] / nb_initial[i];
    }
}

/* break the bond for elastic materials, when bond strain reaches a critical value (we limit the maximum broken number) */
int updateBrittleDamage(const char *dataName, int tstep, int nbreak)
{
    FILE *fpt;
    fpt = fopen(dataName, "a+");
    fprintf(fpt, "TIMESTEP ");
    fprintf(fpt, "%d\n", tstep);

    struct BondStrain b_cr[MAXSMALL * MAXSMALL]; // bonds of which the strain reach the critical value

    // initialization of struct
    for (int i = 0; i < MAXSMALL * MAXSMALL; i++)
    {
        b_cr[i].i_index = -1;
        b_cr[i].j_index = -1;
        b_cr[i].bstrain = 0.0;
    }

    // store bonds that have reached the critical bond strain, if no new crack, k = 0
    int k = 0;
    for (int i = 0; i < nparticle; i++)
    {
        for (int j = 0; j < nb_initial[i]; j++)
        {
            double ave_bstrain = dL[i][j] / distance_initial[i][j]; // note we already make dL to be 0 for broken bonds

            if (ave_bstrain >= critical_bstrain)
            {
                b_cr[k].i_index = i;
                b_cr[k].j_index = j;
                b_cr[k].bstrain = ave_bstrain;
                k++;
            }
        }
    }

    // delete the first nbreak bonds which have largest damage
    int broken_bonds = k;
    if (k > 0 && k <= nbreak)
    {
        for (int i = 0; i < k; i++)
        {
            damage_D[b_cr[i].i_index][b_cr[i].j_index][0] = 1.0;
            damage_w[b_cr[i].i_index][b_cr[i].j_index] = 0.0;
            damage_broken[b_cr[i].i_index][b_cr[i].j_index] = 0.0;
            fprintf(fpt, "%d %d \n", b_cr[i].i_index, neighbors[b_cr[i].i_index][b_cr[i].j_index]);
        }
    }
    else if (k > nbreak)
    {
        // sorting of the struct array, shell sort, from small to large
        int temp_i, temp_j;
        double temp_b;
        for (int r = k / 2; r >= 1; r = r / 2)
        {
            for (int i = r; i < k; ++i)
            {
                temp_i = b_cr[i].i_index;
                temp_j = b_cr[i].j_index;
                temp_b = b_cr[i].bstrain;
                int j = i - r;
                while (j >= 0 && b_cr[j].bstrain > temp_b)
                {
                    b_cr[j + r].bstrain = b_cr[j].bstrain;
                    b_cr[j + r].i_index = b_cr[j].i_index;
                    b_cr[j + r].j_index = b_cr[j].j_index;
                    j = j - r;
                }
                b_cr[j + r].bstrain = temp_b;
                b_cr[j + r].i_index = temp_i;
                b_cr[j + r].j_index = temp_j;
            }
        }

        // // for (int i = 0; i < k; i++)
        // //     printf("%d %d %f\n", b_cr[i].i_index, b_cr[i].j_index, b_cr[i].bstrain);

        // limit the maximum number of broken bonds
        for (int i = k - nbreak; i < k; i++)
        {
            damage_D[b_cr[i].i_index][b_cr[i].j_index][0] = 1.0;
            damage_w[b_cr[i].i_index][b_cr[i].j_index] = 0.0;
            damage_broken[b_cr[i].i_index][b_cr[i].j_index] = 0.0;
            fprintf(fpt, "%d %d \n", b_cr[i].i_index, neighbors[b_cr[i].i_index][b_cr[i].j_index]);
        }
    }

    fclose(fpt);

    return broken_bonds;
}
