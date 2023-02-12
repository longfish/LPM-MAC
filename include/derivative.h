#pragma once
#ifndef DERIVATIVE_H
#define DERIVATIVE_H

/* note: the below code ignores some small perturbation terms (so some values may be different from finite difference) */

template <int nlayer>
std::array<std::array<double, NDIM>, NDIM> fdu2dxyz(Particle<nlayer> *pi) /* The diagonal part */
{
    std::array<std::array<double, NDIM>, NDIM> du2dxyz{0};

    for (int i = 0; i < nlayer; i++)
    {
        for (Bond<nlayer> *bd : pi->bond_layers[i])
        {
            /*d2udxidxi, 1 1*/
            du2dxyz[0][0] += (bd->Kn * bd->csx + 0.5 * pi->Tcs_sumx[i] + 0.5 * bd->Tv * pi->cs_sumx[i]) * bd->csx +
                             (bd->Kn + bd->Tv) * pow(bd->csx, 2);
            /*d2udxidyi, 1 2*/
            du2dxyz[0][1] += (bd->Kn * bd->csy + 0.5 * pi->Tcs_sumy[i] + 0.5 * bd->Tv * pi->cs_sumy[i]) * bd->csx +
                             (bd->Kn + bd->Tv) * bd->csy * bd->csx;
            /*d2udxidzi, 1 3*/
            du2dxyz[0][2] += (bd->Kn * bd->csz + 0.5 * pi->Tcs_sumz[i] + 0.5 * bd->Tv * pi->cs_sumz[i]) * bd->csx +
                             (bd->Kn + bd->Tv) * bd->csz * bd->csx;
            /*d2udyidyi, 2 2*/
            du2dxyz[1][1] += (bd->Kn * bd->csy + 0.5 * pi->Tcs_sumy[i] + 0.5 * bd->Tv * pi->cs_sumy[i]) * bd->csy +
                             (bd->Kn + bd->Tv) * pow(bd->csy, 2);
            /*d2udyidzi, 2 3*/
            du2dxyz[1][2] += (bd->Kn * bd->csz + 0.5 * pi->Tcs_sumz[i] + 0.5 * bd->Tv * pi->cs_sumz[i]) * bd->csy +
                             (bd->Kn + bd->Tv) * bd->csz * bd->csy;
            /*d2udzidzi, 3 3*/
            du2dxyz[2][2] += (bd->Kn * bd->csz + 0.5 * pi->Tcs_sumz[i] + 0.5 * bd->Tv * pi->cs_sumz[i]) * bd->csz +
                             (bd->Kn + bd->Tv) * pow(bd->csz, 2);
        }
    }

    /*d2udyidxi, 2 1*/
    du2dxyz[1][0] = du2dxyz[1][0];
    /*d2udzidxi, 3 1*/
    du2dxyz[2][0] = du2dxyz[0][2];
    /*d2udzidyi, 3 2*/
    du2dxyz[2][1] = du2dxyz[1][2];

    return du2dxyz;
}

template <int nlayer>
std::array<std::array<double, NDIM>, NDIM> fdu2dxyz1(Particle<nlayer> *pi, Particle<nlayer> *pj) /* The off diagonal part */
{
    std::array<std::array<double, NDIM>, NDIM> du2dxyz1{0};

    // find if pj is in pi's first neighbor list
    for (Bond<nlayer> *bd1 : pi->bond_layers[0])
    {
        if (bd1->p2->id == pj->id)
        {
            for (Bond<nlayer> *bd2 : pi->bond_layers[0])
            {
                if (bd1->p2->id == bd2->p2->id)
                {
                    /*d2udxidxj, 1 1*/
                    du2dxyz1[0][0] += (-bd2->Kn * bd2->csx + 0.5 * bd2->p2->Tcs_sumx[0] + 0.5 * bd2->Tv * bd2->p2->cs_sumx[0]) * bd2->csx -
                                      (bd2->Kn + bd2->Tv) * pow(bd2->csx, 2);
                    /*d2udxidyj, 1 2*/
                    du2dxyz1[0][1] += (-bd2->Kn * bd2->csy + 0.5 * bd2->p2->Tcs_sumy[0] + 0.5 * bd2->Tv * bd2->p2->cs_sumy[0]) * bd2->csx -
                                      (bd2->Kn + bd2->Tv) * bd2->csy * bd2->csx;
                    /*d2udxidzj, 1 3*/
                    du2dxyz1[0][2] += (-bd2->Kn * bd2->csz + 0.5 * bd2->p2->Tcs_sumz[0] + 0.5 * bd2->Tv * bd2->p2->cs_sumz[0]) * bd2->csx -
                                      (bd2->Kn + bd2->Tv) * bd2->csz * bd2->csx;
                    /*d2udyidxj, 2 1*/
                    du2dxyz1[1][0] += (-bd2->Kn * bd2->csx + 0.5 * bd2->p2->Tcs_sumx[0] + 0.5 * bd2->Tv * bd2->p2->cs_sumx[0]) * bd2->csy -
                                      (bd2->Kn + bd2->Tv) * bd2->csy * bd2->csx;
                    /*d2udyidyj, 2 2*/
                    du2dxyz1[1][1] += (-bd2->Kn * bd2->csy + 0.5 * bd2->p2->Tcs_sumy[0] + 0.5 * bd2->Tv * bd2->p2->cs_sumy[0]) * bd2->csy -
                                      (bd2->Kn + bd2->Tv) * pow(bd2->csy, 2);
                    /*d2udyidzj, 2 3*/
                    du2dxyz1[1][2] += (-bd2->Kn * bd2->csz + 0.5 * bd2->p2->Tcs_sumz[0] + 0.5 * bd2->Tv * bd2->p2->cs_sumz[0]) * bd2->csy -
                                      (bd2->Kn + bd2->Tv) * bd2->csz * bd2->csy;
                    /*d2udzidxj, 3 1*/
                    du2dxyz1[2][0] += (-bd2->Kn * bd2->csx + 0.5 * bd2->p2->Tcs_sumx[0] + 0.5 * bd2->Tv * bd2->p2->cs_sumx[0]) * bd2->csz -
                                      (bd2->Kn + bd2->Tv) * bd2->csz * bd2->csx;
                    /*d2udzidyj, 3 2*/
                    du2dxyz1[2][1] += (-bd2->Kn * bd2->csy + 0.5 * bd2->p2->Tcs_sumy[0] + 0.5 * bd2->Tv * bd2->p2->cs_sumy[0]) * bd2->csz -
                                      (bd2->Kn + bd2->Tv) * bd2->csz * bd2->csy;
                    /*d2udzidzj, 3 3*/
                    du2dxyz1[2][2] += (-bd2->Kn * bd2->csz + 0.5 * bd2->p2->Tcs_sumz[0] + 0.5 * bd2->Tv * bd2->p2->cs_sumz[0]) * bd2->csz -
                                      (bd2->Kn + bd2->Tv) * pow(bd2->csz, 2);
                }
                else
                {
                    du2dxyz1[0][0] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csx * bd1->csx; /*d2udxidxj, 1 1*/
                    du2dxyz1[0][1] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csx * bd1->csy; /*d2udxidyj, 1 2*/
                    du2dxyz1[0][2] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csx * bd1->csz; /*d2udxidzj, 1 3*/
                    du2dxyz1[1][0] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csy * bd1->csx; /*d2udyidxj, 2 1*/
                    du2dxyz1[1][1] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csy * bd1->csy; /*d2udyidyj, 2 2*/
                    du2dxyz1[1][2] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csy * bd1->csz; /*d2udyidzj, 2 3*/
                    du2dxyz1[2][0] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csz * bd1->csx; /*d2udzidxj, 3 1*/
                    du2dxyz1[2][1] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csz * bd1->csy; /*d2udzidyj, 3 2*/
                    du2dxyz1[2][2] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csz * bd1->csz; /*d2udzidzj, 3 3*/

                    for (Bond<nlayer> *bd3 : bd2->p2->bond_layers[0])
                    {
                        if (bd3->p2->id == pj->id)
                        {
                            du2dxyz1[0][0] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csx * bd3->csx; /*d2udxidxj, 1 1*/
                            du2dxyz1[0][1] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csx * bd3->csy; /*d2udxidyj, 1 2*/
                            du2dxyz1[0][2] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csx * bd3->csz; /*d2udxidzj, 1 3*/
                            du2dxyz1[1][0] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csy * bd3->csx; /*d2udyidxj, 2 1*/
                            du2dxyz1[1][1] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csy * bd3->csy; /*d2udyidyj, 2 2*/
                            du2dxyz1[1][2] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csy * bd3->csz; /*d2udyidzj, 2 3*/
                            du2dxyz1[2][0] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csz * bd3->csx; /*d2udzidxj, 3 1*/
                            du2dxyz1[2][1] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csz * bd3->csy; /*d2udzidyj, 3 2*/
                            du2dxyz1[2][2] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csz * bd3->csz; /*d2udzidzj, 3 3*/
                        }
                    }
                }
            }
            return du2dxyz1;
        }
    }

    for (Bond<nlayer> *bd1 : pi->bond_layers[0])
    {
        for (Bond<nlayer> *bd2 : bd1->p2->bond_layers[0])
        {
            if (bd2->p2->id == pj->id)
            {
                du2dxyz1[0][0] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csx * bd2->csx; /*d2udxidxj, 1 1*/
                du2dxyz1[0][1] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csx * bd2->csy; /*d2udxidyj, 1 2*/
                du2dxyz1[0][2] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csx * bd2->csz; /*d2udxidzj, 1 3*/
                du2dxyz1[1][0] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csy * bd2->csx; /*d2udyidxj, 2 1*/
                du2dxyz1[1][1] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csy * bd2->csy; /*d2udyidyj, 2 2*/
                du2dxyz1[1][2] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csy * bd2->csz; /*d2udyidzj, 2 3*/
                du2dxyz1[2][0] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csz * bd2->csx; /*d2udzidxj, 3 1*/
                du2dxyz1[2][1] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csz * bd2->csy; /*d2udzidyj, 3 2*/
                du2dxyz1[2][2] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csz * bd2->csz; /*d2udzidzj, 3 3*/
            }
        }
    }

    return du2dxyz1;
}

template <int nlayer>
std::array<std::array<double, NDIM>, NDIM> fdu2dxyz2(Particle<nlayer> *pi, Particle<nlayer> *pj) /* The off diagonal part */
{
    std::array<std::array<double, NDIM>, NDIM> du2dxyz2{0};

    // find if pj is in pi's first neighbor list
    for (Bond<nlayer> *bd1 : pi->bond_layers[1])
    {
        if (bd1->p2->id == pj->id)
        {
            for (Bond<nlayer> *bd2 : pi->bond_layers[1])
            {
                if (bd1->p2->id == bd2->p2->id)
                {
                    /*d2udxidxj, 1 1*/
                    du2dxyz2[0][0] += (-bd2->Kn * bd2->csx + 0.5 * bd2->p2->Tcs_sumx[1] + 0.5 * bd2->Tv * bd2->p2->cs_sumx[1]) * bd2->csx -
                                      (bd2->Kn + bd2->Tv) * pow(bd2->csx, 2);
                    /*d2udxidyj, 1 2*/
                    du2dxyz2[0][1] += (-bd2->Kn * bd2->csy + 0.5 * bd2->p2->Tcs_sumy[1] + 0.5 * bd2->Tv * bd2->p2->cs_sumy[1]) * bd2->csx -
                                      (bd2->Kn + bd2->Tv) * bd2->csy * bd2->csx;
                    /*d2udxidzj, 1 3*/
                    du2dxyz2[0][2] += (-bd2->Kn * bd2->csz + 0.5 * bd2->p2->Tcs_sumz[1] + 0.5 * bd2->Tv * bd2->p2->cs_sumz[1]) * bd2->csx -
                                      (bd2->Kn + bd2->Tv) * bd2->csz * bd2->csx;
                    /*d2udyidxj, 2 1*/
                    du2dxyz2[1][0] += (-bd2->Kn * bd2->csx + 0.5 * bd2->p2->Tcs_sumx[1] + 0.5 * bd2->Tv * bd2->p2->cs_sumx[1]) * bd2->csy -
                                      (bd2->Kn + bd2->Tv) * bd2->csy * bd2->csx;
                    /*d2udyidyj, 2 2*/
                    du2dxyz2[1][1] += (-bd2->Kn * bd2->csy + 0.5 * bd2->p2->Tcs_sumy[1] + 0.5 * bd2->Tv * bd2->p2->cs_sumy[1]) * bd2->csy -
                                      (bd2->Kn + bd2->Tv) * pow(bd2->csy, 2);
                    /*d2udyidzj, 2 3*/
                    du2dxyz2[1][2] += (-bd2->Kn * bd2->csz + 0.5 * bd2->p2->Tcs_sumz[1] + 0.5 * bd2->Tv * bd2->p2->cs_sumz[1]) * bd2->csy -
                                      (bd2->Kn + bd2->Tv) * bd2->csz * bd2->csy;
                    /*d2udzidxj, 3 1*/
                    du2dxyz2[2][0] += (-bd2->Kn * bd2->csx + 0.5 * bd2->p2->Tcs_sumx[1] + 0.5 * bd2->Tv * bd2->p2->cs_sumx[1]) * bd2->csz -
                                      (bd2->Kn + bd2->Tv) * bd2->csz * bd2->csx;
                    /*d2udzidyj, 3 2*/
                    du2dxyz2[2][1] += (-bd2->Kn * bd2->csy + 0.5 * bd2->p2->Tcs_sumy[1] + 0.5 * bd2->Tv * bd2->p2->cs_sumy[1]) * bd2->csz -
                                      (bd2->Kn + bd2->Tv) * bd2->csz * bd2->csy;
                    /*d2udzidzj, 3 3*/
                    du2dxyz2[2][2] += (-bd2->Kn * bd2->csz + 0.5 * bd2->p2->Tcs_sumz[1] + 0.5 * bd2->Tv * bd2->p2->cs_sumz[1]) * bd2->csz -
                                      (bd2->Kn + bd2->Tv) * pow(bd2->csz, 2);
                }
                else
                {
                    du2dxyz2[0][0] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csx * bd1->csx; /*d2udxidxj, 1 1*/
                    du2dxyz2[0][1] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csx * bd1->csy; /*d2udxidyj, 1 2*/
                    du2dxyz2[0][2] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csx * bd1->csz; /*d2udxidzj, 1 3*/
                    du2dxyz2[1][0] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csy * bd1->csx; /*d2udyidxj, 2 1*/
                    du2dxyz2[1][1] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csy * bd1->csy; /*d2udyidyj, 2 2*/
                    du2dxyz2[1][2] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csy * bd1->csz; /*d2udyidzj, 2 3*/
                    du2dxyz2[2][0] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csz * bd1->csx; /*d2udzidxj, 3 1*/
                    du2dxyz2[2][1] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csz * bd1->csy; /*d2udzidyj, 3 2*/
                    du2dxyz2[2][2] -= 0.5 * (bd1->Tv + bd2->Tv) * bd2->csz * bd1->csz; /*d2udzidzj, 3 3*/

                    for (Bond<nlayer> *bd3 : bd2->p2->bond_layers[1])
                    {
                        if (bd3->p2->id == pj->id)
                        {
                            du2dxyz2[0][0] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csx * bd3->csx; /*d2udxidxj, 1 1*/
                            du2dxyz2[0][1] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csx * bd3->csy; /*d2udxidyj, 1 2*/
                            du2dxyz2[0][2] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csx * bd3->csz; /*d2udxidzj, 1 3*/
                            du2dxyz2[1][0] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csy * bd3->csx; /*d2udyidxj, 2 1*/
                            du2dxyz2[1][1] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csy * bd3->csy; /*d2udyidyj, 2 2*/
                            du2dxyz2[1][2] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csy * bd3->csz; /*d2udyidzj, 2 3*/
                            du2dxyz2[2][0] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csz * bd3->csx; /*d2udzidxj, 3 1*/
                            du2dxyz2[2][1] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csz * bd3->csy; /*d2udzidyj, 3 2*/
                            du2dxyz2[2][2] -= 0.5 * (bd3->Tv + bd2->Tv) * bd2->csz * bd3->csz; /*d2udzidzj, 3 3*/
                        }
                    }
                }
            }
            return du2dxyz2;
        }
    }

    for (Bond<nlayer> *bd1 : pi->bond_layers[1])
    {
        for (Bond<nlayer> *bd2 : bd1->p2->bond_layers[1])
        {
            if (bd2->p2->id == pj->id)
            {
                du2dxyz2[0][0] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csx * bd2->csx; /*d2udxidxj, 1 1*/
                du2dxyz2[0][1] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csx * bd2->csy; /*d2udxidyj, 1 2*/
                du2dxyz2[0][2] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csx * bd2->csz; /*d2udxidzj, 1 3*/
                du2dxyz2[1][0] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csy * bd2->csx; /*d2udyidxj, 2 1*/
                du2dxyz2[1][1] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csy * bd2->csy; /*d2udyidyj, 2 2*/
                du2dxyz2[1][2] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csy * bd2->csz; /*d2udyidzj, 2 3*/
                du2dxyz2[2][0] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csz * bd2->csx; /*d2udzidxj, 3 1*/
                du2dxyz2[2][1] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csz * bd2->csy; /*d2udzidyj, 3 2*/
                du2dxyz2[2][2] -= 0.5 * (bd2->Tv + bd1->Tv) * bd1->csz * bd2->csz; /*d2udzidzj, 3 3*/
            }
        }
    }

    return du2dxyz2;
}

// template <int nlayer>
// std::array<std::array<double, NDIM>, NDIM> fdu2dxy(int i, double *du2dxy)
// {
//     /* the diagonal part */

//     memset(du2dxy, 0, dim * dim * sizeof(double)); /* set the temporary array as zero-like */

//     int nb1 = countNEqual(neighbors1[i], nneighbors1, -1);
//     for (int k = 0; k < nb1; k++)
//     {
//         /*d2udxidxi, 1 1*/
//         du2dxy[0] += (Kn[i][k][0] * csx[i][k][0] + 0.5 * Tcs_sumx[i][0] + 0.5 * Tv[i][k][0] * cs_sumx[i][0]) * csx[i][k][0] +
//                      (Kn[i][k][0] + Tv[i][k][0]) * pow(csx[i][k][0], 2);
//         /*d2udxidyi, 1 2*/
//         du2dxy[1] += (Kn[i][k][0] * csy[i][k][0] + 0.5 * Tcs_sumy[i][0] + 0.5 * Tv[i][k][0] * cs_sumy[i][0]) * csx[i][k][0] +
//                      (Kn[i][k][0] + Tv[i][k][0]) * csy[i][k][0] * csx[i][k][0];
//         /*d2udyidyi, 2 2*/
//         du2dxy[3] += (Kn[i][k][0] * csy[i][k][0] + 0.5 * Tcs_sumy[i][0] + 0.5 * Tv[i][k][0] * cs_sumy[i][0]) * csy[i][k][0] +
//                      (Kn[i][k][0] + Tv[i][k][0]) * pow(csy[i][k][0], 2);
//     }

//     int nb2 = countNEqual(neighbors2[i], nneighbors2, -1);
//     for (int k = 0; k < nb2; k++)
//     {
//         /*d2udxidxi, 1 1*/
//         du2dxy[0] += (Kn[i][k][1] * csx[i][k][1] + 0.5 * Tcs_sumx[i][1] + 0.5 * Tv[i][k][1] * cs_sumx[i][1]) * csx[i][k][1] +
//                      (Kn[i][k][1] + Tv[i][k][1]) * pow(csx[i][k][1], 2);
//         /*d2udxidyi, 1 2*/
//         du2dxy[1] += (Kn[i][k][1] * csy[i][k][1] + 0.5 * Tcs_sumy[i][1] + 0.5 * Tv[i][k][1] * cs_sumy[i][1]) * csx[i][k][1] +
//                      (Kn[i][k][1] + Tv[i][k][1]) * csy[i][k][1] * csx[i][k][1];
//         /*d2udyidyi, 2 2*/
//         du2dxy[3] += (Kn[i][k][1] * csy[i][k][1] + 0.5 * Tcs_sumy[i][1] + 0.5 * Tv[i][k][1] * cs_sumy[i][1]) * csy[i][k][1] +
//                      (Kn[i][k][1] + Tv[i][k][1]) * pow(csy[i][k][1], 2);
//     }
//     /*d2udyidxi, 2 1*/
//     du2dxy[2] = du2dxy[1];
// }

// template <int nlayer>
// std::array<std::array<double, NDIM>, NDIM> fdu2dxy1(int i, int j, double *du2dxy1)
// {
//     /* set the temporary array as zero-like */
//     memset(du2dxy1, 0, dim * dim * sizeof(double));

//     int nb1 = countNEqual(neighbors1[i], nneighbors1, -1);
//     /* The first and second layers of particles in AFEM, i.e. the neighbors */
//     for (int kk = 0; kk < nb1; kk++)
//     {
//         if (neighbors1[i][kk] == conn[i][j])
//         {
//             for (int k = 0; k < nb1; k++)
//             {
//                 if (k == kk)
//                 {
//                     /*d2udxidxj, 1 1*/
//                     du2dxy1[0] += (-Kn[i][k][0] * csx[i][k][0] + 0.5 * Tcs_sumx[neighbors1[i][k]][0] + 0.5 * Tv[i][k][0] * cs_sumx[neighbors1[i][k]][0]) * csx[i][k][0] -
//                                   (Kn[i][k][0] + Tv[i][k][0]) * pow(csx[i][k][0], 2);
//                     /*d2udxidyj, 1 2*/
//                     du2dxy1[1] += (-Kn[i][k][0] * csy[i][k][0] + 0.5 * Tcs_sumy[neighbors1[i][k]][0] + 0.5 * Tv[i][k][0] * cs_sumy[neighbors1[i][k]][0]) * csx[i][k][0] -
//                                   (Kn[i][k][0] + Tv[i][k][0]) * csy[i][k][0] * csx[i][k][0];
//                     /*d2udyidxj, 2 1*/
//                     du2dxy1[2] += (-Kn[i][k][0] * csx[i][k][0] + 0.5 * Tcs_sumx[neighbors1[i][k]][0] + 0.5 * Tv[i][k][0] * cs_sumx[neighbors1[i][k]][0]) * csy[i][k][0] -
//                                   (Kn[i][k][0] + Tv[i][k][0]) * csy[i][k][0] * csx[i][k][0];
//                     /*d2udyidyj, 2 2*/
//                     du2dxy1[3] += (-Kn[i][k][0] * csy[i][k][0] + 0.5 * Tcs_sumy[neighbors1[i][k]][0] + 0.5 * Tv[i][k][0] * cs_sumy[neighbors1[i][k]][0]) * csy[i][k][0] -
//                                   (Kn[i][k][0] + Tv[i][k][0]) * pow(csy[i][k][0], 2);
//                 }
//                 else
//                 {
//                     du2dxy1[0] -= 0.5 * (Tv[i][kk][0] + Tv[i][k][0]) * csx[i][k][0] * csx[i][kk][0]; /*d2udxidxj, 1 1*/
//                     du2dxy1[1] -= 0.5 * (Tv[i][kk][0] + Tv[i][k][0]) * csx[i][k][0] * csy[i][kk][0]; /*d2udxidyj, 1 2*/
//                     du2dxy1[2] -= 0.5 * (Tv[i][kk][0] + Tv[i][k][0]) * csy[i][k][0] * csx[i][kk][0]; /*d2udyidxj, 2 1*/
//                     du2dxy1[3] -= 0.5 * (Tv[i][kk][0] + Tv[i][k][0]) * csy[i][k][0] * csy[i][kk][0]; /*d2udyidyj, 2 2*/

//                     int nb2 = countNEqual(neighbors1[neighbors1[i][k]], nneighbors1, -1);
//                     for (int n = 0; n < nb2; n++)
//                     {
//                         if (neighbors1[neighbors1[i][k]][n] == conn[i][j])
//                         {
//                             du2dxy1[0] -= 0.5 * (Tv[neighbors1[i][k]][n][0] + Tv[i][k][0]) * csx[i][k][0] * csx[neighbors1[i][k]][n][0]; /*d2udxidxj, 1 1*/
//                             du2dxy1[1] -= 0.5 * (Tv[neighbors1[i][k]][n][0] + Tv[i][k][0]) * csx[i][k][0] * csy[neighbors1[i][k]][n][0]; /*d2udxidyj, 1 2*/
//                             du2dxy1[2] -= 0.5 * (Tv[neighbors1[i][k]][n][0] + Tv[i][k][0]) * csy[i][k][0] * csx[neighbors1[i][k]][n][0]; /*d2udyidxj, 2 1*/
//                             du2dxy1[3] -= 0.5 * (Tv[neighbors1[i][k]][n][0] + Tv[i][k][0]) * csy[i][k][0] * csy[neighbors1[i][k]][n][0]; /*d2udyidyj, 2 2*/
//                         }
//                     }
//                 }
//             }
//             goto end2d1;
//         }
//     }
//     /* The third and fourth nearest neighbor in AFEM */
//     for (int kk = 0; kk < nb1; kk++)
//     {
//         int nb2 = countNEqual(neighbors1[neighbors1[i][kk]], nneighbors1, -1);
//         for (int m = 0; m < nb2; m++)
//         {
//             if (neighbors1[neighbors1[i][kk]][m] == conn[i][j])
//             {
//                 du2dxy1[0] -= 0.5 * (Tv[neighbors1[i][kk]][m][0] + Tv[i][kk][0]) * csx[i][kk][0] * csx[neighbors1[i][kk]][m][0]; /*d2udxidxj, 1 1*/
//                 du2dxy1[1] -= 0.5 * (Tv[neighbors1[i][kk]][m][0] + Tv[i][kk][0]) * csx[i][kk][0] * csy[neighbors1[i][kk]][m][0]; /*d2udxidyj, 1 2*/
//                 du2dxy1[2] -= 0.5 * (Tv[neighbors1[i][kk]][m][0] + Tv[i][kk][0]) * csy[i][kk][0] * csx[neighbors1[i][kk]][m][0]; /*d2udyidxj, 2 1*/
//                 du2dxy1[3] -= 0.5 * (Tv[neighbors1[i][kk]][m][0] + Tv[i][kk][0]) * csy[i][kk][0] * csy[neighbors1[i][kk]][m][0]; /*d2udyidyj, 2 2*/
//             }
//         }
//     }
// end2d1:;
// }

// template <int nlayer>
// std::array<std::array<double, NDIM>, NDIM> fdu2dxy2(int i, int j, double *du2dxy2)
// {
//     /* set the temporary array as zero-like */
//     memset(du2dxy2, 0, dim * dim * sizeof(double));

//     int nb1 = countNEqual(neighbors2[i], nneighbors2, -1);
//     /* The first and second layers of particles in AFEM, i.e. the neighbors */
//     for (int kk = 0; kk < nb1; kk++)
//     {
//         if (neighbors2[i][kk] == conn[i][j])
//         {
//             for (int k = 0; k < nb1; k++)
//             {
//                 if (k == kk)
//                 {
//                     /*d2udxidxj, 1 1*/
//                     du2dxy2[0] += (-Kn[i][k][1] * csx[i][k][1] + 0.5 * Tcs_sumx[neighbors2[i][k]][1] + 0.5 * Tv[i][k][1] * cs_sumx[neighbors2[i][k]][1]) * csx[i][k][1] -
//                                   (Kn[i][k][1] + Tv[i][k][1]) * pow(csx[i][k][1], 2);
//                     /*d2udxidyj, 1 2*/
//                     du2dxy2[1] += (-Kn[i][k][1] * csy[i][k][1] + 0.5 * Tcs_sumy[neighbors2[i][k]][1] + 0.5 * Tv[i][k][1] * cs_sumy[neighbors2[i][k]][1]) * csx[i][k][1] -
//                                   (Kn[i][k][1] + Tv[i][k][1]) * csy[i][k][1] * csx[i][k][1];
//                     /*d2udyidxj, 2 1*/
//                     du2dxy2[2] += (-Kn[i][k][1] * csx[i][k][1] + 0.5 * Tcs_sumx[neighbors2[i][k]][1] + 0.5 * Tv[i][k][1] * cs_sumx[neighbors2[i][k]][1]) * csy[i][k][1] -
//                                   (Kn[i][k][1] + Tv[i][k][1]) * csy[i][k][1] * csx[i][k][1];
//                     /*d2udyidyj, 2 2*/
//                     du2dxy2[3] += (-Kn[i][k][1] * csy[i][k][1] + 0.5 * Tcs_sumy[neighbors2[i][k]][1] + 0.5 * Tv[i][k][1] * cs_sumy[neighbors2[i][k]][1]) * csy[i][k][1] -
//                                   (Kn[i][k][1] + Tv[i][k][1]) * pow(csy[i][k][1], 2);
//                 }
//                 else
//                 {
//                     du2dxy2[0] -= 0.5 * (Tv[i][kk][1] + Tv[i][k][1]) * csx[i][k][1] * csx[i][kk][1]; /*d2udxidxj, 1 1*/
//                     du2dxy2[1] -= 0.5 * (Tv[i][kk][1] + Tv[i][k][1]) * csx[i][k][1] * csy[i][kk][1]; /*d2udxidyj, 1 2*/
//                     du2dxy2[2] -= 0.5 * (Tv[i][kk][1] + Tv[i][k][1]) * csy[i][k][1] * csx[i][kk][1]; /*d2udyidxj, 2 1*/
//                     du2dxy2[3] -= 0.5 * (Tv[i][kk][1] + Tv[i][k][1]) * csy[i][k][1] * csy[i][kk][1]; /*d2udyidyj, 2 2*/

//                     int nb2 = countNEqual(neighbors2[neighbors2[i][k]], nneighbors2, -1);
//                     for (int n = 0; n < nb2; n++)
//                     {
//                         if (neighbors2[neighbors2[i][k]][n] == conn[i][j])
//                         {
//                             du2dxy2[0] -= 0.5 * (Tv[neighbors2[i][k]][n][1] + Tv[i][k][1]) * csx[i][k][1] * csx[neighbors2[i][k]][n][1]; /*d2udxidxj, 1 1*/
//                             du2dxy2[1] -= 0.5 * (Tv[neighbors2[i][k]][n][1] + Tv[i][k][1]) * csx[i][k][1] * csy[neighbors2[i][k]][n][1]; /*d2udxidyj, 1 2*/
//                             du2dxy2[2] -= 0.5 * (Tv[neighbors2[i][k]][n][1] + Tv[i][k][1]) * csy[i][k][1] * csx[neighbors2[i][k]][n][1]; /*d2udyidxj, 2 1*/
//                             du2dxy2[3] -= 0.5 * (Tv[neighbors2[i][k]][n][1] + Tv[i][k][1]) * csy[i][k][1] * csy[neighbors2[i][k]][n][1]; /*d2udyidyj, 2 2*/
//                         }
//                     }
//                 }
//             }
//             goto end2d2;
//         }
//     }
//     /* The third and fourth nearest neighbor in AFEM */
//     for (int kk = 0; kk < nb1; kk++)
//     {
//         int nb2 = countNEqual(neighbors2[neighbors2[i][kk]], nneighbors2, -1);
//         for (int m = 0; m < nb2; m++)
//         {
//             if (neighbors2[neighbors2[i][kk]][m] == conn[i][j])
//             {
//                 du2dxy2[0] -= 0.5 * (Tv[neighbors2[i][kk]][m][1] + Tv[i][kk][1]) * csx[i][kk][1] * csx[neighbors2[i][kk]][m][1]; /*d2udxidxj, 1 1*/
//                 du2dxy2[1] -= 0.5 * (Tv[neighbors2[i][kk]][m][1] + Tv[i][kk][1]) * csx[i][kk][1] * csy[neighbors2[i][kk]][m][1]; /*d2udxidyj, 1 2*/
//                 du2dxy2[2] -= 0.5 * (Tv[neighbors2[i][kk]][m][1] + Tv[i][kk][1]) * csy[i][kk][1] * csx[neighbors2[i][kk]][m][1]; /*d2udyidxj, 2 1*/
//                 du2dxy2[3] -= 0.5 * (Tv[neighbors2[i][kk]][m][1] + Tv[i][kk][1]) * csy[i][kk][1] * csy[neighbors2[i][kk]][m][1]; /*d2udyidyj, 2 2*/
//             }
//         }
//     }

// end2d2:;
// }


#endif