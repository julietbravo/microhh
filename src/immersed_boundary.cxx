/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <iostream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "finite_difference.h"
#include "immersed_boundary.h"

namespace
{
    std::string swib;
    int mblocks;
    int nblocks;
    int iblock;
    int jblock;
    int kblock;

    std::vector<int> ib_pattern;

    void set_no_penetration(double* const restrict ut, double* const restrict vt, double* const restrict wt,
                            double* const restrict u, double* const restrict v, double* const restrict w,
                            const int* const ib,
                            double* const restrict rhoref, double* const restrict rhorefh,
                            double* const restrict dzi, double* const restrict dzhi,
                            const double dxi, const double dyi,
                            const double visc,
                            const int istart, const int iend,
                            const int jstart, const int jend,
                            const int kstart, const int kend,
                            const int jj, const int kk)
    {
        using Finite_difference::O2::interp2;

        const int ii = 1;

        const double dxidxi = dxi*dxi;
        const double dyidyi = dyi*dyi;

        // Set the IB for the east, west, north, and south edges.
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
            {
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + k*kk;

                    // SET NO PENETRATION.
                    // West face.
                    if (ib[ij] & 1)
                    {
                        ut[ijk] = 0.;
                        u [ijk] = 0.;
                    }
                    // East face.
                    if (ib[ij-ii] & 2)
                    {
                        ut[ijk] = 0.;
                        u [ijk] = 0.;
                    }
                    // South face.
                    if (ib[ij] & 4)
                    {
                        vt[ijk] = 0.;
                        v [ijk] = 0.;
                    }
                    // North face.
                    if (ib[ij-jj] & 8)
                    {
                        vt[ijk] = 0.;
                        v [ijk] = 0.;
                    }

                    // SET NO SLIP.
                    // West of west face.
                    if (ib[ij+ii] & 1)
                    {
                        vt[ijk] +=
                                + interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk], v[ijk+ii]) * dxi
                                - visc * ( (v[ijk+ii] - v[ijk]) ) * dxidxi
                                + visc * ( -2.*v[ijk] ) * dxidxi;

                        wt[ijk] +=
                                + interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk], w[ijk+ii]) * dxi
                                - visc * ( (w[ijk+ii] - w[ijk]) ) * dxidxi
                                + visc * ( -2.*w[ijk] ) * dxidxi;
                    }
                    // East of east face.
                    if (ib[ij-ii] & 2)
                    {
                        vt[ijk] +=
                                - interp2(u[ijk-jj], u[ijk]) * interp2(v[ijk-ii], v[ijk]) * dxi
                                + visc * ( (v[ijk] - v[ijk-ii]) ) * dxidxi
                                - visc * ( 2.*v[ijk] ) * dxidxi;

                        wt[ijk] +=
                                - interp2(u[ijk-kk], u[ijk]) * interp2(w[ijk-ii], w[ijk]) * dxi
                                + visc * ( (w[ijk] - w[ijk-ii]) ) * dxidxi
                                - visc * ( 2.*w[ijk] ) * dxidxi;
                    }
                    // South of south face.
                    if (ib[ij+jj] & 4)
                    {
                        ut[ijk] +=
                                + interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk], u[ijk+jj]) * dyi
                                - visc * ( (u[ijk+jj] - u[ijk]) ) * dyidyi
                                + visc * ( -2.*u[ijk] ) * dyidyi;

                        wt[ijk] +=
                                + interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk], w[ijk+jj]) * dyi
                                - visc * ( (w[ijk+jj] - w[ijk]) ) * dyidyi
                                + visc * ( -2.*w[ijk] ) * dyidyi;

                    }
                    // North of north face.
                    if (ib[ij-jj] & 8)
                    {
                        ut[ijk] +=
                                - interp2(v[ijk-ii], v[ijk]) * interp2(u[ijk-jj], u[ijk]) * dyi
                                + visc * ( (u[ijk] - u[ijk-jj]) ) * dyidyi
                                - visc * ( 2.*u[ijk] ) * dyidyi;

                        wt[ijk] +=
                                - interp2(v[ijk-kk], v[ijk]) * interp2(w[ijk-jj], w[ijk]) * dyi
                                + visc * ( (w[ijk] - w[ijk-jj]) ) * dyidyi
                                - visc * ( 2.*w[ijk] ) * dyidyi;
                    }
                }
            }
        }

        // Set the top.
        for (int j=jstart; j<jend; ++j)
        {
            for (int i=istart; i<iend; ++i)
            {
                const int k = kend;
                const int ij = i + j*jj;
                const int ijk = i + j*jj + k*kk;

                const bool at_top = ib[ij] & 16;
                const bool east_of_top = !(ib[ij] & 16) && (ib[ij-ii] & 16);
                const bool north_of_top = !(ib[ij] & 16) && (ib[ij-jj] & 16);

                // Top face.
                if (at_top)
                {
                    wt[ijk] = 0.;
                    w [ijk] = 0.;
                }

                if ( at_top || ( !at_top && east_of_top) )
                {
                    ut[ijk] += - ( rhorefh[k] * interp2(w[ijk-ii], w[ijk]) * interp2(u[ijk-kk], u[ijk]) ) / rhoref[k] * dzi[k]
                               + visc * ( (u[ijk] - u[ijk-kk]) * dzhi[k] ) * dzi[k]
                               - visc * ( 2.*u[ijk] * dzhi[k] ) * dzi[k];
                }
                    
                if ( at_top || ( !at_top && north_of_top) )
                {
                    vt[ijk] += - ( rhorefh[k] * interp2(w[ijk-jj], w[ijk]) * interp2(v[ijk-kk], v[ijk]) ) / rhoref[k] * dzi[k]
                               + visc * ( (v[ijk] - v[ijk-kk]) * dzhi[k] ) * dzi[k]
                               - visc * ( 2.*v[ijk] * dzhi[k] ) * dzi[k];
                }
            }
        }
    }

    void set_scalar(double* const restrict st, const double* const restrict s,
                    const double* const restrict u, const double* const restrict v, const double* const restrict w,
                    const int* const ib,
                    double* const restrict rhoref, double* const restrict rhorefh,
                    double* const restrict dzi, double* const restrict dzhi,
                    const double dxi, const double dyi,
                    const double visc,
                    const double sbot,
                    const int istart, const int iend,
                    const int jstart, const int jend,
                    const int kstart, const int kend,
                    const int jj, const int kk)
    {
        using Finite_difference::O2::interp2;

        const int ii = 1;

        const double dxidxi = dxi*dxi;
        const double dyidyi = dyi*dyi;

        // Enforce no penetration on faces. For simplicity, we
        // also write the ib in the ghost cells on the edge. This
        // is harmless.
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
            {
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + k*kk;

                    const double s_ib = 2.*sbot - s[ijk];

                    // West face.
                    if (ib[ij+ii] & 1)
                    {
                        st[ijk] +=
                                + ( u[ijk+ii] * interp2(s[ijk], s[ijk+ii]) ) * dxi
                                - visc * ( s[ijk+ii] - s[ijk] ) * dxidxi
                                + visc * ( s_ib - s[ijk]   ) * dxidxi;
                    }
                    // East face.
                    if (ib[ij-ii] & 2)
                    {
                        const double s_ib = 2.*sbot - s[ijk];
                        st[ijk] +=
                                - ( u[ijk] * interp2(s[ijk-ii], s[ijk]) ) * dxi
                                + visc * ( s[ijk] - s[ijk-ii] ) * dxidxi
                                - visc * ( s[ijk] - s_ib ) * dxidxi;
                    }
                    // South face.
                    if (ib[ij+jj] & 4)
                    {
                        st[ijk] +=
                                + ( v[ijk+jj] * interp2(s[ijk], s[ijk+jj]) ) * dyi
                                - visc * ( s[ijk+jj] - s[ijk] ) * dyidyi
                                + visc * ( s_ib - s[ijk] ) * dyidyi;
                    }
                    // North face.
                    if (ib[ij-jj] & 8)
                    {
                        st[ijk] +=
                                - ( v[ijk] * interp2(s[ijk-jj], s[ijk]) ) * dyi
                                + visc * ( s[ijk] - s[ijk-jj] ) * dyidyi
                                - visc * ( s[ijk] - s_ib ) * dyidyi;
                    }
                }
            }
        }

        // Set the top and bottom.
        for (int j=jstart; j<jend; ++j)
        {
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*jj;


                if (ib[ij] & 16)
                {
                    // Make sure no scalar flows into the ib from the bottom.
                    int k = kstart;
                    int ijk = i + j*jj + k*kk;
                    double s_ib = 2.*sbot - s[ijk];

                    st[ijk] +=
                             - ( rhorefh[k] * w[ijk] * interp2(s[ijk-kk], s[ijk]) ) / rhoref[k] * dzi[k]
                             + visc * (s[ijk] - s[ijk-kk]) * dzhi[k] * dzi[k];

                    k = kend;
                    ijk = i + j*jj + k*kk;
                    s_ib = 2.*sbot - s[ijk];

                    st[ijk] +=
                             - ( rhorefh[k] * w[ijk] * interp2(s[ijk-kk], s[ijk]) ) / rhoref[k] * dzi[k]
                             + visc * (s[ijk] - s[ijk-kk]) * dzhi[k] * dzi[k]
                             - visc * (s[ijk] - s_ib) * dzhi[k] * dzi[k];
                }
            }
        }
    }
}

Immersed_boundary::Immersed_boundary(Master& masterin, Grid& gridin, Input& input) :
    master(masterin),
    grid(gridin)
{
    int nerror = 0;
    nerror += input.get_item(&swib, "ib", "swib", "");

    if (swib != "1")
    {
        input.flag_as_used("ib", "mblocks");
        input.flag_as_used("ib", "nblocks");
        input.flag_as_used("ib", "iblock");
        input.flag_as_used("ib", "jblock");
        input.flag_as_used("ib", "kblock");
    }
    else
    {
        nerror += input.get_item(&mblocks, "ib", "mblocks", "");
        nerror += input.get_item(&nblocks, "ib", "nblocks", "");
        nerror += input.get_item(&iblock , "ib", "iblock" , "");
        nerror += input.get_item(&jblock , "ib", "jblock" , "");
        nerror += input.get_item(&kblock , "ib", "kblock" , "");
    }

    if (nerror > 0)
        throw 1;
}

Immersed_boundary::~Immersed_boundary()
{
}

void Immersed_boundary::init()
{
    if (swib != "1")
        return;

    // Set the ib patterns to the correct size and initialize to zero.
    ib_pattern.resize(grid.ijcells);
    std::fill(ib_pattern.begin(), ib_pattern.end(), 0);
}

void Immersed_boundary::create(Input& input, Fields& fields)
{
    if (swib != "1")
        return;

    // Get the boundary conditions for each scalar.
    for (auto& s : fields.sp)
    {
        double sbot;
        int nerror = input.get_item(&sbot, "ib", "sbot", s.first);
        if (nerror)
            throw 1;
        sbot_map[s.first] = sbot;
    }

    // Set the west faces
    const int istep = grid.itot/mblocks;
    const int jstep = grid.jtot/nblocks;

    std::vector<int> ib_pattern_tmp(grid.ijcells, 0);
    
    for (int n=0; n<nblocks; ++n)
    {
        for (int m=0; m<mblocks; ++m)
        {
            // Calculate the absolute grid indices of the blocks.
            // Note that we use c-style ranging, thus end is not part of the
            // range, and end-start is the number of points in the range.
            const int iblock_start = m*istep + istep/2 - iblock/2;
            const int iblock_end = m*istep + istep/2 + iblock/2;
            const int jblock_start = n*jstep + jstep/2 - jblock/2;
            const int jblock_end = n*jstep + jstep/2 + jblock/2;

            // const int iblock_start = m*istep;
            // const int iblock_end = iblock_start + iblock;
            // const int jblock_start = n*jstep;
            // const int jblock_end = jblock_start + jblock;
            const int kblock_start = 0;
            const int kblock_end = kblock_start + kblock;

            // Check the ranges of i and j for the specific MPI process.
            const int imin_abs = master.mpicoordx*grid.imax;
            const int imax_abs = (master.mpicoordx+1)*grid.imax;
            const int jmin_abs = master.mpicoordy*grid.jmax;
            const int jmax_abs = (master.mpicoordy+1)*grid.jmax;

            // Check whether there is an edge in range.
            const bool iblock_start_in_range = (iblock_start >= imin_abs) && (iblock_start <  imax_abs);
            const bool iblock_end_in_range   = (iblock_end-1 >= imin_abs) && (iblock_end-1 <  imax_abs);
            const bool iblock_fully_in_range = (iblock_start <  imin_abs) && (iblock_end-1 >= imax_abs);

            const bool jblock_start_in_range = (jblock_start >= jmin_abs) && (jblock_start <  jmax_abs);
            const bool jblock_end_in_range   = (jblock_end-1 >= jmin_abs) && (jblock_end-1 <  jmax_abs);
            const bool jblock_fully_in_range = (jblock_start <  jmin_abs) && (jblock_end-1 >= jmax_abs);

            if ( (iblock_start_in_range || iblock_end_in_range || iblock_fully_in_range) &&
                 (jblock_start_in_range || jblock_end_in_range || jblock_fully_in_range) )
            {
                const int istart = (iblock_start_in_range ? iblock_start : imin_abs) - master.mpicoordx*grid.imax + grid.igc;
                const int iend   = (iblock_end_in_range   ? iblock_end   : imax_abs) - master.mpicoordx*grid.imax + grid.igc;
                const int jstart = (jblock_start_in_range ? jblock_start : jmin_abs) - master.mpicoordy*grid.jmax + grid.jgc;
                const int jend   = (jblock_end_in_range   ? jblock_end   : jmax_abs) - master.mpicoordy*grid.jmax + grid.jgc;

                for (int j=jstart; j<jend; ++j)
                {
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ij = i + j*grid.icells;
                        ib_pattern_tmp[ij] = 1;
                    }
                }
            }
        }
    }

    grid.boundary_cyclic_2d(ib_pattern_tmp.data());

    const int ii = 1;
    const int jj = grid.icells;

    // This can be more elegantly solved with enum. For now, OK.
    for (int j=grid.jstart; j<grid.jend; ++j)
    {
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ij = i + j*jj;
            // Flag west.
            if (ib_pattern_tmp[ij] == 1 && ib_pattern_tmp[ij-ii] == 0)
                ib_pattern[ij] |= 1;
            // Flag east.
            if (ib_pattern_tmp[ij] == 1 && ib_pattern_tmp[ij+ii] == 0)
                ib_pattern[ij] |= 2;
            // Flag south.
            if (ib_pattern_tmp[ij] == 1 && ib_pattern_tmp[ij-jj] == 0)
                ib_pattern[ij] |= 4;
            // Flag north.
            if (ib_pattern_tmp[ij] == 1 && ib_pattern_tmp[ij+jj] == 0)
                ib_pattern[ij] |= 8;
            // Flag top.
            if (ib_pattern_tmp[ij] == 1)
                ib_pattern[ij] |= 16;
        }
    }

    grid.boundary_cyclic_2d(ib_pattern.data());
}

void Immersed_boundary::exec(Fields& fields)
{
    if (swib != "1")
        return;

    set_no_penetration(fields.ut->data, fields.vt->data, fields.wt->data,
                       fields.u->data, fields.v->data, fields.w->data,
                       ib_pattern.data(),
                       fields.rhoref, fields.rhorefh,
                       grid.dzi, grid.dzhi,
                       grid.dxi, grid.dyi,
                       fields.visc,
                       grid.istart, grid.iend,
                       grid.jstart, grid.jend,
                       grid.kstart, grid.kstart+kblock,
                       grid.icells, grid.ijcells);

    grid.boundary_cyclic(fields.u ->data);
    grid.boundary_cyclic(fields.v ->data);

    for (auto& s : fields.st)
    {
        set_scalar(s.second->data, fields.sp[s.first]->data,
                   fields.u->data, fields.v->data, fields.w->data,
                   ib_pattern.data(),
                   fields.rhoref, fields.rhorefh,
                   grid.dzi, grid.dzhi,
                   grid.dxi, grid.dyi,
                   fields.visc,
                   sbot_map[s.first],
                   grid.istart, grid.iend,
                   grid.jstart, grid.jend,
                   grid.kstart, grid.kstart+kblock,
                   grid.icells, grid.ijcells);
    }
}
