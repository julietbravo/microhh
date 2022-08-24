/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#include <algorithm>
#include <cmath>

#include "canopy.h"
#include "input.h"
#include "field3d_operators.h"
#include "grid.h"
#include "netcdf_interface.h"
#include "fast_math.h"
#include "constants.h"
#include "fields.h"

namespace
{
    namespace fm = Fast_math;

    template<typename TF>
    void canopy_drag_u(
            TF* const restrict ut,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict pad,
            const TF utrans,
            const TF vtrans,
            const TF cd,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Interpolate `v` and `w` to `u` locations.
                    const TF u_on_u = u[ijk] + utrans;
                    const TF v_on_u = TF(0.25) * (v[ijk] + v[ijk-ii] + v[ijk-ii+jj] + v[ijk+jj]) + vtrans;
                    const TF w_on_u = TF(0.25) * (w[ijk] + w[ijk-ii] + w[ijk-ii+kk] + w[ijk+kk]);

                    const TF ftau = -cd * pad[k] *
                        std::pow( fm::pow2(u_on_u) +
                                  fm::pow2(v_on_u) +
                                  fm::pow2(w_on_u), TF(0.5) );

                    ut[ijk] += ftau * u_on_u;
                }
    }

    template<typename TF>
    void canopy_drag_v(
            TF* const restrict vt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict pad,
            const TF utrans,
            const TF vtrans,
            const TF cd,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Interpolate `u` and `w` to `v` locations.
                    const TF u_on_v = TF(0.25) * (u[ijk] + u[ijk+ii] + u[ijk+ii-jj] + v[ijk-jj]) + utrans;
                    const TF v_on_v = v[ijk] + vtrans;
                    const TF w_on_v = TF(0.25) * (w[ijk] + w[ijk+kk] + w[ijk+kk-kk] + w[ijk-jj]);

                    const TF ftau = -cd * pad[k] *
                        std::pow( fm::pow2(u_on_v) +
                                  fm::pow2(v_on_v) +
                                  fm::pow2(w_on_v), TF(0.5) );

                    vt[ijk] += ftau * v_on_v;
                }
    }

    template<typename TF>
    void canopy_drag_w(
            TF* const restrict wt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict padh,
            const TF utrans,
            const TF vtrans,
            const TF cd,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Interpolate `u` and `v` to `w` locations.
                    const TF u_on_w = TF(0.25) * (u[ijk] + u[ijk+ii] + u[ijk+ii-kk] + v[ijk-kk]) + utrans;
                    const TF v_on_w = TF(0.25) * (v[ijk] + v[ijk+jj] + v[ijk+jj-kk] + v[ijk-kk]) + vtrans;
                    const TF w_on_w = w[ijk];

                    const TF ftau = -cd * padh[k] *
                        std::pow( fm::pow2(u_on_w) +
                                  fm::pow2(v_on_w) +
                                  fm::pow2(w_on_w), TF(0.5) );

                    wt[ijk] += ftau * w_on_w;
                }
    }
}

template<typename TF>
Canopy<TF>::Canopy(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin), field3d_operators(masterin, gridin, fieldsin)
{
    sw_canopy = inputin.get_item<bool>("canopy", "sw_canopy", "", false);
    cd        = inputin.get_item<TF>("canopy", "cd", "");
}

template <typename TF>
Canopy<TF>::~Canopy()
{
}

template <typename TF>
void Canopy<TF>::init()
{
    if (!sw_canopy)
        return;

    auto& gd = grid.get_grid_data();

    if (sw_canopy)
    {
        pad.resize (gd.kcells);
        padh.resize(gd.kcells);
    }
}

template <typename TF>
void Canopy<TF>::create(
        Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    if (!sw_canopy)
        return;

    auto& gd = grid.get_grid_data();

    // Read plant area density at half levels.
    Netcdf_group& group_nc = input_nc.get_group("init");
    group_nc.get_variable(padh, "padh", {0}, {gd.ktot+1});
    std::rotate(padh.rbegin(), padh.rbegin() + gd.kstart, padh.rend());

    // Interpolate plant area density to full levels.
    for (int k=gd.kstart; k<gd.kend; k++)
        pad[k] = TF(0.5)*(padh[k] + padh[k+1]);

    // Determine end of canopy.
    for (int k=gd.kend-1; k>gd.kstart; k--)
        if (pad[k] < Constants::dsmall && pad[k-1] >= Constants::dsmall)
        {
            kend_canopy = k;
            break;
        }
}

//#ifndef USECUDA
template <typename TF>
void Canopy<TF>::exec()
{
    if (!sw_canopy)
        return;

    auto& gd = grid.get_grid_data();

    // Momentum drag
    canopy_drag_u(
        fields.mt.at("u")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        pad.data(),
        grid.utrans,
        grid.vtrans,
        cd,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, kend_canopy,
        gd.icells, gd.ijcells);

    canopy_drag_v(
        fields.mt.at("v")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        pad.data(),
        grid.utrans,
        grid.vtrans,
        cd,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, kend_canopy,
        gd.icells, gd.ijcells);

    canopy_drag_w(
        fields.mt.at("w")->fld.data(),
        fields.mp.at("u")->fld.data(),
        fields.mp.at("v")->fld.data(),
        fields.mp.at("w")->fld.data(),
        padh.data(),
        grid.utrans,
        grid.vtrans,
        cd,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, kend_canopy,
        gd.icells, gd.ijcells);
}
//#endif

template class Canopy<double>;
template class Canopy<float>;
