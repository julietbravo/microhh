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

#include "canopy.h"
#include "input.h"
#include "field3d_operators.h"
#include "grid.h"

namespace
{
}

template<typename TF>
Canopy<TF>::Canopy(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin), field3d_operators(masterin, gridin, fieldsin)
{
    bool sw_canopy = inputin.get_item<bool>("canopy", "sw_canopy", "", false);
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
        //bla.resize(gd.kcells);
    }
}

template <typename TF>
void Canopy<TF>::create(
        Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    if (!sw_canopy)
        return;

    auto& gd = grid.get_grid_data();
}

//#ifndef USECUDA
template <typename TF>
void Canopy<TF>::exec()
{
    if (!sw_canopy)
        return;

    auto& gd = grid.get_grid_data();
}
//#endif

template class Canopy<double>;
template class Canopy<float>;
