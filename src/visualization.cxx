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

#include <iostream>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "visualization.h"

template<typename TF>
Visualization<TF>::Visualization(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    std::cout << "construct visu" << std::endl;
}

template <typename TF>
Visualization<TF>::~Visualization()
{
}

template <typename TF>
void Visualization<TF>::init()
{
    std::cout << "init visu" << std::endl;
}

template <typename TF>
void Visualization<TF>::exec(Timeloop<TF>& timeloop)
{
    if (!timeloop.in_substep())
        std::cout << "exec visu" << std::endl;
}

template class Visualization<double>;
template class Visualization<float>;
