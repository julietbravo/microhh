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

#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "netcdf_interface.h"
#include "constants.h"
#include "timeloop.h"
#include "trajectory.h"

namespace
{
    template<typename TF> TF netcdf_fp_fillvalue();
    template<> double netcdf_fp_fillvalue<double>() { return NC_FILL_DOUBLE; }
    template<> float  netcdf_fp_fillvalue<float>()  { return NC_FILL_FLOAT; }
}

template<typename TF>
Trajectory<TF>::Trajectory(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_trajectory = inputin.get_item<bool>("trajectory", "swtrajectory", "", false);

    if (sw_trajectory)
    {
        sampletime = inputin.get_item<double>("trajectory", "sampletime", "");
        names = inputin.get_list<std::string>("trajectory", "names", "");
        variables = inputin.get_list<std::string>("trajectory", "variables", "");
    }
}

template<typename TF>
Trajectory<TF>::~Trajectory()
{
}

template<typename TF>
void Trajectory<TF>::init(double ifactor)
{
    if (!sw_trajectory)
        return;

    isampletime = static_cast<unsigned long>(ifactor * sampletime);
    statistics_counter = 0;
}

template<typename TF>
void Trajectory<TF>::create(
        Input& inputin, Netcdf_file& input_nc,
        Timeloop<TF>& timeloop, std::string sim_name)
{
    // Do not create file if column is disabled.
    if (!sw_trajectory)
        return;

    for (auto& name : names)
    {
        // Get time and location from input NetCDF file.
        std::string group_name = "trajectory_" + name;
        Netcdf_group& group_nc = input_nc.get_group(group_name);
        int nt = group_nc.get_dimension_size("time");

        // Create new trajectory.
        auto t = Single_Trajectory<TF>();

        t.time_in.resize(nt);
        t.x_in.resize(nt);
        t.y_in.resize(nt);
        t.z_in.resize(nt);

        group_nc.get_variable(t.time_in, "time", {0}, {nt});
        group_nc.get_variable(t.x_in, "x", {0}, {nt});
        group_nc.get_variable(t.y_in, "y", {0}, {nt});
        group_nc.get_variable(t.z_in, "z", {0}, {nt});

        // Create a NetCDF file for the statistics.
        std::stringstream filename;
        filename << sim_name << "." << "trajectory_" << name << "."
                 << std::setfill('0') << std::setw(7) << timeloop.get_iotime() << ".nc";

        // Create new NetCDF file.
        t.data_file = std::make_unique<Netcdf_file>(master, filename.str(), Netcdf_mode::Create);
        t.data_file->add_dimension("time");

        // Create dimensions and variables.
        t.time_var = std::make_unique<Netcdf_variable<TF>>(
                    t.data_file->template add_variable<TF>("time", {"time"}));
        if (timeloop.has_utc_time())
            t.time_var->add_attribute("units", "seconds since " + timeloop.get_datetime_utc_start_string());
        else
            t.time_var->add_attribute("units", "seconds since start");
        t.time_var->add_attribute("long_name", "Time");

        t.x_var = std::make_unique<Netcdf_variable<TF>>(
                    t.data_file->template add_variable<TF>("x", {"time"}));
        t.x_var->add_attribute("units", "m");
        t.x_var->add_attribute("long_name", "X location trajectory");

        t.y_var = std::make_unique<Netcdf_variable<TF>>(
                    t.data_file->template add_variable<TF>("y", {"time"}));
        t.y_var->add_attribute("units", "m");
        t.y_var->add_attribute("long_name", "Y location trajectory");

        t.z_var = std::make_unique<Netcdf_variable<TF>>(
                    t.data_file->template add_variable<TF>("z", {"time"}));
        t.z_var->add_attribute("units", "m");
        t.z_var->add_attribute("long_name", "Z location trajectory");

        t.start_time = t.time_in[0];
        t.end_time = t.time_in[t.time_in.size()-1];

        for (auto& name : variables)
        {
            // Check if requested field is prognostic. If not, raise exception.
            if (!fields.ap.count(name))
            {
                std::string error = "Trajectory variable \"" + name + "\" is not a prognostic field!";
                throw std::runtime_error(error);
            }

            Time_var<TF> var{t.data_file->template add_variable<TF>(name, {"time"}), TF(0)};
            var.ncvar.add_attribute("units", fields.ap.at(name)->unit);
            var.ncvar.add_attribute("long_name", fields.ap.at(name)->longname);
            //var.ncvar.add_attribute("_FillValue", netcdf_fp_fillvalue<TF>());

            t.time_series.emplace(
                    std::piecewise_construct,
                    std::forward_as_tuple(name),
                    std::forward_as_tuple(std::move(var)));
        }

        // Sync, and append to list of trajectories.
        t.data_file->sync();

        trajectories.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(name),
            std::forward_as_tuple(std::move(t)));
    }
}

template<typename TF>
unsigned long Trajectory<TF>::get_time_limit(unsigned long itime)
{
    if (!sw_trajectory)
        return Constants::ulhuge;

    return isampletime - itime % isampletime;
}

template<typename TF>
bool Trajectory<TF>::do_trajectory(unsigned long itime)
{
    // Check if column are enabled
    if (!sw_trajectory)
        return false;

    // Check if time for execution
    if (itime % isampletime != 0)
        return false;

    // Return true such that column are computed
    return true;
}

template<typename TF>
void Trajectory<TF>::exec(Timeloop<TF>& timeloop, double time, unsigned long itime)
{
    // Write message in trajectory stats is triggered.
    master.print_message("Saving trajectory for time %f\n", time);

    auto& gd=grid.get_grid_data();
    auto& md=master.get_MPI_data();

    const std::vector<int> time_index{statistics_counter};
    const TF fill_value = netcdf_fp_fillvalue<TF>();

    for (auto& t : trajectories)
    {
        auto& tt = t.second;

        bool valid_data = false;
        int mpiid_of_traj;
        TF x_loc, y_loc, z_loc;

        // Always set time in NetCDF. Fill_values for time are not a great thing to deal with..
        tt.time_var->insert(time, time_index);

        // Check on time bounds. Trajectories can start/end during the simulation
        if (time >= tt.start_time && time <= tt.end_time)
        {
            // Get (x,y,z) location by interpolation of the input coordinates in time.
            Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(tt.time_in);

            x_loc = ifac.fac0 * tt.x_in[ifac.index0] + ifac.fac1 * tt.x_in[ifac.index1];
            y_loc = ifac.fac0 * tt.y_in[ifac.index0] + ifac.fac1 * tt.y_in[ifac.index1];
            z_loc = ifac.fac0 * tt.z_in[ifac.index0] + ifac.fac1 * tt.z_in[ifac.index1];

            // Bounds check on domain. NOTE: keep the `>=xsize` instead of `>xsize`...
            if (x_loc >= 0 && x_loc < gd.xsize &&
                y_loc >= 0 && y_loc < gd.ysize &&
                z_loc >= 0 && z_loc < gd.zsize)
            {
                valid_data = true;

                // Get indexes and interpolation factors for full and half levels.
                auto get_index_fac = [&](
                    const std::vector<TF>& x_vec, const TF x_val, const int size)
                {
                    // Get index in `x_vec` left of `x_val`.
                    int n0=0;
                    for (int n=0; n<size; ++n)
                    {
                        if (x_vec[n] <= x_val && x_vec[n+1] > x_val)
                        {
                            n0 = n;
                            break;
                        }
                    }

                    // 2. Calculate interpolation factor.
                    const TF fac = TF(1) - (x_val - x_vec[n0]) / (x_vec[n0+1] - x_vec[n0]);

                    return std::pair(n0, fac);
                };

                // Calculate `mpiid` which contains the current trajectory location.
                // All MPI tasks need to know this ID for the `broadcast` below.
                const TF sub_xsize = gd.imax * gd.dx;
                const TF sub_ysize = gd.jmax * gd.dy;

                const int mpicoordx = x_loc / sub_xsize;
                const int mpicoordy = y_loc / sub_ysize;

                mpiid_of_traj = master.calc_mpiid(mpicoordx, mpicoordy);

                if (md.mpiid == mpiid_of_traj)
                {
                    std::pair<int, TF> ipx  = get_index_fac(gd.x,  x_loc, gd.icells);
                    std::pair<int, TF> ipxh = get_index_fac(gd.xh, x_loc, gd.icells);

                    std::pair<int, TF> jpx  = get_index_fac(gd.y,  y_loc, gd.jcells);
                    std::pair<int, TF> jpxh = get_index_fac(gd.yh, y_loc, gd.jcells);

                    std::pair<int, TF> kpx  = get_index_fac(gd.z,  z_loc, gd.kcells);
                    std::pair<int, TF> kpxh = get_index_fac(gd.zh, z_loc, gd.kcells);

                    // Put in vector for easier access with `field3d->loc`.
                    std::vector<std::pair<int, TF>> ix = {ipx, ipxh};
                    std::vector<std::pair<int, TF>> iy = {jpx, jpxh};
                    std::vector<std::pair<int, TF>> iz = {kpx, kpxh};

                    // Interpolate variables.
                    auto ijk = [&](const int i, const int j, const int k)
                    {
                        return i + j*gd.icells + k*gd.ijcells;
                    };

                    auto interpolate = [&](
                            const std::vector<TF>& fld,
                            std::pair<int, TF> fx,
                            std::pair<int, TF> fy,
                            std::pair<int, TF> fz)
                    {
                        // Short-cuts.
                        const int i0 = fx.first;
                        const int j0 = fy.first;
                        const int k0 = fz.first;

                        const TF fx0 = fx.second;
                        const TF fx1 = TF(1) - fx.second;
                        const TF fy0 = fy.second;
                        const TF fy1 = TF(1) - fy.second;
                        const TF fz0 = fz.second;
                        const TF fz1 = TF(1) - fz.second;

                        // Tri-linear interpolation onto requested location.
                        const TF value =
                            fx0 * fy0 * fz0 * fld[ijk(i0,   j0,   k0  )] +
                            fx1 * fy0 * fz0 * fld[ijk(i0+1, j0,   k0  )] +
                            fx0 * fy1 * fz0 * fld[ijk(i0,   j0+1, k0  )] +
                            fx1 * fy1 * fz0 * fld[ijk(i0+1, j0+1, k0  )] +
                            fx0 * fy0 * fz1 * fld[ijk(i0,   j0,   k0+1)] +
                            fx1 * fy0 * fz1 * fld[ijk(i0+1, j0,   k0+1)] +
                            fx0 * fy1 * fz1 * fld[ijk(i0,   j0+1, k0+1)] +
                            fx1 * fy1 * fz1 * fld[ijk(i0+1, j0+1, k0+1)];

                        return value;
                    };

                    for (auto& name : variables)
                    {
                        const int loc_i = fields.ap.at(name)->loc[0];
                        const int loc_j = fields.ap.at(name)->loc[1];
                        const int loc_k = fields.ap.at(name)->loc[2];

                        tt.time_series.at(name).data = interpolate(
                                fields.ap.at(name)->fld, ix[loc_i], iy[loc_j], iz[loc_k]);
                    }
                } // mpiid == mpiid_of_traj
            } // xyz bounds
        } // time bounds

        if (valid_data)
        {
            tt.x_var->insert(x_loc, time_index);
            tt.y_var->insert(y_loc, time_index);
            tt.z_var->insert(z_loc, time_index);

            // This is a bit wasteful (?), a specific send/recv should be enough,
            // but we are only sending one float/double per variable....
            for (auto& name : variables)
                master.broadcast(&tt.time_series.at(name).data, 1, mpiid_of_traj);

            // Store value in NetCDF file.
            for (auto& name : variables)
                tt.time_series.at(name).ncvar.insert(tt.time_series.at(name).data, time_index);
        }
        else
        {
            // Either time or location out of bounds.
            tt.x_var->insert(fill_value, time_index);
            tt.y_var->insert(fill_value, time_index);
            tt.z_var->insert(fill_value, time_index);

            // Store (fill) value in NetCDF file.
            for (auto& name : variables)
                tt.time_series.at(name).ncvar.insert(fill_value, time_index);
        }

        // Synchronize the NetCDF file
        tt.data_file->sync();
    } // trajectories

    ++statistics_counter;
}

template class Trajectory<double>;
template class Trajectory<float>;
