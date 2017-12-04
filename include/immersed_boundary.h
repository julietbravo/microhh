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

#ifndef IMMERSED_BOUNDARY
#define IMMERSED_BOUNDARY

#include <vector>
#include <bitset>

class Model;
class Grid;
class Fields;
class Input;
class Stats;
struct Mask;

struct Neighbour
{
    int i;
    int j;
    int k;
    int ijk;
    double distance;
};

struct Ghost_cell
{
    int i;      ///< x-index of ghost cell
    int j;      ///< y-index of ghost cell
    int k;      ///< z-index of ghost cell
    int ijk;    ///< combined index

    double xB;  ///< Nearest x location on immersed boundary
    double yB;  ///< Nearest y location on immersed boundary
    double zB;  ///< Nearest z location on immersed boundary

    double xI;  ///< x-Location of image point of ghost cell
    double yI;  ///< x-Location of image point of ghost cell
    double zI;  ///< x-Location of image point of ghost cell

    double dI;  ///< Distance from ghost cell to interpolation point

    std::vector< std::vector<double> > B; ///< Inversed matrix with the locations of all interpolation points
    std::vector<Neighbour> neighbours;    ///< Neighbouring fluid points used in interpolation

    std::vector<double> c_idw; ///< Weights for inverse distance weighted interpolation
    double c_idw_sum;          ///< Sum inverse distance weights
};

class Immersed_boundary
{
    public:
        enum IB_type       {None_type, Poly_type, Dem_type};
        enum Boundary_type {Dirichlet_type, Neumann_type, Flux_type};

        Immersed_boundary(Model*, Input*); ///< Constructor of the class.
        ~Immersed_boundary();              ///< Destructor of the class.

        void init();
        void create();
        void exec_momentum();   ///< Set the immersed boundary ghost cells (velocity components)
        void exec_scalars();    ///< Set the immersed boundary ghost cells (scalars)
        void exec_stats(Mask*); ///< Execute statistics of immersed boundaries
        void get_mask(Field3d*, Field3d*);
        IB_type get_switch() { return ib_type; }

    private:
        Model*  model;  ///< Pointer to model class.
        Fields* fields; ///< Pointer to fields class.
        Grid*   grid;   ///< Pointer to grid class.
        Stats*  stats;  ///< Pointer to grid class.

        ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Ghost_cell> ghost_cells_u;  ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Ghost_cell> ghost_cells_v;  ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Ghost_cell> ghost_cells_w;  ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Ghost_cell> ghost_cells_s;  ///< Vector holding info on all the ghost cells within the boundary

        void read_ghost_cells(std::vector<Ghost_cell>&, std::string, const double*, const double*, const double*, Boundary_type); ///< Function to read user input IB

        void find_ghost_cells(std::vector<Ghost_cell>&, const double*, const double*, const double*, const int, Boundary_type); ///< Function which determines the ghost cells
        bool is_ghost_cell(const double*, const double*, const double*, const int, const int, const int); ///< Function which checks if a cell is a ghost cell
        void find_nearest_location_wall(double&, double&, double&, double&,
                                        const double, const double, const double,
                                        const int, const int, const int); ///< Function which checks if a cell is a ghost cell
        void find_interpolation_points(Ghost_cell&, const double*, const double*, const double*, const int, const int, const int, Boundary_type); ///< Function which searched for the nearest neighbours
        void calc_mask(double*, double*, double*, int*, int*, int*, const double*, const double*, const double*, const double*);
        void zero_ib_tendency(double*, double*, const double*, const double*, const double*);

        double ib_height(const double, const double);

        // General settings IB
        std::string sw_ib;    ///< Namelist IB switch
        IB_type ib_type;      ///< Internal IB switch
        bool zero_ghost_tend; ///< Zero tendencies (velocity/scalar) in ghost cells

        // Interpolation settings
        int n_idw; ///< Number of points used for inverse distance weighted interpolation

        double visc_wall;  ///< Fixed viscosity at wall for LES

        std::vector<double> dem;    // User elevation map (DEM)

        // Boundary conditions
        struct Field3dBc
        {
            double bot;          ///< Value of the bottom boundary.
            Boundary_type bcbot; ///< Switch for the bottom boundary.
        };

        typedef std::map<std::string, Field3dBc*> BcMap;    ///< Map with boundary conditions per scalar
        BcMap sbc;                                          ///< Map with boundary conditions per scalar
};
#endif
