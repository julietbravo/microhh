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

#ifndef VISUALIZATION_H
#define VISUALIZATION_H

#include <GL/glew.h>
#include <GL/freeglut.h>

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timeloop;


template<typename TF>
class Visualization
{
    public:
        Visualization(Master&, Grid<TF>&, Fields<TF>&, Input&); // Constructor of the decay class.
        ~Visualization();                                       // Destructor of the decay class.

        void init();
        void create();
        void exec(Timeloop<TF>&);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        // Size of buffer
        int width;
        int height;

        // Texture and pixel objects
        GLuint pbo = 0;     // OpenGL pixel buffer object
        GLuint tex = 0;     // OpenGL texture object
        struct cudaGraphicsResource *cuda_pbo_resource;
};
#endif
