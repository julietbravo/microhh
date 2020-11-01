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

#include <GL/glew.h>
#include <GL/freeglut.h>

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

namespace
{
    void display_func()
    {
    }
}

template<typename TF>
Visualization<TF>::Visualization(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_visualisation = inputin.get_item<bool>("visualisation", "swvisualisation", "", false);
}

template <typename TF>
void Visualization<TF>::create()
{
    if (!sw_visualisation)
        return;

    auto& gd = grid.get_grid_data();
    width = gd.itot;
    height = gd.ktot;

    char* argv[1];
    int argc=1;
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(width, height);
    glutCreateWindow("MicroHH realtime visualization :-D");
    glewInit();

    // Setup 2D orthographic projection
    gluOrtho2D(0, width, height, 0);

    // DUMMY -> not used as we don't use the
    // glutMainLoop(), but something still needs to
    // be set as the DisplayFunc...
    glutDisplayFunc(display_func);

    // Init pixel buffer
    glGenBuffers(1, &pbo);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo);
    glBufferData(GL_PIXEL_UNPACK_BUFFER, 4*width*height*sizeof(GLubyte), 0, GL_STREAM_DRAW);
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    cudaGraphicsGLRegisterBuffer(
            &cuda_pbo_resource, pbo,
            cudaGraphicsMapFlagsWriteDiscard);

    // Init device (colormaps)
    prepare_device();
}

template <typename TF>
Visualization<TF>::~Visualization()
{
}

template class Visualization<double>;
template class Visualization<float>;
