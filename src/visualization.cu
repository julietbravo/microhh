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
    __device__
    unsigned char clip(int n)
    {
      return n > 255 ? 255 : (n < 0 ? 0 : n);
    }

    template<typename TF> __global__
    void set_values(
            uchar4* image,
            const TF* const __restrict__ fld,
            const int width, const int height,
            const TF vmin, const TF vmax,
            const int j_index,
            const int istart, const int iend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        // Column/row index
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int k = blockIdx.y*blockDim.y + threadIdx.y;

        if ( (i >= width) || (k >= height) )
            return;

        const int ij  = i+(height-k-1)*width;
        const int ijk = (i+istart) + j_index*icells + (k+kstart)*ijcells;

        const TF rel_val = (fld[ijk] - vmin) / (vmax - vmin);
        const unsigned char rel_val_i = clip(rel_val * 255);

        image[ij].x = rel_val_i;
        image[ij].y = rel_val_i;
        image[ij].z = rel_val_i;
        image[ij].w = 255;
    }
}

template <typename TF>
void Visualization<TF>::exec(Timeloop<TF>& timeloop)
{
    if (timeloop.in_substep())
        return;

    auto& gd = grid.get_grid_data();

    // This only works with freeglut, but allows us
    // to use OpenGL without the glutMainLoop()
    glutMainLoopEvent();

    // Map graphics resources so that CUDA can access it
    uchar4 *d_out = 0;
    cudaGraphicsMapResources(1, &cuda_pbo_resource, 0);
    cudaGraphicsResourceGetMappedPointer(
            (void **)&d_out, NULL, cuda_pbo_resource);

    // Launch kernel
    const int blocki = 32;
    const int blockj = 32;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.kmax/blockj + (gd.kmax%blockj > 0);

    dim3 grid_gpu(gridi, gridj, 1);
    dim3 block_gpu(blocki, blockj, 1);

    const TF vmin = 300;
    const TF vmax = 303;

    set_values<TF><<<grid_gpu, block_gpu>>>(
            d_out, fields.ap.at("th")->fld_g,
            width, height, vmin, vmax, gd.jstart,
            gd.istart, gd.iend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    // Unmap graphics resources
    cudaGraphicsUnmapResources(1, &cuda_pbo_resource, 0);

    // Draw texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex2f(0, 0);
    glTexCoord2f(0.0f, 1.0f); glVertex2f(0, height);
    glTexCoord2f(1.0f, 1.0f); glVertex2f(width, height);
    glTexCoord2f(1.0f, 0.0f); glVertex2f(width, 0);
    glEnd();
    glDisable(GL_TEXTURE_2D);

    // Swap buffer to show updated image
    glutSwapBuffers();
}

template class Visualization<double>;
template class Visualization<float>;
