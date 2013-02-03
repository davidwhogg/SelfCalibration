#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes / David Hogg
# August 2012

# This script is used to generate the four survey strategies investigated in
# the self-calibration paper.

from __future__ import division, print_function

import numpy as np
import parameters


def uniform_center_list(sky_limits, FoV, nx, ny):
    assert nx * FoV[0] > sky_limits[1] - sky_limits[0]
    assert ny * FoV[1] > sky_limits[3] - sky_limits[2]
    
    x_step_size = (sky_limits[1] - sky_limits[0] - FoV[0]) / (nx - 1.)
    y_step_size = (sky_limits[3] - sky_limits[2] - FoV[1]) / (ny - 1.)

    x_center_list = sky_limits[0] + 0.5 * FoV[0] + x_step_size * np.arange(nx)
    y_center_list = sky_limits[2] + 0.5 * FoV[1] + y_step_size * np.arange(ny)
    return x_center_list, y_center_list


def generate_uniform_survey(sky_limits, FoV, Ncovering,
                                                rotate=False, offset=False):

    theta = 0.
    nx = 12
    ny = 12

    x_center_list, y_center_list = uniform_center_list(sky_limits, FoV, nx, ny)
    if offset:
        x_center_list_1, y_center_list_1 =\
                    uniform_center_list(sky_limits, FoV, nx + 1, ny - 1)
        x_center_list_2, y_center_list_2 =\
                    uniform_center_list(sky_limits, FoV, nx - 1, ny + 1)

    x = np.zeros((Ncovering * nx * ny, 4))
    indx = 0
    for covering in range(Ncovering):
        xcl = x_center_list
        ycl = y_center_list
        tnx = nx
        tny = ny
        if offset:
            if covering % 3 == 1:
                xcl = x_center_list_1
                ycl = y_center_list_1
                tnx = nx + 1
                tny = ny - 1
            if covering % 3 == 2:
                xcl = x_center_list_2
                ycl = y_center_list_2
                tnx = nx - 1
                tny = ny + 1
        for yy in range(tny):
            x[indx:indx + tnx, 0] = indx + np.arange(tnx)
            x[indx:indx + tnx, 1] = xcl
            x[indx:indx + tnx, 2] = ycl[yy]
            x[indx:indx + tnx, 3] = theta
            indx += tnx
        if rotate:
            theta += 360. / Ncovering

    return x[:indx, :]


def generate_random_survey(srvy, FoV):
    N = len(srvy[:, 0])    
    # Add random offsets
    srvy[:, 1] += np.random.uniform(-0.5 * FoV[0], 0.5 * FoV[0], size=N)
    srvy[:, 2] += np.random.uniform(-0.5 * FoV[1], 0.5 * FoV[1], size=N)
    srvy[:, 3] += np.random.uniform(0., 360, size=N)
    return srvy

if __name__ == "__main__":

    dic = eval(open('parameters.py').read())
    sky_limits = [-4., 4., -4., 4.]
    number_passes = 9
    FoV = dic['FoV']

    xA = generate_uniform_survey(sky_limits, FoV, number_passes)
    np.savetxt('A.txt', xA)

    xB = generate_uniform_survey(sky_limits, FoV, number_passes, rotate=True)
    np.savetxt('B.txt', xB)

    xC = generate_uniform_survey(sky_limits, FoV, number_passes, offset=True)
    np.savetxt('C.txt', xC)

    xD = generate_random_survey(xA, FoV)
    np.savetxt('D.txt', xD)
