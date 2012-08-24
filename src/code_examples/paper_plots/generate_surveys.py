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


def generate_random_survey(sky_limits, FoV, N):
    nx = np.ceil((sky_limits[1] - sky_limits[0]) / FoV[0]).astype(int)
    ny = np.ceil((sky_limits[3] - sky_limits[2]) / FoV[1]).astype(int)
    x_box_size = (sky_limits[1] - sky_limits[0]) / nx
    y_box_size = (sky_limits[3] - sky_limits[2]) / ny
    ntheta = np.floor(N / (nx * ny)).astype(int)
    theta_box_size = 360. / ntheta
    x = np.zeros((N, 4))
    x[:, 0] = range(N)
    # Lay down left edges
    x[:, 1] = sky_limits[0] + x_box_size * np.mod(np.arange(N), nx)
    # Lay down bottom edges
    x[:, 2] = sky_limits[2] + y_box_size \
                        * np.mod(np.floor(np.arange(N) / nx), ny)
    # Lay down orientation edges
    x[:, 3] = 0.0 + theta_box_size \
                        * np.mod(np.floor(np.arange(N) / (nx * ny)), ntheta)
    # Add random offsets
    x[:, 1] += np.random.uniform(0., x_box_size, size=N)
    x[:, 2] += np.random.uniform(0., y_box_size, size=N)
    x[:, 3] += np.random.uniform(0., theta_box_size, size=N)
    return x

if __name__ == "__main__":

    dic = eval(open('parameters.py').read())
    sky_limits = [-4., 4., -4., 4.]
    number_passes = 12
    FoV = dic['FoV']

    xA = generate_uniform_survey(sky_limits, FoV, number_passes)
    np.savetxt('A.txt', xA)

    xB = generate_uniform_survey(sky_limits, FoV, number_passes, rotate=True)
    np.savetxt('B.txt', xB)

    xC = generate_uniform_survey(sky_limits, FoV, number_passes, offset=True)
    np.savetxt('C.txt', xC)

    xD = generate_random_survey(sky_limits, FoV, len(xA))
    np.savetxt('D.txt', xD)
