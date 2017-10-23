#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  posPros.py
#
#  Copyright 2017 Esbel Tomas Valero Orellana <evalero@lisbeth.sabessa>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D




def main(args):
    distNormal = np.loadtxt('norm.dat', usecols=[0])
    denProb = np.loadtxt('norm.dat', usecols=[1])

    plt.figure(1)
    plt.plot(distNormal,denProb,'*r')

    fig = plt.figure(2)
    ax = Axes3D(fig)
    xs = np.loadtxt('grid1.dat', usecols=[0])
    ys = np.loadtxt('grid1.dat', usecols=[1])
    zs = np.loadtxt('grid1.dat', usecols=[2])
    ax.scatter(xs, ys, zs, c='r', marker='o')

    fig = plt.figure(3)
    ax = Axes3D(fig)
    xs = np.loadtxt('grid2.dat', usecols=[0])
    ys = np.loadtxt('grid2.dat', usecols=[1])
    zs = np.loadtxt('grid2.dat', usecols=[2])
    ax.scatter(xs, ys, zs, c='r', marker='o')
    
    plt.show()

    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
