# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2019 Osyris contributors (https://github.com/nvaytet/osyris)
# @author Neil Vaytet

import numpy as np

def create_hash_table(holder):
    """
    Create a hash table for all the cells in the domain
    """

    print("Building hash table")
    holder.hash_table = dict()
    for icell in range(holder.info["ncells"]):
        igrid = int(holder.get("x")[icell]/holder.get("dx")[icell])
        jgrid = int(holder.get("y")[icell]/holder.get("dx")[icell])
        kgrid = int(holder.get("z")[icell]/holder.get("dx")[icell])
        theHash = str(igrid)+','+str(jgrid)+','+str(kgrid)+','+str(int(holder.get("level")[icell]))
        holder.hash_table[theHash] = icell

    return


def interpolate(field, points):
    """
    Interpolate data at any given point in the whole 3D domain
    """

    holder = field.parent

    try:
        hashTable = holder.hash_table
    except AttributeError:
        print("A hash table is needed to perform interpolations")
        holder.create_hash_table()

    points[:, 0] = ((points[:, 0] + holder.info["xc"]) /
                    holder.info["boxsize_scaled"])
    points[:, 1] = ((points[:, 1] + holder.info["yc"]) /
                    holder.info["boxsize_scaled"])
    points[:, 2] = ((points[:, 2] + holder.info["zc"]) /
                    holder.info["boxsize_scaled"])

    npoints = np.shape(points)[0]
    ilevl = holder.info["levelmax"]
    values = np.zeros([npoints])
    for ip in range(npoints):
        not_found = True
        loop_count = 0
        while not_found:
            l = max(min(ilevl+((-1)**loop_count) *
                        int((loop_count+1)/2), holder.info["levelmax"]), 0)
            loop_count += 1
            dxcell = 0.5**l
            igrid = int(points[ip, 0]/dxcell)
            jgrid = int(points[ip, 1]/dxcell)
            kgrid = int(points[ip, 2]/dxcell)
            theHash = str(igrid)+','+str(jgrid)+','+str(kgrid)+','+str(l)
            try:
                icell = holder.hash_table[theHash]
                ilevl = l
                not_found = False
            except KeyError:
                pass

        cube = dict()
        dmax = 0.0
        cube[theHash] = dict()
        cube[theHash]["vars"] = holder.get(field.name, only_leafs=True)[
            holder.hash_table[theHash]]
        cube[theHash]["dist"] = np.sqrt((points[ip, 0]-((igrid+0.5)*dxcell))**2 +
                                        (points[ip, 1]-((jgrid+0.5)*dxcell))**2 +
                                        (points[ip, 2]-((kgrid+0.5)*dxcell))**2)
        dmax = max(dmax, cube[theHash]["dist"])
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if i == j == k == 1:
                        pass
                    else:
                        ii = igrid-1+i
                        jj = jgrid-1+j
                        kk = kgrid-1+k
                        theHash = str(ii)+','+str(jj)+',' + \
                            str(kk)+','+str(ilevl)
                        try:
                            neighbour = holder.hash_table[theHash]
                            cube[theHash] = dict()
                            cube[theHash]["vars"] = holder.get(field.name, only_leafs=True)[
                                holder.hash_table[theHash]]
                            cube[theHash]["dist"] = np.sqrt((points[ip, 0]-((ii+0.5)*dxcell))**2 +
                                                            (points[ip, 1]-((jj+0.5)*dxcell))**2 +
                                                            (points[ip, 2]-((kk+0.5)*dxcell))**2)
                            dmax = max(dmax, cube[theHash]["dist"])
                        except KeyError:
                            theHash = str(2*ii)+','+str(2*jj) + \
                                ','+str(2*kk)+','+str(ilevl+1)
                            try:
                                neighbour = holder.hash_table[theHash]
                                for i1 in range(2):
                                    for j1 in range(2):
                                        for k1 in range(2):
                                            theHash = str(
                                                2*ii+i1)+','+str(2*jj+j1)+','+str(2*kk+k1)+','+str(ilevl+1)
                                            cube[theHash] = dict()
                                            cube[theHash]["vars"] = holder.get(field.name, only_leafs=True)[
                                                holder.hash_table[theHash]]
                                            cube[theHash]["dist"] = np.sqrt((points[ip, 0]-((2*ii+i1+0.5)*dxcell*0.5))**2 +
                                                                            (points[ip, 1]-((2*jj+j1+0.5)*dxcell*0.5))**2 +
                                                                            (points[ip, 2]-((2*kk+k1+0.5)*dxcell*0.5))**2)
                                            dmax = max(
                                                dmax, cube[theHash]["dist"])
                            except KeyError:
                                theHash = str(int(float(
                                    ii)/2.0))+','+str(int(float(jj)/2.0))+','+str(int(float(kk)/2.0))+','+str(ilevl-1)
                                try:
                                    neighbour = holder.hash_table[theHash]
                                    cube[theHash] = dict()
                                    cube[theHash]["vars"] = holder.get(field.name, only_leafs=True)[
                                        holder.hash_table[theHash]]
                                    cube[theHash]["dist"] = np.sqrt((points[ip, 0]-((int(float(ii)/2.0)+0.5)*dxcell*2.0))**2 +
                                                                    (points[ip, 1]-((int(float(jj)/2.0)+0.5)*dxcell*2.0))**2 +
                                                                    (points[ip, 2]-((int(float(kk)/2.0)+0.5)*dxcell*2.0))**2)
                                    dmax = max(dmax, cube[theHash]["dist"])
                                except KeyError:
                                    print("Neighbour not found",
                                          igrid, jgrid, kgrid, i, j, k)

        # Compute inverse distance weighting
        result = 0.0
        weights = 0.0
        for key in cube.keys():
            #w = (0.1-1.0)/dmax * cube[key]["dist"] + 1.0
            w = 1.0 / (np.exp(15.0*(cube[key]["dist"]/dmax-0.5)) + 1.0)
            weights += w
            result += w*cube[key]["vars"]

        values[ip] = result/weights

    return values
