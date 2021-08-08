# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2021 Osyris contributors (https://github.com/nvaytet/osyris)

# flake8: noqa

import sys
import struct
import numpy as np


def to_vtk(holder, fname="osyris_data.vtu", dx=None, dy=None, dz=None):
    """
    Write RAMSES data to VTK file
    """

    try:
        from scipy.spatial import Delaunay
    except ImportError:
        print(
            "Scipy Delaunay library not found. This is needed for VTK output. Exiting.")

    # Print status
    if not fname.endswith(".vtu"):
        fname += ".vtu"
    print("Writing data to VTK file: " + fname)

    if dx is not None:
        if holder.info["ndim"] == 1:
            cube = np.where(np.logical_and(holder.get("x")>=-0.5*dx, \
                                           holder.get("x")<= 0.5*dx))
        elif holder.info["ndim"] == 2:
            cube = np.where(np.logical_and(holder.get("x")>=-0.5*dx, \
                            np.logical_and(holder.get("x")<= 0.5*dx, \
                            np.logical_and(holder.get("y")>=-0.5*dy, \
                                           holder.get("y")<= 0.5*dy))))
        elif holder.info["ndim"] == 3:
            cube = np.where(np.logical_and(holder.get("x")>=-0.5*dx, \
                            np.logical_and(holder.get("x")<= 0.5*dx, \
                            np.logical_and(holder.get("y")>=-0.5*dy, \
                            np.logical_and(holder.get("y")<= 0.5*dy, \
                            np.logical_and(holder.get("z")>=-0.5*dz, \
                                           holder.get("z")<= 0.5*dz))))))
        points = np.array(
            [holder.get("x")[cube],
             holder.get("y")[cube],
             holder.get("z")[cube]]).T
    else:
        # Coordinates ot RAMSES cell centers
        points = np.array([holder.get("x"), holder.get("y"), holder.get("z")]).T

    # Compute Delaunay tetrahedralization from cell nodes
    # Note that this step can take a lot of time!
    ncells = points.shape[0]
    print("Computing Delaunay mesh with %i points." % ncells)
    print("This may take some time...")
    tri = Delaunay(points, qhull_options="QJ Qx Qs Qv")
    ntetra = np.shape(tri.simplices)[0]
    nverts = ntetra * 4
    print("Delaunay mesh with %i tetrahedra complete." % ntetra)

    # Create list of variables by grouping x,y,z components together
    key_list = holder.get_var_list()
    n_components = []
    varlist = []
    for key in key_list:
        thisVar = getattr(holder, key)
        if (not thisVar.vector_component) and (thisVar.group !=
                                               "amr") and (thisVar.group != "part"):
            if thisVar.kind == "vector":
                n_components.append(3)
            elif thisVar.kind == "scalar":
                n_components.append(1)
            else:
                print("Unknown data type: " + thisVar.kind)
                return
            varlist.append(key)

    nvars = len(n_components)

    # Compute byte sizes
    nbytes_xyz = 3 * ncells * 8
    nbytes_cellc = nverts * 4
    nbytes_cello = ntetra * 4
    nbytes_cellt = ntetra * 4
    nbytes_vars = np.zeros([nvars], dtype=np.int32)
    for i in range(nvars):
        nbytes_vars[i] = n_components[i] * ncells * 8

    # Compute byte offsets
    offsets = np.zeros([nvars + 4], dtype=np.int64)
    offsets[0] = 0  # xyz coordinates
    offsets[1] = offsets[0] + 4 + nbytes_xyz  # cell connectivity
    offsets[2] = offsets[1] + 4 + nbytes_cellc  # cell offsets
    offsets[3] = offsets[2] + 4 + nbytes_cello  # cell types
    offsets[4] = offsets[3] + 4 + nbytes_cellt  # first hydro variable
    for i in range(nvars - 1):
        offsets[i + 5] = offsets[i + 4] + 4 + nbytes_vars[i]

    # Open file for binary output
    with open(fname, "wb") as f:
        # Write VTK file header
        f.write(string_to_binary('<?xml version=\"1.0\"?>'))
        f.write(
            string_to_binary(
                '<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">'
            ))
        f.write(string_to_binary('   <UnstructuredGrid>'))
        f.write(
            string_to_binary('   <Piece NumberOfPoints=\"%i\" NumberOfCells=\"%i\">' %
                             (ncells, ntetra)))
        f.write(string_to_binary('      <Points>'))
        f.write(
            string_to_binary(
                '         <DataArray type=\"Float64\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%i\" />'
                % offsets[0]))
        f.write(string_to_binary('      </Points>'))
        f.write(string_to_binary('      <Cells>'))
        f.write(
            string_to_binary(
                '         <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"%i\" />'
                % offsets[1]))
        f.write(
            string_to_binary(
                '         <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"%i\" />'
                % offsets[2]))
        f.write(
            string_to_binary(
                '         <DataArray type=\"Int32\" Name=\"types\" format=\"appended\" offset=\"%i\" />'
                % offsets[3]))
        f.write(string_to_binary('      </Cells>'))
        f.write(string_to_binary('      <PointData>'))
        for i in range(nvars):
            f.write(
                string_to_binary(
                    '         <DataArray type=\"Float64\" Name=\"' + varlist[i] +
                    '\" NumberOfComponents=\"%i\" format=\"appended\" offset=\"%i\" />'
                    % (n_components[i], offsets[i + 4])))
        f.write(string_to_binary('      </PointData>'))
        f.write(string_to_binary('   </Piece>'))
        f.write(string_to_binary('   </UnstructuredGrid>'))
        f.write(string_to_binary('   <AppendedData encoding=\"raw\">'))
        f.write(string_to_binary('_', newline=False))

        # Now write data in binary. Every data field is preceded by its byte size.
        # x,y,z coordinates of the points
        f.write(struct.pack('<i', *[nbytes_xyz]))
        f.write(struct.pack('<%id' % (ncells * 3), *np.ravel(points)))

        # Cell connectivity
        f.write(struct.pack('<i', *[nbytes_cellc]))
        f.write(struct.pack('<%ii' % nverts, *np.ravel(tri.simplices)))

        # Cell offsets
        f.write(struct.pack('<i', *[nbytes_cello]))
        f.write(struct.pack('<%ii' % ntetra, *range(4, ntetra * 4 + 1, 4)))

        # Cell types: number 10 is tetrahedron in VTK file format
        f.write(struct.pack('<i', *[nbytes_cellt]))
        f.write(struct.pack('<%ii' % ntetra, *np.full(ntetra, 10, dtype=np.int32)))

        # Cell variables
        for i in range(nvars):
            if n_components[i] == 3:
                if dx is not None:
                    celldata = np.ravel(
                        np.array([
                            holder.get(varlist[i] + "_x")[cube],
                            holder.get(varlist[i] + "_y")[cube],
                            holder.get(varlist[i] + "_z")[cube]
                        ]).T)
                else:
                    celldata = np.ravel(
                        np.array([
                            holder.get(varlist[i] + "_x"),
                            holder.get(varlist[i] + "_y"),
                            holder.get(varlist[i] + "_z")
                        ]).T)
            else:
                if dx is not None:
                    celldata = holder.get(varlist[i])[cube]
                else:
                    celldata = holder.get(varlist[i])
            f.write(struct.pack('<i', *[nbytes_vars[i]]))
            f.write(struct.pack('<%id' % (ncells * n_components[i]), *celldata))

        # Close file
        f.write(string_to_binary('   </AppendedData>'))
        f.write(string_to_binary('</VTKFile>'))
        # f.close()

    # File size
    fsize_raw = offsets[nvars + 3] + nbytes_vars[nvars - 1]
    if fsize_raw > 1.0e9:
        fsize = float(fsize_raw) / 1.0e9
        funit = "GB"
    elif fsize_raw > 1.0e6:
        fsize = float(fsize_raw) / 1.0e6
        funit = "MB"
    elif fsize_raw > 1000:
        fsize = float(fsize_raw) / 1.0e3
        funit = "KB"
    else:
        fsize = float(fsize_raw)
        funit = "B"

    print("File " + fname + (" of size %.1f" % fsize) + funit + " succesfully written.")

    return


def string_to_binary(s, newline=True):
    if sys.version_info < (3, ):
        return s + ("\n" if newline else "")
    else:
        return str.encode(s) + (b"\n" if newline else b"")
