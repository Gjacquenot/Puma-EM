import os, sys
from mesh_functions_seb import *
from mesh_functions_seb_C import *

if __name__=="__main__":
    path = './geo'
    targetNamesTmp = os.listdir(path)
    targetNames = [elem.split('.')[0] for elem in targetNamesTmp if '.geo' in elem]
    targetNames = ['cubi', 'double_elliptic_cylinder', 'NATO_f16HotSpots', 'closed_waveguide', 'strip', 'corner_reflector', 'tetrahedron', 'EMCC_ogive', 'cube', 'parabolic_reflector', 'cube_crown', 'EMCC_almond', 'EMCC_wedge-plate-cylinder_plate', 'cone-sphere', 'SRR', 'cylinder', 'cone', 'double_ellipsoid', 'corrugated', 'ellipsoid', 'parallelipiped', 'EMCC_wedge-cylinder_plate', 'open_cube', 'test', 'EMCC_cone-gap-sphere', 'cross_plates', 'dihedral_reflector', 'EMCC_double-ogive', 'strip2', 'torus', 'dihedral_with_Y', 'EMCC_cone-sphere', 'wedge', 'cylinder_antenna', 'T', 'sphere', 'parallelipiped_and_plate', 'EMCC_plate-cylinder_plate']
    #targetNames = ['double_elliptic_cylinder']
    for targetName in targetNames:
        print()
        print()
        print()
        print()
        print("####################################")
        print("##", targetName, "##################")
        print("####################################")
        f = 2.12e9
        write_geo(path, targetName, 'lc', c/f/10.0)
        write_geo(path, targetName, 'lx', .05)
        write_geo(path, targetName, 'ly', .05)
        write_geo(path, targetName, 'lz', .05)
        write_geo(path, targetName, 'w', 0.02)
        executeGmsh(path, targetName, 0)
        targetDimensions_scaling_factor = 1.0
        z_offset = 0.0
        ########################################
        ## mesh_functions_seb
        ########################################
        vertexes_coord_1, triangle_vertexes_1, triangles_physicalSurface_1 = read_mesh_GMSH_1(os.path.join(path, targetName + '.msh'), targetDimensions_scaling_factor, z_offset)

        edgeNumber_vertexes, edgeNumber_triangles, triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction = edges_computation(triangle_vertexes_1, vertexes_coord_1)

        triangles_surfaces = reorder_triangle_vertexes(triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction, triangle_vertexes_1, vertexes_coord_1)
        is_closed_surface, connected_surfaces, potential_closed_surfaces = is_surface_closed(triangles_surfaces, edgeNumber_triangles)

        RWGNumber_signedTriangles, RWGNumber_edgeVertexes, N_edges, N_RWG = RWGNumber_signedTriangles_computation(edgeNumber_triangles, edgeNumber_vertexes, triangles_surfaces, is_closed_surface, triangle_vertexes_1, vertexes_coord_1)
        RWGNumber_oppVertexes = RWGNumber_oppVertexes_computation(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, triangle_vertexes_1)
        ########################################


        ########################################
        ## mesh_functions_seb_C
        ########################################
        vertexes_coord_2, triangle_vertexes_2, triangles_physicalSurface_2 = read_mesh_GMSH_1(os.path.join(path, targetName + '.msh'), targetDimensions_scaling_factor, z_offset)

        #triangles_surfaces_C, is_closed_surface_C, RWGNumber_signedTriangles_C, RWGNumber_edgeVertexes_C, RWGNumber_oppVertexes_C = edges_computation_C_old(triangle_vertexes_2, vertexes_coord_2, './geo')
        triangles_surfaces_C, is_closed_surface_C, RWGNumber_signedTriangles_C, RWGNumber_edgeVertexes_C, RWGNumber_oppVertexes_C, triangle_vertexes_2 = edges_computation_C(triangle_vertexes_2, vertexes_coord_2, './geo')
        ########################################

        # comparison mesh_functions_seb vs mesh_functions_seb_C
        diff_triangles_surfaces = sum(abs(triangles_surfaces - triangles_surfaces_C))
        if not diff_triangles_surfaces == 0:
            print("Error in triangles_surfaces for target", targetName)
            sys.stdout.flush()
            sys.exit(1)
        diff_is_closed_surface = sum(abs(is_closed_surface - is_closed_surface_C))
        if not diff_is_closed_surface == 0:
            print("Error in is_closed_surface for target", targetName)
            sys.stdout.flush()
            sys.exit(1)
        diff_triangle_vertexes = sum(abs(triangle_vertexes_1 - triangle_vertexes_2))
        if not diff_triangle_vertexes == 0:
            print("Error in triangle_vertexes for target", targetName)
            sys.stdout.flush()
            sys.exit(1)
        diff_vertexes_coord = sum(abs(vertexes_coord_1 - vertexes_coord_2))
        if not diff_vertexes_coord == 0.0:
            print("Error in vertexes_coord for target", targetName)
            sys.stdout.flush()
            sys.exit(1)
        diff_RWGNumber_signedTriangles = sum(abs(RWGNumber_signedTriangles - RWGNumber_signedTriangles_C))
        if not diff_RWGNumber_signedTriangles == 0:
            print("Error in RWGNumber_signedTriangles for target", targetName)
            sys.stdout.flush()
            sys.exit(1)
        diff_RWGNumber_edgeVertexes = sum(abs(RWGNumber_edgeVertexes - RWGNumber_edgeVertexes_C))
        if not diff_RWGNumber_edgeVertexes == 0:
            print("Error in RWGNumber_edgeVertexes for target", targetName)
            sys.stdout.flush()
            sys.exit(1)
        diff_RWGNumber_oppVertexes = sum(abs(RWGNumber_oppVertexes - RWGNumber_oppVertexes_C))
        if not diff_RWGNumber_oppVertexes == 0:
            print("Error in RWGNumber_oppVertexes for target", targetName)
            sys.stdout.flush()
            sys.exit(1)


