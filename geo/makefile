
OBJECTS_msh = parallelipiped.msh cube.msh strip.msh strip2.msh sphere.msh cubi.msh EMCC_ogive.msh EMCC_double-ogive.msh EMCC_cone-sphere.msh EMCC_cone-gap-sphere.msh EMCC_wedge-cylinder_plate.msh EMCC_plate-cylinder_plate.msh EMCC_wedge-plate-cylinder_plate.msh 

meshes: $(OBJECTS_msh)
	ls -al $(OBJECTS_msh)

parallelipiped.msh: parallelipiped.geo
	gmsh -2 parallelipiped.geo
cube.msh: cube.geo
	gmsh -2 cube.geo
strip.msh: strip.geo
	gmsh -2 strip.geo
strip2.msh: strip2.geo
	gmsh -2 strip2.geo
sphere.msh: sphere.geo
	gmsh -2 sphere.geo
cubi.msh: cubi.geo
	gmsh -2 cubi.geo
EMCC_ogive.msh: EMCC_ogive.geo
	gmsh -2 EMCC_ogive.geo
EMCC_double-ogive.msh: EMCC_double-ogive.geo
	gmsh -2 EMCC_double-ogive.geo
EMCC_cone-sphere.msh: EMCC_cone-sphere.geo
	gmsh -2 EMCC_cone-sphere.geo
EMCC_cone-gap-sphere.msh: EMCC_cone-gap-sphere.geo
	gmsh -2 EMCC_cone-gap-sphere.geo
EMCC_wedge-cylinder_plate.msh: EMCC_wedge-cylinder_plate.geo
	gmsh -2 EMCC_wedge-cylinder_plate.geo
EMCC_plate-cylinder_plate.msh: EMCC_plate-cylinder_plate.geo
	gmsh -2 EMCC_plate-cylinder_plate.geo
EMCC_wedge-plate-cylinder_plate.msh: EMCC_wedge-plate-cylinder_plate.geo
	gmsh -2 EMCC_wedge-plate-cylinder_plate.geo
clean: 
	rm -f *.msh* *~ *.cubes* *.pos* *.unv* *.txt*
