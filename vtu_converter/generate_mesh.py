import mesh_generator

mesh_generator = mesh_generator.mesh_generator()
mesh_generator.set_exfile("heart1.exfile")
mesh_generator.set_outfile("heart111.vtu")
mesh_generator.set_outfile_fibers("fibers220surface.vtu")
mesh_generator.set_refinement_fibers(2,2,0)
mesh_generator.set_refinement(4,4,4)
mesh_generator.set_fibers_surface(True)


mesh_generator.generate_mesh(True)
