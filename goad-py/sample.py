import goad_py as goad
import time

a = goad.sum_as_string(1, 2)
print(a)

goad.goad_py_add()

print("creating shape object")
vertices = [
    (6.025435, 6.025435, -6.025435),
    (6.025435, 6.025435, 6.025435),
    (6.025435, -6.025435, -6.025435),
    (6.025435, -6.025435, 6.025435),
    (-6.025435, 6.025435, -6.025435),
    (-6.025435, 6.025435, 6.025435),
    (-6.025435, -6.025435, -6.025435),
    (-6.025435, -6.025435, 6.025435)
]
faces = [
    (0, 1, 3, 2),
    (2, 3, 7, 6),
    (6, 7, 5, 4),
    (4, 5, 1, 0),
    (2, 6, 4, 0),
    (7, 3, 1, 5)
]
shape_id = 0 # note that we need care to deal with the containment graph
# the shape_id for each face is currently bodged in impl for pymethods on Shape struct
refr_re = 1.31
refr_im = 0.0

shape = goad.Shape(vertices, faces, shape_id, refr_re, refr_im)

shapes = [shape]

print("creating geometry object")

geom = goad.Geom(shapes)

print("creating goad settings")

settings = goad.Settings(
    wavelength=0.532,
    beam_power_threshold=1e-1,
    beam_area_threshold_fac=1e-1,
    total_power_cutoff=0.99,
    medium_refr_index_re=1.0,
    medium_refr_index_im=0.0,
    particle_refr_index_re=1.31,
    particle_refr_index_im=0.0,
    geom_name="hex.obj",
    max_rec=10,
    max_tir=10,
    theta_res=361,
    phi_res=181,
    euler=[0.0, 20.0, 0.0]
)

print("creating goad problem")

problem = goad.Problem(settings)

print("solving goad problem")

from concurrent.futures import ThreadPoolExecutor

executor = ThreadPoolExecutor(max_workers=4)

start_time = time.time()


end_time = time.time()
print(f"Solving took {end_time - start_time:.4f} seconds")

problem.py_print_stats()

# # print the first element of problem.mueller
mueller = problem.mueller

mueller_1d = problem.mueller_1d
