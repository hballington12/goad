import goad_py as goad
import time

a = goad.sum_as_string(1, 2)
print(a)

goad.goad_py_add()

print("creating goad settings")

settings = goad.Settings(
    wavelength=0.532,
    beam_power_threshold=1e-1,
    beam_area_threshold_fac=1e-1,
    cutoff=0.99,
    medium_refr_index_re=1.0,
    medium_refr_index_im=0.0,
    particle_refr_index_re=1.31,
    particle_refr_index_im=0.0,
    geom_name="concave2.obj",
    max_rec=10,
    max_tir=10,
    theta_res=181,
    phi_res=181,
    euler=[0.0, 10.0, 0.0]
)

print("creating goad problem1")

problem = goad.Problem(settings)

print("solving goad problem")

from concurrent.futures import ThreadPoolExecutor

executor = ThreadPoolExecutor(max_workers=4)

start_time = time.time()

problem.py_solve()

end_time = time.time()
print(f"Solving took {end_time - start_time:.4f} seconds")

problem.py_print_stats()

# # print the first element of problem.mueller
mueller = problem.mueller

mueller_1d = problem.mueller_1d

print("mueller_1d shape:", mueller_1d[0][0])

asymmetry = problem.asymmetry

print("asymmetry parameter:", asymmetry)