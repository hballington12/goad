import goad_py as goad
import time

a = goad.sum_as_string(1, 2)
print(a)

goad.goad_py_add()

print("creating goad settings")

settings = goad.Settings(
    geom_name="hex.obj",
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