import goad_py as goad
import time

print("creating goad settings")
settings = goad.Settings(
    geom_name="hex.obj",
    euler=[0.0, 10.0, 0.0]
)

print("creating goad problem")
problem = goad.Problem(settings,geom=None)

print("solving goad problem")
start_time = time.time()
problem.py_solve()
end_time = time.time()
print(f"Solving took {end_time - start_time:.4f} seconds")

problem.py_print_stats()

# print the first element of problem.mueller
mueller = problem.mueller
mueller_1d = problem.mueller_1d
print("mueller_1d s11 first value:", mueller_1d[0][0])

# print the asymmetry parameter
asymmetry = problem.asymmetry
print("asymmetry parameter:", asymmetry)
