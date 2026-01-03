# --8<-- [start:basic]
from goad import Convergence, Param, Settings

settings = Settings(geom_path="../../../examples/data/hex.obj", quiet=True)
convergence = Convergence(settings)
convergence.add_target(Param.Asymmetry, 0.02)  # 2% relative error
convergence.solve()
convergence.save("my_results")  # Save results to disk for later analysis
# --8<-- [end:basic]

# --8<-- [start:multiple]
from goad import Convergence, Param, Settings

settings = Settings(
    geom_path="path/to/geometry.obj",
    quiet=True,
    particle_refr_index_im=0.001,
)
convergence = Convergence(settings)
convergence.add_target(Param.ExtCross, 0.02)  # 2% relative error
convergence.add_target(Param.Asymmetry, 0.02)
convergence.solve()
# --8<-- [end:multiple]

# --8<-- [start:albedo]
from goad import Convergence, Param, Settings

settings = Settings(
    geom_path="path/to/geometry.obj",
    quiet=True,
    particle_refr_index_im=0.001,
)
convergence = Convergence(settings)
convergence.add_target(Param.Albedo, 0.005)  # 0.5% relative error
convergence.solve()
# --8<-- [end:albedo]


# --8<-- [start:extcross]
from goad import Convergence, Param, Settings

settings = Settings(
    geom_path="path/to/geometry.obj",
    quiet=True,
    particle_refr_index_im=0.001,
)
convergence = Convergence(settings)
convergence.add_target(Param.ExtCross, 0.02)  # 2% relative error
convergence.solve()
# --8<-- [end:extcross]

# --8<-- [start:results]
import numpy as np
from goad import Convergence, Param, Settings

np.set_printoptions(threshold=4, precision=3)

settings = Settings(
    geom_path="./data/hex.obj", quiet=True, particle_refr_index_im=0.001
)
convergence = Convergence(settings)
convergence.add_target(Param.Asymmetry, 0.01)  # 1% relative error
convergence.solve()

# Access mean results and standard errors
mean = convergence.mean
sem = convergence.sem

# Get the main integrated parameters
print(f"Asymmetry: {mean.asymmetry:.4f} +/- {sem.asymmetry:.4f}")
print(f"Scattering Cross Section: {mean.scat_cross:.4f}")
print(f"Extinction Cross Section: {mean.ext_cross:.4f}")
print(f"Absorption Cross Section: {mean.ext_cross - mean.scat_cross:.4f}")
print(f"Single Scattering Albedo: {mean.albedo:.4f}")

# Access the Mueller matrix
print(f"Theta bins:\n{mean.bins_1d[:]}")
print(f"[Theta, Phi] bins:\n{mean.bins[:]}")
print(f"Mueller matrix S11: {mean.mueller_1d[:, 0]}")
print(f"Mueller matrix S12: {mean.mueller_1d[:, 1]}")
# --8<-- [end:results]
