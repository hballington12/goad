import goad_py as goad
import time
from scipy.stats import uniform, norm
import monaco as mc
import matplotlib.pyplot as plt
import numpy as np

def goad_run(alpha, beta, gamma):
    """
    Run the GOAD simulation with the given parameters.
    """

    settings = goad.Settings(
        euler=[alpha, beta, gamma],
    )

    problem = goad.Problem(settings)
    problem.py_solve()
    # problem.py_print_stats()

    asymmetry = problem.asymmetry
    # print("asymmetry parameter:", asymmetry)

    return asymmetry

def goad_preprocess(case):
    """
    Pre-process the input parameters for the GOAD simulation.
    """
    # Map the euler values
    alpha = case.invals["alpha"].val * 360
    beta = np.arccos(1 - 2 * case.invals["beta"].val) * 180 / np.pi
    gamma = case.invals["gamma"].val * 360

    return (alpha, beta, gamma)

def goad_postprocess(case,asymmetry):
    """
    Post-process the results of the GOAD simulation.
    """
    case.addOutVal("asymmetry", asymmetry)


ndraws = 2**10
seed = 234234231
fcns = {
    "preprocess" : goad_preprocess,
    "run" : goad_run,
    "postprocess" : goad_postprocess,
}

def goad_sim():
    # Define simulation
    sim = mc.Sim(name="goad", ndraws=ndraws, fcns=fcns, firstcaseismedian=True,
                 seed=seed, singlethreaded=True, verbose=True, debug=True,
                 savecasedata=False, savesimdata=False)

    # Define input variables
    sim.addInVar(name="alpha", dist=uniform, distkwargs={"loc": 0, "scale": 1})
    sim.addInVar(name="beta", dist=uniform, distkwargs={"loc": 0, "scale": 1})
    sim.addInVar(name="gamma", dist=uniform, distkwargs={"loc": 0, "scale": 1})

    # Run sim
    sim.runSim()
    # Get results
    asymmetry_values = sim.outvars["asymmetry"].nums
    print("asymmetry values:", asymmetry_values)

    sim.plot()

    # fig, axs = mc.multi_plot([sim.invars["gamma"], sim.outvars["asymmetry"]],
    #                          title="asymmetry vs gamma", cov_plot=False
    #                          )
    # fig.set_size_inches(8.8, 6.6)
    plt.savefig('launch_angle_vs_landing.png', dpi=100)
    plt.show(block=False) 

    # Optional: Plot histogram of beta values to verify distribution
    # beta_values = sim.invars["beta"].nums
    # plt.hist(beta_values, bins=50, density=True)
    # plt.title("Beta Distribution")
    # plt.xlabel("Beta (degrees)")
    # plt.ylabel("Density")
    plt.show()

    return sim


if __name__ == '__main__':
    sim = goad_sim()