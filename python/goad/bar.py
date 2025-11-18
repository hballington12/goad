import goad

euler = goad.Euler(0, 45, 0)
euler.alpha = 3
print(euler)

convention = goad.EulerConvention("xyz")
print(convention)

scheme = goad.Orientation.uniform(19)
binning = goad.BinningScheme.simple(num_theta=100, num_phi=100)
binning = goad.BinningScheme.interval(
    thetas=[0, 1, 180], theta_spacings=[0.1, 1], phis=[0, 360], phi_spacings=[2]
)
print(binning.thetas())
print(binning.phis())

orientation = goad.Orientation.discrete([euler, euler])
orientation = goad.Orientation.uniform(1)
settings = goad.Settings(
    geom_path="../../examples/data/hex.obj", orientation=orientation
)
print(settings)
mp = goad.MultiProblem(settings)
mp.solve()
print(mp.results.asymmetry)

mueller_2d = mp.results.mueller
mueller_1d = mp.results.mueller_1d

if mueller_1d is not None:
    print(mueller_1d[:, 0])
