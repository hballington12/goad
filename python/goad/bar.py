import goad

euler = goad.Euler(0, 45, 0)
euler.alpha = 3
print(euler)

convention = goad.EulerConvention("xyz")
print(convention)

scheme = goad.Orientation.uniform(19)


orientation = goad.Orientation.discrete([euler, euler])
orientation = goad.Orientation.uniform(1)
settings = goad.Settings(
    geom_path="../../examples/data/hex.obj", orientation=orientation
)
print(settings)
mp = goad.MultiProblem(settings)
mp.solve()
print(mp.results.asymmetry)
