import goad

euler = goad.Euler(0, 45, 0)
euler.alpha = 3
print(euler)

convention = goad.EulerConvention("xyz")
print(convention)

scheme = goad.Orientation.uniform(19)

settings = goad.Settings(geom_path="../../examples/data/hex.obj")
print(settings)
mp = goad.MultiProblem(settings)
mp.solve()
