import goad
from goad.goad import Settings

euler = goad.Euler(0, 45, 0)
euler.alpha = 3
print(euler)

convention = goad.EulerConvention("xyz")
print(convention)

settings = Settings()
