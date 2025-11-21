# --8<-- [start:multiproblem]
from goad import MultiProblem, Settings

# Setup and run a multi-orientation problem with default settings
mp = MultiProblem(Settings(geom_path="path/to/geometry.obj"))
mp.solve()
# --8<-- [end:multiproblem]
