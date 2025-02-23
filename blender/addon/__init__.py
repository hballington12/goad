bl_info = {
    "name": "goad-blender",
    "author": "Harry",
    "version": (1, 0, 5),
    "blender": (3, 0, 0),
    "location": "View3D > Sidebar > Custom Panel",
    "description": "An add-on that imports a Rust-backed Python module",
    "category": "Development",
}

import bpy
import os
import sys
import importlib

# Determine the path to the Rust module
addon_dir = os.path.dirname(__file__)
bin_dir = os.path.join(addon_dir, "goad_py")
so_file = "goad_py.cpython-311-x86_64-linux-gnu.so"
module_path = os.path.join(bin_dir, so_file)

# Load the compiled Rust module
if os.path.exists(module_path):
    spec = importlib.util.spec_from_file_location("goad_py", module_path)
    goad_py = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(goad_py)
    sys.modules["goad_py"] = goad_py
else:
    print(f"Error: Rust module not found at {module_path}")


# Import other add-on parts
from . import panel, operators

# Register functions
def register():
    panel.register()
    operators.register()


def unregister():
    panel.unregister()
    operators.unregister()

if __name__ == "__main__":
    register()
