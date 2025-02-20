bl_info = {
    "name": "Custom Add-on with Rust Module",
    "author": "Your Name",
    "version": (1, 0, 0),
    "blender": (3, 0, 0),
    "location": "View3D > Sidebar > Custom Panel",
    "description": "An add-on that imports a Rust-backed Python module",
    "category": "Development",
}

import bpy
import sys
import os

# Ensure our module can be found
addon_dir = os.path.dirname(__file__)
module_path = os.path.join(addon_dir, "my_module.so")

if module_path not in sys.path:
    sys.path.append(addon_dir)

try:
    import my_module  # Import the Rust module
except ImportError as e:
    print(f"Failed to load my_module: {e}")

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
