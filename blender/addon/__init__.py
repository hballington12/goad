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
import subprocess
import importlib
from io import BytesIO
import site
import threading
import queue

# Function to install required packages
def install_packages():
    # Add user site-packages to Python path
    user_site = site.getusersitepackages()
    if user_site not in sys.path:
        sys.path.append(user_site)
    
    packages = ["matplotlib", "numpy", "PIL", "tornado"]
    python_executable = sys.executable
    for package in packages:
        try:
            importlib.import_module(package)
            module = importlib.import_module(package)
            print(f"Package {package} location: {module.__file__}")
        except ImportError:
            subprocess.check_call([python_executable, "-m", "pip", "install", package])
            module = importlib.import_module(package)
            print(f"Installed {package} at: {module.__file__}")

# Install required packages
install_packages()

import matplotlib
matplotlib.use('WebAgg') # Use the WebAgg backend for plotting

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# Determine the path to the Rust module
addon_dir = os.path.dirname(__file__)
bin_dir = os.path.join(addon_dir, "goad_py")
so_file = "goad_py.cpython-311-x86_64-linux-gnu.so"
module_path = os.path.join(bin_dir, so_file)

# Add binary directory to Python path
if bin_dir not in sys.path:
    sys.path.append(bin_dir)

# Load the compiled Rust module
try:
    if "goad_py" in sys.modules:
        print("Reloading Rust module")
        importlib.reload(sys.modules["goad_py"])
    import goad_py
except ImportError as e:
    print(f"Error loading Rust module: {e}")
    print(f"Expected module path: {module_path}")
    if not os.path.exists(module_path):
        print("Module file does not exist!")


def update_plot(data):
    """
    Update the plot with new data.

    Parameters:
    data (list or numpy array): The data to plot.
    """
    print("Updating plot")
    plt.figure(1)  # Use figure with label 1
    plt.plot(data)  # Plot the new data
    # plt.show()  # Show the plot in a new window

class SCATTERING_PT_Panel(bpy.types.Panel):
    """Light Scattering UI Panel"""
    bl_label = "GOAD"
    bl_idname = "SCATTERING_PT_Panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "GOAD"

    plt.ion()

    def draw(self, context):
        layout = self.layout
        obj = context.object
        
        # Object selection
        if obj:
            layout.label(text=f"Selected: {obj.name}")
        else:
            layout.label(text="No object selected")

        # Light Scattering Controls
        layout.prop(context.scene, "scattering_intensity")
        layout.prop(context.scene, "scattering_angle")

        # Compute Scattering Button
        layout.operator("object.compute_scattering")

class OBJECT_OT_ComputeScattering(bpy.types.Operator):
    """Compute light scattering"""
    bl_idname = "object.compute_scattering"
    bl_label = "Compute Scattering"

    def execute(self, context):
        obj = context.object
        if obj is None:
            self.report({'WARNING'}, "No object selected")
            return {'CANCELLED'}

        # Extract mesh and light properties
        mesh_data = extract_mesh(obj)
        light_params = {
            "wavelength": context.scene.scattering_intensity,
            "angle": context.scene.scattering_angle
        }

        shape_id = 0  # note that we need care to deal with the containment graph
        refr_re = 1.31
        refr_im = 0.0

        vertices = []
        vertex_indices = []
        vertex_counter = 0
        faces = mesh_data["faces"]

        # Process only vertices that are used in faces
        for face in faces:
            face_indices = []
            face_vertices = [mesh_data["vertices"][i] for i in face]
            
            for v in face_vertices:
                vertices.append((v.x, v.y, v.z))
                face_indices.append(vertex_counter)
                vertex_counter += 1
                
            vertex_indices.append(tuple(face_indices))

        # Print vertices in each face
        for i, face in enumerate(vertex_indices):
            face_vertices = [vertices[i] for i in face]
            print(f"Face vertices: {face_vertices}")

        faces = vertex_indices

        shape = goad_py.Shape(vertices, faces, shape_id, refr_re, refr_im)

        shapes = [shape]

        print("creating geometry object")

        geom = goad_py.Geom(shapes)

        print("creating goad settings")

        settings = goad_py.Settings(
            wavelength=0.532,
            beam_power_threshold=1e-1,
            beam_area_threshold_fac=1e-1,
            total_power_cutoff=0.99,
            medium_refr_index_re=1.0,
            medium_refr_index_im=0.0,
            particle_refr_index_re=1.31,
            particle_refr_index_im=0.0,
            geom_name="test",
            max_rec=10,
            max_tir=10,
            theta_res=100,
            phi_res=100
        )

        print("creating goad problem")

        problem = goad_py.Problem(geom, settings)

        print("solving goad problem")

        problem.py_solve()

        problem.py_print_stats()

        # Example data to plot
        data = np.random.rand(100) * 10  # Generate 100 random numbers between 0 and 10
        # update_plot(data)
        update_plot(data)  # Call the function in a new thread

        return {'FINISHED'}

# Extract mesh vertices and faces
def extract_mesh(obj):
    # Get the world matrix which includes all transformations
    world_matrix = obj.matrix_world
    # Transform vertices using the world matrix
    vertices = [(world_matrix @ v.co) for v in obj.data.vertices]
    faces = [f.vertices for f in obj.data.polygons]
    return {"vertices": vertices, "faces": faces}

def register():
    # Reload goad_py module if it exists
    try:
        if "goad_py" in sys.modules:
            print("Reloading Rust module during registration")
            importlib.reload(sys.modules["goad_py"])
        import goad_py
    except ImportError as e:
        print(f"Error loading Rust module during registration: {e}")
        print(f"Expected module path: {module_path}")
        if not os.path.exists(module_path):
            print("Module file does not exist!")
            return

    bpy.utils.register_class(SCATTERING_PT_Panel)
    bpy.utils.register_class(OBJECT_OT_ComputeScattering)
    bpy.types.Scene.scattering_intensity = bpy.props.FloatProperty(name="Wavelength", default=1.0, min=0.1, max=10.0)
    bpy.types.Scene.scattering_angle = bpy.props.FloatProperty(name="Scattering Angle", default=45.0, min=0.0, max=180.0)

def unregister():
    del bpy.types.Scene.scattering_angle
    del bpy.types.Scene.scattering_intensity
    bpy.utils.unregister_class(OBJECT_OT_ComputeScattering)
    bpy.utils.unregister_class(SCATTERING_PT_Panel)

if __name__ == "__main__":
    register()
