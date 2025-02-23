import bpy
import goad_py


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
    bpy.utils.register_class(OBJECT_OT_ComputeScattering)


def unregister():
    bpy.utils.unregister_class(OBJECT_OT_ComputeScattering)
