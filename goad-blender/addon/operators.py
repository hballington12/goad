import bpy
import goad_py  # Your Rust module

class MYADDON_OT_CallRustFunction(bpy.types.Operator):
    bl_idname = "myaddon.call_rust_function"
    bl_label = "Call Rust Function"
    bl_description = "Calls a function from the Rust-backed Python module"

    def execute(self, context):
        try:
            result = my_module.compute_scattering()  # Example Rust function
            self.report({'INFO'}, f"Result: {result}")
        except Exception as e:
            self.report({'ERROR'}, f"Error: {e}")

        return {'FINISHED'}

def register():
    bpy.utils.register_class(MYADDON_OT_CallRustFunction)

def unregister():
    bpy.utils.unregister_class(MYADDON_OT_CallRustFunction)
