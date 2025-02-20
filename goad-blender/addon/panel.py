import bpy

class MYADDON_PT_Panel(bpy.types.Panel):
    bl_label = "My Addon Panel"
    bl_idname = "MYADDON_PT_panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Custom"

    def draw(self, context):
        layout = self.layout
        layout.operator("myaddon.call_rust_function", text="Run Rust Code")

def register():
    bpy.utils.register_class(MYADDON_PT_Panel)

def unregister():
    bpy.utils.unregister_class(MYADDON_PT_Panel)
