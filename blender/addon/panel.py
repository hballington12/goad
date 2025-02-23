import bpy

class SCATTERING_PT_Panel(bpy.types.Panel):
    """Light Scattering UI Panel"""
    bl_label = "GOAD"
    bl_idname = "SCATTERING_PT_Panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "GOAD"

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

def register():
    bpy.utils.register_class(SCATTERING_PT_Panel)
    bpy.types.Scene.scattering_intensity = bpy.props.FloatProperty(name="Wavelength3", default=1.0, min=0.1, max=10.0)
    bpy.types.Scene.scattering_angle = bpy.props.FloatProperty(name="Something Else", default=45.0, min=0.0, max=180.0)


def unregister():
    bpy.utils.unregister_class(SCATTERING_PT_Panel)
    del bpy.types.Scene.scattering_intensity
    del bpy.types.Scene.scattering_angle
