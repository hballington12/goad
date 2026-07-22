"""Triangulate aggregate-2p9mm.obj to satisfy GOAD's planar-face requirement."""
import bpy
import bmesh

from pathlib import Path

output_dir = Path(__file__).parent
src = output_dir / "aggregate-2p9mm.obj"
dst = output_dir / "aggregate-2p9mm_tri.obj"

bpy.ops.wm.read_factory_settings(use_empty=True)
bpy.ops.wm.obj_import(filepath=str(src))
obj = bpy.context.selected_objects[0]

bpy.context.view_layer.objects.active = obj
bpy.ops.object.mode_set(mode="EDIT")
bm = bmesh.from_edit_mesh(obj.data)
bmesh.ops.triangulate(bm, faces=bm.faces[:])
bmesh.update_edit_mesh(obj.data)
bpy.ops.object.mode_set(mode="OBJECT")

bpy.ops.object.select_all(action="SELECT")
bpy.ops.wm.obj_export(
    filepath=str(dst),
    export_selected_objects=True,
    export_materials=False,
    export_uv=False,
)
print(f"Wrote triangulated geometry: {dst}")
