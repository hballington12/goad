"""
Outliner panel widget
"""
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QTreeWidget, QTreeWidgetItem

class OutlinerPanel(QWidget):
    """Widget for displaying scene hierarchy"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()
        
    def initUI(self):
        layout = QVBoxLayout()
        
        # Title
        title_label = QLabel("<b>Scene Outliner</b>")
        
        # Scene objects
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Name", "Type"])
        
        # Add some example items
        scene_root = QTreeWidgetItem(["Scene", "Root"])
        model_item = QTreeWidgetItem(["Cube", "Mesh"])
        camera_item = QTreeWidgetItem(["Main Camera", "Camera"])
        light_item = QTreeWidgetItem(["Light", "Light"])
        
        scene_root.addChild(model_item)
        scene_root.addChild(camera_item)
        scene_root.addChild(light_item)
        
        self.tree.addTopLevelItem(scene_root)
        self.tree.expandAll()
        
        # Add widgets to layout
        layout.addWidget(title_label)
        layout.addWidget(self.tree)
        
        self.setLayout(layout)