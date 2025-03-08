"""
Properties panel widget
"""
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QFrame
from PyQt6.QtCore import Qt

class PropertiesPanel(QWidget):
    """Widget for displaying and editing object properties"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()
        
    def initUI(self):
        layout = QVBoxLayout()
        
        # Title
        title_label = QLabel("<b>Object Properties</b>")
        
        # Properties (example)
        self.vertices_label = QLabel("Vertices: 8")
        self.faces_label = QLabel("Faces: 12")
        self.material_label = QLabel("Material: Default")
        
        # Add widgets to layout
        layout.addWidget(title_label)
        layout.addWidget(QLabel("Statistics:"))
        layout.addWidget(self.vertices_label)
        layout.addWidget(self.faces_label)
        layout.addWidget(self.material_label)
        # Add a separator
        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setFrameShadow(QFrame.Shadow.Sunken)
        layout.addWidget(separator)
        layout.addWidget(separator)
        
        # Display options
        layout.addWidget(QLabel("Display Options:"))
        
        # Wireframe toggle
        self.wireframe_check = QLabel("Wireframe mode: Off")
        layout.addWidget(self.wireframe_check)
        
        # Add stretch to push widgets to the top
        layout.addStretch(1)
        
        self.setLayout(layout)