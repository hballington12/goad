"""
Control panel widget
"""
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QSlider, QLabel, QPushButton

class ControlPanel(QWidget):
    """Widget containing UI controls"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()
        
    def initUI(self):
        layout = QVBoxLayout()
        
        # X-axis Rotation slider
        x_slider_layout = QHBoxLayout()
        self.x_rotation_label = QLabel("X Rotation (degrees):")
        self.x_rotation_slider = QSlider(Qt.Orientation.Horizontal)
        self.x_rotation_slider.setMinimum(0)
        self.x_rotation_slider.setMaximum(360)
        self.x_rotation_slider.setValue(0)
        self.x_rotation_value = QLabel("0°")
        
        x_slider_layout.addWidget(self.x_rotation_label)
        x_slider_layout.addWidget(self.x_rotation_slider)
        x_slider_layout.addWidget(self.x_rotation_value)
        
        # Y-axis Rotation slider (original slider)
        y_slider_layout = QHBoxLayout()
        self.y_rotation_label = QLabel("Y Rotation (degrees):")
        self.y_rotation_slider = QSlider(Qt.Orientation.Horizontal)
        self.y_rotation_slider.setMinimum(0)
        self.y_rotation_slider.setMaximum(360)
        self.y_rotation_slider.setValue(0)
        self.y_rotation_value = QLabel("0°")
        
        y_slider_layout.addWidget(self.y_rotation_label)
        y_slider_layout.addWidget(self.y_rotation_slider)
        y_slider_layout.addWidget(self.y_rotation_value)
        
        # Z-axis Rotation slider
        z_slider_layout = QHBoxLayout()
        self.z_rotation_label = QLabel("Z Rotation (degrees):")
        self.z_rotation_slider = QSlider(Qt.Orientation.Horizontal)
        self.z_rotation_slider.setMinimum(0)
        self.z_rotation_slider.setMaximum(360)
        self.z_rotation_slider.setValue(0)
        self.z_rotation_value = QLabel("0°")
        
        z_slider_layout.addWidget(self.z_rotation_label)
        z_slider_layout.addWidget(self.z_rotation_slider)
        z_slider_layout.addWidget(self.z_rotation_value)
        
        # Update plot button
        self.update_plot_btn = QPushButton("Update Plot")
        
        # Add widgets to layout
        layout.addLayout(x_slider_layout)
        layout.addLayout(y_slider_layout)
        layout.addLayout(z_slider_layout)
        layout.addWidget(self.update_plot_btn)
        layout.addStretch(1)  # Add stretch to push widgets to the top
        
        self.setLayout(layout)
        
        # Connect signals
        self.x_rotation_slider.valueChanged.connect(self.update_x_rotation_label)
        self.y_rotation_slider.valueChanged.connect(self.update_y_rotation_label)
        self.z_rotation_slider.valueChanged.connect(self.update_z_rotation_label)
        
    def update_x_rotation_label(self, value):
        self.x_rotation_value.setText(f"{value}°")
        
    def update_y_rotation_label(self, value):
        self.y_rotation_value.setText(f"{value}°")
        
    def update_z_rotation_label(self, value):
        self.z_rotation_value.setText(f"{value}°")