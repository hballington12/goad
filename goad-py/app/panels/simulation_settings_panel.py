"""
Simulation settings panel widget
"""
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QFormLayout, QDoubleSpinBox, QComboBox, 
                            QCheckBox, QLabel, QPushButton, QSlider, QHBoxLayout)
from PyQt6.QtCore import pyqtSignal, Qt
import goad_py as goad
import numpy as np

class SimulationSettingsPanel(QWidget):
    """Panel for simulation settings"""
    # Define a signal to emit when Mueller data is available
    mueller_data_ready = pyqtSignal(object)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()
        self.setup_goad_configuration()
        self.connect_signals()
        
    def initUI(self):
        layout = QVBoxLayout()
        form_layout = QFormLayout()
        
        # Time step setting
        self.time_step = QDoubleSpinBox()
        self.time_step.setRange(0.001, 1.0)
        self.time_step.setSingleStep(0.001)
        self.time_step.setValue(0.01)
        self.time_step.setDecimals(3)
        form_layout.addRow("Time Step (s):", self.time_step)
        
        # Simulation algorithm
        self.algorithm = QComboBox()
        self.algorithm.addItems(["Euler", "Verlet", "RK4"])
        form_layout.addRow("Algorithm:", self.algorithm)
        
        # Gravity toggle
        self.gravity_enabled = QCheckBox()
        self.gravity_enabled.setChecked(True)
        form_layout.addRow("Gravity Enabled:", self.gravity_enabled)
        
        # Gravity value
        self.gravity_value = QDoubleSpinBox()
        self.gravity_value.setRange(0.0, 100.0)
        self.gravity_value.setSingleStep(0.1)
        self.gravity_value.setValue(9.8)
        self.gravity_value.setDecimals(1)
        form_layout.addRow("Gravity (m/s²):", self.gravity_value)
        
        # Add the form layout to the main layout
        layout.addLayout(form_layout)
        
        # Add Euler angle sliders in separate section
        euler_layout = QVBoxLayout()
        euler_layout.setSpacing(5)
        
        # Euler A slider with value label
        euler_a_layout = QHBoxLayout()
        euler_a_label = QLabel("Euler A (deg):")
        self.euler_a_value = QLabel("0.0")
        self.euler_a_value.setMinimumWidth(40)
        self.euler_a_value.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        
        self.euler_a = QSlider(Qt.Orientation.Horizontal)
        self.euler_a.setRange(0, 360)
        self.euler_a.setValue(0)
        self.euler_a.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.euler_a.setTickInterval(45)
        
        euler_a_layout.addWidget(euler_a_label)
        euler_a_layout.addWidget(self.euler_a)
        euler_a_layout.addWidget(self.euler_a_value)
        euler_layout.addLayout(euler_a_layout)
        
        # Euler B slider with value label
        euler_b_layout = QHBoxLayout()
        euler_b_label = QLabel("Euler B (deg):")
        self.euler_b_value = QLabel("0.0")
        self.euler_b_value.setMinimumWidth(40)
        self.euler_b_value.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        
        self.euler_b = QSlider(Qt.Orientation.Horizontal)
        self.euler_b.setRange(0, 360)
        self.euler_b.setValue(0)
        self.euler_b.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.euler_b.setTickInterval(45)
        
        euler_b_layout.addWidget(euler_b_label)
        euler_b_layout.addWidget(self.euler_b)
        euler_b_layout.addWidget(self.euler_b_value)
        euler_layout.addLayout(euler_b_layout)
        
        # Euler G slider with value label
        euler_g_layout = QHBoxLayout()
        euler_g_label = QLabel("Euler G (deg):")
        self.euler_g_value = QLabel("0.0")
        self.euler_g_value.setMinimumWidth(40)
        self.euler_g_value.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        
        self.euler_g = QSlider(Qt.Orientation.Horizontal)
        self.euler_g.setRange(0, 360)
        self.euler_g.setValue(0)
        self.euler_g.setTickPosition(QSlider.TickPosition.TicksBelow)
        self.euler_g.setTickInterval(45)
        
        euler_g_layout.addWidget(euler_g_label)
        euler_g_layout.addWidget(self.euler_g)
        euler_g_layout.addWidget(self.euler_g_value)
        euler_layout.addLayout(euler_g_layout)
        
        # Add Euler layout to main layout
        layout.addLayout(euler_layout)
        layout.addSpacing(10)
        
        # Add checkbox for live update
        self.live_update = QCheckBox("Live Update")
        self.live_update.setToolTip("Automatically run simulation when settings change")
        self.live_update.setChecked(False)
        layout.addWidget(self.live_update)
        
        # Add Run Simulation button
        self.run_simulation_btn = QPushButton("Run Simulation")
        self.run_simulation_btn.clicked.connect(self.run_simulation)
        layout.addSpacing(10)
        layout.addWidget(self.run_simulation_btn)
        
        # Add some spacing and a help label
        layout.addSpacing(20)
        help_label = QLabel("Adjust simulation parameters above to control physics behavior")
        help_label.setWordWrap(True)
        layout.addWidget(help_label)
        
        # Add stretch to push everything to the top
        layout.addStretch(1)
        
        self.setLayout(layout)
    
    def connect_signals(self):
        """Connect widget signals to update handlers"""
        # Connect Euler angle controls to their update functions
        self.euler_a.valueChanged.connect(self.update_euler_a_display)
        self.euler_b.valueChanged.connect(self.update_euler_b_display)
        self.euler_g.valueChanged.connect(self.update_euler_g_display)
        
        # Connect other controls if they should trigger live updates
        self.time_step.valueChanged.connect(self.on_setting_changed)
        self.algorithm.currentIndexChanged.connect(self.on_setting_changed)
        self.gravity_enabled.stateChanged.connect(self.on_setting_changed)
        self.gravity_value.valueChanged.connect(self.on_setting_changed)
    
    def update_euler_a_display(self, value):
        """Update the Euler A display and GOAD settings"""
        angle = value
        self.euler_a_value.setText(f"{angle:.1f}")
        
        # Update GOAD settings
        current_euler = self.settings.euler
        new_euler = [angle, current_euler[1], current_euler[2]]
        self.settings.euler = new_euler
        
        # Run simulation if live update is enabled
        if self.live_update.isChecked():
            self.run_simulation()
    
    def update_euler_b_display(self, value):
        """Update the Euler B display and GOAD settings"""
        angle = value
        self.euler_b_value.setText(f"{angle:.1f}")
        
        # Update GOAD settings
        current_euler = self.settings.euler
        new_euler = [current_euler[0], angle, current_euler[2]]
        self.settings.euler = new_euler
        self.problem.settings = self.settings
        
        # Run simulation if live update is enabled
        if self.live_update.isChecked():
            self.run_simulation()
    
    def update_euler_g_display(self, value):
        """Update the Euler G display and GOAD settings"""
        angle = value
        self.euler_g_value.setText(f"{angle:.1f}")
        
        # Update GOAD settings
        current_euler = self.settings.euler
        new_euler = [current_euler[0], current_euler[1], angle]
        self.settings.euler = new_euler
        
        # Run simulation if live update is enabled
        if self.live_update.isChecked():
            self.run_simulation()
    
    def on_setting_changed(self):
        """Called when any setting is changed to potentially trigger a simulation run"""
        # Only run if live update is enabled
        if self.live_update.isChecked():
            self.run_simulation()
    
    def setup_goad_configuration(self):
        """Initialize the GOAD problem with default settings"""
        print("Creating GOAD settings and initializing problem...")
        
        self.settings = goad.Settings(
            wavelength=0.532,
            beam_power_threshold=1e-1,
            beam_area_threshold_fac=1e-1,
            total_power_cutoff=0.99,
            medium_refr_index_re=1.0,
            medium_refr_index_im=0.0,
            particle_refr_index_re=1.31,
            particle_refr_index_im=0.0,
            geom_name="hex.obj",
            max_rec=10,
            max_tir=10,
            theta_res=181,
            phi_res=181,
            euler=[0.0, 0.0, 0.0]
        )
        
        # Initialize the problem object so we don't need to recreate it for each run
        print("Creating GOAD problem...")
        self.problem = goad.Problem(self.settings)
        
        print("GOAD settings initialized and ready.")
    
    def run_simulation(self):
        """Handle the Run Simulation button click"""
        print(f"Running simulation with:")
        print(f"- Time step: {self.time_step.value()} seconds")
        print(f"- Algorithm: {self.algorithm.currentText()}")
        print(f"- Gravity enabled: {self.gravity_enabled.isChecked()}")
        print(f"- Gravity value: {self.gravity_value.value()} m/s²")
        print(f"- Euler angles: [{self.euler_a.value()}, {self.euler_b.value()}, {self.euler_g.value()}]")

        self.run_goad_sample()

    def run_goad_sample(self):
        """Run the GOAD simulation using the pre-initialized problem"""

        print("Creating GOAD problem...")
        self.problem = goad.Problem(self.settings)
        
        print("Applying rotation")
        self.problem.py_rerotate()

        print("Solving GOAD problem...")
        
        # Solving the pre-initialized problem
        self.problem.py_solve()
        
        print("Problem solved. Printing statistics...")
        self.problem.py_print_stats()
        
        # Get the Mueller data
        mueller_1d = self.problem.mueller_1d
        
        # Emit the signal with the Mueller data
        self.mueller_data_ready.emit(mueller_1d)
        print("Mueller data emitted.")
