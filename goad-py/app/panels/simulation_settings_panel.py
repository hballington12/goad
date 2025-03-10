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
        
        # Simulation algorithm
        self.algorithm = QComboBox()
        self.algorithm.addItems(["Euler XYZ"])
        form_layout.addRow("Algorithm:", self.algorithm)
        
        # Wavelength control
        self.wavelength = QDoubleSpinBox()
        self.wavelength.setRange(0.1, 10.0)
        self.wavelength.setValue(0.532)
        self.wavelength.setDecimals(3)
        self.wavelength.setSingleStep(0.001)
        self.wavelength.setSuffix(" μm")
        form_layout.addRow("Wavelength:", self.wavelength)
        
        # Beam power threshold
        self.beam_power_threshold = QDoubleSpinBox()
        self.beam_power_threshold.setRange(1e-6, 1.0)
        self.beam_power_threshold.setValue(1e-1)
        self.beam_power_threshold.setDecimals(6)
        self.beam_power_threshold.setSingleStep(0.01)
        form_layout.addRow("Beam Power Threshold:", self.beam_power_threshold)
        
        # Beam area threshold factor
        self.beam_area_threshold_fac = QDoubleSpinBox()
        self.beam_area_threshold_fac.setRange(1e-6, 1.0)
        self.beam_area_threshold_fac.setValue(1e-1)
        self.beam_area_threshold_fac.setDecimals(6)
        self.beam_area_threshold_fac.setSingleStep(0.01)
        form_layout.addRow("Beam Area Threshold:", self.beam_area_threshold_fac)
        
        # Power cutoff
        self.total_power_cutoff = QDoubleSpinBox()
        self.total_power_cutoff.setRange(0.01, 0.9999)
        self.total_power_cutoff.setValue(0.99)
        self.total_power_cutoff.setDecimals(4)
        self.total_power_cutoff.setSingleStep(0.01)
        form_layout.addRow("Power Cutoff:", self.total_power_cutoff)
        
        # Medium refractive index
        self.medium_refr_index_re = QDoubleSpinBox()
        self.medium_refr_index_re.setRange(0.1, 10.0)
        self.medium_refr_index_re.setValue(1.0)
        self.medium_refr_index_re.setDecimals(4)
        self.medium_refr_index_re.setSingleStep(0.01)
        form_layout.addRow("Medium Refr. Index (Real):", self.medium_refr_index_re)
        
        self.medium_refr_index_im = QDoubleSpinBox()
        self.medium_refr_index_im.setRange(0.0, 10.0)
        self.medium_refr_index_im.setValue(0.0)
        self.medium_refr_index_im.setDecimals(4)
        self.medium_refr_index_im.setSingleStep(0.001)
        form_layout.addRow("Medium Refr. Index (Imag):", self.medium_refr_index_im)
        
        # Particle refractive index
        self.particle_refr_index_re = QDoubleSpinBox()
        self.particle_refr_index_re.setRange(0.1, 10.0)
        self.particle_refr_index_re.setValue(1.31)
        self.particle_refr_index_re.setDecimals(4)
        self.particle_refr_index_re.setSingleStep(0.01)
        form_layout.addRow("Particle Refr. Index (Real):", self.particle_refr_index_re)
        
        self.particle_refr_index_im = QDoubleSpinBox()
        self.particle_refr_index_im.setRange(0.0, 10.0)
        self.particle_refr_index_im.setValue(0.0)
        self.particle_refr_index_im.setDecimals(4)
        self.particle_refr_index_im.setSingleStep(0.001)
        form_layout.addRow("Particle Refr. Index (Imag):", self.particle_refr_index_im)
        
        # Geometry selection
        self.geom_name = QComboBox()
        self.geom_name.addItems(["hex.obj", "sphere.obj", "cylinder.obj"])
        form_layout.addRow("Geometry:", self.geom_name)
        
        # Max recursion
        self.max_rec = QDoubleSpinBox()
        self.max_rec.setRange(1, 100)
        self.max_rec.setValue(10)
        self.max_rec.setDecimals(0)
        self.max_rec.setSingleStep(1)
        form_layout.addRow("Max Recursion:", self.max_rec)
        
        # Max TIR
        self.max_tir = QDoubleSpinBox()
        self.max_tir.setRange(1, 100)
        self.max_tir.setValue(10)
        self.max_tir.setDecimals(0)
        self.max_tir.setSingleStep(1)
        form_layout.addRow("Max TIR:", self.max_tir)
        
        # Resolution
        self.theta_res = QDoubleSpinBox()
        self.theta_res.setRange(10, 1000)
        self.theta_res.setValue(181)
        self.theta_res.setDecimals(0)
        self.theta_res.setSingleStep(10)
        form_layout.addRow("Theta Resolution:", self.theta_res)
        
        self.phi_res = QDoubleSpinBox()
        self.phi_res.setRange(10, 1000)
        self.phi_res.setValue(181)
        self.phi_res.setDecimals(0)
        self.phi_res.setSingleStep(10)
        form_layout.addRow("Phi Resolution:", self.phi_res)
        
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
        self.euler_a.valueChanged.connect(self.on_setting_changed)
        self.euler_b.valueChanged.connect(self.on_setting_changed)
        self.euler_g.valueChanged.connect(self.on_setting_changed)
        
        # Connect other controls to update settings
        self.wavelength.valueChanged.connect(self.update_settings)
        self.beam_power_threshold.valueChanged.connect(self.update_settings)
        self.beam_area_threshold_fac.valueChanged.connect(self.update_settings)
        self.total_power_cutoff.valueChanged.connect(self.update_settings)
        self.medium_refr_index_re.valueChanged.connect(self.update_settings)
        self.medium_refr_index_im.valueChanged.connect(self.update_settings)
        self.particle_refr_index_re.valueChanged.connect(self.update_settings)
        self.particle_refr_index_im.valueChanged.connect(self.update_settings)
        self.geom_name.currentTextChanged.connect(self.update_settings)
        self.max_rec.valueChanged.connect(self.update_settings)
        self.max_tir.valueChanged.connect(self.update_settings)
        self.theta_res.valueChanged.connect(self.update_settings)
        self.phi_res.valueChanged.connect(self.update_settings)
        
        # Connect other controls if they should trigger live updates
        self.algorithm.currentIndexChanged.connect(self.on_setting_changed)

    def update_euler_a_display(self, value):
        """Update the Euler A display and GOAD settings"""
        angle = value
        self.euler_a_value.setText(f"{angle:.1f}")
        
        # Update GOAD settings
        current_euler = self.settings.euler
        new_euler = [angle, current_euler[1], current_euler[2]]
        self.on_new_euler_angles(new_euler)
        
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
        self.on_new_euler_angles(new_euler)
        
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
        self.on_new_euler_angles(new_euler)
        
        # Run simulation if live update is enabled
        if self.live_update.isChecked():
            self.run_simulation()
    
    def on_new_euler_angles(self, euler):
        self.settings.euler = euler

        # Update problem settings
        problem_settings = self.problem.settings
        problem_settings.euler = euler
        self.problem.settings = problem_settings

    def on_setting_changed(self):
        """Called when any setting is changed to potentially trigger a simulation run"""
        # Only run if live update is enabled
        if self.live_update.isChecked():
            self.run_simulation()
    
    def update_settings(self):
        """Update all GOAD settings from UI controls"""
        print("Updating GOAD settings...")
        
        # Create new settings object with values from UI controls
        self.settings = goad.Settings(
            wavelength=self.wavelength.value(),
            beam_power_threshold=self.beam_power_threshold.value(),
            beam_area_threshold_fac=self.beam_area_threshold_fac.value(),
            total_power_cutoff=self.total_power_cutoff.value(),
            medium_refr_index_re=self.medium_refr_index_re.value(),
            medium_refr_index_im=self.medium_refr_index_im.value(),
            particle_refr_index_re=self.particle_refr_index_re.value(),
            particle_refr_index_im=self.particle_refr_index_im.value(),
            geom_name=self.geom_name.currentText(),
            max_rec=int(self.max_rec.value()),
            max_tir=int(self.max_tir.value()),
            theta_res=int(self.theta_res.value()),
            phi_res=int(self.phi_res.value()),
            euler=[self.euler_a.value(), self.euler_b.value(), self.euler_g.value()]
        )
        
        # Update the problem with the new settings
        self.problem.settings = self.settings
        
        # Signal that settings have changed
        self.on_setting_changed()

    def setup_goad_configuration(self):
        """Initialize the GOAD problem with default settings"""
        print("Creating GOAD settings and initializing problem...")

        # Create new settings object with values from UI controls
        self.settings = goad.Settings(
            wavelength=self.wavelength.value(),
            beam_power_threshold=self.beam_power_threshold.value(),
            beam_area_threshold_fac=self.beam_area_threshold_fac.value(),
            total_power_cutoff=self.total_power_cutoff.value(),
            medium_refr_index_re=self.medium_refr_index_re.value(),
            medium_refr_index_im=self.medium_refr_index_im.value(),
            particle_refr_index_re=self.particle_refr_index_re.value(),
            particle_refr_index_im=self.particle_refr_index_im.value(),
            geom_name=self.geom_name.currentText(),
            max_rec=int(self.max_rec.value()),
            max_tir=int(self.max_tir.value()),
            theta_res=int(self.theta_res.value()),
            phi_res=int(self.phi_res.value()),
            euler=[self.euler_a.value(), self.euler_b.value(), self.euler_g.value()]
        )
        
        # Initialize the problem object so we don't need to recreate it for each run
        print("Creating GOAD problem...")
        self.problem = goad.Problem(self.settings)
        print("GOAD settings initialized and ready.")
    
    def run_simulation(self):
        """Handle the Run Simulation button click"""

        self.run_goad()

    def run_goad(self):
        """Run the GOAD simulation using the pre-initialized problem"""

        print("Creating GOAD problem...")
        self.problem = goad.Problem(self.settings)
        
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
