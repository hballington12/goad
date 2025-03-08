"""
Simulation settings panel widget
"""
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QFormLayout, QDoubleSpinBox, QComboBox, QCheckBox, QLabel, QPushButton
import goad_py as goad

class SimulationSettingsPanel(QWidget):
    """Panel for simulation settings"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()
        
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
    
    def run_simulation(self):
        """Handle the Run Simulation button click"""
        print(f"Running simulation with:")
        print(f"- Time step: {self.time_step.value()} seconds")
        print(f"- Algorithm: {self.algorithm.currentText()}")
        print(f"- Gravity enabled: {self.gravity_enabled.isChecked()}")
        print(f"- Gravity value: {self.gravity_value.value()} m/s²")

        self.run_goad_sample()

    def run_goad_sample(self):
        print("creating shape object")
        vertices = [
            (6.025435, 6.025435, -6.025435),
            (6.025435, 6.025435, 6.025435),
            (6.025435, -6.025435, -6.025435),
            (6.025435, -6.025435, 6.025435),
            (-6.025435, 6.025435, -6.025435),
            (-6.025435, 6.025435, 6.025435),
            (-6.025435, -6.025435, -6.025435),
            (-6.025435, -6.025435, 6.025435)
        ]
        faces = [
            (0, 1, 3, 2),
            (2, 3, 7, 6),
            (6, 7, 5, 4),
            (4, 5, 1, 0),
            (2, 6, 4, 0),
            (7, 3, 1, 5)
        ]
        shape_id = 0 # note that we need care to deal with the containment graph
        # the shape_id for each face is currently bodged in impl for pymethods on Shape struct
        refr_re = 1.31
        refr_im = 0.0

        shape = goad.Shape(vertices, faces, shape_id, refr_re, refr_im)

        shapes = [shape]

        print("creating geometry object")

        geom = goad.Geom(shapes)

        print("creating goad settings")

        settings = goad.Settings(
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

        problem = goad.Problem(geom, settings)

        print("solving goad problem")

        problem.py_solve()

        problem.py_print_stats()
