"""
Matplotlib plotting widget
"""
import numpy as np
from PyQt6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class MatplotlibWidget(QWidget):
    """Widget to display matplotlib plots"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        
        self.ax = self.figure.add_subplot(111)
        self.update_plot()
        
    def update_plot(self):
        """Generate a sample plot"""
        self.ax.clear()
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        self.ax.plot(x, y)
        self.ax.set_title('Sample Plot')
        self.ax.set_xlabel('X axis')
        self.ax.set_ylabel('Y axis')
        self.ax.grid(True)
        self.canvas.draw()
        
    def create_sine_plot(self):
        """Create a sine wave plot"""
        self.ax.clear()
        x = np.linspace(0, 10, 100)
        y = np.sin(x)
        self.ax.plot(x, y)
        self.ax.set_title('Sine Plot')
        self.ax.set_xlabel('X axis')
        self.ax.set_ylabel('Y axis')
        self.ax.grid(True)
        self.canvas.draw()
        
    def create_cosine_plot(self):
        """Create a cosine wave plot"""
        self.ax.clear()
        x = np.linspace(0, 10, 100)
        y = np.cos(x)
        self.ax.plot(x, 'r-')  # Red line for cosine
        self.ax.set_title('Cosine Plot')
        self.ax.set_xlabel('X axis')
        self.ax.set_ylabel('Y axis')
        self.ax.grid(True)
        self.canvas.draw()
    
    def plot_mueller_data(self, mueller_data):
        """Plot the Mueller matrix data
        
        Args:
            mueller_data: List of lists containing the Mueller matrix elements
        """
        # Convert data to numpy array if it's not already
        mueller_array = np.array(mueller_data)
        
        # Clear the existing plot
        self.ax.clear()
        
        # Create angle array (assuming 0-180 degrees scattering angle)
        angles = np.linspace(0, 180, len(mueller_array))
        
        # Plot Mueller matrix elements
        # S11 element (typically the most significant)
        self.ax.plot(angles, mueller_array[:, 0], 'b-', label='S11')
        
        # Plot other diagonal elements
        # self.ax.plot(angles, mueller_array[:, 5], 'r-', label='S22')
        # self.ax.plot(angles, mueller_array[:, 10], 'g-', label='S33')
        # self.ax.plot(angles, mueller_array[:, 15], 'c-', label='S44')
        
        # Set labels and title
        self.ax.set_xlabel('Scattering Angle (degrees)')
        self.ax.set_ylabel('Mueller Matrix Element')
        self.ax.set_title('GOAD Simulation: Mueller Matrix Elements')
        self.ax.legend()
        self.ax.grid(True)
        self.ax.set_yscale('log')
        
        # Redraw the canvas
        self.canvas.draw()