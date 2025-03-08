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
        self.ax.plot(x, y, 'r-')  # Red line for cosine
        self.ax.set_title('Cosine Plot')
        self.ax.set_xlabel('X axis')
        self.ax.set_ylabel('Y axis')
        self.ax.grid(True)
        self.canvas.draw()