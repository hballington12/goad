"""
Matplotlib plotting widget
"""
import numpy as np
from PyQt6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

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
        
    def plot_mueller_data(self, data):
        """Plot the Mueller matrix data
        
        Args:
            mueller_data: List of lists containing the Mueller matrix elements
        """

        # Unpack the data
        angles, mueller_data = data

        # Convert data to numpy array if it's not already
        mueller_array = np.array(mueller_data)
        
        # Clear the existing plot
        self.ax.clear()
        
        # Plot configuration options (from mueller_plotter.py)
        FS = 16  # font size
        LW = 2   # line width
        GRID = True
        XLABEL = r"$\theta$ / $^\circ$"  # x-axis label
        YLABEL = "Phase Function"  # y-axis label
        LEGEND = True   

        # Enhanced styling - set figure style
        self.figure.patch.set_facecolor('#E8E8E8')  # Light gray figure background
        self.ax.set_facecolor('#F5F5F5')  # Slightly lighter gray plot background
        
        # Apply plot configurations
        plt.rcParams['font.size'] = FS
        plt.rcParams['axes.linewidth'] = LW
        
        # Plot Mueller matrix elements with improved styling
        # Use a color palette for more appealing visuals
        self.ax.plot(angles, mueller_array[:, 0], color='#1f77b4', linewidth=LW+0.5, label='S11', 
                     alpha=0.9)
        
        # Set labels and styling
        self.ax.set_xlabel(XLABEL, fontsize=FS, fontweight='bold')
        self.ax.set_ylabel(YLABEL, fontsize=FS, fontweight='bold')
        self.ax.set_title('GOAD Simulation: Mueller Matrix Elements', 
                          fontsize=FS+2, fontweight='bold', pad=10)
        
        # Enhance ticks and labels
        self.ax.tick_params(axis='both', which='major', labelsize=FS-2, width=1.5, length=6)
        self.ax.tick_params(axis='both', which='minor', width=1, length=3)
        
        # Add a box around the plot
        for spine in self.ax.spines.values():
            spine.set_linewidth(1.5)
            spine.set_color('#555555')
        
        # Apply grid and axis limits (from mueller_plotter.py)
        if GRID:
            self.ax.grid(True, linestyle='--', alpha=0.7, color='#888888')
        self.ax.set_xlim(0, 180)
        self.ax.set_yscale('log')
        
        # Customize the legend with border similar to mueller_plotter.py
        if LEGEND:
            legend = self.ax.legend(loc='upper right', fontsize=FS-2, 
                                   edgecolor='#555555', frameon=True, 
                                   fancybox=True, shadow=True, framealpha=0.9)
            legend.get_frame().set_linewidth(1.5)
        
        # Redraw the canvas
        self.canvas.draw()