"""
3D OpenGL rendering widget
"""
import os
from PyQt6.QtCore import QTimer
from PyQt6.QtOpenGLWidgets import QOpenGLWidget
from pywavefront import Wavefront, visualization
from OpenGL.GL import *
from OpenGL.GLU import *

class OpenGLWidget(QOpenGLWidget):
    """Widget to display 3D content using OpenGL"""
    def __init__(self, parent=None, model_path=None):
        super().__init__(parent)
        # Set initial rotation to look top-down
        self.rotation_y = 0.0
        self.rotation_x = 0.0  
        self.rotation_z = 0.0
        
        # Load the 3D model
        root_path = os.path.dirname(os.path.dirname(__file__))
        if model_path:
            self.model = Wavefront(model_path)
        else:
            self.model = Wavefront(os.path.join(root_path, 'app/hex.obj'))
            
    
    def initializeGL(self):
        """Set up OpenGL state"""
        glClearColor(0.2, 0.2, 0.2, 1.0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        
        # Set up light to match top-down view
        glLightfv(GL_LIGHT0, GL_POSITION, [0.0, 0.0, 1.0, 0.0])
        
    def resizeGL(self, width, height):
        """Handle window resize"""
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        
        # Use orthographic projection for a true flat view
        # Parameters: left, right, bottom, top, near, far
        aspect = float(width) / height if height > 0 else 1.0
        size = 10.0  # Controls how much of the scene is visible
        
        if width >= height:
            glOrtho(-size * aspect, size * aspect, -size, size, -100, 100)
        else:
            glOrtho(-size, size, -size / aspect, size / aspect, -100, 100)
            
        glMatrixMode(GL_MODELVIEW)
        
    def paintGL(self):
        """Draw the 3D scene"""
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        glLoadIdentity()
        
        # Apply rotations - X rotation of 90 degrees gives top-down view
        glRotatef(self.rotation_x, 1.0, 0.0, 0.0)
        glRotatef(self.rotation_y, 0.0, 1.0, 0.0)
        glRotatef(self.rotation_z, 0.0, 0.0, 1.0)
        
        # Adjust scale if needed
        scale_factor = 1.0
        glScalef(scale_factor, scale_factor, scale_factor)
    
        visualization.draw(self.model)
        
    def set_rotation_y(self, value):
        """Update Y rotation value"""
        self.rotation_y = value
        self.update()
        
    def set_rotation_x(self, value):
        """Update X rotation value"""
        self.rotation_x = value
        self.update()
        
    def set_rotation_z(self, value):
        """Update Z rotation value"""
        self.rotation_z = value
        self.update()