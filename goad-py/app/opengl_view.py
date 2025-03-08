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
        self.rotation_y = 0.0
        self.rotation_x = -25.0
        self.rotation_z = 45.0
        
        # Set up a timer for animation
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update)
        self.timer.start(16)  # ~60 fps
        
        # Load the 3D model
        root_path = os.path.dirname(os.path.dirname(__file__))
        if model_path:
            self.model = Wavefront(model_path)
        else:
            self.model = Wavefront(os.path.join(root_path, 'app/cube.obj'))
        
    def initializeGL(self):
        """Set up OpenGL state"""
        glClearColor(0.2, 0.2, 0.2, 1.0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        
        # Set up light
        glLightfv(GL_LIGHT0, GL_POSITION, [-1.0, 1.0, 1.0, 0.0])
        
    def resizeGL(self, width, height):
        """Handle window resize"""
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(90.0, float(width) / height, 1.0, 100.0)
        glMatrixMode(GL_MODELVIEW)
        
    def paintGL(self):
        """Draw the 3D scene"""
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
        glLoadIdentity()
        glTranslatef(-4.0, 2.0, -20.0)
        glRotatef(self.rotation_y, 0.0, 1.0, 0.0)
        glRotatef(self.rotation_x, 1.0, 0.0, 0.0)
        glRotatef(self.rotation_z, 0.0, 0.0, 1.0)
        
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