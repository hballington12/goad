"""
Console panel widget
"""
from PyQt6.QtWidgets import QWidget, QVBoxLayout, QLabel, QTextEdit

class ConsolePanel(QWidget):
    """Widget for displaying console output"""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()
        
    def initUI(self):
        layout = QVBoxLayout()
        
        # Title
        title_label = QLabel("<b>Console Output</b>")
        
        # Text area for console output
        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setStyleSheet("background-color: #1e1e1e; color: #dcdcdc;")
        self.console.append("Application started.")
        self.console.append("OpenGL initialized.")
        self.console.append("Model loaded: cube.obj")
        
        # Command input
        self.command_input = QLabel("Command input would go here")
        
        # Add widgets to layout
        layout.addWidget(title_label)
        layout.addWidget(self.console)
        layout.addWidget(self.command_input)
        
        self.setLayout(layout)
        
    def log(self, message):
        """Add a message to the console output"""
        self.console.append(message)