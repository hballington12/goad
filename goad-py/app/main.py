"""
3D Viewer with integrated plots - Application entry point
"""
import sys
from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QSurfaceFormat
from main_window import MainWindow

def main():
    # Set OpenGL format (important for proper functioning)
    fmt = QSurfaceFormat()
    fmt.setDepthBufferSize(24)
    fmt.setSamples(4)  # Anti-aliasing
    QSurfaceFormat.setDefaultFormat(fmt)
    
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()