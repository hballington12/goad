"""
Main application window
"""
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QMainWindow, QDockWidget, QMessageBox)
from PyQt6.QtGui import QAction

from opengl_view import OpenGLWidget
from plot_view import PlotlyWidget
# from plot_view import MatplotlibWidget
from panels import PropertiesPanel, OutlinerPanel, ConsolePanel, SimulationSettingsPanel

class MainWindow(QMainWindow):
    """Main application window with dockable widgets"""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("3D Viewer with Dockable Panels")
        self.setMinimumSize(1280, 720)
        
        # Set up central widget (3D view)
        self.opengl_widget = OpenGLWidget()
        self.setCentralWidget(self.opengl_widget)
        
        # Create dockable widgets
        self.createDockWidgets()
        
        # Create menu bar and actions
        self.createActions()
        self.createMenus()
        self.createToolBars()
        self.createStatusBar()
        
    def createDockWidgets(self):
        """Create dockable UI panels"""
        # Plot panel (right side)
        self.plot_dock = QDockWidget("Plot View", self)
        self.plot_dock.setAllowedAreas(
            Qt.DockWidgetArea.RightDockWidgetArea | 
            Qt.DockWidgetArea.LeftDockWidgetArea | 
            Qt.DockWidgetArea.BottomDockWidgetArea
        )
        self.plot_widget = PlotlyWidget()
        self.plot_dock.setWidget(self.plot_widget)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.plot_dock)
        self.plot_dock.setVisible(False)
        
        # Properties panel (left side)
        self.properties_dock = QDockWidget("Properties", self)
        self.properties_dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
        self.properties_panel = PropertiesPanel()
        self.properties_dock.setWidget(self.properties_panel)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.properties_dock)
        self.properties_dock.setVisible(False)
        
        # Simulation settings panel (left side)
        self.simulation_dock = QDockWidget("Simulation Settings", self)
        self.simulation_dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
        self.simulation_panel = SimulationSettingsPanel()
        self.simulation_dock.setWidget(self.simulation_panel)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.simulation_dock)
        
        # Outliner panel (left side, below properties)
        self.outliner_dock = QDockWidget("Outliner", self)
        self.outliner_dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
        self.outliner_panel = OutlinerPanel()
        self.outliner_dock.setWidget(self.outliner_panel)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.outliner_dock)
        self.outliner_dock.setVisible(False)
        
        # Console panel (bottom)
        self.console_dock = QDockWidget("Console", self)
        self.console_dock.setAllowedAreas(
            Qt.DockWidgetArea.BottomDockWidgetArea | 
            Qt.DockWidgetArea.TopDockWidgetArea
        )
        self.console_panel = ConsolePanel()
        self.console_dock.setWidget(self.console_panel)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.console_dock)
        
        # Tab the dock widgets
        self.tabifyDockWidget(self.properties_dock, self.simulation_dock)
        self.tabifyDockWidget(self.simulation_dock, self.outliner_dock)
        
        # Connect signals
        self.simulation_panel.euler_a.valueChanged.connect(
            self.opengl_widget.set_rotation_x)
        self.simulation_panel.euler_b.valueChanged.connect(
            self.opengl_widget.set_rotation_y)
        self.simulation_panel.euler_g.valueChanged.connect(
            self.opengl_widget.set_rotation_z)
        self.simulation_panel.mueller_data_ready.connect(
            self.plot_widget.plot_mueller_data)
        self.simulation_panel.mueller_data_ready.connect(
            self.showPlotDock)

    def createActions(self):
        """Create application actions"""
        # File menu actions
        self.newAct = QAction("&New", self, shortcut="Ctrl+N",
                             statusTip="Create a new file",
                             triggered=self.newFile)
        
        self.openAct = QAction("&Open...", self, shortcut="Ctrl+O",
                              statusTip="Open an existing file",
                              triggered=self.openFile)
        
        self.saveAct = QAction("&Save", self, shortcut="Ctrl+S",
                              statusTip="Save the document",
                              triggered=self.saveFile)
        
        self.exitAct = QAction("E&xit", self, shortcut="Ctrl+Q",
                              statusTip="Exit the application",
                              triggered=self.close)
        
        # View menu actions
        self.plotViewAct = QAction("&Plot View", self,
                                  statusTip="Show/hide plot panel",
                                  checkable=True,
                                  checked=False,
                                  triggered=self.togglePlotView)
        
        self.controlsAct = QAction("&Controls", self,
                                 statusTip="Show/hide controls panel",
                                 checkable=True,
                                 checked=True,
                                 triggered=self.toggleControls)
        
        self.propertiesAct = QAction("&Properties", self,
                                    statusTip="Show/hide properties panel",
                                    checkable=True,
                                    checked=False,
                                    triggered=self.toggleProperties)
        
        self.outlinerAct = QAction("&Outliner", self,
                                  statusTip="Show/hide outliner panel",
                                  checkable=True,
                                  checked=False,
                                  triggered=self.toggleOutliner)
        
        self.consoleAct = QAction("&Console", self,
                                 statusTip="Show/hide console panel",
                                 checkable=True,
                                 checked=True,
                                 triggered=self.toggleConsole)
        
        self.simulationSettingsAct = QAction("&Simulation Settings", self,
                                          statusTip="Show/hide simulation settings panel",
                                          checkable=True,
                                          checked=True,
                                          triggered=self.toggleSimulationSettings)
        
        # Layout actions
        self.defaultLayoutAct = QAction("&Reset Layout", self,
                                      statusTip="Reset to default layout",
                                      triggered=self.setDefaultLayout)
        
    def createMenus(self):
        """Create application menus"""
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(self.newAct)
        self.fileMenu.addAction(self.openAct)
        self.fileMenu.addAction(self.saveAct)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAct)
        
        self.viewMenu = self.menuBar().addMenu("&View")
        self.viewMenu.addAction(self.plotViewAct)
        self.viewMenu.addAction(self.controlsAct)
        self.viewMenu.addAction(self.propertiesAct)
        self.viewMenu.addAction(self.simulationSettingsAct)
        self.viewMenu.addAction(self.outlinerAct)
        self.viewMenu.addAction(self.consoleAct)
        self.viewMenu.addSeparator()
        
        self.layoutMenu = self.viewMenu.addMenu("&Layouts")
        self.layoutMenu.addAction(self.defaultLayoutAct)
        
        self.plotMenu = self.menuBar().addMenu("&Plot")
    
    def createToolBars(self):
        """Create application toolbars"""
        self.fileToolBar = self.addToolBar("File")
        self.fileToolBar.addAction(self.newAct)
        self.fileToolBar.addAction(self.openAct)
        self.fileToolBar.addAction(self.saveAct)
        
        self.viewToolBar = self.addToolBar("View")
        self.viewToolBar.addAction(self.defaultLayoutAct)
    
    def createStatusBar(self):
        """Create status bar"""
        self.statusBar().showMessage("Ready")
    
    # File actions
    def newFile(self):
        self.statusBar().showMessage("Creating new file...")
        self.console_panel.log("Creating new file")
    
    def openFile(self):
        self.statusBar().showMessage("Opening file...")
        self.console_panel.log("Opening file dialog")
    
    def saveFile(self):
        self.statusBar().showMessage("Saving file...")
        self.console_panel.log("Saving current file")
    
    # View toggle actions
    def togglePlotView(self):
        self.plot_dock.setVisible(self.plotViewAct.isChecked())
    
    def toggleControls(self):
        self.control_dock.setVisible(self.controlsAct.isChecked())
    
    def toggleProperties(self):
        self.properties_dock.setVisible(self.propertiesAct.isChecked())
    
    def toggleOutliner(self):
        self.outliner_dock.setVisible(self.outlinerAct.isChecked())
    
    def toggleConsole(self):
        self.console_dock.setVisible(self.consoleAct.isChecked())
    
    def toggleSimulationSettings(self):
        self.simulation_dock.setVisible(self.simulationSettingsAct.isChecked())
    
    def showPlotDock(self):
        if not self.plot_dock.isVisible():
            self.plot_dock.setVisible(True)
    
    # Layout actions
    def setDefaultLayout(self):
        self.statusBar().showMessage("Switching to default layout...")
        self.console_panel.log("Switched to default layout")
        
        # Re-add all widgets to their default locations
        for dock in [self.plot_dock, self.properties_dock, 
                    self.simulation_dock, self.outliner_dock, self.console_dock]:
            if dock.isFloating():
                dock.setFloating(False)
                
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.plot_dock)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.properties_dock)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.simulation_dock)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.outliner_dock)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.console_dock)
        self.tabifyDockWidget(self.properties_dock, self.simulation_dock)
        self.tabifyDockWidget(self.simulation_dock, self.outliner_dock)

        self.properties_dock.setVisible(False)
        self.outliner_dock.setVisible(False)
        self.plot_dock.setVisible(False)