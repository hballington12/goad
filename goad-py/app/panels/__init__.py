"""
Import all panel modules
"""
from .control_panel import ControlPanel
from .properties_panel import PropertiesPanel
from .outliner_panel import OutlinerPanel
from .console_panel import ConsolePanel
from .simulation_settings_panel import SimulationSettingsPanel

__all__ = ['ControlPanel', 'PropertiesPanel', 'OutlinerPanel', 'ConsolePanel', 'SimulationSettingsPanel']