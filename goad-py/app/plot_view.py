"""
Plotly plotting widget with persistent figure and data updates - Dark Mode
"""
import numpy as np
from PyQt6.QtWidgets import QWidget, QVBoxLayout
from plotly.graph_objects import Figure
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from PyQt6.QtWebEngineWidgets import QWebEngineView

class PlotlyWidget(QWidget):
    """Widget to display Plotly plots"""
    def __init__(self, parent=None):
        super().__init__(parent)
        # Create a web view to display the Plotly plot
        self.web_view = QWebEngineView()
        
        # Set up the layout
        layout = QVBoxLayout()
        layout.addWidget(self.web_view)
        self.setLayout(layout)
        
        # Initialize an empty figure with proper styling
        self.initialize_figure()
        
        # Render the empty figure initially
        self.update_plot()
        
    def initialize_figure(self):
        """Initialize the figure with dark theme styling but no data"""
        self.fig = go.Figure()
        
        # Add an empty scatter trace that will be updated later
        self.fig.add_trace(
            go.Scatter(
                x=[], 
                y=[],
                mode='lines',
                line=dict(color='#00bfff', width=3),  # Bright blue for dark theme
                name='S11',
                hovertemplate='Angle: %{x:.2f}°<br>S11: %{y:.4e}<extra></extra>'
            )
        )
        
        # Set up the layout once with dark theme
        self.fig.update_layout(
            title={
                'text': 'Mueller Matrix Elements',
                'font': {'size': 18, 'color': '#ffffff'},
                'y': 0.95
            },
            xaxis_title={
                'text': 'Angle (degrees)',
                'font': {'size': 14, 'color': '#e0e0e0'}
            },
            yaxis_title={
                'text': 'Intensity',
                'font': {'size': 14, 'color': '#e0e0e0'}
            },
            yaxis_type='log',
            paper_bgcolor='#1e1e1e',  # Dark background color
            plot_bgcolor='#282828',   # Dark plot area
            hovermode='closest',
            margin=dict(l=60, r=30, t=80, b=60),
            legend=dict(
                x=0.02,
                y=0.98,
                bgcolor='rgba(40, 40, 40, 0.8)',
                bordercolor='#555555',
                font=dict(color='#e0e0e0')
            ),
            template='plotly_dark',
            font=dict(color='#e0e0e0')  # Default font color for all text
        )
        
        # Add grid lines with dark theme colors
        self.fig.update_xaxes(
            showgrid=True, 
            gridwidth=1, 
            gridcolor='#3a3a3a',
            linecolor='#555555',
            tickfont=dict(color='#e0e0e0')
        )
        self.fig.update_yaxes(
            showgrid=True, 
            gridwidth=1, 
            gridcolor='#3a3a3a',
            linecolor='#555555',
            tickfont=dict(color='#e0e0e0')
        )
    
    def update_plot(self):
        """Render the current figure state to the web view"""
        self.web_view.setHtml(self.fig.to_html(
            include_plotlyjs='cdn', 
            config={
                'responsive': True,
                'scrollZoom': True,
                'displayModeBar': True,
                'displaylogo': False,
                'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
                'toImageButtonOptions': {
                    'format': 'png',
                    'filename': 'mueller_plot',
                    'height': 800,
                    'width': 1200,
                    'scale': 2
                }
            }
        ))
        
    def plot_mueller_data(self, data):
        """Update the Mueller matrix data on the existing figure
        
        Args:
            data: Tuple containing (angles, mueller_data)
        """
        # Unpack the data
        angles, mueller_data = data

        # Convert data to numpy array if it's not already
        mueller_array = np.array(mueller_data)
        
        # Update the existing trace data (more efficient than recreating)
        self.fig.update_traces(
            x=angles,
            y=mueller_array[:, 0],
            selector=dict(name='S11')
        )
        
        # Add animation
        self.fig.update_layout(
            transition_duration=500,
            transition={'easing': 'cubic-in-out'}
        )
        
        # Update the plot in the web view
        self.update_plot()