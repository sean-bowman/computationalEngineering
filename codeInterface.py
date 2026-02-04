
# -- Code Interface (Entry-Point) -- #

'''
Main interface for the Surfboard Design and Analysis Toolkit.

Provides a single entry point to:
    - Generate 3D surfboard geometry via the C# SurfboardGeometry module
    - Run hydrodynamic physics analysis via the Python SurfPhysics module
    - Visualize results with interactive Plotly dashboards

Each step is controlled by flags in the JSON config file:
    - geometry.generate    : run C# STL generation (requires .NET + PicoGK)
    - analysis.showDashboard : open interactive Plotly dashboard

To run, execute this script directly or in debug mode.
Swap the configPath variable to analyze different board configurations.

Sean Bowman [02/04/2026]
'''

import json
import os
from SurfboardGeometry.geometryGenerator import GeometryGenerator
from SurfPhysics.analyzer import PhysicsAnalyzer
from SurfPhysics.visualization.visualizer import Visualizer

# Clear terminal
os.system('cls')

# Configuration — swap this path to run different analyses
configPath = 'configs/shortboard_default.json'

# Load config flags
with open(configPath, 'r') as f:
    config = json.load(f)

# Step 1: Generate 3D geometry via C# module
if config.get('geometry', {}).get('generate', False):
    geometryGenerator = GeometryGenerator()
    geometryGenerator.generateFromConfig(configPath)

# Step 2: Run hydrodynamic physics analysis
physicsAnalyzer = PhysicsAnalyzer()
physicsAnalyzer.runAnalysis(configPath)

# Step 3: Visualize results
if config.get('analysis', {}).get('showDashboard', True):
    visualizer = Visualizer()
    visualizer.createDashboard(physicsAnalyzer.results)
