# -- SurfPhysics Package -- #

'''
Surfboard physics simulation module.

Wave theory, hydrodynamic force analysis, and interactive visualization
for parametric surfboard designs.

Sean Bowman [02/03/2026]
'''

__version__ = '0.1.0'

from computationalEngineering.SurfPhysics.analyzer import PhysicsAnalyzer
from computationalEngineering.SurfPhysics.geometry.parameters import SurfboardParameters
from computationalEngineering.SurfPhysics.waves.waveConditions import WaveConditions
from computationalEngineering.SurfPhysics.visualization.visualizer import Visualizer
from computationalEngineering.SurfPhysics.export.viewerExporter import ViewerExporter

# Validation tools
from computationalEngineering.SurfPhysics.validation import (
    MeshComparisonAnalyzer,
    ComparisonResult,
    DistanceStatistics,
)
