# -- SurfPhysics Package -- #

'''
Surfboard physics simulation module.

Wave theory, hydrodynamic force analysis, and interactive visualization
for parametric surfboard designs.

Sean Bowman [02/03/2026]
'''

__version__ = '0.1.0'

from computationalEngineering.Surfboard.SurfPhysics.analyzer import PhysicsAnalyzer
from computationalEngineering.Surfboard.SurfPhysics.geometry.parameters import SurfboardParameters
from computationalEngineering.Surfboard.SurfPhysics.waves.waveConditions import WaveConditions
from computationalEngineering.Surfboard.SurfPhysics.visualization.visualizer import Visualizer
from computationalEngineering.Surfboard.SurfPhysics.export.viewerExporter import ViewerExporter

# Validation tools
from computationalEngineering.Surfboard.SurfPhysics.validation import (
    MeshComparisonAnalyzer,
    ComparisonResult,
    DistanceStatistics,
)
