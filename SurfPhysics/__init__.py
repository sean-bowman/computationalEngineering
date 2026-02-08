# -- SurfPhysics Package -- #

'''
Surfboard physics simulation module.

Wave theory, hydrodynamic force analysis, and interactive visualization
for parametric surfboard designs.

Sean Bowman [02/03/2026]
'''

__version__ = '0.1.0'

from SurfPhysics.analyzer import PhysicsAnalyzer
from SurfPhysics.geometry.parameters import SurfboardParameters
from SurfPhysics.waves.waveConditions import WaveConditions
from SurfPhysics.visualization.visualizer import Visualizer
from SurfPhysics.export.viewerExporter import ViewerExporter

# Validation tools
from SurfPhysics.validation import (
    MeshComparisonAnalyzer,
    ComparisonResult,
    DistanceStatistics,
)
