# -- Surfboard Package -- #

'''
Surfboard design, physics, and visualization toolkit.

Parametric surfboard geometry generation, hydrodynamic analysis,
interactive 3D visualization, and physics animations.

Wildcard import exposes key classes:
    from computationalEngineering.Surfboard import *
    params = SurfboardParameters.shortboard()
    board = BoardGeometry(params)

Sub-modules:
    - SurfPhysics: Parametric geometry, wave theory, hydrodynamics
    - SurfboardGeometry: C# voxel-based geometry generation (PicoGK)
    - SurfAnimations: Manim animation scenes
    - SurfViewer: Three.js interactive web viewer

Sean Bowman [02/08/2026]
'''

# Parametric geometry
from computationalEngineering.Surfboard.SurfPhysics.geometry.parameters import SurfboardParameters
from computationalEngineering.Surfboard.SurfPhysics.geometry.board import BoardGeometry
from computationalEngineering.Surfboard.SurfPhysics.geometry.surfaceMeshGenerator import (
    SurfaceMeshGenerator,
    RESOLUTION_PRESETS,
)
from computationalEngineering.Surfboard.SurfPhysics.geometry.outline import Outline
from computationalEngineering.Surfboard.SurfPhysics.geometry.rocker import RockerProfile
from computationalEngineering.Surfboard.SurfPhysics.geometry.crossSection import CrossSection

# Physics analysis
from computationalEngineering.Surfboard.SurfPhysics.analyzer import PhysicsAnalyzer
from computationalEngineering.Surfboard.SurfPhysics.visualization.visualizer import Visualizer
from computationalEngineering.Surfboard.SurfPhysics.export.viewerExporter import ViewerExporter

# Validation
from computationalEngineering.Surfboard.SurfPhysics.validation.reverseEngineer import ReverseEngineer

# C# geometry wrapper
from computationalEngineering.Surfboard.SurfboardGeometry.geometryGenerator import GeometryGenerator

# Animations
from computationalEngineering.Surfboard.SurfAnimations.render import renderScene, renderAll, SCENES

__all__ = [
    # Geometry
    'SurfboardParameters',
    'BoardGeometry',
    'SurfaceMeshGenerator',
    'RESOLUTION_PRESETS',
    'Outline',
    'RockerProfile',
    'CrossSection',
    # Physics
    'PhysicsAnalyzer',
    'Visualizer',
    'ViewerExporter',
    # Validation
    'ReverseEngineer',
    # C# geometry
    'GeometryGenerator',
    # Animations
    'renderScene',
    'renderAll',
    'SCENES',
]
