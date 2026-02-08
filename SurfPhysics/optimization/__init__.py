# -- Optimization Module -- #

'''
Parameter optimization for matching generated geometry to reference meshes.

Provides two-phase optimization:
1. Direct parameter mapping from bounding box and regional deviations
2. Nelder-Mead simplex fine-tuning for sub-millimeter accuracy

Sean Bowman [02/05/2026]
'''

from SurfPhysics.optimization.parameterOptimizer import (
    ParameterOptimizer,
    OptimizationResult,
    DEFAULT_PARAMS,
)
from SurfPhysics.optimization.directMapper import (
    DirectParameterMapper,
    PARAMETER_BOUNDS,
)
from SurfPhysics.optimization.coordinateTransformer import CoordinateTransformer


__all__ = [
    'ParameterOptimizer',
    'OptimizationResult',
    'DEFAULT_PARAMS',
    'DirectParameterMapper',
    'PARAMETER_BOUNDS',
    'CoordinateTransformer',
]
