# -- WaterSim Package -- #

'''
Water simulation module using Smoothed Particle Hydrodynamics (SPH).

Sloshing tank simulations, dam breaks, and free-surface flow
visualization for computational engineering education.

Sean Bowman [02/05/2026]
'''

__version__ = '0.1.0'

from computationalEngineering.WaterSim.runner import WaterSimRunner
from computationalEngineering.WaterSim.scenarios.sloshingTank import SloshingTankConfig
from computationalEngineering.WaterSim.export.frameExporter import FrameExporter
