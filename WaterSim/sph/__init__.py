# -- SPH Engine Package -- #

'''
Core Smoothed Particle Hydrodynamics (SPH) engine.

Provides kernel functions, particle systems, neighbor search,
boundary handling, time integration, and the WCSPH solver.

Sean Bowman [02/05/2026]
'''

from WaterSim.sph.protocols import SimulationConfig, SimulationState
from WaterSim.sph.kernels import CubicSplineKernel, WendlandC2Kernel, createKernel
from WaterSim.sph.waveMaker import WaveMakerConfig, PistonWaveMaker
from WaterSim.sph.breakingDetection import BreakingEvent, BreakingDetector
