# -- SPH Engine Package -- #

'''
Core Smoothed Particle Hydrodynamics (SPH) engine.

Provides kernel functions, particle systems, neighbor search,
boundary handling, time integration, and the WCSPH solver.

Sean Bowman [02/05/2026]
'''

from computationalEngineering.WaterSim.sph.protocols import SimulationConfig, SimulationState
from computationalEngineering.WaterSim.sph.kernels import CubicSplineKernel, WendlandC2Kernel, createKernel
from computationalEngineering.WaterSim.sph.waveMaker import WaveMakerConfig, PistonWaveMaker
from computationalEngineering.WaterSim.sph.breakingDetection import BreakingEvent, BreakingDetector
