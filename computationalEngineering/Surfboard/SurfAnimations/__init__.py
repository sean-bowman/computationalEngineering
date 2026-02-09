# -- SurfAnimations Package -- #

'''
Manim-based physics animation module.

Provides animated visualizations of wave theory, surfboard
hydrodynamics, performance comparison, and SPH water simulation.

API usage:
    from computationalEngineering.SurfAnimations import renderScene, SCENES
    renderScene('wave_intro', quality='high')

Sean Bowman [02/05/2026]
'''

from computationalEngineering.SurfAnimations.render import renderScene, renderAll, SCENES
from computationalEngineering.SurfAnimations.scenes.waveIntro import WavePhysicsIntro
from computationalEngineering.SurfAnimations.scenes.boardOnWave import BoardOnWaveSideView, BoardOnWavePerspective
from computationalEngineering.SurfAnimations.scenes.performanceComparison import PerformanceComparison
from computationalEngineering.SurfAnimations.scenes.sloshingTank import SloshingTankAnimation
