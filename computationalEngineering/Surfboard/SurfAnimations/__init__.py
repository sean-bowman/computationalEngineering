# -- SurfAnimations Package -- #

'''
Manim-based physics animation module.

Provides animated visualizations of wave theory, surfboard
hydrodynamics, performance comparison, and SPH water simulation.

API usage:
    from computationalEngineering.Surfboard.SurfAnimations import renderScene, SCENES
    renderScene('wave_intro', quality='high')

Sean Bowman [02/05/2026]
'''

from computationalEngineering.Surfboard.SurfAnimations.render import renderScene, renderAll, SCENES
from computationalEngineering.Surfboard.SurfAnimations.scenes.waveIntro import WavePhysicsIntro
from computationalEngineering.Surfboard.SurfAnimations.scenes.boardOnWave import BoardOnWaveSideView, BoardOnWavePerspective
from computationalEngineering.Surfboard.SurfAnimations.scenes.performanceComparison import PerformanceComparison
from computationalEngineering.Surfboard.SurfAnimations.scenes.sloshingTank import SloshingTankAnimation
