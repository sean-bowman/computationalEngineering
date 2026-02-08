# -- SurfAnimations Package -- #

'''
Manim-based physics animation module.

Provides animated visualizations of wave theory, surfboard
hydrodynamics, performance comparison, and SPH water simulation.

API usage:
    from SurfAnimations import renderScene, SCENES
    renderScene('wave_intro', quality='high')

Sean Bowman [02/05/2026]
'''

from SurfAnimations.render import renderScene, renderAll, SCENES
from SurfAnimations.scenes.waveIntro import WavePhysicsIntro
from SurfAnimations.scenes.boardOnWave import BoardOnWaveSideView, BoardOnWavePerspective
from SurfAnimations.scenes.performanceComparison import PerformanceComparison
from SurfAnimations.scenes.sloshingTank import SloshingTankAnimation
