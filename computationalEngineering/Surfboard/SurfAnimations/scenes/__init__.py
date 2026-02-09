# -- Animation Scenes Subpackage -- #

'''
Manim Scene classes for physics animations.

Each scene is a self-contained Manim animation that can be
rendered via the renderScene() API or the CLI.
'''

from computationalEngineering.SurfAnimations.scenes.waveIntro import WavePhysicsIntro
from computationalEngineering.SurfAnimations.scenes.boardOnWave import BoardOnWaveSideView, BoardOnWavePerspective
from computationalEngineering.SurfAnimations.scenes.performanceComparison import PerformanceComparison
from computationalEngineering.SurfAnimations.scenes.sloshingTank import SloshingTankAnimation
