# -- Animation Scenes Subpackage -- #

'''
Manim Scene classes for physics animations.

Each scene is a self-contained Manim animation that can be
rendered via the renderScene() API or the CLI.
'''

from computationalEngineering.Surfboard.SurfAnimations.scenes.waveIntro import WavePhysicsIntro
from computationalEngineering.Surfboard.SurfAnimations.scenes.boardOnWave import BoardOnWaveSideView, BoardOnWavePerspective
from computationalEngineering.Surfboard.SurfAnimations.scenes.performanceComparison import PerformanceComparison
from computationalEngineering.Surfboard.SurfAnimations.scenes.sloshingTank import SloshingTankAnimation
