# -- Animation Scenes Subpackage -- #

'''
Manim Scene classes for physics animations.

Each scene is a self-contained Manim animation that can be
rendered via the renderScene() API or the CLI.
'''

from SurfAnimations.scenes.waveIntro import WavePhysicsIntro
from SurfAnimations.scenes.boardOnWave import BoardOnWaveSideView, BoardOnWavePerspective
from SurfAnimations.scenes.performanceComparison import PerformanceComparison
from SurfAnimations.scenes.sloshingTank import SloshingTankAnimation
