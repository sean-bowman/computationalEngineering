
# -- Code Interface (Entry-Point) -- #

'''
Lightweight entry point for the Computational Engineering Toolkit.

Demonstrates the Surfboard module via wildcard import:
    1. Create board parameters from presets
    2. Generate a parametric surface mesh
    3. Compute volume and export STL
    4. Run reverse-engineering validation
    5. Visualize the 3D mesh in Plotly

Usage:
    python codeInterface.py

Sean Bowman [02/08/2026]
'''

import os

import plotly.graph_objects as go

from computationalEngineering.Surfboard import *

os.system('cls' if os.name == 'nt' else 'clear')

######################################################################
# -- Configuration -- #
######################################################################

boardPreset = 'shortboard'      # 'shortboard', 'longboard', or 'fish'
meshResolution = 'standard'     # 'draft', 'standard', or 'high'

######################################################################
# -- Generate Board -- #
######################################################################

presets = {
    'shortboard': SurfboardParameters.shortboard,
    'longboard': SurfboardParameters.longboard,
    'fish': SurfboardParameters.fish,
}

params = presets[boardPreset]()
board = BoardGeometry(params)

print('=' * 60)
print('  SURFBOARD MODULE SHOWCASE')
print('=' * 60)
params.printSummary()


######################################################################
# -- Surface Mesh -- #
######################################################################

print('\nGenerating surface mesh...')
resSettings = RESOLUTION_PRESETS[meshResolution]
gen = SurfaceMeshGenerator(params, **resSettings)
mesh = gen.generate()

meshVolL = gen.computeVolumeLiters()
paramVolL = board.computeVolume() * 1000

print(f'  Resolution:  {meshResolution} '
      f'({resSettings["nLongitudinal"]} x {resSettings["nLateral"]})')
print(f'  Vertices:    {len(mesh.vertices):,}')
print(f'  Faces:       {len(mesh.faces):,}')
print(f'  Watertight:  {mesh.is_watertight}')
print(f'  Mesh Volume: {meshVolL:.2f} L')
print(f'  Param Volume:{paramVolL:.2f} L (integration)')
print(f'  Difference:  {abs(meshVolL - paramVolL) / paramVolL * 100:.3f}%')

# Export STL
outputDir = 'computationalEngineering/Surfboard/SurfboardGeometry/Output'
os.makedirs(outputDir, exist_ok=True)
stlPath = os.path.join(outputDir, f'{boardPreset}_parametric.stl')
gen.exportStl(stlPath)
print(f'\n  Exported: {stlPath}')


######################################################################
# -- Reverse Engineering Validation -- #
######################################################################

print('\n' + '=' * 60)
print('  REVERSE ENGINEERING VALIDATION')
print('=' * 60)

re = ReverseEngineer(
    stlPath,
    expectedLengthMm=params.length,
    removeFins=False,
    autoTransform=False,
)

profiles = re.extractProfiles(nStations=300, nLateralSamples=50)
fittedParams = re.fitParameters(profiles)

fields = [
    ('Length', 'length'), ('Max Width', 'maxWidth'),
    ('Max Thickness', 'maxThickness'), ('Nose Rocker', 'noseRocker'),
    ('Tail Rocker', 'tailRocker'), ('Deck Crown', 'deckCrown'),
    ('Concave', 'bottomConcave'),
]

print(f'\n  {"Parameter":<16s}  {"Fitted":>8s}  {"Original":>8s}  {"Diff":>6s}')
print('  ' + '-' * 46)
for label, attr in fields:
    fitted = getattr(fittedParams, attr)
    original = getattr(params, attr)
    print(f'  {label:<16s}  {fitted:8.1f}  {original:8.1f}  {fitted - original:+6.1f}')

comparison = re.compareToParametric(fittedParams, profiles)
print(f'\n  Curve RMS Errors:')
print(f'    Outline:    {comparison["outlineRmsMm"]:.2f} mm')
print(f'    Rocker:     {comparison["rockerRmsMm"]:.2f} mm')
print(f'    Thickness:  {comparison["thicknessRmsMm"]:.2f} mm')


######################################################################
# -- 3D Visualization -- #
######################################################################

print('\nOpening 3D visualization...')

plotMesh = mesh
if len(mesh.faces) > 80000:
    try:
        plotMesh = mesh.simplify_quadric_decimation(80000)
    except Exception:
        plotMesh = mesh

verts = plotMesh.vertices
faces = plotMesh.faces

fig = go.Figure(data=[
    go.Mesh3d(
        x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        color='#4fc3f7',
        opacity=0.85,
        flatshading=False,
        lighting=dict(
            ambient=0.3, diffuse=0.7, specular=0.3,
            roughness=0.5, fresnel=0.2,
        ),
        lightposition=dict(x=1000, y=500, z=2000),
        hoverinfo='skip',
    ),
])

fig.update_layout(
    title=f'{boardPreset.title()} — {len(mesh.vertices):,} verts, {meshVolL:.1f}L',
    scene=dict(
        xaxis_title='X — Length (mm)',
        yaxis_title='Y — Width (mm)',
        zaxis_title='Z — Height (mm)',
        aspectmode='data',
        camera=dict(eye=dict(x=1.5, y=0.8, z=0.6), up=dict(x=0, y=0, z=1)),
    ),
    template='plotly_dark',
    height=700,
    width=1000,
)

fig.show()

print('\n' + '=' * 60)
print('  Done.')
print('=' * 60)
