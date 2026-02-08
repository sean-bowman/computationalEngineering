
# -- Code Interface (Entry-Point) -- #

'''
Main interface for the Computational Engineering Toolkit.

Demonstrates the parametric surface mesh generation pipeline:
    1. Generate a surfboard mesh from parametric cross-sections
    2. Visualize the 3D mesh, outline, rocker, and cross-sections in Plotly
    3. Reverse-engineer the generated STL to extract and validate profiles
    4. Compare fitted parameters against the originals

Run directly:
    python codeInterface.py

Sean Bowman [02/07/2026]
'''

import os
import sys
import json

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from computationalEngineering.SurfPhysics.geometry.parameters import SurfboardParameters
from computationalEngineering.SurfPhysics.geometry.board import BoardGeometry
from computationalEngineering.SurfPhysics.geometry.surfaceMeshGenerator import SurfaceMeshGenerator, RESOLUTION_PRESETS
from computationalEngineering.SurfPhysics.validation.reverseEngineer import ReverseEngineer
from computationalEngineering.SurfPhysics.visualization import theme

# Clear terminal
os.system('cls' if os.name == 'nt' else 'clear')


######################################################################
# -- Configuration -- #
######################################################################

# Board to generate — 'shortboard', 'longboard', or 'fish'
boardPreset = 'shortboard'

# Mesh resolution — 'draft', 'standard', or 'high'
meshResolution = 'standard'

# Output directory for generated STL files
outputDir = 'computationalEngineering/SurfboardGeometry/Output'

# Whether to run reverse-engineering validation
runReverseEngineering = True


######################################################################
# -- Board Parameter Selection -- #
######################################################################

presetFactories = {
    'shortboard': SurfboardParameters.shortboard,
    'longboard': SurfboardParameters.longboard,
    'fish': SurfboardParameters.fish,
}

params = presetFactories[boardPreset]()
board = BoardGeometry(params)

print('=' * 60)
print(f'  PARAMETRIC SURFACE MESH GENERATOR')
print('=' * 60)
params.printSummary()


######################################################################
# -- Mesh Generation -- #
######################################################################

print('\nGenerating surface mesh...')
resSettings = RESOLUTION_PRESETS[meshResolution]
gen = SurfaceMeshGenerator(params, **resSettings)
mesh = gen.generate()

print(f'  Resolution:  {meshResolution} '
      f'({resSettings["nLongitudinal"]} x {resSettings["nLateral"]})')
print(f'  Vertices:    {len(mesh.vertices):,}')
print(f'  Faces:       {len(mesh.faces):,}')
print(f'  Watertight:  {mesh.is_watertight}')
print(f'  Volume:      {gen.computeVolumeLiters():.2f} L')

# Compare against parametric volume integration
paramVolL = board.computeVolume() * 1000  # m^3 -> liters
meshVolL = gen.computeVolumeLiters()
volDiffPct = abs(meshVolL - paramVolL) / paramVolL * 100
print(f'  Param Volume:{paramVolL:.2f} L (integration)')
print(f'  Difference:  {volDiffPct:.3f}%')

# Export STL
os.makedirs(outputDir, exist_ok=True)
stlPath = os.path.join(outputDir, f'{boardPreset}_parametric.stl')
gen.exportStl(stlPath)
print(f'\n  Exported: {stlPath}')


######################################################################
# -- 3D Mesh Visualization -- #
######################################################################

print('\nBuilding 3D visualization...')

# Downsample mesh for Plotly performance if face count is high
plotMesh = mesh
if len(mesh.faces) > 80000:
    try:
        plotMesh = mesh.simplify_quadric_decimation(80000)
    except Exception:
        plotMesh = mesh

verts = plotMesh.vertices
faces = plotMesh.faces

# Build a 4-panel figure:
# Row 1: 3D mesh view (spans both columns)
# Row 2: Outline (left), Rocker Profile (right)
# Row 3: Cross-sections at 6 stations

fig = make_subplots(
    rows=3, cols=2,
    specs=[
        [{'type': 'scene', 'colspan': 2}, None],
        [{'type': 'xy'}, {'type': 'xy'}],
        [{'type': 'xy', 'colspan': 2}, None],
    ],
    row_heights=[0.50, 0.25, 0.25],
    subplot_titles=[
        '3D Surface Mesh',
        'Planform Outline', 'Rocker Profile',
        'Cross-Sections (6 Stations)',
    ],
    vertical_spacing=0.08,
    horizontal_spacing=0.08,
)


# -- Panel 1: 3D Mesh (go.Mesh3d) --

fig.add_trace(
    go.Mesh3d(
        x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
        i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
        color=theme.BLUE,
        opacity=0.85,
        flatshading=False,
        lighting=dict(
            ambient=0.3, diffuse=0.7, specular=0.3,
            roughness=0.5, fresnel=0.2,
        ),
        lightposition=dict(x=1000, y=500, z=2000),
        showscale=False,
        hoverinfo='skip',
    ),
    row=1, col=1,
)

# Configure 3D scene
fig.update_scenes(
    dict(
        xaxis=dict(title='X — Length (mm)', showbackground=False),
        yaxis=dict(title='Y — Width (mm)', showbackground=False),
        zaxis=dict(title='Z — Height (mm)', showbackground=False),
        aspectmode='data',
        camera=dict(
            eye=dict(x=1.5, y=0.8, z=0.6),
            up=dict(x=0, y=0, z=1),
        ),
    ),
    row=1, col=1,
)


# -- Panel 2: Planform Outline --

nPts = 200
tValues = np.linspace(0, 1, nPts)
halfWidths = np.array([board.getHalfWidthM(t) * 1000 for t in tValues])
xPositions = tValues * params.length

# Right rail
fig.add_trace(go.Scatter(
    x=xPositions, y=halfWidths, mode='lines',
    line=dict(color=theme.BLUE, width=2),
    showlegend=False,
), row=2, col=1)

# Left rail
fig.add_trace(go.Scatter(
    x=xPositions, y=-halfWidths, mode='lines',
    line=dict(color=theme.BLUE, width=2),
    showlegend=False,
), row=2, col=1)

# Centerline
fig.add_trace(go.Scatter(
    x=xPositions, y=np.zeros(nPts), mode='lines',
    line=dict(color=theme.REFERENCE_LINE, dash='dot', width=0.5),
    showlegend=False,
), row=2, col=1)

fig.update_xaxes(title_text='Length from Nose (mm)', row=2, col=1)
fig.update_yaxes(title_text='Half-Width (mm)', row=2, col=1,
                 scaleanchor='x2', scaleratio=1)


# -- Panel 3: Rocker Profile --

rockerHeights = np.array([board.getRockerHeightM(t) * 1000 for t in tValues])

fig.add_trace(go.Scatter(
    x=xPositions, y=rockerHeights, mode='lines',
    line=dict(color=theme.RED, width=2),
    showlegend=False,
), row=2, col=2)

# Flat reference
fig.add_trace(go.Scatter(
    x=xPositions, y=np.zeros(nPts), mode='lines',
    line=dict(color=theme.REFERENCE_LINE, dash='dot', width=0.5),
    showlegend=False,
), row=2, col=2)

fig.update_xaxes(title_text='Length from Nose (mm)', row=2, col=2)
fig.update_yaxes(title_text='Rocker Height (mm)', row=2, col=2)


# -- Panel 4: Cross-Sections at 6 Stations --

stations = [0.10, 0.25, 0.40, 0.60, 0.80, 0.95]
stationColors = theme.STATION_COLORS

for idx, t in enumerate(stations):
    hw = board.getHalfWidthM(t) * 1000
    if hw < 1.0:
        continue

    nLat = 60
    lfs = np.linspace(0.0, 1.0, nLat)
    deckZ = np.array([board.getDeckHeightM(t, lf) * 1000 for lf in lfs])
    bottomZ = np.array([board.getBottomHeightM(t, lf) * 1000 for lf in lfs])
    yCoords = lfs * hw

    color = stationColors[idx % len(stationColors)]
    label = f't={t:.2f} ({t * params.length:.0f}mm)'

    # Deck (solid)
    fig.add_trace(go.Scatter(
        x=yCoords, y=deckZ, mode='lines',
        name=label,
        line=dict(color=color, width=2),
        legendgroup=label,
    ), row=3, col=1)

    # Bottom (dashed)
    fig.add_trace(go.Scatter(
        x=yCoords, y=bottomZ, mode='lines',
        line=dict(color=color, width=2, dash='dash'),
        showlegend=False,
        legendgroup=label,
    ), row=3, col=1)

fig.update_xaxes(title_text='Distance from Centerline (mm)', row=3, col=1)
fig.update_yaxes(title_text='Height (mm)', row=3, col=1)


# -- Layout --

fig.update_layout(
    title=dict(
        text=f'{boardPreset.title()} — Parametric Surface Mesh '
             f'({meshResolution}, {len(mesh.vertices):,} verts, '
             f'{gen.computeVolumeLiters():.1f}L)',
        font=dict(size=16),
    ),
    template=theme.TEMPLATE,
    height=1200,
    width=1100,
    legend=dict(
        orientation='h', yanchor='bottom', y=-0.02,
        xanchor='center', x=0.5,
    ),
)

fig.show()


######################################################################
# -- Reverse Engineering Validation -- #
######################################################################

if runReverseEngineering:
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

    # Print parameter comparison table
    print(f'\n  {"Parameter":<16s}  {"Fitted":>8s}  {"Original":>8s}  {"Diff":>6s}')
    print('  ' + '-' * 46)
    fields = [
        ('Length', 'length'), ('Max Width', 'maxWidth'),
        ('Max Thickness', 'maxThickness'), ('Nose Width', 'noseWidth'),
        ('Tail Width', 'tailWidth'), ('WP Offset', 'widePointOffset'),
        ('Nose Rocker', 'noseRocker'), ('Tail Rocker', 'tailRocker'),
        ('Deck Crown', 'deckCrown'), ('Concave', 'bottomConcave'),
    ]
    for label, attr in fields:
        fitted = getattr(fittedParams, attr)
        original = getattr(params, attr)
        diff = fitted - original
        print(f'  {label:<16s}  {fitted:8.1f}  {original:8.1f}  {diff:+6.1f}')

    # Curve-level comparison
    comparison = re.compareToParametric(fittedParams, profiles)
    print(f'\n  Curve RMS Errors:')
    print(f'    Outline:    {comparison["outlineRmsMm"]:.2f} mm')
    print(f'    Rocker:     {comparison["rockerRmsMm"]:.2f} mm')
    print(f'    Thickness:  {comparison["thicknessRmsMm"]:.2f} mm')
    print(f'    Deck Prof:  {comparison["deckProfileRmsMm"]:.2f} mm')
    print(f'    Bottom Prof:{comparison["bottomProfileRmsMm"]:.2f} mm')

    # Build comparison plots: extracted vs. parametric curves
    compFig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[
            'Outline: Extracted vs. Parametric',
            'Rocker: Extracted vs. Parametric',
            'Thickness: Extracted vs. Parametric',
            'Cross-Section at Thickest Station',
        ],
        vertical_spacing=0.12,
        horizontal_spacing=0.10,
    )

    from computationalEngineering.SurfPhysics.geometry.outline import Outline
    from computationalEngineering.SurfPhysics.geometry.rocker import RockerProfile
    from computationalEngineering.SurfPhysics.geometry.crossSection import CrossSection

    fitOutline = Outline(fittedParams)
    fitRocker = RockerProfile(fittedParams)
    fitCrossSection = CrossSection(fittedParams)

    profileX = profiles.tValues * profiles.boardLengthMm

    # Outline comparison
    paramHW = np.array([fitOutline.getHalfWidth(float(t)) for t in profiles.tValues])
    compFig.add_trace(go.Scatter(
        x=profileX, y=profiles.halfWidths, mode='lines',
        name='Extracted', line=dict(color=theme.RED, width=2),
    ), row=1, col=1)
    compFig.add_trace(go.Scatter(
        x=profileX, y=paramHW, mode='lines',
        name='Parametric Fit', line=dict(color=theme.BLUE, width=2, dash='dash'),
    ), row=1, col=1)
    compFig.update_xaxes(title_text='Length (mm)', row=1, col=1)
    compFig.update_yaxes(title_text='Half-Width (mm)', row=1, col=1)

    # Rocker comparison
    paramRZ = np.array([fitRocker.getRockerHeight(float(t)) for t in profiles.tValues])
    compFig.add_trace(go.Scatter(
        x=profileX, y=profiles.rockerHeights, mode='lines',
        name='Extracted', line=dict(color=theme.RED, width=2),
        showlegend=False,
    ), row=1, col=2)
    compFig.add_trace(go.Scatter(
        x=profileX, y=paramRZ, mode='lines',
        name='Parametric Fit', line=dict(color=theme.BLUE, width=2, dash='dash'),
        showlegend=False,
    ), row=1, col=2)
    compFig.update_xaxes(title_text='Length (mm)', row=1, col=2)
    compFig.update_yaxes(title_text='Rocker Height (mm)', row=1, col=2)

    # Thickness comparison
    paramThick = np.array([fitCrossSection.getLocalThickness(float(t))
                           for t in profiles.tValues])
    compFig.add_trace(go.Scatter(
        x=profileX, y=profiles.thicknesses, mode='lines',
        name='Extracted', line=dict(color=theme.RED, width=2),
        showlegend=False,
    ), row=2, col=1)
    compFig.add_trace(go.Scatter(
        x=profileX, y=paramThick, mode='lines',
        name='Parametric Fit', line=dict(color=theme.BLUE, width=2, dash='dash'),
        showlegend=False,
    ), row=2, col=1)
    compFig.update_xaxes(title_text='Length (mm)', row=2, col=1)
    compFig.update_yaxes(title_text='Thickness (mm)', row=2, col=1)

    # Cross-section at thickest station
    thickestIdx = int(np.argmax(profiles.thicknesses))
    tThickest = float(profiles.tValues[thickestIdx])
    lfs = profiles.lateralFractions

    compFig.add_trace(go.Scatter(
        x=lfs, y=profiles.deckProfiles[thickestIdx], mode='lines',
        name='Deck (Extracted)', line=dict(color=theme.RED, width=2),
        showlegend=False,
    ), row=2, col=2)
    compFig.add_trace(go.Scatter(
        x=lfs, y=profiles.bottomProfiles[thickestIdx], mode='lines',
        name='Bottom (Extracted)', line=dict(color=theme.ORANGE, width=2),
        showlegend=False,
    ), row=2, col=2)

    paramDeck = np.array([fitCrossSection.getDeckHeight(tThickest, float(lf)) for lf in lfs])
    paramBottom = np.array([fitCrossSection.getBottomHeight(tThickest, float(lf)) for lf in lfs])
    compFig.add_trace(go.Scatter(
        x=lfs, y=paramDeck, mode='lines',
        name='Deck (Param)', line=dict(color=theme.BLUE, width=2, dash='dash'),
        showlegend=False,
    ), row=2, col=2)
    compFig.add_trace(go.Scatter(
        x=lfs, y=paramBottom, mode='lines',
        name='Bottom (Param)', line=dict(color=theme.CYAN, width=2, dash='dash'),
        showlegend=False,
    ), row=2, col=2)
    compFig.update_xaxes(title_text='Lateral Fraction (0=center, 1=rail)', row=2, col=2)
    compFig.update_yaxes(title_text='Height (mm)', row=2, col=2)

    compFig.update_layout(
        title=f'Reverse Engineering Validation — {boardPreset.title()}',
        template=theme.TEMPLATE,
        height=800,
        width=1100,
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5),
    )

    compFig.show()

print('\n' + '=' * 60)
print('  Done.')
print('=' * 60)
