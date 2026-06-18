# Surfboard Domain

**Computational surfboard design, physics simulation, 3D geometry generation, and reference validation.**

This domain covers the full pipeline from parametric board shape definition through physics analysis, animated visualizations, interactive 3D viewing, and reverse-engineering validation against real reference boards.

---

## Quick Start

### Run the reverse-engineering demo

```bash
# From the repo root — fits both Python parametric and C# voxel geometry
# to referenceThruster.stl, then launches the interactive comparison viewer
python codeInterface.py
```

### Generate a board with Python (parametric surface mesh)

```python
from computationalEngineering.Surfboard import *

params = SurfboardParameters.shortboard()  # or .longboard(), .fish()
gen    = SurfaceMeshGenerator.fromPreset(params, 'standard')
gen.exportStl('shortboard.stl')
```

### Generate a board with C# (PicoGK voxel geometry)

```bash
# From the repo root
dotnet run --project computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj -- --type shortboard

# With options
dotnet run --project computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj -- --type fish --voxel 0.25 --fins twin

# All options
dotnet run --project computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj -- --help
```

Output STL files are saved to `SurfboardGeometry/Output/`.

### Launch the interactive 3D web viewer

```bash
python -m http.server 8080
# Open: http://localhost:8080/computationalEngineering/Surfboard/SurfViewer/viewer.html
```

### Render physics animations

```bash
python computationalEngineering/Surfboard/SurfAnimations/render.py --scene wave_intro
python computationalEngineering/Surfboard/SurfAnimations/render.py --all
python computationalEngineering/Surfboard/SurfAnimations/render.py --list
```

---

## Module Guide

```
Surfboard/
├── SurfPhysics/          Python physics, geometry, validation, visualization
├── SurfboardGeometry/    C# / PicoGK voxel geometry generator
├── SurfAnimations/       Manim physics animation scenes
├── SurfViewer/           Three.js interactive 3D web viewer
├── configs/              JSON configuration files
└── documentation/        Technical writeups
```

### SurfPhysics

The core Python package. Import everything with:

```python
from computationalEngineering.Surfboard import *
# or
from computationalEngineering.Surfboard.SurfPhysics.geometry.parameters import SurfboardParameters
```

| Subpackage | Purpose |
|---|---|
| `geometry/` | Parametric board model: outline, rocker, cross-section, mesh generation |
| `waves/` | Linear wave theory: dispersion, phase/group velocity, pressure |
| `hydrodynamics/` | Drag, lift, buoyancy, planing force models |
| `visualization/` | Plotly analysis dashboards and the reference comparison viewer |
| `export/` | Viewer data and STL export utilities |
| `validation/` | Reverse engineering, mesh alignment, and deviation analysis |
| `optimization/` | Two-phase parameter optimizer (direct mapping + Nelder-Mead) |

### SurfboardGeometry (C#)

Voxel-based geometry using [PicoGK](https://github.com/leap71/PicoGK). Geometry is built by spatial painting — placing thousands of overlapping spheres along the board surface, merging into a voxel field, then extracting a triangle mesh via marching cubes.

CLI accepts `--config custom_params.json` for fully custom board shapes (see [Config Reference](#config-reference) below).

### SurfAnimations

Manim scenes for physics education. Available scenes:

| Scene | Description |
|---|---|
| `wave_intro` | Linear wave theory — propagating wave, particle orbits, velocity field |
| `board_side` | 2D side view with force vectors and velocity field |
| `board_perspective` | 3D perspective view with camera orbit |
| `performance` | Shortboard vs longboard vs fish performance comparison |
| `sloshing_tank` | SPH water simulation particle playback |

### SurfViewer

Three.js web viewer with dual geometry modes (STL mesh or parametric surface), solid/wireframe/glass render modes, physics overlays, board comparison, and deviation heatmap. Reads JSON data from `SurfViewer/data/` exported by the Python pipeline.

---

## Board Parameters

All parameters are in millimeters. The `SurfboardParameters` dataclass holds the full board definition:

| Parameter | Description | Default (shortboard) |
|---|---|---|
| `length` | Total board length | 1828 mm (6'0") |
| `maxWidth` | Maximum beam width | 495 mm (19.5") |
| `maxThickness` | Maximum thickness | 62 mm (2.44") |
| `noseWidth` | Width 305 mm from nose | 300 mm |
| `tailWidth` | Width 305 mm from tail | 380 mm |
| `widePointOffset` | Wide point from board center (negative = toward nose) | -25 mm |
| `tailTipHalfWidth` | Half-width at tail tip | 15 mm |
| `noseRocker` | Nose rocker height | 120 mm (4.7") |
| `tailRocker` | Tail rocker height | 40 mm (1.6") |
| `deckCrown` | Deck dome height at centerline | 8 mm |
| `bottomConcave` | Single concave channel depth | 2 mm |

### Presets

| Parameter | Shortboard | Longboard | Fish |
|---|---|---|---|
| Length | 1828 mm (6'0") | 2743 mm (9'0") | 1676 mm (5'6") |
| Width | 495 mm (19.5") | 570 mm (22.5") | 533 mm (21") |
| Thickness | 62 mm (2.44") | 75 mm (3.0") | 65 mm (2.56") |
| Nose Rocker | 120 mm (4.7") | 180 mm (7.0") | 80 mm (3.1") |
| Tail Rocker | 40 mm (1.6") | 25 mm (1.0") | 30 mm (1.2") |
| Default Fins | Thruster | Single | Twin |

### Fin Configurations

| Configuration | Count | Characteristics |
|---|---|---|
| `thruster` | 3 (center + 2 side) | Versatile, balanced drive and pivot |
| `twin` | 2 (side only) | Fast and loose, skatey feel |
| `quad` | 4 (2 front + 2 rear) | Speed and hold in large surf |
| `single` | 1 (center) | Smooth, flowing turns; stable |

---

## Validation and Reverse Engineering

The validation pipeline compares any generated board against a reference STL:

```python
from computationalEngineering.Surfboard.SurfPhysics.validation.reverseEngineer import ReverseEngineer

re       = ReverseEngineer('referenceThruster.stl', expectedLengthMm=1828.0)
profiles = re.extractProfiles(nStations=200)
params   = re.fitParameters(profiles)
errors   = re.compareToParametric(params, profiles)
print(f'Outline RMS: {errors["outlineRmsMm"]:.2f} mm')
```

For C# voxel optimization against a reference:

```python
from computationalEngineering.Surfboard.SurfPhysics.optimization.parameterOptimizer import ParameterOptimizer

optimizer = ParameterOptimizer('referenceThruster.stl', 'SurfboardGeometry/', 'output/')
result    = optimizer.runOptimization('shortboard', targetRmsMm=5.0)
```

The interactive comparison viewer (launched by `codeInterface.py`) overlays both geometry styles on the reference with per-vertex deviation coloring — green inside the reference, red outside.

---

## Config Reference

Configs in `configs/*.json` control the full pipeline via the Python entry point:

```json
{
  "board": {
    "type": "shortboard",
    "finConfiguration": "thruster",
    "foamType": "pu"
  },
  "geometry": {
    "generate": false,
    "voxelSize": 0.5,
    "outputDir": "..."
  },
  "validation": {
    "run": true,
    "referencePath": "...",
    "optimize": true,
    "targetRmsMm": 10.0,
    "maxPhase1Iterations": 3,
    "maxPhase2Iterations": 30
  },
  "waves": { "height": 1.5, "period": 10.0, "depth": 2.5 },
  "rider": { "mass": 75.0 },
  "analysis": { "run": true, "showDashboard": true }
}
```

Preset configs: `shortboard_default.json`, `longboard_small.json`, `fish_overhead.json`.

---

## Documentation

Detailed writeups in [`documentation/`](documentation/):

- **[API Reference](documentation/apiReference.md)** — All importable classes and functions
- **[SurfPhysics Overview](documentation/surfPhysicsOverview.md)** — Wave theory, buoyancy, and force balance
- **[Geometry Implementation](documentation/surfboardGeometryImplementation.md)** — Math details for outline, rocker, and cross-section curves
- **[Parametric Surface Mesh](documentation/parametricSurfaceMeshOverview.md)** — Python STL generation algorithm
- **[SurfboardGeometry Overview](documentation/surfboardGeometryOverview.md)** — Voxel-based C# generation
- **[SurfAnimations Overview](documentation/surfAnimationsOverview.md)** — Manim scene architecture
- **[Fin Dimensions Reference](documentation/finDimensionsReference.md)** — Reference fin measurements
- **[Surfboard Shaping References](documentation/surfboardShapingReferences.md)** — Industry shaping conventions

---

Sean Bowman
