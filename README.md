# Computational Engineering Portfolio

**Generation and analysis of physical models through computer algorithms.**

This repository demonstrates computational engineering -- the process of using software to define, generate, and analyze physical objects and their behavior. Rather than drawing shapes manually in CAD software, geometry is defined programmatically through mathematical functions and parameters, enabling precise control, reproducibility, and design exploration.

The library is organized by domain under `computationalEngineering/`. Each domain (e.g. `Surfboard/`) is a self-contained Python package with its own configs, documentation, and code. Generic tools like `PicoGK` and `WaterSim` live at the top level for reuse across domains.

---

## Quick Start

### Python — Parametric Surface Mesh

```bash
python codeInterface.py
```

The repo-root entry point uses a single wildcard import to access the entire Surfboard module. It generates a parametric surfboard mesh, runs reverse-engineering validation, and opens an interactive 3D Plotly visualization.

**Configuration** (edit variables at the top of `codeInterface.py`):

- `boardPreset` — `'shortboard'`, `'longboard'`, or `'fish'`
- `meshResolution` — `'draft'`, `'standard'`, or `'high'`

**Wildcard import** for programmatic use:

```python
from computationalEngineering.Surfboard import *
params = SurfboardParameters.shortboard()
board = BoardGeometry(params)
```

### C# — Voxel-Based Geometry (PicoGK)

```bash
dotnet run --project computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj
```

Generates voxel-based surfboard geometry using [PicoGK](https://github.com/leap71/PicoGK), LEAP 71's computational geometry kernel. Supports CLI options for board type, fin configuration, and voxel resolution.

```bash
# Generate a longboard with single fin
dotnet run --project computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj -- --type longboard

# Fish board at high resolution
dotnet run --project computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj -- --type fish --voxel 0.25

# Shortboard with quad fins
dotnet run --project computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj -- --fins quad

# Show all options
dotnet run --project computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj -- --help
```

Output STL files are saved to `computationalEngineering/Surfboard/SurfboardGeometry/Output/`.

### Web Viewer — Interactive 3D (Three.js)

```bash
# Serve from the repo root
python -m http.server 8080

# Open in browser:
#   http://localhost:8080/computationalEngineering/Surfboard/SurfViewer/viewer.html
```

The viewer reads JSON data exported by the Python analysis pipeline from `computationalEngineering/Surfboard/SurfViewer/data/`. Features include dual geometry modes (STL mesh or parametric surface), solid/wireframe/glass render modes, physics overlays (waterline, force vectors), and board comparison.

---

## Python API

All modules expose clean package-level imports for programmatic use:

```python
# Wildcard import — all key surfboard classes
from computationalEngineering.Surfboard import *

# Or import specific modules
from computationalEngineering.Surfboard.SurfPhysics.geometry.parameters import SurfboardParameters
from computationalEngineering.Surfboard.SurfPhysics.geometry.board import BoardGeometry
from computationalEngineering.Surfboard.SurfPhysics.geometry.surfaceMeshGenerator import SurfaceMeshGenerator
from computationalEngineering.Surfboard.SurfPhysics.validation.reverseEngineer import ReverseEngineer
from computationalEngineering.Surfboard.SurfboardGeometry.geometryGenerator import GeometryGenerator
from computationalEngineering.Surfboard.SurfAnimations.render import renderScene

# SPH water simulation (generic, top-level)
from computationalEngineering.WaterSim.runner import WaterSimRunner
```

Surfboard configs are in `computationalEngineering/Surfboard/configs/`. Full API documentation is in [computationalEngineering/Surfboard/documentation/apiReference.md](computationalEngineering/Surfboard/documentation/apiReference.md).

---

## Surfboard Geometry (PicoGK / C#)

The C# project generates parametric surfboard geometry using [PicoGK](https://github.com/leap71/PicoGK), LEAP 71's open-source computational geometry kernel. A surfboard's shape is defined by three mathematical profiles that are combined during voxel-based spatial painting:

```
  TOP VIEW (Outline)          SIDE VIEW (Rocker)         CROSS-SECTION
  ─────────────────           ──────────────────         ─────────────

       Nose                                                  Deck Crown
      ╱    ╲                                               ╱‾‾‾‾‾‾‾‾‾‾‾╲
    ╱        ╲        \ Nose Rocker                       │   Volume    │
   │  Wide    │        \                    /              ╲___________╱
   │  Point   │         \_______Flat_______/ Tail Rocker      Concave
    ╲        ╱
      ╲____╱
       Tail
```

### How It Works

1. **Outline** defines the board's planform shape (half-width at each position along the length)
2. **Rocker Profile** defines the bottom curvature from nose to tail
3. **Cross-Section** defines the deck crown, bottom concave, and thickness at each station

These three profiles are combined using **spatial painting** -- placing thousands of overlapping spheres along the board surface to build up the shape as a voxel field. The voxels are then converted to a triangle mesh and exported as an STL file.

### Board Presets

| Parameter   | Shortboard    | Longboard     | Fish          |
| ----------- | ------------- | ------------- | ------------- |
| Length      | 6'0" (1828mm) | 9'0" (2743mm) | 5'6" (1676mm) |
| Width       | 19.5" (495mm) | 22.5" (570mm) | 21" (533mm)   |
| Thickness   | 2.44" (62mm)  | 3.0" (75mm)   | 2.56" (65mm)  |
| Nose Rocker | 4.7" (120mm)  | 7.0" (180mm)  | 3.1" (80mm)   |
| Tail Rocker | 1.6" (40mm)   | 1.0" (25mm)   | 1.2" (30mm)   |

All dimensions are configurable through `SurfboardParameters`.

### Fin Configurations

| Configuration | Fins                 | Default For | Characteristics              |
| ------------- | -------------------- | ----------- | ---------------------------- |
| Thruster      | 3 (center + 2 sides) | Shortboard  | Versatile, good control      |
| Twin          | 2 (sides only)       | Fish        | Fast, loose, skatey feel     |
| Quad          | 4 (2 front + 2 rear) | --          | Speed and hold in large surf |
| Single        | 1 (center only)      | Longboard   | Smooth turns, stability      |

Fin configuration is independent of board type -- any combination can be used via the `--fins` argument.

---

## Prerequisites

- [.NET 10.0 SDK](https://dotnet.microsoft.com/download) or later
- [PicoGK Runtime](https://github.com/leap71/PicoGK) (included as git submodule)
- Python 3.x with `numpy`, `plotly`, `trimesh`
- Windows x64 (for PicoGK native libraries)

## Setup

```bash
# Clone with submodules
git clone --recurse-submodules https://github.com/YOUR_USERNAME/computationalEngineering.git
cd computationalEngineering

# If you already cloned without submodules
git submodule update --init --recursive

# Build the C# surfboard geometry project
dotnet build computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj

# Install Python dependencies
pip install numpy plotly trimesh
```

---

## Physics Animations (Manim)

Animated physics visualizations built with Manim Community Edition, showing wave theory, surfboard hydrodynamics, and multi-board performance comparison.

### Scenes

| Scene                 | Description                                                                                 |
| --------------------- | ------------------------------------------------------------------------------------------- |
| `wave_intro`        | Linear wave theory — propagating wave, labeled properties, particle orbits, velocity field |
| `board_side`        | Shortboard riding a wave (2D side view) with force vectors and velocity field               |
| `board_perspective` | Shortboard on 3D wave surface (2.5D perspective) with force arrows and camera orbit         |
| `performance`       | Three-board comparison — shortboard vs longboard vs fish with force balance and L/D chart  |
| `sloshing_tank`     | SPH water simulation — particle field playback with velocity color mapping                 |

### Usage

```bash
# Render a single scene (1080p)
python computationalEngineering/Surfboard/SurfAnimations/render.py --scene wave_intro

# Render all scenes
python computationalEngineering/Surfboard/SurfAnimations/render.py --all

# Fast preview (480p)
python computationalEngineering/Surfboard/SurfAnimations/render.py --scene board_side --quality low

# List available scenes
python computationalEngineering/Surfboard/SurfAnimations/render.py --list
```

Output MP4 files are saved to `computationalEngineering/Surfboard/SurfAnimations/media/`.

---

## Project Structure

```
repo_root/
├── codeInterface.py                       Python entry point
├── README.md
│
└── computationalEngineering/              Master code directory
    │
    ├── Surfboard/                         Surfboard design domain
    │   ├── SurfboardGeometry/             C# / PicoGK surfboard generator
    │   │   ├── Program.cs                 C# CLI entry point
    │   │   ├── Surfboard/                 Core geometry classes
    │   │   └── Utils/                     Physical constants
    │   ├── SurfPhysics/                   Python physics simulations
    │   │   ├── geometry/                  Parametric board geometry
    │   │   ├── waves/                     Linear wave theory
    │   │   ├── hydrodynamics/             Force models
    │   │   ├── visualization/             Plotly dashboards
    │   │   ├── export/                    Viewer data export
    │   │   ├── validation/                Reverse engineering & mesh comparison
    │   │   └── optimization/              Parameter optimization
    │   ├── SurfAnimations/                Manim physics animations
    │   ├── SurfViewer/                    Three.js interactive 3D viewer
    │   ├── configs/                       Surfboard JSON configs
    │   └── documentation/                 Surfboard technical writeups
    │
    ├── WaterSim/                          SPH water simulation (generic)
    │   ├── sph/                           WCSPH solver engine
    │   ├── scenarios/                     Pre-configured simulations
    │   ├── export/                        Frame data export
    │   ├── configs/                       Water sim JSON configs
    │   └── documentation/                 Water sim technical writeups
    │
    └── PicoGK/                            LEAP 71 geometry kernel (submodule)
```

## Voxel-Based Geometry

This project uses **voxel-based computational geometry** rather than traditional boundary representation (B-rep) CAD. In this approach:

- Geometry is represented as a 3D grid of "on/off" cells (like 3D pixels)
- Complex organic shapes are built by "painting" overlapping spheres into the grid
- Boolean operations (union, subtract, intersect) are trivial and robust
- The result is directly compatible with 3D printing workflows

The voxel size parameter controls resolution: 0.5mm gives good detail for most purposes, while 0.25mm produces smoother surfaces at the cost of computation time.

## Surfboard Design Theory

### Outline

The planform shape uses cosine interpolation between key width stations (nose width, wide point, tail width). This produces smooth, organic curves without the oscillation problems of polynomial interpolation.

### Rocker

The bottom curvature uses power-law functions with different exponents for the nose and tail regions. The nose uses exponent 2.5 (steep kick, flat center) while the tail uses exponent 2.0 (gradual curve). This matches how real boards are shaped -- more aggressive nose entry with a smoother tail transition.

### Cross-Section

At each station along the board, the cross-section defines:

- **Deck crown**: A cosine dome that adds volume at the center while keeping rails thin
- **Bottom concave**: A channel along the centerline that accelerates water flow for speed
- **Thickness distribution**: Follows a cosine envelope that peaks at 40% from the nose

---

## Future Work

- [X] Document all work thoroughly with informational and conversational overviews of all modules
- [X] Prioritize API calls over CLI calls for code library interfacing
- [ ] Update surfboard and fin shapes against reference images and STLs for maximum realism
- [X] Investigate PicoGK geometry generation methodologies: is there a better way than packed spheres?
- [X] Parallel develop surface-mesh based geometry generation through Python
- [X] Consolidate the references and documentation folders
- [ ] Extend WaterSim to 3D with wave-maker boundaries for breaking wave modeling
- [ ] Add new geometry domains (boat hull, etc.) as sibling packages under `computationalEngineering/`
- [ ] Port pytthon utilities and coding standards from existing projects and cross reference progress to github pages portfolio site

## Technology

| Component           | Language              | Purpose                         |
| ------------------- | --------------------- | ------------------------------- |
| Surfboard Geometry  | C# / PicoGK           | Voxel-based 3D model generation |
| Physics Simulations | Python / NumPy        | Wave theory & hydrodynamics     |
| Water Simulation    | Python / NumPy        | SPH sloshing tank simulation    |
| 3D Visualization    | JavaScript / Three.js | Interactive web viewer          |
| Animations          | Python / Manim        | Physics education videos        |

## Documentation

### Surfboard

Detailed documentation is in `computationalEngineering/Surfboard/documentation/`:

- **[API Reference](computationalEngineering/Surfboard/documentation/apiReference.md)** — Quick reference for all importable classes and functions
- **[SurfPhysics Overview](computationalEngineering/Surfboard/documentation/surfPhysicsOverview.md)** — Wave theory, buoyancy, and force balance physics
- **[SurfAnimations Overview](computationalEngineering/Surfboard/documentation/surfAnimationsOverview.md)** — Manim scene architecture and components
- **[SurfboardGeometry Overview](computationalEngineering/Surfboard/documentation/surfboardGeometryOverview.md)** — Voxel-based geometry generation
- **[Geometry Implementation](computationalEngineering/Surfboard/documentation/surfboardGeometryImplementation.md)** — Mathematical details of outline, rocker, and cross-section
- **[Parametric Surface Mesh](computationalEngineering/Surfboard/documentation/parametricSurfaceMeshOverview.md)** — Python parametric surface mesh generation
- **[Fin Dimensions Reference](computationalEngineering/Surfboard/documentation/finDimensionsReference.md)** — Reference fin measurements
- **[Surfboard Shaping References](computationalEngineering/Surfboard/documentation/surfboardShapingReferences.md)** — Industry shaping conventions

### Water Simulation

- **[WaterSim Overview](computationalEngineering/WaterSim/documentation/waterSimOverview.md)** — SPH simulation theory and implementation

---

## References

- [PicoGK by LEAP 71](https://github.com/leap71/PicoGK) -- Computational geometry kernel
- Paine, Frank. *The Science of Surfing*. Cambridge University Press.
- Dean & Dalrymple. *Water Wave Mechanics for Engineers and Scientists*.
- Finney, Ben & Houston, James. *Surfing: A History of the Ancient Hawaiian Sport*.

---

Sean Bowman
