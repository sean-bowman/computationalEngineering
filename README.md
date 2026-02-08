# Computational Engineering Portfolio

**Generation and analysis of physical models through computer algorithms.**

This repository demonstrates computational engineering -- the process of using software to define, generate, and analyze physical objects and their behavior. Rather than drawing shapes manually in CAD software, geometry is defined programmatically through mathematical functions and parameters, enabling precise control, reproducibility, and design exploration.

---

## Python API

All modules expose clean package-level imports for programmatic use. The `codeInterface.py` file demonstrates the full pipeline:

```python
# Core analysis pipeline
from SurfPhysics import PhysicsAnalyzer, SurfboardParameters, Visualizer, ViewerExporter
analyzer = PhysicsAnalyzer()
results = analyzer.runAnalysis('configs/shortboard_default.json')

# SPH water simulation
from WaterSim import WaterSimRunner, SloshingTankConfig
runner = WaterSimRunner()
results = runner.runSloshing(SloshingTankConfig.small2D())

# Manim animations
from SurfAnimations import renderScene
renderScene('wave_intro', quality='high')

# C# geometry generation (via subprocess)
from SurfboardGeometry import GeometryGenerator
generator = GeometryGenerator()
stlPath = generator.generateFromConfig('configs/shortboard_default.json')
```

Each step is config-flag-driven -- see `configs/shortboard_default.json` for all options. Full API documentation is in `documentation/api-reference.md`.

---

## Surfboard Geometry (PicoGK / C#)

The primary project generates parametric surfboard geometry using [PicoGK](https://github.com/leap71/PicoGK), LEAP 71's open-source computational geometry kernel. A surfboard's shape is defined by three mathematical profiles that are combined during voxel-based spatial painting:

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

Fin configuration is independent of board type -- any combination can be used via the `--fins` CLI argument.

---

## Prerequisites

- [.NET 10.0 SDK](https://dotnet.microsoft.com/download) or later
- [PicoGK Runtime](https://github.com/leap71/PicoGK) (included as git submodule)
- Windows x64 (for PicoGK native libraries)

## Setup

```bash
# Clone with submodules
git clone --recurse-submodules https://github.com/YOUR_USERNAME/computationalEngineering.git
cd computationalEngineering

# If you already cloned without submodules
git submodule update --init --recursive

# Build
dotnet build SurfboardGeometry/
```

## Usage

```bash
# Generate a default shortboard (6'0", thruster fins, 0.5mm voxel resolution)
dotnet run --project SurfboardGeometry/

# Generate a longboard (defaults to single fin)
dotnet run --project SurfboardGeometry/ -- --type longboard

# Generate a fish board (defaults to twin fins)
dotnet run --project SurfboardGeometry/ -- --type fish

# Generate a fish board at high resolution
dotnet run --project SurfboardGeometry/ -- --type fish --voxel 0.25

# Shortboard with quad fins
dotnet run --project SurfboardGeometry/ -- --fins quad

# Fast preview (low resolution)
dotnet run --project SurfboardGeometry/ -- --type shortboard --voxel 1

# Generate all 4 fin configurations as separate STL files
dotnet run --project SurfboardGeometry/ -- --all-fins --voxel 2

# Show help
dotnet run --project SurfboardGeometry/ -- --help
```

Output STL files are saved to `SurfboardGeometry/Output/`.

---

## Interactive 3D Viewer (Three.js)

The project includes an interactive Three.js-based 3D viewer for inspecting surfboard geometry and physics data in the browser.

### Features

- **Dual geometry modes**: Load voxel-generated STL meshes or parametric surface reconstructions
- **Render modes**: Solid (fiberglass finish), wireframe, and transparent glass
- **Physics overlays**: Waterline plane, force vectors (weight, buoyancy, drag, lift), draft indicator
- **Board comparison**: Side-by-side or overlay comparison of up to 3 board types

### Usage

```bash
# 1. Export viewer data (run the analysis pipeline with viewer export enabled)
python codeInterface.py

# 2. Serve the viewer locally (Three.js requires HTTP, not file://)
python -m http.server 8080 --directory SurfViewer

# 3. Open in browser
#    http://localhost:8080/viewer.html
```

The viewer reads JSON data exported by the Python analysis pipeline from `SurfViewer/data/`. Each board type produces a `boardData_<type>.json` file containing parametric surface geometry, physics results, and force vector data.

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
python SurfAnimations/render.py --scene wave_intro

# Render all scenes
python SurfAnimations/render.py --all

# Fast preview (480p)
python SurfAnimations/render.py --scene board_side --quality low

# List available scenes
python SurfAnimations/render.py --list
```

Output MP4 files are saved to `SurfAnimations/media/`.

---

## Project Structure

```
computationalEngineering/
├── SurfboardGeometry/                     C# / PicoGK surfboard generator
│   ├── Program.cs                         CLI entry point
│   ├── Surfboard/
│   │   ├── SurfboardParameters.cs         Parametric dimensions & presets
│   │   ├── SurfboardBody.cs               Spatial painting geometry generation
│   │   ├── Outline.cs                     Planform shape (top view)
│   │   ├── RockerProfile.cs               Bottom curvature (side view)
│   │   ├── CrossSection.cs               Deck/bottom profile (cross-section)
│   │   ├── FinConfiguration.cs            Fin setup enum (thruster/twin/quad/single)
│   │   └── FinSystem.cs                   Fin geometry generation
│   ├── Utils/
│   │   └── Constants.cs                   Physical constants
│   └── Output/                            Generated STL files
├── SurfViewer/                            Three.js interactive 3D viewer
│   ├── viewer.html                        Entry point
│   ├── js/                                ES modules (app, scene, renderer, UI)
│   ├── lib/                               Local Three.js dependencies
│   └── data/                              Exported JSON (gitignored)
├── SurfAnimations/                        Manim physics animations
│   ├── scenes/                            Animation scene classes
│   ├── components/                        Reusable Manim mobjects
│   ├── utils/                             Theme colors
│   ├── render.py                          CLI render script
│   └── media/                             Rendered output (gitignored)
├── SurfPhysics/                           Python physics simulations
│   ├── geometry/                          Parametric board geometry
│   ├── waves/                             Linear wave theory
│   ├── hydrodynamics/                     Force models (buoyancy, planing, drag)
│   ├── visualization/                     Plotly dashboards
│   └── export/                            Viewer data export pipeline
├── WaterSim/                              SPH water simulation
│   ├── sph/                               WCSPH solver engine
│   ├── scenarios/                         Pre-configured simulations
│   └── export/                            Frame data export
├── PicoGK/                                LEAP 71 geometry kernel (submodule)
├── configs/                               JSON configuration files
├── documentation/                         Technical writeups & API reference
└── references/                            Source material
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

# Future Work

- [X] Document all work thoroughly with informational and conversational overviews of all modules
- [X] Prioritize API calls over CLI calls for code library interfacing
- [ ] Update surfboard and fin shapes against reference images and STLs for maximum realism
- [X] Investigate PicoGK geometry generation methodologies: is there a better way than packed spheres?
- [ ] Parallel develop surface-mesh based geometry generation through Python
- [X] Consolidate the references and documentation folders
- [ ] Extend WaterSim to 3D with wave-maker boundaries for breaking wave modeling

## Technology

| Component           | Language              | Purpose                         |
| ------------------- | --------------------- | ------------------------------- |
| Surfboard Geometry  | C# / PicoGK           | Voxel-based 3D model generation |
| Physics Simulations | Python / NumPy        | Wave theory & hydrodynamics     |
| Water Simulation    | Python / NumPy        | SPH sloshing tank simulation    |
| 3D Visualization    | JavaScript / Three.js | Interactive web viewer          |
| Animations          | Python / Manim        | Physics education videos        |

## Documentation

Detailed documentation is in the `documentation/` directory:

- **[API Reference](documentation/api-reference.md)** — Quick reference for all importable classes and functions
- **[SurfPhysics Overview](documentation/surf-physics-overview.md)** — Wave theory, buoyancy, and force balance physics
- **[WaterSim Overview](documentation/water-sim-overview.md)** — SPH simulation theory and implementation
- **[SurfAnimations Overview](documentation/surf-animations-overview.md)** — Manim scene architecture and components
- **[SurfboardGeometry Overview](documentation/surfboard-geometry-overview.md)** — Voxel-based geometry generation
- **[Geometry Implementation](documentation/surfboard-geometry-implementation.md)** — Mathematical details of outline, rocker, and cross-section

---

## References

- [PicoGK by LEAP 71](https://github.com/leap71/PicoGK) -- Computational geometry kernel
- Paine, Frank. *The Science of Surfing*. Cambridge University Press.
- Dean & Dalrymple. *Water Wave Mechanics for Engineers and Scientists*.
- Finney, Ben & Houston, James. *Surfing: A History of the Ancient Hawaiian Sport*.

---

Sean Bowman | [GitHub Pages Portfolio](https://YOUR_USERNAME.github.io)
