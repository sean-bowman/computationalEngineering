# Surfboard Geometry -- Module Overview & Python API

A high-level guide to the SurfboardGeometry C# module: what it does, why
it works the way it does, how to use it from Python, and how it connects
to the SurfPhysics simulation toolkit.

For the mathematical details behind each profile class, see
`documentation/surfboard-geometry-implementation.md`.

Sean Bowman [02/05/2026]

---

## 1. What It Does

The SurfboardGeometry module generates parametric surfboard shapes as STL
mesh files. You give it a board type (shortboard, longboard, or fish) and
a voxel resolution, and it produces a watertight triangle mesh ready for
3D printing, CFD simulation, or interactive viewing.

Under the hood, the module uses LEAP 71's PicoGK geometry kernel to build
the board via **spatial painting** -- placing thousands of overlapping
spheres into a 3D lattice, converting that lattice to a voxel field, and
then extracting a triangle mesh surface. The result is a solid, organic
shape with smooth surfaces and no geometric artifacts.

The entry point is `SurfboardGeometry/Program.cs`, which parses CLI
arguments (board type, fin configuration, voxel size), selects a parameter
preset, initializes PicoGK in headless mode, and kicks off geometry
generation in `SurfboardBody.Task()`. From Python, the
`GeometryGenerator` class in `SurfboardGeometry/geometryGenerator.py`
wraps this entire workflow as a subprocess call, so you never need to
touch the C# code directly.

---

## 2. Why Voxels?

Traditional CAD systems use **B-rep** (boundary representation) -- surfaces
defined by NURBS patches or analytic equations. B-rep is great for
machined parts with exact dimensions, but it has a few pain points for
organic shapes like surfboards:

- **Boolean operations are fragile.** Unioning a fin to a board body, or
  subtracting a swallow tail notch, requires the CAD kernel to intersect
  surfaces analytically. If edges are nearly tangent or surfaces are
  coincident, the operation can fail silently or produce corrupt geometry.

- **Organic curvature is hard to control.** Surfboards have continuously
  varying cross-sections, rocker curves, and outline shapes. Defining all
  of this with NURBS surfaces means managing a large number of control
  points and surface patches.

- **Mesh export can produce artifacts.** Converting B-rep surfaces to STL
  requires tessellation, which can introduce gaps, overlapping triangles,
  or non-manifold edges.

PicoGK's voxel approach sidesteps all of these issues:

- **Boolean operations are trivial.** Union is a per-voxel OR. Subtraction
  is a per-voxel AND-NOT. They never fail, regardless of geometry
  complexity. This is how `SurfboardBody.GenerateAndExport()` attaches
  fins (boolean add) and carves the swallow tail notch (boolean subtract)
  without any special-case handling.

- **Organic shapes are natural.** Spatial painting with overlapping spheres
  automatically produces smooth, organic surfaces. The sphere overlap
  density determines the local surface curvature.

- **Output is always watertight.** Marching cubes mesh extraction from
  voxels produces a guaranteed manifold, closed mesh -- exactly what 3D
  printers and simulation tools need.

The trade-off is resolution: voxels discretize space into a grid, so very
fine features (like sharp trailing edges on fins) require a small voxel
size, which increases memory and computation time. The `--voxel` flag lets
you tune this: 2.0 mm for fast previews, 0.25 mm for production quality.

---

## 3. The Three Profiles

A surfboard's shape is fully defined by three profile functions that are
composed during geometry generation. Each one maps a normalized position
along the board (0 = nose tip, 1 = tail tip) to a geometric property.

### Outline -- Planform Shape (`Surfboard/Outline.cs`)

The outline defines the board's half-width at each station along the
length -- the shape you see when looking down at the board from above.
It controls wave-catching ability (wider = more planing area), turning
radius (narrower = tighter arcs), and overall aesthetic.

The outline uses a two-piece power curve split at the wide point, with
exponents calibrated from the nose width and tail width parameters. The
result is a smooth, convex curve with no inflection points -- like a
flexible drafting batten pinned at the key measurement stations.

### Rocker -- Side View Curvature (`Surfboard/RockerProfile.cs`)

Rocker is the curve of the board's bottom when viewed from the side. The
nose kicks up to prevent nosedives on steep takeoffs, while the tail
curves up to enable tight turns. Between them, a flat section (the
planing area) provides speed.

The flat spot sits at 60% from the nose, roughly where your rear foot
goes. The nose rises with a power-law exponent of 2.5 (steep kick), and
the tail rises with an exponent of 2.0 (smoother curve).

### Cross-Section -- Deck, Bottom, and Rails (`Surfboard/CrossSection.cs`)

The cross-section defines the board's thickness at each point: the dome
of the deck (cosine crown), the concave channel on the bottom, and how
they thin toward the nose and tail. The thickness envelope peaks at 40%
from the nose -- right under the surfer's chest for paddle power.

Deck crown adds volume at the center stringer while keeping the rails
thin. Bottom concave channels water flow for speed. Both features fade
to zero at the nose and tail using cosine interpolation for smooth
transitions.

For the full math behind each of these, including the calibration
equations and design rationale, see the
[implementation doc](surfboard-geometry-implementation.md).

---

## 4. Fin System (`Surfboard/FinSystem.cs`)

Fins are generated separately from the board body, then attached via
PicoGK boolean union. The `FinSystem` class supports four configurations:

| Configuration | Fins | Description |
|---------------|------|-------------|
| **Thruster** | 3 | Center + 2 sides. Most versatile setup. |
| **Twin** | 2 | Two large side fins. Fast, loose, classic fish feel. |
| **Quad** | 4 | 2 front + 2 rear trailers. Speed and hold in bigger surf. |
| **Single** | 1 | One large center fin. Traditional longboard stability. |

Each fin is built using the same spatial painting technique as the board
body. The foil cross-section uses a NACA-inspired thickness distribution
(`sin(t * pi)^0.6`) that produces maximum thickness at roughly 30% chord
and tapers to near-zero at the leading and trailing edges. The fin outline
is defined by independent leading and trailing edge curves (power law with
exponent 1.8) that converge at a swept-back tip.

Fin dimensions scale linearly with board length, normalized to a 6'0"
shortboard baseline. The scaling factor is simply `boardLength / 1828mm`.
Side fins are canted outward (5 degrees for thruster, 7 degrees for twin)
to project drive force through turns. Dimensions are calibrated from
real-world FCS II and Futures Fins specifications.

The fish preset also features a **swallow tail** -- a U-shaped notch
carved from the tail via boolean subtraction (`SwallowTailNotch.cs`).
The notch is defined by a signed distance function (PicoGK's
`IBoundedImplicit` interface) with a square-root profile so the channel
opens up quickly from a rounded apex, matching how real fish tails look.

---

## 5. Board Presets

Three presets are built into `SurfboardParameters.cs` as static factory
properties. Each one captures the essential character of that board type:

| Parameter | Shortboard | Longboard | Fish (RNF) |
|-----------|-----------|-----------|------------|
| **Length** | 6'0" (1828 mm) | 9'0" (2743 mm) | 5'6" (1676 mm) |
| **Max Width** | 19.5" (495 mm) | 22.5" (570 mm) | 21.0" (533 mm) |
| **Max Thickness** | 2.44" (62 mm) | 3.0" (75 mm) | 2.56" (65 mm) |
| **Nose Width** | 11.8" (300 mm) | 17.0" (430 mm) | 15.75" (400 mm) |
| **Tail Width** | 15.0" (380 mm) | 14.5" (370 mm) | 16.9" (430 mm) |
| **Nose Rocker** | 4.7" (120 mm) | 7.0" (180 mm) | 3.1" (80 mm) |
| **Tail Rocker** | 1.6" (40 mm) | 1.0" (25 mm) | 1.2" (30 mm) |
| **Bottom Concave** | 2 mm | 1 mm | 3 mm |
| **Tail Shape** | Standard | Standard | Swallow |
| **Default Fins** | Thruster | Single | Twin |
| **Character** | Responsive, steep waves | Glide, noseriding | Speed, small waves |

The shortboard is the default -- a 6'0" performance board with moderate
rocker and a thruster setup. The longboard has a wide nose for
noseriding, minimal tail rocker for speed, and soft round rails. The
fish packs extra volume into a short frame with a flat rocker, deep
concave, and the signature swallow tail for clean water release.

---

## 6. Python API

### Setup

The `SurfboardGeometry` directory is both a C# project and a Python
package. The `__init__.py` at the package root exports the
`GeometryGenerator` class:

```python
from SurfboardGeometry import GeometryGenerator
```

### Generating from a Config File

The simplest way to generate a board is from one of the JSON config files
in `configs/`. These configs contain board type, fin configuration, voxel
size, and output directory settings:

```python
from SurfboardGeometry import GeometryGenerator

generator = GeometryGenerator()
stlPath = generator.generateFromConfig('configs/shortboard_default.json')

print(f'Board type: {generator.boardType}')
print(f'STL file:   {stlPath}')
```

The `generateFromConfig()` method reads the `board` and `geometry`
sections from the JSON, builds a `dotnet run` command, and invokes the
C# CLI as a subprocess. When it finishes, you get back the path to the
generated STL file.

### Config File Structure

The config files contain sections for the full pipeline (geometry, waves,
physics, viewer), but `GeometryGenerator` only reads these keys:

```json
{
    "board": {
        "type": "shortboard",
        "finConfiguration": "thruster"
    },
    "geometry": {
        "voxelSize": 0.5,
        "outputDir": "SurfboardGeometry/Output"
    }
}
```

Valid board types: `'shortboard'`, `'longboard'`, `'fish'`.
Valid fin configurations: `'thruster'`, `'twin'`, `'quad'`, `'single'`,
or `'default'` (uses the board type's default).

### Error Handling

If the C# process fails (missing PicoGK runtime, bad arguments, etc.),
`generateFromConfig()` raises a `subprocess.CalledProcessError` with the
exit code and stderr. If the process succeeds but the STL file is not
found at the expected path, it raises `FileNotFoundError`.

```python
try:
    stlPath = generator.generateFromConfig('configs/fish_overhead.json')
except subprocess.CalledProcessError as e:
    print(f'Geometry generation failed: exit code {e.returncode}')
except FileNotFoundError as e:
    print(f'STL not found: {e}')
```

---

## 7. CLI Usage

You can also invoke the C# generator directly without Python. The project
root is `SurfboardGeometry/`:

```bash
# Default shortboard with thruster fins
dotnet run --project SurfboardGeometry

# 9'0" longboard with single fin
dotnet run --project SurfboardGeometry -- --type longboard

# Fish board at high resolution
dotnet run --project SurfboardGeometry -- --type fish --voxel 0.25

# Shortboard with quad fins, fast preview
dotnet run --project SurfboardGeometry -- --fins quad --voxel 2

# Generate all four fin configurations as separate STL files
dotnet run --project SurfboardGeometry -- --all-fins --voxel 2

# Show help
dotnet run --project SurfboardGeometry -- --help
```

### CLI Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--type <type>` | `shortboard` | Board preset: `shortboard`, `longboard`, or `fish` |
| `--fins <config>` | Per board type | Fin setup: `thruster`, `twin`, `quad`, or `single` |
| `--voxel <size>` | `0.5` | Voxel resolution in mm (smaller = finer, slower) |
| `--all-fins` | off | Export all four fin configs as separate STL files |
| `--help` | -- | Print usage information |

The `--` separator before flags is needed when running through
`dotnet run` so that the arguments are passed to the program rather than
consumed by the dotnet CLI itself.

---

## 8. Relationship to SurfPhysics

The project has two parallel geometry representations that define the same
surfboard shapes:

- **C# SurfboardGeometry** -- Generates voxel-based 3D geometry and
  exports STL mesh files. This is the geometry you can 3D print, render,
  or import into external tools.

- **Python SurfPhysics.geometry** -- A pure-Python re-implementation of
  the same parametric profiles (outline, rocker, cross-section) for
  physics calculations. This is in `SurfPhysics/geometry/` with files
  like `outline.py`, `rocker.py`, `crossSection.py`, and `board.py`.

Both sides share the same parameter definitions. The Python
`SurfboardParameters` dataclass (`SurfPhysics/geometry/parameters.py`)
mirrors the C# `SurfboardParameters` class field-for-field, including the
same factory presets (shortboard, longboard, fish) with identical
dimension values. The Python version adds JSON serialization for
cross-language interop.

**Why two implementations?** The C# side needs PicoGK for spatial
painting and voxel operations -- that is a compiled library with specific
runtime requirements. The Python side needs to evaluate the same profile
curves at arbitrary points for numerical integration (volume, planform
area, wetted surface area, submerged volume) and for the hydrodynamics
force calculations. Keeping the physics in Python means it can use NumPy,
SciPy, and Plotly without any .NET dependency.

The `BoardGeometry` facade (`SurfPhysics/geometry/board.py`) composes the
Python Outline, RockerProfile, and CrossSection classes into a single
interface that the hydrodynamics modules call. It handles the mm-to-meters
conversion at the boundary, so all physics code works in SI units.

In practice, you typically run both together through a config file:

1. `GeometryGenerator` calls `dotnet run` to produce the STL
2. `SurfPhysics` builds a `BoardGeometry` from the same parameter preset
3. The physics modules compute forces, drag, buoyancy, etc.
4. Results are visualized with Plotly dashboards

The config files in `configs/` drive this entire pipeline. A single JSON
file specifies the board type for geometry generation, the wave conditions
for physics, and the analysis/viewer settings.

---

## 9. Output

### STL Files

Generated STL files go to `SurfboardGeometry/Output/` by default (or
wherever the config's `geometry.outputDir` points). File names depend on
the generation mode:

- **Single fin config:** `surfboard.stl`
- **All fin configs (`--all-fins`):** `surfboard_thruster.stl`,
  `surfboard_twin.stl`, `surfboard_quad.stl`, `surfboard_single.stl`

The STL files are binary format, watertight manifold meshes. Triangle
count depends on voxel resolution -- a 0.5 mm voxel board typically
produces a mesh in the hundreds of thousands of triangles.

### Coordinate System

The STL geometry uses the following convention:

- **X** -- Along the board length (0 = nose tip, increasing toward tail)
- **Y** -- Across the board width (0 = centerline, positive = right)
- **Z** -- Vertical (0 = rocker datum at the flat spot, positive = up)

All dimensions in the STL are in **millimeters**.

### Using the Output

The generated STL files can be:

- **3D printed** directly (the mesh is watertight and manifold)
- **Imported into CAD** tools (Fusion 360, SolidWorks) for further
  modification
- **Loaded into CFD** software for hydrodynamic simulation
- **Viewed** in the SurfViewer web component (Three.js-based viewer)
- **Sampled** by `SurfPhysics/export/meshSampler.py` for physics
  validation against the parametric model

---

## Key Source Files

| File | Purpose |
|------|---------|
| `SurfboardGeometry/Program.cs` | CLI entry point, argument parsing |
| `SurfboardGeometry/Surfboard/SurfboardParameters.cs` | Parametric dimensions and presets |
| `SurfboardGeometry/Surfboard/SurfboardBody.cs` | Spatial painting geometry generation |
| `SurfboardGeometry/Surfboard/Outline.cs` | Planform shape (top view half-width) |
| `SurfboardGeometry/Surfboard/RockerProfile.cs` | Bottom curvature (side view) |
| `SurfboardGeometry/Surfboard/CrossSection.cs` | Deck/bottom profile and thickness |
| `SurfboardGeometry/Surfboard/FinSystem.cs` | Fin geometry (4 configurations) |
| `SurfboardGeometry/Surfboard/FinConfiguration.cs` | Fin config enum |
| `SurfboardGeometry/Surfboard/TailShape.cs` | Tail shape enum (Standard, Swallow) |
| `SurfboardGeometry/Surfboard/SwallowTailNotch.cs` | Swallow tail boolean subtraction |
| `SurfboardGeometry/geometryGenerator.py` | Python API wrapper |
| `SurfboardGeometry/__init__.py` | Package export |
| `SurfPhysics/geometry/parameters.py` | Python mirror of SurfboardParameters |
| `SurfPhysics/geometry/board.py` | Python geometry facade for physics |
