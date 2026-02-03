# Computational Engineering Portfolio

**Generation and analysis of physical models through computer algorithms.**

This repository demonstrates computational engineering -- the process of using software to define, generate, and analyze physical objects and their behavior. Rather than drawing shapes manually in CAD software, geometry is defined programmatically through mathematical functions and parameters, enabling precise control, reproducibility, and design exploration.

---

## Surfboard Geometry (PicoGK / C#)

The primary project generates parametric surfboard geometry using [PicoGK](https://github.com/leap71/PicoGK), LEAP 71's open-source computational geometry kernel. A surfboard's shape is defined by three mathematical profiles that are combined during voxel-based spatial painting:

```
  TOP VIEW (Outline)          SIDE VIEW (Rocker)         CROSS-SECTION
  ─────────────────           ──────────────────         ─────────────

       Nose                    Nose Rocker                  Deck Crown
      ╱    ╲                  ╱                           ╱‾‾‾‾‾‾‾‾‾‾╲
    ╱        ╲              ╱                            │     Volume   │
   │  Wide    │           ╱                               ╲___________╱
   │  Point   │         ╱_______Flat_______╲               Concave
    ╲        ╱                              ╲
      ╲____╱                            Tail Rocker
       Tail
```

### How It Works

1. **Outline** defines the board's planform shape (half-width at each position along the length)
2. **Rocker Profile** defines the bottom curvature from nose to tail
3. **Cross-Section** defines the deck crown, bottom concave, and thickness at each station

These three profiles are combined using **spatial painting** -- placing thousands of overlapping spheres along the board surface to build up the shape as a voxel field. The voxels are then converted to a triangle mesh and exported as an STL file.

### Board Presets

| Parameter | Shortboard | Longboard | Fish |
|-----------|-----------|-----------|------|
| Length | 6'0" (1828mm) | 9'0" (2743mm) | 5'6" (1676mm) |
| Width | 19.5" (495mm) | 22.5" (570mm) | 21" (533mm) |
| Thickness | 2.44" (62mm) | 3.0" (75mm) | 2.56" (65mm) |
| Nose Rocker | 4.7" (120mm) | 7.0" (180mm) | 3.1" (80mm) |
| Tail Rocker | 1.6" (40mm) | 1.0" (25mm) | 1.2" (30mm) |

All dimensions are configurable through `SurfboardParameters`.

### Fin Configurations

| Configuration | Fins | Default For | Characteristics |
| ------------- | ---- | ----------- | --------------- |
| Thruster | 3 (center + 2 sides) | Shortboard | Versatile, good control |
| Twin | 2 (sides only) | Fish | Fast, loose, skatey feel |
| Quad | 4 (2 front + 2 rear) | -- | Speed and hold in large surf |
| Single | 1 (center only) | Longboard | Smooth turns, stability |

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

# Show help
dotnet run --project SurfboardGeometry/ -- --help
```

Output STL files are saved to `SurfboardGeometry/Output/`.

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
├── PicoGK/                                LEAP 71 geometry kernel (submodule)
├── SurfPhysics/                           Python physics simulations (planned)
├── documentation/                         Technical writeups
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

## Roadmap

- [x] Parametric surfboard body generation (PicoGK)
- [x] Fin system geometry (thruster, twin, quad, single configurations)
- [ ] Wave physics modeling (linear wave theory, breaking criteria)
- [ ] Hydrodynamic force analysis (buoyancy, drag, planing lift)
- [ ] Interactive 3D visualization (Three.js)
- [ ] Physics animations (Manim)

## Technology

| Component | Language | Purpose |
|-----------|----------|---------|
| Surfboard Geometry | C# / PicoGK | Voxel-based 3D model generation |
| Physics Simulations | Python | Wave theory & hydrodynamics (planned) |
| 3D Visualization | JavaScript / Three.js | Interactive web viewer (planned) |
| Animations | Python / Manim | Physics education videos (planned) |

## References

- [PicoGK by LEAP 71](https://github.com/leap71/PicoGK) -- Computational geometry kernel
- Paine, Frank. *The Science of Surfing*. Cambridge University Press.
- Dean & Dalrymple. *Water Wave Mechanics for Engineers and Scientists*.
- Finney, Ben & Houston, James. *Surfing: A History of the Ancient Hawaiian Sport*.

---

Sean Bowman | [GitHub Pages Portfolio](https://YOUR_USERNAME.github.io)
