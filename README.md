# Computational Engineering Portfolio

**Generation and analysis of physical models through computer algorithms.**

This repository demonstrates computational engineering — using software to define, generate, and analyze physical objects and their behavior. Rather than drawing shapes manually in CAD software, geometry is defined programmatically through mathematical functions and parameters, enabling precise control, reproducibility, and design exploration.

Domains live under `computationalEngineering/`. Each domain is a self-contained Python package with its own configs, documentation, and code. Generic tools like `PicoGK` and `WaterSim` are organized for cross-domain reuse.

---

## Domains

### [Surfboard →](computationalEngineering/Surfboard/README.md)

Parametric surfboard design, hydrodynamic physics, voxel geometry generation, reference validation, and animated visualizations. Includes two geometry generation pipelines (Python parametric surface mesh and C# voxel-based via PicoGK) that are benchmarked against a real reference board.

```bash
python codeInterface.py   # Reverse-engineer reference STL and launch comparison viewer
```

### [WaterSim →](computationalEngineering/WaterSim/README.md)

Weakly Compressible SPH (WCSPH) solver for free-surface water flows. Supports 2D/3D sloshing tanks and numerical wave tanks with breaking wave detection. Used for physics simulations and Manim animations.

```bash
python -m computationalEngineering.WaterSim   # Run default sloshing tank
```

---

## Project Structure

```txt
repo_root/
├── codeInterface.py                 Python entry point (reverse-engineering demo)
├── README.md
│
└── computationalEngineering/
    ├── Surfboard/                   Surfboard design domain → see Surfboard/README.md
    │   ├── SurfPhysics/             Python physics, geometry, validation
    │   ├── SurfboardGeometry/       C# / PicoGK voxel geometry
    │   ├── SurfAnimations/          Manim physics animations
    │   ├── SurfViewer/              Three.js interactive 3D viewer
    │   ├── configs/                 JSON configuration files
    │   └── documentation/           Technical writeups
    │
    ├── WaterSim/                    SPH water simulation → see WaterSim/README.md
    │   ├── sph/                     WCSPH solver engine
    │   ├── scenarios/               Pre-configured simulations
    │   ├── export/                  Frame data export
    │   ├── configs/                 Water sim JSON configs
    │   └── documentation/           Water sim technical writeups
    │
    └── PicoGK/                      LEAP 71 geometry kernel (git submodule)
```

---

## Prerequisites & Setup

```bash
# Clone with submodules
git clone --recurse-submodules https://github.com/YOUR_USERNAME/computationalEngineering.git
cd computationalEngineering

# Initialize submodules if already cloned
git submodule update --init --recursive

# Install Python dependencies
pip install numpy plotly trimesh scipy dash

# Build the C# surfboard geometry project
dotnet build computationalEngineering/Surfboard/SurfboardGeometry/SurfboardGeometry.csproj
```

**Requirements:**

- Python 3.x
- .NET 10.0 SDK or later (for C# voxel geometry)
- Windows x64 (for PicoGK native libraries)

---

## Technology

| Component | Language | Purpose |
| --- | --- | --- |
| Surfboard Geometry | C# / PicoGK | Voxel-based 3D model generation |
| Physics Simulations | Python / NumPy | Wave theory and hydrodynamics |
| Water Simulation | Python / NumPy | WCSPH sloshing and wave tanks |
| 3D Visualization | JavaScript / Three.js | Interactive web viewer |
| Physics Animations | Python / Manim | Education videos |
| Parameter Optimization | Python / SciPy | Reference-board reverse engineering |

---

## Future Work

- [X] Document all work thoroughly with informational overviews of all modules
- [X] Prioritize API calls over CLI calls for code library interfacing
- [X] Update surfboard and fin shapes against reference STLs for maximum realism
- [X] Investigate PicoGK geometry generation methodologies
- [X] Parallel-develop surface-mesh geometry generation through Python
- [X] Consolidate references and documentation folders
- [ ] Extend WaterSim to 3D with wave-maker boundaries for breaking wave modeling
- [ ] Add new geometry domains (boat hull, etc.) as sibling packages

---

## References

- [PicoGK by LEAP 71](https://github.com/leap71/PicoGK) — Computational geometry kernel
- Paine, Frank. *The Science of Surfing*. Cambridge University Press.
- Dean & Dalrymple. *Water Wave Mechanics for Engineers and Scientists*.
- Finney, Ben & Houston, James. *Surfing: A History of the Ancient Hawaiian Sport*.

---

Sean Bowman
