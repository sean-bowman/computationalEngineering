# Surfboard Geometry -- Implementation & Theory Reference

How external references and surfboard shaping theory map to the code
in the SurfboardGeometry project.

See `references/surfboard-shaping-references.md` for full citations.

Sean Bowman [02/03/2026]

---

## 1. Outline -- Planform Shape (`Surfboard/Outline.cs`)

### What It Does
Computes the board's half-width at each position along the length (the
planform shape when viewed from above). The board is symmetric about
its centerline.

### Mathematical Approach
Uses a **two-piece power curve** split at the wide point:

**Nose side** (t = 0 to widePointT):
```
halfWidth = maxHalfWidth * pow(t / widePointT, noseExponent)
```

**Tail side** (t = widePointT to 1):
```
halfWidth = tailTipHalf + (maxHalfWidth - tailTipHalf) * pow(1 - localT, tailExponent)
```

The exponents are **calibrated from measured widths** using logarithmic
inversion. Given a target width at a known station (e.g., nose width at
12" from the tip):
```
exponent = log(targetRatio) / log(positionRatio)
```

This produces smooth, convex curves with no inflection points -- matching
the rail lines of real surfboards.

### Why Power Curves
Power curves (`x^n`) are used instead of polynomials or splines because:
- They guarantee a **convex outline** (no unwanted wiggles or dips)
- A single exponent controls the entire curve character
- They can be calibrated from just two control points (wide point + one measured width)

Traditional shaping uses templates and battens to create similar smooth
curves. The power curve approach is the mathematical equivalent of a
flexible batten constrained at key stations.

### References
- Greenlight Surf Supply -- Outline Design Guide (wide point placement,
  nose/tail width effects on performance)
- surfhydrodynamics.com -- Outline hydrodynamics (Lindsay Lord's optimal
  width/length ratio of ~0.5)
- Natural Curves Surfboards -- Shaper's Journal (station-based measurements)

---

## 2. Rocker Profile -- Bottom Curvature (`Surfboard/RockerProfile.cs`)

### What It Does
Computes the Z-offset (vertical position) of the bottom surface at each
point along the board's length. This is the curvature seen from the side.

### Mathematical Approach
Uses a **composite power-law curve** with different exponents for nose
and tail regions, split at the flat spot (60% from nose):

**Nose region** (t < flatSpot):
```
z = noseRocker * (1 - t/flatSpot)^2.5
```

**Tail region** (t > flatSpot):
```
z = tailRocker * ((t - flatSpot)/(1 - flatSpot))^2.0
```

### Design Rationale
The **nose exponent (2.5)** creates a steep nose kick that transitions
to a flat center section. This matches real shortboard rocker profiles
where the nose has aggressive entry rocker to prevent pearling, while
the center stays flat for speed.

The **tail exponent (2.0)** creates a smoother, more gradual curve.
Tail rocker is typically less aggressive than nose rocker to maintain
drive through turns.

The **flat spot at 60%** from the nose places the lowest point of the
board slightly behind center. This is where the board contacts the water
at rest and initiates planing. The position aligns with the surfer's
rear foot placement.

### Planing Hull Theory
The rocker profile is fundamentally a planing hull design. Bob Simmons
(1949) was the first to apply naval planing hull research to surfboards,
incorporating exaggerated nose lift and a flat-to-slightly-concave
midsection. The key insight: entry rocker creates lift at lower speeds,
but the board planes off the flatter middle section at surfing speed.

From the Greenlight Rocker Design Guide: 'Rocker-induced lift is created
mostly in the entry rocker section -- the first 20-30% of the board
length -- where the bottom is dramatically inclined relative to the
water surface.'

### References
- Greenlight Surf Supply -- Rocker & Foil Design Guide
- surfhydrodynamics.com -- Rocker shape influence on lift and drag
- Surf Simply -- Bob Simmons planing hull theory (MIT naval research)
- Lavery et al. (2018) -- CFD comparison of surfboards with modified rocker

---

## 3. Cross-Section (`Surfboard/CrossSection.cs`)

### What It Does
Defines the board's cross-sectional profile at each station: deck height
(top surface), bottom height (bottom surface), thickness distribution,
deck crown, and bottom concave.

### Mathematical Approach

#### Thickness Distribution
The board thickness varies along its length following a **cosine
envelope** with asymmetric nose/tail falloff:

```
Thickest point: 40% from nose
Nose thickness: 35% of max at tip
Tail thickness: 45% of max at tip
Interpolation: CosineInterpolate(a, b, t) = a + (b-a) * (1-cos(t*pi))/2
```

The 40% thickest-point position places maximum volume slightly ahead of
center, under the surfer's chest when paddling. This provides optimal
paddle power while keeping the performance section (center-to-tail)
thinner for responsiveness.

#### Deck Crown (Cosine Dome)
```
z_deck = halfThickness + crown * cos(lateralFraction * pi/2)
```

The crown follows a **cosine bell curve** along the board length, peaking
at 45% from the nose with falloff exponent 1.5:
```
crownFalloff = DeckCrown * cos(normalizedDistance * pi/2)^1.5
```

This adds volume at the center stringer while keeping rails thin --
exactly how real shapers foil their boards.

#### Bottom Concave (Cosine Valley)
```
z_bottom = -(halfThickness + concave * cos(lateralFraction * pi/2)^2)
```

Concave is active between 25% and 75% of the board length, with cosine
fade-in/fade-out at the boundaries. This creates a water channel along
the centerline in the planing section.

From the Greenlight Bottom Contour Guide: 'Concaves remove material from
the bottom but leave the rail line untouched. Flattening the bottom
rocker reduces curvature, creating flatter planing areas that help the
board plane earlier at lower speeds.'

#### Cosine Interpolation
All smooth transitions use the cosine interpolation function:
```
blend = (1 - cos(t * pi)) / 2
result = a + (b - a) * blend
```

This provides C1-continuous (tangent-continuous) transitions without
the oscillation problems of polynomial interpolation. The cosine easing
function is a standard technique in computer graphics and animation.

### References
- Greenlight Surf Supply -- Bottom Contour Design Guide
- Matveev (2024) -- Windsurfing board bottom modifications (concave effects)
- Simmons planing hull theory (concave channeling water flow)

---

## 4. Fin System (`Surfboard/FinSystem.cs`)

### What It Does
Generates surfboard fin geometry using spatial painting. Supports four
configurations: thruster (3 fins), twin (2 fins), quad (4 fins), and
single (1 fin).

### Fin Foil Shape -- NACA-Inspired Thickness Distribution

The fin cross-section uses a simplified NACA-inspired foil:
```
foilThickness = maxThick * sin(cT * pi)^0.6
```
where `cT` is the normalized chord position (0 = leading edge, 1 = trailing edge).

**Comparison to classical NACA 4-digit:**
The standard NACA symmetric airfoil thickness distribution is:
```
y = (t/0.2) * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)
```

Both produce similar profiles: thickest at approximately 30% chord,
tapering to near-zero at leading and trailing edges. The sin^0.6
formulation is computationally simpler while achieving the same essential
shape characteristics. The NACA equations were themselves derived
empirically (from successful 1930s airfoils), so both approaches are
empirical approximations to real-world foil shapes.

Key differences:
- NACA has a sharper leading edge (sqrt(x) term) and a slightly different
  thickness peak position
- The sin^0.6 version is symmetric about 50% chord, while NACA peaks at
  30% chord
- For surfboard fins at typical Reynolds numbers, the performance
  difference is negligible

### Fin Outline Geometry

Each fin's outline is defined by independent leading and trailing edge
curves that converge at a swept-back tip:

```
Leading edge:  x = baseX + (tipX - baseX) * hT^1.8
Trailing edge: x = baseX + chord + (tipX - (baseX+chord)) * hT^1.8
```

The exponent 1.8 produces a concave leading edge sweep (when viewed from
the water side), matching real fin geometry where the leading edge sweeps
back increasingly toward the tip.

### Cant Angle Mechanics

Side fins are tilted outward by a cant angle using 3D rotation:
```
finZ = baseZ - height * cos(cantAngle)
finYOffset = height * sin(cantAngle)
```

Higher cant angles (7 degrees for twin) project the fin's drive force
outward, generating speed through turns. Lower cant (3 degrees for quad
rear fins) provides straighter tracking.

### Fin Dimensions

All fin dimensions scale linearly with board length:
```
scaleFactor = boardLength / 1828mm  (normalized to 6'0" shortboard)
```

Dimensions are calibrated from industry fin specifications (FCS II and
Futures Fins). See `references/fin-dimensions-reference.md` for the
complete dimension tables and rationale.

### References
- NACA airfoil series (Wikipedia, Stanford, PDAS, Aerospaceweb)
- Gudimetla et al. (2019) -- CFD of 3-fin setup with NACA parametrization
- FCS Fin Data -- Industry fin dimension specifications
- Futures Fins -- F6 Legacy Quad, Performance Single specifications
- Nature Scientific Reports (2025) -- Fin pressure measurements

---

## 5. Spatial Painting (`Surfboard/SurfboardBody.cs`)

### What It Does
Generates the complete surfboard solid body by sweeping cross-sections
along the board length and filling the volume with overlapping spheres.

### Algorithm

```
For each longitudinal station (t = 0 to 1):
  1. Query Outline for half-width at this station
  2. Query RockerProfile for Z-offset (bottom curvature)
  3. For each lateral position across the width:
     a. Query CrossSection for deck height and bottom height
     b. Paint spheres from bottom to deck (vertical fill)
     c. Mirror to negative Y side (symmetry)
  4. Convert sphere lattice to voxel field
  5. Boolean-add fin voxels to body voxels
  6. Convert voxels to triangle mesh and export STL
```

### PicoGK Voxel Approach

The geometry is built using PicoGK's **spatial painting** technique:
overlapping spheres are placed into a `Lattice` object, which is then
converted to a `Voxels` field (3D grid of on/off cells). The voxel
representation has key advantages:

- **Robustness:** Boolean operations (union, subtract, intersect) are
  trivial and never fail, unlike B-rep CAD operations
- **Organic shapes:** Complex curvature is handled naturally by sphere
  overlap density
- **3D printing compatibility:** Voxel fields map directly to additive
  manufacturing processes

**Sphere sizing:** Radius is set to ~2.5x the voxel size for proper
overlap. Adjacent spheres merge when converted to voxels, creating a
continuous solid surface.

**Resolution control:** The voxel size parameter (0.25mm to 2.0mm)
controls the trade-off between surface quality and computation time.

### References
- PicoGK.org -- Coding for Engineers tutorial series (spatial painting)
- Lissner & Kayser (2023) -- Fundamentals of Computational Engineering
- LEAP 71 ShapeKernel -- Computational geometry patterns

---

## 6. Physical Constants (`Utils/Constants.cs`)

### Seawater Properties
- **Density:** 1025 kg/m^3 (standard seawater at 20C)
- **Kinematic viscosity:** 1.05e-6 m^2/s (for Reynolds number calculations)
- **Gravity:** 9.81 m/s^2

These values follow standard oceanographic references and are used for
buoyancy estimation and planned hydrodynamic force calculations.

### Material Properties
- **PU foam density:** 35 kg/m^3 (polyurethane blank, traditional construction)
- **EPS foam density:** 20 kg/m^3 (expanded polystyrene, epoxy construction)
- **Fiberglass density:** 1800 kg/m^3 (resin-saturated glass cloth)
- **Shell thickness:** 1.5mm (typical 4+4oz bottom / 4oz deck glass schedule)

### References
- Standard seawater properties from oceanographic references
- Surfboard material densities from industry construction guides
