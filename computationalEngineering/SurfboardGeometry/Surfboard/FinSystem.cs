// =============================================================================
// FIN SYSTEM - SURFBOARD FIN GEOMETRY GENERATION
// =============================================================================
//
// SURFBOARD CONCEPT: Fin Configurations
// =======================================
//
// Fins provide directional stability and control. Without fins, a surfboard
// would slide sideways on the wave face (like a bodyboard). Different
// configurations affect performance:
//
//   THRUSTER (3 fins):          TWIN (2 fins):
//   ─────────────────           ──────────────
//        │                           │
//        │                           │
//        │                           │
//       ╱│╲  Center fin             │
//      ╱ │ ╲                       ╱ ╲
//     ╱  │  ╲  Side fins          ╱   ╲  Side fins
//    ╱   │   ╲                   ╱     ╲
//
// THRUSTER: Most common setup. Center fin provides hold and control,
//   side fins provide drive and projection through turns.
//   Good all-around performance.
//
// TWIN: Two side fins only. Fast and loose -- the classic fish setup.
//   Larger fins compensate for missing center fin. Higher cant for drive.
//
// QUAD: Four fins in two pairs. Front pair provides drive, rear
//   "trailer" fins add hold and direction. Fast in open-face waves.
//
//   QUAD (4 fins):            SINGLE (1 fin):
//   ──────────────            ────────────────
//        │                           │
//        │                           │
//       ╱ ╲  Front pair             │
//      ╱   ╲                        │
//     ╱ ╲ ╱ ╲  Rear pair           │
//    ╱   ╳   ╲                     │╲  Center fin
//                                   │ ╲
//
// SINGLE: One large center fin. Smooth, drawn-out turns. Traditional
//   longboard setup providing excellent directional stability.
//
// FIN ANATOMY:
// ────────────
//
//   Leading Edge ─╲
//                   ╲
//                    │  ← Foil (thickest ~30% from leading edge)
//                   ╱
//   Trailing Edge ╱
//
//   ↕ Height (base to tip)
//   ↔ Base (chord at bottom)
//   Rake: how far the tip sweeps back from the base
//
// =============================================================================

using System.Numerics;
using PicoGK;

namespace SurfboardGeometry.Surfboard;

/// <summary>
/// Generates surfboard fin geometry using spatial painting.
/// </summary>
/// <remarks>
/// <para>
/// Creates fins as simplified foil shapes and positions them according to
/// the selected configuration (thruster, twin, single). Each fin is built
/// by painting overlapping spheres along the foil profile at each height
/// station from base to tip.
/// </para>
/// <para>
/// PICOGK: Uses the same spatial painting technique as the surfboard body
/// and BellNozzle -- overlapping spheres create a continuous solid when
/// converted to voxels.
/// </para>
/// </remarks>
public class FinSystem
{
    private readonly SurfboardParameters _params;
    private readonly float _voxelSize;

    // Fin dimensions (all in mm)
    private readonly float _centerFinHeight;
    private readonly float _sideFinHeight;
    private readonly float _centerFinBase;    // chord length at base
    private readonly float _sideFinBase;
    private readonly float _finThickness;     // maximum foil thickness
    private readonly float _finRake;          // tip sweep-back distance
    private readonly float _sideCant;         // side fin outward angle (degrees)

    /// <summary>
    /// Create a fin system generator.
    /// </summary>
    /// <param name="parameters">Board parameters (used for positioning)</param>
    /// <param name="voxelSize">Voxel resolution in mm</param>
    public FinSystem(SurfboardParameters parameters, float voxelSize = 0.5f)
    {
        _params = parameters;
        _voxelSize = voxelSize;

        // Fin sizing scales with board length
        float scaleFactor = parameters.Length / 1828f; // normalize to 6'0" shortboard

        _centerFinHeight = 115f * scaleFactor;   // ~4.5" on a 6'0"
        _sideFinHeight = 105f * scaleFactor;      // ~4.1" on a 6'0"
        _centerFinBase = 110f * scaleFactor;      // ~4.3" chord
        _sideFinBase = 100f * scaleFactor;        // ~3.9" chord
        _finThickness = 8f * scaleFactor;         // max foil thickness
        _finRake = 15f * scaleFactor;             // how far tip extends past trailing edge base
        _sideCant = 5f;                           // degrees outward tilt
    }

    // =========================================================================
    // CONFIGURATION DISPATCH
    // =========================================================================

    /// <summary>
    /// Generate fins for the specified configuration.
    /// </summary>
    /// <param name="configuration">Fin setup to generate</param>
    /// <param name="rockerProfile">Rocker profile for positioning fins on the bottom surface</param>
    /// <returns>Voxels containing all fins for the configuration</returns>
    public Voxels Generate(FinConfiguration configuration, RockerProfile rockerProfile)
    {
        return configuration switch
        {
            FinConfiguration.Thruster => GenerateThruster(rockerProfile),
            FinConfiguration.Twin     => GenerateTwin(rockerProfile),
            FinConfiguration.Quad     => GenerateQuad(rockerProfile),
            FinConfiguration.Single   => GenerateSingle(rockerProfile),
            _ => throw new ArgumentOutOfRangeException(
                nameof(configuration),
                configuration,
                "Unknown fin configuration")
        };
    }

    // =========================================================================
    // THRUSTER CONFIGURATION
    // =========================================================================

    /// <summary>
    /// Generate thruster (3-fin) configuration as voxels.
    /// </summary>
    /// <param name="rockerProfile">Rocker profile for positioning fins on the bottom surface</param>
    /// <returns>Voxels containing all three fins</returns>
    /// <remarks>
    /// <para>
    /// THRUSTER POSITIONING:
    /// - Center fin: on the centerline, ~80mm from the tail tip
    /// - Side fins: offset ~120mm from centerline, ~280mm from tail tip
    /// - Side fins are canted outward by 5 degrees for drive
    /// </para>
    /// </remarks>
    public Voxels GenerateThruster(RockerProfile rockerProfile)
    {
        Console.WriteLine("Generating thruster fin system (3 fins)...");

        Lattice lat = new();

        // Fin positions: baseX is the LEADING EDGE of the fin base.
        // The trailing edge extends baseChord mm further toward the tail.
        // Position so trailing edge sits ~20mm inside the tail tip.
        float centerTrailingMargin = 20f;
        float centerFinX = _params.Length - _centerFinBase - centerTrailingMargin;
        float centerFinT = centerFinX / _params.Length;
        float centerFinRockerZ = rockerProfile.GetRockerHeight(centerFinT);

        // Side fins: positioned further forward, offset laterally
        float sideTrailingMargin = 180f;
        float sideFinX = _params.Length - _sideFinBase - sideTrailingMargin;
        float sideFinT = sideFinX / _params.Length;
        float sideFinRockerZ = rockerProfile.GetRockerHeight(sideFinT);
        float sideFinYOffset = 100f; // lateral offset from centerline

        // Paint center fin (vertical, on centerline)
        Console.WriteLine("  Painting center fin...");
        PaintFin(lat,
            baseX: centerFinX,
            baseY: 0f,
            baseZ: centerFinRockerZ,
            height: _centerFinHeight,
            baseChord: _centerFinBase,
            cantAngleDeg: 0f,
            isCenter: true
        );

        // Paint right side fin (canted outward)
        Console.WriteLine("  Painting right side fin...");
        PaintFin(lat,
            baseX: sideFinX,
            baseY: sideFinYOffset,
            baseZ: sideFinRockerZ,
            height: _sideFinHeight,
            baseChord: _sideFinBase,
            cantAngleDeg: _sideCant,
            isCenter: false
        );

        // Paint left side fin (canted outward, mirrored)
        Console.WriteLine("  Painting left side fin...");
        PaintFin(lat,
            baseX: sideFinX,
            baseY: -sideFinYOffset,
            baseZ: sideFinRockerZ,
            height: _sideFinHeight,
            baseChord: _sideFinBase,
            cantAngleDeg: -_sideCant,
            isCenter: false
        );

        Console.WriteLine("  Converting fin lattice to voxels...");
        Voxels finVoxels = new(lat);

        Console.WriteLine("  ✓ Thruster fin system generated");
        return finVoxels;
    }

    // =========================================================================
    // TWIN CONFIGURATION
    // =========================================================================

    /// <summary>
    /// Generate twin (2-fin) configuration as voxels.
    /// </summary>
    /// <param name="rockerProfile">Rocker profile for positioning fins on the bottom surface</param>
    /// <returns>Voxels containing both fins</returns>
    /// <remarks>
    /// <para>
    /// TWIN POSITIONING:
    /// - Two side fins only, no center fin
    /// - Larger than thruster sides to compensate for missing center
    /// - Positioned wider laterally (~120mm from centerline)
    /// - Further forward (~220mm trailing margin) for a looser feel
    /// - Higher cant angle (7 degrees) for drive through turns
    /// </para>
    /// <para>
    /// SURFBOARD: Twin fins are the traditional fish board setup.
    /// Without a center fin, the board is looser and faster but has
    /// less hold in steep, powerful waves. Dimensions referenced from
    /// FCS II Mark Richards Twin fin specs.
    /// </para>
    /// </remarks>
    public Voxels GenerateTwin(RockerProfile rockerProfile)
    {
        Console.WriteLine("Generating twin fin system (2 fins)...");

        Lattice lat = new();

        float scaleFactor = _params.Length / 1828f;

        // Keel-style twin fins: large base and height for drive and hold
        float twinFinHeight = 140f * scaleFactor;    // ~5.5" tall keel fins
        float twinFinBase   = 145f * scaleFactor;    // ~5.7" long base chord (keel shape)
        float twinCant      = 7f;                    // degrees (vs 5 degrees thruster)

        // Positioning: close to tail for classic fish keel placement
        float trailingMargin = 180f;                 // closer to tail for fish look
        float lateralOffset  = 120f;                 // wider than thruster (100mm)

        float finX    = _params.Length - twinFinBase - trailingMargin;
        float finT    = finX / _params.Length;
        float rockerZ = rockerProfile.GetRockerHeight(finT);

        // Paint right side fin
        Console.WriteLine("  Painting right twin fin...");
        PaintFin(lat,
            baseX: finX,
            baseY: lateralOffset,
            baseZ: rockerZ,
            height: twinFinHeight,
            baseChord: twinFinBase,
            cantAngleDeg: twinCant,
            isCenter: false
        );

        // Paint left side fin (mirrored)
        Console.WriteLine("  Painting left twin fin...");
        PaintFin(lat,
            baseX: finX,
            baseY: -lateralOffset,
            baseZ: rockerZ,
            height: twinFinHeight,
            baseChord: twinFinBase,
            cantAngleDeg: -twinCant,
            isCenter: false
        );

        Console.WriteLine("  Converting fin lattice to voxels...");
        Voxels finVoxels = new(lat);

        Console.WriteLine("  ✓ Twin fin system generated");
        return finVoxels;
    }

    // =========================================================================
    // QUAD CONFIGURATION
    // =========================================================================

    /// <summary>
    /// Generate quad (4-fin) configuration as voxels.
    /// </summary>
    /// <param name="rockerProfile">Rocker profile for positioning fins on the bottom surface</param>
    /// <returns>Voxels containing all four fins</returns>
    /// <remarks>
    /// <para>
    /// QUAD POSITIONING:
    /// - Front pair: wider toward rails (~115mm lateral), ~350mm from tail
    /// - Rear pair: smaller "trailer" fins closer to center (~85mm lateral), ~140mm from tail
    /// - Front fins: 5 degree cant for drive
    /// - Rear fins: 3 degree cant for tracking stability
    /// </para>
    /// <para>
    /// SURFBOARD: Quad setups channel water between the fin pairs, generating
    /// speed. Without a center fin, water exits cleanly off the tail. The
    /// front fins provide drive while rear fins add hold and direction.
    /// Dimensions referenced from FCS II Performer Quad fin specs.
    /// </para>
    /// </remarks>
    public Voxels GenerateQuad(RockerProfile rockerProfile)
    {
        Console.WriteLine("Generating quad fin system (4 fins)...");

        Lattice lat = new();

        float scaleFactor = _params.Length / 1828f;

        // Front fins: similar to thruster side fins
        float frontFinHeight = 105f * scaleFactor;   // ~4.1"
        float frontFinBase   = 100f * scaleFactor;   // ~3.9" chord
        float frontCant      = 5f;                   // degrees

        // Rear "trailer" fins: smaller for tracking and hold
        float rearFinHeight  = 75f * scaleFactor;    // ~3.0"
        float rearFinBase    = 70f * scaleFactor;    // ~2.75" chord
        float rearCant       = 3f;                   // degrees (less cant for stability)

        // Front fin positioning: further from tail, wider toward the rails
        // Real quads have the front pair out near the rail line for drive
        float frontTrailingMargin = 250f;
        float frontLateralOffset  = 115f;
        float frontFinX   = _params.Length - frontFinBase - frontTrailingMargin;
        float frontFinT   = frontFinX / _params.Length;
        float frontRockerZ = rockerProfile.GetRockerHeight(frontFinT);

        // Rear fin positioning: close to tail, closer to centerline
        // Rear "trailer" fins tuck inboard to channel water between the pairs
        float rearTrailingMargin = 70f;
        float rearLateralOffset  = 85f;
        float rearFinX   = _params.Length - rearFinBase - rearTrailingMargin;
        float rearFinT   = rearFinX / _params.Length;
        float rearRockerZ = rockerProfile.GetRockerHeight(rearFinT);

        // Paint front right fin
        Console.WriteLine("  Painting front right fin...");
        PaintFin(lat,
            baseX: frontFinX,
            baseY: frontLateralOffset,
            baseZ: frontRockerZ,
            height: frontFinHeight,
            baseChord: frontFinBase,
            cantAngleDeg: frontCant,
            isCenter: false
        );

        // Paint front left fin (mirrored)
        Console.WriteLine("  Painting front left fin...");
        PaintFin(lat,
            baseX: frontFinX,
            baseY: -frontLateralOffset,
            baseZ: frontRockerZ,
            height: frontFinHeight,
            baseChord: frontFinBase,
            cantAngleDeg: -frontCant,
            isCenter: false
        );

        // Paint rear right trailer fin
        Console.WriteLine("  Painting rear right trailer fin...");
        PaintFin(lat,
            baseX: rearFinX,
            baseY: rearLateralOffset,
            baseZ: rearRockerZ,
            height: rearFinHeight,
            baseChord: rearFinBase,
            cantAngleDeg: rearCant,
            isCenter: false
        );

        // Paint rear left trailer fin (mirrored)
        Console.WriteLine("  Painting rear left trailer fin...");
        PaintFin(lat,
            baseX: rearFinX,
            baseY: -rearLateralOffset,
            baseZ: rearRockerZ,
            height: rearFinHeight,
            baseChord: rearFinBase,
            cantAngleDeg: -rearCant,
            isCenter: false
        );

        Console.WriteLine("  Converting fin lattice to voxels...");
        Voxels finVoxels = new(lat);

        Console.WriteLine("  ✓ Quad fin system generated");
        return finVoxels;
    }

    // =========================================================================
    // SINGLE CONFIGURATION
    // =========================================================================

    /// <summary>
    /// Generate single (1-fin) configuration as voxels.
    /// </summary>
    /// <param name="rockerProfile">Rocker profile for positioning fins on the bottom surface</param>
    /// <returns>Voxels containing the single center fin</returns>
    /// <remarks>
    /// <para>
    /// SINGLE FIN POSITIONING:
    /// - One large center fin on the centerline
    /// - Positioned ~80mm trailing margin from the tail
    /// - No cant angle (perfectly vertical)
    /// - Larger and thicker than a thruster center fin for stability
    /// </para>
    /// <para>
    /// SURFBOARD: The original fin setup. A large single fin provides
    /// smooth, drawn-out turns and excellent directional stability.
    /// Traditional for longboards and mid-lengths.
    /// </para>
    /// </remarks>
    public Voxels GenerateSingle(RockerProfile rockerProfile)
    {
        Console.WriteLine("Generating single fin system (1 fin)...");

        Lattice lat = new();

        float scaleFactor = _params.Length / 1828f;

        // Single fin: larger and thicker than thruster center for sole directional control
        float singleFinHeight    = 175f * scaleFactor;  // ~6.9" on a 6'0" (scales up for longboards)
        float singleFinBase      = 165f * scaleFactor;  // ~6.5" chord
        float singleFinThickness = 12f * scaleFactor;   // thicker foil than thruster fins (~8mm)

        // Positioned further from tail than thruster center
        float trailingMargin = 80f;
        float finX    = _params.Length - singleFinBase - trailingMargin;
        float finT    = finX / _params.Length;
        float rockerZ = rockerProfile.GetRockerHeight(finT);

        // Paint the single center fin with thicker foil
        Console.WriteLine("  Painting single center fin...");
        PaintFin(lat,
            baseX: finX,
            baseY: 0f,
            baseZ: rockerZ,
            height: singleFinHeight,
            baseChord: singleFinBase,
            cantAngleDeg: 0f,
            isCenter: true,
            maxThickness: singleFinThickness
        );

        Console.WriteLine("  Converting fin lattice to voxels...");
        Voxels finVoxels = new(lat);

        Console.WriteLine("  ✓ Single fin system generated");
        return finVoxels;
    }

    // =========================================================================
    // FIN PAINTING
    // =========================================================================

    /// <summary>
    /// Paint a single fin using spatial painting.
    /// </summary>
    /// <remarks>
    /// <para>
    /// ALGORITHM: The fin shape is defined by independent leading and trailing
    /// edge curves that converge at a swept-back tip point. This matches how
    /// real fins are shaped:
    ///
    ///     REAL FIN OUTLINE (side view):
    ///
    ///              Tip point (swept back)
    ///              ╱
    ///            ╱   ╲    Trailing edge (relatively straight)
    ///          ╱       │
    ///        ╱         │
    ///      ╱           │
    ///     Leading      │
    ///     edge         │
    ///     (convex      │
    ///      sweep)      │
    ///     ╱             │
    ///    ╱_______________│
    ///    LE base    TE base
    ///    ↔ base chord ↔
    ///
    /// The leading edge sweeps back rapidly (power &lt; 1 for convex curve).
    /// The trailing edge stays relatively straight then curves in at the tip.
    /// Both edges meet at the tip point at X = baseX + finRake.
    /// </para>
    /// </remarks>
    private void PaintFin(
        Lattice lat,
        float baseX,
        float baseY,
        float baseZ,
        float height,
        float baseChord,
        float cantAngleDeg,
        bool isCenter,
        float maxThickness = 0f)
    {
        // Use the provided thickness or fall back to the default
        float effectiveThickness = maxThickness > 0f ? maxThickness : _finThickness;

        float baseSphereRadius = _voxelSize * 2f;

        // Number of height steps (base to tip)
        int heightSteps = (int)(height / (_voxelSize * 3f));
        heightSteps = Math.Max(heightSteps, 20);

        float cantRad = cantAngleDeg * MathF.PI / 180f;

        // The tip point where leading and trailing edges converge.
        // It's at the trailing edge base + a small extra sweep past it.
        // This ensures the tip is BEHIND the trailing edge (toward the tail),
        // matching the swept-back dorsal shape of real fins.
        //
        //   X increases toward tail →
        //
        //   baseX          baseX+chord     tipX
        //   (LE base)      (TE base)       (tip, slightly past TE)
        //     │               │              ╱
        //     │               │            ╱
        //     LE sweeps →→→→→→→→→→→→→→→ ╱
        //                     │ TE stays ╱
        //                     │ straight╱
        //
        float tipX = baseX + baseChord + _finRake;

        for (int i = 0; i <= heightSteps; i++)
        {
            // Normalized height (0 = base, 1 = tip)
            float hT = (float)i / heightSteps;
            float localHeight = hT * height;

            // =====================================================================
            // LEADING EDGE: Concave curve sweeping back from baseX to tipX.
            // Power > 1 keeps the edge near baseX at the base, then sweeps
            // back increasingly toward the tip -- concave when viewed from
            // the water side, matching real fin leading edges.
            // =====================================================================
            float leadingEdgeX = baseX + (tipX - baseX) * MathF.Pow(hT, 1.8f);

            // =====================================================================
            // TRAILING EDGE: Relatively straight, curves in to meet tip.
            // Power > 1 keeps it straight for most of the height, then it
            // curves in to converge with the leading edge at the tip.
            // =====================================================================
            float trailingEdgeBaseX = baseX + baseChord;
            float trailingEdgeX = trailingEdgeBaseX + (tipX - trailingEdgeBaseX) * MathF.Pow(hT, 1.8f);

            // Chord at this height (naturally goes to zero at the tip)
            float chord = trailingEdgeX - leadingEdgeX;

            // Thickness tapers to zero at the tip
            float thicknessTaper = MathF.Pow(1f - hT, 0.7f);
            float maxThick = effectiveThickness * thicknessTaper;

            // Sphere radius decreases near the tip for finer detail
            float sphereRadius = baseSphereRadius * Math.Max(0.4f, 1f - 0.6f * hT);

            // Apply cant angle: height goes downward (-Z), with lateral offset
            float finZ = baseZ - localHeight * MathF.Cos(cantRad);
            float finYOffset = localHeight * MathF.Sin(cantRad);

            // Near the very tip, chord approaches zero --
            // just place a single small sphere to form the point
            if (chord < sphereRadius * 2f)
            {
                lat.AddSphere(
                    new Vector3(tipX, baseY + finYOffset, finZ),
                    sphereRadius
                );
                continue;
            }

            // Number of chord-wise steps
            int chordSteps = (int)(chord / (_voxelSize * 2.5f));
            chordSteps = Math.Max(chordSteps, 10);

            for (int j = 0; j <= chordSteps; j++)
            {
                // Normalized chord position (0 = leading edge, 1 = trailing edge)
                float cT = (float)j / chordSteps;

                // X position along the chord
                float chordX = leadingEdgeX + cT * chord;

                // NACA-like thickness distribution:
                // Thickest at ~30% chord, tapers to near-zero at leading and trailing edges
                float foilThickness = maxThick * MathF.Pow(MathF.Sin(cT * MathF.PI), 0.6f);
                float halfFoil = foilThickness / 2f;

                // For very thin sections near edges, just place a single sphere
                if (halfFoil < sphereRadius * 0.3f)
                {
                    lat.AddSphere(
                        new Vector3(chordX, baseY + finYOffset, finZ),
                        sphereRadius
                    );
                    continue;
                }

                // Paint through the foil thickness (in Y direction)
                int thickSteps = Math.Max((int)(foilThickness / sphereRadius), 1);

                for (int k = 0; k <= thickSteps; k++)
                {
                    float tT = (float)k / thickSteps;
                    float thickOffset = -halfFoil + tT * foilThickness;

                    lat.AddSphere(
                        new Vector3(chordX, baseY + finYOffset + thickOffset, finZ),
                        sphereRadius
                    );
                }
            }
        }
    }
}
