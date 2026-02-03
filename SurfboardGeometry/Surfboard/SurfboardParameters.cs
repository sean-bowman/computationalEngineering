// =============================================================================
// SURFBOARD PARAMETERS - DATA CLASS FOR PARAMETRIC SURFBOARD DIMENSIONS
// =============================================================================
//
// SURFBOARD DESIGN CONCEPT: Parametric Board Geometry
// ====================================================
//
// A surfboard's shape is defined by several key profiles that work together
// to determine how the board interacts with a wave:
//
//   TOP VIEW (PLANFORM / OUTLINE):
//   ──────────────────────────────
//                    Nose
//                   ╱    ╲
//                 ╱        ╲
//               ╱            ╲
//             ╱   Nose Width   ╲          12" from tip
//           ╱                    ╲
//          │                      │
//          │    Wide Point         │       Max width
//          │                      │
//           ╲                    ╱
//            ╲   Tail Width    ╱
//              ╲             ╱
//                ╲_________╱
//                   Tail
//
//   SIDE VIEW (ROCKER PROFILE):
//   ───────────────────────────
//          Nose Rocker
//         ╱
//       ╱
//     ╱
//   ╱___________________________________╲  Tail Rocker
//        Flat section (planing area)      ╲
//
//   CROSS-SECTION (at wide point):
//   ──────────────────────────────
//            Deck Crown
//         ╱‾‾‾‾‾‾‾‾‾‾‾‾╲
//        │                │    Rail
//         ╲_____________ ╱
//          Bottom Concave
//
// DESIGN TRADE-OFFS:
// ------------------
// - More rocker = more maneuverability, less speed
// - Wider board = more stability, less performance
// - More concave = more speed (channeled water flow)
// - Thicker board = more float/paddle, less duck-diving
//
// =============================================================================

namespace SurfboardGeometry.Surfboard;

/// <summary>
/// Parameters defining a parametric surfboard shape.
/// All dimensions in millimeters.
/// </summary>
/// <remarks>
/// <para>
/// SURFBOARD DESIGN: These parameters capture the essential dimensions
/// that define a surfboard's shape. Real shapers use these same measurements
/// (in imperial units) to describe and reproduce board designs.
/// </para>
/// <para>
/// Uses init-only properties for immutability after construction.
/// Factory properties provide common board type presets.
/// </para>
/// </remarks>
public class SurfboardParameters
{
    // =========================================================================
    // PRIMARY DIMENSIONS
    // =========================================================================

    /// <summary>Total board length from nose to tail [mm]</summary>
    /// <remarks>
    /// SURFBOARD: Length is the primary sizing dimension.
    /// Shorter = more maneuverable; longer = more speed and stability.
    /// Typical range: 5'4" (1625mm) to 10'0" (3048mm).
    /// </remarks>
    public float Length { get; init; } = 1828f;         // 6'0"

    /// <summary>Maximum width at the wide point [mm]</summary>
    /// <remarks>
    /// SURFBOARD: Width affects stability and planing speed.
    /// Wider boards plane earlier and are more stable.
    /// Typical range: 18" (457mm) to 23" (584mm).
    /// </remarks>
    public float MaxWidth { get; init; } = 495f;        // 19.5"

    /// <summary>Maximum thickness at the thickest point [mm]</summary>
    /// <remarks>
    /// SURFBOARD: Thickness determines volume (buoyancy).
    /// More thickness = easier paddling but harder to duck dive.
    /// Typical range: 2.0" (51mm) to 3.25" (83mm).
    /// </remarks>
    public float MaxThickness { get; init; } = 62f;      // 2.44"

    // =========================================================================
    // OUTLINE CONTROL
    // =========================================================================

    /// <summary>Width measured 305mm (12") from the nose tip [mm]</summary>
    /// <remarks>
    /// SURFBOARD: Nose width affects paddling entry into waves and
    /// how easily the board catches waves. Wider = more catch.
    /// </remarks>
    public float NoseWidth { get; init; } = 300f;        // ~11.8"

    /// <summary>Width at the tail [mm]</summary>
    /// <remarks>
    /// SURFBOARD: Tail width affects release and hold.
    /// Wider tail = more speed/drive; narrower = more hold in steep waves.
    /// </remarks>
    public float TailWidth { get; init; } = 380f;        // ~15"

    /// <summary>
    /// Wide point offset from board center [mm].
    /// Negative = toward nose, positive = toward tail.
    /// </summary>
    /// <remarks>
    /// SURFBOARD: Forward wide point = more paddle power and wave catching.
    /// Centered or aft wide point = more performance and pivot turning.
    /// </remarks>
    public float WidePointOffset { get; init; } = -25f;

    // =========================================================================
    // ROCKER PROFILE
    // =========================================================================

    /// <summary>Nose rocker height measured at the tip [mm]</summary>
    /// <remarks>
    /// SURFBOARD: More nose rocker prevents pearling (nose diving) and
    /// aids in late takeoffs, but reduces paddle speed.
    /// Typical range: 2.5" (64mm) to 7" (178mm).
    /// </remarks>
    public float NoseRocker { get; init; } = 120f;       // ~4.7"

    /// <summary>Tail rocker height measured at the tail tip [mm]</summary>
    /// <remarks>
    /// SURFBOARD: More tail rocker = tighter turning radius and better
    /// performance in hollow waves. Less = faster down the line.
    /// Typical range: 1" (25mm) to 3.5" (89mm).
    /// </remarks>
    public float TailRocker { get; init; } = 40f;        // ~1.6"

    // =========================================================================
    // CROSS-SECTION PROFILE
    // =========================================================================

    /// <summary>Deck dome height at center (convexity) [mm]</summary>
    /// <remarks>
    /// SURFBOARD: Dome/crown adds volume while keeping rails thinner.
    /// More dome = more volume without adding rail thickness.
    /// </remarks>
    public float DeckCrown { get; init; } = 8f;

    /// <summary>Bottom single concave depth [mm]</summary>
    /// <remarks>
    /// SURFBOARD: Concave channels water flow under the board, creating
    /// lift and speed. Deeper concave = more acceleration but can feel
    /// stiff in turns. 0 = flat bottom.
    /// </remarks>
    public float BottomConcave { get; init; } = 2f;

    /// <summary>Rail edge radius at the thinnest point [mm]</summary>
    /// <remarks>
    /// SURFBOARD: Sharp (hard) rails release water cleanly for speed.
    /// Round (soft) rails are more forgiving and work better in small surf.
    /// </remarks>
    public float RailRadius { get; init; } = 12f;

    // =========================================================================
    // FIN CONFIGURATION
    // =========================================================================

    /// <summary>Default fin configuration for this board type</summary>
    /// <remarks>
    /// SURFBOARD: Each board type has a traditional fin setup.
    /// Shortboards use thruster (3 fins), fish boards use twin (2 fins),
    /// and longboards use a single fin. Users can override via CLI arguments.
    /// </remarks>
    public FinConfiguration DefaultFinConfiguration { get; init; } = FinConfiguration.Thruster;

    // =========================================================================
    // COMPUTED PROPERTIES
    // =========================================================================

    /// <summary>Half the maximum width [mm]</summary>
    public float HalfWidth => MaxWidth / 2f;

    /// <summary>Half the total length [mm]</summary>
    public float HalfLength => Length / 2f;

    /// <summary>X position of the wide point relative to nose [mm]</summary>
    /// <remarks>
    /// The wide point position along the board, measured from the nose.
    /// WidePointOffset shifts it from center: negative = toward nose.
    /// </remarks>
    public float WidePointX => HalfLength + WidePointOffset;

    /// <summary>Half the maximum thickness [mm]</summary>
    public float HalfThickness => MaxThickness / 2f;

    /// <summary>Approximate board volume [liters] using ellipsoid estimation</summary>
    /// <remarks>
    /// SURFBOARD: Volume determines buoyancy and paddle power.
    /// A rough rule: beginner boards ~40-60L, performance shortboards ~22-30L.
    /// This is an approximation; actual volume comes from the voxel model.
    /// V ≈ (4/3) * pi * (L/2) * (W/2) * (T/2) * packing_factor
    /// </remarks>
    public float ApproxVolumeLiters =>
        (4f / 3f) * MathF.PI * (Length / 2f) * (MaxWidth / 2f) * (MaxThickness / 2f)
        * 0.35f     // packing factor (surfboard fills ~35% of bounding ellipsoid)
        / 1_000_000f; // mm^3 to liters

    /// <summary>Planform area estimate [mm^2]</summary>
    public float ApproxPlanformArea =>
        Length * MaxWidth * 0.65f; // surfboard fills ~65% of bounding rectangle

    // =========================================================================
    // DISPLAY HELPERS
    // =========================================================================

    /// <summary>Convert mm to feet-inches string for display</summary>
    private static string MmToFeetInches(float mm)
    {
        float totalInches = mm / 25.4f;
        int feet = (int)MathF.Floor(totalInches / 12f);
        float inches = totalInches - feet * 12f;
        // Handle rounding: if inches rounds to 12, carry over to feet
        if (inches >= 11.95f)
        {
            feet++;
            inches = 0f;
        }
        return $"{feet}'{inches:F1}\"";
    }

    /// <summary>Convert mm to inches string for display</summary>
    private static string MmToInches(float mm)
    {
        return $"{mm / 25.4f:F2}\"";
    }

    // =========================================================================
    // FACTORY PRESETS
    // =========================================================================

    /// <summary>
    /// Standard performance shortboard (6'0" x 19.5" x 2.44").
    /// </summary>
    /// <remarks>
    /// SURFBOARD: The classic high-performance shortboard. Designed for
    /// experienced surfers in good waves. Tight rocker, narrow outline,
    /// moderate volume. Excels in steep, hollow waves.
    /// </remarks>
    public static SurfboardParameters Shortboard => new();

    /// <summary>
    /// Classic longboard (9'0" x 22.5" x 3").
    /// </summary>
    /// <remarks>
    /// SURFBOARD: Traditional longboard for nose riding and glide.
    /// Low rocker for speed, wide outline for stability, thick rails
    /// for smooth turns. Best in small to medium clean waves.
    /// </remarks>
    public static SurfboardParameters Longboard => new()
    {
        Length = 2743f,           // 9'0"
        MaxWidth = 570f,         // 22.5"
        MaxThickness = 75f,      // 3.0"
        NoseWidth = 430f,        // ~17" wide nose for nose riding
        TailWidth = 370f,        // ~14.5"
        WidePointOffset = 0f,    // centered wide point
        NoseRocker = 180f,       // ~7" generous nose rocker
        TailRocker = 25f,        // ~1" minimal tail rocker for speed
        DeckCrown = 10f,
        BottomConcave = 1f,      // subtle concave
        RailRadius = 18f,        // soft, round rails
        DefaultFinConfiguration = FinConfiguration.Single
    };

    /// <summary>
    /// Retro fish (5'6" x 21" x 2.56").
    /// </summary>
    /// <remarks>
    /// SURFBOARD: Wide, flat, and fast. The fish was designed for speed
    /// in small waves. Low rocker, wide nose and tail, extra volume
    /// packed into a short frame. Swallow tail typical.
    /// </remarks>
    public static SurfboardParameters Fish => new()
    {
        Length = 1676f,           // 5'6"
        MaxWidth = 533f,         // 21"
        MaxThickness = 65f,      // 2.56"
        NoseWidth = 380f,        // ~15" wide nose
        TailWidth = 420f,        // ~16.5" wide tail (fish characteristic)
        WidePointOffset = -50f,  // forward wide point for paddle
        NoseRocker = 80f,        // ~3.1" minimal nose rocker
        TailRocker = 30f,        // ~1.2" low tail rocker
        DeckCrown = 6f,
        BottomConcave = 3f,      // deeper concave for speed
        RailRadius = 10f,        // medium-hard rails
        DefaultFinConfiguration = FinConfiguration.Twin
    };

    // =========================================================================
    // UTILITY METHODS
    // =========================================================================

    /// <summary>
    /// Print surfboard parameters to console in a formatted table.
    /// </summary>
    public void PrintSummary()
    {
        Console.WriteLine("╔══════════════════════════════════════════════════════════╗");
        Console.WriteLine("║              SURFBOARD PARAMETERS SUMMARY               ║");
        Console.WriteLine("╠══════════════════════════════════════════════════════════╣");
        Console.WriteLine($"║  Length:              {Length,8:F1} mm  ({MmToFeetInches(Length),8})   ║");
        Console.WriteLine($"║  Max Width:           {MaxWidth,8:F1} mm  ({MmToInches(MaxWidth),8})   ║");
        Console.WriteLine($"║  Max Thickness:       {MaxThickness,8:F1} mm  ({MmToInches(MaxThickness),8})   ║");
        Console.WriteLine("╠══════════════════════════════════════════════════════════╣");
        Console.WriteLine($"║  Nose Width (12\"):    {NoseWidth,8:F1} mm  ({MmToInches(NoseWidth),8})   ║");
        Console.WriteLine($"║  Tail Width:          {TailWidth,8:F1} mm  ({MmToInches(TailWidth),8})   ║");
        Console.WriteLine($"║  Wide Point Offset:   {WidePointOffset,8:F1} mm                    ║");
        Console.WriteLine("╠══════════════════════════════════════════════════════════╣");
        Console.WriteLine($"║  Nose Rocker:         {NoseRocker,8:F1} mm  ({MmToInches(NoseRocker),8})   ║");
        Console.WriteLine($"║  Tail Rocker:         {TailRocker,8:F1} mm  ({MmToInches(TailRocker),8})   ║");
        Console.WriteLine("╠══════════════════════════════════════════════════════════╣");
        Console.WriteLine($"║  Deck Crown:          {DeckCrown,8:F1} mm                    ║");
        Console.WriteLine($"║  Bottom Concave:      {BottomConcave,8:F1} mm                    ║");
        Console.WriteLine($"║  Rail Radius:         {RailRadius,8:F1} mm                    ║");
        Console.WriteLine("╠══════════════════════════════════════════════════════════╣");
        Console.WriteLine($"║  Approx Volume:       {ApproxVolumeLiters,8:F1} L                     ║");
        Console.WriteLine("╚══════════════════════════════════════════════════════════╝");
    }
}
