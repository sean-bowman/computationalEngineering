// =============================================================================
// CROSS-SECTION - DECK AND BOTTOM PROFILE AT EACH LONGITUDINAL STATION
// =============================================================================
//
// SURFBOARD CONCEPT: Cross-Sectional Shape
// ==========================================
//
// The cross-section defines the board's shape when sliced perpendicular
// to the length. It varies along the board:
//
//   THICK CENTER (at wide point):      THIN NOSE/TAIL:
//   ─────────────────────────          ──────────────────
//         Deck Crown                     Thin dome
//      ╱‾‾‾‾‾‾‾‾‾‾‾‾╲                  ╱‾‾‾‾‾╲
//     │    Volume!     │               ╱         ╲
//     │                │              │             │
//      ╲___________  ╱                 ╲___________╱
//       Concave                           Flat/Round
//
// CROSS-SECTION ELEMENTS:
// -----------------------
//
// 1. DECK (top surface): Convex dome that adds volume while keeping
//    rails thin. The crown height decreases toward nose and tail.
//
// 2. BOTTOM (bottom surface): Can be flat, single concave, or double
//    concave. Concave channels water flow for speed and lift.
//
// 3. RAILS (edges): The transition from deck to bottom. Ranges from
//    "soft" (round, 50/50 distribution) to "hard" (sharp, down-turned).
//    Hard rails release water cleanly; soft rails are forgiving.
//
//   Rail profiles:
//   ──────────────
//   Soft (round):     Medium (tucked):    Hard (sharp):
//     ╱‾╲               ╱‾╲                 ╱‾╲
//    │   │              │   │               │   │
//    │   │               ╲  │                ╲  │
//     ╲_╱                 ╲_│                 ╲_│
//
// =============================================================================

namespace SurfboardGeometry.Surfboard;

/// <summary>
/// Computes the cross-sectional profile (deck and bottom heights) at a given
/// longitudinal station.
/// </summary>
/// <remarks>
/// For a given position along the board length and across the width,
/// returns the deck height (top) and bottom height relative to the
/// board's center plane. The thickness envelope varies along the board,
/// being thickest at the center and thinning toward nose and tail.
/// </remarks>
public class CrossSection
{
    private readonly SurfboardParameters _params;

    public CrossSection(SurfboardParameters parameters)
    {
        _params = parameters;
    }

    /// <summary>
    /// Get the deck height (top surface Z) at a given longitudinal and lateral position.
    /// </summary>
    /// <param name="t">Normalized longitudinal position (0 = nose, 1 = tail)</param>
    /// <param name="lateralFraction">
    /// Normalized lateral position (0 = centerline, 1 = rail edge).
    /// </param>
    /// <returns>Deck height in mm above the rocker datum (positive = up)</returns>
    /// <remarks>
    /// The deck profile is a cosine dome:
    ///   z_deck = halfThickness(t) + crown(t) * cos(lateralFraction * pi/2)
    ///
    /// The crown (dome height) peaks at the center and falls off toward
    /// the edges (rails), creating a crowned deck that sheds water and
    /// adds volume without increasing rail thickness.
    /// </remarks>
    public float GetDeckHeight(float t, float lateralFraction)
    {
        t = Math.Clamp(t, 0f, 1f);
        lateralFraction = Math.Clamp(lateralFraction, 0f, 1f);

        float thickness = GetLocalThickness(t);
        float halfThick = thickness / 2f;

        // Crown contribution: cosine dome that peaks at centerline
        // and drops to zero at the rail edge
        float crownAtStation = GetLocalCrown(t);
        float crownContribution = crownAtStation * MathF.Cos(lateralFraction * MathF.PI / 2f);

        return halfThick + crownContribution;
    }

    /// <summary>
    /// Get the bottom height (bottom surface Z) at a given longitudinal and lateral position.
    /// </summary>
    /// <param name="t">Normalized longitudinal position (0 = nose, 1 = tail)</param>
    /// <param name="lateralFraction">
    /// Normalized lateral position (0 = centerline, 1 = rail edge).
    /// </param>
    /// <returns>Bottom height in mm below the rocker datum (negative = down)</returns>
    /// <remarks>
    /// The bottom surface includes optional single concave:
    ///   z_bottom = -halfThickness(t) - concave(t) * cos(lateralFraction * pi)
    ///
    /// The concave creates a channel along the centerline:
    ///
    ///   Flat bottom:           Single concave:
    ///   ───────────            ─────╲  ╱─────
    ///                                ╲╱
    ///
    /// The concave depth is greatest in the middle of the board and
    /// fades to zero at the nose and tail.
    /// </remarks>
    public float GetBottomHeight(float t, float lateralFraction)
    {
        t = Math.Clamp(t, 0f, 1f);
        lateralFraction = Math.Clamp(lateralFraction, 0f, 1f);

        float thickness = GetLocalThickness(t);
        float halfThick = thickness / 2f;

        // Concave contribution: cosine valley at the centerline
        // cos(lateralFraction * PI) ranges from 1 (center) to -1 (edge)
        // We want concave at center, so subtract when positive
        float concaveAtStation = GetLocalConcave(t);
        float concaveContribution = concaveAtStation
            * MathF.Pow(MathF.Cos(lateralFraction * MathF.PI / 2f), 2f);

        return -(halfThick + concaveContribution);
    }

    // =========================================================================
    // THICKNESS DISTRIBUTION ALONG THE BOARD
    // =========================================================================

    /// <summary>
    /// Get the local board thickness at a normalized longitudinal position.
    /// </summary>
    /// <remarks>
    /// SURFBOARD: The thickness distribution follows a smooth curve:
    /// - Thin at the nose (~40% of max thickness)
    /// - Maximum at approximately 40% from nose (just ahead of center)
    /// - Gradually thinning toward the tail (~50% of max)
    ///
    /// This uses a cosine-based envelope centered on the thickest point.
    /// </remarks>
    public float GetLocalThickness(float t)
    {
        // Thickest point is typically ~40% from the nose
        float thickestPointT = 0.40f;

        // Minimum thickness fractions at nose and tail
        float noseThicknessFraction = 0.35f;  // 35% of max at nose tip
        float tailThicknessFraction = 0.45f;  // 45% of max at tail tip

        float fraction;
        if (t <= thickestPointT)
        {
            // Nose side: interpolate from nose fraction up to 1.0 at thickest point
            float localT = t / thickestPointT;
            fraction = CosineInterpolate(noseThicknessFraction, 1f, localT);
        }
        else
        {
            // Tail side: interpolate from 1.0 at thickest point down to tail fraction
            float localT = (t - thickestPointT) / (1f - thickestPointT);
            fraction = CosineInterpolate(1f, tailThicknessFraction, localT);
        }

        return _params.MaxThickness * fraction;
    }

    /// <summary>
    /// Get the deck crown height at a normalized longitudinal position.
    /// </summary>
    /// <remarks>
    /// Crown is greatest at the thick center and fades toward nose and tail.
    /// </remarks>
    private float GetLocalCrown(float t)
    {
        // Crown follows a smooth bell curve, peaking at ~45% from nose
        float peakT = 0.45f;
        float distance = MathF.Abs(t - peakT);
        float maxDistance = MathF.Max(peakT, 1f - peakT);
        float normalized = distance / maxDistance;

        // Cosine falloff from peak
        return _params.DeckCrown * MathF.Pow(MathF.Cos(normalized * MathF.PI / 2f), 1.5f);
    }

    /// <summary>
    /// Get the bottom concave depth at a normalized longitudinal position.
    /// </summary>
    /// <remarks>
    /// SURFBOARD: Concave is strongest in the center of the board (under
    /// the surfer's stance) and fades to zero at nose and tail.
    /// The transition zone is the most important section for speed.
    /// </remarks>
    private float GetLocalConcave(float t)
    {
        // Concave is strongest between 30% and 70% of board length
        float startFade = 0.25f;
        float endFade = 0.75f;

        if (t < startFade)
        {
            float localT = t / startFade;
            return _params.BottomConcave * CosineInterpolate(0f, 1f, localT);
        }
        else if (t > endFade)
        {
            float localT = (t - endFade) / (1f - endFade);
            return _params.BottomConcave * CosineInterpolate(1f, 0f, localT);
        }
        else
        {
            return _params.BottomConcave;
        }
    }

    /// <summary>
    /// Cosine interpolation between two values for smooth transitions.
    /// </summary>
    private static float CosineInterpolate(float a, float b, float t)
    {
        float blend = (1f - MathF.Cos(t * MathF.PI)) / 2f;
        return a + (b - a) * blend;
    }
}
