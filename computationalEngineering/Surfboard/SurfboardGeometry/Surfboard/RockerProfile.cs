// =============================================================================
// ROCKER PROFILE - BOTTOM CURVATURE FROM NOSE TO TAIL
// =============================================================================
//
// SURFBOARD CONCEPT: Rocker
// ==========================
//
// Rocker is the curve of the board's bottom from nose to tail when viewed
// from the side. It's one of the most important design elements:
//
//   SIDE VIEW:
//   ──────────
//          Nose Rocker                              Tail Rocker
//         ╱                                                 ╲
//       ╱                                                     ╲
//     ╱                                                         ╲
//   ╱_____________________________Flat Zone_______________________╲
//   ↑                                                              ↑
//   Nose Tip                                                   Tail Tip
//
// ROCKER EFFECTS:
// ---------------
// - MORE NOSE ROCKER: Prevents pearling (nosedive), helps late takeoffs
//   but adds drag by pushing water, reducing paddle speed.
//
// - MORE TAIL ROCKER: Tighter turning, better in hollow/steep waves,
//   but reduces drive and speed down the line.
//
// - FLAT SECTION: The flatter the middle section, the faster the board
//   will be in a straight line. Performance boards have continuous rocker
//   while speed boards have a flatter center.
//
// ROCKER MEASUREMENT:
// -------------------
// Measured by laying the board on a flat surface and measuring the
// distance from the surface to the nose tip (nose rocker) and tail
// tip (tail rocker). The flat surface touches the board at its lowest
// point (usually around 60% from nose).
//
// =============================================================================

namespace SurfboardGeometry.Surfboard;

/// <summary>
/// Computes the rocker profile (Z-offset of the bottom surface along the board length).
/// </summary>
/// <remarks>
/// Returns the vertical offset (Z) of the bottom surface at each longitudinal
/// position. The lowest point of the board is at Z = 0 (the "flat spot"),
/// with the nose and tail curving upward from there.
///
/// Uses a composite curve: quadratic rise at the nose and tail,
/// with a smooth transition through the flat center section.
/// </remarks>
public class RockerProfile
{
    private readonly SurfboardParameters _params;

    // Normalized position of the lowest point (flat spot)
    // Typically around 55-65% from the nose
    private readonly float _flatSpotT;

    public RockerProfile(SurfboardParameters parameters)
    {
        _params = parameters;

        // The flat spot (lowest point) is typically slightly behind center.
        // On a shortboard it's around 60% from the nose.
        // On a longboard with forward wide point, closer to 55%.
        _flatSpotT = 0.60f;
    }

    /// <summary>
    /// Get the Z-offset (rocker height) at a normalized longitudinal position.
    /// </summary>
    /// <param name="t">Normalized position: 0.0 = nose tip, 1.0 = tail tip</param>
    /// <returns>Z-offset in mm (0 = flat spot, positive = upward)</returns>
    /// <remarks>
    /// MATH: The rocker curve uses power functions with different exponents
    /// for the nose and tail regions:
    ///
    ///   Nose region (t &lt; flatSpot):
    ///     z = noseRocker * (1 - t/flatSpot)^noseExponent
    ///
    ///   Tail region (t &gt; flatSpot):
    ///     z = tailRocker * ((t - flatSpot)/(1 - flatSpot))^tailExponent
    ///
    /// The exponent controls the curve shape:
    ///   - exponent = 1.0: linear (conical)
    ///   - exponent = 2.0: quadratic (gentle curve near flat, steep at tips)
    ///   - exponent = 2.5: more aggressive curve (flatter center)
    ///
    /// Nose typically uses a higher exponent (2.5) for a flatter middle
    /// and steep nose kick. Tail uses ~2.0 for a smoother transition.
    /// </remarks>
    public float GetRockerHeight(float t)
    {
        t = Math.Clamp(t, 0f, 1f);

        if (t <= _flatSpotT)
        {
            // Nose region: rises from flat spot toward nose tip
            // Normalized distance from flat spot toward nose (0 at flat, 1 at nose)
            float noseT = 1f - (t / _flatSpotT);

            // Power curve: steep at the nose tip, flat near the center
            float noseExponent = 2.5f;
            return _params.NoseRocker * MathF.Pow(noseT, noseExponent);
        }
        else
        {
            // Tail region: rises from flat spot toward tail tip
            // Normalized distance from flat spot toward tail (0 at flat, 1 at tail)
            float tailT = (t - _flatSpotT) / (1f - _flatSpotT);

            // Power curve: smoother than nose, gradual kick at the tail
            float tailExponent = 2.0f;
            return _params.TailRocker * MathF.Pow(tailT, tailExponent);
        }
    }
}
