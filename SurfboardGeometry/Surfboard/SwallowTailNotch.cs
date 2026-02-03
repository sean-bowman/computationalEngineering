// =============================================================================
// SWALLOW TAIL NOTCH - V-SHAPED CENTER NOTCH FOR FISH TAILS
// =============================================================================
//
// SURFBOARD CONCEPT: Swallow (Fish) Tail
// =======================================
//
// A swallow tail splits the tail into two lobes via a V-shaped center notch.
// This increases the effective rail length and allows water to release
// cleanly from each lobe independently, improving speed and looseness.
//
//   STANDARD TAIL:              SWALLOW TAIL:
//   ──────────────              ──────────────
//        │                           │
//        │                           │
//         ╲                        ╱   ╲
//          ╲                      ╱     ╲
//           ╲                    │       │
//            ╲                   │ Notch │
//             ╲                  │       │
//              ╲                 ╱╲     ╱╲
//               •               • ╲   ╱ •
//           (point)          (lobe)  V  (lobe)
//
// IMPLEMENTATION:
// ───────────────
// The notch is carved from a wide, blunt tail body using PicoGK's
// boolean subtraction. A V-shaped wedge solid is defined via a signed
// distance function (IBoundedImplicit) and subtracted from the board body.
//
// This follows the same pattern as fin attachment (separate geometry +
// boolean operation), keeping the outline and body painting code unchanged.
//
// =============================================================================

using System.Numerics;
using PicoGK;

namespace SurfboardGeometry.Surfboard;

/// <summary>
/// Signed distance function defining a V-shaped wedge for the swallow tail notch.
/// </summary>
/// <remarks>
/// <para>
/// The wedge is a V-shaped prism centered on the stringer line (Y=0),
/// with its apex pointing forward (toward the nose) and widening toward
/// the tail tip.
/// </para>
/// <para>
/// MATH: At any longitudinal position X in the notch region, the notch
/// half-width is linearly interpolated from 0 (at the apex) to the
/// maximum notch half-width (at the tail tip):
///
///   localHalfWidth = notchHalfWidth * (x - apexX) / notchDepth
///
/// A point is inside the wedge if |Y| &lt; localHalfWidth and X is
/// between the apex and the tail tip. The signed distance is negative
/// inside, positive outside.
/// </para>
/// </remarks>
internal class SwallowNotchImplicit : IBoundedImplicit
{
    private readonly float _apexX;
    private readonly float _tailX;
    private readonly float _notchHalfWidth;
    private readonly float _notchDepth;
    private readonly BBox3 _bounds;

    /// <summary>
    /// Create the implicit function for a swallow tail V-notch.
    /// </summary>
    /// <param name="parameters">Board parameters with swallow tail dimensions</param>
    /// <param name="rockerZAtTail">Rocker Z-offset at the tail tip for bounding box</param>
    public SwallowNotchImplicit(SurfboardParameters parameters, float rockerZAtTail)
    {
        _tailX = parameters.Length;
        _notchDepth = parameters.SwallowNotchDepth;
        _apexX = _tailX - _notchDepth;
        _notchHalfWidth = parameters.SwallowNotchHalfWidth;

        // Bounding box: cover the wedge region with generous Z margins
        // to ensure the notch cuts through the full board thickness
        float zMargin = parameters.MaxThickness + 50f;
        _bounds = new BBox3(
            _apexX - 5f,
            -_notchHalfWidth - 5f,
            rockerZAtTail - zMargin,
            _tailX + 5f,
            _notchHalfWidth + 5f,
            rockerZAtTail + zMargin
        );
    }

    /// <summary>Bounding box of the notch region.</summary>
    public BBox3 oBounds => _bounds;

    /// <summary>
    /// Signed distance to the V-shaped wedge boundary.
    /// </summary>
    /// <param name="vec">Point to evaluate</param>
    /// <returns>Negative inside the wedge, positive outside</returns>
    public float fSignedDistance(in Vector3 vec)
    {
        float x = vec.X;
        float y = vec.Y;

        // Behind the apex (toward nose): outside the notch
        if (x < _apexX)
        {
            // Distance to the apex plane
            return _apexX - x;
        }

        // Forward of apex: compute the notch half-width at this X
        float fraction = (x - _apexX) / _notchDepth;
        fraction = MathF.Min(fraction, 1f);
        float localHalfWidth = _notchHalfWidth * fraction;

        // Signed distance to the angled sides of the V
        // Negative = inside the wedge, positive = outside
        return MathF.Abs(y) - localHalfWidth;
    }
}

/// <summary>
/// Generates the swallow tail notch as a voxel solid for boolean subtraction.
/// </summary>
/// <remarks>
/// Follows the same pattern as FinSystem: takes parameters and voxel size,
/// has a Generate method that returns Voxels. The returned voxels represent
/// the V-shaped wedge that should be subtracted from the board body.
/// </remarks>
public class SwallowTailNotch
{
    private readonly SurfboardParameters _params;
    private readonly float _voxelSize;

    /// <summary>
    /// Create a swallow tail notch generator.
    /// </summary>
    /// <param name="parameters">Board parameters with swallow tail dimensions</param>
    /// <param name="voxelSize">Voxel resolution in mm</param>
    public SwallowTailNotch(SurfboardParameters parameters, float voxelSize)
    {
        _params = parameters;
        _voxelSize = voxelSize;
    }

    /// <summary>
    /// Generate the V-shaped notch as voxels.
    /// </summary>
    /// <param name="rocker">Rocker profile for determining Z position at the tail</param>
    /// <returns>Voxels representing the wedge to be subtracted from the body</returns>
    public Voxels Generate(RockerProfile rocker)
    {
        Console.WriteLine("Generating swallow tail notch...");

        // Get the rocker Z at the tail tip for bounding box positioning
        float tailRockerZ = rocker.GetRockerHeight(1.0f);

        var notchImplicit = new SwallowNotchImplicit(_params, tailRockerZ);
        Voxels voxNotch = new(notchImplicit);

        Console.WriteLine($"  Notch depth:      {_params.SwallowNotchDepth:F0} mm");
        Console.WriteLine($"  Notch half-width: {_params.SwallowNotchHalfWidth:F0} mm at tail tip");
        Console.WriteLine("  ✓ Swallow tail notch generated");

        return voxNotch;
    }
}
