// =============================================================================
// OUTLINE - SURFBOARD PLANFORM SHAPE (TOP VIEW)
// =============================================================================
//
// SURFBOARD CONCEPT: Planform Outline
// ====================================
//
// The outline defines the board's shape when viewed from above. It's the
// most visually distinctive aspect of a surfboard and directly affects:
//
// - WAVE CATCHING: Wider nose/outline = more planing area = easier paddling
// - TURNING: Narrower outline = tighter turning arc
// - HOLD: Curved rail line (outline curvature) grips the wave face
// - SPEED: Straighter outline = less drag at speed
//
// OUTLINE ANATOMY:
// ----------------
//
//   Nose Tip (0.0)
//        │
//   Nose Width ──── measured 12" from tip
//        │
//   Wide Point ──── maximum width (offset from center)
//        │
//   Tail Width ──── width near the tail
//        │
//   Tail Tip (1.0)
//
// This class maps a normalized position along the board (0 = nose, 1 = tail)
// to the half-width at that position. The board is symmetric about its
// centerline, so only the half-width is needed.
//
// CURVE GENERATION:
// -----------------
// Uses a two-piece power curve (nose side + tail side of the wide point).
// Power curves with calibrated exponents produce smooth, convex outlines
// with no inflection points -- matching how real surfboards are shaped.
//
// The NoseWidth and TailWidth parameters are used to calibrate the curve
// exponents so the outline passes through those measured widths at the
// correct stations. This is the same approach a shaper uses: specify a
// few key measurements and let the curve flow smoothly between them.
//
// =============================================================================

namespace SurfboardGeometry.Surfboard;

/// <summary>
/// Computes the surfboard planform outline (half-width at each longitudinal station).
/// </summary>
/// <remarks>
/// <para>
/// The outline is symmetric about the centerline. Given a normalized position
/// along the board (0 = nose tip, 1 = tail tip), returns the half-width.
/// </para>
/// <para>
/// Uses a two-piece power curve split at the wide point. Each side is a single
/// smooth, convex curve with no inflection points. The exponents are calibrated
/// from NoseWidth and TailWidth so the outline passes through those measurements.
/// </para>
/// </remarks>
public class Outline
{
    private readonly SurfboardParameters _params;

    // Normalized position of the wide point (0 = nose, 1 = tail)
    private readonly float _widePointT;

    // Power curve exponents (calibrated from NoseWidth and TailWidth)
    private readonly float _noseExponent;
    private readonly float _tailExponent;

    // Half-widths at key stations
    private readonly float _maxHalfWidth;
    private readonly float _tailTipHalf;

    public Outline(SurfboardParameters parameters)
    {
        _params = parameters;

        _widePointT = _params.WidePointX / _params.Length;
        _maxHalfWidth = _params.HalfWidth;
        _tailTipHalf = parameters.TailTipHalfWidth;

        // =====================================================================
        // CALIBRATE NOSE EXPONENT
        // =====================================================================
        // The nose side power curve is: halfWidth = maxHalfWidth * pow(localT, n)
        // where localT goes from 0 (nose tip) to 1 (wide point).
        //
        // We want the curve to pass through NoseWidth/2 at the nose measurement
        // station (12" = 305mm from the nose tip). Solving for the exponent:
        //
        //   noseWidth/2 = maxHalfWidth * pow(noseStationT / widePointT, n)
        //   n = log(noseWidth/2 / maxHalfWidth) / log(noseStationT / widePointT)
        float noseStationT = 305f / _params.Length;
        float noseRatio = (parameters.NoseWidth / 2f) / _maxHalfWidth;
        float nosePositionRatio = noseStationT / _widePointT;

        // Guard against degenerate cases
        if (noseRatio > 0f && noseRatio < 1f && nosePositionRatio > 0f && nosePositionRatio < 1f)
        {
            _noseExponent = MathF.Log(noseRatio) / MathF.Log(nosePositionRatio);
        }
        else
        {
            _noseExponent = 0.5f; // fallback: square root curve
        }

        // Clamp to reasonable range (0.3 = very full nose, 1.5 = very narrow)
        _noseExponent = Math.Clamp(_noseExponent, 0.3f, 1.5f);

        // =====================================================================
        // CALIBRATE TAIL EXPONENT
        // =====================================================================
        // The tail side power curve is:
        //   halfWidth = tailTipHalf + (maxHalfWidth - tailTipHalf) * pow(1 - localT, n)
        // where localT goes from 0 (wide point) to 1 (tail tip).
        //
        // We want the curve to pass through TailWidth/2 at the tail measurement
        // station (12" = 305mm from the tail tip).
        float tailStationFromTail = 305f / _params.Length;
        float tailStationT = 1f - tailStationFromTail; // normalized position from nose
        float tailLocalT = (tailStationT - _widePointT) / (1f - _widePointT);

        float tailTargetHalf = parameters.TailWidth / 2f;
        // Solve: tailTargetHalf = tailTipHalf + (maxHalfWidth - tailTipHalf) * pow(1 - tailLocalT, n)
        float tailFraction = (tailTargetHalf - _tailTipHalf) / (_maxHalfWidth - _tailTipHalf);
        float tailBlendArg = 1f - tailLocalT;

        if (tailFraction > 0f && tailFraction < 1f && tailBlendArg > 0f && tailBlendArg < 1f)
        {
            _tailExponent = MathF.Log(tailFraction) / MathF.Log(tailBlendArg);
        }
        else
        {
            _tailExponent = 0.6f; // fallback
        }

        // Clamp to reasonable range (0.3 = very full tail, 2.0 = very drawn-in)
        _tailExponent = Math.Clamp(_tailExponent, 0.3f, 2.0f);
    }

    /// <summary>
    /// Get the half-width of the board at a normalized longitudinal position.
    /// </summary>
    /// <param name="t">Normalized position: 0.0 = nose tip, 1.0 = tail tip</param>
    /// <returns>Half-width in mm at the given position</returns>
    /// <remarks>
    /// <para>
    /// Two-piece power curve split at the wide point:
    /// </para>
    /// <para>
    /// NOSE SIDE (t = 0 to widePointT):
    ///   halfWidth = maxHalfWidth * pow(t / widePointT, noseExponent)
    ///   Smoothly expands from 0 at the nose tip to maxHalfWidth.
    /// </para>
    /// <para>
    /// TAIL SIDE (t = widePointT to 1):
    ///   halfWidth = tailTipHalf + (maxHalfWidth - tailTipHalf) * pow(1 - localT, tailExponent)
    ///   Smoothly contracts from maxHalfWidth to tailTipHalf at the tail.
    /// </para>
    /// <para>
    /// Power curves guarantee a convex outline with no inflection points,
    /// which matches the smooth rail lines of real surfboards.
    /// </para>
    /// </remarks>
    public float GetHalfWidth(float t)
    {
        t = Math.Clamp(t, 0f, 1f);

        if (t <= 0f)
        {
            return 0f;
        }
        else if (t < _widePointT)
        {
            // Nose side: power curve from 0 (nose tip) to maxHalfWidth (wide point)
            float localT = t / _widePointT;
            return _maxHalfWidth * MathF.Pow(localT, _noseExponent);
        }
        else if (t < 1f)
        {
            // Tail side: power curve from maxHalfWidth (wide point) to tailTipHalf (tail)
            float localT = (t - _widePointT) / (1f - _widePointT);
            float blend = MathF.Pow(1f - localT, _tailExponent);
            return _tailTipHalf + (_maxHalfWidth - _tailTipHalf) * blend;
        }
        else
        {
            return _tailTipHalf;
        }
    }
}
