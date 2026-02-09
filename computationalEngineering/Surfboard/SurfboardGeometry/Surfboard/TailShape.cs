// =============================================================================
// TAIL SHAPE - SURFBOARD TAIL CONFIGURATION
// =============================================================================

namespace SurfboardGeometry.Surfboard;

/// <summary>
/// Defines the tail shape of the surfboard.
/// </summary>
/// <remarks>
/// The tail shape affects water release, turning characteristics, and
/// overall board feel. Different shapes are suited to different wave
/// conditions and surfing styles.
/// </remarks>
public enum TailShape
{
    /// <summary>
    /// Standard tail that tapers to a narrow point (pin, round, or squash).
    /// Used by shortboards, longboards, and most performance shapes.
    /// </summary>
    Standard,

    /// <summary>
    /// Swallow (fish) tail with a V-shaped center notch creating two lobes.
    /// Classic fish board tail that increases rail length and allows water
    /// to release cleanly from each lobe independently.
    /// </summary>
    Swallow
}
