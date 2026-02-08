// =============================================================================
// FIN CONFIGURATION - SURFBOARD FIN SETUP OPTIONS
// =============================================================================
//
// Different fin configurations change how a surfboard turns, holds, and
// generates speed. The configuration is independent of board type, though
// certain pairings are traditional (e.g., Fish + Twin, Longboard + Single).
//
// =============================================================================

namespace SurfboardGeometry.Surfboard;

/// <summary>
/// Fin configuration options for surfboard setup.
/// </summary>
/// <remarks>
/// <para>
/// SURFBOARD: Fins provide directional stability and control. Different
/// configurations trade off between speed, hold, looseness, and control.
/// The configuration is independent of board type -- any combination can
/// be used, though traditional pairings exist.
/// </para>
/// </remarks>
public enum FinConfiguration
{
    /// <summary>3 fins: center + 2 sides. Most versatile setup.</summary>
    Thruster,

    /// <summary>2 side fins only. Loose, fast, good for small waves.</summary>
    Twin,

    /// <summary>4 fins: 2 front + 2 rear trailers. Speed and hold.</summary>
    Quad,

    /// <summary>1 large center fin. Traditional longboard setup.</summary>
    Single
}
