// =============================================================================
// PHYSICAL AND GEOMETRIC CONSTANTS FOR SURFBOARD DESIGN
// =============================================================================
//
// SURFBOARD DESIGN: Key Physical Properties
// ===========================================
//
// Understanding the physics of surfing requires knowledge of several
// physical constants and material properties:
//
// WATER PROPERTIES:
//   - Seawater density: 1025 kg/m^3 (higher than fresh water due to salt)
//   - Kinematic viscosity: ~1.0e-6 m^2/s (affects drag calculations)
//   - Surface tension: ~0.072 N/m (affects small-scale wave behavior)
//
// SURFBOARD MATERIALS:
//   - PU foam core: ~30-40 kg/m^3 (polyurethane, traditional)
//   - EPS foam core: ~15-25 kg/m^3 (expanded polystyrene, lighter)
//   - Fiberglass laminate: ~1800 kg/m^3 (structural shell)
//   - Typical finished board density: ~100-200 kg/m^3 (depending on glass schedule)
//
// VOXEL RESOLUTION:
//   - 0.25mm: High detail (slow, large files, best for final output)
//   - 0.50mm: Standard (good balance of detail and speed)
//   - 1.00mm: Fast preview (low detail, quick iteration)
//   - 2.00mm: Very fast preview (for parameter tuning)
//
// =============================================================================

namespace SurfboardGeometry.Utils;

/// <summary>
/// Physical and geometric constants for surfboard design and analysis.
/// All SI units unless otherwise noted.
/// </summary>
public static class Constants
{
    // =========================================================================
    // VOXEL RESOLUTION
    // =========================================================================

    /// <summary>Default voxel size in mm for standard-quality output</summary>
    public const float DefaultVoxelSizeMm = 0.5f;

    /// <summary>High-resolution voxel size in mm for detailed output</summary>
    public const float HighResVoxelSizeMm = 0.25f;

    /// <summary>Preview voxel size in mm for fast iteration</summary>
    public const float PreviewVoxelSizeMm = 1.0f;

    // =========================================================================
    // WATER PROPERTIES
    // =========================================================================

    /// <summary>
    /// Seawater density [kg/m^3].
    /// Used for buoyancy and hydrodynamic force calculations.
    /// </summary>
    public const float SeawaterDensity = 1025f;

    /// <summary>
    /// Gravitational acceleration [m/s^2].
    /// Standard gravity used in wave speed and force calculations.
    /// </summary>
    public const float Gravity = 9.81f;

    /// <summary>
    /// Kinematic viscosity of seawater at 20°C [m^2/s].
    /// Used in Reynolds number calculations for drag estimation.
    /// </summary>
    public const float SeawaterKinematicViscosity = 1.05e-6f;

    // =========================================================================
    // SURFBOARD MATERIAL PROPERTIES
    // =========================================================================

    /// <summary>
    /// Typical PU (polyurethane) foam core density [kg/m^3].
    /// Traditional surfboard blank material.
    /// </summary>
    public const float PuFoamDensity = 35f;

    /// <summary>
    /// Typical EPS (expanded polystyrene) foam core density [kg/m^3].
    /// Lighter alternative to PU, used in epoxy boards.
    /// </summary>
    public const float EpsFoamDensity = 20f;

    /// <summary>
    /// Fiberglass laminate density [kg/m^3].
    /// The resin-saturated glass cloth shell.
    /// </summary>
    public const float FiberglassDensity = 1800f;

    /// <summary>
    /// Approximate fiberglass shell thickness [mm].
    /// Typical 4+4oz bottom, 4oz deck glass schedule.
    /// </summary>
    public const float GlassShellThickness = 1.5f;

    // =========================================================================
    // UNIT CONVERSIONS
    // =========================================================================

    /// <summary>Conversion factor: inches to millimeters</summary>
    public const float InchesToMm = 25.4f;

    /// <summary>Conversion factor: feet to millimeters</summary>
    public const float FeetToMm = 304.8f;

    /// <summary>Conversion factor: cubic millimeters to liters</summary>
    public const float CubicMmToLiters = 1e-6f;

    // =========================================================================
    // OUTPUT PATHS
    // =========================================================================

    /// <summary>
    /// Output folder path relative to the build output directory.
    /// </summary>
    public static class Paths
    {
        public static string OutputFolder => Path.Combine(
            AppDomain.CurrentDomain.BaseDirectory,
            "..", "..", "..", "Output"
        );
    }
}
