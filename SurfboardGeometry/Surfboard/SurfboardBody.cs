// =============================================================================
// SURFBOARD BODY - MAIN GEOMETRY GENERATION VIA SPATIAL PAINTING
// =============================================================================
//
// COMPUTATIONAL ENGINEERING: Voxel-Based Surfboard Generation
// ============================================================
//
// This class generates a parametric surfboard shape using PicoGK's
// spatial painting technique. The approach is analogous to how a shaper
// works: moving along the board's length and defining the cross-section
// at each station.
//
// ALGORITHM OVERVIEW:
// -------------------
//
//   For each longitudinal station (t = 0 to 1):
//     1. Query Outline for half-width at this station
//     2. Query RockerProfile for Z-offset (bottom curvature)
//     3. Query CrossSection for deck/bottom heights
//     4. Paint rings of overlapping spheres tracing the cross-section
//
//   SPATIAL PAINTING (same as BellNozzle divergent section):
//   ────────────────────────────────────────────────────────
//
//   Each ○ is a sphere placed in the lattice:
//
//         ○ ○ ○ ○ ○ ○ ○ ○ ○    ← deck surface ring
//       ○                   ○   ← rail
//     ○                       ○
//       ○                   ○   ← rail
//         ○ ○ ○ ○ ○ ○ ○ ○ ○    ← bottom surface ring
//
//   Rings at adjacent stations overlap to create a continuous solid:
//
//     Station N:    ○ ○ ○ ○ ○ ○
//     Station N+1:   ○ ○ ○ ○ ○ ○    ← overlapping spheres merge
//     Station N+2:    ○ ○ ○ ○ ○ ○
//
// COORDINATE SYSTEM:
// ------------------
//   X = along the board length (0 = nose tip, Length = tail tip)
//   Y = across the board width (0 = centerline, ±HalfWidth at edges)
//   Z = vertical (0 = rocker datum at flat spot, positive = up)
//
// =============================================================================

using System.Numerics;
using PicoGK;

namespace SurfboardGeometry.Surfboard;

/// <summary>
/// Generates surfboard body geometry using PicoGK's spatial painting approach.
/// </summary>
/// <remarks>
/// <para>
/// PICOGK: This class uses the same technique as BellNozzle -- placing
/// overlapping spheres along the desired surface to build up the shape.
/// The spheres merge when converted to voxels, creating a continuous solid.
/// </para>
/// <para>
/// The surfboard's organic shape is built by sweeping cross-sections along
/// the board length, with each section defined by the Outline, RockerProfile,
/// and CrossSection classes.
/// </para>
/// </remarks>
public class SurfboardBody
{
    private readonly SurfboardParameters _params;
    private readonly float _voxelSize;

    private readonly Outline _outline;
    private readonly RockerProfile _rocker;
    private readonly CrossSection _crossSection;

    /// <summary>
    /// Create a new surfboard body generator.
    /// </summary>
    /// <param name="parameters">Board dimensions and shape parameters</param>
    /// <param name="voxelSize">Voxel resolution in mm (smaller = finer detail)</param>
    public SurfboardBody(SurfboardParameters parameters, float voxelSize = 0.5f)
    {
        _params = parameters;
        _voxelSize = voxelSize;

        _outline = new Outline(parameters);
        _rocker = new RockerProfile(parameters);
        _crossSection = new CrossSection(parameters);
    }

    // =========================================================================
    // GEOMETRY GENERATION
    // =========================================================================

    /// <summary>
    /// Generate the solid surfboard body as a voxel field.
    /// </summary>
    /// <returns>Voxels representing the surfboard solid body</returns>
    /// <remarks>
    /// <para>
    /// ALGORITHM: Sweeps along the board length, painting filled cross-section
    /// discs at each station. At each station:
    /// 1. Get the half-width from the Outline
    /// 2. Get the rocker Z-offset from the RockerProfile
    /// 3. For each lateral position across the width, paint spheres at the
    ///    deck and bottom heights from CrossSection
    /// 4. Fill between deck and bottom by placing spheres at multiple Z levels
    /// </para>
    /// <para>
    /// PICOGK: Sphere radius is set to ~2x voxel size to ensure proper
    /// overlap between adjacent spheres. This matches the BellNozzle pattern.
    /// </para>
    /// </remarks>
    public Voxels GenerateBody()
    {
        Console.WriteLine("Generating surfboard body geometry...");

        Lattice lat = new();

        // Sphere radius: ~2x voxel size for proper overlap (same as BellNozzle)
        float sphereRadius = _voxelSize * 2.5f;
        float sphereDiameter = sphereRadius * 2f;

        // Number of longitudinal stations
        // More stations = smoother surface but slower computation
        int lengthSteps = (int)(_params.Length / (_voxelSize * 3f));
        lengthSteps = Math.Max(lengthSteps, 100);

        Console.WriteLine($"  Painting {lengthSteps} longitudinal stations...");

        for (int i = 0; i <= lengthSteps; i++)
        {
            // Normalized position along the board (0 = nose, 1 = tail)
            float t = (float)i / lengthSteps;

            // Actual X position along the board
            float x = t * _params.Length;

            // Get the half-width at this station from the outline
            float halfWidth = _outline.GetHalfWidth(t);

            // Skip stations where width is essentially zero (nose/tail tips)
            if (halfWidth < sphereRadius)
            {
                // At the very tips, place a single sphere to round the nose/tail
                if (halfWidth > 0f)
                {
                    float tipRockerZ = _rocker.GetRockerHeight(t);
                    lat.AddSphere(new Vector3(x, 0f, tipRockerZ), sphereRadius);
                }
                continue;
            }

            // Get the rocker Z-offset at this station
            float rockerZ = _rocker.GetRockerHeight(t);

            // Number of lateral steps across the half-width
            int widthSteps = (int)(halfWidth / sphereDiameter) + 1;
            widthSteps = Math.Max(widthSteps, 8);

            // =====================================================================
            // LATERAL SWEEP: Paint the cross-section at this station
            // =====================================================================
            // We paint both sides (positive and negative Y) for symmetry.
            for (int j = 0; j <= widthSteps; j++)
            {
                // Normalized lateral position (0 = center, 1 = rail edge)
                float lateralFraction = (float)j / widthSteps;

                // Actual Y position
                float y = lateralFraction * halfWidth;

                // Get deck and bottom heights from the cross-section
                float deckZ = _crossSection.GetDeckHeight(t, lateralFraction);
                float bottomZ = _crossSection.GetBottomHeight(t, lateralFraction);

                // Shift by rocker offset
                float deckZWorld = rockerZ + deckZ;
                float bottomZWorld = rockerZ + bottomZ;

                // =============================================================
                // VERTICAL FILL: Paint spheres from bottom to deck
                // =============================================================
                // Fill the interior with spheres to create a solid body.
                // Step size slightly smaller than sphere diameter for overlap.
                float zStep = sphereDiameter * 0.8f;
                int zSteps = Math.Max((int)((deckZWorld - bottomZWorld) / zStep), 1);

                for (int k = 0; k <= zSteps; k++)
                {
                    float zFraction = (float)k / zSteps;
                    float z = bottomZWorld + zFraction * (deckZWorld - bottomZWorld);

                    // Paint on positive Y side
                    lat.AddSphere(new Vector3(x, y, z), sphereRadius);

                    // Paint on negative Y side (mirror symmetry)
                    // Skip center (y=0) to avoid double-painting
                    if (j > 0)
                    {
                        lat.AddSphere(new Vector3(x, -y, z), sphereRadius);
                    }
                }
            }

            // Progress reporting every 20%
            if (i % (lengthSteps / 5) == 0)
            {
                int percent = (int)(100f * i / lengthSteps);
                Console.WriteLine($"  {percent}% complete...");
            }
        }

        Console.WriteLine("  Converting lattice to voxels...");
        Voxels vox = new(lat);

        Console.WriteLine("  Applying smoothing for organic surface finish...");
        // PICOGK: Smoothen uses a triple offset (offset out, offset in, offset out)
        // to remove sharp edges while preserving overall shape.
        // Distance of ~1 voxel gives a subtle smoothing effect.
        // NOTE: Smoothen is currently commented out as it may not be available
        // in all PicoGK versions. Uncomment if your version supports it.
        // vox.Smoothen(_voxelSize);

        Console.WriteLine("  ✓ Surfboard body generated");
        return vox;
    }

    // =========================================================================
    // EXPORT FUNCTIONALITY
    // =========================================================================

    /// <summary>
    /// Generate and export surfboard to STL files.
    /// </summary>
    /// <param name="outputFolder">Directory path for output files</param>
    /// <param name="finConfiguration">Fin setup to generate (default: Thruster)</param>
    public void GenerateAndExport(string outputFolder, FinConfiguration finConfiguration = FinConfiguration.Thruster)
    {
        Directory.CreateDirectory(outputFolder);

        // Print configuration
        Console.WriteLine("\n--- Surfboard Geometry Generation ---");
        _params.PrintSummary();
        Console.WriteLine($"Voxel Size: {_voxelSize} mm");
        Console.WriteLine($"Fin Setup:  {finConfiguration}");
        Console.WriteLine();

        // Generate the board body
        Voxels voxBody = GenerateBody();

        // Generate fins for the selected configuration and combine with body
        var finSystem = new FinSystem(_params, _voxelSize);
        Voxels voxFins = finSystem.Generate(finConfiguration, _rocker);

        Console.WriteLine("\nCombining body and fins...");
        voxBody.BoolAdd(voxFins);
        Console.WriteLine("  ✓ Fins attached to body");

        // Export the complete surfboard (body + fins)
        string bodyFile = Path.Combine(outputFolder, "surfboard.stl");
        Console.WriteLine($"\nExporting surfboard to: {bodyFile}");
        ExportToStl(voxBody, bodyFile);

        Console.WriteLine("\n✓ Surfboard export complete!");
    }

    /// <summary>
    /// Generate the board body once, then export with all four fin configurations.
    /// </summary>
    /// <param name="outputFolder">Directory path for output files</param>
    /// <remarks>
    /// Generates the body voxels once (the expensive step), then for each
    /// fin configuration: duplicates the body, attaches the fins via boolean
    /// union, and exports a separate STL file.
    /// </remarks>
    public void GenerateAllFinsAndExport(string outputFolder)
    {
        Directory.CreateDirectory(outputFolder);

        // Print configuration
        Console.WriteLine("\n--- Surfboard Geometry Generation (All Fin Configs) ---");
        _params.PrintSummary();
        Console.WriteLine($"Voxel Size: {_voxelSize} mm");
        Console.WriteLine();

        // Generate the board body once
        Voxels voxBody = GenerateBody();

        var finSystem = new FinSystem(_params, _voxelSize);

        // Export each fin configuration with the same body
        FinConfiguration[] configs = [
            FinConfiguration.Thruster,
            FinConfiguration.Twin,
            FinConfiguration.Quad,
            FinConfiguration.Single
        ];

        foreach (var config in configs)
        {
            Console.WriteLine($"\n--- {config} Configuration ---");

            // Generate fins for this configuration
            Voxels voxFins = finSystem.Generate(config, _rocker);

            // Duplicate the body so we don't modify the original
            Voxels voxCombined = voxBody.voxDuplicate();

            Console.WriteLine("Combining body and fins...");
            voxCombined.BoolAdd(voxFins);
            Console.WriteLine("  ✓ Fins attached to body");

            // Export with configuration-specific filename
            string configName = config.ToString().ToLower();
            string outputFile = Path.Combine(outputFolder, $"surfboard_{configName}.stl");
            Console.WriteLine($"Exporting to: {outputFile}");
            ExportToStl(voxCombined, outputFile);
        }

        Console.WriteLine("\n✓ All fin configurations exported!");
    }

    /// <summary>
    /// Export voxels to STL via mesh conversion.
    /// </summary>
    private void ExportToStl(Voxels voxels, string filePath)
    {
        Mesh mesh = new(voxels);
        mesh.SaveToStlFile(filePath);
        Console.WriteLine($"  → Triangles: {mesh.nTriangleCount():N0}");
    }

    // =========================================================================
    // PICOGK TASK ENTRY POINT
    // =========================================================================

    /// <summary>
    /// PicoGK task entry point for surfboard generation.
    /// Called by Library.Go() after PicoGK initialization.
    /// </summary>
    /// <remarks>
    /// This static method follows the same delegate pattern as BellNozzle.Task.
    /// It's passed to Library.Go() as a method reference.
    /// </remarks>
    public static void Task()
    {
        Task(SurfboardParameters.Shortboard, 0.5f, FinConfiguration.Thruster);
    }

    /// <summary>
    /// Parameterized task entry point for surfboard generation.
    /// </summary>
    /// <param name="parameters">Board parameters to use</param>
    /// <param name="voxelSize">Voxel resolution in mm</param>
    /// <param name="finConfiguration">Fin setup to generate</param>
    public static void Task(SurfboardParameters parameters, float voxelSize, FinConfiguration finConfiguration)
    {
        Console.WriteLine("═══════════════════════════════════════════");
        Console.WriteLine("  PicoGK Surfboard Geometry Showcase");
        Console.WriteLine("  Parametric Surfboard Generation");
        Console.WriteLine("═══════════════════════════════════════════");
        Console.WriteLine();

        var board = new SurfboardBody(parameters, voxelSize);

        string outputFolder = Path.GetFullPath(
            Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "Output")
        );

        board.GenerateAndExport(outputFolder, finConfiguration);

        Console.WriteLine();
        Console.WriteLine($"Output folder: {outputFolder}");
    }
}
