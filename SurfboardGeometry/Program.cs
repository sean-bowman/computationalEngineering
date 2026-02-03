// =============================================================================
// PICOGK SURFBOARD SHOWCASE - MAIN ENTRY POINT
// =============================================================================
//
// COMPUTATIONAL ENGINEERING: Parametric Surfboard Generation
// ===========================================================
//
// This project demonstrates computational engineering applied to surfboard
// design. Using LEAP 71's PicoGK geometry kernel, the workflow is:
//
//   1. DEFINE parameters (dimensions, rocker, outline, cross-section)
//   2. GENERATE geometry using voxel-based spatial painting
//   3. EXPORT to STL format for visualization, 3D printing, or analysis
//
// The surfboard shape is defined by three key profiles that are combined
// during geometry generation:
//
//   - OUTLINE: The planform shape (top view) -- defines width at each station
//   - ROCKER: The bottom curvature (side view) -- defines nose/tail kick
//   - CROSS-SECTION: The deck/bottom/rail profile at each station
//
// Three board presets are included:
//   - SHORTBOARD: 6'0" performance board for experienced surfers
//   - LONGBOARD:  9'0" classic noserider for style and glide
//   - FISH:       5'6" retro speed machine for small waves
//
// =============================================================================

using SurfboardGeometry.Surfboard;
using PicoGK;

namespace SurfboardGeometry;

/// <summary>
/// Program entry point for the surfboard geometry generator.
/// </summary>
class Program
{
    /// <summary>
    /// Application entry point.
    /// </summary>
    /// <param name="args">Command-line arguments</param>
    /// <remarks>
    /// PICOGK: The Library constructor initializes PicoGK in headless mode
    /// (no viewer window). This is ideal for CLI batch processing where we
    /// only need geometry generation and STL export, not interactive viewing.
    /// </remarks>
    static void Main(string[] args)
    {
        // =====================================================================
        // BANNER
        // =====================================================================
        Console.WriteLine();
        Console.WriteLine("╔═══════════════════════════════════════════════════════════╗");
        Console.WriteLine("║                                                           ║");
        Console.WriteLine("║         PicoGK Computational Engineering Showcase         ║");
        Console.WriteLine("║         Parametric Surfboard Geometry Generator            ║");
        Console.WriteLine("║                                                           ║");
        Console.WriteLine("║   Using LEAP 71's computational geometry kernel            ║");
        Console.WriteLine("║   for voxel-based physical model generation                ║");
        Console.WriteLine("║                                                           ║");
        Console.WriteLine("╚═══════════════════════════════════════════════════════════╝");
        Console.WriteLine();

        try
        {
            // =================================================================
            // COMMAND-LINE ARGUMENT PARSING
            // =================================================================
            float voxelSize = 0.5f;
            string boardType = "shortboard";
            string finConfigStr = "default";
            bool allFins = false;

            for (int i = 0; i < args.Length; i++)
            {
                if (args[i] == "--voxel" && i + 1 < args.Length)
                {
                    voxelSize = float.Parse(args[i + 1]);
                }
                else if (args[i] == "--type" && i + 1 < args.Length)
                {
                    boardType = args[i + 1].ToLower();
                }
                else if (args[i] == "--fins" && i + 1 < args.Length)
                {
                    finConfigStr = args[i + 1].ToLower();
                }
                else if (args[i] == "--all-fins")
                {
                    allFins = true;
                }
                else if (args[i] == "--help")
                {
                    PrintUsage();
                    return;
                }
            }

            // =================================================================
            // SELECT BOARD PRESET
            // =================================================================
            SurfboardParameters parameters = boardType switch
            {
                "longboard" => SurfboardParameters.Longboard,
                "fish" => SurfboardParameters.Fish,
                _ => SurfboardParameters.Shortboard
            };

            // =================================================================
            // SELECT FIN CONFIGURATION
            // =================================================================
            // If no --fins argument provided, use the board type's default
            FinConfiguration finConfiguration = finConfigStr switch
            {
                "thruster" => FinConfiguration.Thruster,
                "twin"     => FinConfiguration.Twin,
                "quad"     => FinConfiguration.Quad,
                "single"   => FinConfiguration.Single,
                _          => parameters.DefaultFinConfiguration
            };

            Console.WriteLine($"Board Type:  {boardType}");
            Console.WriteLine($"Fin Setup:   {(allFins ? "All configurations" : finConfiguration.ToString())}");
            Console.WriteLine($"Voxel Size:  {voxelSize} mm");
            Console.WriteLine();

            // =================================================================
            // PICOGK INITIALIZATION AND GEOMETRY GENERATION
            // =================================================================
            // Using PicoGK's headless constructor -- no viewer window.
            // The using block initializes the geometry kernel for the
            // duration of the generation, then cleans up automatically.
            using (new Library(voxelSize))
            {
                if (allFins)
                {
                    SurfboardBody.TaskAllFins(parameters, voxelSize);
                }
                else
                {
                    SurfboardBody.Task(parameters, voxelSize, finConfiguration);
                }
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine();
            Console.WriteLine("═══════════════════════════════════════════");
            Console.WriteLine("ERROR: " + ex.Message);
            Console.WriteLine("═══════════════════════════════════════════");
            Console.WriteLine();
            Console.WriteLine("Make sure you have:");
            Console.WriteLine("  1. Cloned this repo with: git clone --recurse-submodules");
            Console.WriteLine("  2. Or initialized submodules: git submodule update --init --recursive");
            Console.WriteLine("  3. Installed the PicoGK runtime (see PicoGK documentation)");
            Console.WriteLine();
            Console.WriteLine("See README.md for detailed setup instructions.");
            Console.WriteLine();
        }
    }

    /// <summary>
    /// Print command-line usage information.
    /// </summary>
    static void PrintUsage()
    {
        Console.WriteLine("Usage: dotnet run [options]");
        Console.WriteLine();
        Console.WriteLine("Options:");
        Console.WriteLine("  --voxel <size>     Voxel size in mm (default: 0.5)");
        Console.WriteLine("  --type <type>      Board type: shortboard, longboard, or fish (default: shortboard)");
        Console.WriteLine("  --fins <config>    Fin setup: thruster, twin, quad, or single (default: per board type)");
        Console.WriteLine("  --all-fins         Generate all 4 fin configurations as separate STL files");
        Console.WriteLine("  --help             Show this help message");
        Console.WriteLine();
        Console.WriteLine("Examples:");
        Console.WriteLine("  dotnet run                                    # Shortboard with thruster fins");
        Console.WriteLine("  dotnet run -- --voxel 0.25                    # Higher resolution");
        Console.WriteLine("  dotnet run -- --type longboard                # 9'0\" longboard with single fin");
        Console.WriteLine("  dotnet run -- --type fish                     # Fish board with twin fins");
        Console.WriteLine("  dotnet run -- --type fish --voxel 1           # Fast preview of fish board");
        Console.WriteLine("  dotnet run -- --fins quad                     # Shortboard with quad fins");
        Console.WriteLine("  dotnet run -- --type longboard --fins single  # Longboard with single fin");
        Console.WriteLine("  dotnet run -- --all-fins --voxel 2            # All fin configs, fast preview");
        Console.WriteLine();
        Console.WriteLine("Board Presets:");
        Console.WriteLine("  shortboard   6'0\" x 19.5\" x 2.44\"  Performance shortboard  (default fins: thruster)");
        Console.WriteLine("  longboard    9'0\" x 22.5\" x 3.0\"   Classic noserider       (default fins: single)");
        Console.WriteLine("  fish         5'6\" x 21.0\" x 2.56\"  Retro speed board       (default fins: twin)");
        Console.WriteLine();
        Console.WriteLine("Fin Configurations:");
        Console.WriteLine("  thruster   3 fins (center + 2 sides)      Most versatile setup");
        Console.WriteLine("  twin       2 side fins                    Fast, loose, skatey feel");
        Console.WriteLine("  quad       4 fins (2 front + 2 rear)      Speed and hold in large surf");
        Console.WriteLine("  single     1 large center fin              Smooth turns, stability");
    }
}
