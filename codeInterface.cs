// =============================================================================
// CODE INTERFACE — C# ENTRY POINT
// =============================================================================
//
// Repo-root launcher for the PicoGK surfboard geometry generator.
// This thin wrapper delegates to the SurfboardGeometry library, keeping all
// implementation code inside computationalEngineering/SurfboardGeometry/.
//
// Usage:
//   dotnet run --project codeInterface.csproj -- --type shortboard
//   dotnet run --project codeInterface.csproj -- --type longboard --voxel 0.25
//   dotnet run --project codeInterface.csproj -- --type fish --fins twin
//   dotnet run --project codeInterface.csproj -- --all-fins
//   dotnet run --project codeInterface.csproj -- --config configs/shortboard_default.json
//   dotnet run --project codeInterface.csproj -- --help
//
// =============================================================================

using SurfboardGeometry.Surfboard;
using PicoGK;

namespace CodeInterface;

/// <summary>
/// Repo-root entry point for the C# surfboard geometry generator.
/// Delegates to the SurfboardGeometry library for all geometry operations.
/// </summary>
class Program
{
    static void Main(string[] args)
    {
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
            string? configPath = null;

            for (int i = 0; i < args.Length; i++)
            {
                if (args[i] == "--voxel" && i + 1 < args.Length)
                    voxelSize = float.Parse(args[i + 1]);
                else if (args[i] == "--type" && i + 1 < args.Length)
                    boardType = args[i + 1].ToLower();
                else if (args[i] == "--fins" && i + 1 < args.Length)
                    finConfigStr = args[i + 1].ToLower();
                else if (args[i] == "--all-fins")
                    allFins = true;
                else if (args[i] == "--config" && i + 1 < args.Length)
                    configPath = args[i + 1];
                else if (args[i] == "--help")
                {
                    PrintUsage();
                    return;
                }
            }

            // =================================================================
            // SELECT BOARD PARAMETERS
            // =================================================================
            SurfboardParameters parameters;

            if (configPath != null)
            {
                Console.WriteLine($"Loading custom parameters from: {configPath}");
                parameters = SurfboardParameters.FromJson(configPath);
                boardType = "custom";
            }
            else
            {
                parameters = boardType switch
                {
                    "longboard" => SurfboardParameters.Longboard,
                    "fish" => SurfboardParameters.Fish,
                    _ => SurfboardParameters.Shortboard
                };
            }

            // =================================================================
            // SELECT FIN CONFIGURATION
            // =================================================================
            FinConfiguration finConfiguration = finConfigStr switch
            {
                "thruster" => FinConfiguration.Thruster,
                "twin" => FinConfiguration.Twin,
                "quad" => FinConfiguration.Quad,
                "single" => FinConfiguration.Single,
                _ => parameters.DefaultFinConfiguration
            };

            Console.WriteLine($"Board Type:  {boardType}");
            Console.WriteLine($"Fin Setup:   {(allFins ? "All configurations" : finConfiguration.ToString())}");
            Console.WriteLine($"Voxel Size:  {voxelSize} mm");
            Console.WriteLine();

            // =================================================================
            // PICOGK INITIALIZATION AND GEOMETRY GENERATION
            // =================================================================
            using (new Library(voxelSize))
            {
                if (allFins)
                    SurfboardBody.TaskAllFins(parameters, voxelSize);
                else
                    SurfboardBody.Task(parameters, voxelSize, finConfiguration);
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
        }
    }

    static void PrintUsage()
    {
        Console.WriteLine("Usage: dotnet run --project codeInterface.csproj [-- options]");
        Console.WriteLine();
        Console.WriteLine("Options:");
        Console.WriteLine("  --voxel <size>     Voxel size in mm (default: 0.5)");
        Console.WriteLine("  --type <type>      Board type: shortboard, longboard, or fish");
        Console.WriteLine("  --fins <config>    Fin setup: thruster, twin, quad, or single");
        Console.WriteLine("  --all-fins         Generate all 4 fin configurations");
        Console.WriteLine("  --config <path>    Load custom parameters from JSON file");
        Console.WriteLine("  --help             Show this help message");
        Console.WriteLine();
        Console.WriteLine("Examples:");
        Console.WriteLine("  dotnet run --project codeInterface.csproj");
        Console.WriteLine("  dotnet run --project codeInterface.csproj -- --type longboard");
        Console.WriteLine("  dotnet run --project codeInterface.csproj -- --type fish --voxel 1");
        Console.WriteLine("  dotnet run --project codeInterface.csproj -- --all-fins");
    }
}
