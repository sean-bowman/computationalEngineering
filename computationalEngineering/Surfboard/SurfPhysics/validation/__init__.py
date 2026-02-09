# -- Validation Subpackage -- #

'''
Mesh comparison and geometry validation tools.

Provides SDF-based comparison between generated surfboard geometry
and reference STL files to assess shape accuracy and guide parameter
tuning for improved geometry generation.

Sean Bowman [02/05/2026]
'''

# Main orchestrator
from computationalEngineering.Surfboard.SurfPhysics.validation.meshComparison import (
    MeshComparisonAnalyzer,
    ComparisonResult,
)

# Distance metrics
from computationalEngineering.Surfboard.SurfPhysics.validation.distanceMetrics import (
    DistanceStatistics,
    RegionalDeviation,
    DistanceAnalyzer,
)

# Mesh alignment
from computationalEngineering.Surfboard.SurfPhysics.validation.meshAlignment import (
    MeshAligner,
    AlignmentResult,
)

# Fin segmentation
from computationalEngineering.Surfboard.SurfPhysics.validation.finSegmentation import (
    FinSegmenter,
    SegmentedMesh,
)

# Visualization
from computationalEngineering.Surfboard.SurfPhysics.validation.deviationVisualizer import DeviationVisualizer

# Report generation
from computationalEngineering.Surfboard.SurfPhysics.validation.comparisonReport import (
    ComparisonReportGenerator,
    ComparisonReportData,
)

# Reverse engineering
from computationalEngineering.Surfboard.SurfPhysics.validation.reverseEngineer import (
    ReverseEngineer,
    ExtractedProfiles,
)

__all__ = [
    # Orchestrator
    'MeshComparisonAnalyzer',
    'ComparisonResult',
    # Distance metrics
    'DistanceStatistics',
    'RegionalDeviation',
    'DistanceAnalyzer',
    # Alignment
    'MeshAligner',
    'AlignmentResult',
    # Fin segmentation
    'FinSegmenter',
    'SegmentedMesh',
    # Visualization
    'DeviationVisualizer',
    # Reports
    'ComparisonReportGenerator',
    'ComparisonReportData',
    # Reverse engineering
    'ReverseEngineer',
    'ExtractedProfiles',
]
