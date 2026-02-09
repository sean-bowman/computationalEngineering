# -- Validation Module Test -- #

'''
Quick integration test for the mesh comparison module.

Sean Bowman [02/05/2026]
'''

from pathlib import Path
import sys

# Add project root to path
projectRoot = Path(__file__).parent.parent.parent
sys.path.insert(0, str(projectRoot))


def testMeshComparison():
    '''Test the mesh comparison workflow with actual STL files.'''
    from computationalEngineering.Surfboard.SurfPhysics.validation import MeshComparisonAnalyzer

    referencePath = projectRoot / 'referenceThruster.stl'
    generatedPath = projectRoot / 'SurfboardGeometry' / 'Output' / 'surfboard_thruster.stl'

    if not referencePath.exists():
        print(f'Reference STL not found: {referencePath}')
        return False

    if not generatedPath.exists():
        print(f'Generated STL not found: {generatedPath}')
        return False

    print('\n' + '='*60)
    print('MESH COMPARISON TEST')
    print('='*60)

    try:
        # Reference mesh is in different units and axes:
        # - Scale: ~410x (reference units → mm, to match generated board length)
        # - Axes: Reference has length on Z, width on X → Generated has length on X, width on Y
        # Mapping 'zxy': new_X = old_Z (length), new_Y = old_X (width), new_Z = old_Y (thickness)
        analyzer = MeshComparisonAnalyzer(
            referencePath=referencePath,
            generatedPath=generatedPath,
            decimateFactor=0.01,  # 1% for faster test
            forceDecimation=True,  # Also decimate reference mesh
            referenceScaleFactor=410.0,  # Scale reference to mm
            referenceAxisMapping='zxy',  # Remap axes to match generated orientation
        )

        # Run comparison with fewer sample points for speed
        results = analyzer.runFullComparison(
            alignMethod='combined',
            segmentFins=True,
            nSamplePoints=5000,  # Reduced for speed
        )

        # Print summary
        print('\nTest Results:')
        print(f'  Alignment IoU: {results.alignmentResult.boundingBoxIou:.3f}')
        print(f'  Overall RMS: {results.overallStats.rmsMm:.2f} mm')
        print(f'  Hausdorff: {results.hausdorffMm:.2f} mm')

        # Generate report
        reportDir = projectRoot / 'SurfPhysics' / 'validation' / 'test_output'
        paths = analyzer.generateReport(results, reportDir, 'test_comparison')

        print('\n' + '='*60)
        print('TEST PASSED')
        print('='*60)
        return True

    except Exception as e:
        print(f'\nTEST FAILED: {e}')
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = testMeshComparison()
    sys.exit(0 if success else 1)
