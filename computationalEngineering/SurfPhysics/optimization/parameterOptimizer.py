# -- Parameter Optimizer -- #

'''
Main optimization orchestrator for matching generated geometry to reference.

Two-phase optimization:
1. Direct parameter mapping for initial large corrections
2. Nelder-Mead simplex optimization for fine-tuning

Sean Bowman [02/05/2026]
'''

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

try:
    import trimesh
except ImportError:
    trimesh = None

from computationalEngineering.SurfPhysics.optimization.directMapper import DirectParameterMapper, PARAMETER_BOUNDS
from computationalEngineering.SurfPhysics.optimization.coordinateTransformer import CoordinateTransformer


######################################################################
# -- Default Parameters -- #
######################################################################

DEFAULT_PARAMS = {
    'shortboard': {
        'Length': 1828.0,
        'MaxWidth': 495.0,
        'MaxThickness': 62.0,
        'NoseWidth': 300.0,
        'TailWidth': 380.0,
        'WidePointOffset': -25.0,
        'NoseRocker': 120.0,
        'TailRocker': 40.0,
        'DeckCrown': 8.0,
        'BottomConcave': 2.0,
        'RailRadius': 12.0,
    },
    'longboard': {
        'Length': 2743.0,
        'MaxWidth': 570.0,
        'MaxThickness': 75.0,
        'NoseWidth': 430.0,
        'TailWidth': 370.0,
        'WidePointOffset': 0.0,
        'NoseRocker': 180.0,
        'TailRocker': 25.0,
        'DeckCrown': 10.0,
        'BottomConcave': 1.0,
        'RailRadius': 18.0,
    },
    'fish': {
        'Length': 1676.0,
        'MaxWidth': 533.0,
        'MaxThickness': 65.0,
        'NoseWidth': 400.0,
        'TailWidth': 430.0,
        'WidePointOffset': -25.0,
        'NoseRocker': 80.0,
        'TailRocker': 30.0,
        'DeckCrown': 6.0,
        'BottomConcave': 3.0,
        'RailRadius': 10.0,
    },
}


######################################################################
# -- Optimization Result -- #
######################################################################

@dataclass
class OptimizationResult:
    '''Result of optimization run.'''
    success: bool
    finalParams: dict
    finalRmsMm: float
    iterations: int
    history: list[dict]
    message: str


######################################################################
# -- Parameter Optimizer -- #
######################################################################

class ParameterOptimizer:
    '''
    Optimizes surfboard parameters to match reference geometry.

    Two-phase approach:
    1. Direct mapping: Use bounding box and regional deviations for fast initial correction
    2. Nelder-Mead: Fine-tune using gradient-free optimization

    Usage:
        optimizer = ParameterOptimizer(
            referencePath='reference.stl',
            generatorDir='SurfboardGeometry/',
            outputDir='optimization_output/',
        )
        result = optimizer.runOptimization(
            initialBoardType='shortboard',
            targetRmsMm=5.0,
        )
    '''

    # Parameters to optimize (excludes TailShape, fin config, etc.)
    OPTIMIZATION_PARAMS = [
        'Length', 'MaxWidth', 'MaxThickness',
        'NoseWidth', 'TailWidth', 'WidePointOffset',
        'NoseRocker', 'TailRocker',
    ]

    def __init__(
        self,
        referencePath: str | Path,
        generatorDir: str | Path,
        outputDir: str | Path,
        referenceScaleFactor: float | None = None,
        referenceAxisMapping: str | None = None,
    ) -> None:
        '''
        Initialize the optimizer.

        Parameters:
        -----------
        referencePath : str | Path
            Path to reference STL file
        generatorDir : str | Path
            Directory containing SurfboardGeometry C# project
        outputDir : str | Path
            Directory for output files and iteration data
        referenceScaleFactor : float | None
            Scale factor for reference mesh (auto-detect if None)
        referenceAxisMapping : str | None
            Axis mapping for reference mesh (auto-detect if None)
        '''
        if trimesh is None:
            raise ImportError('trimesh is required for optimization')

        self.referencePath = Path(referencePath)
        self.generatorDir = Path(generatorDir)
        self.outputDir = Path(outputDir)
        self.outputDir.mkdir(parents=True, exist_ok=True)

        # Load reference mesh and detect transformation
        print(f'Loading reference: {self.referencePath.name}')
        self._refMesh = trimesh.load(str(self.referencePath))

        if referenceScaleFactor is None or referenceAxisMapping is None:
            detectedScale, detectedMapping = CoordinateTransformer.autoDetect(
                self._refMesh,
                expectedLengthMm=1828.0,  # Assume shortboard by default
            )
            self._scaleFactor = referenceScaleFactor or detectedScale
            self._axisMapping = referenceAxisMapping or detectedMapping
            print(f'  Auto-detected transformation: scale={self._scaleFactor:.1f}, axes={self._axisMapping}')
        else:
            self._scaleFactor = referenceScaleFactor
            self._axisMapping = referenceAxisMapping

        self._history = []

    def runOptimization(
        self,
        initialBoardType: str = 'shortboard',
        targetRmsMm: float = 5.0,
        maxPhase1Iterations: int = 3,
        maxPhase2Iterations: int = 50,
        finConfiguration: str = 'thruster',
    ) -> OptimizationResult:
        '''
        Run complete two-phase optimization.

        Parameters:
        -----------
        initialBoardType : str
            Starting board type ('shortboard', 'longboard', 'fish')
        targetRmsMm : float
            Target RMS deviation to achieve
        maxPhase1Iterations : int
            Maximum iterations for direct mapping phase
        maxPhase2Iterations : int
            Maximum iterations for Nelder-Mead phase
        finConfiguration : str
            Fin configuration to use

        Returns:
        --------
        OptimizationResult : Optimization outcome
        '''
        print('\n' + '='*60)
        print('PARAMETER OPTIMIZATION')
        print('='*60)
        print(f'Reference: {self.referencePath.name}')
        print(f'Initial type: {initialBoardType}')
        print(f'Target RMS: {targetRmsMm:.1f} mm')
        print('='*60 + '\n')

        # Get initial parameters
        currentParams = DEFAULT_PARAMS.get(initialBoardType, DEFAULT_PARAMS['shortboard']).copy()
        self._finConfiguration = finConfiguration

        # =====================================================================
        # PHASE 1: Direct Parameter Mapping
        # =====================================================================
        print('PHASE 1: Direct Parameter Mapping')
        print('-'*40)

        for iteration in range(maxPhase1Iterations):
            print(f'\n--- Phase 1 Iteration {iteration + 1} ---')

            # Generate board with current parameters
            generatedPath = self._generateBoard(currentParams, f'phase1_iter{iteration + 1}')
            if generatedPath is None:
                return OptimizationResult(
                    success=False,
                    finalParams=currentParams,
                    finalRmsMm=float('inf'),
                    iterations=iteration + 1,
                    history=self._history,
                    message='Generation failed',
                )

            # Compare to reference
            results = self._compareToReference(generatedPath)
            if results is None:
                continue

            rmsMm = results.overallStats.rmsMm
            print(f'  Current RMS: {rmsMm:.2f} mm')

            # Record history
            self._history.append({
                'phase': 1,
                'iteration': iteration + 1,
                'params': currentParams.copy(),
                'rmsMm': rmsMm,
            })

            # Check convergence
            if rmsMm < targetRmsMm:
                print(f'\n✓ Target achieved in Phase 1!')
                return OptimizationResult(
                    success=True,
                    finalParams=currentParams,
                    finalRmsMm=rmsMm,
                    iterations=iteration + 1,
                    history=self._history,
                    message=f'Converged in Phase 1, iteration {iteration + 1}',
                )

            # Compute and apply corrections
            mapper = DirectParameterMapper(results)
            corrections = mapper.computeInitialCorrections()

            print(f'  Key corrections:')
            for param, delta in sorted(corrections.items(), key=lambda x: abs(x[1]), reverse=True)[:5]:
                if abs(delta) > 0.5:
                    print(f'    {param}: {delta:+.1f} mm')

            # Apply with decreasing damping
            dampingFactor = 0.7 * (0.8 ** iteration)
            currentParams = mapper.applyCorrections(currentParams, corrections, dampingFactor)

        # =====================================================================
        # Check if Phase 1 was sufficient
        # =====================================================================
        lastRms = self._history[-1]['rmsMm'] if self._history else float('inf')
        if lastRms < targetRmsMm:
            return OptimizationResult(
                success=True,
                finalParams=currentParams,
                finalRmsMm=lastRms,
                iterations=len(self._history),
                history=self._history,
                message='Converged in Phase 1',
            )

        # =====================================================================
        # PHASE 2: Nelder-Mead Fine-Tuning
        # =====================================================================
        print('\n' + 'PHASE 2: Nelder-Mead Optimization')
        print('-'*40)

        try:
            from scipy.optimize import minimize
        except ImportError:
            print('scipy not available, skipping Phase 2')
            return OptimizationResult(
                success=lastRms < targetRmsMm * 2,
                finalParams=currentParams,
                finalRmsMm=lastRms,
                iterations=len(self._history),
                history=self._history,
                message='Phase 1 complete, scipy not available for Phase 2',
            )

        # Build optimization vector
        paramNames = [p for p in self.OPTIMIZATION_PARAMS if p in currentParams]
        x0 = np.array([currentParams[p] for p in paramNames])

        # Parameter scales for better conditioning
        scales = np.array([
            PARAMETER_BOUNDS[p][1] - PARAMETER_BOUNDS[p][0]
            for p in paramNames
        ]) / 10.0  # Scale to ~10 units

        x0Scaled = x0 / scales
        self._nmIterCount = 0

        def objective(xScaled):
            xActual = xScaled * scales
            params = currentParams.copy()
            for name, val in zip(paramNames, xActual):
                params[name] = float(val)

            self._nmIterCount += 1

            # Generate and compare
            genPath = self._generateBoard(params, f'phase2_iter{self._nmIterCount}')
            if genPath is None:
                return 1000.0  # Penalty for failed generation

            results = self._compareToReference(genPath)
            if results is None:
                return 1000.0

            rmsMm = results.overallStats.rmsMm

            self._history.append({
                'phase': 2,
                'iteration': self._nmIterCount,
                'params': params.copy(),
                'rmsMm': rmsMm,
            })

            if self._nmIterCount % 5 == 0:
                print(f'  Phase 2 iteration {self._nmIterCount}: RMS = {rmsMm:.2f} mm')

            # Early stop check
            if rmsMm < targetRmsMm:
                raise StopIteration(f'Target achieved: {rmsMm:.2f} mm')

            return rmsMm

        try:
            result = minimize(
                objective,
                x0Scaled,
                method='Nelder-Mead',
                options={
                    'maxiter': maxPhase2Iterations,
                    'xatol': 0.01,
                    'fatol': 0.5,
                },
            )
            xOptimal = result.x * scales
            finalRms = result.fun

        except StopIteration as e:
            print(f'  Early stop: {e}')
            # Use last evaluated point
            bestHist = min(self._history, key=lambda h: h['rmsMm'])
            xOptimal = np.array([bestHist['params'][p] for p in paramNames])
            finalRms = bestHist['rmsMm']

        # Build final parameters
        finalParams = currentParams.copy()
        for name, val in zip(paramNames, xOptimal):
            finalParams[name] = float(val)

        print('\n' + '='*60)
        print('OPTIMIZATION COMPLETE')
        print('='*60)
        print(f'Final RMS: {finalRms:.2f} mm')
        print(f'Total iterations: {len(self._history)}')
        print('='*60)

        return OptimizationResult(
            success=finalRms < targetRmsMm,
            finalParams=finalParams,
            finalRmsMm=finalRms,
            iterations=len(self._history),
            history=self._history,
            message=f'Optimization complete, final RMS = {finalRms:.2f} mm',
        )

    def _generateBoard(self, params: dict, iterName: str) -> Path | None:
        '''Generate surfboard with given parameters.'''
        # Write parameters to JSON config
        configPath = self.outputDir / f'{iterName}_config.json'
        config = {'customParameters': params}
        configPath.write_text(json.dumps(config, indent=2))

        # Run C# generator
        try:
            result = subprocess.run(
                ['dotnet', 'run', '--', '--config', str(configPath), '--fins', self._finConfiguration],
                cwd=str(self.generatorDir),
                capture_output=True,
                text=True,
                timeout=180,
            )

            if result.returncode != 0:
                print(f'  Generator failed: {result.stderr[:200]}')
                return None

            # Find generated STL
            stlPath = self.generatorDir / 'Output' / 'surfboard.stl'
            if not stlPath.exists():
                # Try with fin config suffix
                stlPath = self.generatorDir / 'Output' / f'surfboard_{self._finConfiguration}.stl'

            if stlPath.exists():
                return stlPath
            else:
                print(f'  STL not found after generation')
                return None

        except subprocess.TimeoutExpired:
            print('  Generator timeout')
            return None
        except Exception as e:
            print(f'  Generator error: {e}')
            return None

    def _compareToReference(self, generatedPath: Path):
        '''Compare generated mesh to reference.'''
        try:
            from computationalEngineering.SurfPhysics.validation import MeshComparisonAnalyzer

            analyzer = MeshComparisonAnalyzer(
                referencePath=self.referencePath,
                generatedPath=generatedPath,
                decimateFactor=0.02,  # Fast comparison
                forceDecimation=True,
                referenceScaleFactor=self._scaleFactor,
                referenceAxisMapping=self._axisMapping,
            )

            return analyzer.runFullComparison(
                alignMethod='combined',
                segmentFins=False,
                nSamplePoints=10000,
            )

        except Exception as e:
            print(f'  Comparison failed: {e}')
            return None

    def saveOptimizedConfig(
        self,
        result: OptimizationResult,
        filename: str = 'optimized_params.json',
    ) -> Path:
        '''Save optimized parameters to JSON file.'''
        outputPath = self.outputDir / filename
        config = {
            'customParameters': result.finalParams,
            'optimization': {
                'success': result.success,
                'finalRmsMm': result.finalRmsMm,
                'iterations': result.iterations,
                'message': result.message,
            },
        }
        outputPath.write_text(json.dumps(config, indent=2))
        print(f'Saved optimized config to: {outputPath}')
        return outputPath
