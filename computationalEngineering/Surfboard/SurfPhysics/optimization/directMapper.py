# -- Direct Parameter Mapper -- #

'''
Maps mesh comparison deviations to geometry parameter adjustments.

Phase 1 of the optimization pipeline: deterministic corrections based
on bounding box differences and regional deviation patterns.

Sean Bowman [02/05/2026]
'''

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from computationalEngineering.SurfPhysics.validation.meshComparison import ComparisonResult
from computationalEngineering.SurfPhysics.validation.distanceMetrics import RegionalDeviation


######################################################################
# -- Parameter Bounds -- #
######################################################################

# Physical bounds for surfboard parameters (all in mm)
PARAMETER_BOUNDS = {
    'Length': (1500.0, 3500.0),           # 5' to 11.5'
    'MaxWidth': (400.0, 650.0),           # 16' to 25.5'
    'MaxThickness': (45.0, 100.0),        # 1.75' to 4'
    'NoseWidth': (150.0, 500.0),          # 6' to 20'
    'TailWidth': (200.0, 500.0),          # 8' to 20'
    'WidePointOffset': (-100.0, 100.0),   # +/- 4'
    'NoseRocker': (40.0, 250.0),          # 1.5' to 10'
    'TailRocker': (15.0, 100.0),          # 0.5' to 4'
    'DeckCrown': (2.0, 20.0),             # 0.1' to 0.8'
    'BottomConcave': (0.0, 8.0),          # 0 to 0.3'
    'RailRadius': (5.0, 25.0),            # 0.2' to 1'
}


######################################################################
# -- Direct Parameter Mapper -- #
######################################################################

class DirectParameterMapper:
    '''
    Computes parameter corrections from mesh comparison results.

    Uses bounding box differences and regional deviations to compute
    deterministic adjustments to geometry parameters.
    '''

    def __init__(self, comparisonResult: ComparisonResult) -> None:
        '''
        Initialize with comparison results.

        Parameters:
        -----------
        comparisonResult : ComparisonResult
            Results from MeshComparisonAnalyzer.runFullComparison()
        '''
        self._result = comparisonResult
        self._refBounds = comparisonResult.referenceMesh.bounds
        self._genBounds = comparisonResult.alignedMesh.bounds

        # Build regional stats lookup
        self._regionalStats = {}
        for rd in comparisonResult.regionalDeviations:
            self._regionalStats[rd.region] = rd

    def computeInitialCorrections(self) -> dict[str, float]:
        '''
        Compute parameter corrections from geometric differences.

        Returns:
        --------
        dict[str, float] : Parameter name to adjustment delta (in mm)
        '''
        corrections = {}

        # Bounding box size differences (reference - generated)
        refSize = self._refBounds[1] - self._refBounds[0]
        genSize = self._genBounds[1] - self._genBounds[0]

        lengthDiff = refSize[0] - genSize[0]    # X-axis
        widthDiff = refSize[1] - genSize[1]     # Y-axis
        thicknessDiff = refSize[2] - genSize[2] # Z-axis

        # ========================================================
        # Primary Dimension Corrections (from bounding box)
        # ========================================================

        # Length: Direct mapping
        corrections['Length'] = lengthDiff

        # MaxWidth: Account for symmetry (half-width * 2)
        corrections['MaxWidth'] = widthDiff

        # MaxThickness: Z-difference with scaling factor
        # Thickness affects deck crown contribution, so apply 0.9 factor
        corrections['MaxThickness'] = thicknessDiff * 0.9

        # ========================================================
        # Regional Corrections (from deviation statistics)
        # ========================================================

        # Nose region corrections
        if 'nose' in self._regionalStats:
            noseStats = self._regionalStats['nose']
            # Negative mean = generated is smaller → increase parameters
            # Positive mean = generated is larger → decrease parameters
            noseBias = -noseStats.stats.meanMm

            # NoseWidth affects lateral extent in nose region
            corrections['NoseWidth'] = noseBias * 1.5

            # NoseRocker affects vertical position at nose tip
            if abs(noseStats.stats.meanMm) > 3.0:
                corrections['NoseRocker'] = noseBias * 0.5

        # Tail region corrections
        if 'tail' in self._regionalStats:
            tailStats = self._regionalStats['tail']
            tailBias = -tailStats.stats.meanMm

            corrections['TailWidth'] = tailBias * 1.5
            if abs(tailStats.stats.meanMm) > 3.0:
                corrections['TailRocker'] = tailBias * 0.3

        # Middle region corrections
        if 'middle' in self._regionalStats:
            middleStats = self._regionalStats['middle']
            middleBias = -middleStats.stats.meanMm

            # Adjust concave/crown if middle region has vertical offset
            if abs(middleBias) > 2.0:
                corrections['BottomConcave'] = middleBias * 0.3
                corrections['DeckCrown'] = middleBias * 0.5

        # WidePointOffset: Based on asymmetry between forward/rear
        if 'forward' in self._regionalStats and 'rear' in self._regionalStats:
            forwardRms = self._regionalStats['forward'].stats.rmsMm
            rearRms = self._regionalStats['rear'].stats.rmsMm

            # If rear has higher deviation, wide point may be too far forward
            if rearRms > forwardRms * 1.3:
                corrections['WidePointOffset'] = 10.0  # Shift toward tail
            elif forwardRms > rearRms * 1.3:
                corrections['WidePointOffset'] = -10.0  # Shift toward nose

        return corrections

    def applyCorrections(
        self,
        currentParams: dict,
        corrections: dict[str, float],
        dampingFactor: float = 0.7,
    ) -> dict:
        '''
        Apply corrections to current parameters with damping.

        Parameters:
        -----------
        currentParams : dict
            Current parameter values
        corrections : dict[str, float]
            Computed corrections
        dampingFactor : float
            Scale factor (0-1) to avoid overshooting

        Returns:
        --------
        dict : Updated parameter values
        '''
        newParams = currentParams.copy()

        for param, delta in corrections.items():
            if param in newParams:
                # Apply damped correction
                newParams[param] = float(newParams[param]) + delta * dampingFactor

        # Enforce parameter bounds
        newParams = self._enforceBounds(newParams)

        return newParams

    def _enforceBounds(self, params: dict) -> dict:
        '''Clamp parameters to physically realistic ranges.'''
        for param, (minVal, maxVal) in PARAMETER_BOUNDS.items():
            if param in params:
                params[param] = max(minVal, min(maxVal, float(params[param])))
        return params

    def getSummary(self) -> str:
        '''Get human-readable summary of the analysis.'''
        refSize = self._refBounds[1] - self._refBounds[0]
        genSize = self._genBounds[1] - self._genBounds[0]

        lines = [
            'Bounding Box Comparison:',
            f'  Reference: {refSize[0]:.1f} x {refSize[1]:.1f} x {refSize[2]:.1f} mm',
            f'  Generated: {genSize[0]:.1f} x {genSize[1]:.1f} x {genSize[2]:.1f} mm',
            f'  Difference: {refSize[0]-genSize[0]:+.1f} x {refSize[1]-genSize[1]:+.1f} x {refSize[2]-genSize[2]:+.1f} mm',
            '',
            'Regional Analysis:',
        ]

        for region, rd in self._regionalStats.items():
            lines.append(f'  {region}: RMS={rd.stats.rmsMm:.1f}mm, Mean={rd.stats.meanMm:.1f}mm')

        return '\n'.join(lines)
