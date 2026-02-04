# -- Surfboard Planform Outline -- #

'''
Computes the surfboard planform outline (half-width at each longitudinal station).

Direct port of SurfboardGeometry/Surfboard/Outline.cs.
Uses a two-piece power curve split at the wide point, with exponents
calibrated from NoseWidth and TailWidth measurements.

Sean Bowman [02/03/2026]
'''

from __future__ import annotations

import math

from SurfPhysics.geometry.parameters import SurfboardParameters


class Outline:
    '''
    Surfboard planform outline — half-width at each longitudinal station.

    The board is symmetric about its centerline, so only the half-width is
    computed. Given a normalized position t (0 = nose tip, 1 = tail tip),
    returns the half-width in mm.

    The outline is a two-piece power curve split at the wide point:
    - Nose side: halfWidth = maxHalfWidth * pow(localT, noseExponent)
    - Tail side: halfWidth = tailTipHalf + (maxHalfWidth - tailTipHalf) * pow(1 - localT, tailExponent)

    Exponents are calibrated so the curve passes through the NoseWidth
    and TailWidth measurements at their respective 305mm stations.
    '''

    def __init__(self, params: SurfboardParameters) -> None:
        '''
        Initialize outline from surfboard parameters.

        Parameters:
        -----------
        params : SurfboardParameters
            Board dimensions including outline control points
        '''
        self._params = params

        # Normalized wide point position (0 = nose, 1 = tail)
        self._widePointT: float = params.widePointX / params.length
        self._maxHalfWidth: float = params.halfWidth
        self._tailTipHalf: float = params.tailTipHalfWidth

        # Calibrate nose exponent from NoseWidth at 305mm station
        # The nose power curve: halfWidth = maxHalfWidth * pow(localT, n)
        # We want it to pass through noseWidth/2 at noseStationT
        noseStationT = 305.0 / params.length
        noseRatio = (params.noseWidth / 2.0) / self._maxHalfWidth
        nosePositionRatio = noseStationT / self._widePointT

        if 0.0 < noseRatio < 1.0 and 0.0 < nosePositionRatio < 1.0:
            self._noseExponent = math.log(noseRatio) / math.log(nosePositionRatio)
        else:
            self._noseExponent = 0.5  # fallback: square root curve

        # Clamp to reasonable range
        self._noseExponent = max(0.3, min(1.5, self._noseExponent))

        # Calibrate tail exponent from TailWidth at 305mm from tail tip
        # The tail power curve:
        #   halfWidth = tailTipHalf + (maxHalfWidth - tailTipHalf) * pow(1 - localT, n)
        tailStationFromTail = 305.0 / params.length
        tailStationT = 1.0 - tailStationFromTail  # position from nose
        tailLocalT = (tailStationT - self._widePointT) / (1.0 - self._widePointT)

        tailTargetHalf = params.tailWidth / 2.0
        tailFraction = (tailTargetHalf - self._tailTipHalf) / (self._maxHalfWidth - self._tailTipHalf)
        tailBlendArg = 1.0 - tailLocalT

        if 0.0 < tailFraction < 1.0 and 0.0 < tailBlendArg < 1.0:
            self._tailExponent = math.log(tailFraction) / math.log(tailBlendArg)
        else:
            self._tailExponent = 0.6  # fallback

        # Clamp to reasonable range
        self._tailExponent = max(0.3, min(2.0, self._tailExponent))

    def getHalfWidth(self, t: float) -> float:
        '''
        Get the half-width of the board at a normalized longitudinal position.

        Parameters:
        -----------
        t : float
            Normalized position: 0.0 = nose tip, 1.0 = tail tip

        Returns:
        --------
        float : Half-width in mm at the given position
        '''
        t = max(0.0, min(1.0, t))

        if t <= 0.0:
            return 0.0
        elif t < self._widePointT:
            # Nose side: power curve from 0 to maxHalfWidth
            localT = t / self._widePointT
            return self._maxHalfWidth * math.pow(localT, self._noseExponent)
        elif t < 1.0:
            # Tail side: power curve from maxHalfWidth to tailTipHalf
            localT = (t - self._widePointT) / (1.0 - self._widePointT)
            blend = math.pow(1.0 - localT, self._tailExponent)
            return self._tailTipHalf + (self._maxHalfWidth - self._tailTipHalf) * blend
        else:
            return self._tailTipHalf
