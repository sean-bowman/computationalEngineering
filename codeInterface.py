
# -- Code Interface (Entry-Point) -- #

'''

Lightweight entry point for the Computational Engineering Toolkit.

Demonstrates the RocketNozzle module:
    1. Load JSON config or pass dict
    2. Generate MOC contour + converging section
    3. Generate regen cooling channels + volutes
    4. Run coupled heat transfer model
    5. Export STL + TXT

Usage:
    python codeInterface.py

Sean Bowman [02/09/2026]

'''

import os

from computationalEngineering.RocketNozzle import RocketNozzle

os.system('cls' if os.name == 'nt' else 'clear')

#--------------------------------------------------------------------#
# -- Configuration -- #
#--------------------------------------------------------------------#

configPath = 'computationalEngineering/RocketNozzle/configs/exampleNozzle.json'

#--------------------------------------------------------------------#
# -- Generate Nozzle -- #
#--------------------------------------------------------------------#

nozzle = RocketNozzle()
nozzle.generate(configPath=configPath)

print('\n' + '=' * 60)
print('  ROCKET NOZZLE DESIGN SUMMARY')
print('=' * 60)

print(f'\n  Gas Properties:')
print(f'    gamma:                {nozzle.chamberGamma}')
print(f'    R (gas constant):     {nozzle.chamberRGasConstant:.2f} J/kg-K')
print(f'    T0 (stagnation):      {nozzle.chamberStagnationTemperature:.0f} K')
print(f'    Pc (chamber):         {nozzle.chamberPressure / 1e6:.2f} MPa')

print(f'\n  Nozzle Geometry:')
print(f'    Thrust:               {nozzle.thrust:.0f} N')
print(f'    Exit Mach:            {nozzle.exitMachNumber:.4f}')
print(f'    Expansion Ratio:      {nozzle.expansionRatio:.1f}')
print(f'    Throat Area:          {nozzle.throatArea * 1e4:.2f} cm^2')
print(f'    Length Fraction:      {nozzle.lengthFraction}')
print(f'    Contour Points:       {len(nozzle.xNozzleWall)}')

if nozzle.makeCoolingChannels == 'on':
    print(f'\n  Regenerative Cooling:')
    print(f'    Channels:             {nozzle.nChannel}')
    print(f'    Channel Type:         {nozzle.channelType}')
    print(f'    Hot Wall Thickness:   {nozzle.hotWallThickness * 1e3:.1f} mm')
    print(f'    Shell Thickness:      {nozzle.shellThickness * 1e3:.1f} mm')
    print(f'    Material:             {nozzle.material}')
    print(f'    Coolant:              {nozzle.coolant}')
    print(f'    Coolant Inlet T:      {nozzle.coolantInitialTemperature:.0f} K')
    print(f'    Coolant Inlet P:      {nozzle.coolantInitialPressure / 1e6:.1f} MPa')
    print(f'    Coolant Mass Flow:    {nozzle.coolantMassFlow:.1f} kg/s')

if nozzle.makeInletVolute == 'on' or nozzle.makeReturnVolute == 'on':
    print(f'\n  Volutes:')
    print(f'    Inlet Volute:         {nozzle.makeInletVolute}')
    print(f'    Return Volute:        {nozzle.makeReturnVolute}')

print('\n' + '=' * 60)
print('  Done.')
print('=' * 60)
