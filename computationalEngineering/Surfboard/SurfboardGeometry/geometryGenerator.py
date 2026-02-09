# -- Geometry Generator -- #

'''
Python wrapper for the C# SurfboardGeometry STL generator.

Invokes the .NET CLI tool as a subprocess to produce voxel-based
surfboard geometry as STL mesh files.

Sean Bowman [02/04/2026]
'''

from __future__ import annotations

import json
import os
import subprocess


class GeometryGenerator:
    '''
    Python wrapper for the C# SurfboardGeometry module.

    Reads board and geometry settings from a JSON config file,
    then invokes `dotnet run --project SurfboardGeometry/` to
    generate an STL mesh file.
    '''

    def __init__(self) -> None:
        self._outputPath: str | None = None
        self._boardType: str | None = None
        self._projectDir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
        )

    @property
    def outputPath(self) -> str | None:
        '''Path to the last generated STL file.'''
        return self._outputPath

    @property
    def boardType(self) -> str | None:
        '''Board type used in the last generation.'''
        return self._boardType

    def generateFromConfig(self, configPath: str) -> str:
        '''
        Generate surfboard STL geometry from a JSON configuration file.

        Reads the "board" and "geometry" sections of the config to
        determine board type, fin configuration, voxel size, and
        output directory. Invokes the C# SurfboardGeometry CLI.

        Parameters:
        -----------
        configPath : str
            Path to the JSON configuration file

        Returns:
        --------
        str : Path to the generated STL file

        Raises:
        -------
        subprocess.CalledProcessError : If dotnet run fails
        FileNotFoundError : If output STL is not found after generation
        '''
        config = self._loadConfig(configPath)

        boardType = config.get('boardType', 'shortboard')
        finConfiguration = config.get('finConfiguration', 'default')
        voxelSize = config.get('voxelSize', 0.5)
        outputDir = config.get('outputDir', 'computationalEngineering/Surfboard/SurfboardGeometry/Output')

        self._boardType = boardType

        # Build the dotnet command
        cmd = [
            'dotnet', 'run',
            '--project', self._projectDir,
            '--',
            '--type', boardType if boardType != 'custom' else 'shortboard',
            '--voxel', str(voxelSize),
        ]

        if finConfiguration != 'default':
            cmd.extend(['--fins', finConfiguration])

        print()
        print('-' * 62)
        print('  GEOMETRY GENERATION (C#)')
        print('-' * 62)
        print(f'  Command: {" ".join(cmd)}')
        print()

        # Run the C# geometry generator
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.stdout:
            print(result.stdout)

        if result.returncode != 0:
            print(f'  Geometry generation failed (exit code {result.returncode})')
            if result.stderr:
                print(f'  stderr: {result.stderr.strip()}')
            raise subprocess.CalledProcessError(result.returncode, cmd)

        # Locate the output STL (go up three levels: SurfboardGeometry -> Surfboard -> computationalEngineering -> repo root)
        rootDir = os.path.dirname(os.path.dirname(os.path.dirname(self._projectDir)))
        stlDir = os.path.join(rootDir, outputDir)
        stlPath = os.path.join(stlDir, 'surfboard.stl')

        if not os.path.exists(stlPath):
            raise FileNotFoundError(f'Expected STL not found at: {stlPath}')

        self._outputPath = stlPath
        print(f'  STL generated: {stlPath}')
        print()

        return stlPath

    def _loadConfig(self, configPath: str) -> dict:
        '''
        Load and flatten relevant config sections from JSON.

        Parameters:
        -----------
        configPath : str
            Path to the JSON configuration file

        Returns:
        --------
        dict : Flattened config with keys: boardType, finConfiguration,
               voxelSize, outputDir
        '''
        with open(configPath, 'r') as f:
            data = json.load(f)

        board = data.get('board', {})
        geometry = data.get('geometry', {})

        return {
            'boardType': board.get('type', 'shortboard'),
            'finConfiguration': board.get('finConfiguration', 'default'),
            'voxelSize': geometry.get('voxelSize', 0.5),
            'outputDir': geometry.get('outputDir', 'computationalEngineering/Surfboard/SurfboardGeometry/Output'),
        }
