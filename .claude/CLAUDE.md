# Coding Guidelines & Workflow Instructions

**Project:** Personal Computational Engineering Portfolio
**Last Updated:** January 13, 2026

This document contains instructions for Claude Code when working on this project.

---

## File Management

### Permissions

Any time that you are searching for and reading files you do not need to ask for permission.

When specifically asked to search the internet, you do not need to ask for permission to submit the subsequent web request.

When running a file as a test of functionality through bash, you have permission always to do so.

### Temporary File Cleanup

**ALWAYS delete temporary files after editing:**

- Pattern: `*.tmp`, `*.tmp.*`, `*.[0-9]+`, `tmpclaude-*`
- Run cleanup after file edits:
  ```bash
  find . -type f \( -name "*.tmp.*" -o -name "*.tmp" -o -name "tmpclaude-*" \) -delete
  ```
- Check for temp files before committing code
- Clean up job files older than 12 hours in `temp/backgroundJobs/`

### File Organization

- Keep related files together in appropriate directories
- Use descriptive file names
- Place documentation in `documentation/` folders
- Store temporary job data in `temp/` directory (gitignored)
- Check for imports in the top of files and do not re-import libraries into methods that are already imported at the top of the class, and any imports that are needed across at least 2 methods can be imported at the top of a class definition file.

---

## Coding Standards

### Naming Conventions

**Python:**

- **Variables & Functions:** `camelCase`
  ```python
  userInput = {}
  def calculateWaveHeight():
  ```
- **Classes:** `PascalCase`
  ```python
  class WaveModel:
  class SurfboardAnalyzer:
  ```
- **Constants:** `camelCase`
  ```python
  maxCpuCores = 16
  sessionTimeout = 300
  ```
- **Private members:** `_leadingUnderscore` (camelCase after underscore)
  ```python
  def _internalMethod():
  _privateVariable = 10
  ```

**JavaScript (for web components):**

- **Variables & Functions:** `camelCase`
  ```javascript
  let sessionId = '';
  function getSessionId() {
  ```
- **Constants:** `camelCase`
  ```javascript
  const sessionKey = 'app_session_id';
  const timeoutMs = 60000;
  ```

### Code Style

**General Principles:**

- Maximum line length: 120 characters (not strict, use judgment)
- Use type hints for function parameters and returns
- Use verbose and descriptive variable names
- Spaces are prefered between variables and their assignments across an equal sigh (=)
- Important section break comments are of the format:
  ```python
  #-*70#
  # -- Important Title -- #
  #-*70#
  ```
- Include docstrings for all classes and public methods
- **ALWAYS use single quotes (') for all strings** - Never use double quotes (")
  - Python strings: `'example'` not `"example"`
  - Docstrings: `'''docstring'''` not `"""docstring"""`
  - F-strings: `f'value: {x}'` not `f"value: {x}"`
  - Exception: Strings containing single quotes may use double quotes to avoid escaping

**Docstring Format:**

```python
def functionName(param1: str, param2: int) -> dict:
    '''
    Brief description of what the function does.

    Parameters:
    -----------
    param1 : str
        Description of param1
    param2 : int
        Description of param2

    Returns:
    --------
    dict : Description of return value

    Examples:
    ---------
    >>> functionName('test', 42)
    {'result': 'success'}
    '''
```

**Import Organization:**

```python
# Standard library imports
import os
import sys
from typing import Optional, Dict

# Third-party imports
import numpy as np
import plotly.graph_objects as go

# Local imports
from utils.geometry import computeOutline
from utils.waveModel import linearWaveTheory
```

---

## Development Workflow

### Before Starting Work

1. Check git status: `git status`
2. Understand current branch and changes
3. Read relevant documentation files
4. Use `TodoWrite` for multi-step tasks

### During Development

**Use TodoWrite for complex tasks:**

```python
# Create todo list for tasks with 3+ steps
TodoWrite([
    {"content": "Task description", "status": "pending", "activeForm": "Doing task"},
    {"content": "Another task", "status": "pending", "activeForm": "Doing another task"}
])
```

**Mark tasks complete immediately:**

- Don't batch completions
- Exactly ONE task should be "in_progress" at a time
- Update status as work progresses

### Code Review Checklist

Before finishing:

- [ ] Delete all temporary files (*.tmp, *.tmp.*, tmpclaude-*)
- [ ] Follow naming conventions (camelCase for Python)
- [ ] Use single quotes (') for all strings, not double quotes (")
- [ ] Include docstrings for new functions/classes
- [ ] Add type hints to function signatures
- [ ] Update relevant documentation

---

# Documentation Standards

### File Headers

**Python files:**

```python

# -- Identifiable Name for File -- #

"""
Brief description of file purpose.

Detailed explanation if needed.

Sean Bowman [MM/DD/YYYY]
"""
```

### Code Comments

- Comment code liberally, especially when calling external functions
- Explain "why" not "what"
- Update comments when code changes
- Use `# TODO:` for future improvements
- Use `# FIXME:` for known issues

### Markdown Documentation

- Use clear headings (##, ###)
- Include code examples in triple backticks with language
- Add table of contents for long documents
- Keep line length reasonable for readability
- Use lists for step-by-step instructions

---

## Error Handling

### General Pattern

```python
try:
    # Main logic here
    result = performCalculation()
except SpecificException as e:
    # Handle specific errors
    logger.error(f'Calculation failed: {e}')
    raise
except Exception as e:
    # Catch-all for unexpected errors
    logger.error(f'Unexpected error: {e}')
    # Clean up resources
finally:
    # Always runs - use for cleanup
    cleanupResources()
```

---

# Git Practices

### Commits

**Only commit when user explicitly requests:**

- Never proactively commit changes
- Ask user before creating commits
- Use descriptive commit messages
- Follow this format:
  ```
  Brief description of change

  Detailed explanation if needed.
  ```

### Branches

- Understand current branch before making changes
- Don't switch branches without user request
- Check if changes belong on current branch

---

## Performance Considerations

### Optimization Priorities

1. **Correctness** - Code must work correctly
2. **Readability** - Code must be maintainable
3. **Performance** - Optimize only when needed

## Security Considerations

### Data Handling

- Never commit sensitive data (.env files, credentials)
- Validate user inputs before processing
- Sanitize file paths to prevent directory traversal
- Validate serialized data before loading

## Project Structure

```
computationalEngineering/
├── SurfboardGeometry/                     # C# / PicoGK surfboard generator
│   ├── Program.cs                         # CLI entry point
│   ├── Surfboard/
│   │   ├── SurfboardParameters.cs         # Parametric dimensions & presets
│   │   ├── SurfboardBody.cs               # Spatial painting geometry generation
│   │   ├── Outline.cs                     # Planform shape (top view)
│   │   ├── RockerProfile.cs               # Bottom curvature (side view)
│   │   ├── CrossSection.cs                # Deck/bottom profile
│   │   ├── FinConfiguration.cs            # Fin setup enum
│   │   └── FinSystem.cs                   # Fin geometry generation
│   ├── Utils/
│   │   └── Constants.cs                   # Physical constants
│   └── Output/                            # Generated STL files
├── PicoGK/                                # LEAP 71 geometry kernel (submodule)
├── SurfPhysics/                           # Python physics simulations (planned)
├── documentation/                         # Technical writeups
└── references/                            # Source material
```

---

## Tools & Dependencies

### Key Libraries

- **PicoGK** - Voxel-based computational geometry kernel (C#)
- **NumPy** - Numerical computations (Python)
- **Plotly** - Interactive plots
- **Three.js** - 3D visualization (planned)

### Development Tools

- **Git** - Version control
- **VSCode** - IDE (Claude Code integration)
- **.NET 10.0** - C# runtime
- **Python 3.x** - Physics simulations

---

## Quick Reference

### Delete Temp Files

```bash
find . -type f \( -name "*.tmp.*" -o -name "*.tmp" -o -name "tmpclaude-*" \) -delete
```

---

## Notes for Claude Code

- **Be proactive** with TodoWrite for complex tasks
- **Ask questions** when requirements are unclear
- **Clean up** temp files automatically after edits
- **Follow conventions** consistently
- **Document** significant changes
- **Never commit** without explicit user request
- **Update documentation** when making architectural changes

---

**Remember:** This is computational engineering software. Accuracy, reliability, and documentation are critical. When in doubt, ask the user for clarification.
