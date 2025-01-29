# Spyffness

A Python library for structural analysis.

## Installation

```bash
pip install Spyffness
```

## Usage

```python
import Spyffness as sf
from Spyffness.Frame import Frame

frame.addMaterial()

frame.addNode(0, [0, 0, 0])
frame.addNode(1, [0, 0, 3])
frame.addBeam(0, 0, 1)
```

## Features

- 3D analysis of frame structures
- Support for different types of materials
- Calculation of stiffness matrices
- Analysis of displacements and reactions