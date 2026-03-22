# nblade - N-dimensional Blade

**High-Performance Geometric Algebra Library**

[![PyPI](https://img.shields.io/pypi/v/nblade.svg)](https://pypi.org/project/nblade/)
[![Crates.io](https://img.shields.io/crates/v/nblade.svg)](https://crates.io/crates/nblade/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview

nblade is a high-performance geometric algebra library implemented in Rust, with Python bindings.

### Features

- **Arbitrary Dimensions** — Up to 64-dimensional vector spaces
- **Arbitrary Signatures** — Support for G(p, q, r) metric signatures
- **High Performance** — Rust backend with parallel computing and SIMD optimization
- **Dual Representation** — Automatic dense/sparse selection
- **Complete Operations** — All standard geometric algebra operations
- **NumPy Integration** — Create vectors directly from NumPy arrays

## Quick Install

```bash
pip install nblade
```

## Quick Start

```python
import nblade

# Create 3D Euclidean geometric algebra
alg = nblade.Algebra.euclidean(3)

# Get basis vectors
e1, e2, e3 = alg.basis_vectors()

# Create a vector
v = alg.vector([1.0, 2.0, 3.0])

# Geometric product
product = e1 * e2  # Bivector e12

# Outer product (wedge)
wedge = e1 ^ e2  # Bivector e12

# Inner product
inner = e1 | e2  # 0 (orthogonal)

# Rotation
import math
rotor = alg.rotor(e1 ^ e2, math.pi / 2)
rotated = e1.rotate_by(rotor)
```

## Documentation

- [Quick Start Guide](en/quickstart.md)
- [Tutorials](en/tutorial/01_basics.md)
- [API Reference](en/api/algebra.md)
- [Examples](en/examples/tutorials.md)

## Links

| Resource | URL |
|----------|-----|
| GitHub | https://github.com/UynajGI/nblade |
| PyPI | https://pypi.org/project/nblade/ |
| Crates.io | https://crates.io/crates/nblade |

## Supported Operations

| Operation | Symbol | Description |
|-----------|--------|-------------|
| Geometric Product | `a * b` | Fundamental GA operation |
| Outer Product | `a ^ b` | Wedge product, spans subspace |
| Inner Product | `a \| b` | Projection relationship |
| Dual | `a.dual()` | Hodge dual |
| Rotation | `a.rotate_by(R)` | Rotate by rotor |
| Projection | `a.project_to(b)` | Project onto subspace |
| Reflection | `a.reflect_in(n)` | Reflect in plane |

## License

MIT License