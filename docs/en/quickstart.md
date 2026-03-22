# Quick Start Guide

## Installation

### Install via pip

```bash
pip install nblade
```

### Install from Source

```bash
# Install maturin
pip install maturin

# Clone the repository
git clone https://github.com/nblade/nblade.git
cd nblade

# Build and install
maturin develop --release
```

## Basic Concepts

### What is Geometric Algebra?

Geometric Algebra (also known as Clifford Algebra) is a generalization of traditional vector algebra. It provides:

1. **Unified operation system**: The geometric product unifies the inner and outer products
2. **Geometric object representation**: Scalars, vectors, planes, volumes, etc. are all represented using the same type of object
3. **Coordinate-free**: Results don't depend on specific coordinate systems

### Core Formula

The geometric product is the most fundamental operation:

```
ab = a·b + a∧b
```

Where:
- `ab` is the geometric product
- `a·b` is the inner product (symmetric part)
- `a∧b` is the outer product (antisymmetric part)

## Quick Start

### Creating a Geometric Algebra

```python
import nblade

# Create a 3D Euclidean geometric algebra G(3,0,0)
alg = nblade.Algebra.euclidean(3)

# View algebra information
print(alg)  # Geometric Algebra G(3, 0, 0) (3D)
```

### Creating Vectors

```python
# Create a vector from a list
v = alg.vector([1.0, 2.0, 3.0])

# Get basis vectors
e1, e2, e3 = alg.basis_vectors()

# Create a scalar
s = alg.scalar(5.0)
```

### Basic Operations

```python
# Create two vectors
a = alg.vector([1.0, 0.0, 0.0])
b = alg.vector([1.0, 1.0, 0.0])

# Geometric product
geom = a * b

# Outer product (wedge product)
outer = a ^ b

# Inner product
inner = a | b
```

### Rotation Operations

```python
import math

# Create rotation plane
plane = e1 ^ e2  # xy plane

# Create a rotor (rotate 45 degrees)
rotor = alg.rotor(plane, math.pi / 4)

# Rotate a vector
rotated = e1.rotate_by(rotor)
```

## Next Steps

- [Basics Tutorial](tutorial/01_basics.md) - Dive deeper into geometric algebra fundamentals
- [API Reference](api/algebra.md) - View the complete API documentation
- [Example Code](https://github.com/nblade/nblade/tree/main/examples) - More practical examples