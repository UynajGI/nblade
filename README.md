# nblade - N-dimensional Blade

[![Crates.io](https://img.shields.io/crates/v/nblade.svg)](https://crates.io/crates/nblade)
[![PyPI](https://img.shields.io/pypi/v/nblade.svg)](https://pypi.org/project/nblade/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![zread](https://img.shields.io/badge/Ask_Zread-_.svg?style=flat&color=00b0aa&labelColor=000000&logo=data%3Aimage%2Fsvg%2Bxml%3Bbase64%2CPHN2ZyB3aWR0aD0iMTYiIGhlaWdodD0iMTYiIHZpZXdCb3g9IjAgMCAxNiAxNiIgZmlsbD0ibm9uZSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj4KPHBhdGggZD0iTTQuOTYxNTYgMS42MDAxSDIuMjQxNTZDMS44ODgxIDEuNjAwMSAxLjYwMTU2IDEuODg2NjQgMS42MDE1NiAyLjI0MDFWNC45NjAxQzEuNjAxNTYgNS4zMTM1NiAxLjg4ODEgNS42MDAxIDIuMjQxNTYgNS42MDAxSDQuOTYxNTZDNS4zMTUwMiA1LjYwMDEgNS42MDE1NiA1LjMxMzU2IDUuNjAxNTYgNC45NjAxVjIuMjQwMUM1LjYwMTU2IDEuODg2NjQgNS4zMTUwMiAxLjYwMDEgNC45NjE1NiAxLjYwMDFaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik00Ljk2MTU2IDEwLjM5OTlIMi4yNDE1NkMxLjg4ODEgMTAuMzk5OSAxLjYwMTU2IDEwLjY4NjQgMS42MDE1NiAxMS4wMzk5VjEzLjc1OTlDMS42MDE1NiAxNC4xMTM0IDEuODg4MSAxNC4zOTk5IDIuMjQxNTYgMTQuMzk5OUg0Ljk2MTU2QzUuMzE1MDIgMTQuMzk5OSA1LjYwMTU2IDE0LjExMzQgNS42MDE1NiAxMy43NTk5VjExLjAzOTlDNS42MDE1NiAxMC42ODY0IDUuMzE1MDIgMTAuMzk5OSA0Ljk2MTU2IDEwLjM5OTlaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik0xMy43NTg0IDEuNjAwMUgxMS4wMzg0QzEwLjY4NSAxLjYwMDEgMTAuMzk4NCAxLjg4NjY0IDEwLjM5ODQgMi4yNDAxVjQuOTYwMUMxMC4zOTg0IDUuMzEzNTYgMTAuNjg1IDUuNjAwMSAxMS4wMzg0IDUuNjAwMUgxMy43NTg0QzE0LjExMTkgNS42MDAxIDE0LjM5ODQgNS4zMTM1NiAxNC4zOTg0IDQuOTYwMVYyLjI0MDFDMTQuMzk4NCAxLjg4NjY0IDE0LjExMTkgMS42MDAxIDEzLjc1ODQgMS42MDAxWiIgZmlsbD0iI2ZmZiIvPgo8cGF0aCBkPSJNNCAxMkwxMiA0TDQgMTJaIiBmaWxsPSIjZmZmIi8%2BCjxwYXRoIGQ9Ik00IDEyTDEyIDQiIHN0cm9rZT0iI2ZmZiIgc3Ryb2tlLXdpZHRoPSIxLjUiIHN0cm9rZS1saW5lY2FwPSJyb3VuZCIvPgo8L3N2Zz4K&logoColor=ffffff)](https://zread.ai/nikolasibalic/ARC-Alkali-Rydberg-Calculator)

**nblade** (N-dimensional Blade) is a high-performance geometric algebra library powered by Rust, with Python bindings. It supports arbitrary dimensions (up to 64D) and arbitrary metric signatures G(p, q, r).

## Features

- **Arbitrary Dimensions**: Up to 64-dimensional vector spaces
- **Arbitrary Signatures**: Support for G(p, q, r) metric signatures (Euclidean, spacetime, conformal, etc.)
- **High Performance**: Rust backend with parallel computing and SIMD optimization
- **Dual Representation**: Automatic selection of dense or sparse representation
- **Complete Operations**: All standard geometric algebra operations
- **NumPy Integration**: Create vectors directly from NumPy arrays

## Installation

### Python

```bash
pip install nblade
```

### From Source

```bash
# Install maturin
pip install maturin

# Build and install
maturin develop --release
```

## Quick Start

```python
import nblade

# Create 3D Euclidean geometric algebra G(3,0,0)
algebra = nblade.Algebra.euclidean(3)

# Create basis vectors
e1 = algebra.basis_vector(0)
e2 = algebra.basis_vector(1)
e3 = algebra.basis_vector(2)

# Create a vector from a list
v = algebra.vector([1, 2, 3])

# Geometric product
product = e1 * e2  # Results in bivector e12

# Outer product (wedge product)
wedge = e1 ^ e2  # Results in bivector e12

# Inner product
inner = e1 | e2  # Results in 0 (orthogonal vectors)

# Dual
I = algebra.config.volume_element()  # Pseudoscalar
v_dual = v.dual()

# Rotation using rotor
import math
plane = e1 ^ e2
rotor = algebra.rotor(plane, math.pi / 2)  # 90-degree rotation
rotated = e1.rotate_by(rotor)
```

## API Reference

### Algebra Class

```python
algebra = nblade.Algebra(dimension, p=0, q=0, r=0)

# Factory methods
algebra = nblade.Algebra.euclidean(dimension)  # G(n, 0, 0)
algebra = nblade.Algebra.spacetime(dimension)  # G(1, n-1, 0)
algebra = nblade.Algebra.cga()                 # G(4, 1, 0)

# Properties
algebra.dimension   # Vector space dimension
algebra.signature   # (p, q, r) tuple

# Methods
algebra.basis_vector(i)      # Create e_i
algebra.vector([x, y, z])    # Create vector from list
algebra.scalar(value)        # Create scalar multivector
algebra.zeros()              # Create zero multivector
algebra.one()                # Create unit scalar
algebra.rotor(plane, angle)  # Create rotor
```

### MultiVector Class

```python
# Construction
mv = nblade.MultiVector.basis_vector(config, i)
mv = nblade.MultiVector.from_scalar(config, value)
mv = nblade.MultiVector.from_coefficients(config, coeffs)
mv = nblade.MultiVector.zeros(config)
mv = nblade.MultiVector.one(config)

# Products
mv.geometric_product(other)  # Geometric product
mv.outer_product(other)      # Outer/wedge product
mv.left_inner(other)         # Left contraction
mv.right_inner(other)        # Right contraction

# Involutions
mv.grade_involution()        # Grade involution (A*)
mv.reversion()               # Reversion (A†)
mv.clifford_conjugate()      # Clifford conjugate (A‡)

# Other operations
mv.dual()                    # Dual (A⊥)
mv.inverse_dual()            # Inverse dual (A⁻⊥)
mv.inverse()                 # Multiplicative inverse
mv.norm()                    # Norm |A|
mv.norm_squared()            # |A|²

# Grade operations
mv.grade(r)                  # r-grade part
mv.even_part()               # Even grades
mv.odd_part()                # Odd grades

# Geometric operations
mv.project_to(blade)         # Project onto blade
mv.reject_from(blade)        # Reject from blade
mv.reflect_in(blade)         # Reflect in blade
mv.rotate_by(rotor)          # Rotate by rotor

# Operator overloads
mv1 + mv2   # Addition
mv1 - mv2   # Subtraction
mv1 * mv2   # Geometric product
mv1 ^ mv2   # Outer product
mv1 | mv2   # Left inner product
~mv        # Grade involution
-mv        # Negation
```

### Module Functions

```python
# Create rotor for rotation in a plane
rotor = nblade.create_rotor(plane, angle)

# Reciprocal frame computation
reciprocal = nblade.reciprocal_frame(vectors)

# Basis expansion
coeffs = nblade.basis_expansion(multivector)
```

## Supported Operations

| Operation | Symbol | Formula |
|-----------|--------|---------|
| Geometric Product | `AB` | `AB = A·B + A∧B` |
| Outer Product | `A∧B` | Antisymmetric part |
| Left Inner | `A⌋B` | Left contraction |
| Right Inner | `A⌊B` | Right contraction |
| Grade Involution | `A*` | `(-1)^r A_r` |
| Reversion | `A†` | `(-1)^(r(r-1)/2) A_r` |
| Clifford Conjugate | `A‡` | `(A*)†` |
| Dual | `A⊥` | `A·I` or `AI` |
| Inverse Dual | `A⁻⊥` | `A⌋I` |

## Performance

nblade is optimized for performance:

- **SIMD acceleration** for 2D-4D operations (with `simd` feature)
- **Adaptive parallelism** for high-dimensional operations (≥6D)
- **Memory pooling** for reduced allocation overhead (with `pool` feature)
- **Dense/Sparse auto-selection** based on coefficient density

### Benchmarks

| Operation | Dimension | Time |
|-----------|-----------|------|
| Geometric Product | 3D | ~5ms |
| Geometric Product | 5D | ~15ms |
| Geometric Product (SIMD) | 3D | ~3ms |

## Examples

### 2D Rotation

```python
import nblade
import math

# Create 2D Euclidean algebra
algebra = nblade.Algebra.euclidean(2)

# Create basis vectors
e1, e2 = algebra.basis_vectors()

# Create rotor for 45-degree rotation
plane = e1 ^ e2
rotor = algebra.rotor(plane, math.pi / 4)

# Rotate e1
rotated = e1.rotate_by(rotor)
print(rotated)  # ~ 0.707*e1 + 0.707*e2
```

### 3D Cross Product via Dual

```python
import nblade

algebra = nblade.Algebra.euclidean(3)
e1, e2, e3 = algebra.basis_vectors()

# Cross product: a × b = (a ∧ b)*
# where * is the dual in 3D

a = e1 + 2*e2
b = 2*e1 + e3

# Compute cross product
a_cross_b = (a ^ b).dual()
print(a_cross_b)
```

## License

MIT License