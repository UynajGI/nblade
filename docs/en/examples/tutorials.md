# Tutorial Examples Overview

This page introduces nblade's tutorial examples to help you quickly master the basic concepts and operations of geometric algebra.

## Example List

### 01_quickstart.py - Quick Start

**Difficulty**: ⭐ Beginner

**Content**:
- Creating geometric algebra
- Basis vectors and vector creation
- Three basic products
- Introduction to rotation operations

**Run**:
```bash
python examples/tutorials/01_quickstart.py
```

**Expected Output**:
```
============================================================
nblade Quick Start - Master Core Features in 5 Minutes
============================================================

【1. Creating Geometric Algebra】
Algebra Type: Geometric Algebra G(3, 0, 0) (3D)
Dimension: 3
Signature (p, q, r): (3, 0, 0)

【2. Basis Vectors and Vector Creation】
Basis vector e1 = 1.000000e1
...
```

**Key Concepts**:
- `nblade.Algebra.euclidean(n)` - Create n-dimensional Euclidean algebra
- `alg.vector([...])` - Create vector
- `*`, `^`, `|` - Geometric product, outer product, inner product

---

### 02_basic_operations.py - Basic Operations

**Difficulty**: ⭐⭐ Basic

**Content**:
- Detailed geometric product explanation
- Detailed outer product explanation
- Detailed inner product explanation
- Relationship between the three products
- Practical application examples

**Run**:
```bash
python examples/tutorials/02_basic_operations.py
```

**Key Concepts**:
- Core formula: `ab = a·b + a∧b`
- Outer product represents the subspace spanned
- Inner product represents projection relationship

---

### 03_vectors.py - Vector Operations

**Difficulty**: ⭐⭐ Basic

**Content**:
- Various methods for creating vectors
- Vector properties (norm, coefficients, etc.)
- Algebraic operations (addition, subtraction, multiplication)
- Geometric operations (projection, reflection)
- Involution operations
- High-dimensional vectors

**Run**:
```bash
python examples/tutorials/03_vectors.py
```

**Key Concepts**:
- `alg.basis_vectors()` - Get all basis vectors
- `v.norm()` - Calculate norm
- `v.project_to(blade)` - Projection
- `v.reflect_in(n)` - Reflection

---

### 05_rotations.py - Rotation Tutorial

**Difficulty**: ⭐⭐⭐ Advanced

**Content**:
- 2D rotation and complex numbers
- 3D rotation and rotors
- Combining multiple rotations
- Rotation interpolation (SLERP)
- Relationship between rotors and quaternions
- Practical application examples

**Run**:
```bash
python examples/tutorials/05_rotations.py
```

**Key Concepts**:
- Rotor: `R = exp(-B·θ/2)`
- Rotation formula: `v' = R·v·R†`
- Combined rotation: `R_combined = R2·R1`

---

## Learning Path

```
01_quickstart.py
       ↓
02_basic_operations.py
       ↓
03_vectors.py
       ↓
05_rotations.py
       ↓
   Application Examples
```

## Frequently Asked Questions

### Q: Why is the outer product result a bivector?

A: The outer product `a ^ b` represents the subspace (plane) spanned by two vectors. In geometric algebra, planes are represented by bivectors, not by normal vectors as in traditional vector analysis.

### Q: What is the relationship between geometric product and inner/outer product?

A: Geometric product is the fundamental operation, decomposed into inner and outer products: `ab = a·b + a∧b`. Inner product is the symmetric part (scalar), outer product is the antisymmetric part (bivector).

### Q: What is the relationship between rotors and quaternions?

A: Quaternions are the even subalgebra of geometric algebra. Quaternion units i, j, k correspond to bivectors e₂∧e₃, e₃∧e₁, e₁∧e₂.

---

## Next Steps

After completing the tutorial examples, you can continue learning:

- [Physics Examples](physics.md) - Rigid body dynamics
- [Graphics Examples](cg.md) - Geometric transformations

Or check the [API Reference](../api/algebra.md) for more details.
