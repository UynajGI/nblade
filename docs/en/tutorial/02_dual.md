# Dual Operations Tutorial

This tutorial explains the concept of dual in geometric algebra and its applications.

## Table of Contents

1. [What is the Dual?](#1-what-is-the-dual)
2. [Mathematical Definition](#2-mathematical-definition)
5. [Dual and Cross Product](#3-dual-and-cross-product)
4. [Inverse Dual](#4-inverse-dual)
5. [Applications](#5-applications)

---

## 1. What is the Dual?

The dual is a fundamental operation in geometric algebra that maps between different grades of multivectors. In 3D Euclidean space:

| Original | Dual |
|----------|------|
| Scalar | Trivector (pseudoscalar) |
| Vector | Bivector |
| Bivector | Vector |
| Trivector | Scalar |

The dual essentially gives us the "orthogonal complement" of a multivector.

### Intuitive Understanding

In 3D:
- The dual of a **vector** is a **bivector** representing the plane perpendicular to that vector
- The dual of a **bivector** is a **vector** perpendicular to that plane
- The dual of a **scalar** is the **volume element** (pseudoscalar)

## 2. Mathematical Definition

### Pseudoscalar

The pseudoscalar `I` is the highest-grade element in an algebra:

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# Pseudoscalar I = e1 ∧ e2 ∧ e3
I = e1 ^ e2 ^ e3
print(f"I = {I}")  # e1∧e2∧e3
print(f"I² = {(I * I).scalar_part()}")  # -1 in 3D Euclidean
```

### Dual Operation

The dual of a multivector `A` is defined as:

```
A* = A · I⁻¹  (right contraction)
```

Or equivalently:

```
A* = A I⁻¹   (geometric product)
```

### Using Dual in nblade

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# Dual of a vector
v = e1
v_dual = v.dual()
print(f"Dual of e1: {v_dual}")  # e2∧e3

# Dual of a bivector
B = e1 ^ e2
B_dual = B.dual()
print(f"Dual of e1∧e2: {B_dual}")  # e3
```

## 3. Dual and Cross Product

In 3D, the cross product can be expressed using the dual:

```
a × b = (a ∧ b)*
```

This is why the cross product only works in 3D - it requires the dual to map a bivector back to a vector.

### Example: Cross Product via Dual

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

a = e1
b = e2

# Cross product via dual: a × b = (a ∧ b)*
wedge = a ^ b        # a ∧ b = e1∧e2
cross = wedge.dual() # dual = e3

print(f"a ∧ b = {wedge}")  # e1∧e2
print(f"(a ∧ b)* = {cross}")  # e3

# This is equivalent to the cross product!
# e1 × e2 = e3
```

### Advantages of the GA Approach

| Cross Product | Geometric Algebra Dual |
|---------------|------------------------|
| Only works in 3D | Works in any dimension |
| Returns a vector | Returns a bivector (more natural) |
| Non-associative | Wedge product is associative |

## 4. Inverse Dual

The inverse dual reverses the dual operation:

```
(A*)* = A (up to sign in some signatures)
```

### Using Inverse Dual in nblade

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

v = e1 + 2*e2 + 3*e3

# Dual
v_dual = v.dual()

# Inverse dual (should recover original)
v_recovered = v_dual.inverse_dual()

print(f"Original: {v}")
print(f"Dual: {v_dual}")
print(f"Inverse dual: {v_recovered}")
```

## 5. Applications

### Normal Vectors from Planes

A plane defined by a bivector B has a normal vector n = B*:

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# Plane defined by e1∧e2 (xy-plane)
plane = e1 ^ e2

# Normal vector (points in z direction)
normal = plane.dual()
print(f"Normal to xy-plane: {normal}")  # e3
```

### Area and Volume

The dual naturally relates areas and volumes:

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# Area element
area = e1 ^ e2  # Bivector representing area

# "Volume" of the area (magnitude)
# |B*| = |B| in 3D Euclidean
area_magnitude = area.norm()
print(f"Area magnitude: {area_magnitude}")  # 1.0
```

### Electromagnetic Field

In physics, the electromagnetic field F = E + iB uses the dual to relate E and B fields:

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()
I = e1 ^ e2 ^ e3  # Pseudoscalar

# Electric field
E = e1 + e2

# Magnetic field as bivector via dual
# B_bivector = I * B_vector = dual(B_vector)
B_vector = e3
B_bivector = B_vector.dual()  # e1∧e2

print(f"Magnetic field bivector: {B_bivector}")
```

---

## Summary

| Concept | Formula | nblade Method |
|---------|---------|---------------|
| Dual | A* = AI⁻¹ | `A.dual()` |
| Inverse Dual | (A*)* | `A.inverse_dual()` |
| Cross Product (3D) | a × b = (a ∧ b)* | `(a ^ b).dual()` |
| Pseudoscalar | I = e₁∧e₂∧...∧eₙ | `alg.config.volume_element()` |

## Further Reading

- Run the example: `python examples/tutorials/04_dual_operations.py`
- See also: [Grade Operations](./03_grade.md)
- API Reference: [MultiVector.dual()](../api/multivector.md#dual)