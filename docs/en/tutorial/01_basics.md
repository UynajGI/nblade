# Geometric Algebra Basics Tutorial

This tutorial introduces the fundamental concepts of geometric algebra and how to use nblade.

## Table of Contents

1. [Introduction to Geometric Algebra](#1-introduction-to-geometric-algebra)
2. [Multivectors](#2-multivectors)
3. [Geometric Product](#3-geometric-product)
4. [Outer Product](#4-outer-product)
5. [Inner Product](#5-inner-product)
6. [Involutions](#6-involutions)
7. [Dual](#7-dual)

---

## 1. Introduction to Geometric Algebra

### What is Geometric Algebra?

Geometric Algebra (also known as Clifford Algebra) is a mathematical framework developed by Hermann Grassmann and William Kingdon Clifford in the 19th century. It unifies:

- **Scalars** (0-dimensional)
- **Vectors** (1-dimensional)
- **Bivectors** (2-dimensional, representing planes)
- **Trivectors** (3-dimensional, representing volumes)
- ... and so on

into a single concept called **Multivector**.

### Why Use Geometric Algebra?

| Traditional Approach | Geometric Algebra Approach |
|---------------------|---------------------------|
| Use matrices for rotation | Use rotors, more elegant |
| Cross product only works in 3D | Outer product works in any dimension |
| Complex numbers, quaternions are separate concepts | Both are subalgebras of geometric algebra |
| Multiple operation rules needed | Geometric product unifies everything |

## 2. Multivectors

A multivector is the fundamental object in geometric algebra, a linear combination of elements of different grades.

### Grade

- **Grade 0**: Scalar, e.g., `5`
- **Grade 1**: Vector, e.g., `e1`, `2e1 + 3e2`
- **Grade 2**: Bivector, e.g., `e1∧e2` (represents a plane)
- **Grade 3**: Trivector, e.g., `e1∧e2∧e3` (represents a volume)

### Creating Multivectors in nblade

```python
import nblade

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# Scalar
s = alg.scalar(5.0)

# Vector
v = alg.vector([1.0, 2.0, 3.0])

# Bivector
B = e1 ^ e2

# Mixed multivector
mixed = s + v + B
```

### Extracting Specific Grades

```python
# Get scalar part
scalar_part = mixed.grade(0)

# Get vector part
vector_part = mixed.grade(1)

# Get bivector part
bivector_part = mixed.grade(2)
```

## 3. Geometric Product

The geometric product is the most fundamental operation in geometric algebra.

### Definition

For two vectors `a` and `b`:

```
ab = a·b + a∧b
```

- `a·b`: Inner product, the symmetric part (scalar)
- `a∧b`: Outer product, the antisymmetric part (bivector)

### Using in nblade

```python
a = alg.vector([1.0, 2.0, 0.0])
b = alg.vector([3.0, 1.0, 0.0])

# Geometric product
product = a * b
print(product)  # scalar + bivector
```

### Properties

1. **Associative**: `(ab)c = a(bc)`
2. **Distributive**: `a(b + c) = ab + ac`
3. **Non-commutative**: `ab ≠ ba` (in general)

### Geometric Product of Basis Vectors

```python
# Orthogonal vectors
print(e1 * e2)  # e12 (bivector)

# Parallel vectors
print(e1 * e1)  # 1 (scalar)
```

## 4. Outer Product

The outer product (wedge product) represents the subspace spanned by vectors.

### Geometric Meaning

- The outer product `a∧b` of vectors `a` and `b` is a **bivector**
- It represents the plane spanned by `a` and `b`
- Its magnitude equals the area of the parallelogram

### Antisymmetry

```python
print(e1 ^ e2)   # e12
print(e2 ^ e1)   # -e12
print(e1 ^ e1)   # 0
```

### Higher-Order Outer Products

```python
# Trivector (volume)
volume = e1 ^ e2 ^ e3
print(volume)  # e123

# In 3D, this is the pseudoscalar
```

### Computing Area

```python
a = alg.vector([2.0, 0.0, 0.0])
b = alg.vector([1.0, 1.0, 0.0])

area = (a ^ b).norm()
print(f"Parallelogram area: {area}")  # 2.0
```

## 5. Inner Product

The inner product represents projection relationships.

### Inner Product of Vectors

```python
a = alg.vector([1.0, 2.0, 3.0])
b = alg.vector([4.0, 5.0, 6.0])

# Inner product (result is scalar)
inner = a | b
print(inner)  # 32
```

### Left and Right Contraction

```python
# Left inner product (left contraction)
left = a | b

# Right inner product (right contraction)
right = a.right_inner(b)
```

### Geometric Meaning

The inner product can be used to:

1. **Compute projection**: `(a|b)/|b|² × b` is the projection of `a` onto `b`
2. **Test orthogonality**: If `a|b = 0`, then `a` and `b` are orthogonal
3. **Compute angle**: `cos(θ) = (a|b)/(|a||b|)`

## 6. Involutions

Involutions are operations that map multivectors to themselves.

### Grade Involution

For an r-grade element, multiply by `(-1)^r`:

```python
v = alg.vector([1, 2, 3])
print(v.grade_involution())  # -v (vectors are grade 1, (-1)^1 = -1)
```

### Reversion

For an r-grade element, multiply by `(-1)^(r(r-1)/2)`:

```python
print(v.reversion())  # v (vectors unchanged)
```

### Clifford Conjugate

Clifford conjugate = grade involution + reversion:

```python
print(v.clifford_conjugate())
```

### Applications

Involutions are mainly used for computing:

- **Norm**: `|A|² = ⟨A†A⟩₀`
- **Inverse**: `A⁻¹ = A†/|A|²` (for invertible elements)

## 7. Dual

The dual maps a k-vector to an (n-k)-vector.

### Definition

In n-dimensional space, the dual of vector `v` is:

```
v* = v·I⁻¹
```

where `I` is the pseudoscalar (highest-grade basis element).

### Using in nblade

```python
v = alg.vector([1, 2, 3])
dual_v = v.dual()
print(dual_v)  # bivector

# Inverse dual
original = dual_v.inverse_dual()
```

### Geometric Meaning

In 3D:

- The dual of a vector is a bivector (plane)
- The dual of a bivector is a vector
- This establishes the connection between **cross product** and **outer product**: `a × b = (a ∧ b)*`

### Pseudoscalar

```python
# Get pseudoscalar
I = alg.config.volume_element()
print(I)  # e123

# Square of pseudoscalar
print(I * I)  # -1 (in 3D Euclidean space)
```

---

## Next Steps

- [Rotation Tutorial](../examples/tutorials.md) - Learn how to use rotors for rotation
- [API Reference](../api/algebra.md) - View complete API documentation
- [Example Code](https://github.com/nblade/nblade/tree/main/examples) - More practical examples