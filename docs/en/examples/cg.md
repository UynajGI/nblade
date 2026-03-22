# Computer Graphics Examples Overview

This page introduces application examples of nblade in computer graphics.

## Example List

### 01_transformations.py - Geometric Transformations

**Difficulty**: ⭐⭐⭐ Advanced

**Content**:
- Reflection transformations
- Projection operations
- Rotation operations
- Scaling transformations
- Combined transformations
- Triangle transformation example
- Introduction to conformal geometric algebra

**Run**:
```bash
python examples/cg/01_transformations.py
```

**Expected Output**:
```
============================================================
nblade Computer Graphics Transformation Example
============================================================

============================================================
Section 1: Reflection
============================================================

Reflection is the most fundamental geometric transformation
Formula: v' = -n·v·n (where n is unit normal vector)

【Vector Reflection】
Original vector v = 1.000000e1 + 2.000000e2 + 3.000000e3
Reflection plane normal n = 1.000000e1 (yz plane)
Reflection result v' = 1.000000e1 + -2.000000e2 + -3.000000e3
...
```

---

## Core Concepts

### Reflection

Reflection is the most fundamental geometric transformation. Other transformations can be obtained through combinations of reflections.

**Formula**: `v' = -n v n`

Where `n` is the unit normal vector of the reflection plane.

```python
v = alg.vector([1, 2, 3])
n = e1  # Normal vector of yz plane

reflected = v.reflect_in(n)
```

### Rotation

Rotation is implemented using rotors.

**Formula**: `v' = R v R†`

```python
import math

# Create rotor
plane = e1 ^ e2  # xy plane
rotor = alg.rotor(plane, math.pi / 4)  # 45°

# Rotate
rotated = v.rotate_by(rotor)
```

### Projection

Project vector onto subspace.

```python
v = alg.vector([1, 2, 3])

# Project to x-axis
proj = v.project_to(e1)

# Orthogonal component
reject = v.reject_from(e1)

# v = proj + reject
```

### Combined Transformations

Multiple transformations are combined through geometric product.

```python
# Rotate first, then reflect
result = v.rotate_by(rotor).reflect_in(n)

# Transformation order matters!
result2 = v.reflect_in(n).rotate_by(rotor)  # Different result
```

---

## Comparison with Traditional Methods

| Operation | Traditional Method | Geometric Algebra Method |
|-----------|-------------------|-------------------------|
| Rotation | 3×3 matrix | Rotor (4 parameters) |
| Reflection | 3×3 matrix | Vector multiplication |
| Combination | Matrix multiplication | Geometric product |
| Interpolation | Euler angles/Quaternions | Rotor SLERP |

**Advantages**:
- Fewer parameters (rotor vs rotation matrix)
- Better numerical stability
- Avoid gimbal lock
- Unified operation rules

---

## Conformal Geometric Algebra (CGA)

Conformal geometric algebra G(4,1,0) extends 3D Euclidean algebra and can represent:

- Points, lines, planes
- Circles, spheres
- Rigid body motions (including translation)

```python
# Create CGA
cga = nblade.Algebra.cga()

# In CGA, translation is also rotor
# Unifies rotation and translation
```

### Advantages of CGA

1. **All geometric objects are multivectors**
2. **Intersection operations through outer product**
3. **Translation and rotation unified as rotors**
4. **Concise formulas, avoid special cases**

---

## Practical Tips

### Batch Transforming Point Sets

```python
# Create multiple points
points = [e1, e2, e3, e1 + e2]

# Create rotor
rotor = alg.rotor(e1 ^ e2, math.pi / 4)

# Batch rotation
rotated_points = [p.rotate_by(rotor) for p in points]
```

### Determine Reflection Plane from Two Points

```python
# Two points determine a line, the perpendicular plane of the line can be used as reflection plane
p1 = alg.vector([1, 0, 0])
p2 = alg.vector([0, 1, 0])

line = p2 - p1
# Reflection plane normal vector needs to be specified externally
```

### Determine Which Side of Plane a Point Is On

```python
point = alg.vector([1, 2, 3])
plane_normal = e3  # xy plane

# Sign of dot product indicates which side
side = (point | plane_normal).scalar_part()
if side > 0:
    print("On positive side of plane")
elif side < 0:
    print("On negative side of plane")
else:
    print("On plane")
```

---

## Frequently Asked Questions

### Q: What is the difference between rotor and quaternion?

A: Quaternions are the even subalgebra of geometric algebra. Rotors are a more general concept that can represent rotations in any dimension. In 3D, they are equivalent.

### Q: How to implement translation?

A: In standard geometric algebra, translation is implemented through vector addition. In conformal geometric algebra, translation is also a rotor, unified with rotation.

### Q: How to do rotation interpolation?

A: Use SLERP (Spherical Linear Interpolation) to interpolate rotors:

```python
def slerp(R1, R2, t):
    # Simplified version
    R = R1.scale(1-t) + R2.scale(t)
    # Normalization...
    return R
```

---

## References

- [Geometric Algebra for Computer Science](https://www.elsevier.com/books/geometric-algebra-for-computer-science/dorst/978-0-12-374942-0) - Dorst, Fontijne, Mann
- [GA Wiki - Conformal Geometric Algebra](https://en.wikipedia.org/wiki/Conformal_geometric_algebra)

---

## Next Steps

- [Tutorial Examples](tutorials.md) - Basic concepts
- [API Reference](../api/algebra.md) - Detailed API
