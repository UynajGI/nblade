# Physics Examples Overview

This page introduces application examples of nblade in physics.

## Example List

### 01_rigid_body.py - Rigid Body Dynamics

**Difficulty**: ⭐⭐⭐ Advanced

**Content**:
- Angular momentum as bivector
- Inertia tensor
- Rigid body rotation dynamics
- Euler equations in geometric algebra form
- Rotational kinetic energy
- Simple rigid body simulation

**Run**:
```bash
python examples/physics/01_rigid_body.py
```

**Expected Output**:
```
============================================================
nblade Rigid Body Physics Example
Application of Geometric Algebra in Classical Mechanics
============================================================

============================================================
Section 1: Angular Momentum as Bivector
============================================================

In traditional vector analysis, angular momentum L = r × p
In geometric algebra, angular momentum naturally represents as bivector L = r ∧ p

Position vector r = 3.000000e1 + 4.000000e2
Momentum vector p = 2.000000e1 + -1.000000e2 + 5.000000e3

Angular momentum bivector L = r∧p = -11.000000e1∧e2 + 15.000000e1∧e3 + 20.000000e2∧e3
...
```

---

## Core Concepts

### Geometric Representation of Angular Momentum

In traditional vector analysis, angular momentum is defined as cross product:

$$\mathbf{L} = \mathbf{r} \times \mathbf{p}$$

In geometric algebra, angular momentum naturally represents as bivector:

$$L = r \wedge p$$

Advantages of this representation:
- Works in arbitrary dimensions (cross product only works in 3D)
- Clear geometric meaning (spanned plane)
- Consistent with other geometric algebra operations

### Relationship with Traditional Cross Product

In 3D, the dual of angular momentum corresponds to traditional cross product:

```python
L_bivector = r ^ p        # Geometric algebra representation
L_vector = L_bivector.dual()  # Corresponds to traditional cross product result
```

### Rigid Body Rotation

Rigid body orientation is described by rotor:

```python
# Angular velocity plane
omega_plane = omega.dual()

# Update rotor
delta_rotor = alg.rotor(omega_plane, omega_magnitude * dt)
rotor = delta_rotor * rotor
```

---

## Theoretical Background

### Euler Equations

Euler equations for rigid body without external torque:

$$\frac{d\mathbf{L}}{dt} + \boldsymbol{\omega} \times \mathbf{L} = 0$$

In geometric algebra:

$$\frac{dL}{dt} + \omega \wedge L = 0$$

### Rotational Kinetic Energy

$$T = \frac{1}{2} \omega \cdot L = \frac{1}{2} \sum_i I_i \omega_i^2$$

In nblade:

```python
T = 0.5 * (omega | L).scalar_part()
```

---

## Frequently Asked Questions

### Q: Why use bivector to represent angular momentum?

A: Bivectors explicitly represent the rotation plane, which is more geometrically intuitive than traditional vectors (which represent rotation axis). Additionally, bivectors work in any dimension.

### Q: How to calculate rotor from angular velocity?

A: The dual of angular velocity is the rotation plane, then use that plane to create rotor:

```python
omega = alg.vector([0, 0, 1])  # Around z-axis
omega_plane = omega.dual()      # e1∧e2
rotor = alg.rotor(omega_plane, angle)
```

---

## References

- [Geometric Algebra for Physicists](https://www.cambridge.org/core/books/geometric-algebra-for-physicists/) - Doran & Lasenby
- [New Foundations for Classical Mechanics](https://www.springer.com/gp/book/9780792355141) - David Hestenes

---

## Next Steps

- [Graphics Examples](cg.md) - Geometric transformations
- [Tutorial Examples](tutorials.md) - Basic concepts
