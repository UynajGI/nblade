# Algebra Class API Reference

The `Algebra` class is the core entry point of nblade, used for creating and configuring geometric algebra spaces.

## Class Definition

```python
class Algebra:
    def __init__(self, dimension: int, p: int = 0, q: int = 0, r: int = 0) -> None
```

### Parameters

| Parameter | Type | Description |
|------|------|------|
| `dimension` | `int` | Vector space dimension (1-64) |
| `p` | `int` | Number of basis vectors with positive square (e_i² = +1) |
| `q` | `int` | Number of basis vectors with negative square (e_i² = -1) |
| `r` | `int` | Number of basis vectors with zero square (e_i² = 0) |

### Example

```python
import nblade

# Create custom signature algebra G(2,1,0)
alg = nblade.Algebra(dimension=3, p=2, q=1, r=0)
```

---

## Factory Methods

### `euclidean(dimension)`

Create Euclidean geometric algebra G(n, 0, 0).

```python
@classmethod
def euclidean(cls, dimension: int) -> Algebra
```

**Parameters**: `dimension` - Space dimension

**Returns**: Algebra instance configured with Euclidean signature

**Example**:
```python
# 3D Euclidean geometric algebra
alg = nblade.Algebra.euclidean(3)

# All basis vectors satisfy e_i² = +1
```

---

### `spacetime(dimension)`

Create spacetime algebra G(1, n-1, 0).

```python
@classmethod
def spacetime(cls, dimension: int) -> Algebra
```

**Parameters**: `dimension` - Spacetime dimension (typically 4)

**Returns**: Algebra instance configured with spacetime signature

**Example**:
```python
# Relativistic spacetime algebra G(1,3,0)
sta = nblade.Algebra.spacetime(4)

# e0² = +1 (time direction)
# e1² = e2² = e3² = -1 (space directions)
```

---

### `cga()`

Create Conformal Geometric Algebra G(4, 1, 0).

```python
@classmethod
def cga(cls) -> Algebra
```

**Returns**: Algebra instance configured for Conformal Geometric Algebra

**Example**:
```python
# Conformal Geometric Algebra, used to represent points, lines, circles, spheres, etc.
cga = nblade.Algebra.cga()
```

---

## Properties

### `dimension`

```python
@property
def dimension(self) -> int
```

Returns the dimension of the vector space.

**Example**:
```python
alg = nblade.Algebra.euclidean(3)
print(alg.dimension)  # 3
```

---

### `basis_count`

```python
@property
def basis_count(self) -> int
```

Returns the number of basis elements (2^dimension).

**Example**:
```python
alg = nblade.Algebra.euclidean(3)
print(alg.basis_count)  # 8 (scalar, 3 vectors, 3 bivectors, 1 trivector)
```

---

### `signature`

```python
@property
def signature(self) -> Tuple[int, int, int]
```

Returns the metric signature (p, q, r).

**Example**:
```python
alg = nblade.Algebra.spacetime(4)
print(alg.signature)  # (1, 3, 0)
```

---

### `config`

```python
@property
def config(self) -> AlgebraConfig
```

Returns the underlying AlgebraConfig object for low-level API calls.

---

## Vector Creation Methods

### `vector(data)`

Create a vector (1-vector) from a list or tuple.

```python
def vector(self, data: Union[List[float], Tuple[float, ...]]) -> MultiVector
```

**Parameters**: `data` - List or tuple of length `dimension`

**Returns**: MultiVector representing the vector

**Exceptions**: `ValueError` - If data length does not match dimension

**Example**:
```python
alg = nblade.Algebra.euclidean(3)

# Create from list
v = alg.vector([1.0, 2.0, 3.0])  # 1*e1 + 2*e2 + 3*e3

# Create from tuple
w = alg.vector((4.0, 5.0, 6.0))
```

---

### `basis_vector(i)`

Create the i-th basis vector e_i.

```python
def basis_vector(self, i: int) -> MultiVector
```

**Parameters**: `i` - Basis vector index (0 to dimension-1)

**Returns**: Basis vector MultiVector

**Example**:
```python
alg = nblade.Algebra.euclidean(3)
e1 = alg.basis_vector(0)  # First basis vector
e2 = alg.basis_vector(1)  # Second basis vector
e3 = alg.basis_vector(2)  # Third basis vector
```

---

### `basis_vectors()`

Get all basis vectors.

```python
def basis_vectors(self) -> List[MultiVector]
```

**Returns**: List of basis vectors

**Example**:
```python
alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()
```

---

### `scalar(value)`

Create a scalar multivector.

```python
def scalar(self, value: float) -> MultiVector
```

**Parameters**: `value` - Scalar value

**Returns**: Scalar MultiVector

**Example**:
```python
s = alg.scalar(5.0)  # Pure scalar
```

---

### `one()`

Create the unit multivector (scalar 1).

```python
def one(self) -> MultiVector
```

**Returns**: Unit scalar MultiVector

---

### `zeros()`

Create the zero multivector.

```python
def zeros(self) -> MultiVector
```

**Returns**: Zero-valued MultiVector

---

### `from_coefficients(coefficients)`

Create a multivector from a coefficient array.

```python
def from_coefficients(self, coefficients: List[float]) -> MultiVector
```

**Parameters**: `coefficients` - Coefficient array for all basis elements, length 2^dimension

**Returns**: MultiVector

**Example**:
```python
# Coefficient order for 3D algebra: [1, e1, e2, e3, e12, e13, e23, e123]
coeffs = [1.0, 2.0, 3.0, 4.0, 0.5, 0.0, 0.0, 0.1]
mv = alg.from_coefficients(coeffs)
```

---

## Rotor Creation

### `rotor(plane, angle)`

Create a rotor (rotator).

```python
def rotor(self, plane: MultiVector, angle: float) -> MultiVector
```

**Parameters**:
- `plane` - Rotation plane (bivector)
- `angle` - Rotation angle (radians)

**Returns**: Rotor MultiVector

**Example**:
```python
import math

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# Rotate 45 degrees in xy plane
plane = e1 ^ e2
rotor = alg.rotor(plane, math.pi / 4)

# Rotate vector using rotor
rotated = e1.rotate_by(rotor)
```

---

## String Representation

### `__repr__()`

Returns detailed algebra information.

```python
alg = nblade.Algebra.euclidean(3)
print(repr(alg))  # Algebra(dimension=3, signature=(3, 0, 0))
```

### `__str__()`

Returns human-readable algebra description.

```python
alg = nblade.Algebra.euclidean(3)
print(str(alg))  # Geometric Algebra G(3, 0, 0) (3D)
```

---

## Complete Example

```python
import nblade
import math

# Create 3D Euclidean algebra
alg = nblade.Algebra.euclidean(3)

# Get basis vectors
e1, e2, e3 = alg.basis_vectors()

# Create vector
v = alg.vector([1.0, 2.0, 3.0])

# Basic operations
geometric = e1 * e2  # Geometric product
outer = e1 ^ e2      # Outer product
inner = e1 | e2      # Inner product

# Create rotor and rotate
rotor = alg.rotor(e1 ^ e2, math.pi / 4)
rotated = v.rotate_by(rotor)

print(f"Original vector: {v}")
print(f"Rotated: {rotated}")
```
