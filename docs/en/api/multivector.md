# MultiVector Class API Reference

`MultiVector` is the core data type in geometric algebra, representing multivectors (linear combinations of scalars, vectors, bivectors, etc.).

## Class Definition

```python
class MultiVector:
    """
    Multivector class, representing multivector elements in geometric algebra.
    Can be a linear combination of scalars, vectors, bivectors, trivectors, etc.
    """
```

---

## Constructor Methods

### `zeros(config)`

Create a zero multivector.

```python
@classmethod
def zeros(cls, config: AlgebraConfig) -> MultiVector
```

**Parameters**: `config` - Algebra configuration object

**Returns**: Zero-valued MultiVector

---

### `one(config)`

Create the unit multivector (scalar 1).

```python
@classmethod
def one(cls, config: AlgebraConfig) -> MultiVector
```

**Returns**: Unit scalar MultiVector

---

### `from_scalar(config, scalar)`

Create a multivector from a scalar value.

```python
@classmethod
def from_scalar(cls, config: AlgebraConfig, scalar: float) -> MultiVector
```

**Parameters**:
- `config` - Algebra configuration
- `scalar` - Scalar value

**Returns**: Scalar MultiVector

**Example**:
```python
alg = nblade.Algebra.euclidean(3)
s = nblade.MultiVector.from_scalar(alg.config, 5.0)
# Or use the convenience method
s = alg.scalar(5.0)
```

---

### `basis_vector(config, i)`

Create basis vector e_i.

```python
@classmethod
def basis_vector(cls, config: AlgebraConfig, i: int) -> MultiVector
```

**Parameters**:
- `config` - Algebra configuration
- `i` - Basis vector index (0-based)

**Returns**: Basis vector MultiVector

---

### `from_coefficients(config, coefficients)`

Create a multivector from a coefficient array.

```python
@classmethod
def from_coefficients(cls, config: AlgebraConfig, coefficients: List[float]) -> MultiVector
```

**Parameters**:
- `config` - Algebra configuration
- `coefficients` - Coefficient array for all basis elements

**Returns**: MultiVector

**Example**:
```python
# Coefficient order in 3D: [scalar, e1, e2, e3, e12, e13, e23, e123]
coeffs = [1.0, 2.0, 3.0, 4.0, 0.5, 0.0, 0.0, 0.1]
mv = nblade.MultiVector.from_coefficients(alg.config, coeffs)
```

---

### `from_numpy(config, array)`

Create a multivector from a NumPy array.

```python
@classmethod
def from_numpy(cls, config: AlgebraConfig, array: numpy.ndarray) -> MultiVector
```

**Parameters**:
- `config` - Algebra configuration
- `array` - NumPy array

**Returns**: MultiVector

---

## Properties

### `config`

```python
@property
def config(self) -> AlgebraConfig
```

Returns the associated algebra configuration.

---

## Basic Operations

### Geometric Product `geometric_product()`

```python
def geometric_product(self, other: MultiVector) -> MultiVector
```

Compute the geometric product with another multivector.

**Operator**: `*`

**Example**:
```python
a = alg.vector([1, 0, 0])
b = alg.vector([0, 1, 0])

# Method call
product = a.geometric_product(b)

# Or use operator
product = a * b
```

---

### Outer Product `outer_product()`

```python
def outer_product(self, other: MultiVector) -> MultiVector
```

Compute the outer product (wedge product).

**Operator**: `^`

**Example**:
```python
a = alg.vector([1, 0, 0])
b = alg.vector([0, 1, 0])

wedge = a ^ b  # Bivector
```

---

### Left Inner Product `left_inner()`

```python
def left_inner(self, other: MultiVector) -> MultiVector
```

Compute the left inner product (left contraction).

**Operator**: `|`

**Example**:
```python
a = alg.vector([1, 2, 3])
b = alg.vector([4, 5, 6])

inner = a | b  # Scalar
```

---

### Right Inner Product `right_inner()`

```python
def right_inner(self, other: MultiVector) -> MultiVector
```

Compute the right inner product (right contraction).

**Example**:
```python
a = alg.vector([1, 2, 3])
b = alg.vector([4, 5, 6])

right = a.right_inner(b)
```

---

## Grade Operations

### `grade(r)`

Extract the r-grade part.

```python
def grade(self, r: int) -> MultiVector
```

**Parameters**: `r` - Grade (0=scalar, 1=vector, 2=bivector, ...)

**Returns**: MultiVector of specified grade

**Example**:
```python
# Mixed multivector
mv = alg.scalar(1.0) + e1 + (e1 ^ e2)

scalar = mv.grade(0)    # Scalar part
vector = mv.grade(1)    # Vector part
bivector = mv.grade(2)  # Bivector part
```

---

### `even_part()`

Extract the even-grade part.

```python
def even_part(self) -> MultiVector
```

**Returns**: Even-grade components (scalar + bivector + quadvector + ...)

---

### `odd_part()`

Extract the odd-grade part.

```python
def odd_part(self) -> MultiVector
```

**Returns**: Odd-grade components (vector + trivector + ...)

---

### `scalar_part()`

Extract the scalar part.

```python
def scalar_part(self) -> float
```

**Returns**: Scalar value

**Example**:
```python
a = alg.vector([3, 4, 0])
b = alg.vector([3, 4, 0])

inner = a | b
print(inner.scalar_part())  # 25
```

---

### `coefficients()`

Get coefficients for all basis elements.

```python
def coefficients(self) -> List[float]
```

**Returns**: List of coefficients

---

## Involution Operations

### `grade_involution()`

Grade involution.

```python
def grade_involution(self) -> MultiVector
```

Multiply r-grade elements by (-1)^r.

**Operator**: `~` (Note: In Python ~ is bitwise NOT, so use the method)

**Example**:
```python
v = alg.vector([1, 2, 3])
involution = v.grade_involution()  # -v (vector is grade 1)
```

---

### `reversion()`

Reversion.

```python
def reversion(self) -> MultiVector
```

Multiply r-grade elements by (-1)^(r(r-1)/2).

**Example**:
```python
v = alg.vector([1, 2, 3])
rev = v.reversion()  # v (vectors are unchanged)

B = e1 ^ e2
B_rev = B.reversion()  # -B (bivectors flip sign)
```

---

### `clifford_conjugate()`

Clifford conjugate.

```python
def clifford_conjugate(self) -> MultiVector
```

Equal to grade involution + reversion.

**Example**:
```python
mv = alg.scalar(1) + e1 + (e1 ^ e2)
conj = mv.clifford_conjugate()
```

---

## Dual Operations

### `dual()`

Hodge dual.

```python
def dual(self) -> MultiVector
```

Maps k-vectors to (n-k)-vectors.

**Example**:
```python
v = alg.vector([1, 2, 3])
dual_v = v.dual()  # Bivector

# In 3D, the dual establishes the connection between cross product and outer product
# a × b = (a ∧ b).dual()
```

---

### `inverse_dual()`

Inverse dual.

```python
def inverse_dual(self) -> MultiVector
```

**Example**:
```python
v = alg.vector([1, 2, 3])
original = v.dual().inverse_dual()  # Equals v
```

---

## Norm and Inverse

### `norm()`

Compute the norm.

```python
def norm(self) -> float
```

**Returns**: |A| = sqrt(<A†A>₀)

**Example**:
```python
v = alg.vector([3, 4, 0])
print(v.norm())  # 5.0
```

---

### `norm_squared()`

Compute the squared norm.

```python
def norm_squared(self) -> float
```

**Returns**: |A|²

---

### `inverse()`

Compute the multiplicative inverse.

```python
def inverse(self) -> MultiVector
```

**Returns**: A⁻¹ such that A * A⁻¹ = 1

**Exceptions**: Raises exception if not invertible

**Example**:
```python
v = alg.vector([1, 2, 3])
v_inv = v.inverse()

# Verify
product = v * v_inv
print(product.scalar_part())  # Approximately 1
```

---

### `is_invertible()`

Check if invertible.

```python
def is_invertible(self) -> bool
```

---

## Geometric Operations

### `project_to(blade)`

Project onto the subspace represented by a blade.

```python
def project_to(self, blade: MultiVector) -> MultiVector
```

**Parameters**: `blade` - Blade representing the target subspace

**Returns**: Projected MultiVector

**Example**:
```python
v = alg.vector([1, 2, 3])
e1 = alg.basis_vector(0)

proj = v.project_to(e1)  # Project onto x-axis
print(proj)  # e1
```

---

### `reject_from(blade)`

Reject from the subspace represented by a blade.

```python
def reject_from(self, blade: MultiVector) -> MultiVector
```

**Returns**: Orthogonal component

**Example**:
```python
v = alg.vector([1, 2, 3])
e1 = alg.basis_vector(0)

reject = v.reject_from(e1)  # Component perpendicular to x-axis
# v = proj + reject
```

---

### `reflect_in(blade)`

Reflect in the hyperplane represented by a blade.

```python
def reflect_in(self, blade: MultiVector) -> MultiVector
```

**Formula**: v' = -n v n (where n is the normal vector)

**Example**:
```python
v = alg.vector([1, 2, 3])
e1 = alg.basis_vector(0)

reflected = v.reflect_in(e1)  # Reflect in yz plane
# x component is negated
```

---

### `rotate_by(rotor)`

Rotate using a rotor.

```python
def rotate_by(self, rotor: MultiVector) -> MultiVector
```

**Parameters**: `rotor` - Rotation rotor

**Returns**: Rotated MultiVector

**Formula**: v' = R v R†

**Example**:
```python
import math

alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# Create rotor
plane = e1 ^ e2
rotor = alg.rotor(plane, math.pi / 4)

# Rotate
rotated = e1.rotate_by(rotor)
```

---

## Other Operations

### `scale(factor)`

Scalar multiplication.

```python
def scale(self, factor: float) -> MultiVector
```

**Parameters**: `factor` - Scale factor

**Returns**: Scaled MultiVector

---

### `add(other)`

Addition.

```python
def add(self, other: MultiVector) -> MultiVector
```

**Operator**: `+`

---

### `commutator(other)`

Commutator [A, B] = (AB - BA) / 2.

```python
def commutator(self, other: MultiVector) -> MultiVector
```

---

### `scalar_product(other)`

Scalar product <AB>₀.

```python
def scalar_product(self, other: MultiVector) -> float
```

---

### `is_zero()`

Check if zero.

```python
def is_zero(self) -> bool
```

---

## Operator Overloads

| Operator | Method | Description |
|--------|------|------|
| `a + b` | `add` | Addition |
| `a - b` | `__sub__` | Subtraction |
| `a * b` | `geometric_product` | Geometric product |
| `a ^ b` | `outer_product` | Outer product |
| `a \| b` | `left_inner` | Left inner product |
| `-a` | `__neg__` | Negation |
| `~a` | `grade_involution` | Grade involution |

---

## Complete Example

```python
import nblade
import math

# Create algebra
alg = nblade.Algebra.euclidean(3)
e1, e2, e3 = alg.basis_vectors()

# Create vectors
v = alg.vector([1.0, 2.0, 3.0])
w = alg.vector([4.0, 5.0, 6.0])

# Basic operations
print("Geometric product:", v * w)
print("Outer product:", v ^ w)
print("Inner product:", v | w)

# Grade operations
mixed = alg.scalar(1.0) + v + (e1 ^ e2)
print("Scalar part:", mixed.grade(0))
print("Vector part:", mixed.grade(1))

# Rotation
rotor = alg.rotor(e1 ^ e2, math.pi / 4)
rotated = v.rotate_by(rotor)
print("Rotated:", rotated)

# Norm and inverse
print("Norm:", v.norm())
print("Inverse:", v.inverse())

# Projection and rejection
proj = v.project_to(e1)
reject = v.reject_from(e1)
print("Projection:", proj)
print("Rejection:", reject)
```
