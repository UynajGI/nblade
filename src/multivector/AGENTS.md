# MULTIVECTOR MODULE

**Parent:** [../AGENTS.md](../AGENTS.md)

Core data representation for Geometric Algebra. Auto-selects between dense (ndarray) and sparse (HashMap) storage based on coefficient density.

## FILES

```
multivector/
├── mod.rs      # MultiVector enum + 40+ operation dispatches
├── dense.rs    # DenseMultiVector (ndarray::Array1<f64>)
└── sparse.rs   # SparseMultiVector (HashMap<BasisIndex, f64>)
```

## API SURFACE

**Construction:**
```rust
MultiVector::zeros(config)           // Zero multivector
MultiVector::one(config)             // Scalar 1
MultiVector::from_scalar(config, s)  // Scalar value
MultiVector::basis_vector(config, i) // e_i (0-indexed)
MultiVector::from_coefficients(config, vec)  // Auto-selects dense/sparse
```

**Products (4 variants each):**
- `geometric_product(&other)` — Clifford product
- `outer_product(&other)` — Wedge product (⋀)
- `left_inner(&other)` — Contraction (A⌋B)
- `right_inner(&other)` — Right contraction (A⌊B)

**Operations:**
- `grade_involution()`, `reversion()`, `clifford_conjugate()`
- `dual()`, `inverse()`
- `norm()`, `norm_squared()`, `scalar_product(&other)`
- `commutator(&other)` — [A,B] = (AB-BA)/2

**Geometry:**
- `project_to(&blade)`, `reject_from(&blade)`
- `reflect_in(&blade)`, `rotate_by(&rotor)`

## INTERNAL CONVENTIONS

**Dense/Sparse Threshold:**
```rust
// sparse.rs
const SPARSE_THRESHOLD: f64 = 0.1;  // <10% non-zero → sparse
```

**Dispatch Pattern (every operation):**
```rust
match (self, other) {
    (Dense(a), Dense(b)) => Dense(a.op(b)),
    (Sparse(a), Sparse(b)) => Sparse(a.op(b)),
    (Dense(a), Sparse(b)) => Dense(a.op(&b.to_dense())),
    (Sparse(a), Dense(b)) => Dense(a.to_dense().op(b)),
}
```

**Zero Check:**
```rust
const EPSILON: f64 = 1e-15;
coef.abs() < EPSILON  // Considered zero
```

## WHERE TO LOOK

| Task | File | Function |
|------|------|----------|
| Add operation | `dense.rs` / `sparse.rs` | `add()`, `sub()` |
| Geometric product impl | `dense.rs` | Uses `MultiplicationTables` |
| Grade projection | `dense.rs:111` | `grade_projection(r)` |
| Display formatting | `dense.rs:202`, `sparse.rs` | `impl Display` |
| Sparse → Dense conversion | `sparse.rs` | `to_dense()` |

## NOTES

- **No traits** — All operations on concrete types
- **Arc-shared config** — `config: AlgebraConfigRef` field in both variants
- **Mixed ops → Dense** — Always converts to dense for mixed dense/sparse operations