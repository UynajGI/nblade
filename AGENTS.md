# NBLADE PROJECT KNOWLEDGE BASE

**Generated:** 2026-03-22T12:00:00+08:00
**Stack:** Rust 2021 Edition + Python 3.8+ (PyO3/maturin)

## OVERVIEW

nblade (N-dimensional Blade) - High-performance Geometric Algebra library for arbitrary dimensions (up to 64D). Supports any metric signature G(p,q,r) with dense/sparse multivector representations and parallel computation via rayon.

## STRUCTURE

```
nblade/
├── src/
│   ├── lib.rs              # Entry point, public API, prelude
│   ├── signature.rs        # Metric signature (p,q,r)
│   ├── algebra_config.rs   # Algebra configuration (Arc-shared)
│   ├── multiplication_table.rs  # Precomputed tables (rayon parallel)
│   ├── basis/              # Bitwise basis index operations
│   ├── multivector/        # Dense/Sparse multivector types (CORE)
│   ├── products/           # Geometric, outer, inner products
│   ├── operations/         # Involution, dual, inverse, norm
│   ├── geometry/           # Projection, rotation, reflection
│   └── python/             # PyO3 bindings (feature-gated)
├── tests/                  # Integration tests
├── benches/                # Criterion benchmarks
├── python/                 # Python package (nblade)
├── examples/               # Tutorial and application examples
│   ├── tutorials/          # Basic tutorials
│   ├── physics/            # Physics applications
│   └── cg/                 # Computer graphics examples
└── docs/                   # MkDocs documentation
    ├── zh/                 # Chinese documentation
    └── en/                 # English documentation
```

## WHERE TO LOOK

| Task | Location | Notes |
|------|----------|-------|
| Create multivector | `src/multivector/mod.rs` | `MultiVector::basis_vector()`, `from_coefficients()` |
| Geometric product | `src/products/geometric.rs` | Uses precomputed tables |
| Basis operations | `src/basis/index.rs` | Bitwise encoding, sign calculations |
| Metric signatures | `src/signature.rs` | Euclidean, spacetime, custom |
| Thread sharing | `src/algebra_config.rs` | `AlgebraConfigRef = Arc<AlgebraConfig>` |
| Python bindings | `src/python/mod.rs` | Feature `python` required |
| Tests | `tests/test_arbitrary_dimension.rs` | 19 tests, 1D-16D coverage |
| Benchmarks | `benches/geometric_product.rs` | Criterion, 3D and 5D |
| Volume element | `src/algebra_config.rs` | `volume_element()`, `volume_element_squared()`, `volume_element_inverse()` |
| Reciprocal frame | `src/basis/frame.rs` | `reciprocal_frame()`, `verify_reciprocal_frame()`, `metric_tensor()` |
| Basis expansion | `src/basis/frame.rs` | `basis_expansion()`, `basis_reconstruction()` |
| Dual operations | `src/operations/dual.rs` | `dual()`, `inverse_dual()` |

## CODE MAP

| Symbol | Type | Location | Role |
|--------|------|----------|------|
| `MultiVector` | Enum | `src/multivector/mod.rs` | Core type, 40+ operations |
| `DenseMultiVector` | Struct | `src/multivector/dense.rs` | ndarray-backed, all 2^n coeffs |
| `SparseMultiVector` | Struct | `src/multivector/sparse.rs` | HashMap-backed, <10% density |
| `Signature` | Struct | `src/signature.rs` | Metric (p,q,r), basis squares |
| `AlgebraConfig` | Struct | `src/algebra_config.rs` | Dimension, signature, metric factors |
| `BasisIndex` | Type | `src/basis/index.rs` | `u64` for up to 64 dimensions |
| `MultiplicationTables` | Struct | `src/multiplication_table.rs` | Precomputed lookup tables |
| `reciprocal_frame` | Function | `src/basis/frame.rs` | Compute reciprocal frame vectors |
| `verify_reciprocal_frame` | Function | `src/basis/frame.rs` | Verify aⁱ⌋aⱼ = δⁱⱼ |
| `metric_tensor` | Function | `src/basis/frame.rs` | Compute metric tensor g_ij |
| `basis_expansion` | Function | `src/basis/frame.rs` | Expand multivector in basis |
| `basis_reconstruction` | Function | `src/basis/frame.rs` | Reconstruct from coefficients |

## CONVENTIONS

**Basis Index Encoding:**
```rust
// Bitwise: bit i set = includes e_i
// 0b00101 = e1∧e3, 0b111 = e1∧e2∧e3
let index: BasisIndex = 0b101;  // e1∧e3
```

**Auto Dense/Sparse Selection:**
```rust
// Density < 10% → Sparse, else Dense
const SPARSE_THRESHOLD: f64 = 0.1;
```

**Zero Tolerance:**
```rust
const EPSILON: f64 = 1e-15;  // Used throughout for is_zero(), filtering
```

**Test Config Pattern (repeated in 13 files):**
```rust
fn test_config() -> AlgebraConfigRef {
    Arc::new(AlgebraConfig::euclidean(3))
}
```

**Float Assertions:**
```rust
assert!((value - expected).abs() < 1e-10);  // Never assert_eq! for floats
```

## ANTI-PATTERNS (THIS PROJECT)

- **No traits** — Uses concrete types and enum dispatch (MultiVector)
- **Chinese comments** — Docstrings in Chinese throughout
- **Manual Send/Sync** — `unsafe impl Send/Sync for AlgebraConfig` (lines 144-145 in algebra_config.rs)
- **No CI/CD** — Zero workflow files, Makefile, or justfile

## UNIQUE STYLES

**Module Re-exports:**
```rust
// basis/mod.rs - wildcard re-export pattern
pub use index::*;
pub use grade::*;
```

**Prelude Module:**
```rust
pub mod prelude {
    pub use crate::Signature;
    pub use crate::AlgebraConfig;
    pub use crate::MultiVector;
    pub use crate::basis::grade_of_index;
    pub use crate::basis::index::index_to_string;
}
```

**Dual Crate Type:**
```toml
# Cargo.toml
crate-type = ["cdylib", "rlib"]  # Python extension + Rust library
```

## COMMANDS

```bash
# Build
cargo build [--release]

# Test
cargo test --all-features

# Benchmark
cargo bench

# Python wheel (requires maturin)
maturin develop
maturin build --release
```

## NOTES

- **Max dimension: 64** (BasisIndex = u64)
- **No git history** — Fresh repo or squashed
- **Empty `src/utils/`** directory exists but unused
- **Python feature default ON** — `default = ["python"]`
- **Release profile:** LTO enabled, single codegen unit for optimization