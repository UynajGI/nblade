# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.6] - 2026-06-14

### Features
- Added mypy type checking to CI
- Added comprehensive input validation with descriptive error messages
  (dimension bounds, signature consistency, basis index, coefficient length)

### Documentation
- Clarified Python operators vs math notation in README tables

## [0.1.5] - 2026-06-14

### Features
- Added LaTeX rendering support (`_repr_latex_`) for Jupyter notebook display
- Added Python property-based tests (22 hypothesis tests, 2D-5D)
- Added coverage tracking (cargo-llvm-cov + pytest-cov) to CI

### Bug Fixes
- Fixed CI: merged Python build+test into single venv step
- Fixed clippy `manual_strip` and ruff violations across Python code

### Documentation
- Documented `inverse()` limitation (pure-grade only) in API docs
- Documented `|` (left_inner) metric signature caveat in API docs
- Updated dimension limit docs (64D → practical ~12D)
- Split README into English-only and Chinese-only versions

## [0.1.3] - 2026-06-14

### Features
- Added `__eq__`, `__ne__`, `get_coefficient()` to Python MultiVector bindings

### Bug Fixes
- Fixed 5 Python test bugs (basis_count assertion, basis_vectors unpacking, commutator formula)
- Fixed `^` operator precedence bugs in examples (Python parses `e1 ^ e2 + 2*(e2 ^ e3)` as `e1 ^ (e2 + ...)`)
- Fixed dual-dual expectation in advanced examples (library convention: dual(dual(v)) = v)
- Fixed inverse example to use pure bivector (rev(A)/|A|^2 formula only valid for pure-grade elements)
- Fixed non-Euclidean norm bug: `(A|B).scalar_part()` ignores metric signature, replaced with `norm_squared()`
- Fixed Cramer's rule example in applications (`.grade(1)` → `.scalar_part()`)
- Fixed dual-dual roundtrip assertions in tutorial (dual(dual(v))=v, inverse_dual(dual(v))=-v in 3D)

### Documentation
- Added `get_coefficient`, `__eq__`/`__ne__` to MultiVector API docs (en + zh)
- Fixed `coefficients()` return type in docs (List[float] → numpy.ndarray)
- Renumbered example files to consistent 01-starting numbering
- Updated examples/README.md paths to match renumbered files
- Fixed i18n: replaced mkdocs-i18n with mkdocs-static-i18n, renamed zh-cn → zh

## [0.1.2] - 2026-03-23

### Features
- Added comprehensive examples: spacetime algebra, conformal GA, electromagnetism, rotor interpolation
- Added Chinese documentation (zh-cn)
- Added CONTRIBUTING.md with bilingual guidelines
- Added separate CI workflows for TestPyPI, PyPI, and Crates.io publishing

### Documentation
- Reorganized docs structure for i18n support
- Added zh-cn/index.md as Chinese landing page

### Bug Fixes
- Fixed ReadTheDocs build configuration
- Fixed mkdocs i18n language code (zh → zh-cn)

## [0.1.1] - 2026-03-22

### Features
- Initial release with core geometric algebra operations
- Python bindings via PyO3
- Support for arbitrary dimensions (up to 64D; practical limit ~12D for dense operations)
- Support for arbitrary metric signatures G(p, q, r)
- Dense and sparse multivector representations
- NumPy integration

[Unreleased]: https://github.com/UynajGI/nblade/compare/v0.1.6...HEAD
[0.1.6]: https://github.com/UynajGI/nblade/compare/v0.1.5...v0.1.6
[0.1.5]: https://github.com/UynajGI/nblade/compare/v0.1.4...v0.1.5
[0.1.3]: https://github.com/UynajGI/nblade/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/UynajGI/nblade/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/UynajGI/nblade/releases/tag/v0.1.1
