# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Features
- Added `__eq__`, `__ne__`, `get_coefficient()` to Python MultiVector bindings

### Bug Fixes
- Fixed 5 Python test bugs (basis_count assertion, basis_vectors unpacking, commutator formula)
- Fixed `^` operator precedence bugs in examples (Python parses `e1 ^ e2 + 2*(e2 ^ e3)` as `e1 ^ (e2 + ...)`)
- Fixed dual-dual expectation in advanced example (library convention: dual(dual(v)) = v)
- Fixed inverse example to use pure bivector (rev(A)/|A|^2 formula only valid for pure-grade elements)

### Documentation
- Added `get_coefficient`, `__eq__`/`__ne__` to MultiVector API docs (en + zh-cn)
- Fixed `coefficients()` return type in docs (List[float] → numpy.ndarray)

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
- Support for arbitrary dimensions (up to 64D)
- Support for arbitrary metric signatures G(p, q, r)
- Dense and sparse multivector representations
- NumPy integration

[Unreleased]: https://github.com/UynajGI/nblade/compare/v0.1.2...HEAD
[0.1.2]: https://github.com/UynajGI/nblade/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/UynajGI/nblade/releases/tag/v0.1.1