# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0](https://github.com/UynajGI/nblade/compare/v0.1.1...v0.2.0) (2026-03-22)


### Features

* add release workflows, contributing guide, and reorganize docs ([21ede49](https://github.com/UynajGI/nblade/commit/21ede49d067efa64f0c022141ce3db23509bef11))
* **ci:** add release-please for automatic changelog generation ([2e36ba6](https://github.com/UynajGI/nblade/commit/2e36ba6c0aec39457f5a99da047b740781d69d98))


### Bug Fixes

* **ci:** add virtualenv for maturin develop + cargo fmt ([ec294aa](https://github.com/UynajGI/nblade/commit/ec294aadd5ff190fd013d429a0a27834d92f293b))
* **clippy:** add clippy allows for common patterns ([777c8f5](https://github.com/UynajGI/nblade/commit/777c8f5fc1018ca9a74a79199f7d1d5b5c228a82))
* **i18n:** use zh-cn language code for ReadTheDocs compatibility ([d4bff43](https://github.com/UynajGI/nblade/commit/d4bff43504d6de1d01b52c63e685dcde2b874cad))
* **mkdocs:** fix i18n plugin config format for ReadTheDocs ([3c2c0a4](https://github.com/UynajGI/nblade/commit/3c2c0a4682a87a70bfeee49383883ee476f3ae1c))
* **mkdocs:** fix i18n plugin languages config format ([7d8eceb](https://github.com/UynajGI/nblade/commit/7d8eceb3b9a49eaae2e5e75083aff7b6c66998b3))

## [Unreleased]

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
