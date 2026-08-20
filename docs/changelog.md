# CHANGELOG

Notable changes to lanctools (starting with v1.0.0) will be documented here.

## [v1.0.0] - 2026-08-20

[v1.0.0]: https://github.com/frankp-0/lanctools/compare/v0.8.0...v1.0.0>

### Added

- LancData now supports close() and context-manager usage.
- Added validation for public query indices and FlatLanc structure.
- Added robust validation for FLARE/RFMix conversion inputs.
- Added validation and safer handling for .lanc file merging.

### Changed

- Invalid indices, malformed .lanc files, mismatched PLINK metadata, samples, chromosomes, or variants now raise explicit errors.
- Native FLARE/RFMix parser errors are surfaced clearly.
- get_lanc_dosage typing now reflects its actual int32 result.
- Documentation now demonstrates LancData resource management.
