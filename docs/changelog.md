# CHANGELOG

Notable changes to lanctools (starting with v1.0.0) will be documented here.

## [v1.0.3] - 2026-08-30

[v1.0.3]: https://github.com/frankp-0/lanctools/compare/v1.0.1...v1.0.3>

### Fixed

- Bug closing tracts on chromosome switches
- Use per-chromosome min/max positions to clip tracts

## [v1.0.1] - 2026-08-20

[v1.0.1]: https://github.com/frankp-0/lanctools/compare/v0.1.0...v1.0.1>

### Fixed

- Match breakpoints with .pvar file by CHR,POS not just POS

## [v1.0.0] - 2026-08-20

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
