# Changelog

All notable changes to GROOT are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). An entry that begins with **Breaking.** asks a user to
change a build.

## [Unreleased]

## [1.0.0] - 2026-11-14

### Added

- `vtk_io` test: VTK gas field read, solution and wall files written and read back, ascii and raw (run in CI).
- A workflow that generates the contributors list (`authors.md`) from the git history.
- `install.sh build --update-deps` selects the newest ORION 2.x.
- [Versioning](docs/development/versioning.md) page.
- The configure stops when ORION or FiNeR is not the release GROOT is tested with (`GROOT_PIN_CHECK=WARNING` to build
  it on purpose).
- The documentation is built with Zensical (strict in CI), with a new home page.

### Changed

- **Breaking.** `ORION_PATH` and `FINER_PATH` are removed: use `FETCHCONTENT_SOURCE_DIR_ORION` and
  `FETCHCONTENT_SOURCE_DIR_FINER` (`install.sh --include-*` as before). `install.sh update` is removed.
- **Breaking.** FiNeR (v2.3.0) moves to `lib/FiNeR`; delete an old `lib/third_party/FiNeR`. Its five libraries are no
  longer submodules of GROOT: they follow their default branches (`GROOT_FINER_DEPS=BACKUP` for fixed releases).
- `vtk raw` output is raw binary, as the wall file; it was base64 in the solution file.
- The first reference point of `source_j25` in the furnace test is not compared (see its V&V page).
- `install.sh` is executable.

### Fixed

- `vtk ascii` and `vtk raw` input: the field is read from `OUTPUT/field.vtm`, with its time.
- VTK output: variable names without quotes (the block files were not valid XML), and the wall file no longer crashes.
- GROOT builds as a part of a larger project (hydra): its sources and `data/snb` are found from GROOT's own folder.

## [0.1.0] - 2026-05-21

First public release. See the [README](README.md) for the features.

[Unreleased]: https://github.com/open-hydra/GROOT/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/open-hydra/GROOT/compare/v0.1.0...v1.0.0
[0.1.0]: https://github.com/open-hydra/GROOT/releases/tag/v0.1.0
