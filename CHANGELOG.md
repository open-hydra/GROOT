# Changelog

All notable changes to GROOT are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). An entry that begins with **Breaking.** asks a user to
change a build.

## [Unreleased]

### Added

- Add workflow to automayically source contributors list.
- `install.sh build --update-deps` selects the newest ORION 2.x.
- [Versioning](docs/development/versioning.md) page.
- The configure stops when ORION or FiNeR is not the release GROOT is tested with (`GROOT_PIN_CHECK=WARNING` to build it on purpose).

### Changed

- **Breaking.** `ORION_PATH` and `FINER_PATH` are removed: use `FETCHCONTENT_SOURCE_DIR_ORION` and
  `FETCHCONTENT_SOURCE_DIR_FINER` (`install.sh --include-*` as before). `install.sh update` is removed.
- **Breaking.** FiNeR (v2.3.0) moves to `lib/FiNeR`; delete an old `lib/third_party/FiNeR`. Its five libraries are no
  longer submodules of GROOT: they follow their default branches (`GROOT_FINER_DEPS=BACKUP` for fixed releases).
- `install.sh` is executable.

## [0.1.0] - 2026-05-21

First public release. See the [README](README.md) for the features.

[Unreleased]: https://github.com/open-hydra/GROOT/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/open-hydra/GROOT/releases/tag/v0.1.0
