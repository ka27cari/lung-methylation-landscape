# Code release notes (2026-11-03)

This revision is a **directory/name refactor**, not a new analysis or revision of manuscript statistics.

- Flattened numbered Python stage folders into descriptive names under `pipeline/`.
- Retained the explicit stage order in `run.py` and corrected relative project-root paths for the moved scripts.
- Consolidated original R analyses in `r/analysis/` and figure code in `r/figures/`, with descriptive file names.
- Updated shared R setup references to the new locations; original R calculations were not intentionally changed.
- Fixed the external PBMC Python test import path, which originally depended on other tests importing the source package first.
- Added a repository-layout test. The source-data tooling and scientific data are unchanged.

**Checks in packaging environment:**
- Python synthetic and contract tests: 14 passing.
- Standalone external PBMC projection tests: 4 passing.
- `python run.py run --execute --allow-clean`: completed using authorized private staging and the saved recurrent-program discovery checkpoint (expensive de novo discovery was not rerun).
- Internal analytical output check: `status=pass`, `problems=[]`.
- `python tools/verify_science.py --root rebuild`: `Aggregate scientific contracts: PASS`.
- `python tools/build_source_data.py --root rebuild --out <outside-repo>`: 29 aggregate CSV tables exported.
- **R runtime was unavailable.** Renamed R source references were inspected, but R scripts were not executed or syntax-parsed with R.

No participant-level data or local `rebuild/` outputs belong in the release archive.
