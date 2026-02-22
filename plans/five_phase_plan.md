# VLE Final Project Five-Phase Plan

## Phase 1 - Baseline Lock and Data Governance

Goal:
- Freeze the upstream fork boundary and keep upstream code unchanged.
- Make all project parameters traceable and auditable.

Implementation:
- Use `vle_consulting/data/parameter_db.json` as the single parameter source.
- Require source metadata for every parameter block.
- Enforce explicit `tau = a + b/T` terms for NRTL.

Exit criteria:
- Fork area has no new edits.
- Parameter loader rejects missing-source datasets.

## Phase 2 - Core Thermodynamic Engine

Goal:
- Build stable and testable calculation modules.

Implementation:
- `src/models.py`: Antoine, Van Laar, NRTL models.
- `src/solver.py`: bubble-point solver and T-x-y generation.
- `src/analysis.py`: GE/HE and azeotrope detection.

Exit criteria:
- Pipeline can run from CLI without report-generation code.
- Results are reproducible from database values only.

## Phase 3 - Interface and Output Standardization

Goal:
- Standardize naming and execution flow.

Implementation:
- Keep simple entrypoint `src/main.py`.
- Keep orchestration in `src/pipeline.py`.
- Keep output naming simple and lowercase.

Exit criteria:
- One command runs full analysis.
- Output files are deterministic and easy to parse.

## Phase 4 - Verification and Quality Gates

Goal:
- Prevent silent regressions.

Implementation:
- `tests/test_parameter_store.py` validates data-source and NRTL constraints.
- `tests/test_pipeline.py` validates runtime shape and numerical sanity.

Exit criteria:
- Unit tests pass in local environment.
- Basic numerical checks remain stable across reruns.

## Phase 5 - Expansion for Final Submission

Goal:
- Extend from single-system analysis to submission-ready engineering package.

Implementation backlog:
- Add uncertainty/sensitivity study (pressure, parameter perturbation).
- Add parameter-set versioning for alternative literature sources.
- Add numeric result export (`csv`/`json`) for appendix tables.
- Add experiment-to-model comparison hooks if new lab data are available.
- Add CI job for `vle_consulting` standalone test run.

Exit criteria:
- Final package contains reproducible computation, traceable data, and validation logs.
