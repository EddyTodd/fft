## Scope

Describe the FFT mechanism, planner/API change, kernel work, correctness change, documentation update, or migration work in this PR.

## Library contract

- [ ] Public algorithm/plan/planner/kernel API changes are intentional and documented.
- [ ] Forward/inverse normalization and workspace/scratch semantics are unchanged or explicitly described.
- [ ] Planner selection/capability behavior is unchanged or explicitly described.
- [ ] Architecture-specific kernels do not silently substitute unsupported implementations under the same identity.
- [ ] Numerical claims identify the error norm/reference and whether they are guarantees or regression contracts.

## Correctness

- [ ] Deterministic coverage includes the transform-size families affected by the change (power-of-two/composite/coprime/prime/arbitrary as applicable).
- [ ] Reusable-plan and round-trip behavior was checked when planning/workspace code changed.
- [ ] Changed ISA kernels were exercised on each architecture claimed below, or are explicitly marked unvalidated.
- [ ] All installed public headers remain independently includable and the external package consumer remains valid.

## Research boundary

- [ ] Generic benchmark timing/statistics/campaign/provenance/reporting code remains in `EddyTodd/bench`.
- [ ] Any deletion of vendor/kernel/arbitrary/real-plan research assets has executed treatment/correctness/evidence parity.

## Validation performed

List exact commands, compiler/platform/ISA, transform sizes/families, and results. Distinguish compilation from real-target execution.

## Release impact

- [ ] No release-note entry required.
- [ ] `CHANGELOG.md` updated for a release-impacting change.
- [ ] Version/citation changes, if any, follow `docs/release-process.md`.
