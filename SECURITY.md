# Security policy

## Supported versions

Until multiple maintained release lines exist, security fixes target current `main` and the latest stable `1.x` line. Older snapshots are not maintained independently.

## Reporting a vulnerability

Do not publish exploit details, proof-of-concept inputs, or sensitive crash traces in a public issue before coordinated disclosure.

Use GitHub private vulnerability reporting if the repository exposes it. If no private channel is available, open a minimal public issue requesting a private security contact channel and omit exploit details until one is established.

Useful reports include:

- affected commit/version, architecture, OS, compiler, and transform size;
- minimal call sequence or triggering data shape;
- sanitizer/crash diagnostics;
- affected algorithm/plan/kernel if known;
- concrete impact and attacker-controlled prerequisites;
- proposed mitigation if available.

## Security-relevant scope

Examples include memory-safety defects, workspace/scratch-size mistakes, integer overflow in size/index calculations, undefined behavior reachable through public APIs, unsupported-ISA execution, malformed-size inputs causing unintended resource exhaustion, and correctness defects with a concrete downstream security consequence.

Numerical-accuracy disagreements and performance regressions normally follow the ordinary bug/research process unless they create a specific security impact.

## Disclosure

Please allow time to reproduce, fix, validate representative transform families and affected architectures, and prepare a release before public disclosure. Release notes must state validation limitations explicitly.
