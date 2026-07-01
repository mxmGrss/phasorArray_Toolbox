# Contributing to PhasorArray Toolbox


## Reporting Issues

Open a GitHub issue with a minimal reproducible example. Include MATLAB version,
toolbox version (`ver` output), and the full error message.

## Pull Requests

PRs are welcome. Please ensure:
- New functions follow the existing style (`arguments` blocks, English comments).
- Mathematical claims are accompanied by a reference or a derivation sketch in the
  header comment.

## Development History and AI Assistance

**Up to November 2025**, all code in this toolbox was written manually by the authors.
The original implementation (approximately 70% of the codebase) reflects several years
of research work at CRAN.

**From November 2025 onwards**, the toolbox is maintained as a side project alongside
active postdoctoral research. New developments are AI co-authored (Claude, Anthropic)
under the following workflow:

- Mathematical derivations and algorithmic decisions are developed and reviewed by the
  authors before implementation.
- Every line of generated code is manually reviewed and validated against theory.
- Numerical correctness is verified by residual checks and, where applicable,
  cross-validation against independent implementations.
- Commits involving AI assistance carry the `Co-Authored-By: Claude` trailer.

This disclosure follows the spirit of emerging norms in open-source and academic
software projects regarding AI-assisted development.

