# Contributing to PhasorArray Toolbox

Contributions are welcome — this toolbox is maintained as a side project
alongside active research, so a little process keeps things manageable.

## Reporting Issues

Search existing issues first, then open a GitHub issue with a minimal
reproducible example. Include MATLAB version, toolbox version (`ver` output),
and the full error message.

## Pull Requests

1. For anything larger than a bug fix, open an issue first so we can discuss
   scope before you invest time.
2. Fork, branch from `main`, run `installToolbox.m`, make your changes.
3. Pre-validate: both test suites must pass without new failures:

       test_PhasorArray_basic();
       test_PhasorArray_advanced();

   As of today the suites use a custom lightweight runner (struct-based, not
   `matlab.unittest`); new tests should follow the same style. A migration to
   `matlab.unittest` is planned to enable CI on pull requests.
4. Keep commits atomic (one change per commit) and PRs focused (one feature
   or fix per PR).
5. Purely cosmetic changes (reformatting, renaming) are generally not
   accepted on their own.

New functions follow the existing style:
`arguments` blocks with validators, English comments, no global variables,
`warning('PhasorArray:...')` identifiers for diagnostics. Mathematical claims
are accompanied by a reference or a derivation sketch in the header comment.

## Licensing of Contributions

The toolbox is MIT-licensed. There is no Contributor License Agreement: by
submitting a pull request you agree that your contribution is licensed under
the same MIT license as the project ("inbound = outbound"). You retain the
copyright on your contribution.

## Development History and AI Assistance

**Up to November 2025**, all code in this toolbox was written manually by the
authors. The original implementation (approximately 70% of the codebase)
reflects several years of research work at CRAN.

**From November 2025 onwards**, the toolbox is maintained as a side project
alongside active postdoctoral research. New developments are AI co-authored
(Claude, Anthropic) under the following workflow:

- Mathematical derivations and algorithmic decisions are developed and
  reviewed by the authors before implementation.
- Every line of generated code is manually reviewed and validated against
  theory.
- Numerical correctness is verified by residual checks and, where applicable,
  cross-validation against independent implementations.
- Commits involving AI assistance carry the `Co-Authored-By: Claude` trailer.

This disclosure follows the spirit of emerging norms in open-source and
academic software projects regarding AI-assisted development.

Thanks!

Maxime Grosso (CRAN, Université de Lorraine)
