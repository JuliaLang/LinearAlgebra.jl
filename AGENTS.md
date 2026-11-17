# AGENTS.md

Guidance for agents and other contributors to LinearAlgebra.jl, distilled from the
project's recent pull requests and reviews. Follow it unless a maintainer says otherwise.

## What this repo is

The Julia standard library for dense and structured linear algebra: the BLAS and LAPACK
wrappers (`src/blas.jl`, `src/lapack.jl`, loaded through libblastrampoline), the
structured matrix types (`Diagonal`, `Bidiagonal`, `Tridiagonal`, `SymTridiagonal`,
triangular, `Symmetric`/`Hermitian`, `UpperHessenberg`), the factorizations and their
`AbstractQ` types, `UniformScaling`, and generic fallbacks for all of it.

Being a stdlib, and one that is baked into Julia's system image, shapes everything below:

- Develop and test on Julia nightly. Each Julia release ships a fixed copy of this
  package, so there are no `@static if VERSION` branches; to reach an older Julia,
  backport instead (see Backports).
- Code reaches users only when JuliaLang/julia bumps its pinned LinearAlgebra version
  (`stdlib/LinearAlgebra.version`).
- The sysimage already contains a LinearAlgebra, and this package pirates Base methods.
  A plain `julia --project` silently runs the **bundled** copy, not this checkout. Read
  "Setting up and running tests" before running anything.
- SparseArrays is pinned in lockstep and builds on internal hooks of this package.
  Changes to those hooks need a paired PR there.
- Base runs this package's tests in its own CI, serially, without ParallelTestRunner.
  Keep `test/runtests.jl`'s serial fallback, and keep test time down.
- No new dependencies, neither for the package nor for the tests. Test-only types live
  in `test/testhelpers/`; if none fits, do without. Do not mention external packages in
  code or docstrings.
- Julia's manual includes `docs/src/index.md` and the docstrings it lists. Every exported
  or public name needs a docstring listed there (`runtests.jl` checks
  `Docs.undocumented_names`).
- Every source and test file starts with the Julia MIT license banner.
- The default branch is `master`. Release branches are `release-1.x`; backports collect
  on `backports-release-1.x`.

## Setting up and running tests

The environment needs only Julia nightly (`julia +nightly` with juliaup) and the
package's test dependencies (`ParallelTestRunner`, which `Pkg.test` installs). Two
things are required for every run, or the bundled LinearAlgebra is tested instead of
your changes:

1. `JULIA_PRUNE_OLD_LA=true` in the environment, together with `test/prune_old_LA.jl`
   being included before `using LinearAlgebra`. It deletes the sysimage copy's methods so
   that the checkout can replace it. `runtests.jl`, every test file and `docs/make.jl`
   include it themselves; a script of your own must include it first.
2. `--compiled-modules=existing`. Without it you get "Failed to precompile
   LinearAlgebra", "Replacing module `LinearAlgebra`", or stack traces pointing into the
   stdlib copy.

Confirm in a new script that `pathof(LinearAlgebra)` points into the checkout. Run
everything from the repository root:

```sh
# one or more test files, through ParallelTestRunner as CI does (names from test/testgroups)
JULIA_PRUNE_OLD_LA=true JULIA_NUM_THREADS=1 julia +nightly --project --compiled-modules=existing -e 'using Pkg; Pkg.test(julia_args=`--compiled-modules=existing`, test_args=`givens triangular`)'

# the whole suite; add `--comprehensive` to test_args for comprehensive mode
JULIA_PRUNE_OLD_LA=true JULIA_NUM_THREADS=1 julia +nightly --project --compiled-modules=existing -e 'using Pkg; Pkg.test(julia_args=`--compiled-modules=existing`)'

# one file directly, without Pkg.test: the fastest edit-test loop
JULIA_PRUNE_OLD_LA=true julia +nightly --project --compiled-modules=existing -e 'include("test/givens.jl")'
JULIA_PRUNE_OLD_LA=true LINALG_TEST_COMPREHENSIVE=true julia +nightly --project --compiled-modules=existing -e 'include("test/givens.jl")'

# checks CI runs besides the tests
julia .ci/check-whitespace.jl
JULIA_PRUNE_OLD_LA=true julia +nightly --project --compiled-modules=existing test/ambiguous_exec.jl   # after changing signatures
JULIA_PRUNE_OLD_LA=true julia +nightly --project=docs --compiled-modules=existing -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'   # docs and doctests
```

Notes:

- A single file takes about 20 s directly and about 90 s through `Pkg.test`; almost all
  of it is compilation. The full comprehensive suite takes 20–30 minutes, so run long
  jobs in the background with a generous timeout. Output written to a file is buffered,
  so an empty log does not mean that nothing is happening.
- A top-level `@testset` that fails aborts the rest of an `include`d file. To see every
  result, wrap the code you are checking in one outer testset.
- If the package server is unreachable, `JULIA_PKG_SERVER=""` makes Pkg fall back to
  git.
- For coverage of a change, run the relevant file with
  `--code-coverage=@$(pwd)/src` and inspect the `.cov` files next to the sources (delete
  them afterwards).
- Nightly sometimes breaks tests for reasons unrelated to this package. Before blaming a
  change for a failure, run the same test on `master`.

### Checking that a test fails on master

Make a worktree of `master` once, and run the same probe script against both checkouts
by switching `--project`:

```sh
git fetch origin master && git worktree add ../LA-master origin/master   # once; later: git -C ../LA-master checkout --detach origin/master
```

```julia
# probe.jl: the new testset(s), copied verbatim
include(joinpath(ENV["LA_DIR"], "test", "prune_old_LA.jl"))
using Test, LinearAlgebra
@testset "probe" begin
    # ...
end
```

```sh
LA_DIR=$PWD JULIA_PRUNE_OLD_LA=true julia +nightly --project=. --compiled-modules=existing probe.jl
LA_DIR=../LA-master JULIA_PRUNE_OLD_LA=true julia +nightly --project=../LA-master --compiled-modules=existing probe.jl
```

If the probe needs a test helper, `include` it from `test/testhelpers/` after
`prune_old_LA.jl`. Do not load `Random` or other stdlibs that depend on LinearAlgebra
before `prune_old_LA.jl`.

## Layout

- Source files in `src/` are named by type or area (`diagonal.jl`, `bidiag.jl`,
  `tridiag.jl`, `triangular.jl`, `symmetric.jl`, `hessenberg.jl`, one file per
  factorization, `matmul.jl`, `dense.jl` for dense matrix functions, `generic.jl` for
  `AbstractArray` fallbacks, `special.jl` for operations between different structured
  types, `structuredbroadcast.jl`). Test files mirror them; the triangular tests are
  split over `triangular*.jl` with shared helpers in `testtriag.jl`.
- Every test file is a module that includes `prune_old_LA.jl` and, when it needs them,
  `testhelpers/testhelpers.jl` (`COMPREHENSIVE` and the test-only types: `Furlongs`
  for dimensionful numbers, `Quaternions` for a non-commutative eltype, `DualNumbers`,
  `OffsetArrays`, `SizedArrays`, `ImmutableArrays`, `StructArrays`, `FillArrays`,
  `InfiniteArrays`, `StridedArrays`). Add a type there, not to a test file.
- `test/testgroups` lists the test files by decreasing runtime; `runtests.jl` reads it.
  A new test file must be added there.
- `test/README.md` describes standard and comprehensive mode and sharding.

## Coding conventions

- **Stay generic.** Code reachable by generic element types must not call BLAS or
  LAPACK directly; specialize for `BlasFloat` by dispatch and keep a generic fallback
  for everything else. Aim for maximum element-type independence: `BigFloat`,
  `Rational{BigInt}`, quaternions, dual numbers, and dimensionful numbers (units).
  Do not assume `T * T` or `one(T)` is of type `T`; infer result eltypes the way the
  existing `promote_op`-based machinery does, and check that the result is inferred
  with a concrete eltype. A `BlasFloat` restriction that is a leftover from when only
  the LAPACK path existed should go.
- **Use dispatch, not type branches** in a method body (`T <: X ? ... : ...`).
  Specializing for concrete types is good when it avoids runtime work.
- **Do not throw on working cases.** When a fast path (BLAS strides, LAPACK
  requirements) does not apply, fall back to the generic kernel instead of erroring. Inf
  and NaN propagate silently; only add checks that prevent hangs, crashes or wrong
  results, and say so in the PR when a change rejects input that used to be accepted.
  Change behaviour in a separate step from fixing a bug.
- **Follow dense semantics** for structured matrices: the same return values (including
  NaN/Inf propagation) as `Matrix(A)`, with the structured result type where one fits.
- **Allocations.** Do not add allocations to in-place (`!`) paths; reuse memory where
  possible.
- **Specialize on functions.** Write `f::F ... where {F}` rather than `f::Function` for a
  function argument that is only passed on, so the method specializes (and trims).
- **Narrowing signatures** risks ambiguities with downstream packages (SparseArrays,
  GPUArrays and the GPU vendor packages); check them, and run the ambiguity test.
- **Keep code minimal.** Prefer one-line method definitions where they fit, forward
  keyword arguments instead of repeating defaults, inline needless indirections,
  factor repeated checks into a helper, and reuse the existing internal helpers
  (e.g. `_sym_uplo`, `_valeltype`, `_tohalf`) instead of reimplementing them. A
  refactor moves code; it does not rewrite it.
- **Error messages** are `LazyString`s built from comma-separated parts, with
  interpolation done lazily. `ArgumentError` for bad values, `DimensionMismatch` for
  shapes.
- **Comments** only for genuinely tricky code, and they say why. History belongs in the
  PR, not in the code. Internal helpers may have short docstrings.
- When adding a `@static` guard, do not reindent the guarded code; keep the diff small.
- Style-only changes go in their own PR.

## Tests

Coverage of this package is near perfect, so a new test is the exception, not the
default. The rules:

- **Test only changed behaviour.** Add a test only for what the PR changes: a fixed bug,
  a new method, a new branch. Do not test functionality that existing tests already
  cover, and do not add an assertion that an existing loop already makes. If the
  changed code is already exercised (check with coverage, or by planting an error in
  the branch and running the existing tests), the PR may need no new test at all;
  say which existing tests cover it.
- **New tests should fail on the previous master.** Check every new assertion against
  `master` (see above) and drop the ones that pass there, unless one is needed to cover
  a new branch. Keep only the minimal reproducer.
- **Do reach every new branch.** A test that only goes through the BLAS/LAPACK route is
  uninteresting when the change is in the generic path; use an element type or wrapper
  that actually reaches the code that changed.
- **Standard vs comprehensive mode** (`test/README.md`): standard mode tests each
  feature once, over `Float64` and `ComplexF64` plus one generic eltype where the
  generic path is the point. Variations (more eltypes, wrappers, `uplo`, sizes,
  promotion pairs) go under `@static if COMPREHENSIVE`, with the guard on the values a
  test runs over, not around a copy of its body. A guard is worth it only when it saves
  real compile time; a single quick check needs none.
- **Make tests able to fail.** Use inputs that expose the bug: a non-symmetric parent
  for `Symmetric`/`Hermitian` (so that reading the wrong triangle shows), both `:U` and
  `:L`, Hermitian data with nonzero imaginary diagonal parts, complex values that differ
  from their conjugates, a destination filled with nonzeros (`ones`) so that a missing
  write shows. Use the public interface (`F.R0`), not internal fields.
- **Use an independent reference.** Compare against `Matrix(A)` or another code path
  that does not run through the code under test; `lu(Matrix(A))` or a LAPACK result for
  a generic-eltype computation.
- **Inputs are fixed or seeded.** Julia's CI fails the stdlib bump on a test that fails
  for some draws. Prefer literals and constructions with a guaranteed property.
  `Rational{Int}` overflows easily; use `BigInt`.
- **Write plain tests.** Plain testsets with meaningful names rather than deep nesting;
  a `for (A, B, ...) in (...)` loop computing in the body rather than closures; hoist
  repeated conversions (`MA = Matrix(A)`). Use `==` rather than `isequal` unless NaN is
  involved. Do not qualify exported names with `LinearAlgebra.`.
- No allocation tests, except where zero allocations are expected. No wall-clock
  assertions. No tests of docstrings.
- Name issue regressions `"issue #NNNN"` and keep existing issue references.

## Pull requests and GitHub

- Report before acting outward: do not open a PR, push to a PR branch, or post a
  comment, review or reply on GitHub unless asked. When asked to look at a PR or issue,
  report back rather than posting. A review, when asked for, says it was written by an
  agent.
- Work on a branch off current `master`. Maintainers often push to PR branches
  themselves (fixes, conflict resolutions, simplifications), so fetch the branch before
  you push to it, and build on top of it.
- One logical change per PR; split anything separable (an unrelated fix, an
  optimization, a cleanup) into its own PR, stacked if needed. Do not carry commits from
  another PR. Within a PR, order commits so the fix applies cleanly to `master` on its
  own: fix, then tests, then cleanup. When asked to amend, fold fixes into the earlier
  commits instead of adding new ones.
- PRs are squash-merged with the title as the subject. Titles are imperative and name
  the API involved (see `git log --oneline`).
- Descriptions are short: one or two sentences, `Fixes #NNNN` where it applies, and for
  a behaviour change the current and the fixed behaviour (a reproducer with output).
  Mention deliberate behaviour changes (e.g. newly rejected input) and whether a
  backport is needed. Do not report what CI checks (tests run, pass counts, whitespace)
  or where the tests went; the diff shows that. Keep the description in sync when the
  PR changes.
- Credit others: when a PR builds on someone's code or ports someone's PR, add a
  `Co-authored-by:` trailer for them.
- The repository has no PR template.

## Backports and release notes

- Maintained releases are 1.10 (LTS) and the most recent ones; the maintainer adds the
  `backport 1.x` labels. When asked whether a fix needs a backport, check the code on
  `origin/release-1.x` and whether the fix applies there.
- Backports are bug and regression fixes only, never new methods, behaviour changes or
  master-only refactors. A backport PR targets `backports-release-1.x`, and its tests
  must not use `COMPREHENSIVE` (it only exists from 1.14 on); spell out the types.
  Test it with that release's Julia.
- For 1.10, LinearAlgebra still lived in JuliaLang/julia under `stdlib/`, so backports go
  to that repository's `backports-release-1.10`. Reference this repository's PRs there
  as `JuliaLang/LinearAlgebra.jl#NNNN`; a bare `#NNNN` links to julia's.
- `HISTORY.md` records new or newly public API, user-visible behaviour changes,
  deprecations and notable speedups per release, not bug fixes, tests or refactors.
  Follow the existing format: one bullet per change, PR link at the end.
