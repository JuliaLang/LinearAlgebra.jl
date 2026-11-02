# LinearAlgebra.jl tests

## Standard and comprehensive mode

The test suite runs in one of two modes:

- **Standard** mode is what `Pkg.test()`, Julia's own CI and all CI jobs except the
  coverage job run. It should test each feature once for correctness, typically over
  `Float64` and `ComplexF64`, plus one generic (non-BLAS) element type where the generic
  fallback is what is being tested. It must stay fast: Julia's CI runs all test files
  serially in one process, and almost all test time is compilation.
- **Comprehensive** mode additionally runs the tests guarded with
  `@static if COMPREHENSIVE`: more element types, wrappers, `uplo` variants, sizes,
  promotion pairs and corner cases. It is enabled with
  `Pkg.test("LinearAlgebra"; test_args=["--comprehensive"])` or by setting
  `LINALG_TEST_COMPREHENSIVE=true`. In CI it runs on the ubuntu x64 job, which also
  collects coverage, split into shards via `LINALG_TEST_SHARD` (e.g. `1/2`).

`COMPREHENSIVE` is exported from `LinearAlgebraTestHelpers` (`testhelpers/testhelpers.jl`).
Guard with `@static`, so that a standard run does not lower or compile the guarded code,
and put the guard on the values a test runs over rather than around a copy of its body:

```julia
@testset for elty in (Float64, ComplexF64, (@static COMPREHENSIVE ? (Float32, ComplexF32, BigFloat, Int) : ())...)
    ...
    @static if COMPREHENSIVE
        @test_throws DimensionMismatch ...
    end
end

@static if COMPREHENSIVE
@testset "issue #1234" begin ... end
end
```

New variations of existing tests and regression tests for reported issues belong in
comprehensive mode; an unguarded test should be one representative case of a feature.

## Sharding

`test/testgroups` lists the test files in order of decreasing runtime. With
`LINALG_TEST_SHARD=i/n`, `runtests.jl` runs only the `i`-th of `n` shards, dealing the
files out in snake order (1, 2, 2, 1, 1, 2, ...) to balance the shards. Keep the order
roughly up to date when adding test files or when runtimes change significantly.
