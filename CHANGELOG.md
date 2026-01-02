LinearAlgebra v1.12 Release Notes
=======================

- `rank` can now take a `QRPivoted` matrix to allow rank estimation via QR factorization ([#54283](https://github.com/JuliaLang/julia/issues/54283))
- Added keyword argument `alg` to `eigen`, `eigen!`, `eigvals` and `eigvals!` for self-adjoint matrix types (`RealHermSymComplexHerm`) to switch between different eigendecomposition algorithms ([#49355](https://github.com/JuliaLang/julia/pull/49355))
- Added a generic version of the (unblocked) pivoted Cholesky decomposition (callable via `cholesky[!](A, RowMaximum())`) ([#54619](https://github.com/JuliaLang/julia/pull/54619))
- The number of default BLAS threads now respects process affinity, instead of using the total number of logical threads available on the system ([#55574](https://github.com/JuliaLang/julia/pull/55574))
- Added a new function `zeroslike` that generates zero elements for matrix-valued banded matrices. Custom array types may specialize this function ([#55252](https://github.com/JuliaLang/julia/pull/55252))
- `A * B` matrix multiplication now calls `matprod_dest(A, B, T::Type)` to generate the destination; this function is now public ([#55537](https://github.com/JuliaLang/julia/pull/55537))
- The function `haszero(T::Type)` to check if a type has a unique zero element is now public ([#56223](https://github.com/JuliaLang/julia/pull/56223))
- Added `diagview` to return a view into a specific band of an `AbstractMatrix` ([#56175](https://github.com/JuliaLang/julia/pull/56175))

LinearAlgebra v1.11 Release Notes
=======================

- `cbrt(::AbstractMatrix{<:Real})` is now defined and returns real-valued matrix cube roots of real-valued matrices ([#50661](https://github.com/JuliaLang/julia/pull/50661))
- `eigvals`/`eigen(A, bunchkaufman(B))` and `eigvals`/`eigen(A, lu(B))`, which utilize the Bunchkaufman (LDL) and LU decomposition of `B`, now efficiently compute the generalized eigenvalues (and eigenvectors) of `A` and `B` ([#50471](https://github.com/JuliaLang/julia/pull/50471))
- Specialized dispatch for `eigvals`/`eigen(::Hermitian{<:Tridiagonal})` now performs a similarity transformation to create a real symmetric tridiagonal matrix and solves it using LAPACK routines ([#49546](https://github.com/JuliaLang/julia/pull/49546))
- Structured matrices now retain either the axes of the parent (for `Symmetric`/`Hermitian`/`AbstractTriangular`/`UpperHessenberg`) or that of the principal diagonal (for banded matrices) ([#52480](https://github.com/JuliaLang/julia/pull/52480))
- `bunchkaufman` and `bunchkaufman!` now work for any `AbstractFloat`, `Rational`, and their complex variants. `bunchkaufman` now also supports `Integer` types via internal conversion to `Rational{BigInt}`. Added new function `inertia` to compute the inertia of the diagonal factor in the Bunchkaufman factorization; for complex symmetric matrices, `inertia` only computes the number of zero eigenvalues of the diagonal factor ([#51487](https://github.com/JuliaLang/julia/pull/51487))
- Packages that specialize `mul!` with a method signature like `mul!(::AbstractMatrix, ::MyMatrix, ::AbstractMatrix, ::Number, ::Number)` no longer encounter method ambiguities with LinearAlgebra’s structured matrix types (e.g., `AbstractTriangular`) ([#52837](https://github.com/JuliaLang/julia/pull/52837))
- `lu` and `issuccess(::LU)` now accept an `allowsingular` keyword argument. When set to true, a valid factorization with rank-deficient U factor will be treated as success instead of throwing an error. Such factorizations are now shown by printing the factors together with a "rank-deficient" note rather than printing a "Failed Factorization" message  ([#52957](https://github.com/JuliaLang/julia/pull/52957))

LinearAlgebra v1.10 Release Notes
=======================

- `AbstractQ` no longer subtypes `AbstractMatrix`. Additionally, `adjoint(Q::AbstractQ)` now wraps `Q` in an `AdjointQ` instead of `Adjoint`. `AdjointQ` itself subtypes `AbstractQ`. This change clarifies that `AbstractQ` instances generally behave like function-based, matrix-backed linear operators, which may not allow efficient indexing. Many `AbstractQ` types can act on vectors/matrices of varying sizes, behaving like a matrix with context-dependent dimensions. The full API is described in the Julia documentation ([#46196](https://github.com/JuliaLang/julia/pull/46196)).
- Adjoints and transposes of `Factorization` objects are now wrapped in `AdjointFactorization` and `TransposeFactorization` types, respectively, instead of the generic `Adjoint` and `Transpose` wrappers. These new types subtype `Factorization` ([#46874](https://github.com/JuliaLang/julia/pull/46874)).
- Added new functions `hermitianpart` and `hermitianpart!` to extract the Hermitian (real symmetric) part of a matrix ([#31836](https://github.com/JuliaLang/julia/pull/31836)).
- The norm of the adjoint or transpose of an `AbstractMatrix` now returns the norm of the parent matrix by default, consistent with the behavior for `AbstractVector`s ([#49020](https://github.com/JuliaLang/julia/pull/49020)).
- `eigen(A, B)` and `eigvals(A, B)`, where one of `A` or `B` is symmetric or Hermitian, are now fully supported ([#49533](https://github.com/JuliaLang/julia/pull/49533)).
- `eigvals` and `eigen(A, cholesky(B))` now compute the generalized eigenvalues (and eigenvectors) of `A` and `B` via Cholesky decomposition for positive definite `B`. The second argument is the output of `cholesky` ([#49533](https://github.com/JuliaLang/julia/pull/49533)).

LinearAlgebra v1.9 Release Notes
=======================

- The methods `a / b` and `b \ a`, where `a` is a scalar and `b` a vector, have been removed. These were previously equivalent to `a * pinv(b)` but could be confusing ([#44358](https://github.com/JuliaLang/julia/issues/44358)).
- LinearAlgebra is now fully reliant on `libblastrampoline` (LBT) for calling BLAS and LAPACK. OpenBLAS is shipped by default. Building the system image with other BLAS/LAPACK libraries is no longer supported; instead, the LBT mechanism should be used to swap BLAS/LAPACK with vendor-provided libraries ([#44360](https://github.com/JuliaLang/julia/pull/44360)).
- `lu` now supports a new pivoting strategy `RowNonZero()` which chooses the first non-zero pivot element. This is useful for new arithmetic types and educational purposes ([#44571](https://github.com/JuliaLang/julia/pull/44571)).
- `normalize(x, p=2)` now supports any normed vector space `x`, including scalars ([#44925](https://github.com/JuliaLang/julia/pull/44925)).
- The default number of BLAS threads is now set to the number of CPU threads on ARM CPUs, and half the number of CPU threads on other architectures ([#45412](https://github.com/JuliaLang/julia/pull/45412), [#46085](https://github.com/JuliaLang/julia/pull/46085)).

LinearAlgebra v1.8 Release Notes
=======================

- The `BLAS` submodule now supports the level-2 BLAS subroutine `spr!` ([#42830](https://github.com/JuliaLang/julia/pull/42830)).
- `cholesky[!]` now supports `LinearAlgebra.PivotingStrategy` (singleton type) values as its optional pivot argument. The default remains `cholesky(A, NoPivot())` instead of `cholesky(A, RowMaximum())`; the previous `Val{true/false}`-based calls are now deprecated ([#41640](https://github.com/JuliaLang/julia/pull/41640)).
- `LinearAlgebra.jl` is now completely independent of `SparseArrays.jl`, both in source code and unit testing ([#43127](https://github.com/JuliaLang/julia/pull/43127)).  
  As a result, sparse arrays are no longer silently returned by LinearAlgebra methods applied to Base or LinearAlgebra objects. Specifically, this introduces the following breaking changes:
  - Concatenations involving special "sparse" matrices (e.g., `*diagonal`) now return dense matrices. Consequently, the `D1` and `D2` fields of `SVD` objects constructed via `getproperty` are now dense matrices.
  - `3-arg similar(::SpecialSparseMatrix, ::Type, ::Dims)` now returns a dense zero matrix. As a consequence, products of bi-, tri-, and symmetric tridiagonal matrices with each other produce dense output. Constructing 3-arg similar matrices for special "sparse" (nonstatic) matrices now fails due to the lack of `zero(::Type{Matrix{T}})`.

LinearAlgebra v1.7 Release Notes
=======================

- Use `Libblastrampoline` to pick a BLAS and LAPACK at runtime. By default it forwards to OpenBLAS in the Julia distribution. The forwarding mechanism can be used by packages to replace the BLAS and LAPACK with user preferences ([#39455](https://github.com/JuliaLang/julia/pull/39455)).
- On `aarch64`, OpenBLAS now uses an ILP64 BLAS like all other 64-bit platforms ([#39436](https://github.com/JuliaLang/julia/pull/39436)).
- OpenBLAS is updated to 0.3.13 ([#39216](https://github.com/JuliaLang/julia/pull/39216)).
- SuiteSparse is updated to 5.8.1 ([#39455](https://github.com/JuliaLang/julia/pull/39455)).
- The shape of an `UpperHessenberg` matrix is preserved under certain arithmetic operations, e.g. when multiplying or dividing by an `UpperTriangular` matrix ([#40039](https://github.com/JuliaLang/julia/pull/40039)).
- Real quasitriangular Schur factorizations `S` can now be efficiently converted to complex upper-triangular form with `Schur{Complex}(S)` ([#40573](https://github.com/JuliaLang/julia/pull/40573)).
- `cis(A)` now supports matrix arguments ([#40194](https://github.com/JuliaLang/julia/pull/40194)).
- `dot` now supports `UniformScaling` with `AbstractMatrix` ([#40250](https://github.com/JuliaLang/julia/pull/40250)).
- `qr[!]` and `lu[!]` now support `LinearAlgebra.PivotingStrategy` (singleton type) values as their optional pivot argument: defaults are `qr(A, NoPivot())` (vs. `qr(A, ColumnNorm())` for pivoting) and `lu(A, RowMaximum())` (vs. `lu(A, NoPivot())` without pivoting); the former `Val{true/false}`-based calls are deprecated ([#40623](https://github.com/JuliaLang/julia/pull/40623)).
- `det(M::AbstractMatrix{BigInt})` now calls `det_bareiss(M)`, which uses the Bareiss algorithm to calculate precise values ([#40868](https://github.com/JuliaLang/julia/pull/40868)).

LinearAlgebra v1.6 Release Notes
=======================

- New method `LinearAlgebra.issuccess(::CholeskyPivoted)` for checking whether pivoted Cholesky factorization was successful ([#36002](https://github.com/JuliaLang/julia/pull/36002))
- `UniformScaling` can now be indexed using ranges to return dense matrices and vectors ([#24359](https://github.com/JuliaLang/julia/pull/24359))
- New function `LinearAlgebra.BLAS.get_num_threads()` for getting the number of BLAS threads ([#36360](https://github.com/JuliaLang/julia/pull/36360))
- `(+)(::UniformScaling)` is now defined, making `+I` a valid unary operation ([#36784](https://github.com/JuliaLang/julia/pull/36784))
- Instances of `UniformScaling` are no longer `isequal` to matrices; previous behavior violated the rule that `isequal(x, y)` implies `hash(x) == hash(y)`
- Transposing `*Triangular` matrices now returns matrices of the opposite triangular type, consistent with `adjoint!(::*Triangular)` and `transpose!(::*Triangular)`. Packages containing methods with, e.g., `Adjoint{<:Any,<:LowerTriangular{<:Any,<:OwnMatrixType}}` should replace that by `UpperTriangular{<:Any,<:Adjoint{<:Any,<:OwnMatrixType}}` in the method signature ([#38168](https://github.com/JuliaLang/julia/pull/38168))
