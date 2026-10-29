LinearAlgebra v1.14 Release Notes
=======================

- Eigenvalues are now sorted by default for all matrix types, including `Diagonal` and triangular matrices, and all `eigen`, `eigvals` and `eigvecs` methods accept the `sortby` keyword argument ([#1477](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1477), [#1595](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1595), [#1626](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1626), [#1645](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1645), [#1608](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1608))
- `eigen` and `eigvals` of a complex matrix now always return complex eigenvalues, even if the matrix happens to be Hermitian, which makes them type-stable. Wrap the matrix in `Hermitian` to obtain real eigenvalues ([#1717](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1717))
- `eigvals` and `svdvals` now return `Float16` results for `Float16` input, consistent with `eigen` and `svd` ([#1718](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1718))
- `eigvecs` of triangular matrices now returns eigenvectors normalized to unit length, like for all other matrix types ([#1603](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1603))
- For arrays with an abstract element type (e.g., `Matrix{Real}`), factorizations, eigen- and singular value decompositions, `lu` and `\` now determine the working element type from the stored values instead of from the abstract element type alone ([#1710](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1710), [#1744](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1744))
- Added the matrix absolute value `abs(A) = sqrt(A'A)` for matrices and `UniformScaling` ([#1680](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1680), [#1699](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1699))
- `lq` now works for generic element types, and multiplication with the `Q` factors of `QRCompactWY`, `LQ` and `Hessenberg` factorizations falls back to generic implementations for non-BLAS element types ([#1701](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1701), [#1702](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1702), [#1704](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1704))
- `AbstractQ` objects are now displayed with their entries like matrices, computing only the columns that are shown ([#1705](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1705))
- `convert(::Type{<:AbstractArray}, F)` for `Factorization`s and `AbstractQ`s is deprecated in favor of the explicit `Matrix(F)`, `Array(F)` or `AbstractMatrix(F)` ([#1712](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1712))
- `lowrankupdate[!]` and `lowrankdowndate[!]` now also accept a matrix `V` to compute rank-k up- and downdates of a `Cholesky` factorization of `A ± V*V'` in a single pass ([#1731](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1731))
- `nullspace` now accepts an `SVD` object, which allows reusing an existing factorization or choosing the SVD algorithm ([#1576](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1576))
- `sqrt(A)` accepts a `check` keyword argument; by default, it now throws an error if the computed square root of a non-Hermitian matrix with multiple zero eigenvalues does not satisfy `X^2 ≈ A` ([#1558](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1558))
- Left division now supports right-hand sides whose elements are themselves vectors (e.g., a `Vector` of static vectors) ([#1446](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1446))
- `\` with LAPACK-based factorizations now preserves the array type of strided right-hand sides instead of always returning an `Array` ([#1746](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1746))
- The public function `matprod_dest` has been removed. The destination arrays of `*`, `\` and `/` are now uniformly allocated via `LinearAlgebra.matop_dest(op, A, B)`, which packages may specialize, also for `UpperHessenberg` matrices and `Factorization`s ([#1578](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1578), [#1587](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1587), [#1591](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1591), [#1606](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1606))
- The new public functions `nzrows`, `nzcols` and `nonzeroinds` return the indices of the nonzero entries in a column, row or vector, respectively. Sparse array types may specialize them to speed up generic matrix multiplication ([#1539](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1539))
- Structured broadcasting has been extended to `Symmetric` and `Hermitian` matrices, and broadcasting with zero-absorbing operations such as `*` now preserves the structure of banded matrices ([#1607](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1607), [#1636](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1636), [#1549](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1549))
- Triangular solves with `Diagonal` matrices now return triangular matrices ([#1758](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1758))
- `I * X` now works for arrays of any dimension, and `factorize(::UniformScaling)` is now defined ([#1657](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1657), [#1759](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1759))
- `hermitianpart` now supports numbers ([#1445](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1445))
- Constructing an `Adjoint` of an `Adjoint`, or a `Transpose` of a `Transpose`, now throws an `ArgumentError`; use `adjoint` and `transpose` instead ([#1739](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1739))
- `SymTridiagonal` now requires the off-diagonal to have exactly one element fewer than the diagonal ([#1569](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1569))
- `checksquare` now throws a `DimensionMismatch` for arguments that are not matrices ([#1476](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1476))
- LinearAlgebra no longer re-exports names that are already exported by Base, such as `copyto!`, `kron`, `transpose`, `\` and `/` ([#1537](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1537), [#1540](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1540))
- `eigvals` of `UpperHessenberg` matrices now uses LAPACK's Hessenberg eigensolver and is about twice as fast ([#1557](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1557))
- The eigensolvers for `SymTridiagonal` matrices now use the same strategy as LAPACK's `syevr`, falling back to bisection and inverse iteration if the default algorithm fails ([#1493](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1493))
- `norm(A, 2)` and `norm(A, 1)` for arrays of numbers now use faster native implementations ([#1719](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1719), [#1724](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1724))
- `lu` now allocates the pivot vector with the same array type as the input ([#1647](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1647))
- `LAPACK.gelsd!` and `LAPACK.gelsy!` gained non-allocating methods that accept preallocated workspace arrays ([#1714](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1714))
- BLAS calls into OpenBLAS can now be cancelled via Julia's task cancellation mechanism ([#1672](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1672))

LinearAlgebra v1.13 Release Notes
=======================

- `lu!(F::LU, A)` now reuses the storage of an existing `LU` factorization `F` to factorize a matrix `A` of the same size, for general dense matrices ([#1131](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1131))
- `svd` and `svd!` accept the keyword arguments `atol` and `rtol` to truncate small singular values, and `ldiv!(::SVD, B; atol, rtol)` allows stable SVD-based least-squares solves. `pinv` now also accepts an `SVD` object and returns the `SVD` of the pseudoinverse ([#1387](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1387), [#1398](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1398))
- New lower-level methods `mul!(C, tA, A, B, α, β)` and `mul!(C, tA, tB, A, B, α, β)`, where `tA` and `tB` are the characters `'N'`, `'T'` or `'C'`, supersede the non-public `generic_matvecmul!` and `generic_matmatmul!`. Packages providing their own storage types are advised to overload these instead of 5-arg `mul!` ([#1671](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1671))
- `eigvecs(A, eigvals)` now computes the eigenvectors of a `Symmetric` or `Hermitian` matrix `A` corresponding to the given eigenvalues ([#1268](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1268))
- `eigmin` and `eigmax` now accept `Eigen` and `GeneralizedEigen` objects ([#1468](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1468))
- `eigen` of `UpperHessenberg` matrices now uses LAPACK and is about 50% faster than the dense eigensolver ([#1599](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1599))
- `cholesky` now supports `Tridiagonal` matrices as well as `Hermitian`/`Symmetric` wrappers of `Bidiagonal`, `Tridiagonal` and `SymTridiagonal` matrices, returning `Bidiagonal` factors ([#1409](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1409))
- The functions `fillstored!` to fill the stored entries of a structured matrix, `fillband!` to fill a range of bands of a matrix, and `uplo` to access the stored triangle of `Symmetric`, `Hermitian` and `Bidiagonal` matrices are now public ([#1333](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1333), [#1345](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1345), [#1440](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1440))
- `BLAS.set_num_threads` and `BLAS.get_num_threads` are now public ([#1581](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1581))
- On Apple Silicon, the default number of BLAS threads now excludes the efficiency cores, and is one less than the total number of cores on chips without efficiency cores ([#1686](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1686))
- `UpperHessenberg` matrices now support structured broadcasting, which preserves their structure where possible ([#1325](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1325))
- `adjoint` and `transpose` of block `Diagonal` and `Bidiagonal` matrices are now lazy, so the result shares memory with the original matrix ([#1368](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1368))
- Comparing a matrix with a `UniformScaling` using `==` now checks the elements for equality, which fixes comparisons of block matrices ([#1369](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1369))
- `AbstractQ` objects can now be converted with `AbstractMatrix(Q)` and `AbstractArray(Q)` ([#1470](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1470))
- `convert` from an `AbstractMatrix` to `Bidiagonal`, `Tridiagonal` and `SymTridiagonal` is now supported ([#1212](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1212))
- The 3-arg `dot(x, A, y)` now has an efficient method for all `Symmetric` and `Hermitian` matrices `A` ([#1410](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1410))
- Generic (non-BLAS) implementations of `syrk` and `herk` speed up products like `A'A` for non-BLAS element types ([#1249](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1249))
- `SingularException` now shows a more descriptive error message ([#1204](https://github.com/JuliaLang/LinearAlgebra.jl/pull/1204))

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

LinearAlgebra v1.5 Release Notes
=======================

- The BLAS submodule now supports the level-2 BLAS subroutine `hpmv!` ([#34211](https://github.com/JuliaLang/julia/issues/34211)).
- `normalize` now supports multidimensional arrays ([#34239](https://github.com/JuliaLang/julia/issues/34239)).
- `lq` factorizations can now be used to compute the minimum-norm solution to under-determined systems ([#34350](https://github.com/JuliaLang/julia/issues/34350)).
- `sqrt(::Hermitian)` now treats slightly negative eigenvalues as zero for nearly semidefinite matrices, and accepts a new `rtol` keyword argument for this tolerance ([#35057](https://github.com/JuliaLang/julia/issues/35057)).
- The BLAS submodule now supports the level-2 BLAS subroutine `spmv!` ([#34320](https://github.com/JuliaLang/julia/issues/34320)).
- The BLAS submodule now supports the level-1 BLAS subroutine `rot!` ([#35124](https://github.com/JuliaLang/julia/issues/35124)).
- New generic `rotate!(x, y, c, s)` and `reflect!(x, y, c, s)` functions ([#35124](https://github.com/JuliaLang/julia/issues/35124)).

LinearAlgebra v1.4 Release Notes
=======================

- `qr` and `qr!` functions support `blocksize` keyword argument ([#33053](https://github.com/JuliaLang/julia/issues/33053)).
- `dot` now admits a 3-argument method `dot(x, A, y)` to compute generalized dot products `dot(x, A*y)` without computing and storing an intermediate result `A*y` ([#32739](https://github.com/JuliaLang/julia/issues/32739)).
- `ldlt` and non-pivoted `lu` now throw a new `ZeroPivotException` type ([#33372](https://github.com/JuliaLang/julia/issues/33372)).
- `cond(A, p)` with `p=1` or `p=Inf` now computes the exact condition number instead of an estimate ([#33547](https://github.com/JuliaLang/julia/issues/33547)).
- `UniformScaling` objects may now be exponentiated such that `(a*I)^x = a^x * I`.

LinearAlgebra v1.3 Release Notes
=======================

- The BLAS submodule no longer exports `dot`, which conflicted with that in LinearAlgebra ([#31838](https://github.com/JuliaLang/julia/issues/31838)).
- `diagm` and `spdiagm` now accept optional `m,n` initial arguments to specify a size ([#31654](https://github.com/JuliaLang/julia/issues/31654)).
- `Hessenberg` factorizations `H` now support efficient shifted solves `(H+µI) \ b` and determinants, and use a specialized tridiagonal factorization for Hermitian matrices. There is also a new `UpperHessenberg` matrix type ([#31853](https://github.com/JuliaLang/julia/issues/31853)).
- Added keyword argument `alg` to `svd` and `svd!` to switch between different SVD algorithms ([#31057](https://github.com/JuliaLang/julia/issues/31057)).
- Five-argument `mul!(C, A, B, α, β)` now implements in-place multiplication fused with addition C = A B α + C β ([#23919](https://github.com/JuliaLang/julia/issues/23919)).

LinearAlgebra v1.2 Release Notes
=======================

- Added keyword arguments `rtol` and `atol` to `pinv` and `nullspace` ([#29998](https://github.com/JuliaLang/julia/issues/29998)).
- `UniformScaling` instances are now callable such that e.g. `I(3)` will produce a `Diagonal` matrix ([#30298](https://github.com/JuliaLang/julia/issues/30298)).
- Eigenvalues λ of general matrices are now sorted lexicographically by (Re λ, Im λ) ([#21598](https://github.com/JuliaLang/julia/issues/21598)).
- `one` for structured matrices (`Diagonal`, `Bidiagonal`, `Tridiagonal`, `Symtridiagonal`) now preserves structure and type ([#29777](https://github.com/JuliaLang/julia/issues/29777)).
- `diagm(v)` is now shorthand for `diagm(0 => v)` ([#31125](https://github.com/JuliaLang/julia/issues/31125)).

LinearAlgebra v1.1 Release Notes
=======================

- `isdiag` and `isposdef` now support `Diagonal` and `UniformScaling` ([#29638](https://github.com/JuliaLang/julia/issues/29638)).
- Added `mul!`, `rmul!`, and `lmul!` methods for `UniformScaling` ([#29506](https://github.com/JuliaLang/julia/issues/29506)).
- `Symmetric` and `Hermitian` matrices now preserve the wrapper when scaled with a number ([#29469](https://github.com/JuliaLang/julia/issues/29469)).
- Exponentiation operator `^` now supports raising an `Irrational` to an `AbstractMatrix` power ([#29782](https://github.com/JuliaLang/julia/issues/29782)).
- Added keyword arguments `rtol` and `atol` to `rank` ([#29926](https://github.com/JuliaLang/julia/issues/29926)).
