# This file is a part of Julia. License is MIT: https://julialang.org/license

module LAPACK
@doc """
Interfaces to LAPACK subroutines.
""" LAPACK

using ..LinearAlgebra.BLAS: @blasfunc, chkuplo

using ..LinearAlgebra: libblastrampoline, BlasFloat, BlasInt, LAPACKException, DimensionMismatch,
    SingularException, PosDefException, chkstride1, checksquare, triu, tril, dot

using Base: iszero, require_one_based_indexing


# Legacy binding maintained for backwards-compatibility but new packages
# should not look at this, instead preferring to parse the output
# of BLAS.get_config()
const liblapack = libblastrampoline

#Generic LAPACK error handlers
"""
Handle only negative LAPACK error codes

*NOTE* use only if the positive error code is useful.
"""
function chkargsok(ret::BlasInt)
    if ret < 0
        throw(ArgumentError(lazy"invalid argument #$(-ret) to LAPACK call"))
    end
end

"Handle all nonzero info codes"
function chklapackerror(ret::BlasInt, f...)
    if ret == 0
        return
    elseif ret < 0
        throw(ArgumentError(lazy"invalid argument #$(-ret) to LAPACK call"))
    else # ret > 0
        chklapackerror_positive(ret, f...)
    end
end

chklapackerror_positive(ret, f...) = throw(LAPACKException(ret))

function chknonsingular(ret::BlasInt)
    if ret > 0
        throw(SingularException(ret))
    end
end

function chkposdef(ret::BlasInt)
    if ret > 0
        throw(PosDefException(ret))
    end
end

# Generic fallback function to assert that parameters are valid
# In specific cases, the following functions may be more useful
macro chkvalidparam(position::Int, param, validvalues)
    :(chkvalidparam($position, $(string(param)), $(esc(param)), $validvalues))
end
function chkvalidparam(position::Int, var::String, val, validvals)
    # mimic `repr` for chars without explicitly calling it
    # This is because `repr` introduces dynamic dispatch
    _repr(c::AbstractChar) = "'$c'"
    _repr(c) = c
    if val ∉ validvals
        throw(ArgumentError(
            lazy"argument #$position: $var must be one of $validvals, but $(_repr(val)) was passed"))
    end
    return val
end

"Check that {c}transpose is correctly specified"
function chktrans(trans::AbstractChar)
    if !(trans == 'N' || trans == 'C' || trans == 'T')
        throw(ArgumentError(lazy"trans argument must be 'N' (no transpose), 'T' (transpose), or 'C' (conjugate transpose), got '$trans'"))
    end
    trans
end

"Check that left/right hand side multiply is correctly specified"
function chkside(side::AbstractChar)
    if !(side == 'L' || side == 'R')
        throw(ArgumentError(lazy"side argument must be 'L' (left hand multiply) or 'R' (right hand multiply), got '$side'"))
    end
    side
end

"Check that unit diagonal flag is correctly specified"
function chkdiag(diag::AbstractChar)
    if !(diag == 'U' || diag =='N')
        throw(ArgumentError(lazy"diag argument must be 'U' (unit diagonal) or 'N' (non-unit diagonal), got '$diag'"))
    end
    diag
end

subsetrows(X::AbstractVector, Y::AbstractArray, k) = Y[1:k]
subsetrows(X::AbstractMatrix, Y::AbstractArray, k) = Y[1:k, :]

function chkfinite(A::AbstractMatrix)
    for a in A
        if !isfinite(a)
            throw(ArgumentError("matrix contains Infs or NaNs"))
        end
    end
    return true
end

function chkuplofinite(A::AbstractMatrix, uplo::AbstractChar)
    require_one_based_indexing(A)
    chkuplo(uplo)
    m, n = size(A)
    if uplo == 'U'
        @inbounds for j in 1:n, i in 1:j
            if !isfinite(A[i,j])
                throw(ArgumentError("matrix contains Infs or NaNs"))
            end
        end
    else
        @inbounds for j in 1:n, i in j:m
            if !isfinite(A[i,j])
                throw(ArgumentError("matrix contains Infs or NaNs"))
            end
        end
    end
end

# LAPACK version number
function version()
    major = Ref{BlasInt}(0)
    minor = Ref{BlasInt}(0)
    patch = Ref{BlasInt}(0)
    ccall((@blasfunc(ilaver_), libblastrampoline), Cvoid,
          (Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
          major, minor, patch)
    return VersionNumber(major[], minor[], patch[])
end

# (GB) general banded matrices, LU decomposition and solver
for (gbtrf, gbtrs, elty) in
    ((:dgbtrf_,:dgbtrs_,:Float64),
     (:sgbtrf_,:sgbtrs_,:Float32),
     (:zgbtrf_,:zgbtrs_,:ComplexF64),
     (:cgbtrf_,:cgbtrs_,:ComplexF32))
    @eval begin
        # SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, KL, KU, LDAB, M, N
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   AB( LDAB, * )
        function gbtrf!(kl::Integer, ku::Integer, m::Integer, AB::AbstractMatrix{$elty})
            require_one_based_indexing(AB)
            chkstride1(AB)
            n    = size(AB, 2)
            mnmn = min(m, n)
            ipiv = similar(AB, BlasInt, mnmn)
            info = Ref{BlasInt}()
            ccall((@blasfunc($gbtrf), libblastrampoline), Cvoid,
                  (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                  m, n, kl, ku, AB, max(1,stride(AB,2)), ipiv, info)
            chklapackerror(info[])
            AB, ipiv
        end

        # SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)
        # *     .. Scalar Arguments ..
        #       CHARACTER          TRANS
        #       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
        function gbtrs!(trans::AbstractChar, kl::Integer, ku::Integer, m::Integer,
                        AB::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt},
                        B::AbstractVecOrMat{$elty})
            require_one_based_indexing(AB, B)
            chkstride1(AB, B, ipiv)
            chktrans(trans)
            info = Ref{BlasInt}()
            n    = size(AB,2)
            if m != n || m != size(B,1)
                throw(DimensionMismatch(lazy"matrix AB has dimensions $(size(AB)), but right hand side matrix B has dimensions $(size(B))"))
            end
            ccall((@blasfunc($gbtrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ref{BlasInt}, Clong),
                  trans, n, kl, ku, size(B,2), AB, max(1,stride(AB,2)), ipiv,
                  B, max(1,stride(B,2)), info, 1)
            chklapackerror(info[])
            B
        end
    end
end

"""
    gbtrf!(kl, ku, m, AB) -> (AB, ipiv)

Compute the LU factorization of a banded matrix `AB`. `kl` is the first
subdiagonal containing a nonzero band, `ku` is the last superdiagonal
containing one, and `m` is the first dimension of the matrix `AB`. Returns
the LU factorization in-place and `ipiv`, the vector of pivots used.
"""
gbtrf!(kl::Integer, ku::Integer, m::Integer, AB::AbstractMatrix)

"""
    gbtrs!(trans, kl, ku, m, AB, ipiv, B)

Solve the equation `AB * X = B`. `trans` determines the orientation of `AB`. It may
be `N` (no transpose), `T` (transpose), or `C` (conjugate transpose). `kl` is the
first subdiagonal containing a nonzero band, `ku` is the last superdiagonal
containing one, and `m` is the first dimension of the matrix `AB`. `ipiv` is the vector
of pivots returned from `gbtrf!`. Returns the vector or matrix `X`, overwriting `B` in-place.
"""
gbtrs!(trans::AbstractChar, kl::Integer, ku::Integer, m::Integer, AB::AbstractMatrix, ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat)

## (GE) general matrices: balancing and back-transforming
for (gebal, gebak, elty, relty) in
    ((:dgebal_, :dgebak_, :Float64, :Float64),
     (:sgebal_, :sgebak_, :Float32, :Float32),
     (:zgebal_, :zgebak_, :ComplexF64, :Float64),
     (:cgebal_, :cgebak_, :ComplexF32, :Float32))
    @eval begin
        #     SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          JOB
        #      INTEGER            IHI, ILP, INFO, LDA, N
        #     .. Array Arguments ..
        #      DOUBLE PRECISION   A( LDA, * ), SCALE( * )
        function gebal!(job::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            @chkvalidparam 1 job ('N', 'P', 'S', 'B')
            n = checksquare(A)
            chkfinite(A) # balancing routines don't support NaNs and Infs
            ihi = Ref{BlasInt}()
            ilo = Ref{BlasInt}()
            scale = similar(A, $relty, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($gebal), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$relty}, Ref{BlasInt}, Clong),
                  job, n, A, max(1,stride(A,2)), ilo, ihi, scale, info, 1)
            chklapackerror(info[])
            ilo[], ihi[], scale
        end

        #     SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          JOB, SIDE
        #      INTEGER            IHI, ILP, INFO, LDV, M, N
        #     .. Array Arguments ..
        #      DOUBLE PRECISION   SCALE( * ), V( LDV, * )
        function gebak!(job::AbstractChar, side::AbstractChar,
                        ilo::BlasInt, ihi::BlasInt, scale::AbstractVector{$relty},
                        V::AbstractMatrix{$elty})
            require_one_based_indexing(scale, V)
            @chkvalidparam 1 job ('N', 'P', 'S', 'B')
            chkstride1(scale, V)
            chkside(side)
            chkfinite(V) # balancing routines don't support NaNs and Infs
            n = checksquare(V)
            info = Ref{BlasInt}()
            ccall((@blasfunc($gebak), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                   Ptr{$relty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                   Clong, Clong),
                  job, side, size(V,1), ilo, ihi, scale, n, V, max(1,stride(V,2)), info,
                  1, 1)
            chklapackerror(info[])
            V
        end
    end
end

"""
    gebal!(job, A) -> (ilo, ihi, scale)

Balance the matrix `A` before computing its eigensystem or Schur factorization.
`job` can be one of `N` (`A` will not be permuted or scaled), `P` (`A` will only
be permuted), `S` (`A` will only be scaled), or `B` (`A` will be both permuted
and scaled). Modifies `A` in-place and returns `ilo`, `ihi`, and `scale`. If
permuting was turned on, `A[i,j] = 0` if `j > i` and `1 < j < ilo` or `j > ihi`.
`scale` contains information about the scaling/permutations performed.
"""
gebal!(job::AbstractChar, A::AbstractMatrix)

"""
    gebak!(job, side, ilo, ihi, scale, V)

Transform the eigenvectors `V` of a matrix balanced using `gebal!` to
the unscaled/unpermuted eigenvectors of the original matrix. Modifies `V`
in-place. `side` can be `L` (left eigenvectors are transformed) or `R`
(right eigenvectors are transformed).
"""
gebak!(job::AbstractChar, side::AbstractChar, ilo::BlasInt, ihi::BlasInt, scale::AbstractVector, V::AbstractMatrix)

# (GE) general matrices, direct decompositions
#
# These mutating functions take as arguments all the values they
# return, even if the value of the function does not depend on them
# (e.g. the tau argument).  This is so that a factorization can be
# updated in place.  The condensed mutating functions, usually a
# function of A only, are defined after this block.
for (gebrd, gelqf, geqlf, geqrf, geqp3, geqrt, geqrt3, gerqf, getrf, elty, relty) in
    ((:dgebrd_,:dgelqf_,:dgeqlf_,:dgeqrf_,:dgeqp3_,:dgeqrt_,:dgeqrt3_,:dgerqf_,:dgetrf_,:Float64,:Float64),
     (:sgebrd_,:sgelqf_,:sgeqlf_,:sgeqrf_,:sgeqp3_,:sgeqrt_,:sgeqrt3_,:sgerqf_,:sgetrf_,:Float32,:Float32),
     (:zgebrd_,:zgelqf_,:zgeqlf_,:zgeqrf_,:zgeqp3_,:zgeqrt_,:zgeqrt3_,:zgerqf_,:zgetrf_,:ComplexF64,:Float64),
     (:cgebrd_,:cgelqf_,:cgeqlf_,:cgeqrf_,:cgeqp3_,:cgeqrt_,:cgeqrt3_,:cgerqf_,:cgetrf_,:ComplexF32,:Float32))
    @eval begin
        # SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,
        #                    INFO )
        # .. Scalar Arguments ..
        # INTEGER            INFO, LDA, LWORK, M, N
        # .. Array Arguments ..
        #  DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAUP( * ),
        #           TAUQ( * ), WORK( * )
        function gebrd!(A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            m, n  = size(A)
            k     = min(m, n)
            d     = similar(A, $relty, k)
            e     = similar(A, $relty, k)
            tauq  = similar(A, $elty, k)
            taup  = similar(A, $elty, k)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gebrd), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$relty}, Ptr{$relty}, Ptr{$elty}, Ptr{$elty},
                     Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                     m, n, A, max(1,stride(A,2)),
                     d, e, tauq, taup,
                     work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, d, e, tauq, taup
        end

        # SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LWORK, M, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function gelqf!(A::AbstractMatrix{$elty}, tau::AbstractVector{$elty})
            require_one_based_indexing(A, tau)
            chkstride1(A,tau)
            m     = BlasInt(size(A, 1))
            n     = BlasInt(size(A, 2))
            lda   = BlasInt(max(1,stride(A, 2)))
            if length(tau) != min(m,n)
                throw(DimensionMismatch(lazy"tau has length $(length(tau)), but needs length $(min(m,n))"))
            end
            lwork = BlasInt(-1)
            work  = Vector{$elty}(undef, 1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gelqf), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      m, n, A, lda, tau, work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, tau
        end

        # SUBROUTINE DGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LWORK, M, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function geqlf!(A::AbstractMatrix{$elty}, tau::AbstractVector{$elty})
            require_one_based_indexing(A, tau)
            chkstride1(A,tau)
            m     = BlasInt(size(A, 1))
            n     = BlasInt(size(A, 2))
            lda   = BlasInt(max(1,stride(A, 2)))
            if length(tau) != min(m,n)
                throw(DimensionMismatch(lazy"tau has length $(length(tau)), but needs length $(min(m,n))"))
            end
            lwork = BlasInt(-1)
            work  = Vector{$elty}(undef, 1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($geqlf), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      m, n, A, lda, tau, work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, tau
        end

        # SUBROUTINE DGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LWORK, M, N
        # *     .. Array Arguments ..
        #       INTEGER            JPVT( * )
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function geqp3!(A::AbstractMatrix{$elty}, jpvt::AbstractVector{BlasInt}, tau::AbstractVector{$elty})
            require_one_based_indexing(A, jpvt, tau)
            chkstride1(A,jpvt,tau)
            m,n = size(A)
            if length(tau) != min(m,n)
                throw(DimensionMismatch(lazy"tau has length $(length(tau)), but needs length $(min(m,n))"))
            end
            if length(jpvt) != n
                throw(DimensionMismatch(lazy"jpvt has length $(length(jpvt)), but needs length $n"))
            end
            lda = stride(A,2)
            if lda == 0
                return A, tau, jpvt
            end # Early exit
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            cmplx = eltype(A)<:Complex
            if cmplx
                rwork = Vector{$relty}(undef, 2n)
            end
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                if cmplx
                    ccall((@blasfunc($geqp3), libblastrampoline), Cvoid,
                          (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{$relty}, Ptr{BlasInt}),
                          m, n, A, lda,
                          jpvt, tau, work, lwork,
                          rwork, info)
                else
                    ccall((@blasfunc($geqp3), libblastrampoline), Cvoid,
                          (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{BlasInt}),
                          m, n, A, lda,
                          jpvt, tau, work,
                          lwork, info)
                end
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            return A, tau, jpvt
        end

        function geqrt!(A::AbstractMatrix{$elty}, T::AbstractMatrix{$elty})
            require_one_based_indexing(A, T)
            chkstride1(A)
            m, n = size(A)
            minmn = min(m, n)
            nb = size(T, 1)
            if nb > minmn
                throw(ArgumentError(lazy"block size $nb > $minmn too large"))
            end
            lda = max(1, stride(A,2))
            work = Vector{$elty}(undef, nb*n)
            if minmn > 0
                info = Ref{BlasInt}()
                ccall((@blasfunc($geqrt), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ptr{BlasInt}),
                     m, n, nb, A,
                     lda, T, max(1,stride(T,2)), work,
                     info)
                chklapackerror(info[])
            end
            A, T
        end

        function geqrt3!(A::AbstractMatrix{$elty}, T::AbstractMatrix{$elty})
            require_one_based_indexing(A, T)
            chkstride1(A)
            chkstride1(T)
            m, n = size(A)
            p, q = size(T)
            if m < n
                throw(DimensionMismatch(lazy"input matrix A has dimensions ($m,$n), but should have more rows than columns"))
            end
            if p != n || q != n
                throw(DimensionMismatch(lazy"block reflector T has dimensions ($p,$q), but should have dimensions ($n,$n)"))
            end
            if n > 0    # this implies `m > 0` because of `m >= n`
                info = Ref{BlasInt}()
                ccall((@blasfunc($geqrt3), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                     m, n, A, max(1, stride(A, 2)),
                     T, max(1,stride(T,2)), info)
                chklapackerror(info[])
            end
            A, T
        end

        ## geqrfp! - positive elements on diagonal of R - not defined yet
        # SUBROUTINE DGEQRFP( M, N, A, LDA, TAU, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LWORK, M, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function geqrf!(A::AbstractMatrix{$elty}, tau::AbstractVector{$elty})
            require_one_based_indexing(A, tau)
            chkstride1(A,tau)
            m, n  = size(A)
            if length(tau) != min(m,n)
                throw(DimensionMismatch(lazy"tau has length $(length(tau)), but needs length $(min(m,n))"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2                # first call returns lwork as work[1]
                ccall((@blasfunc($geqrf), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      m, n, A, max(1,stride(A,2)), tau, work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = max(BlasInt(1),BlasInt(real(work[1])))
                    resize!(work, lwork)
                end
            end
            A, tau
        end

        # SUBROUTINE DGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LWORK, M, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function gerqf!(A::AbstractMatrix{$elty},tau::AbstractVector{$elty})
            require_one_based_indexing(A, tau)
            chkstride1(A,tau)
            m, n  = size(A)
            if length(tau) != min(m,n)
                throw(DimensionMismatch(lazy"tau has length $(length(tau)), but needs length $(min(m,n))"))
            end
            lwork = BlasInt(-1)
            work  = Vector{$elty}(undef, 1)
            info  = Ref{BlasInt}()
            for i = 1:2                # first call returns lwork as work[1]
                ccall((@blasfunc($gerqf), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      m, n, A, max(1,stride(A,2)), tau, work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = max(BlasInt(m), BlasInt(real(work[1])))
                    resize!(work, lwork)
                end
            end
            A, tau
        end

        # SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, M, N
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * )
        function getrf!(A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt}; check::Bool=true)
            require_one_based_indexing(A)
            check && chkfinite(A)
            chkstride1(A)
            m, n = size(A)
            lda  = max(1,stride(A, 2))
            info = Ref{BlasInt}()
            ccall((@blasfunc($getrf), libblastrampoline), Cvoid,
                  (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                   Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                  m, n, A, lda, ipiv, info)
            chkargsok(info[])
            A, ipiv, info[] #Error code is stored in LU factorization type
        end
    end
end

"""
    gebrd!(A) -> (A, d, e, tauq, taup)

Reduce `A` in-place to bidiagonal form `A = QBP'`. Returns `A`, containing the
bidiagonal matrix `B`; `d`, containing the diagonal elements of `B`; `e`,
containing the off-diagonal elements of `B`; `tauq`, containing the
elementary reflectors representing `Q`; and `taup`, containing the
elementary reflectors representing `P`.
"""
gebrd!(A::AbstractMatrix)

"""
    gelqf!(A, tau)

Compute the `LQ` factorization of `A`, `A = LQ`. `tau` contains scalars
which parameterize the elementary reflectors of the factorization. `tau`
must have length greater than or equal to the smallest dimension of `A`.

Returns
`A` and `tau` modified in-place.
"""
gelqf!(A::AbstractMatrix, tau::AbstractVector)

"""
    geqlf!(A, tau)

Compute the `QL` factorization of `A`, `A = QL`. `tau` contains scalars
which parameterize the elementary reflectors of the factorization. `tau`
must have length greater than or equal to the smallest dimension of `A`.

Returns `A` and `tau` modified in-place.
"""
geqlf!(A::AbstractMatrix, tau::AbstractVector)

"""
    geqp3!(A, [jpvt, tau]) -> (A, tau, jpvt)

Compute the pivoted `QR` factorization of `A`, `AP = QR` using BLAS level 3.
`P` is a pivoting matrix, represented by `jpvt`. `tau` stores the elementary
reflectors. The arguments `jpvt` and `tau` are optional and allow
for passing preallocated arrays. When passed, `jpvt` must have length greater
than or equal to `n` if `A` is an `(m x n)` matrix and `tau` must have length
greater than or equal to the smallest dimension of `A`. On entry, if `jpvt[j]`
does not equal zero then the `j`th column of `A` is permuted to the front of
`AP`.

`A`, `jpvt`, and `tau` are modified in-place.
"""
geqp3!(A::AbstractMatrix, jpvt::AbstractVector{BlasInt}, tau::AbstractVector)

function geqp3!(A::AbstractMatrix{<:BlasFloat}, jpvt::AbstractVector{BlasInt})
    require_one_based_indexing(A, jpvt)
    m, n = size(A)
    geqp3!(A, jpvt, similar(A, min(m, n)))
end

function geqp3!(A::AbstractMatrix{<:BlasFloat})
    require_one_based_indexing(A)
    m, n = size(A)
    geqp3!(A, zeros(BlasInt, n), similar(A, min(m, n)))
end

"""
    geqrt!(A, T)

Compute the blocked `QR` factorization of `A`, `A = QR`. `T` contains upper
triangular block reflectors which parameterize the elementary reflectors of
the factorization. The first dimension of `T` sets the block size and it must
be between 1 and `n`. The second dimension of `T` must equal the smallest
dimension of `A`.

Returns `A` and `T` modified in-place.
"""
geqrt!(A::AbstractMatrix, T::AbstractMatrix)

"""
    geqrt3!(A, T)

Recursively computes the blocked `QR` factorization of `A`, `A = QR`. `T`
contains upper triangular block reflectors which parameterize the
elementary reflectors of the factorization.  The first dimension of `T` sets the
block size and it must be between 1 and `n`. The second dimension of `T` must
equal the smallest dimension of `A`.

Returns `A` and `T` modified in-place.
"""
geqrt3!(A::AbstractMatrix, T::AbstractMatrix)

"""
    geqrf!(A, tau)

Compute the `QR` factorization of `A`, `A = QR`. `tau` contains scalars
which parameterize the elementary reflectors of the factorization. `tau`
must have length greater than or equal to the smallest dimension of `A`.

Returns `A` and `tau` modified in-place.
"""
geqrf!(A::AbstractMatrix, tau::AbstractVector)

"""
    gerqf!(A, tau)

Compute the `RQ` factorization of `A`, `A = RQ`. `tau` contains scalars
which parameterize the elementary reflectors of the factorization. `tau`
must have length greater than or equal to the smallest dimension of `A`.

Returns `A` and `tau` modified in-place.
"""
gerqf!(A::AbstractMatrix, tau::AbstractVector)

"""
    getrf!(A, ipiv) -> (A, ipiv, info)

Compute the pivoted `LU` factorization of `A`, `A = LU`. `ipiv` contains the pivoting
information and `info` a code which indicates success (`info = 0`), a singular value
in `U` (`info = i`, in which case `U[i,i]` is singular), or an error code (`info < 0`).
"""
getrf!(A::AbstractMatrix, ipiv::AbstractVector; check::Bool=true)

"""
    gelqf!(A) -> (A, tau)

Compute the `LQ` factorization of `A`, `A = LQ`.

Returns `A`, modified in-place, and `tau`, which contains scalars
which parameterize the elementary reflectors of the factorization.
"""
gelqf!(A::AbstractMatrix{<:BlasFloat}) = ((m,n) = size(A); gelqf!(A, similar(A, min(m, n))))

"""
    geqlf!(A) -> (A, tau)

Compute the `QL` factorization of `A`, `A = QL`.

Returns `A`, modified in-place, and `tau`, which contains scalars
which parameterize the elementary reflectors of the factorization.
"""
geqlf!(A::AbstractMatrix{<:BlasFloat}) = ((m,n) = size(A); geqlf!(A, similar(A, min(m, n))))

"""
    geqrt!(A, nb) -> (A, T)

Compute the blocked `QR` factorization of `A`, `A = QR`. `nb` sets the block size
and it must be between 1 and `n`, the second dimension of `A`.

Returns `A`, modified in-place, and `T`, which contains upper
triangular block reflectors which parameterize the elementary reflectors of
the factorization.
"""
geqrt!(A::AbstractMatrix{<:BlasFloat}, nb::Integer) = geqrt!(A, similar(A, nb, minimum(size(A))))

"""
    geqrt3!(A) -> (A, T)

Recursively computes the blocked `QR` factorization of `A`, `A = QR`.

Returns `A`, modified in-place, and `T`, which contains upper triangular block
reflectors which parameterize the elementary reflectors of the factorization.
"""
geqrt3!(A::AbstractMatrix{<:BlasFloat}) = (n = size(A, 2); geqrt3!(A, similar(A, n, n)))

"""
    geqrf!(A) -> (A, tau)

Compute the `QR` factorization of `A`, `A = QR`.

Returns `A`, modified in-place, and `tau`, which contains scalars
which parameterize the elementary reflectors of the factorization.
"""
geqrf!(A::AbstractMatrix{<:BlasFloat}) = ((m,n) = size(A); geqrf!(A, similar(A, min(m, n))))

"""
    gerqf!(A) -> (A, tau)

Compute the `RQ` factorization of `A`, `A = RQ`.

Returns `A`, modified in-place, and `tau`, which contains scalars
which parameterize the elementary reflectors of the factorization.
"""
gerqf!(A::AbstractMatrix{<:BlasFloat}) = ((m,n) = size(A); gerqf!(A, similar(A, min(m, n))))

"""
    getrf!(A) -> (A, ipiv, info)

Compute the pivoted `LU` factorization of `A`, `A = LU`.

Returns `A`, modified in-place, `ipiv`, the pivoting information, and an `info`
code which indicates success (`info = 0`), a singular value in `U`
(`info = i`, in which case `U[i,i]` is singular), or an error code (`info < 0`).
"""
getrf!(A::AbstractMatrix{T}; check::Bool=true) where {T <: BlasFloat} = ((m,n) = size(A); getrf!(A, similar(A, BlasInt, min(m, n)); check))

## Tools to compute and apply elementary reflectors
for (larfg, elty) in
    ((:dlarfg_, Float64),
     (:slarfg_, Float32),
     (:zlarfg_, ComplexF64),
     (:clarfg_, ComplexF32))
    @eval begin
        #        .. Scalar Arguments ..
        #        INTEGER            incx, n
        #        DOUBLE PRECISION   alpha, tau
        #        ..
        #        .. Array Arguments ..
        #        DOUBLE PRECISION   x( * )
        function larfg!(x::AbstractVector{$elty})
            require_one_based_indexing(x)
            N    = BlasInt(length(x))
            α    = Ref{$elty}(x[1])
            incx = BlasInt(1)
            τ    = Ref{$elty}(0)
            ccall((@blasfunc($larfg), libblastrampoline), Cvoid,
                (Ref{BlasInt}, Ref{$elty}, Ptr{$elty}, Ref{BlasInt}, Ref{$elty}),
                N, α, pointer(x, 2), incx, τ)
            @inbounds x[1] = one($elty)
            return τ[]
        end
    end
end

for (larf, elty) in
    ((:dlarf_, Float64),
     (:slarf_, Float32),
     (:zlarf_, ComplexF64),
     (:clarf_, ComplexF32))
    @eval begin
        #        .. Scalar Arguments ..
        #        CHARACTER          side
        #        INTEGER            incv, ldc, m, n
        #        DOUBLE PRECISION   tau
        #        ..
        #        .. Array Arguments ..
        #        DOUBLE PRECISION   c( ldc, * ), v( * ), work( * )
        function larf!(side::AbstractChar, v::AbstractVector{$elty},
                       τ::$elty, C::AbstractMatrix{$elty}, work::AbstractVector{$elty})
            require_one_based_indexing(v, C, work)
            m, n = size(C)
            chkside(side)
            ldc = max(1, stride(C, 2))
            l = side == 'L' ? n : m
            incv  = BlasInt(1)
            ccall((@blasfunc($larf), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                 Ref{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Clong),
                side, m, n, v, incv,
                τ, C, ldc, work, 1)
            return C
        end

        function larf!(side::AbstractChar, v::AbstractVector{$elty},
                       τ::$elty, C::AbstractMatrix{$elty})
            require_one_based_indexing(v, C)
            m, n = size(C)
            chkside(side)
            lwork = side == 'L' ? n : m
            return larf!(side, v, τ, C, Vector{$elty}(undef,lwork))
        end
    end
end

## Complete orthogonaliztion tools
for (tzrzf, ormrz, elty) in
    ((:dtzrzf_,:dormrz_,:Float64),
     (:stzrzf_,:sormrz_,:Float32),
     (:ztzrzf_,:zunmrz_,:ComplexF64),
     (:ctzrzf_,:cunmrz_,:ComplexF32))
    @eval begin
         #       SUBROUTINE ZTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         #
         #       .. Scalar Arguments ..
         #       INTEGER            INFO, LDA, LWORK, M, N
         #       ..
         #       .. Array Arguments ..
         #       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
        function tzrzf!(A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            m, n = size(A)
            if n < m
                throw(DimensionMismatch(lazy"input matrix A has dimensions ($m,$n), but cannot have fewer columns than rows"))
            end
            lda = max(1, stride(A,2))
            tau = similar(A, $elty, m)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($tzrzf), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                    m, n, A, lda,
                    tau, work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, tau
        end

        #       SUBROUTINE ZUNMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,
        #                          WORK, LWORK, INFO )
        #
        #       .. Scalar Arguments ..
        #       CHARACTER          SIDE, TRANS
        #       INTEGER            INFO, K, L, LDA, LDC, LWORK, M, N
        #       ..
        #       .. Array Arguments ..
        #       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
        function ormrz!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix{$elty},
                        tau::AbstractVector{$elty}, C::AbstractMatrix{$elty})
            require_one_based_indexing(A, tau, C)
            chktrans(trans)
            chkside(side)
            chkstride1(A, tau, C)
            m, n = size(C)
            k = length(tau)
            l = size(A, 2) - size(A, 1)
            lda = max(1, stride(A,2))
            ldc = max(1, stride(C,2))
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ormrz), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                     Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                    side, trans, m, n,
                    k, l, A, lda,
                    tau, C, ldc, work,
                    lwork, info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            C
        end
    end
end

"""
    ormrz!(side, trans, A, tau, C)

Multiplies the matrix `C` by `Q` from the transformation supplied by
`tzrzf!`. Depending on `side` or `trans` the multiplication can be
left-sided (`side = L, Q*C`) or right-sided (`side = R, C*Q`) and `Q`
can be unmodified (`trans = N`), transposed (`trans = T`), or conjugate
transposed (`trans = C`). Returns matrix `C` which is modified in-place
with the result of the multiplication.
"""
ormrz!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix, tau::AbstractVector, C::AbstractMatrix)

"""
    tzrzf!(A) -> (A, tau)

Transforms the upper trapezoidal matrix `A` to upper triangular form in-place.
Returns `A` and `tau`, the scalar parameters for the elementary reflectors
of the transformation.
"""
tzrzf!(A::AbstractMatrix)

## (GE) general matrices, solvers with factorization, solver and inverse
for (gels, gesv, getrs, getri, elty) in
    ((:dgels_,:dgesv_,:dgetrs_,:dgetri_,:Float64),
     (:sgels_,:sgesv_,:sgetrs_,:sgetri_,:Float32),
     (:zgels_,:zgesv_,:zgetrs_,:zgetri_,:ComplexF64),
     (:cgels_,:cgesv_,:cgetrs_,:cgetri_,:ComplexF32))
    @eval begin
        #      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,INFO)
        # *     .. Scalar Arguments ..
        #       CHARACTER          TRANS
        #       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
        function gels!(trans::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chktrans(trans)
            chkstride1(A, B)
            btrn  = trans == 'T'
            m, n  = size(A)
            if size(B,1) != (btrn ? n : m)
                throw(DimensionMismatch(lazy"matrix A has dimensions ($m,$n), transposed: $btrn, but leading dimension of B is $(size(B,1))"))
            end
            info  = Ref{BlasInt}()
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gels), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      (btrn ? 'T' : 'N'), m, n, size(B,2), A, max(1,stride(A,2)),
                      B, max(1,stride(B,2)), work, lwork, info, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            k   = min(m, n)
            F   = m < n ? tril(A[1:k, 1:k]) : triu(A[1:k, 1:k])
            ssr = Vector{$elty}(undef, size(B, 2))
            for i = 1:size(B,2)
                x = zero($elty)
                for j = k+1:size(B,1)
                    x += abs2(B[j,i])
                end
                ssr[i] = x
            end
            F, subsetrows(B, B, k), ssr
        end

        # SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LDB, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function gesv!(A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            n = checksquare(A)
            if size(B,1) != n
                throw(DimensionMismatch(lazy"B has leading dimension $(size(B,1)), but needs $n"))
            end
            ipiv = similar(A, BlasInt, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($gesv), libblastrampoline), Cvoid,
                  (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                  n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info)
            chklapackerror(info[])
            B, A, ipiv
        end

        #     SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          TRANS
        #      INTEGER            INFO, LDA, LDB, N, NRHS
        #     .. Array Arguments ..
        #      INTEGER            IPIV( * )
        #      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function getrs!(trans::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, ipiv, B)
            chktrans(trans)
            chkstride1(A, B, ipiv)
            n = checksquare(A)
            if n != size(B, 1)
                throw(DimensionMismatch(lazy"B has leading dimension $(size(B,1)), but needs $n"))
            end
            if n != length(ipiv)
                throw(DimensionMismatch(lazy"ipiv has length $(length(ipiv)), but needs to be $n"))
            end
            nrhs = size(B, 2)
            info = Ref{BlasInt}()
            ccall((@blasfunc($getrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  trans, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info, 1)
            chklapackerror(info[])
            B
        end

        #     SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
        #*     .. Scalar Arguments ..
        #      INTEGER            INFO, LDA, LWORK, N
        #*     .. Array Arguments ..
        #      INTEGER            IPIV( * )
        #      DOUBLE PRECISION   A( LDA, * ), WORK( * )
        function getri!(A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A, ipiv)
            chkstride1(A, ipiv)
            n = checksquare(A)
            if n != length(ipiv)
                throw(DimensionMismatch(lazy"ipiv has length $(length(ipiv)), but needs $n"))
            end
            lda = max(1,stride(A, 2))
            lwork = BlasInt(-1)
            work  = Vector{$elty}(undef, 1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($getri), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      n, A, lda, ipiv, work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A
        end
    end
end

"""
    gels!(trans, A, B) -> (F, B, ssr)

Solves the linear equation `A * X = B`, `transpose(A) * X = B`, or `adjoint(A) * X = B` using
a QR or LQ factorization. Modifies the matrix/vector `B` in place with the
solution. `A` is overwritten with its `QR` or `LQ` factorization. `trans`
may be one of `N` (no modification), `T` (transpose), or `C` (conjugate
transpose). `gels!` searches for the minimum norm/least squares solution.
`A` may be under or over determined. The solution is returned in `B`.
"""
gels!(trans::AbstractChar, A::AbstractMatrix, B::AbstractVecOrMat)

"""
    gesv!(A, B) -> (B, A, ipiv)

Solves the linear equation `A * X = B` where `A` is a square matrix using
the `LU` factorization of `A`. `A` is overwritten with its `LU`
factorization and `B` is overwritten with the solution `X`. `ipiv` contains the
pivoting information for the `LU` factorization of `A`.
"""
gesv!(A::AbstractMatrix, B::AbstractVecOrMat)

"""
    getrs!(trans, A, ipiv, B)

Solves the linear equation `A * X = B`, `transpose(A) * X = B`, or `adjoint(A) * X = B` for
square `A`. Modifies the matrix/vector `B` in place with the solution. `A`
is the `LU` factorization from `getrf!`, with `ipiv` the pivoting
information. `trans` may be one of `N` (no modification), `T` (transpose),
or `C` (conjugate transpose).
"""
getrs!(trans::AbstractChar, A::AbstractMatrix, ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat)

"""
    getri!(A, ipiv)

Computes the inverse of `A`, using its `LU` factorization found by
`getrf!`. `ipiv` is the pivot information output and `A`
contains the `LU` factorization of `getrf!`. `A` is overwritten with
its inverse.
"""
getri!(A::AbstractMatrix, ipiv::AbstractVector{BlasInt})

for (gesvx, elty) in
    ((:dgesvx_,:Float64),
     (:sgesvx_,:Float32))
    @eval begin
        #     SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
        #                        EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
        #                        WORK, IWORK, INFO )
        #
        #     .. Scalar Arguments ..
        #     CHARACTER          EQUED, FACT, TRANS
        #     INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
        #     DOUBLE PRECISION   RCOND
        #     ..
        #     .. Array Arguments ..
        #     INTEGER            IPIV( * ), IWORK( * )
        #     DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
        #    $                   BERR( * ), C( * ), FERR( * ), R( * ),
        #    $                   WORK( * ), X( LDX, *
        #
        function gesvx!(fact::AbstractChar, trans::AbstractChar, A::AbstractMatrix{$elty},
                        AF::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt}, equed::AbstractChar,
                        R::AbstractVector{$elty}, C::AbstractVector{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, AF, ipiv, R, C, B)
            @chkvalidparam 1 fact ('F', 'N', 'E')
            chktrans(trans)
            chkstride1(ipiv, R, C, B)
            n    = checksquare(A)
            lda  = stride(A,2)
            n    = checksquare(AF)
            ldaf = stride(AF,2)
            nrhs = size(B,2)
            ldb  = stride(B,2)
            rcond = Ref{$elty}()
            ferr  = similar(A, $elty, nrhs)
            berr  = similar(A, $elty, nrhs)
            work  = Vector{$elty}(undef, 4n)
            iwork = Vector{BlasInt}(undef, n)
            info  = Ref{BlasInt}()
            X = similar(A, $elty, n, nrhs)
            ccall((@blasfunc($gesvx), libblastrampoline), Cvoid,
              (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
               Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
               Ref{UInt8}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
               Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty},
               Ptr{$elty}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
              fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B,
              ldb, X, n, rcond, ferr, berr, work, iwork, info, 1, 1, 1)
            chklapackerror(info[])
            if info[] == n + 1
                @warn "Matrix is singular to working precision"
            else
                chknonsingular(info[])
            end
            #WORK(1) contains the reciprocal pivot growth factor norm(A)/norm(U)
            X, equed, R, C, B, rcond[], ferr, berr, work[1]
        end

        function gesvx!(A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            n = size(A,1)
            X, equed, R, C, B, rcond, ferr, berr, rpgf =
                gesvx!('N', 'N', A,
                       similar(A, $elty, n, n),
                       similar(A, BlasInt, n),
                       'N',
                       similar(A, $elty, n),
                       similar(A, $elty, n),
                       B)
            X, rcond, ferr, berr, rpgf
        end
    end
end
for (gesvx, elty, relty) in
    ((:zgesvx_,:ComplexF64,:Float64),
     (:cgesvx_,:ComplexF32 ,:Float32))
    @eval begin
        #     SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
        #                        EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
        #                        WORK, RWORK, INFO )
        #
        #     .. Scalar Arguments ..
        #     CHARACTER          EQUED, FACT, TRANS
        #     INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
        #     DOUBLE PRECISION   RCOND
        #     ..
        #     .. Array Arguments ..
        #     INTEGER            IPIV( * )
        #     DOUBLE PRECISION   BERR( * ), C( * ), FERR( * ), R( * ),
        #    $                   RWORK( * )
        #     COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
        #    $                   WORK( * ), X( LDX, * )
        function gesvx!(fact::AbstractChar, trans::AbstractChar, A::AbstractMatrix{$elty},
                        AF::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt}, equed::AbstractChar,
                        R::AbstractVector{$relty}, C::AbstractVector{$relty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, AF, ipiv, R, C, B)
            @chkvalidparam 1 fact ('F', 'N', 'E')
            chktrans(trans)
            chkstride1(A, AF, ipiv, R, C, B)
            n   = checksquare(A)
            lda = stride(A,2)
            n   = checksquare(AF)
            ldaf = stride(AF,2)
            nrhs = size(B,2)
            ldb = stride(B,2)
            rcond = Ref{$relty}()
            ferr  = similar(A, $relty, nrhs)
            berr  = similar(A, $relty, nrhs)
            work  = Vector{$elty}(undef, 2n)
            rwork = Vector{$relty}(undef, 2n)
            info  = Ref{BlasInt}()
            X = similar(A, $elty, n, nrhs)
            ccall((@blasfunc($gesvx), libblastrampoline), Cvoid,
              (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
               Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
               Ref{UInt8}, Ptr{$relty}, Ptr{$relty}, Ptr{$elty}, Ref{BlasInt},
               Ptr{$elty}, Ref{BlasInt}, Ptr{$relty}, Ptr{$relty}, Ptr{$relty},
               Ptr{$elty}, Ptr{$relty}, Ref{BlasInt}, Clong, Clong, Clong),
              fact, trans, n, nrhs, A, lda, AF, ldaf, ipiv, equed, R, C, B,
              ldb, X, n, rcond, ferr, berr, work, rwork, info, 1, 1, 1)
            chklapackerror(info[])
            if info[] == n + 1
                @warn "Matrix is singular to working precision"
            else
                chknonsingular(info[])
            end
            #RWORK(1) contains the reciprocal pivot growth factor norm(A)/norm(U)
            X, equed, R, C, B, rcond[], ferr, berr, rwork[1]
        end

        #Wrapper for the no-equilibration, no-transpose calculation
        function gesvx!(A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            n = size(A,1)
            X, equed, R, C, B, rcond, ferr, berr, rpgf =
                gesvx!('N', 'N', A,
                       similar(A, $elty, n, n),
                       similar(A, BlasInt, n),
                       'N',
                       similar(A, $relty, n),
                       similar(A, $relty, n),
                       B)
            X, rcond, ferr, berr, rpgf
        end
    end
end

"""
    gesvx!(fact, trans, A, AF, ipiv, equed, R, C, B) -> (X, equed, R, C, B, rcond, ferr, berr, work)

Solves the linear equation `A * X = B` (`trans = N`), `transpose(A) * X = B`
(`trans = T`), or `adjoint(A) * X = B` (`trans = C`) using the `LU` factorization
of `A`. `fact` may be `E`, in which case `A` will be equilibrated and copied
to `AF`; `F`, in which case `AF` and `ipiv` from a previous `LU` factorization
are inputs; or `N`, in which case `A` will be copied to `AF` and then
factored. If `fact = F`, `equed` may be `N`, meaning `A` has not been
equilibrated; `R`, meaning `A` was multiplied by `Diagonal(R)` from the left;
`C`, meaning `A` was multiplied by `Diagonal(C)` from the right; or `B`, meaning
`A` was multiplied by `Diagonal(R)` from the left and `Diagonal(C)` from the right.
If `fact = F` and `equed = R` or `B` the elements of `R` must all be positive.
If `fact = F` and `equed = C` or `B` the elements of `C` must all be positive.

Returns the solution `X`; `equed`, which is an output if `fact` is not `N`,
and describes the equilibration that was performed; `R`, the row equilibration
diagonal; `C`, the column equilibration diagonal; `B`, which may be overwritten
with its equilibrated form `Diagonal(R)*B` (if `trans = N` and `equed = R,B`) or
`Diagonal(C)*B` (if `trans = T,C` and `equed = C,B`); `rcond`, the reciprocal
condition number of `A` after equilbrating; `ferr`, the forward error bound for
each solution vector in `X`; `berr`, the forward error bound for each solution
vector in `X`; and `work`, the reciprocal pivot growth factor.
"""
gesvx!(fact::AbstractChar, trans::AbstractChar, A::AbstractMatrix, AF::AbstractMatrix,
    ipiv::AbstractVector{BlasInt}, equed::AbstractChar, R::AbstractVector, C::AbstractVector, B::AbstractVecOrMat)

"""
    gesvx!(A, B)

The no-equilibration, no-transpose simplification of `gesvx!`.
"""
gesvx!(A::AbstractMatrix, B::AbstractVecOrMat)

for (gelsd, gelsy, elty) in
    ((:dgelsd_,:dgelsy_,:Float64),
     (:sgelsd_,:sgelsy_,:Float32))
    @eval begin
        # SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
        #      $                   WORK, LWORK, IWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
        #       DOUBLE PRECISION   RCOND
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
        function gelsd!(A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty}, rcond::Real=-one($elty))
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            m, n  = size(A)
            if size(B, 1) != m
                throw(DimensionMismatch(lazy"B has leading dimension $(size(B,1)) but needs $m"))
            end
            newB = [B; zeros($elty, max(0, n - size(B, 1)), size(B, 2))]
            s     = similar(A, $elty, min(m, n))
            rnk   = Ref{BlasInt}()
            info  = Ref{BlasInt}()
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            iwork = Vector{BlasInt}(undef, 1)
            for i = 1:2  # first call returns lwork as work[1] and iwork length as iwork[1]
                ccall((@blasfunc($gelsd), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$elty}, Ref{$elty}, Ref{BlasInt}, Ptr{$elty},
                       Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
                      m, n, size(B,2),
                      A, max(1,stride(A,2)), newB, max(1,stride(B,2),n),
                      s, $elty(rcond), rnk, work,
                      lwork, iwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    resize!(iwork, iwork[1])
                end
            end
            subsetrows(B, newB, n), rnk[]
        end

        #       SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,
        #      $                   WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
        #       DOUBLE PRECISION   RCOND
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            JPVT( * )
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
        function gelsy!(A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty}, rcond::Real=eps($elty))
            require_one_based_indexing(A, B)
            chkstride1(A)
            m = size(A, 1)
            n = size(A, 2)
            nrhs = size(B, 2)
            if size(B, 1) != m
                throw(DimensionMismatch(lazy"B has leading dimension $(size(B,1)) but needs $m"))
            end
            newB = [B; zeros($elty, max(0, n - size(B, 1)), size(B, 2))]
            lda = max(1, stride(A,2))
            ldb = max(1, stride(newB,2))
            jpvt = zeros(BlasInt, n)
            rnk = Ref{BlasInt}()
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gelsy), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                     Ref{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{BlasInt}),
                    m, n, nrhs, A,
                    lda, newB, ldb, jpvt,
                    $elty(rcond), rnk, work, lwork,
                    info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                end
            end
            subsetrows(B, newB, n), rnk[]
        end
    end
end

for (gelsd, gelsy, elty, relty) in
    ((:zgelsd_,:zgelsy_,:ComplexF64,:Float64),
     (:cgelsd_,:cgelsy_,:ComplexF32,:Float32))
    @eval begin
        # SUBROUTINE ZGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
        #      $                   WORK, LWORK, RWORK, IWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
        #       DOUBLE PRECISION   RCOND
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   RWORK( * ), S( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
        function gelsd!(A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty}, rcond::Real=-one($relty))
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            m, n  = size(A)
            if size(B, 1) != m
                throw(DimensionMismatch(lazy"B has leading dimension $(size(B,1)) but needs $m"))
            end
            newB = [B; zeros($elty, max(0, n - size(B, 1)), size(B, 2))]
            s     = similar(A, $relty, min(m, n))
            rnk   = Ref{BlasInt}()
            info  = Ref{BlasInt}()
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, 1)
            iwork = Vector{BlasInt}(undef, 1)
            for i = 1:2  # first call returns lwork as work[1], rwork length as rwork[1] and iwork length as iwork[1]
                ccall((@blasfunc($gelsd), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                       Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$relty},
                       Ref{$relty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$relty}, Ref{BlasInt}, Ref{BlasInt}),
                      m, n, size(B,2), A,
                      max(1,stride(A,2)), newB, max(1,stride(B,2),n), s,
                      $relty(rcond), rnk, work, lwork,
                      rwork, iwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    resize!(rwork, BlasInt(rwork[1]))
                    resize!(iwork, iwork[1])
                end
            end
            subsetrows(B, newB, n), rnk[]
        end

        #       SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,
        #      $                   WORK, LWORK, RWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
        #       DOUBLE PRECISION   RCOND
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            JPVT( * )
        #       DOUBLE PRECISION   RWORK( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
        function gelsy!(A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty}, rcond::Real=eps($relty))
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            m, n = size(A)
            nrhs = size(B, 2)
            if size(B, 1) != m
                throw(DimensionMismatch(lazy"B has leading dimension $(size(B,1)) but needs $m"))
            end
            newB = [B; zeros($elty, max(0, n - size(B, 1)), size(B, 2))]
            lda = max(1, m)
            ldb = max(1, m, n)
            jpvt = zeros(BlasInt, n)
            rnk = Ref{BlasInt}(1)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, 2n)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gelsy), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                     Ref{$relty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$relty}, Ptr{BlasInt}),
                    m, n, nrhs, A,
                    lda, newB, ldb, jpvt,
                    $relty(rcond), rnk, work, lwork,
                    rwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            subsetrows(B, newB, n), rnk[]
        end
    end
end

"""
    gelsd!(A, B, rcond) -> (B, rnk)

Computes the least norm solution of `A * X = B` by finding the `SVD`
factorization of `A`, then dividing-and-conquering the problem. `B`
is overwritten with the solution `X`. Singular values below `rcond`
will be treated as zero. Returns the solution in `B` and the effective rank
of `A` in `rnk`.
"""
gelsd!(A::AbstractMatrix, B::AbstractVecOrMat, rcond::Real)

"""
    gelsy!(A, B, rcond) -> (B, rnk)

Computes the least norm solution of `A * X = B` by finding the full `QR`
factorization of `A`, then dividing-and-conquering the problem. `B`
is overwritten with the solution `X`. Singular values below `rcond`
will be treated as zero. Returns the solution in `B` and the effective rank
of `A` in `rnk`.
"""
gelsy!(A::AbstractMatrix, B::AbstractVecOrMat, rcond::Real)

for (gglse, elty) in ((:dgglse_, :Float64),
                      (:sgglse_, :Float32),
                      (:zgglse_, :ComplexF64),
                      (:cgglse_, :ComplexF32))
    @eval begin
        # SUBROUTINE DGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK,
        #      $                   INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, LDA, LDB, LWORK, M, N, P
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( * ), D( * ),
        #      $                   WORK( * ), X( * )
        function gglse!(A::AbstractMatrix{$elty}, c::AbstractVector{$elty},
                        B::AbstractMatrix{$elty}, d::AbstractVector{$elty})
            require_one_based_indexing(A, c, B, d)
            chkstride1(A, c, B, d)
            m, n = size(A)
            p = size(B, 1)
            if size(B, 2) != n
                throw(DimensionMismatch(lazy"B has second dimension $(size(B,2)), needs $n"))
            end
            if length(c) != m
                throw(DimensionMismatch(lazy"c has length $(length(c)), needs $m"))
            end
            if length(d) != p
                throw(DimensionMismatch(lazy"d has length $(length(d)), needs $p"))
            end
            X = zeros($elty, n)
            info  = Ref{BlasInt}()
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gglse), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                       Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                       Ptr{$elty}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{BlasInt}),
                      m, n, p, A, max(1,stride(A,2)), B, max(1,stride(B,2)), c, d, X,
                      work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            X, dot(view(c, n - p + 1:m), view(c, n - p + 1:m))
        end
    end
end

"""
    gglse!(A, c, B, d) -> (X,res)

Solves the equation `A * x = c` where `x` is subject to the equality
constraint `B * x = d`. Uses the formula `||c - A*x||^2 = 0` to solve.
Returns `X` and the residual sum-of-squares.
"""
gglse!(A::AbstractMatrix, c::AbstractVector, B::AbstractMatrix, d::AbstractVector)

# (GE) general matrices eigenvalue-eigenvector and singular value decompositions
for (geev, gesvd, gesdd, ggsvd, elty, relty) in
    ((:dgeev_,:dgesvd_,:dgesdd_,:dggsvd_,:Float64,:Float64),
     (:sgeev_,:sgesvd_,:sgesdd_,:sggsvd_,:Float32,:Float32),
     (:zgeev_,:zgesvd_,:zgesdd_,:zggsvd_,:ComplexF64,:Float64),
     (:cgeev_,:cgesvd_,:cgesdd_,:cggsvd_,:ComplexF32,:Float32))
    @eval begin
        #      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
        #      $                  LDVR, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVL, JOBVR
        #       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
        #      $                   WI( * ), WORK( * ), WR( * )
        function geev!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            @chkvalidparam 1 jobvl ('N', 'V')
            @chkvalidparam 2 jobvr ('N', 'V')
            chkfinite(A) # balancing routines don't support NaNs and Infs
            lvecs = jobvl == 'V'
            rvecs = jobvr == 'V'
            VL    = similar(A, $elty, (n, lvecs ? n : 0))
            VR    = similar(A, $elty, (n, rvecs ? n : 0))
            cmplx = eltype(A) <: Complex
            if cmplx
                W     = similar(A, $elty, n)
                rwork = similar(A, $relty, 2n)
            else
                WR    = similar(A, $elty, n)
                WI    = similar(A, $elty, n)
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                if cmplx
                    ccall((@blasfunc($geev), libblastrampoline), Cvoid,
                          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{$relty}, Ref{BlasInt}, Clong, Clong),
                          jobvl, jobvr, n, A, max(1,stride(A,2)), W, VL, n, VR, n,
                          work, lwork, rwork, info, 1, 1)
                else
                    ccall((@blasfunc($geev), libblastrampoline), Cvoid,
                          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                           Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                          jobvl, jobvr, n, A, max(1,stride(A,2)), WR, WI, VL, n,
                          VR, n, work, lwork, info, 1, 1)
                end
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            cmplx ? (W, VL, VR) : (WR, WI, VL, VR)
        end

        #    SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK,
        #                   LWORK, IWORK, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          JOBZ
        #      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
        #*     ..
        #*     .. Array Arguments ..
        #      INTEGER            IWORK( * )
        #      DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
        #                        VT( LDVT, * ), WORK( * )
        function gesdd!(job::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            @chkvalidparam 1 job ('A', 'S', 'O', 'N')
            m, n   = size(A)
            minmn  = min(m, n)
            if job == 'A'
                U  = similar(A, $elty, (m, m))
                VT = similar(A, $elty, (n, n))
            elseif job == 'S'
                U  = similar(A, $elty, (m, minmn))
                VT = similar(A, $elty, (minmn, n))
            elseif job == 'O'
                U  = similar(A, $elty, (m, m >= n ? 0 : m))
                VT = similar(A, $elty, (n, m >= n ? n : 0))
            else
                U  = similar(A, $elty, (m, 0))
                VT = similar(A, $elty, (n, 0))
            end
            work   = Vector{$elty}(undef, 1)
            lwork  = BlasInt(-1)
            S      = similar(A, $relty, minmn)
            cmplx  = eltype(A)<:Complex
            if cmplx
                rwork = Vector{$relty}(undef, job == 'N' ? 7*minmn : minmn*max(5*minmn+7, 2*max(m,n)+2*minmn+1))
            end
            iwork  = Vector{BlasInt}(undef, 8*minmn)
            info   = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                if cmplx
                    ccall((@blasfunc($gesdd), libblastrampoline), Cvoid,
                          (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$relty}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{$relty}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                          job, m, n, A, max(1,stride(A,2)), S, U, max(1,stride(U,2)), VT, max(1,stride(VT,2)),
                          work, lwork, rwork, iwork, info, 1)
                else
                    ccall((@blasfunc($gesdd), libblastrampoline), Cvoid,
                          (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                           Ptr{BlasInt}, Ref{BlasInt}, Clong),
                          job, m, n, A, max(1,stride(A,2)), S, U, max(1,stride(U,2)), VT, max(1,stride(VT,2)),
                          work, lwork, iwork, info, 1)
                end
                chklapackerror(info[])
                if i == 1
                    # Work around issue with truncated Float32 representation of lwork in
                    # sgesdd by using nextfloat. See
                    # http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=13&t=4587&p=11036&hilit=sgesdd#p11036
                    # and
                    # https://github.com/scipy/scipy/issues/5401
                    lwork = round(BlasInt, nextfloat(real(work[1])))
                    resize!(work, lwork)
                end
            end
            if job == 'O'
                if m >= n
                    return (A, S, VT)
                else
                    # ()__
                    # ||::Z__
                    # ||::|:::Z____
                    # ||::|:::|====|
                    # ||==|===|====|
                    # ||""|===|====|
                    # ||  `"""|====|
                    # ||      `""""`
                    return (U, S, A)
                end
            end
            return (U, S, VT)
        end

        # SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBU, JOBVT
        #       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
        #      $                   VT( LDVT, * ), WORK( * )
        function gesvd!(jobu::AbstractChar, jobvt::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            @chkvalidparam 1 jobu ('A', 'S', 'O', 'N')
            @chkvalidparam 2 jobvt ('A', 'S', 'O', 'N')
            (jobu == jobvt == 'O') && throw(ArgumentError("jobu and jobvt cannot both be O"))
            m, n   = size(A)
            minmn  = min(m, n)
            S      = similar(A, $relty, minmn)
            U      = similar(A, $elty, jobu  == 'A' ? (m, m) : (jobu  == 'S' ? (m, minmn) : (m, 0)))
            VT     = similar(A, $elty, jobvt == 'A' ? (n, n) : (jobvt == 'S' ? (minmn, n) : (n, 0)))
            work   = Vector{$elty}(undef, 1)
            cmplx  = eltype(A) <: Complex
            if cmplx
                rwork = Vector{$relty}(undef, 5minmn)
            end
            lwork  = BlasInt(-1)
            info   = Ref{BlasInt}()
            for i in 1:2  # first call returns lwork as work[1]
                if cmplx
                    ccall((@blasfunc($gesvd), libblastrampoline), Cvoid,
                          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                           Ptr{$elty}, Ref{BlasInt}, Ptr{$relty}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$relty}, Ref{BlasInt}, Clong, Clong),
                          jobu, jobvt, m, n, A, max(1,stride(A,2)), S, U, max(1,stride(U,2)), VT, max(1,stride(VT,2)),
                          work, lwork, rwork, info, 1, 1)
                else
                    ccall((@blasfunc($gesvd), libblastrampoline), Cvoid,
                          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                           Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty},
                           Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                           Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                          jobu, jobvt, m, n, A, max(1,stride(A,2)), S, U, max(1,stride(U,2)), VT, max(1,stride(VT,2)),
                          work, lwork, info, 1, 1)
                end
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            if jobu == 'O'
                return (A, S, VT)
            elseif jobvt == 'O'
                    # =============|===========|()
                                   # # # #::::::
                                   # # # #::::::
                                   # # # #::::::
                                   # # # #::::::
                                   # # # # # # #
                                   # # # # # # #
                                   # # # # # # #
                return (U, S, A)   # # # # # # #
            else                   # # # # # # #
                return (U, S, VT)  # # # # # # #

            end
        end

        #       SUBROUTINE ZGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,
        #      $                   LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK,
        #      $                   RWORK, IWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBQ, JOBU, JOBV
        #       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   ALPHA( * ), BETA( * ), RWORK( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
        #      $                   U( LDU, * ), V( LDV, * ), WORK( * )
        function ggsvd!(jobu::AbstractChar, jobv::AbstractChar, jobq::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            @chkvalidparam 1 jobu ('U', 'N')
            @chkvalidparam 2 jobv ('V', 'N')
            @chkvalidparam 3 jobq ('Q', 'N')
            m, n = size(A)
            if size(B, 2) != n
                throw(DimensionMismatch(lazy"B has second dimension $(size(B,2)) but needs $n"))
            end
            p = size(B, 1)
            k = Vector{BlasInt}(undef, 1)
            l = Vector{BlasInt}(undef, 1)
            lda = max(1,stride(A, 2))
            ldb = max(1,stride(B, 2))
            alpha = similar(A, $relty, n)
            beta = similar(A, $relty, n)
            ldu = max(1, m)
            U = jobu == 'U' ? similar(A, $elty, ldu, m) : similar(A, $elty, 0)
            ldv = max(1, p)
            V = jobv == 'V' ? similar(A, $elty, ldv, p) : similar(A, $elty, 0)
            ldq = max(1, n)
            Q = jobq == 'Q' ? similar(A, $elty, ldq, n) : similar(A, $elty, 0)
            work = Vector{$elty}(undef, max(3n, m, p) + n)
            cmplx = eltype(A) <: Complex
            if cmplx
                rwork = Vector{$relty}(undef, 2n)
            end
            iwork = Vector{BlasInt}(undef, n)
            info = Ref{BlasInt}()
            if cmplx
                ccall((@blasfunc($ggsvd), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                    Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$relty}, Ptr{$relty}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ptr{$relty}, Ptr{BlasInt}, Ptr{BlasInt},
                    Clong, Clong, Clong),
                    jobu, jobv, jobq, m,
                    n, p, k, l,
                    A, lda, B, ldb,
                    alpha, beta, U, ldu,
                    V, ldv, Q, ldq,
                    work, rwork, iwork, info,
                    1, 1, 1)
            else
                ccall((@blasfunc($ggsvd), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                    Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$relty}, Ptr{$relty}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt},
                    Clong, Clong, Clong),
                    jobu, jobv, jobq, m,
                    n, p, k, l,
                    A, lda, B, ldb,
                    alpha, beta, U, ldu,
                    V, ldv, Q, ldq,
                    work, iwork, info,
                    1, 1, 1)
            end
            chklapackerror(info[])
            if m - k[1] - l[1] >= 0
                R = triu(A[1:k[1] + l[1],n - k[1] - l[1] + 1:n])
            else
                R = triu([A[1:m, n - k[1] - l[1] + 1:n]; B[m - k[1] + 1:l[1], n - k[1] - l[1] + 1:n]])
            end
            U, V, Q, alpha, beta, k[1], l[1], R
        end
    end
end

"""
    geev!(jobvl, jobvr, A) -> (W, VL, VR)

Finds the eigensystem of `A`. If `jobvl = N`, the left eigenvectors of
`A` aren't computed. If `jobvr = N`, the right eigenvectors of `A`
aren't computed. If `jobvl = V` or `jobvr = V`, the corresponding
eigenvectors are computed. Returns the eigenvalues in `W`, the right
eigenvectors in `VR`, and the left eigenvectors in `VL`.
"""
geev!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix)

"""
    gesdd!(job, A) -> (U, S, VT)

Finds the singular value decomposition of `A`, `A = U * S * V'`,
using a divide and conquer approach. If `job = A`, all the columns of `U` and
the rows of `V'` are computed. If `job = N`, no columns of `U` or rows of `V'`
are computed. If `job = O`, `A` is overwritten with the columns of (thin) `U`
and the rows of (thin) `V'`. If `job = S`, the columns of (thin) `U` and the
rows of (thin) `V'` are computed and returned separately.
"""
gesdd!(job::AbstractChar, A::AbstractMatrix)

"""
    gesvd!(jobu, jobvt, A) -> (U, S, VT)

Finds the singular value decomposition of `A`, `A = U * S * V'`.
If `jobu = A`, all the columns of `U` are computed. If `jobvt = A` all the rows
of `V'` are computed. If `jobu = N`, no columns of `U` are computed. If
`jobvt = N` no rows of `V'` are computed. If `jobu = O`, `A` is overwritten with
the columns of (thin) `U`. If `jobvt = O`, `A` is overwritten with the rows
of (thin) `V'`. If `jobu = S`, the columns of (thin) `U` are computed
and returned separately. If `jobvt = S` the rows of (thin) `V'` are
computed and returned separately. `jobu` and `jobvt` can't both be `O`.

Returns `U`, `S`, and `Vt`, where `S` are the singular values of `A`.
"""
gesvd!(jobu::AbstractChar, jobvt::AbstractChar, A::AbstractMatrix)

"""
    ggsvd!(jobu, jobv, jobq, A, B) -> (U, V, Q, alpha, beta, k, l, R)

Finds the generalized singular value decomposition of `A` and `B`, `U'*A*Q = D1*R`
and `V'*B*Q = D2*R`. `D1` has `alpha` on its diagonal and `D2` has `beta` on its
diagonal. If `jobu = U`, the orthogonal/unitary matrix `U` is computed. If
`jobv = V` the orthogonal/unitary matrix `V` is computed. If `jobq = Q`,
the orthogonal/unitary matrix `Q` is computed. If `jobu`, `jobv` or `jobq` is
`N`, that matrix is not computed. This function is only available in LAPACK
versions prior to 3.6.0.
"""
ggsvd!(jobu::AbstractChar, jobv::AbstractChar, jobq::AbstractChar, A::AbstractMatrix, B::AbstractMatrix)


for (f, elty) in ((:dggsvd3_, :Float64),
                  (:sggsvd3_, :Float32))
    @eval begin
        function ggsvd3!(jobu::AbstractChar, jobv::AbstractChar, jobq::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            @chkvalidparam 1 jobu ('U', 'N')
            @chkvalidparam 2 jobv ('V', 'N')
            @chkvalidparam 3 jobq ('Q', 'N')
            m, n = size(A)
            if size(B, 2) != n
                throw(DimensionMismatch(lazy"B has second dimension $(size(B,2)) but needs $n"))
            end
            p = size(B, 1)
            k = Ref{BlasInt}()
            l = Ref{BlasInt}()
            lda = max(1, stride(A, 2))
            ldb = max(1, stride(B, 2))
            alpha = similar(A, $elty, n)
            beta = similar(A, $elty, n)
            ldu = max(1, m)
            U = jobu == 'U' ? similar(A, $elty, ldu, m) : similar(A, $elty, 0)
            ldv = max(1, p)
            V = jobv == 'V' ? similar(A, $elty, ldv, p) : similar(A, $elty, 0)
            ldq = max(1, n)
            Q = jobq == 'Q' ? similar(A, $elty, ldq, n) : similar(A, $elty, 0)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            iwork = Vector{BlasInt}(undef, n)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($f), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                    Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                    Clong, Clong, Clong),
                    jobu, jobv, jobq, m,
                    n, p, k, l,
                    A, lda, B, ldb,
                    alpha, beta, U, ldu,
                    V, ldv, Q, ldq,
                    work, lwork, iwork, info,
                    1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                end
            end
            if m - k[] - l[] >= 0
                R = triu(A[1:k[] + l[],n - k[] - l[] + 1:n])
            else
                R = triu([A[1:m, n - k[] - l[] + 1:n]; B[m - k[] + 1:l[], n - k[] - l[] + 1:n]])
            end
            return U, V, Q, alpha, beta, k[], l[], R
        end
    end
end

for (f, elty, relty) in ((:zggsvd3_, :ComplexF64, :Float64),
                         (:cggsvd3_, :ComplexF32, :Float32))
    @eval begin
        function ggsvd3!(jobu::AbstractChar, jobv::AbstractChar, jobq::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            @chkvalidparam 1 jobu ('U', 'N')
            @chkvalidparam 2 jobv ('V', 'N')
            @chkvalidparam 3 jobq ('Q', 'N')
            m, n = size(A)
            if size(B, 2) != n
                throw(DimensionMismatch(lazy"B has second dimension $(size(B,2)) but needs $n"))
            end
            p = size(B, 1)
            k = Vector{BlasInt}(undef, 1)
            l = Vector{BlasInt}(undef, 1)
            lda = max(1,stride(A, 2))
            ldb = max(1,stride(B, 2))
            alpha = similar(A, $relty, n)
            beta = similar(A, $relty, n)
            ldu = max(1, m)
            U = jobu == 'U' ? similar(A, $elty, ldu, m) : similar(A, $elty, 0)
            ldv = max(1, p)
            V = jobv == 'V' ? similar(A, $elty, ldv, p) : similar(A, $elty, 0)
            ldq = max(1, n)
            Q = jobq == 'Q' ? similar(A, $elty, ldq, n) : similar(A, $elty, 0)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, 2n)
            iwork = Vector{BlasInt}(undef, n)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($f), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                    Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$relty}, Ptr{$relty}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$relty}, Ptr{BlasInt},
                    Ref{BlasInt}, Clong, Clong, Clong),
                    jobu, jobv, jobq, m,
                    n, p, k, l,
                    A, lda, B, ldb,
                    alpha, beta, U, ldu,
                    V, ldv, Q, ldq,
                    work, lwork, rwork, iwork,
                    info, 1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                end
            end
            if m - k[1] - l[1] >= 0
                R = triu(A[1:k[1] + l[1],n - k[1] - l[1] + 1:n])
            else
                R = triu([A[1:m, n - k[1] - l[1] + 1:n]; B[m - k[1] + 1:l[1], n - k[1] - l[1] + 1:n]])
            end
            return U, V, Q, alpha, beta, k[1], l[1], R
        end
    end
end

"""
    ggsvd3!(jobu, jobv, jobq, A, B) -> (U, V, Q, alpha, beta, k, l, R)

Finds the generalized singular value decomposition of `A` and `B`, `U'*A*Q = D1*R`
and `V'*B*Q = D2*R`. `D1` has `alpha` on its diagonal and `D2` has `beta` on its
diagonal. If `jobu = U`, the orthogonal/unitary matrix `U` is computed. If
`jobv = V` the orthogonal/unitary matrix `V` is computed. If `jobq = Q`,
the orthogonal/unitary matrix `Q` is computed. If `jobu`, `jobv`, or `jobq` is
`N`, that matrix is not computed. This function requires LAPACK 3.6.0.
"""
ggsvd3!

## Expert driver and generalized eigenvalue problem
for (geevx, ggev, ggev3, elty) in
    ((:dgeevx_,:dggev_,:dggev3_,:Float64),
     (:sgeevx_,:sggev_,:sggev3_,:Float32))
    @eval begin
        #     SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,
        #                          VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,
        #                          RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
        #
        #       .. Scalar Arguments ..
        #       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
        #       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
        #       DOUBLE PRECISION   ABNRM
        #       ..
        #       .. Array Arguments ..
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   A( LDA, * ), RCONDE( * ), RCONDV( * ),
        #      $                   SCALE( * ), VL( LDVL, * ), VR( LDVR, * ),
        #      $                   WI( * ), WORK( * ), WR( * )
        function geevx!(balanc::AbstractChar, jobvl::AbstractChar, jobvr::AbstractChar, sense::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            @chkvalidparam 1 balanc ('N', 'P', 'S', 'B')
            @chkvalidparam 4 sense ('N', 'E', 'V', 'B')
            if sense ∈ ('E', 'B') && !(jobvl == jobvr == 'V')
                throw(ArgumentError(lazy"sense = '$sense' requires jobvl = 'V' and jobvr = 'V'"))
            end
            n = checksquare(A)
            ldvl = 0
            if jobvl == 'V'
                ldvl = n
            elseif jobvl == 'N'
                ldvl = 0
            else
                throw(ArgumentError(lazy"jobvl must be 'V' or 'N', but $jobvl was passed"))
            end
            ldvr = 0
            if jobvr == 'V'
                ldvr = n
            elseif jobvr == 'N'
                ldvr = 0
            else
                throw(ArgumentError(lazy"jobvr must be 'V' or 'N', but $jobvr was passed"))
            end
            chkfinite(A) # balancing routines don't support NaNs and Infs
            lda = max(1,stride(A,2))
            wr = similar(A, $elty, n)
            wi = similar(A, $elty, n)
            VL = similar(A, $elty, ldvl, n)
            VR = similar(A, $elty, ldvr, n)
            ilo = Ref{BlasInt}()
            ihi = Ref{BlasInt}()
            scale = similar(A, $elty, n)
            abnrm = Ref{$elty}()
            rconde = similar(A, $elty, n)
            rcondv = similar(A, $elty, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            iworksize = 0
            if sense == 'N' || sense == 'E'
                iworksize = 0
            elseif sense == 'V' || sense == 'B'
                iworksize = 2*n - 2
            else
                throw(ArgumentError(lazy"sense must be 'N', 'E', 'V' or 'B', but $sense was passed"))
            end
            iwork = Vector{BlasInt}(undef, iworksize)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($geevx), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                       Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                       Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                       Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty},
                       Ptr{$elty}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty},
                       Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                       Clong, Clong, Clong, Clong),
                       balanc, jobvl, jobvr, sense,
                       n, A, lda, wr,
                       wi, VL, max(1,ldvl), VR,
                       max(1,ldvr), ilo, ihi, scale,
                       abnrm, rconde, rcondv, work,
                       lwork, iwork, info,
                       1, 1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                end
            end
            A, wr, wi, VL, VR, ilo[], ihi[], scale, abnrm[], rconde, rcondv
        end

        #       SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
        #      $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVL, JOBVR
        #       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
        #      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
        #      $                   VR( LDVR, * ), WORK( * )
        function ggev!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A,B)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"A has dimensions $(size(A)), and B has dimensions $(size(B)), but A and B must have the same size"))
            end
            ldvl = 0
            if jobvl == 'V'
                ldvl = n
            elseif jobvl == 'N'
                ldvl = 1
            else
                throw(ArgumentError(lazy"jobvl must be 'V' or 'N', but $jobvl was passed"))
            end
            ldvr = 0
            if jobvr == 'V'
                ldvr = n
            elseif jobvr == 'N'
                ldvr = 1
            else
                throw(ArgumentError(lazy"jobvr must be 'V' or 'N', but $jobvr was passed"))
            end
            lda = max(1, stride(A, 2))
            ldb = max(1, stride(B, 2))
            alphar = similar(A, $elty, n)
            alphai = similar(A, $elty, n)
            beta = similar(A, $elty, n)
            vl = similar(A, $elty, ldvl, n)
            vr = similar(A, $elty, ldvr, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ggev), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ptr{$elty}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ref{BlasInt}, Clong, Clong),
                    jobvl, jobvr, n, A,
                    lda, B, ldb, alphar,
                    alphai, beta, vl, ldvl,
                    vr, ldvr, work, lwork,
                    info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                end
            end
            alphar, alphai, beta, vl, vr
        end

        #       SUBROUTINE DGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
        #      $                   BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVL, JOBVR
        #       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
        #      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
        #      $                   VR( LDVR, * ), WORK( * )
        function ggev3!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A,B)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"A has dimensions $(size(A)), and B has dimensions $(size(B)), but A and B must have the same size"))
            end
            ldvl = 0
            if jobvl == 'V'
                ldvl = n
            elseif jobvl == 'N'
                ldvl = 1
            else
                throw(ArgumentError(lazy"jobvl must be 'V' or 'N', but $jobvl was passed"))
            end
            ldvr = 0
            if jobvr == 'V'
                ldvr = n
            elseif jobvr == 'N'
                ldvr = 1
            else
                throw(ArgumentError(lazy"jobvr must be 'V' or 'N', but $jobvr was passed"))
            end
            lda = max(1, stride(A, 2))
            ldb = max(1, stride(B, 2))
            alphar = similar(A, $elty, n)
            alphai = similar(A, $elty, n)
            beta = similar(A, $elty, n)
            vl = similar(A, $elty, ldvl, n)
            vr = similar(A, $elty, ldvr, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ggev3), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ptr{$elty}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ref{BlasInt}, Clong, Clong),
                    jobvl, jobvr, n, A,
                    lda, B, ldb, alphar,
                    alphai, beta, vl, ldvl,
                    vr, ldvr, work, lwork,
                    info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                end
            end
            alphar, alphai, beta, vl, vr
        end
    end
end

for (geevx, ggev, ggev3, elty, relty) in
    ((:zgeevx_,:zggev_,:zggev3_,:ComplexF64,:Float64),
     (:cgeevx_,:cggev_,:cggev3_,:ComplexF32,:Float32))
    @eval begin
        #     SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL,
        #                          LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE,
        #                          RCONDV, WORK, LWORK, RWORK, INFO )
        #
        #       .. Scalar Arguments ..
        #       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
        #       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
        #       DOUBLE PRECISION   ABNRM
        #       ..
        #       .. Array Arguments ..
        #       DOUBLE PRECISION   RCONDE( * ), RCONDV( * ), RWORK( * ),
        #      $                   SCALE( * )
        #       COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
        #      $                   W( * ), WORK( * )
        function geevx!(balanc::AbstractChar, jobvl::AbstractChar, jobvr::AbstractChar, sense::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            if balanc ∉ ('N', 'P', 'S', 'B')
                throw(ArgumentError(lazy"balanc must be 'N', 'P', 'S', or 'B', but $balanc was passed"))
            end
            if sense ∉ ('N','E','V','B')
                throw(ArgumentError(lazy"sense must be 'N', 'E', 'V' or 'B', but $sense was passed"))
            end
            if sense ∈ ('E', 'B') && !(jobvl == jobvr == 'V')
                throw(ArgumentError(lazy"sense = '$sense' requires jobvl = 'V' and jobvr = 'V'"))
            end
            n = checksquare(A)
            ldvl = 0
            if jobvl == 'V'
                ldvl = n
            elseif jobvl == 'N'
                ldvl = 0
            else
                throw(ArgumentError(lazy"jobvl must be 'V' or 'N', but $jobvl was passed"))
            end
            ldvr = 0
            if jobvr == 'V'
                ldvr = n
            elseif jobvr == 'N'
                ldvr = 0
            else
                throw(ArgumentError(lazy"jobvr must be 'V' or 'N', but $jobvr was passed"))
            end
            chkfinite(A) # balancing routines don't support NaNs and Infs
            lda = max(1,stride(A,2))
            w = similar(A, $elty, n)
            VL = similar(A, $elty, ldvl, n)
            VR = similar(A, $elty, ldvr, n)
            ilo = Ref{BlasInt}()
            ihi = Ref{BlasInt}()
            scale = similar(A, $relty, n)
            abnrm = Ref{$relty}()
            rconde = similar(A, $relty, n)
            rcondv = similar(A, $relty, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, 2n)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($geevx), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                       Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$relty}, Ptr{$relty},
                       Ptr{$relty}, Ptr{$relty}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$relty}, Ref{BlasInt}, Clong, Clong, Clong, Clong),
                       balanc, jobvl, jobvr, sense,
                       n, A, lda, w,
                       VL, max(1,ldvl), VR, max(1,ldvr),
                       ilo, ihi, scale, abnrm,
                       rconde, rcondv, work, lwork,
                       rwork, info, 1, 1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                end
            end
            A, w, VL, VR, ilo[], ihi[], scale, abnrm[], rconde, rcondv
        end

        # SUBROUTINE ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
        #      $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVL, JOBVR
        #       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   RWORK( * )
        #       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
        #      $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
        #      $                   WORK( * )
        function ggev!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"A has dimensions $(size(A)), and B has dimensions $(size(B)), but A and B must have the same size"))
            end
            ldvl = 0
            if jobvl == 'V'
                ldvl = n
            elseif jobvl == 'N'
                ldvl = 1
            else
                throw(ArgumentError(lazy"jobvl must be 'V' or 'N', but $jobvl was passed"))
            end
            ldvr = 0
            if jobvr == 'V'
                ldvr = n
            elseif jobvr == 'N'
                ldvr = 1
            else
                throw(ArgumentError(lazy"jobvr must be 'V' or 'N', but $jobvr was passed"))
            end
            lda = max(1, stride(A, 2))
            ldb = max(1, stride(B, 2))
            alpha = similar(A, $elty, n)
            beta = similar(A, $elty, n)
            vl = similar(A, $elty, ldvl, n)
            vr = similar(A, $elty, ldvr, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, 8n)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ggev), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$relty},
                     Ref{BlasInt}, Clong, Clong),
                    jobvl, jobvr, n, A,
                    lda, B, ldb, alpha,
                    beta, vl, ldvl, vr,
                    ldvr, work, lwork, rwork,
                    info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                end
            end
            alpha, beta, vl, vr
        end

        # SUBROUTINE ZGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
        #      $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVL, JOBVR
        #       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   RWORK( * )
        #       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
        #      $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
        #      $                   WORK( * )
        function ggev3!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"A has dimensions $(size(A)), and B has dimensions $(size(B)), but A and B must have the same size"))
            end
            ldvl = 0
            if jobvl == 'V'
                ldvl = n
            elseif jobvl == 'N'
                ldvl = 1
            else
                throw(ArgumentError(lazy"jobvl must be 'V' or 'N', but $jobvl was passed"))
            end
            ldvr = 0
            if jobvr == 'V'
                ldvr = n
            elseif jobvr == 'N'
                ldvr = 1
            else
                throw(ArgumentError(lazy"jobvr must be 'V' or 'N', but $jobvr was passed"))
            end
            lda = max(1, stride(A, 2))
            ldb = max(1, stride(B, 2))
            alpha = similar(A, $elty, n)
            beta = similar(A, $elty, n)
            vl = similar(A, $elty, ldvl, n)
            vr = similar(A, $elty, ldvr, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, 8n)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ggev3), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$relty},
                     Ref{BlasInt}, Clong, Clong),
                    jobvl, jobvr, n, A,
                    lda, B, ldb, alpha,
                    beta, vl, ldvl, vr,
                    ldvr, work, lwork, rwork,
                    info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                end
            end
            alpha, beta, vl, vr
        end
    end
end

"""
    geevx!(balanc, jobvl, jobvr, sense, A) -> (A, w, VL, VR, ilo, ihi, scale, abnrm, rconde, rcondv)

Finds the eigensystem of `A` with matrix balancing. If `jobvl = N`, the
left eigenvectors of `A` aren't computed. If `jobvr = N`, the right
eigenvectors of `A` aren't computed. If `jobvl = V` or `jobvr = V`, the
corresponding eigenvectors are computed. If `balanc = N`, no balancing is
performed. If `balanc = P`, `A` is permuted but not scaled. If
`balanc = S`, `A` is scaled but not permuted. If `balanc = B`, `A` is
permuted and scaled. If `sense = N`, no reciprocal condition numbers are
computed. If `sense = E`, reciprocal condition numbers are computed for
the eigenvalues only. If `sense = V`, reciprocal condition numbers are
computed for the right eigenvectors only. If `sense = B`, reciprocal
condition numbers are computed for the right eigenvectors and the
eigenvectors. If `sense = E,B`, the right and left eigenvectors must be
computed.
"""
geevx!(balanc::AbstractChar, jobvl::AbstractChar, jobvr::AbstractChar, sense::AbstractChar, A::AbstractMatrix)

"""
    ggev!(jobvl, jobvr, A, B) -> (alpha, beta, vl, vr)

Finds the generalized eigendecomposition of `A` and `B`. If `jobvl = N`,
the left eigenvectors aren't computed. If `jobvr = N`, the right
eigenvectors aren't computed. If `jobvl = V` or `jobvr = V`, the
corresponding eigenvectors are computed.
"""
ggev!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix, B::AbstractMatrix)

"""
    ggev3!(jobvl, jobvr, A, B) -> (alpha, beta, vl, vr)

Finds the generalized eigendecomposition of `A` and `B` using a blocked
algorithm. If `jobvl = N`, the left eigenvectors aren't computed. If
`jobvr = N`, the right eigenvectors aren't computed. If `jobvl = V` or
`jobvr = V`, the corresponding eigenvectors are computed.  This function
requires LAPACK 3.6.0.
"""
ggev3!(jobvl::AbstractChar, jobvr::AbstractChar, A::AbstractMatrix, B::AbstractMatrix)

# One step incremental condition estimation of max/min singular values
for (laic1, elty) in
    ((:dlaic1_,:Float64),
     (:slaic1_,:Float32))
    @eval begin
        #  SUBROUTINE DLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )
        #
        #  .. Scalar Arguments ..
        #  INTEGER            J, JOB
        #  DOUBLE PRECISION   C, GAMMA, S, SEST, SESTPR
        #  ..
        #  .. Array Arguments ..
        #  DOUBLE PRECISION   W( J ), X( J )
        function laic1!(job::Integer, x::AbstractVector{$elty},
                        sest::$elty, w::AbstractVector{$elty}, gamma::$elty)
            require_one_based_indexing(x, w)
            @chkvalidparam 1 job (1,2)
            j = length(x)
            if j != length(w)
                throw(DimensionMismatch(lazy"vectors must have same length, but length of x is $j and length of w is $(length(w))"))
            end
            sestpr = Ref{$elty}()
            s = Ref{$elty}()
            c = Ref{$elty}()
            ccall((@blasfunc($laic1), libblastrampoline), Cvoid,
                (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{$elty},
                 Ptr{$elty}, Ref{$elty}, Ref{$elty}, Ref{$elty},
                 Ref{$elty}),
                job, j, x, sest,
                w, gamma, sestpr, s,
                c)
            sestpr[], s[], c[]
        end
    end
end
for (laic1, elty, relty) in
    ((:zlaic1_,:ComplexF64,:Float64),
     (:claic1_,:ComplexF32,:Float32))
    @eval begin
       #  SUBROUTINE ZLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )
       #
       #  .. Scalar Arguments ..
       #  INTEGER            J, JOB
       #  DOUBLE PRECISION   SEST, SESTPR
       #  COMPLEX*16         C, GAMMA, S
       #  ..
       #  .. Array Arguments ..
       #  COMPLEX*16         W( J ), X( J )
        function laic1!(job::Integer, x::AbstractVector{$elty},
                        sest::$relty, w::AbstractVector{$elty}, gamma::$elty)
            require_one_based_indexing(x, w)
            @chkvalidparam 1 job (1,2)
            j = length(x)
            if j != length(w)
                throw(DimensionMismatch(lazy"vectors must have same length, but length of x is $j and length of w is $(length(w))"))
            end
            sestpr = Ref{$relty}()
            s = Ref{$elty}()
            c = Ref{$elty}()
            ccall((@blasfunc($laic1), libblastrampoline), Cvoid,
                (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{$relty},
                 Ptr{$elty}, Ref{$elty}, Ref{$relty}, Ref{$elty},
                 Ref{$elty}),
                job, j, x, sest,
                w, gamma, sestpr, s,
                c)
            sestpr[], s[], c[]
        end
    end
end

# (GT) General tridiagonal, decomposition, solver and direct solver
for (gtsv, gttrf, gttrs, elty) in
    ((:dgtsv_,:dgttrf_,:dgttrs_,:Float64),
     (:sgtsv_,:sgttrf_,:sgttrs_,:Float32),
     (:zgtsv_,:zgttrf_,:zgttrs_,:ComplexF64),
     (:cgtsv_,:cgttrf_,:cgttrs_,:ComplexF32))
    @eval begin
        #       SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
        #       .. Scalar Arguments ..
        #       INTEGER            INFO, LDB, N, NRHS
        #       .. Array Arguments ..
        #       DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
        function gtsv!(dl::AbstractVector{$elty}, d::AbstractVector{$elty}, du::AbstractVector{$elty},
                       B::AbstractVecOrMat{$elty})
            require_one_based_indexing(dl, d, du, B)
            chkstride1(B, dl, d, du)
            n = length(d)
            if !(n >= length(dl) >= n - 1)
                throw(DimensionMismatch(lazy"subdiagonal has length $(length(dl)), but should be $n or $(n - 1)"))
            end
            if !(n >= length(du) >= n - 1)
                throw(DimensionMismatch(lazy"superdiagonal has length $(length(du)), but should be $n or $(n - 1)"))
            end
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has leading dimension $(size(B,1)), but should have $n"))
            end
            if n == 0
                return B # Early exit if possible
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($gtsv), libblastrampoline), Cvoid,
                  (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                  n, size(B,2), dl, d, du, B, max(1,stride(B,2)), info)
            chklapackerror(info[])
            B
        end

        #       SUBROUTINE DGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
        #       .. Scalar Arguments ..
        #       INTEGER            INFO, N
        #       .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   D( * ), DL( * ), DU( * ), DU2( * )
        function gttrf!(dl::AbstractVector{$elty}, d::AbstractVector{$elty}, du::AbstractVector{$elty})
            require_one_based_indexing(dl, d, du)
            chkstride1(dl,d,du)
            n    = length(d)
            if length(dl) != n - 1
                throw(DimensionMismatch(lazy"subdiagonal has length $(length(dl)), but should be $(n - 1)"))
            end
            if length(du) != n - 1
                throw(DimensionMismatch(lazy"superdiagonal has length $(length(du)), but should be $(n - 1)"))
            end
            du2  = similar(d, $elty, n-2)
            ipiv = similar(d, BlasInt, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($gttrf), libblastrampoline), Cvoid,
                  (Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty},
                   Ptr{BlasInt}, Ref{BlasInt}),
                  n, dl, d, du, du2, ipiv, info)
            chklapackerror(info[])
            dl, d, du, du2, ipiv
        end

        #       SUBROUTINE DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )
        #       .. Scalar Arguments ..
        #       CHARACTER          TRANS
        #       INTEGER            INFO, LDB, N, NRHS
        #       .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
        function gttrs!(trans::AbstractChar, dl::AbstractVector{$elty}, d::AbstractVector{$elty},
                        du::AbstractVector{$elty}, du2::AbstractVector{$elty}, ipiv::AbstractVector{BlasInt},
                        B::AbstractVecOrMat{$elty})
            require_one_based_indexing(dl, d, du, du2, ipiv, B)
            chktrans(trans)
            chkstride1(B, ipiv, dl, d, du, du2)
            n = length(d)
            if length(dl) != n - 1
                throw(DimensionMismatch(lazy"subdiagonal has length $(length(dl)), but should be $(n - 1)"))
            end
            if length(du) != n - 1
                throw(DimensionMismatch(lazy"superdiagonal has length $(length(du)), but should be $(n - 1)"))
            end
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has leading dimension $(size(B,1)), but should have $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($gttrs), libblastrampoline), Cvoid,
                   (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                    Ptr{$elty}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty},
                    Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                   trans, n, size(B,2), dl, d, du, du2, ipiv, B, max(1,stride(B,2)), info, 1)
            chklapackerror(info[])
            B
         end
    end
end

"""
    gtsv!(dl, d, du, B)

Solves the equation `A * X = B` where `A` is a tridiagonal matrix with
`dl` on the subdiagonal, `d` on the diagonal, and `du` on the
superdiagonal.

Overwrites `B` with the solution `X` and returns it.
"""
gtsv!(dl::AbstractVector, d::AbstractVector, du::AbstractVector, B::AbstractVecOrMat)

"""
    gttrf!(dl, d, du) -> (dl, d, du, du2, ipiv)

Finds the `LU` factorization of a tridiagonal matrix with `dl` on the
subdiagonal, `d` on the diagonal, and `du` on the superdiagonal.

Modifies `dl`, `d`, and `du` in-place and returns them and the second
superdiagonal `du2` and the pivoting vector `ipiv`.
"""
gttrf!(dl::AbstractVector, d::AbstractVector, du::AbstractVector)

"""
    gttrs!(trans, dl, d, du, du2, ipiv, B)

Solves the equation `A * X = B` (`trans = N`), `transpose(A) * X = B` (`trans = T`),
or `adjoint(A) * X = B` (`trans = C`) using the `LU` factorization computed by
`gttrf!`. `B` is overwritten with the solution `X`.
"""
gttrs!(trans::AbstractChar, dl::AbstractVector, d::AbstractVector, du::AbstractVector, du2::AbstractVector,
       ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat)

## (OR) orthogonal (or UN, unitary) matrices, extractors and multiplication
for (orglq, orgqr, orgql, orgrq, ormlq, ormqr, ormql, ormrq, gemqrt, elty) in
    ((:dorglq_,:dorgqr_,:dorgql_,:dorgrq_,:dormlq_,:dormqr_,:dormql_,:dormrq_,:dgemqrt_,:Float64),
     (:sorglq_,:sorgqr_,:sorgql_,:sorgrq_,:sormlq_,:sormqr_,:sormql_,:sormrq_,:sgemqrt_,:Float32),
     (:zunglq_,:zungqr_,:zungql_,:zungrq_,:zunmlq_,:zunmqr_,:zunmql_,:zunmrq_,:zgemqrt_,:ComplexF64),
     (:cunglq_,:cungqr_,:cungql_,:cungrq_,:cunmlq_,:cunmqr_,:cunmql_,:cunmrq_,:cgemqrt_,:ComplexF32))
    @eval begin
        # SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, K, LDA, LWORK, M, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function orglq!(A::AbstractMatrix{$elty}, tau::AbstractVector{$elty}, k::Integer = length(tau))
            require_one_based_indexing(A, tau)
            chkstride1(A,tau)
            n = size(A, 2)
            m = min(n, size(A, 1))
            if k > m
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= m = $m"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($orglq), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                       Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      m, n, k, A, max(1,stride(A,2)), tau, work, lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            if m < size(A,1)
                A[1:m,:]
            else
                A
            end
        end

        # SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, K, LDA, LWORK, M, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function orgqr!(A::AbstractMatrix{$elty}, tau::AbstractVector{$elty}, k::Integer = length(tau))
            require_one_based_indexing(A, tau)
            chkstride1(A,tau)
            m = size(A, 1)
            n = min(m, size(A, 2))
            if k > n
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= n = $n"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($orgqr), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                       Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      m, n, k, A,
                      max(1,stride(A,2)), tau, work, lwork,
                      info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            if n < size(A,2)
                reshape(resize!(vec(A), m * n), m, n)
            else
                A
            end
        end

        # SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, K, LDA, LWORK, M, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function orgql!(A::AbstractMatrix{$elty}, tau::AbstractVector{$elty}, k::Integer = length(tau))
            require_one_based_indexing(A, tau)
            chkstride1(A,tau)
            m = size(A, 1)
            n = min(m, size(A, 2))
            if k > n
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= n = $n"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($orgql), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                       Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      m, n, k, A,
                      max(1,stride(A,2)), tau, work, lwork,
                      info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            if n < size(A,2)
                reshape(resize!(vec(A), m * n), m, n)
            else
                A
            end
        end

        # SUBROUTINE DORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       INTEGER            INFO, K, LDA, LWORK, M, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function orgrq!(A::AbstractMatrix{$elty}, tau::AbstractVector{$elty}, k::Integer = length(tau))
            require_one_based_indexing(A, tau)
            chkstride1(A,tau)
            m, n = size(A)
            if n < m
                throw(DimensionMismatch(lazy"input matrix A has dimensions ($m,$n), but cannot have fewer columns than rows"))
            end
            if k > n
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= n = $n"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($orgrq), libblastrampoline), Cvoid,
                      (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                       Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                      m, n, k, A,
                      max(1,stride(A,2)), tau, work, lwork,
                      info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A
        end

        #      SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
        #                         WORK, LWORK, INFO )
        #      .. Scalar Arguments ..
        #      CHARACTER          SIDE, TRANS
        #      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
        #      .. Array Arguments ..
        #      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
        function ormlq!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix{$elty},
                        tau::AbstractVector{$elty}, C::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, tau, C)
            chktrans(trans)
            chkside(side)
            chkstride1(A, C, tau)
            m,n = ndims(C) == 2 ? size(C) : (size(C, 1), 1)
            nA = size(A, 2)
            k   = length(tau)
            if side == 'L' && m != nA
                throw(DimensionMismatch(lazy"for a left-sided multiplication, the first dimension of C, $m, must equal the second dimension of A, $nA"))
            end
            if side == 'R' && n != nA
                throw(DimensionMismatch(lazy"for a right-sided multiplication, the second dimension of C, $n, must equal the second dimension of A, $nA"))
            end
            if side == 'L' && k > m
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= m = $m"))
            end
            if side == 'R' && k > n
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= n = $n"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ormlq), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                      side, trans, m, n, k, A, max(1,stride(A,2)), tau,
                      C, max(1,stride(C,2)), work, lwork, info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            C
        end

        #      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
        #                         WORK, INFO )
        #      .. Scalar Arguments ..
        #      CHARACTER          SIDE, TRANS
        #      INTEGER            INFO, K, LDA, LDC, M, N
        #      .. Array Arguments ..
        #      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
        function ormqr!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix{$elty},
                        tau::AbstractVector{$elty}, C::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, tau, C)
            chktrans(trans)
            chkside(side)
            chkstride1(A, C, tau)
            m,n = ndims(C) == 2 ? size(C) : (size(C, 1), 1)
            mA  = size(A, 1)
            k   = length(tau)
            if side == 'L' && m != mA
                throw(DimensionMismatch(lazy"for a left-sided multiplication, the first dimension of C, $m, must equal the first dimension of A, $mA"))
            end
            if side == 'R' && n != mA
                throw(DimensionMismatch(lazy"for a right-sided multiplication, the second dimension of C, $n, must equal the first dimension of A, $mA"))
            end
            if side == 'L' && k > m
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= m = $m"))
            end
            if side == 'R' && k > n
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= n = $n"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ormqr), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                       Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ref{BlasInt}, Clong, Clong),
                      side, trans, m, n,
                      k, A, max(1,stride(A,2)), tau,
                      C, max(1, stride(C,2)), work, lwork,
                      info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            C
        end

        #      SUBROUTINE DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
        #                         WORK, INFO )
        #      .. Scalar Arguments ..
        #      CHARACTER          SIDE, TRANS
        #      INTEGER            INFO, K, LDA, LDC, M, N
        #      .. Array Arguments ..
        #      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
        function ormql!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix{$elty},
                        tau::AbstractVector{$elty}, C::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, tau, C)
            chktrans(trans)
            chkside(side)
            chkstride1(A, C, tau)
            m,n = ndims(C) == 2 ? size(C) : (size(C, 1), 1)
            mA  = size(A, 1)
            k   = length(tau)
            if side == 'L' && m != mA
                throw(DimensionMismatch(lazy"for a left-sided multiplication, the first dimension of C, $m, must equal the first dimension of A, $mA"))
            end
            if side == 'R' && n != mA
                throw(DimensionMismatch(lazy"for a right-sided multiplication, the second dimension of C, $n, must equal the first dimension of A, $mA"))
            end
            if side == 'L' && k > m
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= m = $m"))
            end
            if side == 'R' && k > n
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= n = $n"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ormql), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                       Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ref{BlasInt}, Clong, Clong),
                      side, trans, m, n,
                      k, A, max(1,stride(A,2)), tau,
                      C, max(1, stride(C,2)), work, lwork,
                      info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            C
        end

        #      SUBROUTINE DORMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
        #                         WORK, LWORK, INFO )
        #      .. Scalar Arguments ..
        #      CHARACTER          SIDE, TRANS
        #      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
        #      .. Array Arguments ..
        #      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
        function ormrq!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix{$elty},
                        tau::AbstractVector{$elty}, C::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, tau, C)
            chktrans(trans)
            chkside(side)
            chkstride1(A, C, tau)
            m,n = ndims(C) == 2 ? size(C) : (size(C, 1), 1)
            nA  = size(A, 2)
            k   = length(tau)
            if side == 'L' && m != nA
                throw(DimensionMismatch(lazy"for a left-sided multiplication, the first dimension of C, $m, must equal the second dimension of A, $nA"))
            end
            if side == 'R' && n != nA
                throw(DimensionMismatch(lazy"for a right-sided multiplication, the second dimension of C, $n, must equal the second dimension of A, $nA"))
            end
            if side == 'L' && k > m
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= m = $m"))
            end
            if side == 'R' && k > n
                throw(DimensionMismatch(lazy"invalid number of reflectors: k = $k should be <= n = $n"))
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ormrq), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                      side, trans, m, n, k, A, max(1,stride(A,2)), tau,
                      C, max(1,stride(C,2)), work, lwork, info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            C
        end

        function gemqrt!(side::AbstractChar, trans::AbstractChar, V::AbstractMatrix{$elty}, T::AbstractMatrix{$elty}, C::AbstractVecOrMat{$elty})
            require_one_based_indexing(V, T, C)
            chktrans(trans)
            chkside(side)
            chkstride1(V, T, C)
            m,n = ndims(C) == 2 ? size(C) : (size(C, 1), 1)
            nb, k = size(T)
            if k == 0
                return C
            end
            if side == 'L'
                if !(0 <= k <= m)
                    throw(DimensionMismatch(lazy"wrong value for k = $k: must be between 0 and $m"))
                end
                if m != size(V,1)
                    throw(DimensionMismatch(lazy"first dimensions of C, $m, and V, $(size(V,1)) must match"))
                end
                ldv = stride(V,2)
                if ldv < max(1, m)
                    throw(DimensionMismatch(lazy"Q and C don't fit! The stride of V, $ldv, is too small"))
                end
                wss = n*k
            elseif side == 'R'
                if !(0 <= k <= n)
                    throw(DimensionMismatch(lazy"wrong value for k = $k: must be between 0 and $n"))
                end
                if n != size(V,1)
                    throw(DimensionMismatch(lazy"second dimension of C, $n, and first dimension of V, $(size(V,1)) must match"))
                end
                ldv = stride(V,2)
                if ldv < max(1, n)
                    throw(DimensionMismatch(lazy"Q and C don't fit! The stride of V, $ldv, is too small"))
                end
                wss = m*k
            end
            if !(1 <= nb <= k)
                throw(DimensionMismatch(lazy"wrong value for nb = $nb, which must be between 1 and $k"))
            end
            ldc = stride(C, 2)
            work = Vector{$elty}(undef, wss)
            info = Ref{BlasInt}()
            ccall((@blasfunc($gemqrt), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                 Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                 Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                 Ptr{$elty}, Ref{BlasInt}, Clong, Clong),
                side, trans, m, n,
                k, nb, V, ldv,
                T, max(1,stride(T,2)), C, max(1,ldc),
                work, info, 1, 1)
            chklapackerror(info[])
            return C
        end
    end
end

"""
    orglq!(A, tau, k = length(tau))

Explicitly finds the matrix `Q` of a `LQ` factorization after calling
`gelqf!` on `A`. Uses the output of `gelqf!`. `A` is overwritten by `Q`.
"""
orglq!(A::AbstractMatrix, tau::AbstractVector, k::Integer = length(tau))

"""
    orgqr!(A, tau, k = length(tau))

Explicitly finds the matrix `Q` of a `QR` factorization after calling
`geqrf!` on `A`. Uses the output of `geqrf!`. `A` is overwritten by `Q`.
"""
orgqr!(A::AbstractMatrix, tau::AbstractVector, k::Integer = length(tau))

"""
    orgql!(A, tau, k = length(tau))

Explicitly finds the matrix `Q` of a `QL` factorization after calling
`geqlf!` on `A`. Uses the output of `geqlf!`. `A` is overwritten by `Q`.
"""
orgql!(A::AbstractMatrix, tau::AbstractVector, k::Integer = length(tau))

"""
    orgrq!(A, tau, k = length(tau))

Explicitly finds the matrix `Q` of a `RQ` factorization after calling
`gerqf!` on `A`. Uses the output of `gerqf!`. `A` is overwritten by `Q`.
"""
orgrq!(A::AbstractMatrix, tau::AbstractVector, k::Integer = length(tau))

"""
    ormlq!(side, trans, A, tau, C)

Computes `Q * C` (`trans = N`), `transpose(Q) * C` (`trans = T`), `adjoint(Q) * C`
(`trans = C`) for `side = L` or the equivalent right-sided multiplication
for `side = R` using `Q` from a `LQ` factorization of `A` computed using
`gelqf!`. `C` is overwritten.
"""
ormlq!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix, tau::AbstractVector, C::AbstractVecOrMat)

"""
    ormqr!(side, trans, A, tau, C)

Computes `Q * C` (`trans = N`), `transpose(Q) * C` (`trans = T`), `adjoint(Q) * C`
(`trans = C`) for `side = L` or the equivalent right-sided multiplication
for `side = R` using `Q` from a `QR` factorization of `A` computed using
`geqrf!`. `C` is overwritten.
"""
ormqr!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix, tau::AbstractVector, C::AbstractVecOrMat)

"""
    ormql!(side, trans, A, tau, C)

Computes `Q * C` (`trans = N`), `transpose(Q) * C` (`trans = T`), `adjoint(Q) * C`
(`trans = C`) for `side = L` or the equivalent right-sided multiplication
for `side = R` using `Q` from a `QL` factorization of `A` computed using
`geqlf!`. `C` is overwritten.
"""
ormql!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix, tau::AbstractVector, C::AbstractVecOrMat)

"""
    ormrq!(side, trans, A, tau, C)

Computes `Q * C` (`trans = N`), `transpose(Q) * C` (`trans = T`), `adjoint(Q) * C`
(`trans = C`) for `side = L` or the equivalent right-sided multiplication
for `side = R` using `Q` from a `RQ` factorization of `A` computed using
`gerqf!`. `C` is overwritten.
"""
ormrq!(side::AbstractChar, trans::AbstractChar, A::AbstractMatrix, tau::AbstractVector, C::AbstractVecOrMat)

"""
    gemqrt!(side, trans, V, T, C)

Computes `Q * C` (`trans = N`), `transpose(Q) * C` (`trans = T`), `adjoint(Q) * C`
(`trans = C`) for `side = L` or the equivalent right-sided multiplication
for `side = R` using `Q` from a `QR` factorization of `A` computed using
`geqrt!`. `C` is overwritten.
"""
gemqrt!(side::AbstractChar, trans::AbstractChar, V::AbstractMatrix, T::AbstractMatrix, C::AbstractVecOrMat)

# (PO) positive-definite symmetric matrices,
for (posv, potrf, potri, potrs, pstrf, elty, rtyp) in
    ((:dposv_,:dpotrf_,:dpotri_,:dpotrs_,:dpstrf_,:Float64,:Float64),
     (:sposv_,:spotrf_,:spotri_,:spotrs_,:spstrf_,:Float32,:Float32),
     (:zposv_,:zpotrf_,:zpotri_,:zpotrs_,:zpstrf_,:ComplexF64,:Float64),
     (:cposv_,:cpotrf_,:cpotri_,:cpotrs_,:cpstrf_,:ComplexF32,:Float32))
    @eval begin
        #     SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          UPLO
        #      INTEGER            INFO, LDA, LDB, N, NRHS
        #     .. Array Arguments ..
        #      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function posv!(uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            n = checksquare(A)
            chkuplo(uplo)
            if size(B,1) != n
                throw(DimensionMismatch(lazy"first dimension of B, $(size(B,1)), and size of A, ($n,$n), must match!"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($posv), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, n, size(B,2), A, max(1,stride(A,2)), B, max(1,stride(B,2)), info, 1)
            chkargsok(info[])
            chkposdef(info[])
            A, B
        end

        # SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * )
        function potrf!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            checksquare(A)
            chkuplo(uplo)
            lda = max(1,stride(A,2))
            if lda == 0
                return A, 0
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($potrf), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, size(A,1), A, lda, info, 1)
            chkargsok(info[])
            #info[] > 0 means the leading minor of order info[] is not positive definite
            #ordinarily, throw Exception here, but return error code here
            #this simplifies isposdef! and factorize
            return A, info[] # info stored in Cholesky
        end

        #       SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
        #       .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, N
        #       .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * )
        function potri!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            chkuplo(uplo)
            info = Ref{BlasInt}()
            ccall((@blasfunc($potri), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, size(A,1), A, max(1,stride(A,2)), info, 1)
            chkargsok(info[])
            chknonsingular(info[])
            A
        end

        #     SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
        #     .. Scalar Arguments ..
        #      CHARACTER          UPLO
        #      INTEGER            INFO, LDA, LDB, N, NRHS
        #     .. Array Arguments ..
        #      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function potrs!(uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            n = checksquare(A)
            chkuplo(uplo)
            nrhs = size(B,2)
            if size(B,1) != n
                throw(DimensionMismatch(lazy"first dimension of B, $(size(B,1)), and size of A, ($n,$n), must match!"))
            end
            lda = max(1,stride(A,2))
            if lda == 0 || nrhs == 0
                return B
            end
            ldb = max(1,stride(B,2))
            info = Ref{BlasInt}()
            ccall((@blasfunc($potrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                    Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                   uplo, n, nrhs, A,
                   lda, B, ldb, info, 1)
            chklapackerror(info[])
            return B
        end

        #       SUBROUTINE DPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )
        #       .. Scalar Arguments ..
        #       DOUBLE PRECISION   TOL
        #       INTEGER            INFO, LDA, N, RANK
        #       CHARACTER          UPLO
        #       .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), WORK( 2*N )
        #       INTEGER            PIV( N )
        function pstrf!(uplo::AbstractChar, A::AbstractMatrix{$elty}, tol::Real)
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkuplo(uplo)
            piv  = similar(A, BlasInt, n)
            rank = Vector{BlasInt}(undef, 1)
            work = Vector{$rtyp}(undef, 2n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($pstrf), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                   Ptr{BlasInt}, Ref{$rtyp}, Ptr{$rtyp}, Ref{BlasInt}, Clong),
                  uplo, n, A, max(1,stride(A,2)), piv, rank, tol, work, info, 1)
            chkargsok(info[])
            A, piv, rank[1], info[] #Stored in CholeskyPivoted
        end
    end
end

"""
    posv!(uplo, A, B) -> (A, B)

Finds the solution to `A * X = B` where `A` is a symmetric or Hermitian
positive definite matrix. If `uplo = U` the upper Cholesky decomposition
of `A` is computed. If `uplo = L` the lower Cholesky decomposition of `A`
is computed. `A` is overwritten by its Cholesky decomposition. `B` is
overwritten with the solution `X`.
"""
posv!(uplo::AbstractChar, A::AbstractMatrix, B::AbstractVecOrMat)

"""
    potrf!(uplo, A)

Computes the Cholesky (upper if `uplo = U`, lower if `uplo = L`)
decomposition of positive-definite matrix `A`. `A` is overwritten and
returned with an info code.
"""
potrf!(uplo::AbstractChar, A::AbstractMatrix)

"""
    potri!(uplo, A)

Computes the inverse of positive-definite matrix `A` after calling
`potrf!` to find its (upper if `uplo = U`, lower if `uplo = L`) Cholesky
decomposition.

`A` is overwritten by its inverse and returned.
"""
potri!(uplo::AbstractChar, A::AbstractMatrix)

"""
    potrs!(uplo, A, B)

Finds the solution to `A * X = B` where `A` is a symmetric or Hermitian
positive definite matrix whose Cholesky decomposition was computed by
`potrf!`. If `uplo = U` the upper Cholesky decomposition of `A` was
computed. If `uplo = L` the lower Cholesky decomposition of `A` was
computed. `B` is overwritten with the solution `X`.
"""
potrs!(uplo::AbstractChar, A::AbstractMatrix, B::AbstractVecOrMat)

"""
    pstrf!(uplo, A, tol) -> (A, piv, rank, info)

Computes the (upper if `uplo = U`, lower if `uplo = L`) pivoted Cholesky
decomposition of positive-definite matrix `A` with a user-set tolerance
`tol`. `A` is overwritten by its Cholesky decomposition.

Returns `A`, the pivots `piv`, the rank of `A`, and an `info` code. If `info = 0`,
the factorization succeeded. If `info = i > 0 `, then `A` is indefinite or
rank-deficient.
"""
pstrf!(uplo::AbstractChar, A::AbstractMatrix, tol::Real)

# (PT) positive-definite, symmetric, tri-diagonal matrices
# Direct solvers for general tridiagonal and symmetric positive-definite tridiagonal
for (ptsv, pttrf, elty, relty) in
    ((:dptsv_,:dpttrf_,:Float64,:Float64),
     (:sptsv_,:spttrf_,:Float32,:Float32),
     (:zptsv_,:zpttrf_,:ComplexF64,:Float64),
     (:cptsv_,:cpttrf_,:ComplexF32,:Float32))
    @eval begin
        #       SUBROUTINE DPTSV( N, NRHS, D, E, B, LDB, INFO )
        #       .. Scalar Arguments ..
        #       INTEGER            INFO, LDB, N, NRHS
        #       .. Array Arguments ..
        #       DOUBLE PRECISION   B( LDB, * ), D( * ), E( * )
        function ptsv!(D::AbstractVector{$relty}, E::AbstractVector{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(D, E, B)
            chkstride1(B, D, E)
            n = length(D)
            if length(E) != n - 1
                throw(DimensionMismatch(lazy"E has length $(length(E)), but needs $(n - 1)"))
            end
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)) but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($ptsv), libblastrampoline), Cvoid,
                  (Ref{BlasInt}, Ref{BlasInt}, Ptr{$relty}, Ptr{$elty},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                  n, size(B,2), D, E, B, max(1,stride(B,2)), info)
            chklapackerror(info[])
            B
        end

        #       SUBROUTINE DPTTRF( N, D, E, INFO )
        #       .. Scalar Arguments ..
        #       INTEGER            INFO, N
        #       .. Array Arguments ..
        #       DOUBLE PRECISION   D( * ), E( * )
        function pttrf!(D::AbstractVector{$relty}, E::AbstractVector{$elty})
            require_one_based_indexing(D, E)
            chkstride1(D, E)
            n = length(D)
            if length(E) != n - 1
                throw(DimensionMismatch(lazy"E has length $(length(E)), but needs $(n - 1)"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($pttrf), libblastrampoline), Cvoid,
                  (Ref{BlasInt}, Ptr{$relty}, Ptr{$elty}, Ptr{BlasInt}),
                  n, D, E, info)
            chklapackerror(info[])
            D, E
        end
    end
end

"""
    ptsv!(D, E, B)

Solves `A * X = B` for positive-definite tridiagonal `A`. `D` is the
diagonal of `A` and `E` is the off-diagonal. `B` is overwritten with the
solution `X` and returned.
"""
ptsv!(D::AbstractVector, E::AbstractVector, B::AbstractVecOrMat)

"""
    pttrf!(D, E)

Computes the LDLt factorization of a positive-definite tridiagonal matrix
with `D` as diagonal and `E` as off-diagonal. `D` and `E` are overwritten
and returned.
"""
pttrf!(D::AbstractVector, E::AbstractVector)

for (pttrs, elty, relty) in
    ((:dpttrs_,:Float64,:Float64),
     (:spttrs_,:Float32,:Float32))
    @eval begin
        #       SUBROUTINE DPTTRS( N, NRHS, D, E, B, LDB, INFO )
        #       .. Scalar Arguments ..
        #       INTEGER            INFO, LDB, N, NRHS
        #       .. Array Arguments ..
        #       DOUBLE PRECISION   B( LDB, * ), D( * ), E( * )
        function pttrs!(D::AbstractVector{$relty}, E::AbstractVector{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(D, E, B)
            chkstride1(B, D, E)
            n = length(D)
            if length(E) != n - 1
                throw(DimensionMismatch(lazy"E has length $(length(E)), but needs $(n - 1)"))
            end
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)) but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($pttrs), libblastrampoline), Cvoid,
                  (Ref{BlasInt}, Ref{BlasInt}, Ptr{$relty}, Ptr{$elty},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}),
                  n, size(B,2), D, E, B, max(1,stride(B,2)), info)
            chklapackerror(info[])
            B
        end
    end
end

for (pttrs, elty, relty) in
    ((:zpttrs_,:ComplexF64,:Float64),
     (:cpttrs_,:ComplexF32,:Float32))
    @eval begin
        #       SUBROUTINE ZPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDB, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   D( * )
        #       COMPLEX*16         B( LDB, * ), E( * )
        function pttrs!(uplo::AbstractChar, D::AbstractVector{$relty}, E::AbstractVector{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(D, E, B)
            chkstride1(B, D, E)
            chkuplo(uplo)
            n = length(D)
            if length(E) != n - 1
                throw(DimensionMismatch(lazy"E has length $(length(E)), but needs $(n - 1)"))
            end
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)) but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($pttrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$relty}, Ptr{$elty},
                   Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, n, size(B,2), D, E, B, max(1,stride(B,2)), info, 1)
            chklapackerror(info[])
            B
        end
    end
end

"""
    pttrs!(D, E, B)

Solves `A * X = B` for positive-definite tridiagonal `A` with diagonal
`D` and off-diagonal `E` after computing `A`'s LDLt factorization using
`pttrf!`. `B` is overwritten with the solution `X`.
"""
pttrs!(D::AbstractVector, E::AbstractVector, B::AbstractVecOrMat)

## (TR) triangular matrices: solver and inverse
for (trtri, trtrs, elty) in
    ((:dtrtri_,:dtrtrs_,:Float64),
     (:strtri_,:strtrs_,:Float32),
     (:ztrtri_,:ztrtrs_,:ComplexF64),
     (:ctrtri_,:ctrtrs_,:ComplexF32))
    @eval begin
        #     SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
        #*     .. Scalar Arguments ..
        #      CHARACTER          DIAG, UPLO
        #      INTEGER            INFO, LDA, N
        #     .. Array Arguments ..
        #      DOUBLE PRECISION   A( LDA, * )
        function trtri!(uplo::AbstractChar, diag::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkuplo(uplo)
            chkdiag(diag)
            lda = max(1,stride(A, 2))
            info = Ref{BlasInt}()
            ccall((@blasfunc($trtri), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ref{BlasInt}, Clong, Clong),
                  uplo, diag, n, A, lda, info, 1, 1)
            chklapackerror(info[])
            A
        end

        #      SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          DIAG, TRANS, UPLO
        #       INTEGER            INFO, LDA, LDB, N, NRHS
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function trtrs!(uplo::AbstractChar, trans::AbstractChar, diag::AbstractChar,
                        A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chktrans(trans)
            chkdiag(diag)
            chkstride1(A)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)) but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($trtrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                   Clong, Clong, Clong),
                  uplo, trans, diag, n, size(B,2), A, max(1,stride(A,2)),
                  B, max(1,stride(B,2)), info,
                  1, 1, 1)
            chklapackerror(info[], trtrs!)
            B
        end
    end
end
chklapackerror_positive(ret, ::typeof(trtrs!)) = chknonsingular(ret)

"""
    trtri!(uplo, diag, A)

Finds the inverse of (upper if `uplo = U`, lower if `uplo = L`)
triangular matrix `A`. If `diag = N`, `A` has non-unit diagonal elements.
If `diag = U`, all diagonal elements of `A` are one. `A` is overwritten
with its inverse.
"""
trtri!(uplo::AbstractChar, diag::AbstractChar, A::AbstractMatrix)

"""
    trtrs!(uplo, trans, diag, A, B)

Solves `A * X = B` (`trans = N`), `transpose(A) * X = B` (`trans = T`), or
`adjoint(A) * X = B` (`trans = C`) for (upper if `uplo = U`, lower if `uplo = L`)
triangular matrix `A`. If `diag = N`, `A` has non-unit diagonal elements.
If `diag = U`, all diagonal elements of `A` are one. `B` is overwritten
with the solution `X`.
"""
trtrs!(uplo::AbstractChar, trans::AbstractChar, diag::AbstractChar, A::AbstractMatrix, B::AbstractVecOrMat)

#Eigenvector computation and condition number estimation
for (trcon, trevc, trrfs, elty) in
    ((:dtrcon_,:dtrevc_,:dtrrfs_,:Float64),
     (:strcon_,:strevc_,:strrfs_,:Float32))
    @eval begin
        # SUBROUTINE DTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK,
        #                  IWORK, INFO )
        # .. Scalar Arguments ..
        # CHARACTER          DIAG, NORM, UPLO
        # INTEGER            INFO, LDA, N
        # DOUBLE PRECISION   RCOND
        # .. Array Arguments ..
        # INTEGER            IWORK( * )
        # DOUBLE PRECISION   A( LDA, * ), WORK( * )
        function trcon!(norm::AbstractChar, uplo::AbstractChar, diag::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            chkdiag(diag)
            n = checksquare(A)
            chkuplo(uplo)
            @chkvalidparam 1 norm ('O', '1', 'I')
            rcond = Ref{$elty}()
            work  = Vector{$elty}(undef, 3n)
            iwork = Vector{BlasInt}(undef, n)
            info  = Ref{BlasInt}()
            ccall((@blasfunc($trcon), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ref{$elty}, Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt},
                   Clong, Clong, Clong),
                  norm, uplo, diag, n,
                  A, max(1,stride(A,2)), rcond, work, iwork, info,
                  1, 1, 1)
            chklapackerror(info[])
            rcond[]
        end

        # SUBROUTINE DTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
        #                    LDVR, MM, M, WORK, INFO )
        #
        # .. Scalar Arguments ..
        # CHARACTER          HOWMNY, SIDE
        # INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
        # ..
        # .. Array Arguments ..
        # LOGICAL            SELECT( * )
        # DOUBLE PRECISION   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
        #$                   WORK( * )
        Base.@constprop :aggressive function trevc!(side::AbstractChar, howmny::AbstractChar, select::AbstractVector{BlasInt}, T::AbstractMatrix{$elty},
                        VL::AbstractMatrix{$elty} = similar(T),
                        VR::AbstractMatrix{$elty} = similar(T))
            require_one_based_indexing(select, T, VL, VR)
            # Extract
            if side ∉ ('L','R','B')
                throw(ArgumentError(lazy"side argument must be 'L' (left eigenvectors), 'R' (right eigenvectors), or 'B' (both), got $side"))
            end
            @chkvalidparam 2 howmny ('A', 'B', 'S')
            n, mm = checksquare(T), size(VL, 2)
            ldt, ldvl, ldvr = stride(T, 2), stride(VL, 2), stride(VR, 2)

            # Check
            chkstride1(T, select, VL, VR)

            # Allocate
            m = Ref{BlasInt}()
            work = Vector{$elty}(undef, 3n)
            info = Ref{BlasInt}()

            ccall((@blasfunc($trevc), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt},
                 Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                 Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                 Ptr{$elty}, Ref{BlasInt}, Clong, Clong),
                side, howmny, select, n,
                T, ldt, VL, ldvl,
                VR, ldvr, mm, m,
                work, info, 1, 1)
            chklapackerror(info[])

            VLn = size(VL, 1)
            VRn = size(VR, 1)
            #Decide what exactly to return
            if howmny == 'S' #compute selected eigenvectors
                if side == 'L' #left eigenvectors only
                    return select, reshape(resize!(vec(VL), VLn * m[]), VLn, m[])
                elseif side == 'R' #right eigenvectors only
                    return select, reshape(resize!(vec(VR), VRn * m[]), VRn, m[])
                else #side == 'B' #both eigenvectors
                    return select, reshape(resize!(vec(VL), VLn * m[]), VLn, m[]), reshape(resize!(vec(VR), VRn * m[]), VRn, m[])
                end
            else #compute all eigenvectors
                if side == 'L' #left eigenvectors only
                    return reshape(resize!(vec(VL), VLn * m[]), VLn, m[])
                elseif side == 'R' #right eigenvectors only
                    return reshape(resize!(vec(VR), VRn * m[]), VRn, m[])
                else #side == 'B' #both eigenvectors
                    return reshape(resize!(vec(VL), VLn * m[]), VLn, m[]), reshape(resize!(vec(VR), VRn * m[]), VRn, m[])
                end
            end
        end

        # SUBROUTINE DTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X,
        #                    LDX, FERR, BERR, WORK, IWORK, INFO )
        # .. Scalar Arguments ..
        # CHARACTER          DIAG, TRANS, UPLO
        # INTEGER            INFO, LDA, LDB, LDX, N, NRHS
        # .. Array Arguments ..
        # INTEGER            IWORK( * )
        # DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), BERR( * ), FERR( * ),
        #$                   WORK( * ), X( LDX, * )
        function trrfs!(uplo::AbstractChar, trans::AbstractChar, diag::AbstractChar,
                A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty}, X::AbstractVecOrMat{$elty},
                Ferr::AbstractVector{$elty} = similar(B, $elty, size(B,2)),
                Berr::AbstractVector{$elty} = similar(B, $elty, size(B,2)))
            require_one_based_indexing(A, B, X, Ferr, Berr)
            chkstride1(A, B, X, Ferr, Berr)
            chktrans(trans)
            chkuplo(uplo)
            chkdiag(diag)
            n = size(A,2)
            nrhs = size(B,2)
            if nrhs != size(X,2)
                throw(DimensionMismatch(lazy"second dimensions of B, $nrhs, and X, $(size(X,2)), must match"))
            end
            work = Vector{$elty}(undef, 3n)
            iwork = Vector{BlasInt}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($trrfs), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                 Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                 Ptr{$elty}, Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                uplo, trans, diag, n,
                nrhs, A, max(1,stride(A,2)), B, max(1,stride(B,2)), X, max(1,stride(X,2)),
                Ferr, Berr, work, iwork, info, 1, 1, 1)
            chklapackerror(info[])
            Ferr, Berr
        end
    end
end

for (trcon, trevc, trrfs, elty, relty) in
    ((:ztrcon_,:ztrevc_,:ztrrfs_,:ComplexF64,:Float64),
     (:ctrcon_,:ctrevc_,:ctrrfs_,:ComplexF32, :Float32))
    @eval begin
        # SUBROUTINE ZTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK,
        #                   RWORK, INFO )
        # .. Scalar Arguments ..
        # CHARACTER          DIAG, NORM, UPLO
        # INTEGER            INFO, LDA, N
        # DOUBLE PRECISION   RCOND
        # .. Array Arguments ..
        # DOUBLE PRECISION   RWORK( * )
        # COMPLEX*16         A( LDA, * ), WORK( * )
        function trcon!(norm::AbstractChar, uplo::AbstractChar, diag::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            @chkvalidparam 1 norm ('O', '1', 'I')
            chkuplo(uplo)
            chkdiag(diag)
            rcond = Ref{$relty}(1)
            work  = Vector{$elty}(undef, 2n)
            rwork = Vector{$relty}(undef, n)
            info  = Ref{BlasInt}()
            ccall((@blasfunc($trcon), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ref{$relty}, Ptr{$elty}, Ptr{$relty}, Ptr{BlasInt},
                   Clong, Clong, Clong),
                  norm, uplo, diag, n,
                  A, max(1,stride(A,2)), rcond, work, rwork, info,
                  1, 1, 1)
            chklapackerror(info[])
            rcond[]
        end

        # SUBROUTINE ZTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
        #                    LDVR, MM, M, WORK, RWORK, INFO )
        #
        # .. Scalar Arguments ..
        # CHARACTER          HOWMNY, SIDE
        # INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
        # ..
        # .. Array Arguments ..
        # LOGICAL            SELECT( * )
        # DOUBLE PRECISION   RWORK( * )
        # COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
        #$                   WORK( * )
        function trevc!(side::AbstractChar, howmny::AbstractChar, select::AbstractVector{BlasInt}, T::AbstractMatrix{$elty},
                        VL::AbstractMatrix{$elty} = similar(T),
                        VR::AbstractMatrix{$elty} = similar(T))
            require_one_based_indexing(select, T, VL, VR)
            # Extract
            n, mm = checksquare(T), size(VL, 2)
            ldt, ldvl, ldvr = stride(T, 2), stride(VL, 2), stride(VR, 2)

            # Check
            chkstride1(T, select, VL, VR)
            if side ∉ ('L','R','B')
                throw(ArgumentError(lazy"side argument must be 'L' (left eigenvectors), 'R' (right eigenvectors), or 'B' (both), got $side"))
            end
            @chkvalidparam 2 howmny ('A', 'B', 'S')

            # Allocate
            m = Ref{BlasInt}()
            work = Vector{$elty}(undef, 2n)
            rwork = Vector{$relty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($trevc), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt},
                 Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                 Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                 Ptr{$elty}, Ptr{$relty}, Ref{BlasInt}, Clong, Clong),
                side, howmny, select, n,
                T, ldt, VL, ldvl,
                VR, ldvr, mm, m,
                work, rwork, info, 1, 1)
            chklapackerror(info[])

            VLn = size(VL, 1)
            VRn = size(VR, 1)
            #Decide what exactly to return
            if howmny == 'S' #compute selected eigenvectors
                if side == 'L' #left eigenvectors only
                    return select, reshape(resize!(vec(VL), VLn * m[]), VLn, m[])
                elseif side == 'R' #right eigenvectors only
                    return select, reshape(resize!(vec(VR), VRn * m[]), VRn, m[])
                else #side == 'B' #both eigenvectors
                    return select, reshape(resize!(vec(VL), VLn * m[]), VLn, m[]), reshape(resize!(vec(VR), VRn * m[]), VRn, m[])
                end
            else #compute all eigenvectors
                if side == 'L' #left eigenvectors only
                    return reshape(resize!(vec(VL), VLn * m[]), VLn, m[])
                elseif side == 'R' #right eigenvectors only
                    return reshape(resize!(vec(VR), VRn * m[]), VRn, m[])
                else #side == 'B' #both eigenvectors
                    return reshape(resize!(vec(VL), VLn * m[]), VLn, m[]), reshape(resize!(vec(VR), VRn * m[]), VRn, m[])
                end
            end
        end

        # SUBROUTINE ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X,
        #                    LDX, FERR, BERR, WORK, IWORK, INFO )
        # .. Scalar Arguments ..
        # CHARACTER          DIAG, TRANS, UPLO
        # INTEGER            INFO, LDA, LDB, LDX, N, NRHS
        # .. Array Arguments ..
        # INTEGER            IWORK( * )
        # DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), BERR( * ), FERR( * ),
        #$                   WORK( * ), X( LDX, * )
        function trrfs!(uplo::AbstractChar, trans::AbstractChar, diag::AbstractChar,
                        A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty}, X::AbstractVecOrMat{$elty},
                        Ferr::AbstractVector{$relty} = similar(B, $relty, size(B,2)),
                        Berr::AbstractVector{$relty} = similar(B, $relty, size(B,2)))
            require_one_based_indexing(A, B, X, Ferr, Berr)
            chkstride1(A, B, X, Ferr, Berr)
            chktrans(trans)
            chkuplo(uplo)
            chkdiag(diag)
            n = size(A,2)
            nrhs = size(B,2)
            if nrhs != size(X,2)
                throw(DimensionMismatch(lazy"second dimensions of B, $nrhs, and X, $(size(X,2)), must match"))
            end
            work  = Vector{$elty}(undef, 2n)
            rwork = Vector{$relty}(undef, n)
            info  = Ref{BlasInt}()
            ccall((@blasfunc($trrfs), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                 Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                 Ptr{$relty}, Ptr{$relty}, Ptr{$elty}, Ptr{$relty}, Ref{BlasInt}, Clong, Clong, Clong),
                uplo, trans, diag, n,
                nrhs, A, max(1,stride(A,2)), B, max(1,stride(B,2)), X, max(1,stride(X,2)),
                Ferr, Berr, work, rwork, info, 1, 1, 1)
            chklapackerror(info[])
            Ferr, Berr
        end
    end
end

"""
    trcon!(norm, uplo, diag, A)

Finds the reciprocal condition number of (upper if `uplo = U`, lower if
`uplo = L`) triangular matrix `A`. If `diag = N`, `A` has non-unit
diagonal elements. If `diag = U`, all diagonal elements of `A` are one.
If `norm = I`, the condition number is found in the infinity norm. If
`norm = O` or `1`, the condition number is found in the one norm.
"""
trcon!(norm::AbstractChar, uplo::AbstractChar, diag::AbstractChar, A::AbstractMatrix)

"""
    trevc!(side, howmny, select, T, VL = similar(T), VR = similar(T))

Finds the eigensystem of an upper triangular matrix `T`. If `side = R`,
the right eigenvectors are computed. If `side = L`, the left
eigenvectors are computed. If `side = B`, both sets are computed. If
`howmny = A`, all eigenvectors are found. If `howmny = B`, all
eigenvectors are found and backtransformed using `VL` and `VR`. If
`howmny = S`, only the eigenvectors corresponding to the values in
`select` are computed.
"""
trevc!(side::AbstractChar, howmny::AbstractChar, select::AbstractVector{BlasInt}, T::AbstractMatrix,
        VL::AbstractMatrix = similar(T), VR::AbstractMatrix = similar(T))

"""
    trrfs!(uplo, trans, diag, A, B, X, Ferr, Berr) -> (Ferr, Berr)

Estimates the error in the solution to `A * X = B` (`trans = N`),
`transpose(A) * X = B` (`trans = T`), `adjoint(A) * X = B` (`trans = C`) for `side = L`,
or the equivalent equations a right-handed `side = R` `X * A` after
computing `X` using `trtrs!`. If `uplo = U`, `A` is upper triangular.
If `uplo = L`, `A` is lower triangular. If `diag = N`, `A` has non-unit
diagonal elements. If `diag = U`, all diagonal elements of `A` are one.
`Ferr` and `Berr` are optional inputs. `Ferr` is the forward error and
`Berr` is the backward error, each component-wise.
"""
trrfs!(uplo::AbstractChar, trans::AbstractChar, diag::AbstractChar, A::AbstractMatrix, B::AbstractVecOrMat,
       X::AbstractVecOrMat, Ferr::AbstractVector, Berr::AbstractVector)

## (ST) Symmetric tridiagonal - eigendecomposition
for (stev, stebz, stegr, stein, elty) in
    ((:dstev_,:dstebz_,:dstegr_,:dstein_,:Float64),
     (:sstev_,:sstebz_,:sstegr_,:sstein_,:Float32)
#     , (:zstev_,:ComplexF64)  Need to rewrite for ZHEEV, rwork, etc.
#     , (:cstev_,:ComplexF32)
     )
    @eval begin
        function stev!(job::AbstractChar, dv::AbstractVector{$elty}, ev::AbstractVector{$elty})
            require_one_based_indexing(dv, ev)
            @chkvalidparam 1 job ('N', 'V')
            chkstride1(dv, ev)
            n = length(dv)
            if length(ev) != n - 1 && length(ev) != n
                throw(DimensionMismatch(lazy"ev has length $(length(ev)) but needs one less than or equal to dv's length, $n)"))
            end
            Zmat = similar(dv, $elty, (n, job != 'N' ? n : 0))
            work = Vector{$elty}(undef, max(1, 2n-2))
            info = Ref{BlasInt}()
            ccall((@blasfunc($stev), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty},
                   Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong),
                  job, n, dv, ev, Zmat, n, work, info, 1)
            chklapackerror(info[])
            dv, Zmat
        end

        #*  DSTEBZ computes the eigenvalues of a symmetric tridiagonal
        #*  matrix T.  The user may ask for all eigenvalues, all eigenvalues
        #*  in the half-open interval (VL, VU], or the IL-th through IU-th
        #*  eigenvalues.
        function stebz!(range::AbstractChar, order::AbstractChar, vl::$elty, vu::$elty, il::Integer, iu::Integer, abstol::Real, dv::AbstractVector{$elty}, ev::AbstractVector{$elty})
            require_one_based_indexing(dv, ev)
            @chkvalidparam 1 range ('A', 'V', 'I')
            @chkvalidparam 2 order ('B', 'E')
            chkstride1(dv, ev)
            n = length(dv)
            if length(ev) != n - 1
                throw(DimensionMismatch(lazy"ev has length $(length(ev)) but needs one less than dv's length, $n)"))
            end
            m = Ref{BlasInt}()
            nsplit = Vector{BlasInt}(undef, 1)
            w = similar(dv, $elty, n)
            tmp = 0.0
            iblock = similar(dv, BlasInt,n)
            isplit = similar(dv, BlasInt,n)
            work = Vector{$elty}(undef, 4*n)
            iwork = Vector{BlasInt}(undef, 3*n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($stebz), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{$elty},
                Ref{$elty}, Ref{BlasInt}, Ref{BlasInt}, Ref{$elty},
                Ptr{$elty}, Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt},
                Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty},
                Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong),
                range, order, n, vl,
                vu, il, iu, abstol,
                dv, ev, m, nsplit,
                w, iblock, isplit, work,
                iwork, info, 1, 1)
            chklapackerror(info[])
            resize!(w, m[]), resize!(iblock, m[]), resize!(isplit, nsplit[1])
        end

        function stegr!(jobz::AbstractChar, range::AbstractChar, dv::AbstractVector{$elty}, ev::AbstractVector{$elty}, vl::Real, vu::Real, il::Integer, iu::Integer)
            require_one_based_indexing(dv, ev)
            @chkvalidparam 1 jobz ('N', 'V')
            @chkvalidparam 2 range ('A', 'V', 'I')
            chkstride1(dv, ev)
            n = length(dv)
            ne = length(ev)
            if ne == n - 1
                eev = [ev; zero($elty)]
            elseif ne == n
                eev = copy(ev)
                eev[n] = zero($elty)
            else
                throw(DimensionMismatch(lazy"ev has length $ne but needs one less than or equal to dv's length, $n)"))
            end

            abstol = Vector{$elty}(undef, 1)
            m = Ref{BlasInt}()
            w = similar(dv, $elty, n)
            ldz = jobz == 'N' ? 1 : n
            Zn = range == 'I' ? iu-il+1 : n
            Z = similar(dv, $elty, ldz * Zn)
            isuppz = similar(dv, BlasInt, 2 * Zn)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            iwork = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1] and liwork as iwork[1]
                ccall((@blasfunc($stegr), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                    Ptr{$elty}, Ref{$elty}, Ref{$elty}, Ref{BlasInt},
                    Ref{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{$elty},
                    Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                    Clong, Clong),
                    jobz, range, n, dv,
                    eev, vl, vu, il,
                    iu, abstol, m, w,
                    Z, ldz, isuppz, work,
                    lwork, iwork, liwork, info,
                    1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                    liwork = iwork[1]
                    resize!(iwork, liwork)
                end
            end
            return resize!(w, m[]), reshape(resize!(Z, ldz * m[]), ldz, m[])
        end

        function stein!(dv::AbstractVector{$elty}, ev_in::AbstractVector{$elty}, w_in::AbstractVector{$elty}, iblock_in::AbstractVector{BlasInt}, isplit_in::AbstractVector{BlasInt})
            require_one_based_indexing(dv, ev_in, w_in, iblock_in, isplit_in)
            chkstride1(dv, ev_in, w_in, iblock_in, isplit_in)
            n = length(dv)
            ne = length(ev_in)
            if ne == n - 1
                ev = [ev_in; zero($elty)]
            elseif ne == n
                ev = copy(ev_in)
                ev[n] = zero($elty)
            else
                throw(DimensionMismatch(lazy"ev_in has length $ne but needs one less than or equal to dv's length, $n)"))
            end
            ldz = n #Leading dimension
            #Number of eigenvalues to find
            if !(1 <= length(w_in) <= n)
                throw(DimensionMismatch(lazy"w_in has length $(length(w_in)), but needs to be between 1 and $n"))
            end
            m = length(w_in)
            #If iblock and isplit are invalid input, assume worst-case block partitioning,
            # i.e. set the block scheme to be the entire matrix
            iblock = similar(dv, BlasInt,n)
            isplit = similar(dv, BlasInt,n)
            w = similar(dv, $elty,n)
            if length(iblock_in) < m #Not enough block specifications
                iblock[1:m] = fill(BlasInt(1), m)
                w[1:m] = sort(w_in)
            else
                iblock[1:m] = iblock_in
                w[1:m] = w_in #Assume user has sorted the eigenvalues properly
            end
            if length(isplit_in) < 1 #Not enough block specifications
                isplit[1] = n
            else
                isplit[1:length(isplit_in)] = isplit_in
            end
            z = similar(dv, $elty,(n,m))
            work  = Vector{$elty}(undef, 5*n)
            iwork = Vector{BlasInt}(undef, n)
            ifail = Vector{BlasInt}(undef, m)
            info  = Ref{BlasInt}()
            ccall((@blasfunc($stein), libblastrampoline), Cvoid,
                (Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{$elty},
                Ref{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt},
                Ptr{BlasInt}),
                n, dv, ev, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info)
            chklapackerror(info[])
            if any(ifail .!= 0)
                # TODO: better error message / type
                error("failed to converge eigenvectors:\n$(findall(!iszero, ifail))")
            end
            z
        end
    end
end
stegr!(jobz::AbstractChar, dv::AbstractVector, ev::AbstractVector) = stegr!(jobz, 'A', dv, ev, 0.0, 0.0, 0, 0)

# Allow user to skip specification of iblock and isplit
stein!(dv::AbstractVector, ev::AbstractVector, w_in::AbstractVector) = stein!(dv, ev, w_in, zeros(BlasInt,0), zeros(BlasInt,0))
# Allow user to specify just one eigenvector to get in stein!
stein!(dv::AbstractVector, ev::AbstractVector, eval::Real) = stein!(dv, ev, [eval], zeros(BlasInt,0), zeros(BlasInt,0))

"""
    stev!(job, dv, ev) -> (dv, Zmat)

Computes the eigensystem for a symmetric tridiagonal matrix with `dv` as
diagonal and `ev` as off-diagonal. If `job = N` only the eigenvalues are
found and returned in `dv`. If `job = V` then the eigenvectors are also found
and returned in `Zmat`.
"""
stev!(job::AbstractChar, dv::AbstractVector, ev::AbstractVector)

"""
    stebz!(range, order, vl, vu, il, iu, abstol, dv, ev) -> (dv, iblock, isplit)

Computes the eigenvalues for a symmetric tridiagonal matrix with `dv` as
diagonal and `ev` as off-diagonal. If `range = A`, all the eigenvalues
are found. If `range = V`, the eigenvalues in the half-open interval
`(vl, vu]` are found. If `range = I`, the eigenvalues with indices between
`il` and `iu` are found. If `order = B`, eigvalues are ordered within a
block. If `order = E`, they are ordered across all the blocks.
`abstol` can be set as a tolerance for convergence.
"""
stebz!(range::AbstractChar, order::AbstractChar, vl, vu, il::Integer, iu::Integer, abstol::Real, dv::AbstractVector, ev::AbstractVector)

"""
    stegr!(jobz, range, dv, ev, vl, vu, il, iu) -> (w, Z)

Computes the eigenvalues (`jobz = N`) or eigenvalues and eigenvectors
(`jobz = V`) for a symmetric tridiagonal matrix with `dv` as diagonal
and `ev` as off-diagonal. If `range = A`, all the eigenvalues
are found. If `range = V`, the eigenvalues in the half-open interval
`(vl, vu]` are found. If `range = I`, the eigenvalues with indices between
`il` and `iu` are found. The eigenvalues are returned in `w` and the eigenvectors
in `Z`.
"""
stegr!(jobz::AbstractChar, range::AbstractChar, dv::AbstractVector, ev::AbstractVector, vl::Real, vu::Real, il::Integer, iu::Integer)

"""
    stein!(dv, ev_in, w_in, iblock_in, isplit_in)

Computes the eigenvectors for a symmetric tridiagonal matrix with `dv`
as diagonal and `ev_in` as off-diagonal. `w_in` specifies the input
eigenvalues for which to find corresponding eigenvectors. `iblock_in`
specifies the submatrices corresponding to the eigenvalues in `w_in`.
`isplit_in` specifies the splitting points between the submatrix blocks.
"""
stein!(dv::AbstractVector, ev_in::AbstractVector, w_in::AbstractVector, iblock_in::AbstractVector{BlasInt}, isplit_in::AbstractVector{BlasInt})

## (SY) symmetric real matrices - Bunch-Kaufman decomposition,
## solvers (direct and factored) and inverse.
for (syconv, sysv, sytrf, sytri, sytrs, elty) in
    ((:dsyconv_,:dsysv_,:dsytrf_,:dsytri_,:dsytrs_,:Float64),
     (:ssyconv_,:ssysv_,:ssytrf_,:ssytri_,:ssytrs_,:Float32))
    @eval begin
        #       SUBROUTINE DSYCONV( UPLO, WAY, N, A, LDA, IPIV, WORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO, WAY
        #       INTEGER            INFO, LDA, N
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), WORK( * )
        function syconv!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A, ipiv)
            chkstride1(A, ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            work = Vector{$elty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($syconv), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong, Clong),
                  uplo, 'C', n, A, max(1,stride(A,2)), ipiv, work, info, 1, 1)
            chklapackerror(info[])
            A, work
        end

        #       SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
        #                         LWORK, INFO )
        #       .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
        #       .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
        function sysv!(uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A,B)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            ipiv  = similar(A, BlasInt, n)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($sysv), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)),
                      work, lwork, info,  1)
                chkargsok(info[])
                chknonsingular(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            B, A, ipiv
        end

        #       SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), WORK( * )
        function sytrf!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkuplo(uplo)
            if n == 0
                return A, ipiv, zero(BlasInt)
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($sytrf), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, A, stride(A,2), ipiv, work, lwork, info, 1)
                chkargsok(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            return A, ipiv, info[]
        end

        function sytrf!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkuplo(uplo)
            n = checksquare(A)
            ipiv = similar(A, BlasInt, n)
            sytrf!(uplo, A, ipiv)
        end

        #       SUBROUTINE DSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), WORK( * )
#         function sytri!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::Vector{BlasInt})
#             chkstride1(A)
#             n = checksquare(A)
#             chkuplo(uplo)
#             work  = Vector{$elty}(undef, 1)
#             lwork = BlasInt(-1)
#             info  = Ref{BlasInt}()
#             for i in 1:2
#                 ccall((@blasfunc($sytri), libblastrampoline), Cvoid,
#                       (Ptr{UInt8}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
#                        Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
#                       &uplo, &n, A, &max(1,stride(A,2)), ipiv, work, &lwork, info, 1)
#                 @assertargsok
#                 chknonsingular(info[])
#                 if lwork < 0
#                     lwork = BlasInt(real(work[1]))
#                     work = Vector{$elty}(undef, lwork)
#                 end
#             end
#             A
#         end

        #      SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
        #     .. Scalar Arguments ..
        #      CHARACTER          UPLO
        #      INTEGER            INFO, LDA, N
        #     .. Array Arguments ..
        #      INTEGER            IPIV( * )
        #      DOUBLE PRECISION   A( LDA, * ), WORK( * )
        function sytri!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A, ipiv)
            chkstride1(A, ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            work = Vector{$elty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($sytri), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong),
                  uplo, n, A, max(1,stride(A,2)), ipiv, work, info, 1)
            chkargsok(info[])
            chknonsingular(info[])
            A
        end

        #       SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        #
        #       .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, N, NRHS
        #       .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function sytrs!(uplo::AbstractChar, A::AbstractMatrix{$elty},
                       ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, ipiv, B)
            chkstride1(A,B,ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($sytrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info, 1)
            chklapackerror(info[])
            B
        end
    end
end

# Rook-pivoting variants of symmetric-matrix algorithms
for (sysv, sytrf, sytri, sytrs, syconvf, elty) in
    ((:dsysv_rook_,:dsytrf_rook_,:dsytri_rook_,:dsytrs_rook_,:dsyconvf_rook_,:Float64),
     (:ssysv_rook_,:ssytrf_rook_,:ssytri_rook_,:ssytrs_rook_,:ssyconvf_rook_,:Float32))
    @eval begin
        #       SUBROUTINE DSYSV_ROOK(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
        #                             LWORK, INFO )
        #       .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
        #       .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
        function sysv_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A,B)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            ipiv  = similar(A, BlasInt, n)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($sysv), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)),
                      work, lwork, info, 1)
                chkargsok(info[])
                chknonsingular(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            B, A, ipiv
        end

        #       SUBROUTINE DSYTRF_ROOK(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), WORK( * )
        function sytrf_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkuplo(uplo)
            ipiv  = similar(A, BlasInt, n)
            if n == 0
                return A, ipiv, zero(BlasInt)
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($sytrf), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, A, stride(A,2), ipiv, work, lwork, info, 1)
                chkargsok(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            return A, ipiv, info[]
        end

        #      SUBROUTINE DSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )
        #     .. Scalar Arguments ..
        #      CHARACTER          UPLO
        #      INTEGER            INFO, LDA, N
        #     .. Array Arguments ..
        #      INTEGER            IPIV( * )
        #      DOUBLE PRECISION   A( LDA, * ), WORK( * )
        function sytri_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A, ipiv)
            chkstride1(A, ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            work = Vector{$elty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($sytri), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong),
                  uplo, n, A, max(1,stride(A,2)), ipiv, work, info, 1)
            chkargsok(info[])
            chknonsingular(info[])
            A
        end

        #       SUBROUTINE DSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        #
        #       .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, N, NRHS
        #       .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        function sytrs_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty},
                       ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, ipiv, B)
            chkstride1(A,B,ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($sytrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info, 1)
            chklapackerror(info[])
            B
        end

        # SUBROUTINE DSYCONVF_ROOK( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
        #
        # .. Scalar Arguments ..
        # CHARACTER          UPLO, WAY
        # INTEGER            INFO, LDA, N
        # ..
        # .. Array Arguments ..
        # INTEGER            IPIV( * )
        # DOUBLE PRECISION   A( LDA, * ), E( * )
        function syconvf_rook!(uplo::AbstractChar, way::AbstractChar,
                                A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt},
                                e::AbstractVector{$elty} = Vector{$elty}(undef, length(ipiv)))
            require_one_based_indexing(A, ipiv, e)
            # extract
            n = checksquare(A)
            lda = max(1, stride(A, 2))

            # check
            chkuplo(uplo)
            if way != 'C' && way != 'R'
                throw(ArgumentError("way must be C or R"))
            end
            if length(ipiv) != n
                throw(ArgumentError(lazy"length of pivot vector was $(length(ipiv)) but should have been $n"))
            end
            if length(e) != n
                throw(ArgumentError(lazy"length of e vector was $(length(e)) but should have been $n"))
            end

            # allocate
            info = Ref{BlasInt}()

            ccall((@blasfunc($syconvf), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                 Ref{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt},
                 Clong, Clong),
                uplo, way, n, A,
                lda, e, ipiv, info,
                1, 1)

            chklapackerror(info[])
            return A, e
        end
    end
end

## (SY) hermitian matrices - eigendecomposition, Bunch-Kaufman decomposition,
## solvers (direct and factored) and inverse.
for (syconv, hesv, hetrf, hetri, hetrs, elty, relty) in
    ((:zsyconv_,:zhesv_,:zhetrf_,:zhetri_,:zhetrs_,:ComplexF64, :Float64),
     (:csyconv_,:chesv_,:chetrf_,:chetri_,:chetrs_,:ComplexF32, :Float32))
    @eval begin
       #   SUBROUTINE ZSYCONV( UPLO, WAY, N, A, LDA, IPIV, WORK, INFO )
       #
       #        .. Scalar Arguments ..
       #        CHARACTER          UPLO, WAY
       #        INTEGER            INFO, LDA, N
       #        ..
       #        .. Array Arguments ..
       #        INTEGER            IPIV( * )
       #        COMPLEX*16         A( LDA, * ), WORK( * )
        function syconv!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A,ipiv)
            chkstride1(A,ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            work = Vector{$elty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($syconv), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong, Clong),
                  uplo, 'C', n, A, max(1,stride(A,2)), ipiv, work, info, 1, 1)
            chklapackerror(info[])
            A, work
        end

        #       SUBROUTINE ZHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
        function hesv!(uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A,B)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            ipiv  = similar(A, BlasInt, n)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($hesv), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)),
                      work, lwork, info, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            B, A, ipiv
        end

        #       SUBROUTINE ZHETRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function hetrf!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkuplo(uplo)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i in 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($hetrf), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, A, max(1,stride(A,2)), ipiv, work, lwork, info, 1)
                chkargsok(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, ipiv, info[]
        end

        function hetrf!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkuplo(uplo)
            n = checksquare(A)
            ipiv = similar(A, BlasInt, n)
            hetrf!(uplo, A, ipiv)
        end

#       SUBROUTINE ZHETRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
# *     .. Scalar Arguments ..
#       CHARACTER          UPLO
#       INTEGER            INFO, LDA, LWORK, N
# *     ..
# *     .. Array Arguments ..
#       INTEGER            IPIV( * )
#       COMPLEX*16         A( LDA, * ), WORK( * )
#         function hetri!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::Vector{BlasInt})
#             chkstride1(A)
#             n = checksquare(A)
#             chkuplo(uplo)
#             work  = Vector{$elty}(undef, 1)
#             lwork = BlasInt(-1)
#             info  = Ref{BlasInt}()
#             for i in 1:2
#                 ccall((@blasfunc($hetri), libblastrampoline), Cvoid,
#                       (Ptr{UInt8}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
#                        Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
#                       &uplo, &n, A, &max(1,stride(A,2)), ipiv, work, &lwork, info, 1)
#                 chklapackerror(info[])
#                 if lwork < 0
#                     lwork = BlasInt(real(work[1]))
#                     work = Vector{$elty}(undef, lwork)
#                 end
#             end
#             A
#         end


        #       SUBROUTINE ZHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function hetri!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A, ipiv)
            chkstride1(A, ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            work = Vector{$elty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($hetri), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong),
                  uplo, n, A, max(1,stride(A,2)), ipiv, work, info, 1)
            chklapackerror(info[])
            A
        end

        #       SUBROUTINE ZHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * )
        function hetrs!(uplo::AbstractChar, A::AbstractMatrix{$elty},
                       ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, ipiv, B)
            chkuplo(uplo)
            chkstride1(A,B,ipiv)
            n = checksquare(A)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($hetrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info, 1)
            chklapackerror(info[])
            B
        end
    end
end

for (hesv, hetrf, hetri, hetrs, elty, relty) in
    ((:zhesv_rook_,:zhetrf_rook_,:zhetri_rook_,:zhetrs_rook_,:ComplexF64, :Float64),
     (:chesv_rook_,:chetrf_rook_,:chetri_rook_,:chetrs_rook_,:ComplexF32, :Float32))
    @eval begin
        #       SUBROUTINE ZHESV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
        function hesv_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A,B)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            ipiv  = similar(A, BlasInt, n)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($hesv), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)),
                      work, lwork, info, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            B, A, ipiv
        end

        #       SUBROUTINE ZHETRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function hetrf_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkuplo(uplo)
            ipiv  = similar(A, BlasInt, n)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i in 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($hetrf), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, A, max(1,stride(A,2)), ipiv, work, lwork, info, 1)
                chkargsok(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, ipiv, info[]
        end

        #       SUBROUTINE ZHETRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function hetri_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A,ipiv)
            chkstride1(A,ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            work = Vector{$elty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($hetri), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong),
                  uplo, n, A, max(1,stride(A,2)), ipiv, work, info, 1)
            chklapackerror(info[])
            A
        end

        #       SUBROUTINE ZHETRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * )
        function hetrs_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty},
                             ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, ipiv, B)
            chkstride1(A,B,ipiv)
            chkuplo(uplo)
            n = checksquare(A)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($hetrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info, 1)
            chklapackerror(info[])
            B
        end
    end
end

for (sysv, sytrf, sytri, sytrs, elty, relty) in
    ((:zsysv_,:zsytrf_,:zsytri_,:zsytrs_,:ComplexF64, :Float64),
     (:csysv_,:csytrf_,:csytri_,:csytrs_,:ComplexF32, :Float32))
    @eval begin
        #       SUBROUTINE ZSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
        #      $                  LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
        function sysv!(uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A,B)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            ipiv  = similar(A, BlasInt, n)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($sysv), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)),
                      work, lwork, info, 1)
                chkargsok(info[])
                chknonsingular(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            B, A, ipiv
        end

        #       SUBROUTINE ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function sytrf!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkuplo(uplo)
            if n == 0
                return A, ipiv, zero(BlasInt)
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($sytrf), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, A, max(1,stride(A,2)), ipiv, work, lwork, info, 1)
                chkargsok(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, ipiv, info[]
        end

        function sytrf!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkuplo(uplo)
            n = checksquare(A)
            ipiv = similar(A, BlasInt, n)
            sytrf!(uplo, A, ipiv)
        end

#       SUBROUTINE ZSYTRI2( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
# *     .. Scalar Arguments ..
#       CHARACTER          UPLO
#       INTEGER            INFO, LDA, LWORK, N
# *     ..
# *     .. Array Arguments ..
#       INTEGER            IPIV( * )
#       COMPLEX*16         A( LDA, * ), WORK( * )
#         function sytri!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::Vector{BlasInt})
#             chkstride1(A)
#             n = checksquare(A)
#             chkuplo(uplo)
#             work  = Vector{$elty}(undef, 1)
#             lwork = BlasInt(-1)
#             info  = Ref{BlasInt}()
#             for i in 1:2
#                 ccall((@blasfunc($sytri), libblastrampoline), Cvoid,
#                       (Ptr{UInt8}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt},
#                        Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
#                       &uplo, &n, A, &max(1,stride(A,2)), ipiv, work, &lwork, info, 1)
#                 chklapackerror(info[])
#                 if lwork < 0
#                     lwork = BlasInt(real(work[1]))
#                     work = Vector{$elty}(undef, lwork)
#                 end
#             end
#             A
#         end

        #       SUBROUTINE ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function sytri!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A, ipiv)
            chkstride1(A, ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            work = Vector{$elty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($sytri), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong),
                  uplo, n, A, max(1,stride(A,2)), ipiv, work, info, 1)
            chklapackerror(info[])
            A
        end

        #       SUBROUTINE ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * )
        function sytrs!(uplo::AbstractChar, A::AbstractMatrix{$elty},
                       ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, ipiv, B)
            chkstride1(A,B,ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($sytrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info,  1)
            chklapackerror(info[])
            B
        end
    end
end

for (sysv, sytrf, sytri, sytrs, syconvf, elty, relty) in
    ((:zsysv_rook_,:zsytrf_rook_,:zsytri_rook_,:zsytrs_rook_,:zsyconvf_rook_,:ComplexF64, :Float64),
     (:csysv_rook_,:csytrf_rook_,:csytri_rook_,:csytrs_rook_,:csyconvf_rook_,:ComplexF32, :Float32))
    @eval begin
        #       SUBROUTINE ZSYSV_ROOK(UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
        #      $                      LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
        function sysv_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, B)
            chkstride1(A,B)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            ipiv  = similar(A, BlasInt, n)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($sysv), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)),
                      work, lwork, info, 1)
                chkargsok(info[])
                chknonsingular(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            B, A, ipiv
        end

        #       SUBROUTINE ZSYTRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function sytrf_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkuplo(uplo)
            ipiv = similar(A, BlasInt, n)
            if n == 0
                return A, ipiv, zero(BlasInt)
            end
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($sytrf), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                       Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                      uplo, n, A, max(1,stride(A,2)), ipiv, work, lwork, info, 1)
                chkargsok(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, ipiv, info[]
        end

        #       SUBROUTINE ZSYTRI_ROOK( UPLO, N, A, LDA, IPIV, WORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function sytri_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
            require_one_based_indexing(A, ipiv)
            chkstride1(A, ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            work = Vector{$elty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($sytri), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong),
                  uplo, n, A, max(1,stride(A,2)), ipiv, work, info, 1)
            chklapackerror(info[])
            A
        end

        #       SUBROUTINE ZSYTRS_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LDB, N, NRHS
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IPIV( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * )
        function sytrs_rook!(uplo::AbstractChar, A::AbstractMatrix{$elty},
                             ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat{$elty})
            require_one_based_indexing(A, ipiv, B)
            chkstride1(A,B,ipiv)
            n = checksquare(A)
            chkuplo(uplo)
            if n != size(B,1)
                throw(DimensionMismatch(lazy"B has first dimension $(size(B,1)), but needs $n"))
            end
            info = Ref{BlasInt}()
            ccall((@blasfunc($sytrs), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ptr{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                  uplo, n, size(B,2), A, max(1,stride(A,2)), ipiv, B, max(1,stride(B,2)), info, 1)
            chklapackerror(info[])
            B
        end

        # SUBROUTINE ZSYCONVF_ROOK( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
        #
        # .. Scalar Arguments ..
        # CHARACTER          UPLO, WAY
        # INTEGER            INFO, LDA, N
        # ..
        # .. Array Arguments ..
        # INTEGER            IPIV( * )
        # COMPLEX*16         A( LDA, * ), E( * )
        function syconvf_rook!(uplo::AbstractChar, way::AbstractChar,
                                A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt},
                                e::AbstractVector{$elty} = Vector{$elty}(undef, length(ipiv)))
            require_one_based_indexing(A, ipiv, e)
            chkstride1(A, ipiv, e)

            # extract
            n   = checksquare(A)
            lda = stride(A, 2)

            # check
            chkuplo(uplo)
            if way != 'C' && way != 'R'
                throw(ArgumentError(lazy"way must be 'C' or 'R'"))
            end
            if length(ipiv) != n
                throw(ArgumentError(lazy"length of pivot vector was $(length(ipiv)) but should have been $n"))
            end
            if length(e) != n
                throw(ArgumentError(lazy"length of e vector was $(length(e)) but should have been $n"))
            end

            # allocate
            info = Ref{BlasInt}()

            ccall((@blasfunc($syconvf), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                 Ref{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt},
                 Clong, Clong),
                uplo, way, n, A,
                max(1, lda), e, ipiv, info,
                1, 1)

            chklapackerror(info[])
            return A, e
        end
    end
end

"""
    syconv!(uplo, A, ipiv) -> (A, work)

Converts a symmetric matrix `A` (which has been factorized into a
triangular matrix) into two matrices `L` and `D`. If `uplo = U`, `A`
is upper triangular. If `uplo = L`, it is lower triangular. `ipiv` is
the pivot vector from the triangular factorization. `A` is overwritten
by `L` and `D`.
"""
syconv!(uplo::AbstractChar, A::AbstractMatrix, ipiv::AbstractVector{BlasInt})

"""
    sysv!(uplo, A, B) -> (B, A, ipiv)

Finds the solution to `A * X = B` for symmetric matrix `A`. If `uplo = U`,
the upper half of `A` is stored. If `uplo = L`, the lower half is stored.
`B` is overwritten by the solution `X`. `A` is overwritten by its
Bunch-Kaufman factorization. `ipiv` contains pivoting information about the
factorization.
"""
sysv!(uplo::AbstractChar, A::AbstractMatrix, B::AbstractVecOrMat)

"""
    sytrf!(uplo, A) -> (A, ipiv, info)

Computes the Bunch-Kaufman factorization of a symmetric matrix `A`. If
`uplo = U`, the upper half of `A` is stored. If `uplo = L`, the lower
half is stored.

Returns `A`, overwritten by the factorization, a pivot vector `ipiv`, and
the error code `info` which is a non-negative integer. If `info` is positive
the matrix is singular and the diagonal part of the factorization is exactly
zero at position `info`.
"""
sytrf!(uplo::AbstractChar, A::AbstractMatrix)

"""
    sytrf!(uplo, A, ipiv) -> (A, ipiv, info)

Computes the Bunch-Kaufman factorization of a symmetric matrix `A`. If
`uplo = U`, the upper half of `A` is stored. If `uplo = L`, the lower
half is stored.

Returns `A`, overwritten by the factorization, the pivot vector `ipiv`, and
the error code `info` which is a non-negative integer. If `info` is positive
the matrix is singular and the diagonal part of the factorization is exactly
zero at position `info`.
"""
sytrf!(uplo::AbstractChar, A::AbstractMatrix, ipiv::AbstractVector{BlasInt})

"""
    sytri!(uplo, A, ipiv)

Computes the inverse of a symmetric matrix `A` using the results of
`sytrf!`. If `uplo = U`, the upper half of `A` is stored. If `uplo = L`,
the lower half is stored. `A` is overwritten by its inverse.
"""
sytri!(uplo::AbstractChar, A::AbstractMatrix, ipiv::AbstractVector{BlasInt})

"""
    sytrs!(uplo, A, ipiv, B)

Solves the equation `A * X = B` for a symmetric matrix `A` using the
results of `sytrf!`. If `uplo = U`, the upper half of `A` is stored.
If `uplo = L`, the lower half is stored. `B` is overwritten by the
solution `X`.
"""
sytrs!(uplo::AbstractChar, A::AbstractMatrix, ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat)


"""
    hesv!(uplo, A, B) -> (B, A, ipiv)

Finds the solution to `A * X = B` for Hermitian matrix `A`. If `uplo = U`,
the upper half of `A` is stored. If `uplo = L`, the lower half is stored.
`B` is overwritten by the solution `X`. `A` is overwritten by its
Bunch-Kaufman factorization. `ipiv` contains pivoting information about the
factorization.
"""
hesv!(uplo::AbstractChar, A::AbstractMatrix, B::AbstractVecOrMat)

"""
    hetrf!(uplo, A) -> (A, ipiv, info)

Computes the Bunch-Kaufman factorization of a Hermitian matrix `A`. If
`uplo = U`, the upper half of `A` is stored. If `uplo = L`, the lower
half is stored.

Returns `A`, overwritten by the factorization, a pivot vector `ipiv`, and
the error code `info` which is a non-negative integer. If `info` is positive
the matrix is singular and the diagonal part of the factorization is exactly
zero at position `info`.
"""
hetrf!(uplo::AbstractChar, A::AbstractMatrix)

"""
    hetrf!(uplo, A, ipiv) -> (A, ipiv, info)

Computes the Bunch-Kaufman factorization of a Hermitian matrix `A`. If
`uplo = U`, the upper half of `A` is stored. If `uplo = L`, the lower
half is stored.

Returns `A`, overwritten by the factorization, the pivot vector `ipiv`, and
the error code `info` which is a non-negative integer. If `info` is positive
the matrix is singular and the diagonal part of the factorization is exactly
zero at position `info`.
"""
hetrf!(uplo::AbstractChar, A::AbstractMatrix, ipiv::AbstractVector{BlasInt})

"""
    hetri!(uplo, A, ipiv)

Computes the inverse of a Hermitian matrix `A` using the results of
`sytrf!`. If `uplo = U`, the upper half of `A` is stored. If `uplo = L`,
the lower half is stored. `A` is overwritten by its inverse.
"""
hetri!(uplo::AbstractChar, A::AbstractMatrix, ipiv::AbstractVector{BlasInt})

"""
    hetrs!(uplo, A, ipiv, B)

Solves the equation `A * X = B` for a Hermitian matrix `A` using the
results of `sytrf!`. If `uplo = U`, the upper half of `A` is stored.
If `uplo = L`, the lower half is stored. `B` is overwritten by the
solution `X`.
"""
hetrs!(uplo::AbstractChar, A::AbstractMatrix, ipiv::AbstractVector{BlasInt}, B::AbstractVecOrMat)

for f in (:syevd!, :syev!)
    _f = Symbol(:_, f)
    @eval function $f(jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix)
        W, A = $_f(jobz, uplo, A)
        jobz == 'V' ? (W, A) : W
    end
end

# Symmetric (real) eigensolvers
for (syev, syevr, syevd, sygvd, elty) in
    ((:dsyev_,:dsyevr_,:dsyevd_,:dsygvd_,:Float64),
     (:ssyev_,:ssyevr_,:ssyevd_,:ssygvd_,:Float32))
    @eval begin
        #       SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
        Base.@constprop :none function _syev!(jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            @chkvalidparam 1 jobz ('N', 'V')
            chkuplo(uplo)
            chkstride1(A)
            n = checksquare(A)
            W     = similar(A, $elty, n)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($syev), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                      Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                      jobz, uplo, n, A, max(1,stride(A,2)), W, work, lwork, info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            W, A
        end

        #       SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
        #      $                   ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
        #      $                   IWORK, LIWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, RANGE, UPLO
        #       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
        #       DOUBLE PRECISION   ABSTOL, VL, VU
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            ISUPPZ( * ), IWORK( * )
        #       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
        function syevr!(jobz::AbstractChar, range::AbstractChar, uplo::AbstractChar, A::AbstractMatrix{$elty},
                        vl::AbstractFloat, vu::AbstractFloat, il::Integer, iu::Integer, abstol::AbstractFloat)
            require_one_based_indexing(A)
            @chkvalidparam 1 jobz ('N', 'V')
            @chkvalidparam 2 range ('A', 'V', 'I')
            chkstride1(A)
            n = checksquare(A)
            if range == 'I' && !(1 <= il <= iu <= n)
                throw(ArgumentError(lazy"illegal choice of eigenvalue indices (il = $il, iu = $iu), which must be between 1 and n = $n"))
            end
            if range == 'V' && vl >= vu
                throw(ArgumentError(lazy"lower boundary, $vl, must be less than upper boundary, $vu"))
            end
            chkuplofinite(A, uplo)
            lda = stride(A,2)
            m = Ref{BlasInt}()
            W = similar(A, $elty, n)
            ldz = n
            if jobz == 'N'
                Z = similar(A, $elty, 0)
            elseif jobz == 'V'
                Z = similar(A, $elty, ldz * n)
            end
            isuppz = similar(A, BlasInt, 2*n)
            work   = Vector{$elty}(undef, 1)
            lwork  = BlasInt(-1)
            iwork  = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            info   = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1] and liwork as iwork[1]
                ccall((@blasfunc($syevr), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                        Ptr{$elty}, Ref{BlasInt}, Ref{$elty}, Ref{$elty},
                        Ref{BlasInt}, Ref{BlasInt}, Ref{$elty}, Ptr{BlasInt},
                        Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                        Ref{BlasInt}, Clong, Clong, Clong),
                    jobz, range, uplo, n,
                    A, max(1,lda), vl, vu,
                    il, iu, abstol, m,
                    W, Z, max(1,ldz), isuppz,
                    work, lwork, iwork, liwork,
                    info, 1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    liwork = iwork[1]
                    resize!(iwork, liwork)
                end
            end
            zm = jobz == 'V' ? m[] : 0
            resize!(W, m[]), reshape(resize!(Z, ldz * zm), ldz, zm)
        end
        syevr!(jobz::AbstractChar, A::AbstractMatrix{$elty}) =
            syevr!(jobz, 'A', 'U', A, 0.0, 0.0, 0, 0, -1.0)

        #       SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK,
        #      $                   IWORK, LIWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, UPLO
        #       INTEGER            INFO, LDA, LIWORK, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
        Base.@constprop :none function _syevd!(jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            @chkvalidparam 1 jobz ('N', 'V')
            chkstride1(A)
            n = checksquare(A)
            chkuplofinite(A, uplo)
            lda = stride(A,2)
            m = Ref{BlasInt}()
            W = similar(A, $elty, n)
            work   = Vector{$elty}(undef, 1)
            lwork  = BlasInt(-1)
            iwork  = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            info   = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1] and liwork as iwork[1]
                ccall((@blasfunc($syevd), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                        Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                        Ref{BlasInt}, Clong, Clong),
                    jobz, uplo, n, A, max(1,lda),
                    W, work, lwork, iwork, liwork,
                    info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    liwork = iwork[1]
                    resize!(iwork, liwork)
                end
            end
            W, A
        end

        # Generalized eigenproblem
        #           SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
        #      $                   LWORK, IWORK, LIWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, UPLO
        #       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
        function sygvd!(itype::Integer, jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            @chkvalidparam 1 itype 1:3
            @chkvalidparam 2 jobz ('N', 'V')
            chkuplo(uplo)
            chkstride1(A, B)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"dimensions of A, ($n,$n), and B, ($m,$m), must match"))
            end
            lda = max(1, stride(A, 2))
            ldb = max(1, stride(B, 2))
            w = similar(A, $elty, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            iwork = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1] and liwork as iwork[1]
                ccall((@blasfunc($sygvd), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                     Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                     Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                    itype, jobz, uplo, n,
                    A, lda, B, ldb,
                    w, work, lwork, iwork,
                    liwork, info, 1, 1)
                chkargsok(info[])
                if i == 1
                    lwork = BlasInt(work[1])
                    resize!(work, lwork)
                    liwork = iwork[1]
                    resize!(iwork, liwork)
                end
            end
            chkposdef(info[])
            w, A, B
        end
    end
end
# Hermitian eigensolvers
for (syev, syevr, syevd, sygvd, elty, relty) in
    ((:zheev_,:zheevr_,:zheevd_,:zhegvd_,:ComplexF64,:Float64),
     (:cheev_,:cheevr_,:cheevd_,:chegvd_,:ComplexF32,:Float32))
    @eval begin
        # SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   RWORK( * ), W( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        Base.@constprop :none function _syev!(jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            @chkvalidparam 1 jobz ('N', 'V')
            chkstride1(A)
            chkuplofinite(A, uplo)
            n = checksquare(A)
            W     = similar(A, $relty, n)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, max(1, 3n-2))
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($syev), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                      Ptr{$relty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$relty}, Ptr{BlasInt},
                      Clong, Clong),
                      jobz, uplo, n, A, stride(A,2), W, work, lwork, rwork, info,
                      1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            W, A
        end

        #       SUBROUTINE ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
        #      $                   ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
        #      $                   RWORK, LRWORK, IWORK, LIWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, RANGE, UPLO
        #       INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK,
        #      $                   M, N
        #       DOUBLE PRECISION   ABSTOL, VL, VU
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            ISUPPZ( * ), IWORK( * )
        #       DOUBLE PRECISION   RWORK( * ), W( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )
        function syevr!(jobz::AbstractChar, range::AbstractChar, uplo::AbstractChar, A::AbstractMatrix{$elty},
                        vl::AbstractFloat, vu::AbstractFloat, il::Integer, iu::Integer, abstol::AbstractFloat)
            require_one_based_indexing(A)
            @chkvalidparam 1 jobz ('N', 'V')
            @chkvalidparam 2 range ('A', 'V', 'I')
            chkstride1(A)
            chkuplofinite(A, uplo)
            n = checksquare(A)
            if range == 'I' && !(1 <= il <= iu <= n)
                throw(ArgumentError(lazy"illegal choice of eigenvalue indices (il = $il, iu=$iu), which must be between 1 and n = $n"))
            end
            if range == 'V' && vl >= vu
                throw(ArgumentError(lazy"lower boundary, $vl, must be less than upper boundary, $vu"))
            end
            lda = max(1,stride(A,2))
            m = Ref{BlasInt}()
            W = similar(A, $relty, n)
            if jobz == 'N'
                ldz = 1
                Z = similar(A, $elty, 0)
            elseif jobz == 'V'
                ldz = n
                Z = similar(A, $elty, ldz * n)
            end
            isuppz = similar(A, BlasInt, 2*n)
            work   = Vector{$elty}(undef, 1)
            lwork  = BlasInt(-1)
            rwork  = Vector{$relty}(undef, 1)
            lrwork = BlasInt(-1)
            iwork  = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            info   = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
                ccall((@blasfunc($syevr), libblastrampoline), Cvoid,
                      (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ref{$elty}, Ref{$elty},
                       Ref{BlasInt}, Ref{BlasInt}, Ref{$elty}, Ptr{BlasInt},
                       Ptr{$relty}, Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt},
                       Ptr{$elty}, Ref{BlasInt}, Ptr{$relty}, Ref{BlasInt},
                       Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                       Clong, Clong, Clong),
                      jobz, range, uplo, n,
                      A, lda, vl, vu,
                      il, iu, abstol, m,
                      W, Z, ldz, isuppz,
                      work, lwork, rwork, lrwork,
                      iwork, liwork, info,
                      1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    lrwork = BlasInt(rwork[1])
                    resize!(rwork, lrwork)
                    liwork = iwork[1]
                    resize!(iwork, liwork)
                end
            end
            zm = jobz == 'V' ? m[] : 0
            resize!(W, m[]), reshape(resize!(Z, ldz * zm), ldz, zm)
        end
        syevr!(jobz::AbstractChar, A::AbstractMatrix{$elty}) =
            syevr!(jobz, 'A', 'U', A, 0.0, 0.0, 0, 0, -1.0)

        #       SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
        #      $                   LRWORK, IWORK, LIWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, UPLO
        #       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   RWORK( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        Base.@constprop :none function _syevd!(jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            @chkvalidparam 1 jobz ('N', 'V')
            chkstride1(A)
            chkuplofinite(A, uplo)
            n = checksquare(A)
            lda = max(1, stride(A,2))
            m = Ref{BlasInt}()
            W = similar(A, $relty, n)
            work   = Vector{$elty}(undef, 1)
            lwork  = BlasInt(-1)
            rwork  = Vector{$relty}(undef, 1)
            lrwork = BlasInt(-1)
            iwork  = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            info   = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
                ccall((@blasfunc($syevd), liblapack), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$relty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$relty}, Ref{BlasInt},
                    Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                    jobz, uplo, n, A, stride(A,2),
                    W, work, lwork, rwork, lrwork,
                    iwork, liwork, info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    lrwork = BlasInt(rwork[1])
                    resize!(rwork, lrwork)
                    liwork = iwork[1]
                    resize!(iwork, liwork)
                end
            end
            W, A
        end

        #       SUBROUTINE ZHEGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
        #      $                   LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBZ, UPLO
        #       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   RWORK( * ), W( * )
        #       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
        function sygvd!(itype::Integer, jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            @chkvalidparam 1 itype 1:3
            @chkvalidparam 2 jobz ('N', 'V')
            chkstride1(A, B)
            chkuplofinite(A, uplo)
            chkuplofinite(B, uplo)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"dimensions of A, ($n,$n), and B, ($m,$m), must match"))
            end
            lda = max(1, stride(A, 2))
            ldb = max(1, stride(B, 2))
            w = similar(A, $relty, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            iwork = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, 1)
            lrwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1], lrwork as rwork[1] and liwork as iwork[1]
                ccall((@blasfunc($sygvd), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                     Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$relty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$relty},
                     Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                     Clong, Clong),
                    itype, jobz, uplo, n,
                    A, lda, B, ldb,
                    w, work, lwork, rwork,
                    lrwork, iwork, liwork, info,
                    1, 1)
                chkargsok(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    liwork = iwork[1]
                    resize!(iwork, liwork)
                    lrwork = BlasInt(rwork[1])
                    resize!(rwork, lrwork)
                end
            end
            chkposdef(info[])
            w, A, B
        end
    end
end

"""
    syev!(jobz, uplo, A)

Finds the eigenvalues (`jobz = N`) or eigenvalues and eigenvectors
(`jobz = V`) of a symmetric matrix `A`. If `uplo = U`, the upper triangle
of `A` is used. If `uplo = L`, the lower triangle of `A` is used.
"""
syev!(jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix)

"""
    syevr!(jobz, range, uplo, A, vl, vu, il, iu, abstol) -> (W, Z)

Finds the eigenvalues (`jobz = N`) or eigenvalues and eigenvectors
(`jobz = V`) of a symmetric matrix `A`. If `uplo = U`, the upper triangle
of `A` is used. If `uplo = L`, the lower triangle of `A` is used. If
`range = A`, all the eigenvalues are found. If `range = V`, the
eigenvalues in the half-open interval `(vl, vu]` are found.
If `range = I`, the eigenvalues with indices between `il` and `iu` are
found. `abstol` can be set as a tolerance for convergence.

The eigenvalues are returned in `W` and the eigenvectors in `Z`.
"""
syevr!(jobz::AbstractChar, range::AbstractChar, uplo::AbstractChar, A::AbstractMatrix,
       vl::AbstractFloat, vu::AbstractFloat, il::Integer, iu::Integer, abstol::AbstractFloat)

"""
    syevd!(jobz, uplo, A)

Finds the eigenvalues (`jobz = N`) or eigenvalues and eigenvectors
(`jobz = V`) of a symmetric matrix `A`. If `uplo = U`, the upper triangle
of `A` is used. If `uplo = L`, the lower triangle of `A` is used.
"""
syevd!(jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix)

"""
    sygvd!(itype, jobz, uplo, A, B) -> (w, A, B)

Finds the generalized eigenvalues (`jobz = N`) or eigenvalues and
eigenvectors (`jobz = V`) of a symmetric matrix `A` and symmetric
positive-definite matrix `B`. If `uplo = U`, the upper triangles
of `A` and `B` are used. If `uplo = L`, the lower triangles of `A` and
`B` are used. If `itype = 1`, the problem to solve is
`A * x = lambda * B * x`. If `itype = 2`, the problem to solve is
`A * B * x = lambda * x`. If `itype = 3`, the problem to solve is
`B * A * x = lambda * x`.
"""
sygvd!(itype::Integer, jobz::AbstractChar, uplo::AbstractChar, A::AbstractMatrix, B::AbstractMatrix)

## (BD) Bidiagonal matrices - singular value decomposition
for (bdsqr, relty, elty) in
    ((:dbdsqr_,:Float64,:Float64),
     (:sbdsqr_,:Float32,:Float32),
     (:zbdsqr_,:Float64,:ComplexF64),
     (:cbdsqr_,:Float32,:ComplexF32))
    @eval begin
        function bdsqr!(uplo::AbstractChar, d::AbstractVector{$relty}, e_::AbstractVector{$relty},
                        Vt::AbstractMatrix{$elty}, U::AbstractMatrix{$elty}, C::AbstractMatrix{$elty})
            require_one_based_indexing(d, e_, Vt, U, C)
            chkstride1(d, e_, Vt, U, C)
            # Extract number
            n = length(d)
            ncvt, nru, ncc = size(Vt, 2), size(U, 1), size(C, 2)
            ldvt, ldu, ldc = max(1, stride(Vt,2)), max(1, stride(U, 2)), max(1, stride(C,2))
            # Do checks
            chkuplo(uplo)
            if length(e_) != n - 1
                throw(DimensionMismatch(lazy"off-diagonal has length $(length(e_)) but should have length $(n - 1)"))
            end
            if ncvt > 0 && ldvt < n
                throw(DimensionMismatch(lazy"leading dimension of Vt, $ldvt, must be at least $n"))
            end
            if ldu < nru
                throw(DimensionMismatch(lazy"leading dimension of U, $ldu, must be at least $nru"))
            end
            if size(U, 2) != n
                throw(DimensionMismatch(lazy"U must have $n columns but has $(size(U, 2))"))
            end
            if ncc > 0 && ldc < n
                throw(DimensionMismatch(lazy"leading dimension of C, $ldc, must be at least $n"))
            end
            # Allocate
            work = Vector{$relty}(undef, 4n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($bdsqr), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                 Ref{BlasInt}, Ptr{$relty}, Ptr{$relty}, Ptr{$elty},
                 Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                 Ref{BlasInt}, Ptr{$relty}, Ref{BlasInt}, Clong),
                uplo, n, ncvt, nru,
                ncc, d, e_, Vt,
                ldvt, U, ldu, C,
                ldc, work, info, 1)
            chklapackerror(info[])
            d, Vt, U, C #singular values in descending order, P**T * VT, U * Q, Q**T * C
        end
    end
end

"""
    bdsqr!(uplo, d, e_, Vt, U, C) -> (d, Vt, U, C)

Computes the singular value decomposition of a bidiagonal matrix with
`d` on the diagonal and `e_` on the off-diagonal. If `uplo = U`, `e_` is
the superdiagonal. If `uplo = L`, `e_` is the subdiagonal. Can optionally also
compute the product `Q' * C`.

Returns the singular values in `d`, and the matrix `C` overwritten with `Q' * C`.
"""
bdsqr!(uplo::AbstractChar, d::AbstractVector, e_::AbstractVector, Vt::AbstractMatrix, U::AbstractMatrix, C::AbstractMatrix)

#Defined only for real types
for (bdsdc, elty) in
    ((:dbdsdc_,:Float64),
     (:sbdsdc_,:Float32))
    @eval begin
        #*  DBDSDC computes the singular value decomposition (SVD) of a real
        #*  N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
        #*  using a divide and conquer method
        #*     .. Scalar Arguments ..
        #      CHARACTER          COMPQ, UPLO
        #      INTEGER            INFO, LDU, LDVT, N
        #*     ..
        #*     .. Array Arguments ..
        #      INTEGER            IQ( * ), IWORK( * )
        #      DOUBLE PRECISION   D( * ), E( * ), Q( * ), U( LDU, * ),
        #     $                   VT( LDVT, * ), WORK( * )
        function bdsdc!(uplo::AbstractChar, compq::AbstractChar, d::AbstractVector{$elty}, e_::AbstractVector{$elty})
            require_one_based_indexing(d, e_)
            chkstride1(d, e_)
            n, ldiq, ldq, ldu, ldvt = length(d), 1, 1, 1, 1
            chkuplo(uplo)
            if compq == 'N'
                lwork = 6*n
            elseif compq == 'P'
                @warn "COMPQ='P' is not tested"
                #TODO turn this into an actual LAPACK call
                #smlsiz=ilaenv(9, $elty === :Float64 ? 'dbdsqr' : 'sbdsqr', string(uplo, compq), n,n,n,n)
                smlsiz=100 #For now, completely overkill
                ldq = n*(11+2*smlsiz+8*round(Int,log((n/(smlsiz+1)))/log(2)))
                ldiq = n*(3+3*round(Int,log(n/(smlsiz+1))/log(2)))
                lwork = 6*n
            elseif compq == 'I'
                ldvt=ldu=max(1, n)
                lwork=3*n^2 + 4*n
            else
                throw(ArgumentError(lazy"COMPQ argument must be 'N', 'P' or 'I', got $(repr(compq))"))
            end
            u  = similar(d, $elty, (ldu,  n))
            vt = similar(d, $elty, (ldvt, n))
            q  = similar(d, $elty, ldq)
            iq = similar(d, BlasInt, ldiq)
            work  = Vector{$elty}(undef, lwork)
            iwork = Vector{BlasInt}(undef, 8n)
            info  = Ref{BlasInt}()
            ccall((@blasfunc($bdsdc), libblastrampoline), Cvoid,
               (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty},
                Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                Ptr{$elty}, Ptr{BlasInt}, Ptr{$elty}, Ptr{BlasInt}, Ptr{BlasInt},
                Clong, Clong),
                uplo, compq, n, d, e_,
                u, ldu, vt, ldvt,
                q, iq, work, iwork, info,
                1, 1)
            chklapackerror(info[])
            d, e_, u, vt, q, iq
        end
    end
end

"""
    bdsdc!(uplo, compq, d, e_) -> (d, e, u, vt, q, iq)

Computes the singular value decomposition of a bidiagonal matrix with `d` on the
diagonal and `e_` on the off-diagonal using a divide and conqueq method.
If `uplo = U`, `e_` is the superdiagonal. If `uplo = L`, `e_` is the subdiagonal.
If `compq = N`, only the singular values are found. If `compq = I`, the singular
values and vectors are found. If `compq = P`, the singular values
and vectors are found in compact form. Only works for real types.

Returns the singular values in `d`, and if `compq = P`, the compact singular
vectors in `iq`.
"""
bdsdc!(uplo::AbstractChar, compq::AbstractChar, d::AbstractVector, e_::AbstractVector)

for (gecon, elty) in
    ((:dgecon_,:Float64),
     (:sgecon_,:Float32))
    @eval begin
        #  SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,
        #      $                   INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          NORM
        #       INTEGER            INFO, LDA, N
        #       DOUBLE PRECISION   ANORM, RCOND
        # *     ..
        # *     .. Array Arguments ..
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   A( LDA, * ), WORK( * )
        function gecon!(normtype::AbstractChar, A::AbstractMatrix{$elty}, anorm::$elty)
            require_one_based_indexing(A)
            @chkvalidparam 1 normtype ('0', '1', 'I')
            chkstride1(A)
            n = checksquare(A)
            lda = max(1, stride(A, 2))
            rcond = Ref{$elty}()
            work = Vector{$elty}(undef, 4n)
            iwork = Vector{BlasInt}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($gecon), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ref{$elty}, Ref{$elty}, Ptr{$elty}, Ptr{BlasInt},
                   Ref{BlasInt}, Clong),
                  normtype, n, A, lda, anorm, rcond, work, iwork,
                  info, 1)
            chklapackerror(info[])
            rcond[]
        end
    end
end

for (gecon, elty, relty) in
    ((:zgecon_,:ComplexF64,:Float64),
     (:cgecon_,:ComplexF32,:Float32))
    @eval begin
        #       SUBROUTINE ZGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK,
        #      $                   INFO )
        # *     .. Scalar Arguments ..
        #       CHARACTER          NORM
        #       INTEGER            INFO, LDA, N
        #       DOUBLE PRECISION   ANORM, RCOND
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   RWORK( * )
        #       COMPLEX*16         A( LDA, * ), WORK( * )
        function gecon!(normtype::AbstractChar, A::AbstractMatrix{$elty}, anorm::$relty)
            require_one_based_indexing(A)
            @chkvalidparam 1 normtype ('0', '1', 'I')
            chkstride1(A)
            n = checksquare(A)
            lda = max(1, stride(A, 2))
            rcond = Ref{$relty}()
            work = Vector{$elty}(undef, 2n)
            rwork = Vector{$relty}(undef, 2n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($gecon), libblastrampoline), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ref{$relty}, Ref{$relty}, Ptr{$elty}, Ptr{$relty},
                   Ref{BlasInt}, Clong),
                  normtype, n, A, lda, anorm, rcond, work, rwork,
                  info, 1)
            chklapackerror(info[])
            rcond[]
        end
    end
end

"""
    gecon!(normtype, A, anorm)

Finds the reciprocal condition number of matrix `A`. If `normtype = I`,
the condition number is found in the infinity norm. If `normtype = O` or
`1`, the condition number is found in the one norm. `A` must be the
result of `getrf!` and `anorm` is the norm of `A` in the relevant norm.
"""
gecon!(normtype::AbstractChar, A::AbstractMatrix, anorm)

for (gehrd, elty) in
    ((:dgehrd_,:Float64),
     (:sgehrd_,:Float32),
     (:zgehrd_,:ComplexF64),
     (:cgehrd_,:ComplexF32))
    @eval begin

        #                 .. Scalar Arguments ..
        #       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION  A( LDA, * ), TAU( * ), WORK( * )
        function gehrd!(ilo::Integer, ihi::Integer, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkfinite(A) # balancing routines don't support NaNs and Infs
            tau = similar(A, $elty, max(0,n - 1))
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gehrd), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{BlasInt}),
                    n, ilo, ihi, A,
                    max(1, stride(A, 2)), tau, work, lwork,
                    info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, tau
        end
    end
end
gehrd!(A::AbstractMatrix) = gehrd!(1, size(A, 1), A)

"""
    gehrd!(ilo, ihi, A) -> (A, tau)

Converts a matrix `A` to Hessenberg form. If `A` is balanced with `gebal!`
then `ilo` and `ihi` are the outputs of `gebal!`. Otherwise they should be
`ilo = 1` and `ihi = size(A,2)`. `tau` contains the elementary reflectors of
the factorization.
"""
gehrd!(ilo::Integer, ihi::Integer, A::AbstractMatrix)

for (orghr, elty) in
    ((:dorghr_,:Float64),
     (:sorghr_,:Float32),
     (:zunghr_,:ComplexF64),
     (:cunghr_,:ComplexF32))
    @eval begin
        # *     .. Scalar Arguments ..
        #       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function orghr!(ilo::Integer, ihi::Integer, A::AbstractMatrix{$elty}, tau::AbstractVector{$elty})
            require_one_based_indexing(A, tau)
            chkstride1(A, tau)
            n = checksquare(A)
            if n - length(tau) != 1
                throw(DimensionMismatch(lazy"tau has length $(length(tau)), needs $(n - 1)"))
            end
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($orghr), libblastrampoline), Cvoid,
                    (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{BlasInt}),
                    n, ilo, ihi, A,
                    max(1, stride(A, 2)), tau, work, lwork,
                    info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A
        end
    end
end

"""
    orghr!(ilo, ihi, A, tau)

Explicitly finds `Q`, the orthogonal/unitary matrix from `gehrd!`. `ilo`,
`ihi`, `A`, and `tau` must correspond to the input/output to `gehrd!`.
"""
orghr!(ilo::Integer, ihi::Integer, A::AbstractMatrix, tau::AbstractVector)

for (ormhr, elty) in
    ((:dormhr_,:Float64),
     (:sormhr_,:Float32),
     (:zunmhr_,:ComplexF64),
     (:cunmhr_,:ComplexF32))
    @eval begin
        # .. Scalar Arguments ..
        # CHARACTER          side, trans
        # INTEGER            ihi, ilo, info, lda, ldc, lwork, m, n
        # ..
        # .. Array Arguments ..
        # DOUBLE PRECISION   a( lda, * ), c( ldc, * ), tau( * ), work( * )
        function ormhr!(side::AbstractChar, trans::AbstractChar, ilo::Integer, ihi::Integer, A::AbstractMatrix{$elty},
            tau::AbstractVector{$elty}, C::AbstractVecOrMat{$elty})

            require_one_based_indexing(A, tau, C)
            chkstride1(A, tau, C)
            chkside(side)
            chktrans(trans)
            n = checksquare(A)
            ntau = length(tau)
            mC, nC = size(C, 1), size(C, 2)

            if max(n - 1, 0) != ntau
                throw(DimensionMismatch(lazy"tau has length $ntau, needs $(n - 1)"))
            end
            if (side == 'L' && mC != n) || (side == 'R' && nC != n)
                throw(DimensionMismatch("A and C matrices are not conformable"))
            end

            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ormhr), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                     Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                     Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                    side, trans, mC, nC,
                    ilo, ihi, A, max(1, stride(A, 2)),
                    tau, C, max(1, stride(C, 2)), work,
                    lwork, info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            C
        end
    end
end

for (hseqr, elty) in
    ((:zhseqr_,:ComplexF64),
     (:chseqr_,:ComplexF32))
    @eval begin
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOB, COMPZ
        #       INTEGER            N, ILO, IHI, LWORK, LDH, LDZ, INFO
        # *     ..
        # *     .. Array Arguments ..
        #       COMPLEX*16         H( LDH, * ), Z( LDZ, * ), WORK( * )
        function hseqr!(job::AbstractChar, compz::AbstractChar, ilo::Integer, ihi::Integer,
                        H::AbstractMatrix{$elty}, Z::AbstractMatrix{$elty})
            require_one_based_indexing(H, Z)
            @chkvalidparam 1 job ('E', 'S')
            @chkvalidparam 2 compz ('N', 'I', 'V')
            chkstride1(H)
            n = checksquare(H)
            checksquare(Z) == n || throw(DimensionMismatch())
            ldh = max(1, stride(H, 2))
            ldz = max(1, stride(Z, 2))
            w = similar(H, $elty, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($hseqr), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                    Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{BlasInt}),
                    job, compz, n, ilo, ihi,
                    H, ldh, w, Z, ldz, work,
                    lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            H, Z, w
        end
    end
end

for (hseqr, elty) in
    ((:dhseqr_,:Float64),
     (:shseqr_,:Float32))
    @eval begin
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOB, COMPZ
        #       INTEGER            N, ILO, IHI, LWORK, LDH, LDZ, INFO
        # *     ..
        # *     .. Array Arguments ..
        #       COMPLEX*16         H( LDH, * ), Z( LDZ, * ), WORK( * )
        function hseqr!(job::AbstractChar, compz::AbstractChar, ilo::Integer, ihi::Integer,
                        H::AbstractMatrix{$elty}, Z::AbstractMatrix{$elty})
            require_one_based_indexing(H, Z)
            @chkvalidparam 1 job ('E', 'S')
            @chkvalidparam 2 compz ('N', 'I', 'V')
            chkstride1(H)
            n = checksquare(H)
            checksquare(Z) == n || throw(DimensionMismatch())
            ldh = max(1, stride(H, 2))
            ldz = max(1, stride(Z, 2))
            wr = similar(H, $elty, n)
            wi = similar(H, $elty, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($hseqr), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                    Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{BlasInt}),
                    job, compz, n, ilo, ihi,
                    H, ldh, wr, wi, Z, ldz, work,
                    lwork, info)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            H, Z, complex.(wr, wi)
        end
    end
end
hseqr!(H::StridedMatrix{T}, Z::StridedMatrix{T}) where {T<:BlasFloat} = hseqr!('S', 'V', 1, size(H, 1), H, Z)
hseqr!(H::StridedMatrix{T}) where {T<:BlasFloat} = hseqr!('S', 'I', 1, size(H, 1), H, similar(H))

"""
    hseqr!(job, compz, ilo, ihi, H, Z) -> (H, Z, w)

Computes all eigenvalues and (optionally) the Schur factorization of a matrix
reduced to Hessenberg form. If `H` is balanced with `gebal!`
then `ilo` and `ihi` are the outputs of `gebal!`. Otherwise they should be
`ilo = 1` and `ihi = size(H,2)`. `tau` contains the elementary reflectors of
the factorization.
"""
hseqr!(job::AbstractChar, compz::AbstractChar, ilo::Integer, ihi::Integer, H::AbstractMatrix, Z::AbstractMatrix)

for (hetrd, elty) in
    ((:dsytrd_,Float64),
     (:ssytrd_,Float32),
     (:zhetrd_,ComplexF64),
     (:chetrd_,ComplexF32))
    relty = real(elty)
    @eval begin

        #                 .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION  A( LDA, * ), D( * ), E( * ), TAU( * ), WORK( * )
        function hetrd!(uplo::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            chkstride1(A)
            n = checksquare(A)
            chkuplofinite(A, uplo) # balancing routines don't support NaNs and Infs
            tau = similar(A, $elty, max(0,n - 1))
            d = Vector{$relty}(undef, n)
            e = Vector{$relty}(undef, max(0,n - 1))
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($hetrd), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$relty}, Ptr{$relty},
                    Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                    uplo, n, A, max(1, stride(A, 2)), d, e, tau, work, lwork, info, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, tau, d, e
        end
    end
end

"""
    hetrd!(uplo, A) -> (A, tau, d, e)

Converts a Hermitian matrix `A` to real-symmetric tridiagonal Hessenberg form.
If `uplo = U`, the upper half of `A` is stored; if `uplo = L`, the lower half is stored.
`tau` contains the elementary reflectors of the factorization, `d` contains the
diagonal and `e` contains the upper/lower diagonal.
"""
hetrd!(uplo::AbstractChar, A::AbstractMatrix)

for (orgtr, elty) in
    ((:dorgtr_,:Float64),
     (:sorgtr_,:Float32),
     (:zungtr_,:ComplexF64),
     (:cungtr_,:ComplexF32))
    @eval begin
        # *     .. Scalar Arguments ..
        #       CHARACTER          UPLO
        #       INTEGER            INFO, LDA, LWORK, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
        function orgtr!(uplo::AbstractChar, A::AbstractMatrix{$elty}, tau::AbstractVector{$elty})
            require_one_based_indexing(A, tau)
            chkstride1(A, tau)
            n = checksquare(A)
            if n - length(tau) != 1
                throw(DimensionMismatch(lazy"tau has length $(length(tau)), needs $(n - 1)"))
            end
            chkuplo(uplo)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($orgtr), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ref{BlasInt},
                     Ref{BlasInt}, Clong),
                    uplo, n, A,
                    max(1, stride(A, 2)), tau, work, lwork,
                    info, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A
        end
    end
end

"""
    orgtr!(uplo, A, tau)

Explicitly finds `Q`, the orthogonal/unitary matrix from `hetrd!`. `uplo`,
`A`, and `tau` must correspond to the input/output to `hetrd!`.
"""
orgtr!(uplo::AbstractChar, A::AbstractMatrix, tau::AbstractVector)

for (ormtr, elty) in
    ((:dormtr_,:Float64),
     (:sormtr_,:Float32),
     (:zunmtr_,:ComplexF64),
     (:cunmtr_,:ComplexF32))
    @eval begin
        # .. Scalar Arguments ..
        # CHARACTER          side, trans, uplo
        # INTEGER            info, lda, ldc, lwork, m, n
        # ..
        # .. Array Arguments ..
        # DOUBLE PRECISION   a( lda, * ), c( ldc, * ), tau( * ), work( * )
        function ormtr!(side::AbstractChar, uplo::AbstractChar, trans::AbstractChar, A::AbstractMatrix{$elty},
                        tau::AbstractVector{$elty}, C::AbstractVecOrMat{$elty})

            require_one_based_indexing(A, tau, C)
            chkstride1(A, tau, C)
            n = checksquare(A)
            chkside(side)
            chkuplo(uplo)
            chktrans(trans)
            mC, nC = size(C, 1), size(C, 2)

            if n - length(tau) != 1
                throw(DimensionMismatch(lazy"tau has length $(length(tau)), needs $(n - 1)"))
            end
            if (side == 'L' && mC != n) || (side == 'R' && nC != n)
                throw(DimensionMismatch(lazy"A and C matrices are not conformable"))
            end

            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($ormtr), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                     Ptr{$elty}, Ref{BlasInt},
                     Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                     Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                    side, uplo, trans, mC, nC,
                    A, max(1, stride(A, 2)),
                    tau, C, max(1, stride(C, 2)), work,
                    lwork, info, 1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            C
        end
    end
end

for (gees, gges, gges3, elty) in
    ((:dgees_,:dgges_,:dgges3_,:Float64),
     (:sgees_,:sgges_,:sgges3_,:Float32))
    @eval begin
        #     .. Scalar Arguments ..
        #     CHARACTER          JOBVS, SORT
        #     INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
        #     ..
        #     .. Array Arguments ..
        #     LOGICAL            BWORK( * )
        #     DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),
        #    $                   WR( * )
        function gees!(jobvs::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            @chkvalidparam 1 jobvs ('N', 'V')
            chkstride1(A)
            n     = checksquare(A)
            sdim  = Vector{BlasInt}(undef, 1)
            wr    = similar(A, $elty, n)
            wi    = similar(A, $elty, n)
            vs    = similar(A, $elty, jobvs == 'V' ? n : 0, n)
            ldvs  = max(size(vs, 1), 1)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gees), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ptr{Cvoid}, Ref{BlasInt},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{$elty},
                        Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                        Ref{BlasInt}, Ptr{Cvoid}, Ref{BlasInt}, Clong, Clong),
                    jobvs, 'N', C_NULL, n,
                        A, max(1, stride(A, 2)), sdim, wr,
                        wi, vs, ldvs, work,
                        lwork, C_NULL, info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            iszero(wi) ? (A, vs, wr) : (A, vs, complex.(wr, wi))
        end

        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVSL, JOBVSR, SORT
        #       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
        # *     ..
        # *     .. Array Arguments ..
        #       LOGICAL            BWORK( * )
        #       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
        #      $                   B( LDB, * ), BETA( * ), VSL( LDVSL, * ),
        #      $                   VSR( LDVSR, * ), WORK( * )
        function gges!(jobvsl::AbstractChar, jobvsr::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            @chkvalidparam 1 jobvsl ('N', 'V')
            @chkvalidparam 2 jobvsr ('N', 'V')
            chkstride1(A, B)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"dimensions of A, ($n,$n), and B, ($m,$m), must match"))
            end
            sdim = BlasInt(0)
            alphar = similar(A, $elty, n)
            alphai = similar(A, $elty, n)
            beta = similar(A, $elty, n)
            ldvsl = jobvsl == 'V' ? max(1, n) : 1
            vsl = similar(A, $elty, ldvsl, n)
            ldvsr = jobvsr == 'V' ? max(1, n) : 1
            vsr = similar(A, $elty, ldvsr, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gges), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{Cvoid},
                        Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                        Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty},
                        Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                        Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{Cvoid},
                        Ref{BlasInt}, Clong, Clong, Clong),
                    jobvsl, jobvsr, 'N', C_NULL,
                    n, A, max(1,stride(A, 2)), B,
                    max(1,stride(B, 2)), sdim, alphar, alphai,
                    beta, vsl, ldvsl, vsr,
                    ldvsr, work, lwork, C_NULL,
                    info, 1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, B, complex.(alphar, alphai), beta, vsl[1:(jobvsl == 'V' ? n : 0),:], vsr[1:(jobvsr == 'V' ? n : 0),:]
        end

        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVSL, JOBVSR, SORT
        #       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
        # *     ..
        # *     .. Array Arguments ..
        #       LOGICAL            BWORK( * )
        #       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
        #      $                   B( LDB, * ), BETA( * ), VSL( LDVSL, * ),
        #      $                   VSR( LDVSR, * ), WORK( * )
        function gges3!(jobvsl::AbstractChar, jobvsr::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            @chkvalidparam 1 jobvsl ('N', 'V')
            @chkvalidparam 2 jobvsr ('N', 'V')
            chkstride1(A, B)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"dimensions of A, ($n,$n), and B, ($m,$m), must match"))
            end
            sdim = BlasInt(0)
            alphar = similar(A, $elty, n)
            alphai = similar(A, $elty, n)
            beta = similar(A, $elty, n)
            ldvsl = jobvsl == 'V' ? max(1, n) : 1
            vsl = similar(A, $elty, ldvsl, n)
            ldvsr = jobvsr == 'V' ? max(1, n) : 1
            vsr = similar(A, $elty, ldvsr, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gges3), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{Cvoid},
                        Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                        Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty},
                        Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                        Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{Cvoid},
                        Ref{BlasInt}, Clong, Clong, Clong),
                    jobvsl, jobvsr, 'N', C_NULL,
                    n, A, max(1,stride(A, 2)), B,
                    max(1,stride(B, 2)), sdim, alphar, alphai,
                    beta, vsl, ldvsl, vsr,
                    ldvsr, work, lwork, C_NULL,
                    info, 1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, B, complex.(alphar, alphai), beta, vsl[1:(jobvsl == 'V' ? n : 0),:], vsr[1:(jobvsr == 'V' ? n : 0),:]
        end
    end
end

for (gees, gges, gges3, elty, relty) in
    ((:zgees_,:zgges_,:zgges3_,:ComplexF64,:Float64),
     (:cgees_,:cgges_,:cgges3_,:ComplexF32,:Float32))
    @eval begin
        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVS, SORT
        #       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
        # *     ..
        # *     .. Array Arguments ..
        #       LOGICAL            BWORK( * )
        #       DOUBLE PRECISION   RWORK( * )
        #       COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
        function gees!(jobvs::AbstractChar, A::AbstractMatrix{$elty})
            require_one_based_indexing(A)
            @chkvalidparam 1 jobvs ('N', 'V')
            chkstride1(A)
            n     = checksquare(A)
            sort  = 'N'
            sdim  = BlasInt(0)
            w     = similar(A, $elty, n)
            vs    = similar(A, $elty, jobvs == 'V' ? n : 1, n)
            ldvs  = max(size(vs, 1), 1)
            work  = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, n)
            info  = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gees), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ptr{Cvoid}, Ref{BlasInt},
                        Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                        Ptr{$relty}, Ptr{Cvoid}, Ref{BlasInt}, Clong, Clong),
                    jobvs, sort, C_NULL, n,
                        A, max(1, stride(A, 2)), sdim, w,
                        vs, ldvs, work, lwork,
                        rwork, C_NULL, info, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, vs, w
        end

        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVSL, JOBVSR, SORT
        #       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
        # *     ..
        # *     .. Array Arguments ..
        #       LOGICAL            BWORK( * )
        #       DOUBLE PRECISION   RWORK( * )
        #       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
        #      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
        #      $                   WORK( * )
        function gges!(jobvsl::AbstractChar, jobvsr::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            @chkvalidparam 1 jobvsl ('N', 'V')
            @chkvalidparam 2 jobvsr ('N', 'V')
            chkstride1(A, B)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"dimensions of A, ($n,$n), and B, ($m,$m), must match"))
            end
            sdim = BlasInt(0)
            alpha = similar(A, $elty, n)
            beta = similar(A, $elty, n)
            ldvsl = jobvsl == 'V' ? max(1, n) : 1
            vsl = similar(A, $elty, ldvsl, n)
            ldvsr = jobvsr == 'V' ? max(1, n) : 1
            vsr = similar(A, $elty, ldvsr, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, 8n)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gges), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{Cvoid},
                        Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                        Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{$relty}, Ptr{Cvoid},
                        Ref{BlasInt}, Clong, Clong, Clong),
                    jobvsl, jobvsr, 'N', C_NULL,
                    n, A, max(1, stride(A, 2)), B,
                    max(1, stride(B, 2)), sdim, alpha, beta,
                    vsl, ldvsl, vsr, ldvsr,
                    work, lwork, rwork, C_NULL,
                    info, 1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, B, alpha, beta, vsl[1:(jobvsl == 'V' ? n : 0),:], vsr[1:(jobvsr == 'V' ? n : 0),:]
        end

        # *     .. Scalar Arguments ..
        #       CHARACTER          JOBVSL, JOBVSR, SORT
        #       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
        # *     ..
        # *     .. Array Arguments ..
        #       LOGICAL            BWORK( * )
        #       DOUBLE PRECISION   RWORK( * )
        #       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
        #      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
        #      $                   WORK( * )
        function gges3!(jobvsl::AbstractChar, jobvsr::AbstractChar, A::AbstractMatrix{$elty}, B::AbstractMatrix{$elty})
            require_one_based_indexing(A, B)
            @chkvalidparam 1 jobvsl ('N', 'V')
            @chkvalidparam 2 jobvsr ('N', 'V')
            chkstride1(A, B)
            n, m = map(checksquare, (A, B))
            if n != m
                throw(DimensionMismatch(lazy"dimensions of A, ($n,$n), and B, ($m,$m), must match"))
            end
            sdim = BlasInt(0)
            alpha = similar(A, $elty, n)
            beta = similar(A, $elty, n)
            ldvsl = jobvsl == 'V' ? max(1, n) : 1
            vsl = similar(A, $elty, ldvsl, n)
            ldvsr = jobvsr == 'V' ? max(1, n) : 1
            vsr = similar(A, $elty, ldvsr, n)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            rwork = Vector{$relty}(undef, 8n)
            info = Ref{BlasInt}()
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($gges3), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{Cvoid},
                        Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                        Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty}, Ptr{$elty},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{$relty}, Ptr{Cvoid},
                        Ref{BlasInt}, Clong, Clong, Clong),
                    jobvsl, jobvsr, 'N', C_NULL,
                    n, A, max(1, stride(A, 2)), B,
                    max(1, stride(B, 2)), sdim, alpha, beta,
                    vsl, ldvsl, vsr, ldvsr,
                    work, lwork, rwork, C_NULL,
                    info, 1, 1, 1)
                chklapackerror(info[])
                if i == 1
                    lwork = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            A, B, alpha, beta, vsl[1:(jobvsl == 'V' ? n : 0),:], vsr[1:(jobvsr == 'V' ? n : 0),:]
        end
    end
end

"""
    gees!(jobvs, A) -> (A, vs, w)

Computes the eigenvalues (`jobvs = N`) or the eigenvalues and Schur
vectors (`jobvs = V`) of matrix `A`. `A` is overwritten by its Schur form.

Returns `A`, `vs` containing the Schur vectors, and `w`, containing the
eigenvalues.
"""
gees!(jobvs::AbstractChar, A::AbstractMatrix)


"""
    gges!(jobvsl, jobvsr, A, B) -> (A, B, alpha, beta, vsl, vsr)

Computes the generalized eigenvalues, generalized Schur form, left Schur
vectors (`jobsvl = V`), or right Schur vectors (`jobvsr = V`) of `A` and
`B`.

The generalized eigenvalues are returned in `alpha` and `beta`. The left Schur
vectors are returned in `vsl` and the right Schur vectors are returned in `vsr`.
"""
gges!(jobvsl::AbstractChar, jobvsr::AbstractChar, A::AbstractMatrix, B::AbstractMatrix)

"""
    gges3!(jobvsl, jobvsr, A, B) -> (A, B, alpha, beta, vsl, vsr)

Computes the generalized eigenvalues, generalized Schur form, left Schur
vectors (`jobsvl = V`), or right Schur vectors (`jobvsr = V`) of `A` and
`B` using a blocked algorithm. This function requires LAPACK 3.6.0.

The generalized eigenvalues are returned in `alpha` and `beta`. The left Schur
vectors are returned in `vsl` and the right Schur vectors are returned in `vsr`.
"""
gges3!(jobvsl::AbstractChar, jobvsr::AbstractChar, A::AbstractMatrix, B::AbstractMatrix)

for (trexc, trsen, tgsen, elty) in
    ((:dtrexc_, :dtrsen_, :dtgsen_, :Float64),
     (:strexc_, :strsen_, :stgsen_, :Float32))
    @eval begin
        # *     .. Scalar Arguments ..
        #       CHARACTER          COMPQ
        #       INTEGER            IFST, ILST, INFO, LDQ, LDT, N
        # *     ..
        # *     .. Array Arguments ..
        #       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
        function trexc!(compq::AbstractChar, ifst::BlasInt, ilst::BlasInt, T::AbstractMatrix{$elty}, Q::AbstractMatrix{$elty})
            require_one_based_indexing(T, Q)
            @chkvalidparam 1 compq ('V', 'N')
            chkstride1(T, Q)
            n = checksquare(T)
            ldt = max(1, stride(T, 2))
            ldq = max(1, stride(Q, 2))
            work = Vector{$elty}(undef, n)
            info = Ref{BlasInt}()
            ccall((@blasfunc($trexc), libblastrampoline), Cvoid,
                  (Ref{UInt8},  Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ref{BlasInt}, Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Clong),
                  compq, n,
                  T, ldt, Q, ldq,
                  ifst, ilst,
                  work, info, 1)
            chklapackerror(info[])
            T, Q
        end
        trexc!(ifst::BlasInt, ilst::BlasInt, T::AbstractMatrix{$elty}, Q::AbstractMatrix{$elty}) =
            trexc!('V', ifst, ilst, T, Q)

        # *     .. Scalar Arguments ..
        #       CHARACTER          COMPQ, JOB
        #       INTEGER            INFO, LDQ, LDT, LIWORK, LWORK, M, N
        #       DOUBLE PRECISION   S, SEP
        # *     ..
        # *     .. Array Arguments ..
        #       LOGICAL            SELECT( * )
        #       INTEGER            IWORK( * )
        #       DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ), WR( * )
        function trsen!(job::AbstractChar, compq::AbstractChar, select::AbstractVector{BlasInt},
                        T::AbstractMatrix{$elty}, Q::AbstractMatrix{$elty})
            require_one_based_indexing(T, Q, select)
            @chkvalidparam 1 job ('N', 'E', 'V', 'B')
            @chkvalidparam 2 compq ('V', 'N')
            chkstride1(T, Q, select)
            n = checksquare(T)
            ldt = max(1, stride(T, 2))
            ldq = max(1, stride(Q, 2))
            wr = similar(T, $elty, n)
            wi = similar(T, $elty, n)
            m = sum(select)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            iwork = Vector{BlasInt}(undef, 1)
            liwork = BlasInt(-1)
            info = Ref{BlasInt}()
            select = convert(Array{BlasInt}, select)
            s = Ref{$elty}(zero($elty))
            sep = Ref{$elty}(zero($elty))
            for i = 1:2  # first call returns lwork as work[1] and liwork as iwork[1]
                ccall((@blasfunc($trsen), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ptr{$elty}, Ref{BlasInt}, Ref{$elty}, Ref{$elty},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                    Ref{BlasInt}, Clong, Clong),
                    job, compq, select, n,
                    T, ldt, Q, ldq,
                    wr, wi, m, s, sep,
                    work, lwork, iwork, liwork,
                    info, 1, 1)
                chklapackerror(info[])
                if i == 1 # only estimated optimal lwork, liwork
                    lwork  = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    liwork = BlasInt(real(iwork[1]))
                    resize!(iwork, liwork)
                end
            end
            iszero(wi) ? (T, Q, wr, s[], sep[]) : (T, Q, complex.(wr, wi), s[], sep[])
        end
        trsen!(select::AbstractVector{BlasInt}, T::AbstractMatrix{$elty}, Q::AbstractMatrix{$elty}) =
            trsen!('N', 'V', select, T, Q)

        #        .. Scalar Arguments ..
        #        LOGICAL            WANTQ, WANTZ
        #        INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,
        #       $                   M, N
        #        DOUBLE PRECISION   PL, PR
        #        ..
        #        .. Array Arguments ..
        #        LOGICAL            SELECT( * )
        #        INTEGER            IWORK( * )
        #        DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
        #       $                   B( LDB, * ), BETA( * ), DIF( * ), Q( LDQ, * ),
        #       $                   WORK( * ), Z( LDZ, * )
        #        ..
        function tgsen!(select::AbstractVector{BlasInt}, S::AbstractMatrix{$elty}, T::AbstractMatrix{$elty},
                        Q::AbstractMatrix{$elty}, Z::AbstractMatrix{$elty})
            require_one_based_indexing(select, S, T, Q, Z)
            chkstride1(select, S, T, Q, Z)
            n, nt, nq, nz = checksquare(S, T, Q, Z)
            if n != nt
                throw(DimensionMismatch(lazy"dimensions of S, ($n,$n), and T, ($nt,$nt), must match"))
            end
            if n != nq
                throw(DimensionMismatch(lazy"dimensions of S, ($n,$n), and Q, ($nq,$nq), must match"))
            end
            if n != nz
                throw(DimensionMismatch(lazy"dimensions of S, ($n,$n), and Z, ($nz,$nz), must match"))
            end
            lds = max(1, stride(S, 2))
            ldt = max(1, stride(T, 2))
            ldq = max(1, stride(Q, 2))
            ldz = max(1, stride(Z, 2))
            m = sum(select)
            alphai = similar(T, $elty, n)
            alphar = similar(T, $elty, n)
            beta = similar(T, $elty, n)
            lwork = BlasInt(-1)
            work = Vector{$elty}(undef, 1)
            liwork = BlasInt(-1)
            iwork = Vector{BlasInt}(undef, 1)
            info = Ref{BlasInt}()
            select = convert(Array{BlasInt}, select)
            for i = 1:2  # first call returns lwork as work[1] and liwork as iwork[1]
                ccall((@blasfunc($tgsen), libblastrampoline), Cvoid,
                       (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                        Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                        Ref{BlasInt}, Ptr{$elty}, Ptr{$elty}, Ptr{$elty},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                        Ref{BlasInt}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                        Ptr{BlasInt}),
                    0, 1, 1, select,
                    n, S, lds, T,
                    ldt, alphar, alphai, beta,
                    Q, ldq, Z, ldz,
                    m, C_NULL, C_NULL, C_NULL,
                    work, lwork, iwork, liwork,
                    info)
                chklapackerror(info[])
                if i == 1 # only estimated optimal lwork, liwork
                    lwork  = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    liwork = BlasInt(real(iwork[1]))
                    resize!(iwork, liwork)
                end
            end
            S, T, complex.(alphar, alphai), beta, Q, Z
        end
    end
end

for (trexc, trsen, tgsen, elty, relty) in
    ((:ztrexc_, :ztrsen_, :ztgsen_, :ComplexF64, :Float64),
     (:ctrexc_, :ctrsen_, :ctgsen_, :ComplexF32, :Float32))
    @eval begin
        #      .. Scalar Arguments ..
        #      CHARACTER          COMPQ
        #      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
        #      ..
        #      .. Array Arguments ..
        #      DOUBLE PRECISION   Q( LDQ, * ), T( LDT, * ), WORK( * )
        function trexc!(compq::AbstractChar, ifst::BlasInt, ilst::BlasInt, T::AbstractMatrix{$elty}, Q::AbstractMatrix{$elty})
            require_one_based_indexing(T, Q)
            @chkvalidparam 1 compq ('V', 'N')
            chkstride1(T, Q)
            n = checksquare(T)
            ldt = max(1, stride(T, 2))
            ldq = max(1, stride(Q, 2))
            info = Ref{BlasInt}()
            ccall((@blasfunc($trexc), libblastrampoline), Cvoid,
                  (Ref{UInt8},  Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ref{BlasInt}, Ref{BlasInt},
                   Ref{BlasInt}, Clong),
                  compq, n,
                  T, ldt, Q, ldq,
                  ifst, ilst,
                  info,  1)
            chklapackerror(info[])
            T, Q
        end
        trexc!(ifst::BlasInt, ilst::BlasInt, T::AbstractMatrix{$elty}, Q::AbstractMatrix{$elty}) =
            trexc!('V', ifst, ilst, T, Q)

        #      .. Scalar Arguments ..
        #      CHARACTER          COMPQ, JOB
        #      INTEGER            INFO, LDQ, LDT, LWORK, M, N
        #      DOUBLE PRECISION   S, SEP
        #      ..
        #      .. Array Arguments ..
        #      LOGICAL            SELECT( * )
        #      COMPLEX            Q( LDQ, * ), T( LDT, * ), W( * ), WORK( * )
        function trsen!(job::AbstractChar, compq::AbstractChar, select::AbstractVector{BlasInt},
                        T::AbstractMatrix{$elty}, Q::AbstractMatrix{$elty})
            require_one_based_indexing(select, T, Q)
            @chkvalidparam 1 job ('N', 'E', 'V', 'B')
            @chkvalidparam 2 compq ('N', 'V')
            chkstride1(select, T, Q)
            n = checksquare(T)
            ldt = max(1, stride(T, 2))
            ldq = max(1, stride(Q, 2))
            w = similar(T, $elty, n)
            m = sum(select)
            work = Vector{$elty}(undef, 1)
            lwork = BlasInt(-1)
            info = Ref{BlasInt}()
            select = convert(Array{BlasInt}, select)
            s = Ref{$relty}(zero($relty))
            sep = Ref{$relty}(zero($relty))
            for i = 1:2  # first call returns lwork as work[1]
                ccall((@blasfunc($trsen), libblastrampoline), Cvoid,
                    (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                    Ptr{$elty}, Ref{BlasInt}, Ref{$relty}, Ref{$relty},
                    Ptr{$elty}, Ref{BlasInt},
                    Ref{BlasInt}, Clong, Clong),
                    job, compq, select, n,
                    T, ldt, Q, ldq,
                    w, m, s, sep,
                    work, lwork,
                    info, 1, 1)
                chklapackerror(info[])
                if i == 1 # only estimated optimal lwork, liwork
                    lwork  = BlasInt(real(work[1]))
                    resize!(work, lwork)
                end
            end
            T, Q, w, s[], sep[]
        end
        trsen!(select::AbstractVector{BlasInt}, T::AbstractMatrix{$elty}, Q::AbstractMatrix{$elty}) =
            trsen!('N', 'V', select, T, Q)

        #        .. Scalar Arguments ..
        #        LOGICAL            WANTQ, WANTZ
        #        INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,
        #       $                   M, N
        #        DOUBLE PRECISION   PL, PR
        #        ..
        #        .. Array Arguments ..
        #        LOGICAL            SELECT( * )
        #        INTEGER            IWORK( * )
        #        DOUBLE PRECISION   DIF( * )
        #        COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
        #       $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
        #        ..
        function tgsen!(select::AbstractVector{BlasInt}, S::AbstractMatrix{$elty}, T::AbstractMatrix{$elty},
                        Q::AbstractMatrix{$elty}, Z::AbstractMatrix{$elty})
            require_one_based_indexing(select, S, T, Q, Z)
            chkstride1(select, S, T, Q, Z)
            n, nt, nq, nz = checksquare(S, T, Q, Z)
            if n != nt
                throw(DimensionMismatch(lazy"dimensions of S, ($n,$n), and T, ($nt,$nt), must match"))
            end
            if n != nq
                throw(DimensionMismatch(lazy"dimensions of S, ($n,$n), and Q, ($nq,$nq), must match"))
            end
            if n != nz
                throw(DimensionMismatch(lazy"dimensions of S, ($n,$n), and Z, ($nz,$nz), must match"))
            end
            lds = max(1, stride(S, 2))
            ldt = max(1, stride(T, 2))
            ldq = max(1, stride(Q, 2))
            ldz = max(1, stride(Z, 2))
            m = sum(select)
            alpha = similar(T, $elty, n)
            beta = similar(T, $elty, n)
            lwork = BlasInt(-1)
            work = Vector{$elty}(undef, 1)
            liwork = BlasInt(-1)
            iwork = Vector{BlasInt}(undef, 1)
            info = Ref{BlasInt}()
            select = convert(Array{BlasInt}, select)
            for i = 1:2  # first call returns lwork as work[1] and liwork as iwork[1]
                ccall((@blasfunc($tgsen), libblastrampoline), Cvoid,
                       (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                        Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty},
                        Ref{BlasInt}, Ptr{$elty}, Ptr{$elty},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                        Ref{BlasInt}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid},
                        Ptr{$elty}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                        Ptr{BlasInt}),
                    0, 1, 1, select,
                    n, S, lds, T,
                    ldt, alpha, beta,
                    Q, ldq, Z, ldz,
                    m, C_NULL, C_NULL, C_NULL,
                    work, lwork, iwork, liwork,
                    info)
                chklapackerror(info[])
                if i == 1 # only estimated optimal lwork, liwork
                    lwork  = BlasInt(real(work[1]))
                    resize!(work, lwork)
                    liwork = BlasInt(real(iwork[1]))
                    resize!(iwork, liwork)
                end
            end
            S, T, alpha, beta, Q, Z
        end
    end
end

"""
    trexc!(compq, ifst, ilst, T, Q) -> (T, Q)
    trexc!(ifst, ilst, T, Q) -> (T, Q)

Reorder the Schur factorization `T` of a matrix, such that the diagonal block
of `T` with row index `ifst` is moved to row index `ilst`. If `compq = V`, the Schur
vectors `Q` are reordered. If `compq = N` they are not modified. The 4-arg method
calls the 5-arg method with `compq = V`.
"""
trexc!(compq::AbstractChar, ifst::BlasInt, ilst::BlasInt, T::AbstractMatrix, Q::AbstractMatrix)

"""
    trsen!(job, compq, select, T, Q) -> (T, Q, w, s, sep)
    trsen!(select, T, Q) -> (T, Q, w, s, sep)

Reorder the Schur factorization of a matrix and optionally finds reciprocal
condition numbers. If `job = N`, no condition numbers are found. If `job = E`,
only the condition number for this cluster of eigenvalues is found. If
`job = V`, only the condition number for the invariant subspace is found.
If `job = B` then the condition numbers for the cluster and subspace are
found. If `compq = V` the Schur vectors `Q` are updated. If `compq = N`
the Schur vectors are not modified. `select` determines which
eigenvalues are in the cluster. The 3-arg method calls the 5-arg method
with `job = N` and `compq = V`.

Returns `T`, `Q`, reordered eigenvalues in `w`, the condition number of the
cluster of eigenvalues `s`, and the condition number of the invariant subspace
`sep`.
"""
trsen!(compq::AbstractChar, job::AbstractChar, select::AbstractVector{BlasInt}, T::AbstractMatrix, Q::AbstractMatrix)

"""
    tgsen!(select, S, T, Q, Z) -> (S, T, alpha, beta, Q, Z)

Reorders the vectors of a generalized Schur decomposition. `select` specifies
the eigenvalues in each cluster.
"""
tgsen!(select::AbstractVector{BlasInt}, S::AbstractMatrix, T::AbstractMatrix, Q::AbstractMatrix, Z::AbstractMatrix)

for (fn, elty, relty) in ((:dtrsyl_, :Float64, :Float64),
                   (:strsyl_, :Float32, :Float32),
                   (:ztrsyl_, :ComplexF64, :Float64),
                   (:ctrsyl_, :ComplexF32, :Float32))
    @eval begin
        function trsyl!(transa::AbstractChar, transb::AbstractChar, A::AbstractMatrix{$elty},
                        B::AbstractMatrix{$elty}, C::AbstractMatrix{$elty}, isgn::Int=1)
            require_one_based_indexing(A, B, C)
            chktrans(transa)
            chktrans(transb)
            chkstride1(A, B, C)
            m, n = checksquare(A), checksquare(B)
            lda = max(1, stride(A, 2))
            ldb = max(1, stride(B, 2))
            m1, n1 = size(C)
            if m != m1 || n != n1
                throw(DimensionMismatch(lazy"dimensions of A, ($m,$n), and C, ($m1,$n1), must match"))
            end
            ldc = max(1, stride(C, 2))
            scale = Ref{$relty}()
            info  = Ref{BlasInt}()
            ccall((@blasfunc($fn), libblastrampoline), Cvoid,
                (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                 Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                 Ptr{$relty}, Ref{BlasInt}, Clong, Clong),
                transa, transb, isgn, m, n,
                A, lda, B, ldb, C, ldc,
                scale, info, 1, 1)
            chklapackerror(info[])
            C, scale[]
        end
    end
end

"""
    trsyl!(transa, transb, A, B, C, isgn=1) -> (C, scale)

Solves the Sylvester matrix equation `A * X +/- X * B = scale*C` where `A` and
`B` are both quasi-upper triangular. If `transa = N`, `A` is not modified.
If `transa = T`, `A` is transposed. If `transa = C`, `A` is conjugate
transposed. Similarly for `transb` and `B`. If `isgn = 1`, the equation
`A * X + X * B = scale * C` is solved. If `isgn = -1`, the equation
`A * X - X * B = scale * C` is solved.

Returns `X` (overwriting `C`) and `scale`.
"""
trsyl!(transa::AbstractChar, transb::AbstractChar, A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, isgn::Int=1)

for (fn, elty) in ((:dlacpy_, :Float64),
                   (:slacpy_, :Float32),
                   (:zlacpy_, :ComplexF64),
                   (:clacpy_, :ComplexF32))
    @eval begin
        # SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
        #     .. Scalar Arguments ..
        #      CHARACTER          UPLO
        #      INTEGER            LDA, LDB, M, N
        #     ..
        #     .. Array Arguments ..
        #     DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        #     ..
        function lacpy!(B::AbstractMatrix{$elty}, A::AbstractMatrix{$elty}, uplo::AbstractChar)
            require_one_based_indexing(A, B)
            chkstride1(A, B)
            m, n = size(A)
            m1, n1 = size(B)
            if uplo == 'U'
                lacpy_size_check((m1, n1), (n < m ? n : m, n))
            elseif uplo == 'L'
                lacpy_size_check((m1, n1), (m, m < n ? m : n))
            else
                lacpy_size_check((m1, n1), (m, n))
            end
            lda = max(1, stride(A, 2))
            ldb = max(1, stride(B, 2))
            ccall((@blasfunc($fn), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
                  Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt}, Clong),
                  uplo, m, n, A, lda, B, ldb, 1)
            B
        end
    end
end

# The noinline annotation reduces latency
@noinline lacpy_size_check((m1, n1), (m, n)) = (m1 < m || n1 < n) && throw(DimensionMismatch(lazy"B of size ($m1,$n1) should have at least size ($m,$n)"))

"""
    lacpy!(B, A, uplo) -> B

Copies all or part of a matrix `A` to another matrix `B`.
uplo specifies the part of the matrix `A` to be copied to `B`.
Set `uplo = 'L'` for the lower triangular part, `uplo = 'U'`
for the upper triangular part, any other character for all
the matrix `A`.

# Examples
```jldoctest
julia> A = [1. 2. ; 3. 4.];

julia> B = [0. 0. ; 0. 0.];

julia> LAPACK.lacpy!(B, A, 'U')
2×2 Matrix{Float64}:
 1.0  2.0
 0.0  4.0
```
"""
lacpy!(B::AbstractMatrix, A::AbstractMatrix, uplo::AbstractChar)

end # module
