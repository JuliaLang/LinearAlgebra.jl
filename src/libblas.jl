# function cdotc(n, cx, incx, cy, incy)
#     @ccall libblastrampoline.cdotc_(n::Ref{BlasInt}, cx::Ptr{ComplexF32},
#                                     incx::Ref{BlasInt},
#                                     cy::Ptr{ComplexF32}, incy::Ref{BlasInt})::ComplexF32
# end

# function cdotu(n, cx, incx, cy, incy)
#     @ccall libblastrampoline.cdotu_(n::Ref{BlasInt}, cx::Ptr{ComplexF32},
#                                     incx::Ref{BlasInt},
#                                     cy::Ptr{ComplexF32}, incy::Ref{BlasInt})::ComplexF32
# end

# function zdotc(n, zx, incx, zy, incy)
#     @ccall libblastrampoline.zdotc_(n::Ref{BlasInt}, zx::Ptr{ComplexF64},
#                                     incx::Ref{BlasInt},
#                                     zy::Ptr{ComplexF64}, incy::Ref{BlasInt})::ComplexF64
# end

# function zdotu(n, zx, incx, zy, incy)
#     @ccall libblastrampoline.zdotu_(n::Ref{BlasInt}, zx::Ptr{ComplexF64},
#                                     incx::Ref{BlasInt},
#                                     zy::Ptr{ComplexF64}, incy::Ref{BlasInt})::ComplexF64
# end

function srotg(a, b, c, s)
    @ccall libblastrampoline.srotg_(a::Ref{Float32}, b::Ref{Float32}, c::Ref{Float32},
                                    s::Ref{Float32})::Cvoid
end

function drotg(a, b, c, s)
    @ccall libblastrampoline.drotg_(a::Ref{Float64}, b::Ref{Float64}, c::Ref{Float64},
                                    s::Ref{Float64})::Cvoid
end

function crotg(a, b, c, s)
    @ccall libblastrampoline.crotg_(a::Ref{ComplexF32}, b::Ref{ComplexF32}, c::Ref{Float32},
                                    s::Ref{ComplexF32})::Cvoid
end

function zrotg(a, b, c, s)
    @ccall libblastrampoline.zrotg_(a::Ref{ComplexF64}, b::Ref{ComplexF64}, c::Ref{Float64},
                                    s::Ref{ComplexF64})::Cvoid
end

function snrm2(n, x, incx)
    @ccall libblastrampoline.snrm2_(n::Ref{BlasInt}, x::Ptr{Float32},
                                    incx::Ref{BlasInt})::Float32
end

function dnrm2(n, x, incx)
    @ccall libblastrampoline.dnrm2_(n::Ref{BlasInt}, x::Ptr{Float64},
                                    incx::Ref{BlasInt})::Float64
end

function scnrm2(n, x, incx)
    @ccall libblastrampoline.scnrm2_(n::Ref{BlasInt}, x::Ptr{ComplexF32},
                                     incx::Ref{BlasInt})::Float32
end

function dznrm2(n, x, incx)
    @ccall libblastrampoline.dznrm2_(n::Ref{BlasInt}, x::Ptr{ComplexF64},
                                     incx::Ref{BlasInt})::Float64
end

function caxpy(n, ca, cx, incx, cy, incy)
    @ccall libblastrampoline.caxpy_(n::Ref{BlasInt}, ca::Ref{ComplexF32},
                                    cx::Ptr{ComplexF32},
                                    incx::Ref{BlasInt}, cy::Ptr{ComplexF32},
                                    incy::Ref{BlasInt})::Cvoid
end

function ccopy(n, cx, incx, cy, incy)
    @ccall libblastrampoline.ccopy_(n::Ref{BlasInt}, cx::Ptr{ComplexF32},
                                    incx::Ref{BlasInt},
                                    cy::Ptr{ComplexF32}, incy::Ref{BlasInt})::Cvoid
end

function cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.cgbmv_(trans::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    kl::Ref{BlasInt},
                                    ku::Ref{BlasInt}, alpha::Ref{ComplexF32},
                                    a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt}, x::Ptr{ComplexF32},
                                    incx::Ref{BlasInt},
                                    beta::Ref{ComplexF32}, y::Ptr{ComplexF32},
                                    incy::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.cgemm_(transa::Ref{UInt8}, transb::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt}, k::Ref{BlasInt},
                                    alpha::Ref{ComplexF32},
                                    a::Ptr{ComplexF32}, lda::Ref{BlasInt},
                                    b::Ptr{ComplexF32},
                                    ldb::Ref{BlasInt}, beta::Ref{ComplexF32},
                                    c::Ptr{ComplexF32},
                                    ldc::Ref{BlasInt}, 1::Clong, 1::Clong)::Cvoid
end

function cgemmt(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.cgemmt_(uplo::Ref{UInt8}, transa::Ref{UInt8},
                                     transb::Ref{UInt8},
                                     n::Ref{BlasInt}, k::Ref{BlasInt},
                                     alpha::Ref{ComplexF32},
                                     a::Ptr{ComplexF32}, lda::Ref{BlasInt},
                                     b::Ptr{ComplexF32},
                                     ldb::Ref{BlasInt}, beta::Ref{ComplexF32},
                                     c::Ptr{ComplexF32},
                                     ldc::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function cgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.cgemmtr_(uplo::Ref{UInt8}, transa::Ref{UInt8},
                                      transb::Ref{UInt8},
                                      n::Ref{BlasInt}, k::Ref{BlasInt},
                                      alpha::Ref{ComplexF32},
                                      a::Ptr{ComplexF32}, lda::Ref{BlasInt},
                                      b::Ptr{ComplexF32},
                                      ldb::Ref{BlasInt}, beta::Ref{ComplexF32},
                                      c::Ptr{ComplexF32},
                                      ldc::Ref{BlasInt}, 1::Clong, 1::Clong,
                                      1::Clong)::Cvoid
end

function cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.cgemv_(trans::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt},
                                    beta::Ref{ComplexF32},
                                    y::Ptr{ComplexF32}, incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function cgerc(m, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.cgerc_(m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt},
                                    y::Ptr{ComplexF32},
                                    incy::Ref{BlasInt}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt})::Cvoid
end

function cgeru(m, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.cgeru_(m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt},
                                    y::Ptr{ComplexF32},
                                    incy::Ref{BlasInt}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt})::Cvoid
end

function chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.chbmv_(uplo::Ref{UInt8}, n::Ref{BlasInt}, k::Ref{BlasInt},
                                    alpha::Ref{ComplexF32}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt},
                                    beta::Ref{ComplexF32},
                                    y::Ptr{ComplexF32}, incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.chemm_(side::Ref{UInt8}, uplo::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    b::Ptr{ComplexF32}, ldb::Ref{BlasInt},
                                    beta::Ref{ComplexF32},
                                    c::Ptr{ComplexF32}, ldc::Ref{BlasInt}, 1::Clong,
                                    1::Clong)::Cvoid
end

function chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.chemv_(uplo::Ref{UInt8}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32},
                                    a::Ptr{ComplexF32}, lda::Ref{BlasInt},
                                    x::Ptr{ComplexF32},
                                    incx::Ref{BlasInt}, beta::Ref{ComplexF32},
                                    y::Ptr{ComplexF32},
                                    incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function cher(uplo, n, alpha, x, incx, a, lda)
    @ccall libblastrampoline.cher_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float32},
                                   x::Ptr{ComplexF32}, incx::Ref{BlasInt},
                                   a::Ptr{ComplexF32},
                                   lda::Ref{BlasInt}, 1::Clong)::Cvoid
end

function cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.cher2_(uplo::Ref{UInt8}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt},
                                    y::Ptr{ComplexF32},
                                    incy::Ref{BlasInt}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.cher2k_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                     k::Ref{BlasInt}, alpha::Ref{ComplexF32},
                                     a::Ptr{ComplexF32},
                                     lda::Ref{BlasInt}, b::Ptr{ComplexF32},
                                     ldb::Ref{BlasInt},
                                     beta::Ref{Float32}, c::Ptr{ComplexF32},
                                     ldc::Ref{BlasInt},
                                     1::Clong, 1::Clong)::Cvoid
end

function cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    @ccall libblastrampoline.cherk_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                    k::Ref{BlasInt},
                                    alpha::Ref{Float32}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    beta::Ref{Float32}, c::Ptr{ComplexF32},
                                    ldc::Ref{BlasInt},
                                    1::Clong, 1::Clong)::Cvoid
end

function chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    @ccall libblastrampoline.chpmv_(uplo::Ref{UInt8}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32},
                                    ap::Ptr{ComplexF32}, x::Ptr{ComplexF32},
                                    incx::Ref{BlasInt},
                                    beta::Ref{ComplexF32}, y::Ptr{ComplexF32},
                                    incy::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function chpr(uplo, n, alpha, x, incx, ap)
    @ccall libblastrampoline.chpr_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float32},
                                   x::Ptr{ComplexF32}, incx::Ref{BlasInt},
                                   ap::Ptr{ComplexF32},
                                   1::Clong)::Cvoid
end

function chpr2(uplo, n, alpha, x, incx, y, incy, ap)
    @ccall libblastrampoline.chpr2_(uplo::Ref{UInt8}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt},
                                    y::Ptr{ComplexF32},
                                    incy::Ref{BlasInt}, ap::Ptr{ComplexF32},
                                    1::Clong)::Cvoid
end

function cscal(n, ca, cx, incx)
    @ccall libblastrampoline.cscal_(n::Ref{BlasInt}, ca::Ref{ComplexF32},
                                    cx::Ptr{ComplexF32},
                                    incx::Ref{BlasInt})::Cvoid
end

function csrot(n, cx, incx, cy, incy, c, s)
    @ccall libblastrampoline.csrot_(n::Ref{BlasInt}, cx::Ptr{ComplexF32},
                                    incx::Ref{BlasInt},
                                    cy::Ptr{ComplexF32}, incy::Ref{BlasInt},
                                    c::Ref{Float32},
                                    s::Ref{Float32})::Cvoid
end

function csscal(n, sa, cx, incx)
    @ccall libblastrampoline.csscal_(n::Ref{BlasInt}, sa::Ref{Float32}, cx::Ptr{ComplexF32},
                                     incx::Ref{BlasInt})::Cvoid
end

function cswap(n, cx, incx, cy, incy)
    @ccall libblastrampoline.cswap_(n::Ref{BlasInt}, cx::Ptr{ComplexF32},
                                    incx::Ref{BlasInt},
                                    cy::Ptr{ComplexF32}, incy::Ref{BlasInt})::Cvoid
end

function csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.csymm_(side::Ref{UInt8}, uplo::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    b::Ptr{ComplexF32}, ldb::Ref{BlasInt},
                                    beta::Ref{ComplexF32},
                                    c::Ptr{ComplexF32}, ldc::Ref{BlasInt}, 1::Clong,
                                    1::Clong)::Cvoid
end

function csyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.csyr2k_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                     k::Ref{BlasInt}, alpha::Ref{ComplexF32},
                                     a::Ptr{ComplexF32},
                                     lda::Ref{BlasInt}, b::Ptr{ComplexF32},
                                     ldb::Ref{BlasInt},
                                     beta::Ref{ComplexF32}, c::Ptr{ComplexF32},
                                     ldc::Ref{BlasInt},
                                     1::Clong, 1::Clong)::Cvoid
end

function csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    @ccall libblastrampoline.csyrk_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                    k::Ref{BlasInt},
                                    alpha::Ref{ComplexF32}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    beta::Ref{ComplexF32}, c::Ptr{ComplexF32},
                                    ldc::Ref{BlasInt},
                                    1::Clong, 1::Clong)::Cvoid
end

function ctbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    @ccall libblastrampoline.ctbmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong)::Cvoid
end

function ctbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    @ccall libblastrampoline.ctbsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong)::Cvoid
end

function ctpmv(uplo, trans, diag, n, ap, x, incx)
    @ccall libblastrampoline.ctpmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, ap::Ptr{ComplexF32},
                                    x::Ptr{ComplexF32},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function ctpsv(uplo, trans, diag, n, ap, x, incx)
    @ccall libblastrampoline.ctpsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, ap::Ptr{ComplexF32},
                                    x::Ptr{ComplexF32},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    @ccall libblastrampoline.ctrmm_(side::Ref{UInt8}, uplo::Ref{UInt8}, transa::Ref{UInt8},
                                    diag::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    b::Ptr{ComplexF32}, ldb::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong, 1::Clong)::Cvoid
end

function ctrmv(uplo, trans, diag, n, a, lda, x, incx)
    @ccall libblastrampoline.ctrmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, a::Ptr{ComplexF32}, lda::Ref{BlasInt},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong)::Cvoid
end

function ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    @ccall libblastrampoline.ctrsm_(side::Ref{UInt8}, uplo::Ref{UInt8}, transa::Ref{UInt8},
                                    diag::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF32}, a::Ptr{ComplexF32},
                                    lda::Ref{BlasInt},
                                    b::Ptr{ComplexF32}, ldb::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong, 1::Clong)::Cvoid
end

function ctrsv(uplo, trans, diag, n, a, lda, x, incx)
    @ccall libblastrampoline.ctrsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, a::Ptr{ComplexF32}, lda::Ref{BlasInt},
                                    x::Ptr{ComplexF32}, incx::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong)::Cvoid
end

function dasum(n, dx, incx)
    @ccall libblastrampoline.dasum_(n::Ref{BlasInt}, dx::Ptr{Float64},
                                    incx::Ref{BlasInt})::Float64
end

function daxpy(n, da, dx, incx, dy, incy)
    @ccall libblastrampoline.daxpy_(n::Ref{BlasInt}, da::Ref{Float64}, dx::Ptr{Float64},
                                    incx::Ref{BlasInt}, dy::Ptr{Float64},
                                    incy::Ref{BlasInt})::Cvoid
end

function dcabs1(z)
    @ccall libblastrampoline.dcabs1_(z::Ref{ComplexF64})::Float64
end

function dcopy(n, dx, incx, dy, incy)
    @ccall libblastrampoline.dcopy_(n::Ref{BlasInt}, dx::Ptr{Float64}, incx::Ref{BlasInt},
                                    dy::Ptr{Float64}, incy::Ref{BlasInt})::Cvoid
end

function ddot(n, dx, incx, dy, incy)
    @ccall libblastrampoline.ddot_(n::Ref{BlasInt}, dx::Ptr{Float64}, incx::Ref{BlasInt},
                                   dy::Ptr{Float64}, incy::Ref{BlasInt})::Float64
end

function dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.dgbmv_(trans::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    kl::Ref{BlasInt},
                                    ku::Ref{BlasInt}, alpha::Ref{Float64}, a::Ptr{Float64},
                                    lda::Ref{BlasInt}, x::Ptr{Float64}, incx::Ref{BlasInt},
                                    beta::Ref{Float64}, y::Ptr{Float64}, incy::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.dgemm_(transa::Ref{UInt8}, transb::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, alpha::Ref{Float64},
                                    a::Ptr{Float64}, lda::Ref{BlasInt}, b::Ptr{Float64},
                                    ldb::Ref{BlasInt}, beta::Ref{Float64}, c::Ptr{Float64},
                                    ldc::Ref{BlasInt}, 1::Clong, 1::Clong)::Cvoid
end

function dgemmt(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.dgemmt_(uplo::Ref{UInt8}, transa::Ref{UInt8},
                                     transb::Ref{UInt8},
                                     n::Ref{BlasInt}, k::Ref{BlasInt}, alpha::Ref{Float64},
                                     a::Ptr{Float64}, lda::Ref{BlasInt}, b::Ptr{Float64},
                                     ldb::Ref{BlasInt}, beta::Ref{Float64}, c::Ptr{Float64},
                                     ldc::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function dgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.dgemmtr_(uplo::Ref{UInt8}, transa::Ref{UInt8},
                                      transb::Ref{UInt8},
                                      n::Ref{BlasInt}, k::Ref{BlasInt}, alpha::Ref{Float64},
                                      a::Ptr{Float64}, lda::Ref{BlasInt}, b::Ptr{Float64},
                                      ldb::Ref{BlasInt}, beta::Ref{Float64},
                                      c::Ptr{Float64},
                                      ldc::Ref{BlasInt}, 1::Clong, 1::Clong,
                                      1::Clong)::Cvoid
end

function dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.dgemv_(trans::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{Float64}, a::Ptr{Float64}, lda::Ref{BlasInt},
                                    x::Ptr{Float64}, incx::Ref{BlasInt}, beta::Ref{Float64},
                                    y::Ptr{Float64}, incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function dger(m, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.dger_(m::Ref{BlasInt}, n::Ref{BlasInt}, alpha::Ref{Float64},
                                   x::Ptr{Float64},
                                   incx::Ref{BlasInt}, y::Ptr{Float64}, incy::Ref{BlasInt},
                                   a::Ptr{Float64}, lda::Ref{BlasInt})::Cvoid
end

function drot(n, dx, incx, dy, incy, c, s)
    @ccall libblastrampoline.drot_(n::Ref{BlasInt}, dx::Ptr{Float64}, incx::Ref{BlasInt},
                                   dy::Ptr{Float64}, incy::Ref{BlasInt}, c::Ref{Float64},
                                   s::Ref{Float64})::Cvoid
end

function drotm(n, dx, incx, dy, incy, dparam)
    @ccall libblastrampoline.drotm_(n::Ref{BlasInt}, dx::Ptr{Float64}, incx::Ref{BlasInt},
                                    dy::Ptr{Float64}, incy::Ref{BlasInt},
                                    dparam::Ptr{Float64})::Cvoid
end

function drotmg(dd1, dd2, dx1, dy1, dparam)
    @ccall libblastrampoline.drotmg_(dd1::Ref{Float64}, dd2::Ref{Float64},
                                     dx1::Ref{Float64},
                                     dy1::Ref{Float64}, dparam::Ptr{Float64})::Cvoid
end

function dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.dsbmv_(uplo::Ref{UInt8}, n::Ref{BlasInt}, k::Ref{BlasInt},
                                    alpha::Ref{Float64}, a::Ptr{Float64}, lda::Ref{BlasInt},
                                    x::Ptr{Float64}, incx::Ref{BlasInt}, beta::Ref{Float64},
                                    y::Ptr{Float64}, incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function dscal(n, da, dx, incx)
    @ccall libblastrampoline.dscal_(n::Ref{BlasInt}, da::Ref{Float64}, dx::Ptr{Float64},
                                    incx::Ref{BlasInt})::Cvoid
end

function dsdot(n, sx, incx, sy, incy)
    @ccall libblastrampoline.dsdot_(n::Ref{BlasInt}, sx::Ptr{Float32}, incx::Ref{BlasInt},
                                    sy::Ptr{Float32}, incy::Ref{BlasInt})::Float64
end

function dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    @ccall libblastrampoline.dspmv_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float64},
                                    ap::Ptr{Float64}, x::Ptr{Float64}, incx::Ref{BlasInt},
                                    beta::Ref{Float64}, y::Ptr{Float64}, incy::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function dspr(uplo, n, alpha, x, incx, ap)
    @ccall libblastrampoline.dspr_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float64},
                                   x::Ptr{Float64}, incx::Ref{BlasInt}, ap::Ptr{Float64},
                                   1::Clong)::Cvoid
end

function dspr2(uplo, n, alpha, x, incx, y, incy, ap)
    @ccall libblastrampoline.dspr2_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float64},
                                    x::Ptr{Float64}, incx::Ref{BlasInt}, y::Ptr{Float64},
                                    incy::Ref{BlasInt}, ap::Ptr{Float64}, 1::Clong)::Cvoid
end

function dswap(n, dx, incx, dy, incy)
    @ccall libblastrampoline.dswap_(n::Ref{BlasInt}, dx::Ptr{Float64}, incx::Ref{BlasInt},
                                    dy::Ptr{Float64}, incy::Ref{BlasInt})::Cvoid
end

function dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.dsymm_(side::Ref{UInt8}, uplo::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt},
                                    alpha::Ref{Float64}, a::Ptr{Float64}, lda::Ref{BlasInt},
                                    b::Ptr{Float64}, ldb::Ref{BlasInt}, beta::Ref{Float64},
                                    c::Ptr{Float64}, ldc::Ref{BlasInt}, 1::Clong,
                                    1::Clong)::Cvoid
end

function dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.dsymv_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float64},
                                    a::Ptr{Float64}, lda::Ref{BlasInt}, x::Ptr{Float64},
                                    incx::Ref{BlasInt}, beta::Ref{Float64}, y::Ptr{Float64},
                                    incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function dsyr(uplo, n, alpha, x, incx, a, lda)
    @ccall libblastrampoline.dsyr_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float64},
                                   x::Ptr{Float64}, incx::Ref{BlasInt}, a::Ptr{Float64},
                                   lda::Ref{BlasInt}, 1::Clong)::Cvoid
end

function dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.dsyr2_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float64},
                                    x::Ptr{Float64}, incx::Ref{BlasInt}, y::Ptr{Float64},
                                    incy::Ref{BlasInt}, a::Ptr{Float64}, lda::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.dsyr2k_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                     k::Ref{BlasInt}, alpha::Ref{Float64}, a::Ptr{Float64},
                                     lda::Ref{BlasInt}, b::Ptr{Float64}, ldb::Ref{BlasInt},
                                     beta::Ref{Float64}, c::Ptr{Float64}, ldc::Ref{BlasInt},
                                     1::Clong, 1::Clong)::Cvoid
end

function dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    @ccall libblastrampoline.dsyrk_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                    k::Ref{BlasInt},
                                    alpha::Ref{Float64}, a::Ptr{Float64}, lda::Ref{BlasInt},
                                    beta::Ref{Float64}, c::Ptr{Float64}, ldc::Ref{BlasInt},
                                    1::Clong,
                                    1::Clong)::Cvoid
end

function dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    @ccall libblastrampoline.dtbmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, a::Ptr{Float64},
                                    lda::Ref{BlasInt},
                                    x::Ptr{Float64}, incx::Ref{BlasInt}, 1::Clong, 1::Clong,
                                    1::Clong)::Cvoid
end

function dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    @ccall libblastrampoline.dtbsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, a::Ptr{Float64},
                                    lda::Ref{BlasInt},
                                    x::Ptr{Float64}, incx::Ref{BlasInt}, 1::Clong, 1::Clong,
                                    1::Clong)::Cvoid
end

function dtpmv(uplo, trans, diag, n, ap, x, incx)
    @ccall libblastrampoline.dtpmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, ap::Ptr{Float64}, x::Ptr{Float64},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function dtpsv(uplo, trans, diag, n, ap, x, incx)
    @ccall libblastrampoline.dtpsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, ap::Ptr{Float64}, x::Ptr{Float64},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    @ccall libblastrampoline.dtrmm_(side::Ref{UInt8}, uplo::Ref{UInt8}, transa::Ref{UInt8},
                                    diag::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{Float64}, a::Ptr{Float64}, lda::Ref{BlasInt},
                                    b::Ptr{Float64}, ldb::Ref{BlasInt}, 1::Clong, 1::Clong,
                                    1::Clong, 1::Clong)::Cvoid
end

function dtrmv(uplo, trans, diag, n, a, lda, x, incx)
    @ccall libblastrampoline.dtrmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, a::Ptr{Float64}, lda::Ref{BlasInt},
                                    x::Ptr{Float64},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    @ccall libblastrampoline.dtrsm_(side::Ref{UInt8}, uplo::Ref{UInt8}, transa::Ref{UInt8},
                                    diag::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{Float64}, a::Ptr{Float64}, lda::Ref{BlasInt},
                                    b::Ptr{Float64}, ldb::Ref{BlasInt}, 1::Clong, 1::Clong,
                                    1::Clong, 1::Clong)::Cvoid
end

function dtrsv(uplo, trans, diag, n, a, lda, x, incx)
    @ccall libblastrampoline.dtrsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, a::Ptr{Float64}, lda::Ref{BlasInt},
                                    x::Ptr{Float64},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function dzasum(n, zx, incx)
    @ccall libblastrampoline.dzasum_(n::Ref{BlasInt}, zx::Ptr{ComplexF64},
                                     incx::Ref{BlasInt})::Float64
end

function icamax(n, cx, incx)
    @ccall libblastrampoline.icamax_(n::Ref{BlasInt}, cx::Ptr{ComplexF32},
                                     incx::Ref{BlasInt})::BlasInt
end

function idamax(n, dx, incx)
    @ccall libblastrampoline.idamax_(n::Ref{BlasInt}, dx::Ptr{Float64},
                                     incx::Ref{BlasInt})::BlasInt
end

function isamax(n, sx, incx)
    @ccall libblastrampoline.isamax_(n::Ref{BlasInt}, sx::Ptr{Float32},
                                     incx::Ref{BlasInt})::BlasInt
end

function izamax(n, zx, incx)
    @ccall libblastrampoline.izamax_(n::Ref{BlasInt}, zx::Ptr{ComplexF64},
                                     incx::Ref{BlasInt})::BlasInt
end

function sasum(n, sx, incx)
    @ccall libblastrampoline.sasum_(n::Ref{BlasInt}, sx::Ptr{Float32},
                                    incx::Ref{BlasInt})::Float32
end

function saxpy(n, sa, sx, incx, sy, incy)
    @ccall libblastrampoline.saxpy_(n::Ref{BlasInt}, sa::Ref{Float32}, sx::Ptr{Float32},
                                    incx::Ref{BlasInt}, sy::Ptr{Float32},
                                    incy::Ref{BlasInt})::Cvoid
end

function scabs1(z)
    @ccall libblastrampoline.scabs1_(z::Ref{ComplexF32})::Float32
end

function scasum(n, cx, incx)
    @ccall libblastrampoline.scasum_(n::Ref{BlasInt}, cx::Ptr{ComplexF32},
                                     incx::Ref{BlasInt})::Float32
end

function scopy(n, sx, incx, sy, incy)
    @ccall libblastrampoline.scopy_(n::Ref{BlasInt}, sx::Ptr{Float32}, incx::Ref{BlasInt},
                                    sy::Ptr{Float32}, incy::Ref{BlasInt})::Cvoid
end

function sdot(n, sx, incx, sy, incy)
    @ccall libblastrampoline.sdot_(n::Ref{BlasInt}, sx::Ptr{Float32}, incx::Ref{BlasInt},
                                   sy::Ptr{Float32}, incy::Ref{BlasInt})::Float32
end

function sdsdot(n, sb, sx, incx, sy, incy)
    @ccall libblastrampoline.sdsdot_(n::Ref{BlasInt}, sb::Ref{Float32}, sx::Ptr{Float32},
                                     incx::Ref{BlasInt}, sy::Ptr{Float32},
                                     incy::Ref{BlasInt})::Float32
end

function sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.sgbmv_(trans::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    kl::Ref{BlasInt},
                                    ku::Ref{BlasInt}, alpha::Ref{Float32}, a::Ptr{Float32},
                                    lda::Ref{BlasInt}, x::Ptr{Float32}, incx::Ref{BlasInt},
                                    beta::Ref{Float32}, y::Ptr{Float32}, incy::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.sgemm_(transa::Ref{UInt8}, transb::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, alpha::Ref{Float32},
                                    a::Ptr{Float32}, lda::Ref{BlasInt}, b::Ptr{Float32},
                                    ldb::Ref{BlasInt}, beta::Ref{Float32}, c::Ptr{Float32},
                                    ldc::Ref{BlasInt}, 1::Clong, 1::Clong)::Cvoid
end

function sgemmt(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.sgemmt_(uplo::Ref{UInt8}, transa::Ref{UInt8},
                                     transb::Ref{UInt8},
                                     n::Ref{BlasInt}, k::Ref{BlasInt}, alpha::Ref{Float32},
                                     a::Ptr{Float32}, lda::Ref{BlasInt}, b::Ptr{Float32},
                                     ldb::Ref{BlasInt}, beta::Ref{Float32}, c::Ptr{Float32},
                                     ldc::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function sgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.sgemmtr_(uplo::Ref{UInt8}, transa::Ref{UInt8},
                                      transb::Ref{UInt8},
                                      n::Ref{BlasInt}, k::Ref{BlasInt}, alpha::Ref{Float32},
                                      a::Ptr{Float32}, lda::Ref{BlasInt}, b::Ptr{Float32},
                                      ldb::Ref{BlasInt}, beta::Ref{Float32},
                                      c::Ptr{Float32},
                                      ldc::Ref{BlasInt}, 1::Clong, 1::Clong,
                                      1::Clong)::Cvoid
end

function sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.sgemv_(trans::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{Float32}, a::Ptr{Float32}, lda::Ref{BlasInt},
                                    x::Ptr{Float32}, incx::Ref{BlasInt}, beta::Ref{Float32},
                                    y::Ptr{Float32}, incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function sger(m, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.sger_(m::Ref{BlasInt}, n::Ref{BlasInt}, alpha::Ref{Float32},
                                   x::Ptr{Float32},
                                   incx::Ref{BlasInt}, y::Ptr{Float32}, incy::Ref{BlasInt},
                                   a::Ptr{Float32}, lda::Ref{BlasInt})::Cvoid
end

function srot(n, sx, incx, sy, incy, c, s)
    @ccall libblastrampoline.srot_(n::Ref{BlasInt}, sx::Ptr{Float32}, incx::Ref{BlasInt},
                                   sy::Ptr{Float32}, incy::Ref{BlasInt}, c::Ref{Float32},
                                   s::Ref{Float32})::Cvoid
end

function srotm(n, sx, incx, sy, incy, sparam)
    @ccall libblastrampoline.srotm_(n::Ref{BlasInt}, sx::Ptr{Float32}, incx::Ref{BlasInt},
                                    sy::Ptr{Float32}, incy::Ref{BlasInt},
                                    sparam::Ptr{Float32})::Cvoid
end

function srotmg(sd1, sd2, sx1, sy1, sparam)
    @ccall libblastrampoline.srotmg_(sd1::Ref{Float32}, sd2::Ref{Float32},
                                     sx1::Ref{Float32},
                                     sy1::Ref{Float32}, sparam::Ptr{Float32})::Cvoid
end

function ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.ssbmv_(uplo::Ref{UInt8}, n::Ref{BlasInt}, k::Ref{BlasInt},
                                    alpha::Ref{Float32}, a::Ptr{Float32}, lda::Ref{BlasInt},
                                    x::Ptr{Float32}, incx::Ref{BlasInt}, beta::Ref{Float32},
                                    y::Ptr{Float32}, incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function sscal(n, sa, sx, incx)
    @ccall libblastrampoline.sscal_(n::Ref{BlasInt}, sa::Ref{Float32}, sx::Ptr{Float32},
                                    incx::Ref{BlasInt})::Cvoid
end

function sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    @ccall libblastrampoline.sspmv_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float32},
                                    ap::Ptr{Float32}, x::Ptr{Float32}, incx::Ref{BlasInt},
                                    beta::Ref{Float32}, y::Ptr{Float32}, incy::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function sspr(uplo, n, alpha, x, incx, ap)
    @ccall libblastrampoline.sspr_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float32},
                                   x::Ptr{Float32}, incx::Ref{BlasInt}, ap::Ptr{Float32},
                                   1::Clong)::Cvoid
end

function sspr2(uplo, n, alpha, x, incx, y, incy, ap)
    @ccall libblastrampoline.sspr2_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float32},
                                    x::Ptr{Float32}, incx::Ref{BlasInt}, y::Ptr{Float32},
                                    incy::Ref{BlasInt}, ap::Ptr{Float32}, 1::Clong)::Cvoid
end

function sswap(n, sx, incx, sy, incy)
    @ccall libblastrampoline.sswap_(n::Ref{BlasInt}, sx::Ptr{Float32}, incx::Ref{BlasInt},
                                    sy::Ptr{Float32}, incy::Ref{BlasInt})::Cvoid
end

function ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.ssymm_(side::Ref{UInt8}, uplo::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt},
                                    alpha::Ref{Float32}, a::Ptr{Float32}, lda::Ref{BlasInt},
                                    b::Ptr{Float32}, ldb::Ref{BlasInt}, beta::Ref{Float32},
                                    c::Ptr{Float32}, ldc::Ref{BlasInt}, 1::Clong,
                                    1::Clong)::Cvoid
end

function ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.ssymv_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float32},
                                    a::Ptr{Float32}, lda::Ref{BlasInt}, x::Ptr{Float32},
                                    incx::Ref{BlasInt}, beta::Ref{Float32}, y::Ptr{Float32},
                                    incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function ssyr(uplo, n, alpha, x, incx, a, lda)
    @ccall libblastrampoline.ssyr_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float32},
                                   x::Ptr{Float32}, incx::Ref{BlasInt}, a::Ptr{Float32},
                                   lda::Ref{BlasInt}, 1::Clong)::Cvoid
end

function ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.ssyr2_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float32},
                                    x::Ptr{Float32}, incx::Ref{BlasInt}, y::Ptr{Float32},
                                    incy::Ref{BlasInt}, a::Ptr{Float32}, lda::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.ssyr2k_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                     k::Ref{BlasInt}, alpha::Ref{Float32}, a::Ptr{Float32},
                                     lda::Ref{BlasInt}, b::Ptr{Float32}, ldb::Ref{BlasInt},
                                     beta::Ref{Float32}, c::Ptr{Float32}, ldc::Ref{BlasInt},
                                     1::Clong, 1::Clong)::Cvoid
end

function ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    @ccall libblastrampoline.ssyrk_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                    k::Ref{BlasInt},
                                    alpha::Ref{Float32}, a::Ptr{Float32}, lda::Ref{BlasInt},
                                    beta::Ref{Float32}, c::Ptr{Float32}, ldc::Ref{BlasInt},
                                    1::Clong,
                                    1::Clong)::Cvoid
end

function stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    @ccall libblastrampoline.stbmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, a::Ptr{Float32},
                                    lda::Ref{BlasInt},
                                    x::Ptr{Float32}, incx::Ref{BlasInt}, 1::Clong, 1::Clong,
                                    1::Clong)::Cvoid
end

function stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    @ccall libblastrampoline.stbsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, a::Ptr{Float32},
                                    lda::Ref{BlasInt},
                                    x::Ptr{Float32}, incx::Ref{BlasInt}, 1::Clong, 1::Clong,
                                    1::Clong)::Cvoid
end

function stpmv(uplo, trans, diag, n, ap, x, incx)
    @ccall libblastrampoline.stpmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, ap::Ptr{Float32}, x::Ptr{Float32},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function stpsv(uplo, trans, diag, n, ap, x, incx)
    @ccall libblastrampoline.stpsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, ap::Ptr{Float32}, x::Ptr{Float32},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    @ccall libblastrampoline.strmm_(side::Ref{UInt8}, uplo::Ref{UInt8}, transa::Ref{UInt8},
                                    diag::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{Float32}, a::Ptr{Float32}, lda::Ref{BlasInt},
                                    b::Ptr{Float32}, ldb::Ref{BlasInt}, 1::Clong, 1::Clong,
                                    1::Clong, 1::Clong)::Cvoid
end

function strmv(uplo, trans, diag, n, a, lda, x, incx)
    @ccall libblastrampoline.strmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, a::Ptr{Float32}, lda::Ref{BlasInt},
                                    x::Ptr{Float32},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    @ccall libblastrampoline.strsm_(side::Ref{UInt8}, uplo::Ref{UInt8}, transa::Ref{UInt8},
                                    diag::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{Float32}, a::Ptr{Float32}, lda::Ref{BlasInt},
                                    b::Ptr{Float32}, ldb::Ref{BlasInt}, 1::Clong, 1::Clong,
                                    1::Clong, 1::Clong)::Cvoid
end

function strsv(uplo, trans, diag, n, a, lda, x, incx)
    @ccall libblastrampoline.strsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, a::Ptr{Float32}, lda::Ref{BlasInt},
                                    x::Ptr{Float32},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function zaxpy(n, za, zx, incx, zy, incy)
    @ccall libblastrampoline.zaxpy_(n::Ref{BlasInt}, za::Ref{ComplexF64},
                                    zx::Ptr{ComplexF64},
                                    incx::Ref{BlasInt}, zy::Ptr{ComplexF64},
                                    incy::Ref{BlasInt})::Cvoid
end

function zcopy(n, zx, incx, zy, incy)
    @ccall libblastrampoline.zcopy_(n::Ref{BlasInt}, zx::Ptr{ComplexF64},
                                    incx::Ref{BlasInt},
                                    zy::Ptr{ComplexF64}, incy::Ref{BlasInt})::Cvoid
end

function zdrot(n, zx, incx, zy, incy, c, s)
    @ccall libblastrampoline.zdrot_(n::Ref{BlasInt}, zx::Ptr{ComplexF64},
                                    incx::Ref{BlasInt},
                                    zy::Ptr{ComplexF64}, incy::Ref{BlasInt},
                                    c::Ref{Float64},
                                    s::Ref{Float64})::Cvoid
end

function zdscal(n, da, zx, incx)
    @ccall libblastrampoline.zdscal_(n::Ref{BlasInt}, da::Ref{Float64}, zx::Ptr{ComplexF64},
                                     incx::Ref{BlasInt})::Cvoid
end

function zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.zgbmv_(trans::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    kl::Ref{BlasInt},
                                    ku::Ref{BlasInt}, alpha::Ref{ComplexF64},
                                    a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt}, x::Ptr{ComplexF64},
                                    incx::Ref{BlasInt},
                                    beta::Ref{ComplexF64}, y::Ptr{ComplexF64},
                                    incy::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.zgemm_(transa::Ref{UInt8}, transb::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt}, k::Ref{BlasInt},
                                    alpha::Ref{ComplexF64},
                                    a::Ptr{ComplexF64}, lda::Ref{BlasInt},
                                    b::Ptr{ComplexF64},
                                    ldb::Ref{BlasInt}, beta::Ref{ComplexF64},
                                    c::Ptr{ComplexF64},
                                    ldc::Ref{BlasInt}, 1::Clong, 1::Clong)::Cvoid
end

function zgemmt(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.zgemmt_(uplo::Ref{UInt8}, transa::Ref{UInt8},
                                     transb::Ref{UInt8},
                                     n::Ref{BlasInt}, k::Ref{BlasInt},
                                     alpha::Ref{ComplexF64},
                                     a::Ptr{ComplexF64}, lda::Ref{BlasInt},
                                     b::Ptr{ComplexF64},
                                     ldb::Ref{BlasInt}, beta::Ref{ComplexF64},
                                     c::Ptr{ComplexF64},
                                     ldc::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function zgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.zgemmtr_(uplo::Ref{UInt8}, transa::Ref{UInt8},
                                      transb::Ref{UInt8},
                                      n::Ref{BlasInt}, k::Ref{BlasInt},
                                      alpha::Ref{ComplexF64},
                                      a::Ptr{ComplexF64}, lda::Ref{BlasInt},
                                      b::Ptr{ComplexF64},
                                      ldb::Ref{BlasInt}, beta::Ref{ComplexF64},
                                      c::Ptr{ComplexF64},
                                      ldc::Ref{BlasInt}, 1::Clong, 1::Clong,
                                      1::Clong)::Cvoid
end

function zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.zgemv_(trans::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt},
                                    beta::Ref{ComplexF64},
                                    y::Ptr{ComplexF64}, incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function zgerc(m, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.zgerc_(m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt},
                                    y::Ptr{ComplexF64},
                                    incy::Ref{BlasInt}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt})::Cvoid
end

function zgeru(m, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.zgeru_(m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt},
                                    y::Ptr{ComplexF64},
                                    incy::Ref{BlasInt}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt})::Cvoid
end

function zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.zhbmv_(uplo::Ref{UInt8}, n::Ref{BlasInt}, k::Ref{BlasInt},
                                    alpha::Ref{ComplexF64}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt},
                                    beta::Ref{ComplexF64},
                                    y::Ptr{ComplexF64}, incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.zhemm_(side::Ref{UInt8}, uplo::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    b::Ptr{ComplexF64}, ldb::Ref{BlasInt},
                                    beta::Ref{ComplexF64},
                                    c::Ptr{ComplexF64}, ldc::Ref{BlasInt}, 1::Clong,
                                    1::Clong)::Cvoid
end

function zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    @ccall libblastrampoline.zhemv_(uplo::Ref{UInt8}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64},
                                    a::Ptr{ComplexF64}, lda::Ref{BlasInt},
                                    x::Ptr{ComplexF64},
                                    incx::Ref{BlasInt}, beta::Ref{ComplexF64},
                                    y::Ptr{ComplexF64},
                                    incy::Ref{BlasInt}, 1::Clong)::Cvoid
end

function zher(uplo, n, alpha, x, incx, a, lda)
    @ccall libblastrampoline.zher_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float64},
                                   x::Ptr{ComplexF64}, incx::Ref{BlasInt},
                                   a::Ptr{ComplexF64},
                                   lda::Ref{BlasInt}, 1::Clong)::Cvoid
end

function zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    @ccall libblastrampoline.zher2_(uplo::Ref{UInt8}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt},
                                    y::Ptr{ComplexF64},
                                    incy::Ref{BlasInt}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.zher2k_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                     k::Ref{BlasInt}, alpha::Ref{ComplexF64},
                                     a::Ptr{ComplexF64},
                                     lda::Ref{BlasInt}, b::Ptr{ComplexF64},
                                     ldb::Ref{BlasInt},
                                     beta::Ref{Float64}, c::Ptr{ComplexF64},
                                     ldc::Ref{BlasInt},
                                     1::Clong, 1::Clong)::Cvoid
end

function zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    @ccall libblastrampoline.zherk_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                    k::Ref{BlasInt},
                                    alpha::Ref{Float64}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    beta::Ref{Float64}, c::Ptr{ComplexF64},
                                    ldc::Ref{BlasInt},
                                    1::Clong, 1::Clong)::Cvoid
end

function zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    @ccall libblastrampoline.zhpmv_(uplo::Ref{UInt8}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64},
                                    ap::Ptr{ComplexF64}, x::Ptr{ComplexF64},
                                    incx::Ref{BlasInt},
                                    beta::Ref{ComplexF64}, y::Ptr{ComplexF64},
                                    incy::Ref{BlasInt},
                                    1::Clong)::Cvoid
end

function zhpr(uplo, n, alpha, x, incx, ap)
    @ccall libblastrampoline.zhpr_(uplo::Ref{UInt8}, n::Ref{BlasInt}, alpha::Ref{Float64},
                                   x::Ptr{ComplexF64}, incx::Ref{BlasInt},
                                   ap::Ptr{ComplexF64},
                                   1::Clong)::Cvoid
end

function zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
    @ccall libblastrampoline.zhpr2_(uplo::Ref{UInt8}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt},
                                    y::Ptr{ComplexF64},
                                    incy::Ref{BlasInt}, ap::Ptr{ComplexF64},
                                    1::Clong)::Cvoid
end

function zscal(n, za, zx, incx)
    @ccall libblastrampoline.zscal_(n::Ref{BlasInt}, za::Ref{ComplexF64},
                                    zx::Ptr{ComplexF64},
                                    incx::Ref{BlasInt})::Cvoid
end

function zswap(n, zx, incx, zy, incy)
    @ccall libblastrampoline.zswap_(n::Ref{BlasInt}, zx::Ptr{ComplexF64},
                                    incx::Ref{BlasInt},
                                    zy::Ptr{ComplexF64}, incy::Ref{BlasInt})::Cvoid
end

function zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.zsymm_(side::Ref{UInt8}, uplo::Ref{UInt8}, m::Ref{BlasInt},
                                    n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    b::Ptr{ComplexF64}, ldb::Ref{BlasInt},
                                    beta::Ref{ComplexF64},
                                    c::Ptr{ComplexF64}, ldc::Ref{BlasInt}, 1::Clong,
                                    1::Clong)::Cvoid
end

function zsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    @ccall libblastrampoline.zsyr2k_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                     k::Ref{BlasInt}, alpha::Ref{ComplexF64},
                                     a::Ptr{ComplexF64},
                                     lda::Ref{BlasInt}, b::Ptr{ComplexF64},
                                     ldb::Ref{BlasInt},
                                     beta::Ref{ComplexF64}, c::Ptr{ComplexF64},
                                     ldc::Ref{BlasInt},
                                     1::Clong, 1::Clong)::Cvoid
end

function zsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    @ccall libblastrampoline.zsyrk_(uplo::Ref{UInt8}, trans::Ref{UInt8}, n::Ref{BlasInt},
                                    k::Ref{BlasInt},
                                    alpha::Ref{ComplexF64}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    beta::Ref{ComplexF64}, c::Ptr{ComplexF64},
                                    ldc::Ref{BlasInt},
                                    1::Clong, 1::Clong)::Cvoid
end

function ztbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    @ccall libblastrampoline.ztbmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong)::Cvoid
end

function ztbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    @ccall libblastrampoline.ztbsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, k::Ref{BlasInt}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong)::Cvoid
end

function ztpmv(uplo, trans, diag, n, ap, x, incx)
    @ccall libblastrampoline.ztpmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, ap::Ptr{ComplexF64},
                                    x::Ptr{ComplexF64},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function ztpsv(uplo, trans, diag, n, ap, x, incx)
    @ccall libblastrampoline.ztpsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, ap::Ptr{ComplexF64},
                                    x::Ptr{ComplexF64},
                                    incx::Ref{BlasInt}, 1::Clong, 1::Clong, 1::Clong)::Cvoid
end

function ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    @ccall libblastrampoline.ztrmm_(side::Ref{UInt8}, uplo::Ref{UInt8}, transa::Ref{UInt8},
                                    diag::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    b::Ptr{ComplexF64}, ldb::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong, 1::Clong)::Cvoid
end

function ztrmv(uplo, trans, diag, n, a, lda, x, incx)
    @ccall libblastrampoline.ztrmv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, a::Ptr{ComplexF64}, lda::Ref{BlasInt},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong)::Cvoid
end

function ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    @ccall libblastrampoline.ztrsm_(side::Ref{UInt8}, uplo::Ref{UInt8}, transa::Ref{UInt8},
                                    diag::Ref{UInt8}, m::Ref{BlasInt}, n::Ref{BlasInt},
                                    alpha::Ref{ComplexF64}, a::Ptr{ComplexF64},
                                    lda::Ref{BlasInt},
                                    b::Ptr{ComplexF64}, ldb::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong, 1::Clong)::Cvoid
end

function ztrsv(uplo, trans, diag, n, a, lda, x, incx)
    @ccall libblastrampoline.ztrsv_(uplo::Ref{UInt8}, trans::Ref{UInt8}, diag::Ref{UInt8},
                                    n::Ref{BlasInt}, a::Ptr{ComplexF64}, lda::Ref{BlasInt},
                                    x::Ptr{ComplexF64}, incx::Ref{BlasInt}, 1::Clong,
                                    1::Clong,
                                    1::Clong)::Cvoid
end
