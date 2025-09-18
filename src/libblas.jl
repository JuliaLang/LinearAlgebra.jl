function ssymv(uplo, n, alpha, A, lda, px, stx, beta, py, sty)
    return ccall((@blasfunc(ssymv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Clong),
                 uplo, n, alpha, A, lda, px, stx, beta, py, sty, 1)
end

function dsymv(uplo, n, alpha, A, lda, px, stx, beta, py, sty)
    return ccall((@blasfunc(dsymv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Clong),
                 uplo, n, alpha, A, lda, px, stx, beta, py, sty, 1)
end

function csymv(uplo, n, alpha, A, lda, px, stx, beta, py, sty)
    return ccall((@blasfunc(csymv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong),
                 uplo, n, alpha, A, lda, px, stx, beta, py, sty, 1)
end

function zsymv(uplo, n, alpha, A, lda, px, stx, beta, py, sty)
    return ccall((@blasfunc(zsymv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong),
                 uplo, n, alpha, A, lda, px, stx, beta, py, sty, 1)
end

function ssyr(uplo, n, a, x, stx, A, lda)
    return ccall((@blasfunc(ssyr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong),
                 uplo, n, a, x, stx, A, lda, 1)
end

function dsyr(uplo, n, a, x, stx, A, lda)
    return ccall((@blasfunc(dsyr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong),
                 uplo, n, a, x, stx, A, lda, 1)
end

function csyr(uplo, n, a, x, stx, A, lda)
    return ccall((@blasfunc(csyr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong),
                 uplo, n, a, x, stx, A, lda, 1)
end

function zsyr(uplo, n, a, x, stx, A, lda)
    return ccall((@blasfunc(zsyr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong),
                 uplo, n, a, x, stx, A, lda, 1)
end

function drot(n, dx, incx, dy, incy, c, s)
    return ccall((@blasfunc(drot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}),
                 n, dx, incx, dy, incy, c, s)
end

function srot(n, sx, incx, sy, incy, c, s)
    return ccall((@blasfunc(srot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}),
                 n, sx, incx, sy, incy, c, s)
end

function zdrot(n, zx, incx, zy, incy, c, s)
    return ccall((@blasfunc(zdrot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}),
                 n, zx, incx, zy, incy, c, s)
end

function csrot(n, cx, incx, cy, incy, c, s)
    return ccall((@blasfunc(csrot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}),
                 n, cx, incx, cy, incy, c, s)
end

function zrot(n, zx, incx, zy, incy, c, s)
    return ccall((@blasfunc(zrot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{ComplexF64}),
                 n, zx, incx, zy, incy, c, s)
end

function crot(n, cx, incx, cy, incy, c, s)
    return ccall((@blasfunc(crot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{ComplexF32}),
                 n, cx, incx, cy, incy, c, s)
end

function srotg(a, b, c, s)
    return ccall((@blasfunc(srotg_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), a, b, c, s)
end

function drotg(a, b, c, s)
    return ccall((@blasfunc(drotg_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), a, b, c, s)
end

function crotg(a, b, c, s)
    return ccall((@blasfunc(crotg_), libblastrampoline), Cvoid,
                 (Ref{ComplexF32}, Ref{ComplexF32}, Ref{Float32}, Ref{ComplexF32}), a, b, c, s)
end

function zrotg(a, b, c, s)
    return ccall((@blasfunc(zrotg_), libblastrampoline), Cvoid,
                 (Ref{ComplexF64}, Ref{ComplexF64}, Ref{Float64}, Ref{ComplexF64}), a, b, c, s)
end

function snrm2(n, x, incx)
    return ccall((@blasfunc(snrm2_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), n, x, incx)
end

function dnrm2(n, x, incx)
    return ccall((@blasfunc(dnrm2_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, x, incx)
end

function scnrm2(n, x, incx)
    return ccall((@blasfunc(scnrm2_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, x, incx)
end

function dznrm2(n, x, incx)
    return ccall((@blasfunc(dznrm2_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, x, incx)
end

function saxpby(n, sa, sx, incx, sb, sy, incy)
    return ccall((@blasfunc(saxpby_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}), n, sa, sx, incx, sb, sy, incy)
end

function daxpby(n, da, dx, incx, db, dy, incy)
    return ccall((@blasfunc(daxpby_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}), n, da, dx, incx, db, dy, incy)
end

function caxpby(n, ca, cx, incx, cb, cy, incy)
    return ccall((@blasfunc(caxpby_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}), n, ca, cx, incx, cb, cy, incy)
end

function zaxpby(n, za, zx, incx, zb, zy, incy)
    return ccall((@blasfunc(zaxpby_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}), n, za, zx, incx, zb, zy, incy)
end

function sgemmt(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(sgemmt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1, 1)
end

function dgemmt(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(dgemmt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1, 1)
end

function cgemmt(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(cgemmt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1, 1)
end

function zgemmt(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(zgemmt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1, 1)
end

function caxpy(n, ca, cx, incx, cy, incy)
    return ccall((@blasfunc(caxpy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}), n, ca, cx, incx, cy, incy)
end

function ccopy(n, cx, incx, cy, incy)
    return ccall((@blasfunc(ccopy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, cx,
                 incx, cy, incy)
end

function cdotc(n, cx, incx, cy, incy)
    return ccall((@blasfunc(cdotc_), libblastrampoline), ComplexF32,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, cx,
                 incx, cy, incy)
end

function cdotu(n, cx, incx, cy, incy)
    return ccall((@blasfunc(cdotu_), libblastrampoline), ComplexF32,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, cx,
                 incx, cy, incy)
end

function cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(cgbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), trans, m, n, kl, ku, alpha, a, lda, x,
                 incx, beta, y, incy, 1)
end

function cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(cgemm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), transa, transb, m, n, k, alpha,
                 a, lda, b, ldb, beta, c, ldc, 1, 1)
end

function cgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(cgemmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, transa, transb, n,
                 k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1, 1)
end

function cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(cgemv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong), trans, m, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function cgerc(m, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(cgerc_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, alpha, x,
                 incx, y, incy, a, lda)
end

function cgeru(m, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(cgeru_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, alpha, x,
                 incx, y, incy, a, lda)
end

function chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(chbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong), uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(chemm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), side, uplo, m, n, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(chemv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Clong), uplo, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function cher(uplo, n, alpha, x, incx, a, lda)
    return ccall((@blasfunc(cher_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, alpha, x, incx, a, lda, 1)
end

function cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(cher2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 alpha, x, incx, y, incy, a, lda, 1)
end

function cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(cher2k_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), uplo, trans, n, k, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    return ccall((@blasfunc(cherk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Clong, Clong), uplo, trans, n, k, alpha, a, lda, beta, c, ldc, 1, 1)
end

function chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    return ccall((@blasfunc(chpmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 alpha, ap, x, incx, beta, y, incy, 1)
end

function chpr(uplo, n, alpha, x, incx, ap)
    return ccall((@blasfunc(chpr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Clong), uplo, n, alpha, x, incx, ap, 1)
end

function chpr2(uplo, n, alpha, x, incx, y, incy, ap)
    return ccall((@blasfunc(chpr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong), uplo, n, alpha, x,
                 incx, y, incy, ap, 1)
end

function cscal(n, ca, cx, incx)
    return ccall((@blasfunc(cscal_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}), n, ca, cx, incx)
end

function csrot(n, cx, incx, cy, incy, c, s)
    return ccall((@blasfunc(csrot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}), n, cx, incx, cy, incy, c, s)
end

function csscal(n, sa, cx, incx)
    return ccall((@blasfunc(csscal_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}), n, sa, cx, incx)
end

function cswap(n, cx, incx, cy, incy)
    return ccall((@blasfunc(cswap_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, cx,
                 incx, cy, incy)
end

function csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(csymm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), side, uplo, m, n, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function csyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(csyr2k_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), uplo, trans, n, k, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    return ccall((@blasfunc(csyrk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Clong, Clong), uplo, trans, n, k, alpha, a, lda, beta, c, ldc, 1, 1)
end

function ctbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    return ccall((@blasfunc(ctbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, k, a, lda, x, incx, 1, 1, 1)
end

function ctbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    return ccall((@blasfunc(ctbsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, k, a, lda, x, incx, 1, 1, 1)
end

function ctpmv(uplo, trans, diag, n, ap, x, incx)
    return ccall((@blasfunc(ctpmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n,
                 ap, x, incx, 1, 1, 1)
end

function ctpsv(uplo, trans, diag, n, ap, x, incx)
    return ccall((@blasfunc(ctpsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n,
                 ap, x, incx, 1, 1, 1)
end

function ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    return ccall((@blasfunc(ctrmm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), side, uplo, transa, diag, m, n, alpha, a,
                 lda, b, ldb, 1, 1, 1, 1)
end

function ctrmv(uplo, trans, diag, n, a, lda, x, incx)
    return ccall((@blasfunc(ctrmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, a,
                 lda, x, incx, 1, 1, 1)
end

function ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    return ccall((@blasfunc(ctrsm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), side, uplo, transa, diag, m, n, alpha, a,
                 lda, b, ldb, 1, 1, 1, 1)
end

function ctrsv(uplo, trans, diag, n, a, lda, x, incx)
    return ccall((@blasfunc(ctrsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, a,
                 lda, x, incx, 1, 1, 1)
end

function dasum(n, dx, incx)
    return ccall((@blasfunc(dasum_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, dx, incx)
end

function daxpy(n, da, dx, incx, dy, incy)
    return ccall((@blasfunc(daxpy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}), n, da, dx, incx, dy, incy)
end

function dcabs1(z)
    return ccall((@blasfunc(dcabs1_), libblastrampoline), Float64, (Ref{ComplexF64},), z)
end

function dcopy(n, dx, incx, dy, incy)
    return ccall((@blasfunc(dcopy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, dx, incx,
                 dy, incy)
end

function ddot(n, dx, incx, dy, incy)
    return ccall((@blasfunc(ddot_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, dx, incx,
                 dy, incy)
end

function dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(dgbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Clong), trans, m, n, kl, ku, alpha, a, lda, x,
                 incx, beta, y, incy, 1)
end

function dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(dgemm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), transa, transb, m, n, k, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function dgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(dgemmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, transa, transb, n, k,
                 alpha, a, lda, b, ldb, beta, c, ldc, 1, 1, 1)
end

function dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(dgemv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong),
                 trans, m, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function dger(m, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(dger_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), m, n, alpha, x, incx, y, incy, a,
                 lda)
end

function drot(n, dx, incx, dy, incy, c, s)
    return ccall((@blasfunc(drot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}), n, dx, incx, dy, incy, c, s)
end

function drotm(n, dx, incx, dy, incy, dparam)
    return ccall((@blasfunc(drotm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}), n, dx, incx, dy, incy, dparam)
end

function drotmg(dd1, dd2, dx1, dy1, dparam)
    return ccall((@blasfunc(drotmg_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64}),
                 dd1, dd2, dx1, dy1, dparam)
end

function dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(dsbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong),
                 uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function dscal(n, da, dx, incx)
    return ccall((@blasfunc(dscal_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}), n, da, dx, incx)
end

function dsdot(n, sx, incx, sy, incy)
    return ccall((@blasfunc(dsdot_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), n, sx, incx,
                 sy, incy)
end

function dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    return ccall((@blasfunc(dspmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, alpha,
                 ap, x, incx, beta, y, incy, 1)
end

function dspr(uplo, n, alpha, x, incx, ap)
    return ccall((@blasfunc(dspr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Clong), uplo, n, alpha, x, incx, ap, 1)
end

function dspr2(uplo, n, alpha, x, incx, y, incy, ap)
    return ccall((@blasfunc(dspr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong), uplo, n, alpha, x, incx, y,
                 incy, ap, 1)
end

function dswap(n, dx, incx, dy, incy)
    return ccall((@blasfunc(dswap_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, dx, incx,
                 dy, incy)
end

function dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(dsymm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, 1,
                 1)
end

function dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(dsymv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong),
                 uplo, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function dsyr(uplo, n, alpha, x, incx, a, lda)
    return ccall((@blasfunc(dsyr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, alpha, x, incx, a, lda, 1)
end

function dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(dsyr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, alpha,
                 x, incx, y, incy, a, lda, 1)
end

function dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(dsyr2k_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1,
                 1)
end

function dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    return ccall((@blasfunc(dsyrk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), uplo,
                 trans, n, k, alpha, a, lda, beta, c, ldc, 1, 1)
end

function dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    return ccall((@blasfunc(dtbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, k, a, lda, x, incx, 1, 1, 1)
end

function dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    return ccall((@blasfunc(dtbsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, k, a, lda, x, incx, 1, 1, 1)
end

function dtpmv(uplo, trans, diag, n, ap, x, incx)
    return ccall((@blasfunc(dtpmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, ap, x, incx, 1, 1,
                 1)
end

function dtpsv(uplo, trans, diag, n, ap, x, incx)
    return ccall((@blasfunc(dtpsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, ap, x, incx, 1, 1,
                 1)
end

function dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    return ccall((@blasfunc(dtrmm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong), side, uplo, transa, diag, m, n, alpha, a, lda, b,
                 ldb, 1, 1, 1, 1)
end

function dtrmv(uplo, trans, diag, n, a, lda, x, incx)
    return ccall((@blasfunc(dtrmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, a,
                 lda, x, incx, 1, 1, 1)
end

function dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    return ccall((@blasfunc(dtrsm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong), side, uplo, transa, diag, m, n, alpha, a, lda, b,
                 ldb, 1, 1, 1, 1)
end

function dtrsv(uplo, trans, diag, n, a, lda, x, incx)
    return ccall((@blasfunc(dtrsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, a,
                 lda, x, incx, 1, 1, 1)
end

function dzasum(n, zx, incx)
    return ccall((@blasfunc(dzasum_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, zx, incx)
end

function icamax(n, cx, incx)
    return ccall((@blasfunc(icamax_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, cx, incx)
end

function idamax(n, dx, incx)
    return ccall((@blasfunc(idamax_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, dx, incx)
end

function isamax(n, sx, incx)
    return ccall((@blasfunc(isamax_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), n, sx, incx)
end

function izamax(n, zx, incx)
    return ccall((@blasfunc(izamax_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, zx, incx)
end

function sasum(n, sx, incx)
    return ccall((@blasfunc(sasum_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), n, sx, incx)
end

function saxpy(n, sa, sx, incx, sy, incy)
    return ccall((@blasfunc(saxpy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}), n, sa, sx, incx, sy, incy)
end

function scabs1(z)
    return ccall((@blasfunc(scabs1_), libblastrampoline), Float32, (Ref{ComplexF32},), z)
end

function scasum(n, cx, incx)
    return ccall((@blasfunc(scasum_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, cx, incx)
end

function scopy(n, sx, incx, sy, incy)
    return ccall((@blasfunc(scopy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), n, sx, incx,
                 sy, incy)
end

function sdot(n, sx, incx, sy, incy)
    return ccall((@blasfunc(sdot_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), n, sx, incx,
                 sy, incy)
end

function sdsdot(n, sb, sx, incx, sy, incy)
    return ccall((@blasfunc(sdsdot_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}), n, sb, sx, incx, sy, incy)
end

function sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(sgbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Clong), trans, m, n, kl, ku, alpha, a, lda, x,
                 incx, beta, y, incy, 1)
end

function sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(sgemm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), transa, transb, m, n, k, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function sgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(sgemmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, transa, transb, n, k,
                 alpha, a, lda, b, ldb, beta, c, ldc, 1, 1, 1)
end

function sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(sgemv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong),
                 trans, m, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function sger(m, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(sger_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), m, n, alpha, x, incx, y, incy, a,
                 lda)
end

function srot(n, sx, incx, sy, incy, c, s)
    return ccall((@blasfunc(srot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}), n, sx, incx, sy, incy, c, s)
end

function srotm(n, sx, incx, sy, incy, sparam)
    return ccall((@blasfunc(srotm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}), n, sx, incx, sy, incy, sparam)
end

function srotmg(sd1, sd2, sx1, sy1, sparam)
    return ccall((@blasfunc(srotmg_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ptr{Float32}),
                 sd1, sd2, sx1, sy1, sparam)
end

function ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(ssbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong),
                 uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function sscal(n, sa, sx, incx)
    return ccall((@blasfunc(sscal_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}), n, sa, sx, incx)
end

function sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    return ccall((@blasfunc(sspmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, alpha,
                 ap, x, incx, beta, y, incy, 1)
end

function sspr(uplo, n, alpha, x, incx, ap)
    return ccall((@blasfunc(sspr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Clong), uplo, n, alpha, x, incx, ap, 1)
end

function sspr2(uplo, n, alpha, x, incx, y, incy, ap)
    return ccall((@blasfunc(sspr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong), uplo, n, alpha, x, incx, y,
                 incy, ap, 1)
end

function sswap(n, sx, incx, sy, incy)
    return ccall((@blasfunc(sswap_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), n, sx, incx,
                 sy, incy)
end

function ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(ssymm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, 1,
                 1)
end

function ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(ssymv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong),
                 uplo, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function ssyr(uplo, n, alpha, x, incx, a, lda)
    return ccall((@blasfunc(ssyr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, alpha, x, incx, a, lda, 1)
end

function ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(ssyr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, alpha,
                 x, incx, y, incy, a, lda, 1)
end

function ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(ssyr2k_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1,
                 1)
end

function ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    return ccall((@blasfunc(ssyrk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), uplo,
                 trans, n, k, alpha, a, lda, beta, c, ldc, 1, 1)
end

function stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    return ccall((@blasfunc(stbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, k, a, lda, x, incx, 1, 1, 1)
end

function stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    return ccall((@blasfunc(stbsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, k, a, lda, x, incx, 1, 1, 1)
end

function stpmv(uplo, trans, diag, n, ap, x, incx)
    return ccall((@blasfunc(stpmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, ap, x, incx, 1, 1,
                 1)
end

function stpsv(uplo, trans, diag, n, ap, x, incx)
    return ccall((@blasfunc(stpsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, ap, x, incx, 1, 1,
                 1)
end

function strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    return ccall((@blasfunc(strmm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong), side, uplo, transa, diag, m, n, alpha, a, lda, b,
                 ldb, 1, 1, 1, 1)
end

function strmv(uplo, trans, diag, n, a, lda, x, incx)
    return ccall((@blasfunc(strmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, a,
                 lda, x, incx, 1, 1, 1)
end

function strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    return ccall((@blasfunc(strsm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong), side, uplo, transa, diag, m, n, alpha, a, lda, b,
                 ldb, 1, 1, 1, 1)
end

function strsv(uplo, trans, diag, n, a, lda, x, incx)
    return ccall((@blasfunc(strsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, a,
                 lda, x, incx, 1, 1, 1)
end

function zaxpy(n, za, zx, incx, zy, incy)
    return ccall((@blasfunc(zaxpy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}), n, za, zx, incx, zy, incy)
end

function zcopy(n, zx, incx, zy, incy)
    return ccall((@blasfunc(zcopy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, zx,
                 incx, zy, incy)
end

function zdotc(n, zx, incx, zy, incy)
    return ccall((@blasfunc(zdotc_), libblastrampoline), ComplexF64,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, zx,
                 incx, zy, incy)
end

function zdotu(n, zx, incx, zy, incy)
    return ccall((@blasfunc(zdotu_), libblastrampoline), ComplexF64,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, zx,
                 incx, zy, incy)
end

function zdrot(n, zx, incx, zy, incy, c, s)
    return ccall((@blasfunc(zdrot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}), n, zx, incx, zy, incy, c, s)
end

function zdscal(n, da, zx, incx)
    return ccall((@blasfunc(zdscal_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}), n, da, zx, incx)
end

function zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(zgbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), trans, m, n, kl, ku, alpha, a, lda, x,
                 incx, beta, y, incy, 1)
end

function zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(zgemm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), transa, transb, m, n, k, alpha,
                 a, lda, b, ldb, beta, c, ldc, 1, 1)
end

function zgemmtr(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(zgemmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, transa, transb, n,
                 k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1, 1)
end

function zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(zgemv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong), trans, m, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function zgerc(m, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(zgerc_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, alpha, x,
                 incx, y, incy, a, lda)
end

function zgeru(m, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(zgeru_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, alpha, x,
                 incx, y, incy, a, lda)
end

function zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(zhbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong), uplo, n, k, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(zhemm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), side, uplo, m, n, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(zhemv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Clong), uplo, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function zher(uplo, n, alpha, x, incx, a, lda)
    return ccall((@blasfunc(zher_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, alpha, x, incx, a, lda, 1)
end

function zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    return ccall((@blasfunc(zher2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 alpha, x, incx, y, incy, a, lda, 1)
end

function zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(zher2k_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), uplo, trans, n, k, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    return ccall((@blasfunc(zherk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Clong, Clong), uplo, trans, n, k, alpha, a, lda, beta, c, ldc, 1, 1)
end

function zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    return ccall((@blasfunc(zhpmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 alpha, ap, x, incx, beta, y, incy, 1)
end

function zhpr(uplo, n, alpha, x, incx, ap)
    return ccall((@blasfunc(zhpr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Clong), uplo, n, alpha, x, incx, ap, 1)
end

function zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
    return ccall((@blasfunc(zhpr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong), uplo, n, alpha, x,
                 incx, y, incy, ap, 1)
end

function zscal(n, za, zx, incx)
    return ccall((@blasfunc(zscal_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}), n, za, zx, incx)
end

function zswap(n, zx, incx, zy, incy)
    return ccall((@blasfunc(zswap_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, zx,
                 incx, zy, incy)
end

function zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(zsymm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), side, uplo, m, n, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function zsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    return ccall((@blasfunc(zsyr2k_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), uplo, trans, n, k, alpha, a,
                 lda, b, ldb, beta, c, ldc, 1, 1)
end

function zsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    return ccall((@blasfunc(zsyrk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Clong, Clong), uplo, trans, n, k, alpha, a, lda, beta, c, ldc, 1, 1)
end

function ztbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    return ccall((@blasfunc(ztbmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, k, a, lda, x, incx, 1, 1, 1)
end

function ztbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    return ccall((@blasfunc(ztbsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, k, a, lda, x, incx, 1, 1, 1)
end

function ztpmv(uplo, trans, diag, n, ap, x, incx)
    return ccall((@blasfunc(ztpmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n,
                 ap, x, incx, 1, 1, 1)
end

function ztpsv(uplo, trans, diag, n, ap, x, incx)
    return ccall((@blasfunc(ztpsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n,
                 ap, x, incx, 1, 1, 1)
end

function ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    return ccall((@blasfunc(ztrmm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), side, uplo, transa, diag, m, n, alpha, a,
                 lda, b, ldb, 1, 1, 1, 1)
end

function ztrmv(uplo, trans, diag, n, a, lda, x, incx)
    return ccall((@blasfunc(ztrmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, a,
                 lda, x, incx, 1, 1, 1)
end

function ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    return ccall((@blasfunc(ztrsm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), side, uplo, transa, diag, m, n, alpha, a,
                 lda, b, ldb, 1, 1, 1, 1)
end

function ztrsv(uplo, trans, diag, n, a, lda, x, incx)
    return ccall((@blasfunc(ztrsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, a,
                 lda, x, incx, 1, 1, 1)
end
