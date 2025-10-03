function sgedmd(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
                k, reig, imeig, z, ldz, res, b, ldb, v, ldv, work, lwork, iwork, liwork, info)
    ccall((@blasfunc(sgedmd_), libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
           Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
           Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32},
           Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
           Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong, Clong),
          jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
          k, reig, imeig, z, ldz, res, b, ldb, v, ldv, work, lwork, iwork, liwork, info,
          1, 1, 1, 1, 1, 1, 1)
end

function dgedmd(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
                k, reig, imeig, z, ldz, res, b, ldb, v, ldv, work, lwork, iwork, liwork, info)
    ccall((@blasfunc(dgedmd_), libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
           Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong, Clong),
          jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
          k, reig, imeig, z, ldz, res, b, ldb, v, ldv, work, lwork, iwork, liwork, info,
          1, 1, 1, 1, 1, 1, 1)
end

function cgedmd(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
                k, eig, z, ldz, res, b, ldb, v, ldv, work, lwork, rwork, iwork, liwork, info)
    ccall((@blasfunc(cgedmd_), libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
           Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
           Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{ComplexF32},
           Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
           Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong, Clong),
          jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
          k, eig, z, ldz, res, b, ldb, v, ldv, work, lwork, rwork, iwork, liwork, info,
          1, 1, 1, 1, 1, 1, 1)
end

function zgedmd(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
                k, eig, z, ldz, res, b, ldb, v, ldv, work, lwork, rwork, iwork, liwork, info)
    ccall((@blasfunc(zgedmd_), libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
           Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
           Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{ComplexF64},
           Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
           Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong, Clong),
          jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
          k, eig, z, ldz, res, b, ldb, v, ldv, work, lwork, rwork, iwork, liwork, info,
          1, 1, 1, 1, 1, 1, 1)
end

function sgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
                 k, reig, imeig, z, ldz, res, b, ldb, v, ldv, s, lds, work, lwork, iwork, liwork, info)
    ccall((@blasfunc(sgedmdq_), libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
           Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
           Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32},
           Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
           Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
           Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong,
           Clong),
          jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
          k, reig, imeig, z, ldz, res, b, ldb, v, ldv, s, lds, work, lwork, iwork, liwork, info,
          1, 1, 1, 1, 1, 1, 1)
end

function dgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
                 k, reig, imeig, z, ldz, res, b, ldb, v, ldv, s, lds, work, lwork, iwork, liwork, info)
    ccall((@blasfunc(dgedmdq_), libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
           Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
           Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong,
           Clong),
          jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
          k, reig, imeig, z, ldz, res, b, ldb, v, ldv, s, lds, work, lwork, iwork, liwork, info,
          1, 1, 1, 1, 1, 1, 1)
end

function cgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
                 k, eig, z, ldz, res, b, ldb, v, ldv, s, lds, work, lwork, rwork, iwork, liwork, info)
    ccall((@blasfunc(cgedmdq_), libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
           Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
           Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{ComplexF32},
           Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
           Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
           Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong, Clong),
          jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol, k, eig,
          z, ldz, res, b, ldb, v, ldv, s, lds, work, lwork, rwork, iwork, liwork, info, 1, 1, 1, 1, 1, 1, 1)
end

function zgedmdq(jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol,
                 k, eig, z, ldz, res, b, ldb, v, ldv, s, lds, work, lwork, rwork, iwork, liwork, info)
    ccall((@blasfunc(zgedmdq_), libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
           Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
           Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{ComplexF64},
           Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
           Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
           Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong, Clong),
          jobs, jobz, jobr, jobq, jobt, jobf, whtsvd, m, n, f, ldf, x, ldx, y, ldy, nrnk, tol, k, eig,
          z, ldz, res, b, ldb, v, ldv, s, lds, work, lwork, rwork, iwork, liwork, info, 1, 1, 1, 1, 1, 1, 1)
end

function slartg(f, g, cs, sn, r)
    ccall((@blasfunc(slartg_), libblastrampoline), Cvoid,
          (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}),
          f, g, cs, sn, r)
end

function dlartg(f, g, cs, sn, r)
    ccall((@blasfunc(dlartg_), libblastrampoline), Cvoid,
          (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
          f, g, cs, sn, r)
end

function clartg(f, g, cs, sn, r)
    ccall((@blasfunc(clartg_), libblastrampoline), Cvoid,
          (Ref{ComplexF32}, Ref{ComplexF32}, Ref{Float32}, Ref{ComplexF32}, Ref{ComplexF32}),
          f, g, cs, sn, r)
end

function zlartg(f, g, cs, sn, r)
    ccall((@blasfunc(zlartg_), libblastrampoline), Cvoid,
          (Ref{ComplexF64}, Ref{ComplexF64}, Ref{Float64}, Ref{ComplexF64}, Ref{ComplexF64}),
          f, g, cs, sn, r)
end

function slassq(n, x, incx, scale, sumsq)
    ccall((@blasfunc(slassq_), libblastrampoline), Cvoid,
          (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}),
          n, x, incx, scale, sumsq)
end

function dlassq(n, x, incx, scale, sumsq)
    ccall((@blasfunc(dlassq_), libblastrampoline), Cvoid,
          (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}),
          n, x, incx, scale, sumsq)
end

function classq(n, x, incx, scale, sumsq)
    ccall((@blasfunc(classq_), libblastrampoline), Cvoid,
          (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}),
          n, x, incx, scale, sumsq)
end

function zlassq(n, x, incx, scale, sumsq)
    ccall((@blasfunc(zlassq_), libblastrampoline), Cvoid,
          (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}),
          n, x, incx, scale, sumsq)
end

function ilaver(major, minor, patch)
    return ccall((@blasfunc(ilaver_), libblastrampoline), Cvoid,
                 (Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
                 major, minor, patch)
end

function cbbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2,
                ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d,
                b22e, rwork, lrwork, info)
    return ccall((@blasfunc(cbbcsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong),
                 jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2,
                 ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d,
                 b22e, rwork, lrwork, info, 1, 1, 1, 1, 1)
end

function cbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info)
    return ccall((@blasfunc(cbdsqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n,
                 ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info, 1)
end

function cgbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work,
                rwork, info)
    return ccall((@blasfunc(cgbbrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong), vect, m, n, ncc, kl, ku,
                 ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, rwork, info, 1)
end

function cgbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(cgbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ref{BlasInt}, Clong), norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work,
                 rwork, info, 1)
end

function cgbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(cgbequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{BlasInt}), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
end

function cgbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(cgbequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{BlasInt}), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
end

function cgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr,
                berr, work, rwork, info)
    return ccall((@blasfunc(cgbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong), trans, n, kl, ku, nrhs,
                 ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info,
                 1)
end

function cgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x,
                 ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info)
    return ccall((@blasfunc(cgbrfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong), trans, equed, n, kl, ku, nrhs, ab, ldab, afb,
                 ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, rwork, info, 1, 1)
end

function cgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    return ccall((@blasfunc(cgbsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), n, kl, ku, nrhs, ab,
                 ldab, ipiv, b, ldb, info)
end

function cgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
                ldb, x, ldx, rcond, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(cgbsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{UInt8}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong), fact,
                 trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x,
                 ldx, rcond, ferr, berr, work, rwork, info, 1, 1, 1)
end

function cgbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
                 ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, rwork, info)
    return ccall((@blasfunc(cgbsvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{UInt8}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong), fact,
                 trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x,
                 ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, rwork, info, 1, 1, 1)
end

function cgbtf2(m, n, kl, ku, ab, ldab, ipiv, info)
    return ccall((@blasfunc(cgbtf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, kl, ku, ab, ldab, ipiv, info)
end

function cgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
    return ccall((@blasfunc(cgbtrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, kl, ku, ab, ldab, ipiv, info)
end

function cgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    return ccall((@blasfunc(cgbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info, 1)
end

function cgebak(job, side, n, ilo, ihi, scale, m, v, ldv, info)
    return ccall((@blasfunc(cgebak_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job,
                 side, n, ilo, ihi, scale, m, v, ldv, info, 1, 1)
end

function cgebal(job, n, a, lda, ilo, ihi, scale, info)
    return ccall((@blasfunc(cgebal_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong), job, n, a, lda, ilo, ihi, scale, info, 1)
end

function cgebd2(m, n, a, lda, d, e, tauq, taup, work, info)
    return ccall((@blasfunc(cgebd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}), m, n, a, lda, d, e, tauq, taup, work, info)
end

function cgebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info)
    return ccall((@blasfunc(cgebrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, d, e, tauq, taup, work, lwork, info)
end

function cgecon(norm, n, a, lda, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(cgecon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong), norm, n,
                 a, lda, anorm, rcond, work, rwork, info, 1)
end

function cgeequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(cgeequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}), m, n,
                 a, lda, r, c, rowcnd, colcnd, amax, info)
end

function cgeequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(cgeequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}), m, n,
                 a, lda, r, c, rowcnd, colcnd, amax, info)
end

function cgees(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork,
               info)
    return ccall((@blasfunc(cgees_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), jobvs, sort,
                 select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info, 1,
                 1)
end

function cgeesx(jobvs, sort, select, sense, n, a, lda, sdim, w, vs, ldvs, rconde, rcondv,
                work, lwork, rwork, bwork, info)
    return ccall((@blasfunc(cgeesx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobvs, sort, select, sense, n,
                 a, lda, sdim, w, vs, ldvs, rconde, rcondv, work, lwork, rwork, bwork, info,
                 1, 1, 1)
end

function cgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    return ccall((@blasfunc(cgeev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), jobvl,
                 jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info, 1, 1)
end

function cgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi,
                scale, abnrm, rconde, rcondv, work, lwork, rwork, info)
    return ccall((@blasfunc(cgeevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong, Clong), balanc, jobvl,
                 jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm,
                 rconde, rcondv, work, lwork, rwork, info, 1, 1, 1, 1)
end

function cgehd2(n, ilo, ihi, a, lda, tau, work, info)
    return ccall((@blasfunc(cgehd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}), n, ilo, ihi, a, lda, tau,
                 work, info)
end

function cgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cgehrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), n, ilo, ihi, a,
                 lda, tau, work, lwork, info)
end

function cgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv,
                cwork, lwork, rwork, lrwork, iwork, info)
    return ccall((@blasfunc(cgejsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong, Clong, Clong, Clong), joba, jobu, jobv, jobr, jobt, jobp, m, n, a,
                 lda, sva, u, ldu, v, ldv, cwork, lwork, rwork, lrwork, iwork, info, 1, 1,
                 1, 1, 1, 1)
end

function cgelq(m, n, a, lda, t, tsize, work, lwork, info)
    return ccall((@blasfunc(cgelq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, tsize,
                 work, lwork, info)
end

function cgelq2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(cgelq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function cgelqf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cgelqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function cgelqt(m, n, mb, a, lda, t, ldt, work, info)
    return ccall((@blasfunc(cgelqt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, mb, a, lda,
                 t, ldt, work, info)
end

function cgelqt3(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(cgelqt3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function cgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(cgels_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function cgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(cgelsd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), m, n,
                 nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info)
end

function cgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info)
    return ccall((@blasfunc(cgelss_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), m, n, nrhs, a, lda,
                 b, ldb, s, rcond, rank, work, lwork, rwork, info)
end

function cgelst(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(cgelst_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function cgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info)
    return ccall((@blasfunc(cgelsy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), m, n, nrhs, a, lda,
                 b, ldb, jpvt, rcond, rank, work, lwork, rwork, info)
end

function cgemlq(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cgemlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, a, lda, t, tsize, c, ldc, work, lwork, info, 1, 1)
end

function cgemlqt(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info)
    return ccall((@blasfunc(cgemlqt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, mb, v, ldv, t, ldt, c, ldc, work, info, 1, 1)
end

function cgemqr(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cgemqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, a, lda, t, tsize, c, ldc, work, lwork, info, 1, 1)
end

function cgemqrt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info)
    return ccall((@blasfunc(cgemqrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, nb, v, ldv, t, ldt, c, ldc, work, info, 1, 1)
end

function cgeql2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(cgeql2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function cgeqlf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cgeqlf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function cgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
    return ccall((@blasfunc(cgeqp3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), m,
                 n, a, lda, jpvt, tau, work, lwork, rwork, info)
end

function cgeqp3rk(m, n, nrhs, kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk,
                  jpiv, tau, work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(cgeqp3rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs, kmax, abstol, reltol, a, lda, k,
                 maxc2nrmk, relmaxc2nrmk, jpiv, tau, work, lwork, rwork, iwork, info)
end

function cgeqr(m, n, a, lda, t, tsize, work, lwork, info)
    return ccall((@blasfunc(cgeqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, tsize,
                 work, lwork, info)
end

function cgeqr2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(cgeqr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function cgeqr2p(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(cgeqr2p_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function cgeqrf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cgeqrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function cgeqrfp(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cgeqrfp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function cgeqrt(m, n, nb, a, lda, t, ldt, work, info)
    return ccall((@blasfunc(cgeqrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, nb, a, lda,
                 t, ldt, work, info)
end

function cgeqrt2(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(cgeqrt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function cgeqrt3(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(cgeqrt3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function cgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(cgerfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{Float32}, Ref{BlasInt}, Clong), trans, n, nrhs, a, lda, af, ldaf, ipiv,
                 b, ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function cgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 rwork, info)
    return ccall((@blasfunc(cgerfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong),
                 trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 rwork, info, 1, 1)
end

function cgerq2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(cgerq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function cgerqf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cgerqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function cgesc2(n, a, lda, rhs, ipiv, jpiv, scale)
    return ccall((@blasfunc(cgesc2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{Float32}), n, a, lda, rhs, ipiv, jpiv, scale)
end

function cgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(cgesdd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info,
                 1)
end

function cgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(cgesv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), n, nrhs, a, lda, ipiv, b, ldb,
                 info)
end

function cgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
    return ccall((@blasfunc(cgesvd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), jobu,
                 jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info, 1, 1)
end

function cgesvdq(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank,
                 iwork, liwork, cwork, lcwork, rwork, lrwork, info)
    return ccall((@blasfunc(cgesvdq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong, Clong), joba, jobp, jobr, jobu, jobv, m, n, a, lda,
                 s, u, ldu, v, ldv, numrank, iwork, liwork, cwork, lcwork, rwork, lrwork,
                 info, 1, 1, 1, 1, 1)
end

function cgesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt,
                 work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(cgesvdx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u,
                 ldu, vt, ldvt, work, lwork, rwork, iwork, info, 1, 1, 1)
end

function cgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, cwork, lwork, rwork,
                lrwork, info)
    return ccall((@blasfunc(cgesvj_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, cwork,
                 lwork, rwork, lrwork, info, 1, 1, 1)
end

function cgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                rcond, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(cgesvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, nrhs, a, lda, af, ldaf,
                 ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1,
                 1, 1)
end

function cgesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info)
    return ccall((@blasfunc(cgesvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, nrhs, a, lda, af, ldaf,
                 ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds,
                 err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, 1, 1, 1)
end

function cgetc2(n, a, lda, ipiv, jpiv, info)
    return ccall((@blasfunc(cgetc2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                 n, a, lda, ipiv, jpiv, info)
end

function cgetf2(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(cgetf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                 m, n, a, lda, ipiv, info)
end

function cgetrf(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(cgetrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                 m, n, a, lda, ipiv, info)
end

function cgetrf2(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(cgetrf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                 m, n, a, lda, ipiv, info)
end

function cgetri(n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(cgetri_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}), n, a, lda, ipiv, work, lwork, info)
end

function cgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(cgetrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), trans, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function cgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(cgetsls_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function cgetsqrhrt(m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(cgetsqrhrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info)
end

function cggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info)
    return ccall((@blasfunc(cggbak_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info, 1, 1)
end

function cggbal(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info)
    return ccall((@blasfunc(cggbal_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Clong), job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work,
                 info, 1)
end

function cgges(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl,
               ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info)
    return ccall((@blasfunc(cgges_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim,
                 alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info, 1, 1,
                 1)
end

function cgges3(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl,
                ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info)
    return ccall((@blasfunc(cgges3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim,
                 alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info, 1, 1,
                 1)
end

function cggesx(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta,
                vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork,
                bwork, info)
    return ccall((@blasfunc(cggesx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), jobvsl, jobvsr, sort, selctg, sense, n, a,
                 lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv,
                 work, lwork, rwork, iwork, liwork, bwork, info, 1, 1, 1, 1)
end

function cggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work,
               lwork, rwork, info)
    return ccall((@blasfunc(cggev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), jobvl, jobvr, n, a,
                 lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info, 1,
                 1)
end

function cggev3(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work,
                lwork, rwork, info)
    return ccall((@blasfunc(cggev3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), jobvl, jobvr, n, a,
                 lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info, 1,
                 1)
end

function cggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr,
                ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork,
                rwork, iwork, bwork, info)
    return ccall((@blasfunc(cggevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl,
                 ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv,
                 work, lwork, rwork, iwork, bwork, info, 1, 1, 1, 1)
end

function cggglm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info)
    return ccall((@blasfunc(cggglm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), n, m, p, a, lda,
                 b, ldb, d, x, y, work, lwork, info)
end

function cgghd3(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork,
                info)
    return ccall((@blasfunc(cgghd3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work,
                 lwork, info, 1, 1)
end

function cgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info)
    return ccall((@blasfunc(cgghrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), compq, compz, n,
                 ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info, 1, 1)
end

function cgglse(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    return ccall((@blasfunc(cgglse_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, p, a, lda,
                 b, ldb, c, d, x, work, lwork, info)
end

function cggqrf(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info)
    return ccall((@blasfunc(cggqrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), n, m, p, a, lda, taua, b, ldb,
                 taub, work, lwork, info)
end

function cggrqf(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info)
    return ccall((@blasfunc(cggrqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, p, n, a, lda, taua, b, ldb,
                 taub, work, lwork, info)
end

function cggsvd3(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v,
                 ldv, q, ldq, work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(cggsvd3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobu,
                 jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q,
                 ldq, work, lwork, rwork, iwork, info, 1, 1, 1)
end

function cggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v,
                 ldv, q, ldq, iwork, rwork, tau, work, lwork, info)
    return ccall((@blasfunc(cggsvp3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola,
                 tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info,
                 1, 1, 1)
end

function cgsvj0(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work,
                lwork, info)
    return ccall((@blasfunc(cgsvj0_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), jobv, m, n, a, lda, d, sva, mv, v, ldv, eps,
                 sfmin, tol, nsweep, work, lwork, info, 1)
end

function cgsvj1(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work,
                lwork, info)
    return ccall((@blasfunc(cgsvj1_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), jobv, m, n, n1, a, lda, d, sva, mv, v, ldv,
                 eps, sfmin, tol, nsweep, work, lwork, info, 1)
end

function cgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(cgtcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong), norm, n, dl, d, du, du2, ipiv, anorm, rcond, work,
                 info, 1)
end

function cgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr,
                berr, work, rwork, info)
    return ccall((@blasfunc(cgtrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ref{BlasInt}, Clong), trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b,
                 ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function cgtsv(n, nrhs, dl, d, du, b, ldb, info)
    return ccall((@blasfunc(cgtsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), n, nrhs, dl, d, du, b, ldb, info)
end

function cgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx,
                rcond, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(cgtsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), fact, trans, n,
                 nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr,
                 berr, work, rwork, info, 1, 1)
end

function cgttrf(n, dl, d, du, du2, ipiv, info)
    return ccall((@blasfunc(cgttrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ref{BlasInt}), n, dl, d, du, du2, ipiv, info)
end

function cgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info)
    return ccall((@blasfunc(cgttrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info, 1)
end

function cgtts2(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb)
    return ccall((@blasfunc(cgtts2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}),
                 itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb)
end

function chb2st_kernels(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt,
                        work)
    return ccall((@blasfunc(chb2st_kernels_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong),
                 uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work,
                 1)
end

function chbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, rwork, info)
    return ccall((@blasfunc(chbev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z, ldz, work,
                 rwork, info, 1, 1)
end

function chbev_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, info)
    return ccall((@blasfunc(chbev_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z,
                 ldz, work, lwork, rwork, info, 1, 1)
end

function chbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork,
                liwork, info)
    return ccall((@blasfunc(chbevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork,
                 liwork, info, 1, 1)
end

function chbevd_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork,
                       iwork, liwork, info)
    return ccall((@blasfunc(chbevd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork,
                 liwork, info, 1, 1)
end

function chbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z,
                ldz, work, rwork, iwork, ifail, info)
    return ccall((@blasfunc(chbevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo, n, kd, ab,
                 ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork,
                 ifail, info, 1, 1, 1)
end

function chbevx_2stage(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol,
                       m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
    return ccall((@blasfunc(chbevx_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo,
                 n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                 rwork, iwork, ifail, info, 1, 1, 1)
end

function chbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, rwork, info)
    return ccall((@blasfunc(chbgst_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), vect, uplo, n,
                 ka, kb, ab, ldab, bb, ldbb, x, ldx, work, rwork, info, 1, 1)
end

function chbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, rwork, info)
    return ccall((@blasfunc(chbgv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), jobz,
                 uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, rwork, info, 1, 1)
end

function chbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, rwork,
                lrwork, iwork, liwork, info)
    return ccall((@blasfunc(chbgvd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, ka, kb, ab, ldab, bb,
                 ldbb, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function chbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu,
                abstol, m, w, z, ldz, work, rwork, iwork, ifail, info)
    return ccall((@blasfunc(chbgvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz,
                 range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol,
                 m, w, z, ldz, work, rwork, iwork, ifail, info, 1, 1, 1)
end

function chbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info)
    return ccall((@blasfunc(chbtrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong, Clong), vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work,
                 info, 1, 1)
end

function checon(uplo, n, a, lda, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(checon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, ipiv, anorm, rcond, work, info, 1)
end

function checon_3(uplo, n, a, lda, e, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(checon_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, e, ipiv, anorm, rcond, work, info, 1)
end

function checon_rook(uplo, n, a, lda, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(checon_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, ipiv, anorm, rcond, work, info, 1)
end

function cheequb(uplo, n, a, lda, s, scond, amax, work, info)
    return ccall((@blasfunc(cheequb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, s, scond, amax, work, info, 1)
end

function cheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
    return ccall((@blasfunc(cheev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong), jobz, uplo, n, a, lda, w, work, lwork, rwork, info, 1, 1)
end

function cheev_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
    return ccall((@blasfunc(cheev_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong), jobz, uplo, n, a, lda, w, work, lwork, rwork, info, 1, 1)
end

function cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    return ccall((@blasfunc(cheevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, a, lda, w,
                 work, lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function cheevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork,
                       info)
    return ccall((@blasfunc(cheevd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, a, lda, w,
                 work, lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function cheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz,
                work, lwork, rwork, lrwork, iwork, liwork, info)
    return ccall((@blasfunc(cheevr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), jobz, range, uplo, n, a, lda, vl, vu, il, iu,
                 abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork,
                 info, 1, 1, 1)
end

function cheevr_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                       isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    return ccall((@blasfunc(cheevr_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), jobz, range, uplo, n, a, lda, vl, vu, il, iu,
                 abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork,
                 info, 1, 1, 1)
end

function cheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work,
                lwork, rwork, iwork, ifail, info)
    return ccall((@blasfunc(cheevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz,
                 range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                 rwork, iwork, ifail, info, 1, 1, 1)
end

function cheevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                       work, lwork, rwork, iwork, ifail, info)
    return ccall((@blasfunc(cheevx_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz,
                 range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                 rwork, iwork, ifail, info, 1, 1, 1)
end

function chegs2(itype, uplo, n, a, lda, b, ldb, info)
    return ccall((@blasfunc(chegs2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), itype, uplo, n, a, lda, b,
                 ldb, info, 1)
end

function chegst(itype, uplo, n, a, lda, b, ldb, info)
    return ccall((@blasfunc(chegst_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), itype, uplo, n, a, lda, b,
                 ldb, info, 1)
end

function chegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    return ccall((@blasfunc(chegv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, a, lda, b,
                 ldb, w, work, lwork, rwork, info, 1, 1)
end

function chegv_2stage(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    return ccall((@blasfunc(chegv_2stage_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, a, lda, b,
                 ldb, w, work, lwork, rwork, info, 1, 1)
end

function chegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork,
                liwork, info)
    return ccall((@blasfunc(chegvd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork,
                 liwork, info, 1, 1)
end

function chegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w,
                z, ldz, work, lwork, rwork, iwork, ifail, info)
    return ccall((@blasfunc(chegvx_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), itype, jobz, range,
                 uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                 rwork, iwork, ifail, info, 1, 1, 1)
end

function cherfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(cherfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, af, ldaf, ipiv, b,
                 ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function cherfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 rwork, info)
    return ccall((@blasfunc(cherfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), uplo, equed, n,
                 nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds,
                 err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, 1, 1)
end

function chesv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(chesv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function chesv_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(chesv_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function chesv_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork,
                         info)
    return ccall((@blasfunc(chesv_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info, 1)
end

function chesv_rk(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(chesv_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, e, ipiv, b, ldb,
                 work, lwork, info, 1)
end

function chesv_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(chesv_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function chesvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr,
                berr, work, lwork, rwork, info)
    return ccall((@blasfunc(chesvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), fact,
                 uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr,
                 work, lwork, rwork, info, 1, 1)
end

function chesvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info)
    return ccall((@blasfunc(chesvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong, Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b,
                 ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, rwork, info, 1, 1, 1)
end

function cheswapr(uplo, n, a, lda, i1, i2)
    return ccall((@blasfunc(cheswapr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, i1, i2, 1)
end

function chetd2(uplo, n, a, lda, d, e, tau, info)
    return ccall((@blasfunc(chetd2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, a, lda, d, e,
                 tau, info, 1)
end

function chetf2(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(chetf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function chetf2_rk(uplo, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(chetf2_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv, info, 1)
end

function chetf2_rook(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(chetf2_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function chetrd(uplo, n, a, lda, d, e, tau, work, lwork, info)
    return ccall((@blasfunc(chetrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, d, e, tau, work, lwork, info, 1)
end

function chetrd_2stage(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info)
    return ccall((@blasfunc(chetrd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), vect, uplo, n, a,
                 lda, d, e, tau, hous2, lhous2, work, lwork, info, 1, 1)
end

function chetrd_hb2st(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork,
                      info)
    return ccall((@blasfunc(chetrd_hb2st_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), stage1, vect,
                 uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info, 1, 1, 1)
end

function chetrd_he2hb(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info)
    return ccall((@blasfunc(chetrd_he2hb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info,
                 1)
end

function chetrf(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(chetrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function chetrf_aa(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(chetrf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function chetrf_aa_2stage(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info)
    return ccall((@blasfunc(chetrf_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info, 1)
end

function chetrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(chetrf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, lwork, info, 1)
end

function chetrf_rook(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(chetrf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function chetri(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(chetri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function chetri2(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(chetri2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function chetri2x(uplo, n, a, lda, ipiv, work, nb, info)
    return ccall((@blasfunc(chetri2x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, nb, info, 1)
end

function chetri_3(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(chetri_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, lwork, info, 1)
end

function chetri_3x(uplo, n, a, lda, e, ipiv, work, nb, info)
    return ccall((@blasfunc(chetri_3x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, nb, info, 1)
end

function chetri_rook(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(chetri_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function chetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(chetrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function chetrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info)
    return ccall((@blasfunc(chetrs2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 nrhs, a, lda, ipiv, b, ldb, work, info, 1)
end

function chetrs_3(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info)
    return ccall((@blasfunc(chetrs_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info, 1)
end

function chetrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(chetrs_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function chetrs_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info)
    return ccall((@blasfunc(chetrs_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2,
                 b, ldb, info, 1)
end

function chetrs_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(chetrs_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function chfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c)
    return ccall((@blasfunc(chfrk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Clong, Clong,
                  Clong), transr, uplo, trans, n, k, alpha, a, lda, beta, c, 1, 1, 1)
end

function chgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz,
                work, lwork, rwork, info)
    return ccall((@blasfunc(chgeqz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong),
                 job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z,
                 ldz, work, lwork, rwork, info, 1, 1, 1)
end

function chpcon(uplo, n, ap, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(chpcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, ap, ipiv,
                 anorm, rcond, work, info, 1)
end

function chpev(jobz, uplo, n, ap, w, z, ldz, work, rwork, info)
    return ccall((@blasfunc(chpev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), jobz, uplo, n, ap, w, z, ldz, work, rwork, info, 1, 1)
end

function chpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork,
                info)
    return ccall((@blasfunc(chpevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n,
                 ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function chpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork,
                iwork, ifail, info)
    return ccall((@blasfunc(chpevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo,
                 n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail,
                 info, 1, 1, 1)
end

function chpgst(itype, uplo, n, ap, bp, info)
    return ccall((@blasfunc(chpgst_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong), itype, uplo, n, ap, bp, info, 1)
end

function chpgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, info)
    return ccall((@blasfunc(chpgv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), itype, jobz,
                 uplo, n, ap, bp, w, z, ldz, work, rwork, info, 1, 1)
end

function chpgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, rwork, lrwork, iwork,
                liwork, info)
    return ccall((@blasfunc(chpgvd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, ap, bp, w, z, ldz, work,
                 lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function chpgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz,
                work, rwork, iwork, ifail, info)
    return ccall((@blasfunc(chpgvx_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu,
                 abstol, m, w, z, ldz, work, rwork, iwork, ifail, info, 1, 1, 1)
end

function chprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(chprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt},
                  Clong), uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 rwork, info, 1)
end

function chpsv(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(chpsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b,
                 ldb, info, 1)
end

function chpsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(chpsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), fact, uplo, n, nrhs, ap, afp,
                 ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1, 1)
end

function chptrd(uplo, n, ap, d, e, tau, info)
    return ccall((@blasfunc(chptrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, ap, d, e, tau, info, 1)
end

function chptrf(uplo, n, ap, ipiv, info)
    return ccall((@blasfunc(chptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, ap, ipiv, info, 1)
end

function chptri(uplo, n, ap, ipiv, work, info)
    return ccall((@blasfunc(chptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong), uplo, n, ap, ipiv, work, info, 1)
end

function chptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(chptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b,
                 ldb, info, 1)
end

function chsein(side, eigsrc, initv, select, n, h, ldh, w, vl, ldvl, vr, ldvr, mm, m, work,
                rwork, ifaill, ifailr, info)
    return ccall((@blasfunc(chsein_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), side, eigsrc, initv, select,
                 n, h, ldh, w, vl, ldvl, vr, ldvr, mm, m, work, rwork, ifaill, ifailr, info,
                 1, 1, 1)
end

function chseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info)
    return ccall((@blasfunc(chseqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job, compz, n, ilo, ihi, h, ldh, w,
                 z, ldz, work, lwork, info, 1, 1)
end

function cla_gbamv(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy)
    return ccall((@blasfunc(cla_gbamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}), trans, m, n, kl, ku, alpha, ab, ldab, x, incx,
                 beta, y, incy)
end

function cla_gbrcond_c(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, capply, info, work,
                       rwork)
    return ccall((@blasfunc(cla_gbrcond_c_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Clong), trans, n, kl, ku, ab, ldab, afb,
                 ldafb, ipiv, c, capply, info, work, rwork, 1)
end

function cla_gbrcond_x(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, x, info, work, rwork)
    return ccall((@blasfunc(cla_gbrcond_x_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Clong), trans, n, kl, ku, ab, ldab, afb,
                 ldafb, ipiv, x, info, work, rwork, 1)
end

function cla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb,
                             ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms,
                             err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond,
                             ithresh, rthresh, dz_ub, ignore_cwise, info)
    return ccall((@blasfunc(cla_gbrfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb,
                 ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm,
                 err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                 ignore_cwise, info)
end

function cla_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb)
    return ccall((@blasfunc(cla_gbrpvgrw_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}), n, kl, ku, ncols, ab, ldab, afb, ldafb)
end

function cla_geamv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(cla_geamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}), trans,
                 m, n, alpha, a, lda, x, incx, beta, y, incy)
end

function cla_gercond_c(trans, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork)
    return ccall((@blasfunc(cla_gercond_c_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Clong), trans, n, a, lda, af, ldaf, ipiv, c, capply, info,
                 work, rwork, 1)
end

function cla_gercond_x(trans, n, a, lda, af, ldaf, ipiv, x, info, work, rwork)
    return ccall((@blasfunc(cla_gercond_x_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Clong), trans, n, a, lda, af, ldaf, ipiv, x, info, work,
                 rwork, 1)
end

function cla_gerfsx_extended(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ,
                             c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb,
                             dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info)
    return ccall((@blasfunc(cla_gerfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}), prec_type, trans_type, n, nrhs, a,
                 lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n,
                 errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise,
                 info)
end

function cla_gerpvgrw(n, ncols, a, lda, af, ldaf)
    return ccall((@blasfunc(cla_gerpvgrw_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}), n, ncols, a, lda, af, ldaf)
end

function cla_heamv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(cla_heamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}), uplo,
                 n, alpha, a, lda, x, incx, beta, y, incy)
end

function cla_hercond_c(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork)
    return ccall((@blasfunc(cla_hercond_c_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Clong), uplo, n, a, lda, af, ldaf, ipiv, c, capply, info,
                 work, rwork, 1)
end

function cla_hercond_x(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork)
    return ccall((@blasfunc(cla_hercond_x_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Clong), uplo, n, a, lda, af, ldaf, ipiv, x, info, work,
                 rwork, 1)
end

function cla_herfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                             ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp,
                             res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                             ignore_cwise, info)
    return ccall((@blasfunc(cla_herfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), prec_type, uplo, n, nrhs, a,
                 lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms,
                 err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh,
                 rthresh, dz_ub, ignore_cwise, info, 1)
end

function cla_herpvgrw(uplo, n, info, a, lda, af, ldaf, ipiv, work)
    return ccall((@blasfunc(cla_herpvgrw_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Clong), uplo, n,
                 info, a, lda, af, ldaf, ipiv, work, 1)
end

function cla_lin_berr(n, nz, nrhs, res, ayb, berr)
    return ccall((@blasfunc(cla_lin_berr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32},
                  Ptr{Float32}), n, nz, nrhs, res, ayb, berr)
end

function cla_porcond_c(uplo, n, a, lda, af, ldaf, c, capply, info, work, rwork)
    return ccall((@blasfunc(cla_porcond_c_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Clong), uplo, n, a, lda, af, ldaf, c, capply, info, work,
                 rwork, 1)
end

function cla_porcond_x(uplo, n, a, lda, af, ldaf, x, info, work, rwork)
    return ccall((@blasfunc(cla_porcond_x_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32},
                  Clong), uplo, n, a, lda, af, ldaf, x, info, work, rwork, 1)
end

function cla_porfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb,
                             y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res,
                             ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise,
                             info)
    return ccall((@blasfunc(cla_porfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), prec_type, uplo, n, nrhs, a,
                 lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm,
                 err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                 ignore_cwise, info, 1)
end

function cla_porpvgrw(uplo, ncols, a, lda, af, ldaf, work)
    return ccall((@blasfunc(cla_porpvgrw_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Clong), uplo, ncols, a, lda, af, ldaf, work, 1)
end

function cla_syamv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(cla_syamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}), uplo,
                 n, alpha, a, lda, x, incx, beta, y, incy)
end

function cla_syrcond_c(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork)
    return ccall((@blasfunc(cla_syrcond_c_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Clong), uplo, n, a, lda, af, ldaf, ipiv, c, capply, info,
                 work, rwork, 1)
end

function cla_syrcond_x(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork)
    return ccall((@blasfunc(cla_syrcond_x_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Clong), uplo, n, a, lda, af, ldaf, ipiv, x, info, work,
                 rwork, 1)
end

function cla_syrfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                             ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp,
                             res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                             ignore_cwise, info)
    return ccall((@blasfunc(cla_syrfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), prec_type, uplo, n, nrhs, a,
                 lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms,
                 err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh,
                 rthresh, dz_ub, ignore_cwise, info, 1)
end

function cla_syrpvgrw(uplo, n, info, a, lda, af, ldaf, ipiv, work)
    return ccall((@blasfunc(cla_syrpvgrw_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Clong), uplo, n,
                 info, a, lda, af, ldaf, ipiv, work, 1)
end

function cla_wwaddw(n, x, y, w)
    return ccall((@blasfunc(cla_wwaddw_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}), n, x, y, w)
end

function clabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy)
    return ccall((@blasfunc(clabrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, nb, a, lda, d, e, tauq,
                 taup, x, ldx, y, ldy)
end

function clacgv(n, x, incx)
    return ccall((@blasfunc(clacgv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, x, incx)
end

function clacn2(n, v, x, est, kase, isave)
    return ccall((@blasfunc(clacn2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}), n, v, x, est, kase, isave)
end

function clacon(n, v, x, est, kase)
    return ccall((@blasfunc(clacon_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ref{BlasInt}), n,
                 v, x, est, kase)
end

function clacp2(uplo, m, n, a, lda, b, ldb)
    return ccall((@blasfunc(clacp2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, m, n, a, lda, b, ldb, 1)
end

function clacpy(uplo, m, n, a, lda, b, ldb)
    return ccall((@blasfunc(clacpy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, m, n, a, lda, b, ldb, 1)
end

function clacrm(m, n, a, lda, b, ldb, c, ldc, rwork)
    return ccall((@blasfunc(clacrm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}), m, n, a, lda, b, ldb, c, ldc,
                 rwork)
end

function clacrt(n, cx, incx, cy, incy, c, s)
    return ccall((@blasfunc(clacrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{ComplexF32}, Ref{ComplexF32}), n, cx, incx, cy, incy, c, s)
end

function cladiv(x, y)
    return ccall((@blasfunc(cladiv_), libblastrampoline), ComplexF32,
                 (Ref{ComplexF32}, Ref{ComplexF32}), x, y)
end

function claed0(qsiz, n, d, e, q, ldq, qstore, ldqs, rwork, iwork, info)
    return ccall((@blasfunc(claed0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ref{BlasInt}), qsiz, n, d, e, q, ldq, qstore, ldqs, rwork, iwork, info)
end

function claed7(n, cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, ldq, rho, indxq, qstore, qptr,
                prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, info)
    return ccall((@blasfunc(claed7_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), n,
                 cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, ldq, rho, indxq, qstore, qptr,
                 prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, info)
end

function claed8(k, n, qsiz, q, ldq, d, rho, cutpnt, z, dlambda, q2, ldq2, w, indxp, indx,
                indxq, perm, givptr, givcol, givnum, info)
    return ccall((@blasfunc(claed8_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}), k, n, qsiz, q, ldq, d,
                 rho, cutpnt, z, dlambda, q2, ldq2, w, indxp, indx, indxq, perm, givptr,
                 givcol, givnum, info)
end

function claein(rightv, noinit, n, h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info)
    return ccall((@blasfunc(claein_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}), rightv, noinit, n,
                 h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info)
end

function claesy(a, b, c, rt1, rt2, evscal, cs1, sn1)
    return ccall((@blasfunc(claesy_), libblastrampoline), Cvoid,
                 (Ref{ComplexF32}, Ref{ComplexF32}, Ref{ComplexF32}, Ref{ComplexF32},
                  Ref{ComplexF32}, Ref{ComplexF32}, Ref{ComplexF32}, Ref{ComplexF32}), a, b,
                 c, rt1, rt2, evscal, cs1, sn1)
end

function claev2(a, b, c, rt1, rt2, cs1, sn1)
    return ccall((@blasfunc(claev2_), libblastrampoline), Cvoid,
                 (Ref{ComplexF32}, Ref{ComplexF32}, Ref{ComplexF32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{ComplexF32}), a, b, c, rt1, rt2, cs1, sn1)
end

function clag2z(m, n, sa, ldsa, a, lda, info)
    return ccall((@blasfunc(clag2z_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, sa, ldsa, a, lda, info)
end

function clags2(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq)
    return ccall((@blasfunc(clags2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ref{ComplexF32}, Ref{Float32}, Ref{Float32},
                  Ref{ComplexF32}, Ref{Float32}, Ref{Float32}, Ref{ComplexF32},
                  Ref{Float32}, Ref{ComplexF32}, Ref{Float32}, Ref{ComplexF32}), upper, a1,
                 a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq)
end

function clagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb)
    return ccall((@blasfunc(clagtm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), trans, n, nrhs, alpha,
                 dl, d, du, x, ldx, beta, b, ldb, 1)
end

function clahef(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(clahef_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb,
                 a, lda, ipiv, w, ldw, info, 1)
end

function clahef_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work)
    return ccall((@blasfunc(clahef_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong), uplo, j1,
                 m, nb, a, lda, ipiv, h, ldh, work, 1)
end

function clahef_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info)
    return ccall((@blasfunc(clahef_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info, 1)
end

function clahef_rook(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(clahef_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb,
                 a, lda, ipiv, w, ldw, info, 1)
end

function clahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info)
    return ccall((@blasfunc(clahqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz,
                 z, ldz, info)
end

function clahr2(n, k, nb, a, lda, tau, t, ldt, y, ldy)
    return ccall((@blasfunc(clahr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}),
                 n, k, nb, a, lda, tau, t, ldt, y, ldy)
end

function claic1(job, j, x, sest, w, gamma, sestpr, s, c)
    return ccall((@blasfunc(claic1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{Float32}, Ptr{ComplexF32},
                  Ref{ComplexF32}, Ref{Float32}, Ref{ComplexF32}, Ref{ComplexF32}), job, j,
                 x, sest, w, gamma, sestpr, s, c)
end

function clals0(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol,
                givnum, ldgnum, poles, difl, difr, z, k, c, s, rwork, info)
    return ccall((@blasfunc(clals0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}), icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx,
                 perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c,
                 s, rwork, info)
end

function clalsa(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z,
                poles, givptr, givcol, ldgcol, perm, givnum, c, s, rwork, iwork, info)
    return ccall((@blasfunc(clalsa_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), icompq, smlsiz, n,
                 nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr,
                 givcol, ldgcol, perm, givnum, c, s, rwork, iwork, info)
end

function clalsd(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, rwork, iwork, info)
    return ccall((@blasfunc(clalsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, smlsiz, n, nrhs, d, e,
                 b, ldb, rcond, rank, work, rwork, iwork, info, 1)
end

function clamswlq(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info)
    return ccall((@blasfunc(clamswlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork,
                 info, 1, 1)
end

function clamtsqr(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info)
    return ccall((@blasfunc(clamtsqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork,
                 info, 1, 1)
end

function clangb(norm, n, kl, ku, ab, ldab, work)
    return ccall((@blasfunc(clangb_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Clong), norm, n, kl, ku, ab, ldab, work, 1)
end

function clange(norm, m, n, a, lda, work)
    return ccall((@blasfunc(clange_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Clong), norm, m, n, a, lda, work, 1)
end

function clangt(norm, n, dl, d, du)
    return ccall((@blasfunc(clangt_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Clong), norm, n, dl, d, du, 1)
end

function clanhb(norm, uplo, n, k, ab, ldab, work)
    return ccall((@blasfunc(clanhb_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Clong, Clong), norm, uplo, n, k, ab, ldab, work, 1, 1)
end

function clanhe(norm, uplo, n, a, lda, work)
    return ccall((@blasfunc(clanhe_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Clong, Clong), norm, uplo, n, a, lda, work, 1, 1)
end

function clanhf(norm, transr, uplo, n, a, work)
    return ccall((@blasfunc(clanhf_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Clong, Clong, Clong), norm, transr, uplo, n, a, work, 1, 1,
                 1)
end

function clanhp(norm, uplo, n, ap, work)
    return ccall((@blasfunc(clanhp_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Clong,
                  Clong), norm, uplo, n, ap, work, 1, 1)
end

function clanhs(norm, n, a, lda, work)
    return ccall((@blasfunc(clanhs_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Clong),
                 norm, n, a, lda, work, 1)
end

function clanht(norm, n, d, e)
    return ccall((@blasfunc(clanht_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Clong), norm, n, d,
                 e, 1)
end

function clansb(norm, uplo, n, k, ab, ldab, work)
    return ccall((@blasfunc(clansb_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Clong, Clong), norm, uplo, n, k, ab, ldab, work, 1, 1)
end

function clansp(norm, uplo, n, ap, work)
    return ccall((@blasfunc(clansp_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Clong,
                  Clong), norm, uplo, n, ap, work, 1, 1)
end

function clansy(norm, uplo, n, a, lda, work)
    return ccall((@blasfunc(clansy_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Clong, Clong), norm, uplo, n, a, lda, work, 1, 1)
end

function clantb(norm, uplo, diag, n, k, ab, ldab, work)
    return ccall((@blasfunc(clantb_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Clong, Clong, Clong), norm, uplo, diag, n, k, ab,
                 ldab, work, 1, 1, 1)
end

function clantp(norm, uplo, diag, n, ap, work)
    return ccall((@blasfunc(clantp_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Clong, Clong, Clong), norm, uplo, diag, n, ap, work, 1, 1,
                 1)
end

function clantr(norm, uplo, diag, m, n, a, lda, work)
    return ccall((@blasfunc(clantr_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Clong, Clong, Clong), norm, uplo, diag, m, n, a,
                 lda, work, 1, 1, 1)
end

function clapll(n, x, incx, y, incy, ssmin)
    return ccall((@blasfunc(clapll_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}), n, x, incx, y, incy, ssmin)
end

function clapmr(forwrd, m, n, x, ldx, k)
    return ccall((@blasfunc(clapmr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}),
                 forwrd, m, n, x, ldx, k)
end

function clapmt(forwrd, m, n, x, ldx, k)
    return ccall((@blasfunc(clapmt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}),
                 forwrd, m, n, x, ldx, k)
end

function claqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed)
    return ccall((@blasfunc(claqgb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{UInt8}, Clong), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax,
                 equed, 1)
end

function claqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed)
    return ccall((@blasfunc(claqge_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{UInt8},
                  Clong), m, n, a, lda, r, c, rowcnd, colcnd, amax, equed, 1)
end

function claqhb(uplo, n, kd, ab, ldab, s, scond, amax, equed)
    return ccall((@blasfunc(claqhb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{UInt8}, Clong, Clong), uplo,
                 n, kd, ab, ldab, s, scond, amax, equed, 1, 1)
end

function claqhe(uplo, n, a, lda, s, scond, amax, equed)
    return ccall((@blasfunc(claqhe_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{UInt8}, Clong, Clong), uplo, n, a, lda, s,
                 scond, amax, equed, 1, 1)
end

function claqhp(uplo, n, ap, s, scond, amax, equed)
    return ccall((@blasfunc(claqhp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{UInt8}, Clong, Clong), uplo, n, ap, s, scond, amax,
                 equed, 1, 1)
end

function claqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work)
    return ccall((@blasfunc(claqp2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}), m, n,
                 offset, a, lda, jpvt, tau, vn1, vn2, work)
end

function claqp2rk(m, n, nrhs, ioffset, kmax, abstol, reltol, kp1, maxc2nrm, a, lda, k,
                  maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2, work, info)
    return ccall((@blasfunc(claqp2rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, nrhs,
                 ioffset, kmax, abstol, reltol, kp1, maxc2nrm, a, lda, k, maxc2nrmk,
                 relmaxc2nrmk, jpiv, tau, vn1, vn2, work, info)
end

function claqp3rk(m, n, nrhs, ioffset, nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb,
                  maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info)
    return ccall((@blasfunc(claqp3rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs, ioffset,
                 nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb, maxc2nrmk,
                 relmaxc2nrmk, jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info)
end

function claqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf)
    return ccall((@blasfunc(claqps_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, offset, nb, kb, a,
                 lda, jpvt, tau, vn1, vn2, auxv, f, ldf)
end

function claqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info)
    return ccall((@blasfunc(claqr0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo,
                 ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info)
end

function claqr1(n, h, ldh, s1, s2, v)
    return ccall((@blasfunc(claqr1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32}, Ref{ComplexF32},
                  Ptr{ComplexF32}), n, h, ldh, s1, s2, v)
end

function claqr2(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v,
                ldv, nh, t, ldt, nv, wv, ldwv, work, lwork)
    return ccall((@blasfunc(claqr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), wantt, wantz, n,
                 ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt,
                 nv, wv, ldwv, work, lwork)
end

function claqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v,
                ldv, nh, t, ldt, nv, wv, ldwv, work, lwork)
    return ccall((@blasfunc(claqr3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), wantt, wantz, n,
                 ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt,
                 nv, wv, ldwv, work, lwork)
end

function claqr4(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info)
    return ccall((@blasfunc(claqr4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo,
                 ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info)
end

function claqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz,
                v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh)
    return ccall((@blasfunc(claqr5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), wantt, wantz, kacc22, n, ktop,
                 kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv,
                 nh, wh, ldwh)
end

function claqsb(uplo, n, kd, ab, ldab, s, scond, amax, equed)
    return ccall((@blasfunc(claqsb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{UInt8}, Clong, Clong), uplo,
                 n, kd, ab, ldab, s, scond, amax, equed, 1, 1)
end

function claqsp(uplo, n, ap, s, scond, amax, equed)
    return ccall((@blasfunc(claqsp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{UInt8}, Clong, Clong), uplo, n, ap, s, scond, amax,
                 equed, 1, 1)
end

function claqsy(uplo, n, a, lda, s, scond, amax, equed)
    return ccall((@blasfunc(claqsy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{UInt8}, Clong, Clong), uplo, n, a, lda, s,
                 scond, amax, equed, 1, 1)
end

function clar1v(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz,
                mingma, r, isuppz, nrminv, resid, rqcorr, work)
    return ccall((@blasfunc(clar1v_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ptr{Float32}), n, b1, bn,
                 lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r,
                 isuppz, nrminv, resid, rqcorr, work)
end

function clar2v(n, x, y, z, incx, c, s, incc)
    return ccall((@blasfunc(clar2v_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}), n, x, y, z, incx, c, s, incc)
end

function clarcm(m, n, a, lda, b, ldb, c, ldc, rwork)
    return ccall((@blasfunc(clarcm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}), m, n, a, lda, b, ldb, c, ldc,
                 rwork)
end

function clarf(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(clarf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong),
                 side, m, n, v, incv, tau, c, ldc, work, 1)
end

function clarf1f(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(clarf1f_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong),
                 side, m, n, v, incv, tau, c, ldc, work, 1)
end

function clarf1l(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(clarf1l_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong),
                 side, m, n, v, incv, tau, c, ldc, work, 1)
end

function clarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork)
    return ccall((@blasfunc(clarfb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong,
                  Clong, Clong), side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c,
                 ldc, work, ldwork, 1, 1, 1, 1)
end

function clarfb_gett(ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork)
    return ccall((@blasfunc(clarfb_gett_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong), ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork,
                 1)
end

function clarfg(n, alpha, x, incx, tau)
    return ccall((@blasfunc(clarfg_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32}),
                 n, alpha, x, incx, tau)
end

function clarfgp(n, alpha, x, incx, tau)
    return ccall((@blasfunc(clarfgp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32}),
                 n, alpha, x, incx, tau)
end

function clarft(direct, storev, n, k, v, ldv, tau, t, ldt)
    return ccall((@blasfunc(clarft_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), direct,
                 storev, n, k, v, ldv, tau, t, ldt, 1, 1)
end

function clarfx(side, m, n, v, tau, c, ldc, work)
    return ccall((@blasfunc(clarfx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong), side, m, n, v, tau,
                 c, ldc, work, 1)
end

function clarfy(uplo, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(clarfy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong), uplo, n, v, incv,
                 tau, c, ldc, work, 1)
end

function clargv(n, x, incx, y, incy, c, incc)
    return ccall((@blasfunc(clargv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}), n, x, incx, y, incy, c, incc)
end

function clarnv(idist, iseed, n, x)
    return ccall((@blasfunc(clarnv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}), idist, iseed, n, x)
end

function clarrv(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr,
                wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info)
    return ccall((@blasfunc(clarrv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), n, vl, vu, d, l, pivmin, isplit, m,
                 dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z,
                 ldz, isuppz, work, iwork, info)
end

function clarscl2(m, n, d, x, ldx)
    return ccall((@blasfunc(clarscl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, d,
                 x, ldx)
end

function clartv(n, x, incx, y, incy, c, s, incc)
    return ccall((@blasfunc(clartv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}), n, x, incx, y, incy, c, s,
                 incc)
end

function clarz(side, m, n, l, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(clarz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong),
                 side, m, n, l, v, incv, tau, c, ldc, work, 1)
end

function clarzb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work,
                ldwork)
    return ccall((@blasfunc(clarzb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong), side, trans, direct, storev, m, n, k, l, v, ldv, t,
                 ldt, c, ldc, work, ldwork, 1, 1, 1, 1)
end

function clarzt(direct, storev, n, k, v, ldv, tau, t, ldt)
    return ccall((@blasfunc(clarzt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), direct,
                 storev, n, k, v, ldv, tau, t, ldt, 1, 1)
end

function clascl(type, kl, ku, cfrom, cto, m, n, a, lda, info)
    return ccall((@blasfunc(clascl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), type, kl, ku,
                 cfrom, cto, m, n, a, lda, info, 1)
end

function clascl2(m, n, d, x, ldx)
    return ccall((@blasfunc(clascl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, d,
                 x, ldx)
end

function claset(uplo, m, n, alpha, beta, a, lda)
    return ccall((@blasfunc(claset_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF32}, Ref{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, m, n, alpha, beta, a, lda, 1)
end

function clasr(side, pivot, direct, m, n, c, s, a, lda)
    return ccall((@blasfunc(clasr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong), side,
                 pivot, direct, m, n, c, s, a, lda, 1, 1, 1)
end

function claswlq(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(claswlq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 mb, nb, a, lda, t, ldt, work, lwork, info)
end

function claswp(n, a, lda, k1, k2, ipiv, incx)
    return ccall((@blasfunc(claswp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, a, lda, k1, k2, ipiv, incx)
end

function clasyf(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(clasyf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb,
                 a, lda, ipiv, w, ldw, info, 1)
end

function clasyf_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work)
    return ccall((@blasfunc(clasyf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Clong), uplo, j1,
                 m, nb, a, lda, ipiv, h, ldh, work, 1)
end

function clasyf_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info)
    return ccall((@blasfunc(clasyf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info, 1)
end

function clasyf_rook(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(clasyf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb,
                 a, lda, ipiv, w, ldw, info, 1)
end

function clatbs(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info)
    return ccall((@blasfunc(clatbs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong), uplo, trans, diag, normin, n, kd,
                 ab, ldab, x, scale, cnorm, info, 1, 1, 1, 1)
end

function clatdf(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv)
    return ccall((@blasfunc(clatdf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{Float32}, Ref{Float32}, Ptr{BlasInt}, Ptr{BlasInt}), ijob, n, z, ldz, rhs,
                 rdsum, rdscal, ipiv, jpiv)
end

function clatps(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info)
    return ccall((@blasfunc(clatps_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), uplo, trans, diag, normin, n, ap, x, scale,
                 cnorm, info, 1, 1, 1, 1)
end

function clatrd(uplo, n, nb, a, lda, e, tau, w, ldw)
    return ccall((@blasfunc(clatrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo,
                 n, nb, a, lda, e, tau, w, ldw, 1)
end

function clatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info)
    return ccall((@blasfunc(clatrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong), uplo, trans, diag, normin, n, a,
                 lda, x, scale, cnorm, info, 1, 1, 1, 1)
end

function clatrs3(uplo, trans, diag, normin, n, nrhs, a, lda, x, ldx, scale, cnorm, work,
                 lwork, info)
    return ccall((@blasfunc(clatrs3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), uplo, trans, diag, normin, n, nrhs, a, lda, x, ldx, scale, cnorm,
                 work, lwork, info, 1, 1, 1, 1)
end

function clatrz(m, n, l, a, lda, tau, work)
    return ccall((@blasfunc(clatrz_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}), m, n, l, a, lda, tau, work)
end

function clatsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(clatsqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 mb, nb, a, lda, t, ldt, work, lwork, info)
end

function claunhr_col_getrfnp(m, n, a, lda, d, info)
    return ccall((@blasfunc(claunhr_col_getrfnp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}), m, n, a, lda, d, info)
end

function claunhr_col_getrfnp2(m, n, a, lda, d, info)
    return ccall((@blasfunc(claunhr_col_getrfnp2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}), m, n, a, lda, d, info)
end

function clauu2(uplo, n, a, lda, info)
    return ccall((@blasfunc(clauu2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function clauum(uplo, n, a, lda, info)
    return ccall((@blasfunc(clauum_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function cpbcon(uplo, n, kd, ab, ldab, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(cpbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, anorm, rcond, work, rwork, info, 1)
end

function cpbequ(uplo, n, kd, ab, ldab, s, scond, amax, info)
    return ccall((@blasfunc(cpbequ_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Clong), uplo, n, kd,
                 ab, ldab, s, scond, amax, info, 1)
end

function cpbrfs(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(cpbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x,
                 ldx, ferr, berr, work, rwork, info, 1)
end

function cpbstf(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(cpbstf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function cpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(cpbsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab,
                 ldab, b, ldb, info, 1)
end

function cpbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx,
                rcond, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(cpbsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong, Clong), fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb,
                 equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1, 1, 1)
end

function cpbtf2(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(cpbtf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function cpbtrf(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(cpbtrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function cpbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(cpbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab,
                 ldab, b, ldb, info, 1)
end

function cpftrf(transr, uplo, n, a, info)
    return ccall((@blasfunc(cpftrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong,
                  Clong), transr, uplo, n, a, info, 1, 1)
end

function cpftri(transr, uplo, n, a, info)
    return ccall((@blasfunc(cpftri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong,
                  Clong), transr, uplo, n, a, info, 1, 1)
end

function cpftrs(transr, uplo, n, nrhs, a, b, ldb, info)
    return ccall((@blasfunc(cpftrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), transr, uplo, n,
                 nrhs, a, b, ldb, info, 1, 1)
end

function cpocon(uplo, n, a, lda, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(cpocon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, anorm, rcond, work, rwork, info, 1)
end

function cpoequ(n, a, lda, s, scond, amax, info)
    return ccall((@blasfunc(cpoequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}), n, a, lda, s, scond, amax, info)
end

function cpoequb(n, a, lda, s, scond, amax, info)
    return ccall((@blasfunc(cpoequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}), n, a, lda, s, scond, amax, info)
end

function cporfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, rwork,
                info)
    return ccall((@blasfunc(cporfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr,
                 berr, work, rwork, info, 1)
end

function cporfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr,
                 n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork,
                 info)
    return ccall((@blasfunc(cporfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), uplo, equed, n, nrhs, a, lda, af,
                 ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, rwork, info, 1, 1)
end

function cposv(uplo, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(cposv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, b,
                 ldb, info, 1)
end

function cposvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                ferr, berr, work, rwork, info)
    return ccall((@blasfunc(cposvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong,
                  Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx,
                 rcond, ferr, berr, work, rwork, info, 1, 1, 1)
end

function cposvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                 rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params,
                 work, rwork, info)
    return ccall((@blasfunc(cposvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong,
                  Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info, 1, 1, 1)
end

function cpotf2(uplo, n, a, lda, info)
    return ccall((@blasfunc(cpotf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function cpotrf(uplo, n, a, lda, info)
    return ccall((@blasfunc(cpotrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function cpotrf2(uplo, n, a, lda, info)
    return ccall((@blasfunc(cpotrf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function cpotri(uplo, n, a, lda, info)
    return ccall((@blasfunc(cpotri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function cpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(cpotrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, b,
                 ldb, info, 1)
end

function cppcon(uplo, n, ap, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(cppcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{Float32}, Ref{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, ap, anorm,
                 rcond, work, rwork, info, 1)
end

function cppequ(uplo, n, ap, s, scond, amax, info)
    return ccall((@blasfunc(cppequ_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Clong), uplo, n, ap, s, scond, amax, info, 1)
end

function cpprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(cpprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n,
                 nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function cppsv(uplo, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(cppsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, b, ldb, info, 1)
end

function cppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr,
                work, rwork, info)
    return ccall((@blasfunc(cppsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{UInt8}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong), fact,
                 uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work,
                 rwork, info, 1, 1, 1)
end

function cpptrf(uplo, n, ap, info)
    return ccall((@blasfunc(cpptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, ap,
                 info, 1)
end

function cpptri(uplo, n, ap, info)
    return ccall((@blasfunc(cpptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, ap,
                 info, 1)
end

function cpptrs(uplo, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(cpptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, b, ldb, info, 1)
end

function cpstf2(uplo, n, a, lda, piv, rank, tol, work, info)
    return ccall((@blasfunc(cpstf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, a, lda, piv, rank,
                 tol, work, info, 1)
end

function cpstrf(uplo, n, a, lda, piv, rank, tol, work, info)
    return ccall((@blasfunc(cpstrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, a, lda, piv, rank,
                 tol, work, info, 1)
end

function cptcon(n, d, e, anorm, rcond, rwork, info)
    return ccall((@blasfunc(cptcon_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}), n, d, e, anorm, rcond, rwork, info)
end

function cpteqr(compz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(cpteqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work,
                 info, 1)
end

function cptrfs(uplo, n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(cptrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, nrhs, d, e, df, ef, b, ldb, x,
                 ldx, ferr, berr, work, rwork, info, 1)
end

function cptsv(n, nrhs, d, e, b, ldb, info)
    return ccall((@blasfunc(cptsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}), n, nrhs, d, e, b, ldb, info)
end

function cptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, rwork,
                info)
    return ccall((@blasfunc(cptsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong), fact, n, nrhs, d, e, df,
                 ef, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1)
end

function cpttrf(n, d, e, info)
    return ccall((@blasfunc(cpttrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}), n, d, e, info)
end

function cpttrs(uplo, n, nrhs, d, e, b, ldb, info)
    return ccall((@blasfunc(cpttrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, d, e, b,
                 ldb, info, 1)
end

function cptts2(iuplo, n, nrhs, d, e, b, ldb)
    return ccall((@blasfunc(cptts2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}), iuplo, n, nrhs, d, e, b, ldb)
end

function crot(n, cx, incx, cy, incy, c, s)
    return ccall((@blasfunc(crot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{ComplexF32}), n, cx, incx, cy, incy, c, s)
end

function crscl(n, a, x, incx)
    return ccall((@blasfunc(crscl_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}), n, a, x, incx)
end

function cspcon(uplo, n, ap, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(cspcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, ap, ipiv,
                 anorm, rcond, work, info, 1)
end

function cspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    return ccall((@blasfunc(cspmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 alpha, ap, x, incx, beta, y, incy, 1)
end

function cspr(uplo, n, alpha, x, incx, ap)
    return ccall((@blasfunc(cspr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Clong), uplo, n, alpha, x, incx, ap, 1)
end

function csprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(csprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt},
                  Clong), uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 rwork, info, 1)
end

function cspsv(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(cspsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b,
                 ldb, info, 1)
end

function cspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(cspsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), fact, uplo, n, nrhs, ap, afp,
                 ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1, 1)
end

function csptrf(uplo, n, ap, ipiv, info)
    return ccall((@blasfunc(csptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, ap, ipiv, info, 1)
end

function csptri(uplo, n, ap, ipiv, work, info)
    return ccall((@blasfunc(csptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong), uplo, n, ap, ipiv, work, info, 1)
end

function csptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(csptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b,
                 ldb, info, 1)
end

function csrscl(n, sa, sx, incx)
    return ccall((@blasfunc(csrscl_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}), n, sa, sx, incx)
end

function cstedc(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info)
    return ccall((@blasfunc(cstedc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work, lwork, rwork,
                 lrwork, iwork, liwork, info, 1)
end

function cstegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
                lwork, iwork, liwork, info)
    return ccall((@blasfunc(cstegr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n,
                 d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork,
                 liwork, info, 1, 1)
end

function cstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info)
    return ccall((@blasfunc(cstein_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail,
                 info)
end

function cstemr(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac,
                work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(cstemr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n,
                 d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork,
                 iwork, liwork, info, 1, 1)
end

function csteqr(compz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(csteqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work,
                 info, 1)
end

function csycon(uplo, n, a, lda, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(csycon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, ipiv, anorm, rcond, work, info, 1)
end

function csycon_3(uplo, n, a, lda, e, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(csycon_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, e, ipiv, anorm, rcond, work, info, 1)
end

function csycon_rook(uplo, n, a, lda, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(csycon_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, ipiv, anorm, rcond, work, info, 1)
end

function csyconv(uplo, way, n, a, lda, ipiv, e, info)
    return ccall((@blasfunc(csyconv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a, lda, ipiv, e,
                 info, 1, 1)
end

function csyconvf(uplo, way, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(csyconvf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a,
                 lda, e, ipiv, info, 1, 1)
end

function csyconvf_rook(uplo, way, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(csyconvf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a,
                 lda, e, ipiv, info, 1, 1)
end

function csyequb(uplo, n, a, lda, s, scond, amax, work, info)
    return ccall((@blasfunc(csyequb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, s, scond, amax, work, info, 1)
end

function csymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(csymv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Clong), uplo, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function csyr(uplo, n, alpha, x, incx, a, lda)
    return ccall((@blasfunc(csyr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, alpha, x, incx, a, lda, 1)
end

function csyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(csyrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, af, ldaf, ipiv, b,
                 ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function csyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 rwork, info)
    return ccall((@blasfunc(csyrfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), uplo, equed, n,
                 nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds,
                 err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, 1, 1)
end

function csysv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(csysv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function csysv_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(csysv_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function csysv_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork,
                         info)
    return ccall((@blasfunc(csysv_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info, 1)
end

function csysv_rk(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(csysv_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, e, ipiv, b, ldb,
                 work, lwork, info, 1)
end

function csysv_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(csysv_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function csysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr,
                berr, work, lwork, rwork, info)
    return ccall((@blasfunc(csysvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), fact,
                 uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr,
                 work, lwork, rwork, info, 1, 1)
end

function csysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info)
    return ccall((@blasfunc(csysvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong, Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b,
                 ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, rwork, info, 1, 1, 1)
end

function csyswapr(uplo, n, a, lda, i1, i2)
    return ccall((@blasfunc(csyswapr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, i1, i2, 1)
end

function csytf2(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(csytf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function csytf2_rk(uplo, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(csytf2_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv, info, 1)
end

function csytf2_rook(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(csytf2_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function csytrf(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(csytrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function csytrf_aa(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(csytrf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function csytrf_aa_2stage(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info)
    return ccall((@blasfunc(csytrf_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info, 1)
end

function csytrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(csytrf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, lwork, info, 1)
end

function csytrf_rook(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(csytrf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function csytri(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(csytri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function csytri2(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(csytri2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function csytri2x(uplo, n, a, lda, ipiv, work, nb, info)
    return ccall((@blasfunc(csytri2x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, nb, info, 1)
end

function csytri_3(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(csytri_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, lwork, info, 1)
end

function csytri_3x(uplo, n, a, lda, e, ipiv, work, nb, info)
    return ccall((@blasfunc(csytri_3x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, nb, info, 1)
end

function csytri_rook(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(csytri_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function csytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(csytrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function csytrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info)
    return ccall((@blasfunc(csytrs2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n,
                 nrhs, a, lda, ipiv, b, ldb, work, info, 1)
end

function csytrs_3(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info)
    return ccall((@blasfunc(csytrs_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info, 1)
end

function csytrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(csytrs_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function csytrs_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info)
    return ccall((@blasfunc(csytrs_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2,
                 b, ldb, info, 1)
end

function csytrs_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(csytrs_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function ctbcon(norm, uplo, diag, n, kd, ab, ldab, rcond, work, rwork, info)
    return ccall((@blasfunc(ctbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong, Clong), norm, uplo, diag, n, kd, ab, ldab, rcond, work, rwork,
                 info, 1, 1, 1)
end

function ctbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(ctbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, kd, nrhs, ab, ldab,
                 b, ldb, x, ldx, ferr, berr, work, rwork, info, 1, 1, 1)
end

function ctbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(ctbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info, 1,
                 1, 1)
end

function ctfsm(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb)
    return ccall((@blasfunc(ctfsm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong, Clong), transr, side, uplo, trans, diag, m, n,
                 alpha, a, b, ldb, 1, 1, 1, 1, 1)
end

function ctftri(transr, uplo, diag, n, a, info)
    return ccall((@blasfunc(ctftri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Clong, Clong, Clong), transr, uplo, diag, n, a, info, 1, 1, 1)
end

function ctfttp(transr, uplo, n, arf, ap, info)
    return ccall((@blasfunc(ctfttp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong, Clong), transr, uplo, n, arf, ap, info, 1, 1)
end

function ctfttr(transr, uplo, n, arf, a, lda, info)
    return ccall((@blasfunc(ctfttr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), transr, uplo, n, arf, a, lda, info,
                 1, 1)
end

function ctgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, m, work,
                rwork, info)
    return ccall((@blasfunc(ctgevc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr,
                 ldvr, mm, m, work, rwork, info, 1, 1)
end

function ctgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, info)
    return ccall((@blasfunc(ctgex2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}), wantq, wantz, n, a, lda, b, ldb, q, ldq,
                 z, ldz, j1, info)
end

function ctgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, info)
    return ccall((@blasfunc(ctgexc_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}), wantq, wantz, n, a, lda, b,
                 ldb, q, ldq, z, ldz, ifst, ilst, info)
end

function ctgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz,
                m, pl, pr, dif, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(ctgsen_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}), ijob, wantq, wantz, select, n, a, lda,
                 b, ldb, alpha, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork,
                 liwork, info)
end

function ctgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u,
                ldu, v, ldv, q, ldq, work, ncycle, info)
    return ccall((@blasfunc(ctgsja_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta,
                 u, ldu, v, ldv, q, ldq, work, ncycle, info, 1, 1, 1)
end

function ctgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m,
                work, lwork, iwork, info)
    return ccall((@blasfunc(ctgsna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), job,
                 howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work,
                 lwork, iwork, info, 1, 1)
end

function ctgsy2(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                rdsum, rdscal, info)
    return ccall((@blasfunc(ctgsy2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Clong), trans, ijob,
                 m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal,
                 info, 1)
end

function ctgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                dif, work, lwork, iwork, info)
    return ccall((@blasfunc(ctgsyl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong), trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e,
                 lde, f, ldf, scale, dif, work, lwork, iwork, info, 1)
end

function ctpcon(norm, uplo, diag, n, ap, rcond, work, rwork, info)
    return ccall((@blasfunc(ctpcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong,
                  Clong), norm, uplo, diag, n, ap, rcond, work, rwork, info, 1, 1, 1)
end

function ctplqt(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info)
    return ccall((@blasfunc(ctplqt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}), m, n, l, mb, a, lda, b, ldb, t, ldt, work, info)
end

function ctplqt2(m, n, l, a, lda, b, ldb, t, ldt, info)
    return ccall((@blasfunc(ctplqt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 l, a, lda, b, ldb, t, ldt, info)
end

function ctpmlqt(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info)
    return ccall((@blasfunc(ctpmlqt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, l, mb, v, ldv, t, ldt, a,
                 lda, b, ldb, work, info, 1, 1)
end

function ctpmqrt(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info)
    return ccall((@blasfunc(ctpmqrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, l, nb, v, ldv, t, ldt, a,
                 lda, b, ldb, work, info, 1, 1)
end

function ctpqrt(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info)
    return ccall((@blasfunc(ctpqrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}), m, n, l, nb, a, lda, b, ldb, t, ldt, work, info)
end

function ctpqrt2(m, n, l, a, lda, b, ldb, t, ldt, info)
    return ccall((@blasfunc(ctpqrt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 l, a, lda, b, ldb, t, ldt, info)
end

function ctprfb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb,
                work, ldwork)
    return ccall((@blasfunc(ctprfb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong, Clong, Clong), side, trans,
                 direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork,
                 1, 1, 1, 1)
end

function ctprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work, rwork,
                info)
    return ccall((@blasfunc(ctprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong,
                  Clong), uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work,
                 rwork, info, 1, 1, 1)
end

function ctptri(uplo, diag, n, ap, info)
    return ccall((@blasfunc(ctptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong,
                  Clong), uplo, diag, n, ap, info, 1, 1)
end

function ctptrs(uplo, trans, diag, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(ctptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, nrhs, ap, b, ldb, info, 1, 1, 1)
end

function ctpttf(transr, uplo, n, ap, arf, info)
    return ccall((@blasfunc(ctpttf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong, Clong), transr, uplo, n, ap, arf, info, 1, 1)
end

function ctpttr(uplo, n, ap, a, lda, info)
    return ccall((@blasfunc(ctpttr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, ap, a, lda, info, 1)
end

function ctrcon(norm, uplo, diag, n, a, lda, rcond, work, rwork, info)
    return ccall((@blasfunc(ctrcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong,
                  Clong), norm, uplo, diag, n, a, lda, rcond, work, rwork, info, 1, 1, 1)
end

function ctrevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, rwork,
                info)
    return ccall((@blasfunc(ctrevc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), side,
                 howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, rwork, info, 1,
                 1)
end

function ctrevc3(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork,
                 rwork, lrwork, info)
    return ccall((@blasfunc(ctrevc3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m,
                 work, lwork, rwork, lrwork, info, 1, 1)
end

function ctrexc(compq, n, t, ldt, q, ldq, ifst, ilst, info)
    return ccall((@blasfunc(ctrexc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong), compq, n, t, ldt, q,
                 ldq, ifst, ilst, info, 1)
end

function ctrrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr, work, rwork,
                info)
    return ccall((@blasfunc(ctrrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong, Clong), uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx,
                 ferr, berr, work, rwork, info, 1, 1, 1)
end

function ctrsen(job, compq, select, n, t, ldt, q, ldq, w, m, s, sep, work, lwork, info)
    return ccall((@blasfunc(ctrsen_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job,
                 compq, select, n, t, ldt, q, ldq, w, m, s, sep, work, lwork, info, 1, 1)
end

function ctrsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, m, work,
                ldwork, rwork, info)
    return ccall((@blasfunc(ctrsna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), job, howmny, select, n, t, ldt,
                 vl, ldvl, vr, ldvr, s, sep, mm, m, work, ldwork, rwork, info, 1, 1)
end

function ctrsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info)
    return ccall((@blasfunc(ctrsyl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ref{BlasInt}, Clong, Clong), trana, tranb, isgn, m, n, a, lda,
                 b, ldb, c, ldc, scale, info, 1, 1)
end

function ctrsyl3(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, swork, ldswork,
                 info)
    return ccall((@blasfunc(ctrsyl3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), trana,
                 tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, swork, ldswork, info, 1,
                 1)
end

function ctrti2(uplo, diag, n, a, lda, info)
    return ccall((@blasfunc(ctrti2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), uplo, diag, n, a, lda, info, 1, 1)
end

function ctrtri(uplo, diag, n, a, lda, info)
    return ccall((@blasfunc(ctrtri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), uplo, diag, n, a, lda, info, 1, 1)
end

function ctrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(ctrtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, trans, diag, n, nrhs, a, lda, b, ldb, info, 1, 1, 1)
end

function ctrttf(transr, uplo, n, a, lda, arf, info)
    return ccall((@blasfunc(ctrttf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), transr, uplo, n, a, lda, arf,
                 info, 1, 1)
end

function ctrttp(uplo, n, a, lda, ap, info)
    return ccall((@blasfunc(ctrttp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong), uplo, n, a, lda, ap, info, 1)
end

function ctzrzf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(ctzrzf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function cunbdb(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22,
                theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info)
    return ccall((@blasfunc(cunbdb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), trans, signs, m, p, q, x11, ldx11,
                 x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2,
                 work, lwork, info, 1, 1)
end

function cunbdb1(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(cunbdb1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}),
                 m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
end

function cunbdb2(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(cunbdb2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}),
                 m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
end

function cunbdb3(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(cunbdb3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}),
                 m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
end

function cunbdb4(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom,
                 work, lwork, info)
    return ccall((@blasfunc(cunbdb4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}), m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1,
                 taup2, tauq1, phantom, work, lwork, info)
end

function cunbdb5(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
    return ccall((@blasfunc(cunbdb5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m1, m2, n, x1, incx1,
                 x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
end

function cunbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
    return ccall((@blasfunc(cunbdb6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m1, m2, n, x1, incx1,
                 x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
end

function cuncsd(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12,
                x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t,
                work, lwork, rwork, lrwork, iwork, info)
    return ccall((@blasfunc(cuncsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong, Clong, Clong), jobu1, jobu2, jobv1t, jobv2t,
                 trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22,
                 theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, rwork,
                 lrwork, iwork, info, 1, 1, 1, 1, 1, 1)
end

function cuncsd2by1(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1,
                    u2, ldu2, v1t, ldv1t, work, lwork, rwork, lrwork, iwork, info)
    return ccall((@blasfunc(cuncsd2by1_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), jobu1, jobu2, jobv1t, m, p, q, x11,
                 ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork,
                 rwork, lrwork, iwork, info, 1, 1, 1)
end

function cung2l(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(cung2l_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, k, a, lda, tau, work,
                 info)
end

function cung2r(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(cung2r_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, k, a, lda, tau, work,
                 info)
end

function cungbr(vect, m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cungbr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), vect, m,
                 n, k, a, lda, tau, work, lwork, info, 1)
end

function cunghr(n, ilo, ihi, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cunghr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), n, ilo, ihi, a,
                 lda, tau, work, lwork, info)
end

function cungl2(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(cungl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, k, a, lda, tau, work,
                 info)
end

function cunglq(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cunglq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda,
                 tau, work, lwork, info)
end

function cungql(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cungql_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda,
                 tau, work, lwork, info)
end

function cungqr(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cungqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda,
                 tau, work, lwork, info)
end

function cungr2(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(cungr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, k, a, lda, tau, work,
                 info)
end

function cungrq(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cungrq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda,
                 tau, work, lwork, info)
end

function cungtr(uplo, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(cungtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, tau, work,
                 lwork, info, 1)
end

function cungtsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(cungtsqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 mb, nb, a, lda, t, ldt, work, lwork, info)
end

function cungtsqr_row(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(cungtsqr_row_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 mb, nb, a, lda, t, ldt, work, lwork, info)
end

function cunhr_col(m, n, nb, a, lda, t, ldt, d, info)
    return ccall((@blasfunc(cunhr_col_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, nb, a, lda,
                 t, ldt, d, info)
end

function cunm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cunm22_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, n1, n2, q, ldq, c,
                 ldc, work, lwork, info, 1, 1)
end

function cunm2l(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(cunm2l_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 info, 1, 1)
end

function cunm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(cunm2r_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 info, 1, 1)
end

function cunmbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cunmbr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), vect, side,
                 trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info, 1, 1, 1)
end

function cunmhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cunmhr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 ilo, ihi, a, lda, tau, c, ldc, work, lwork, info, 1, 1)
end

function cunml2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(cunml2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 info, 1, 1)
end

function cunmlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cunmlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1)
end

function cunmql(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cunmql_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1)
end

function cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cunmqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1)
end

function cunmr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(cunmr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 info, 1, 1)
end

function cunmr3(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(cunmr3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, l, a,
                 lda, tau, c, ldc, work, info, 1, 1)
end

function cunmrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cunmrq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1)
end

function cunmrz(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cunmrz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt},
                  Ptr{ComplexF32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, l, a, lda, tau, c, ldc, work, lwork, info, 1, 1)
end

function cunmtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(cunmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), side, uplo, trans, m, n, a,
                 lda, tau, c, ldc, work, lwork, info, 1, 1, 1)
end

function cupgtr(uplo, n, ap, tau, q, ldq, work, info)
    return ccall((@blasfunc(cupgtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF32}, Ptr{ComplexF32}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}, Clong), uplo, n, ap, tau, q, ldq,
                 work, info, 1)
end

function cupmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info)
    return ccall((@blasfunc(cupmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ptr{ComplexF32}, Ptr{ComplexF32}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt},
                  Clong, Clong, Clong), side, uplo, trans, m, n, ap, tau, c, ldc, work,
                 info, 1, 1, 1)
end

function dbbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2,
                ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d,
                b22e, work, lwork, info)
    return ccall((@blasfunc(dbbcsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong), jobu1, jobu2, jobv1t,
                 jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t,
                 ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, work, lwork, info,
                 1, 1, 1, 1, 1)
end

function dbdsdc(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info)
    return ccall((@blasfunc(dbdsdc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, compq, n, d, e,
                 u, ldu, vt, ldvt, q, iq, work, iwork, info, 1, 1)
end

function dbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info)
    return ccall((@blasfunc(dbdsqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, ncvt,
                 nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info, 1)
end

function dbdsvdx(uplo, jobz, range, n, d, e, vl, vu, il, iu, ns, s, z, ldz, work, iwork,
                 info)
    return ccall((@blasfunc(dbdsvdx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), uplo, jobz, range, n, d, e, vl, vu, il, iu, ns, s, z, ldz, work,
                 iwork, info, 1, 1, 1)
end

function ddisna(job, m, n, d, sep, info)
    return ccall((@blasfunc(ddisna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Clong), job, m, n, d, sep, info, 1)
end

function dgbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work,
                info)
    return ccall((@blasfunc(dgbbrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Clong), vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt,
                 ldpt, c, ldc, work, info, 1)
end

function dgbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dgbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, iwork, info,
                 1)
end

function dgbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(dgbequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{BlasInt}), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
end

function dgbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(dgbequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{BlasInt}), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
end

function dgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr,
                berr, work, iwork, info)
    return ccall((@blasfunc(dgbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), trans, n, kl, ku, nrhs, ab, ldab, afb,
                 ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, 1)
end

function dgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x,
                 ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, iwork, info)
    return ccall((@blasfunc(dgbrfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c,
                 b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, iwork, info, 1, 1)
end

function dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    return ccall((@blasfunc(dgbsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), n, kl, ku, nrhs, ab, ldab,
                 ipiv, b, ldb, info)
end

function dgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
                ldb, x, ldx, rcond, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(dgbsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, kl, ku, nrhs,
                 ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr,
                 work, iwork, info, 1, 1, 1)
end

function dgbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
                 ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, iwork, info)
    return ccall((@blasfunc(dgbsvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, kl, ku, nrhs,
                 ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 iwork, info, 1, 1, 1)
end

function dgbtf2(m, n, kl, ku, ab, ldab, ipiv, info)
    return ccall((@blasfunc(dgbtf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, kl, ku, ab, ldab, ipiv, info)
end

function dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
    return ccall((@blasfunc(dgbtrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, kl, ku, ab, ldab, ipiv, info)
end

function dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    return ccall((@blasfunc(dgbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), trans,
                 n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info, 1)
end

function dgebak(job, side, n, ilo, ihi, scale, m, v, ldv, info)
    return ccall((@blasfunc(dgebak_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job, side,
                 n, ilo, ihi, scale, m, v, ldv, info, 1, 1)
end

function dgebal(job, n, a, lda, ilo, ihi, scale, info)
    return ccall((@blasfunc(dgebal_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong), job, n, a, lda, ilo, ihi, scale, info, 1)
end

function dgebd2(m, n, a, lda, d, e, tauq, taup, work, info)
    return ccall((@blasfunc(dgebd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), m, n, a, lda, d, e,
                 tauq, taup, work, info)
end

function dgebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info)
    return ccall((@blasfunc(dgebrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a,
                 lda, d, e, tauq, taup, work, lwork, info)
end

function dgecon(norm, n, a, lda, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dgecon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), norm, n, a, lda,
                 anorm, rcond, work, iwork, info, 1)
end

function dgeequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(dgeequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}), m, n, a, lda, r, c,
                 rowcnd, colcnd, amax, info)
end

function dgeequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(dgeequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}), m, n, a, lda, r, c,
                 rowcnd, colcnd, amax, info)
end

function dgees(jobvs, sort, select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork,
               info)
    return ccall((@blasfunc(dgees_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), jobvs, sort,
                 select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info, 1, 1)
end

function dgeesx(jobvs, sort, select, sense, n, a, lda, sdim, wr, wi, vs, ldvs, rconde,
                rcondv, work, lwork, iwork, liwork, bwork, info)
    return ccall((@blasfunc(dgeesx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobvs, sort, select, sense, n,
                 a, lda, sdim, wr, wi, vs, ldvs, rconde, rcondv, work, lwork, iwork, liwork,
                 bwork, info, 1, 1, 1)
end

function dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    return ccall((@blasfunc(dgeev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobvl, jobvr, n, a,
                 lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info, 1, 1)
end

function dgeevx(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo,
                ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info)
    return ccall((@blasfunc(dgeevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), balanc, jobvl, jobvr, sense, n, a, lda, wr,
                 wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work,
                 lwork, iwork, info, 1, 1, 1, 1)
end

function dgehd2(n, ilo, ihi, a, lda, tau, work, info)
    return ccall((@blasfunc(dgehd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}), n, ilo, ihi, a, lda, tau, work, info)
end

function dgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dgehrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), n, ilo, ihi, a, lda, tau, work,
                 lwork, info)
end

function dgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work,
                lwork, iwork, info)
    return ccall((@blasfunc(dgejsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong), joba, jobu, jobv,
                 jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work, lwork, iwork,
                 info, 1, 1, 1, 1, 1, 1)
end

function dgelq(m, n, a, lda, t, tsize, work, lwork, info)
    return ccall((@blasfunc(dgelq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, tsize, work, lwork,
                 info)
end

function dgelq2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(dgelq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function dgelqf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dgelqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function dgelqt(m, n, mb, a, lda, t, ldt, work, info)
    return ccall((@blasfunc(dgelqt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), m, n, mb, a, lda, t, ldt, work, info)
end

function dgelqt3(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(dgelqt3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(dgels_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
    return ccall((@blasfunc(dgelsd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work,
                 lwork, iwork, info)
end

function dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
    return ccall((@blasfunc(dgelss_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
end

function dgelst(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(dgelst_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
    return ccall((@blasfunc(dgelsy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork,
                 info)
end

function dgemlq(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dgemlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, t,
                 tsize, c, ldc, work, lwork, info, 1, 1)
end

function dgemlqt(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info)
    return ccall((@blasfunc(dgemlqt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, mb, v, ldv,
                 t, ldt, c, ldc, work, info, 1, 1)
end

function dgemqr(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dgemqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, t,
                 tsize, c, ldc, work, lwork, info, 1, 1)
end

function dgemqrt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info)
    return ccall((@blasfunc(dgemqrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, nb, v, ldv,
                 t, ldt, c, ldc, work, info, 1, 1)
end

function dgeql2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(dgeql2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function dgeqlf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dgeqlf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function dgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info)
    return ccall((@blasfunc(dgeqp3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, jpvt, tau, work, lwork,
                 info)
end

function dgeqp3rk(m, n, nrhs, kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk,
                  jpiv, tau, work, lwork, iwork, info)
    return ccall((@blasfunc(dgeqp3rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs,
                 kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk, jpiv, tau, work,
                 lwork, iwork, info)
end

function dgeqr(m, n, a, lda, t, tsize, work, lwork, info)
    return ccall((@blasfunc(dgeqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, tsize, work, lwork,
                 info)
end

function dgeqr2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(dgeqr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function dgeqr2p(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(dgeqr2p_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function dgeqrf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dgeqrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function dgeqrfp(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dgeqrfp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function dgeqrt(m, n, nb, a, lda, t, ldt, work, info)
    return ccall((@blasfunc(dgeqrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), m, n, nb, a, lda, t, ldt, work, info)
end

function dgeqrt2(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(dgeqrt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function dgeqrt3(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(dgeqrt3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function dgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(dgerfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1)
end

function dgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 iwork, info)
    return ccall((@blasfunc(dgerfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), trans,
                 equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr,
                 n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork,
                 info, 1, 1)
end

function dgerq2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(dgerq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function dgerqf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dgerqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function dgesc2(n, a, lda, rhs, ipiv, jpiv, scale)
    return ccall((@blasfunc(dgesc2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{Float64}), n, a, lda, rhs, ipiv, jpiv, scale)
end

function dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
    return ccall((@blasfunc(dgesdd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), jobz, m, n, a, lda, s, u, ldu, vt, ldvt,
                 work, lwork, iwork, info, 1)
end

function dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(dgesv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), n, nrhs, a, lda, ipiv, b, ldb, info)
end

function dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    return ccall((@blasfunc(dgesvd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobu, jobvt, m, n, a,
                 lda, s, u, ldu, vt, ldvt, work, lwork, info, 1, 1)
end

function dgesvdq(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank,
                 iwork, liwork, work, lwork, rwork, lrwork, info)
    return ccall((@blasfunc(dgesvdq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong,
                  Clong), joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv,
                 numrank, iwork, liwork, work, lwork, rwork, lrwork, info, 1, 1, 1, 1, 1)
end

function dgesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt,
                 work, lwork, iwork, info)
    return ccall((@blasfunc(dgesvdx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobu,
                 jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work,
                 lwork, iwork, info, 1, 1, 1)
end

function dgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info)
    return ccall((@blasfunc(dgesvj_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), joba, jobu, jobv, m, n, a,
                 lda, sva, mv, v, ldv, work, lwork, info, 1, 1, 1)
end

function dgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                rcond, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(dgesvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, nrhs, a, lda, af, ldaf,
                 ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, 1,
                 1, 1)
end

function dgesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, iwork, info)
    return ccall((@blasfunc(dgesvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed,
                 r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, iwork, info, 1, 1, 1)
end

function dgetc2(n, a, lda, ipiv, jpiv, info)
    return ccall((@blasfunc(dgetc2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), n,
                 a, lda, ipiv, jpiv, info)
end

function dgetf2(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(dgetf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m,
                 n, a, lda, ipiv, info)
end

function dgetrf(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(dgetrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m,
                 n, a, lda, ipiv, info)
end

function dgetrf2(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(dgetrf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m,
                 n, a, lda, ipiv, info)
end

function dgetri(n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(dgetri_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), n, a, lda, ipiv, work, lwork, info)
end

function dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(dgetrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), trans, n, nrhs, a, lda, ipiv,
                 b, ldb, info, 1)
end

function dgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(dgetsls_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function dgetsqrhrt(m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(dgetsqrhrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}),
                 m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info)
end

function dggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info)
    return ccall((@blasfunc(dggbak_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info, 1, 1)
end

function dggbal(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info)
    return ccall((@blasfunc(dggbal_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Clong), job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info, 1)
end

function dgges(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta,
               vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info)
    return ccall((@blasfunc(dgges_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai,
                 beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info, 1, 1, 1)
end

function dgges3(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta,
                vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info)
    return ccall((@blasfunc(dgges3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai,
                 beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info, 1, 1, 1)
end

function dggesx(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alphar,
                alphai, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork,
                liwork, bwork, info)
    return ccall((@blasfunc(dggesx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong), jobvsl,
                 jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alphar, alphai, beta,
                 vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, bwork,
                 info, 1, 1, 1, 1)
end

function dggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr,
               work, lwork, info)
    return ccall((@blasfunc(dggev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                 ldvl, vr, ldvr, work, lwork, info, 1, 1)
end

function dggev3(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr,
                work, lwork, info)
    return ccall((@blasfunc(dggev3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                 ldvl, vr, ldvr, work, lwork, info, 1, 1)
end

function dggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv,
                work, lwork, iwork, bwork, info)
    return ccall((@blasfunc(dggevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong), balanc, jobvl, jobvr, sense, n, a,
                 lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale,
                 rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info, 1,
                 1, 1, 1)
end

function dggglm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info)
    return ccall((@blasfunc(dggglm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), n, m, p, a, lda, b, ldb, d, x, y, work, lwork,
                 info)
end

function dgghd3(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork,
                info)
    return ccall((@blasfunc(dgghd3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), compq,
                 compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info, 1,
                 1)
end

function dgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info)
    return ccall((@blasfunc(dgghrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), compq, compz, n, ilo, ihi, a, lda, b,
                 ldb, q, ldq, z, ldz, info, 1, 1)
end

function dgglse(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    return ccall((@blasfunc(dgglse_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, p, a, lda, b, ldb, c, d, x, work, lwork,
                 info)
end

function dggqrf(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info)
    return ccall((@blasfunc(dggqrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info)
end

function dggrqf(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info)
    return ccall((@blasfunc(dggrqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info)
end

function dggsvd3(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v,
                 ldv, q, ldq, work, lwork, iwork, info)
    return ccall((@blasfunc(dggsvd3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), jobu, jobv, jobq, m, n, p, k, l, a, lda,
                 b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info, 1,
                 1, 1)
end

function dggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v,
                 ldv, q, ldq, iwork, tau, work, lwork, info)
    return ccall((@blasfunc(dggsvp3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobu, jobv, jobq, m, p, n, a,
                 lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work,
                 lwork, info, 1, 1, 1)
end

function dgsvj0(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work,
                lwork, info)
    return ccall((@blasfunc(dgsvj0_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep,
                 work, lwork, info, 1)
end

function dgsvj1(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work,
                lwork, info)
    return ccall((@blasfunc(dgsvj1_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), jobv, m, n, n1, a, lda, d, sva, mv, v, ldv,
                 eps, sfmin, tol, nsweep, work, lwork, info, 1)
end

function dgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dgtcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), norm, n, dl, d, du, du2, ipiv, anorm, rcond,
                 work, iwork, info, 1)
end

function dgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr,
                berr, work, iwork, info)
    return ccall((@blasfunc(dgtrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), trans, n, nrhs,
                 dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1)
end

function dgtsv(n, nrhs, dl, d, du, b, ldb, info)
    return ccall((@blasfunc(dgtsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), n, nrhs, dl, d, du, b, ldb, info)
end

function dgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx,
                rcond, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(dgtsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb,
                 x, ldx, rcond, ferr, berr, work, iwork, info, 1, 1)
end

function dgttrf(n, dl, d, du, du2, ipiv, info)
    return ccall((@blasfunc(dgttrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}), n, dl, d, du, du2, ipiv, info)
end

function dgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info)
    return ccall((@blasfunc(dgttrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info, 1)
end

function dgtts2(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb)
    return ccall((@blasfunc(dgtts2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}), itrans, n, nrhs, dl, d,
                 du, du2, ipiv, b, ldb)
end

function dhgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q,
                ldq, z, ldz, work, lwork, info)
    return ccall((@blasfunc(dhgeqz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), job,
                 compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z,
                 ldz, work, lwork, info, 1, 1, 1)
end

function dhsein(side, eigsrc, initv, select, n, h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m,
                work, ifaill, ifailr, info)
    return ccall((@blasfunc(dhsein_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), side, eigsrc, initv, select,
                 n, h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, work, ifaill, ifailr, info,
                 1, 1, 1)
end

function dhseqr(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info)
    return ccall((@blasfunc(dhseqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job, compz, n, ilo,
                 ihi, h, ldh, wr, wi, z, ldz, work, lwork, info, 1, 1)
end

function dla_gbamv(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy)
    return ccall((@blasfunc(dla_gbamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}), trans, m, n, kl, ku, alpha, ab, ldab, x, incx,
                 beta, y, incy)
end

function dla_gbrcond(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, cmode, c, info, work,
                     iwork)
    return ccall((@blasfunc(dla_gbrcond_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Clong), trans, n, kl, ku, ab, ldab, afb, ldafb,
                 ipiv, cmode, c, info, work, iwork, 1)
end

function dla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb,
                             ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms,
                             err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond,
                             ithresh, rthresh, dz_ub, ignore_cwise, info)
    return ccall((@blasfunc(dla_gbrfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}), prec_type, trans_type,
                 n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy,
                 berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail,
                 rcond, ithresh, rthresh, dz_ub, ignore_cwise, info)
end

function dla_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb)
    return ccall((@blasfunc(dla_gbrpvgrw_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}), n, kl, ku, ncols, ab, ldab, afb, ldafb)
end

function dla_geamv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(dla_geamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}), trans, m,
                 n, alpha, a, lda, x, incx, beta, y, incy)
end

function dla_gercond(trans, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork)
    return ccall((@blasfunc(dla_gercond_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Clong), trans, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork, 1)
end

function dla_gerfsx_extended(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ,
                             c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb,
                             dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info)
    return ccall((@blasfunc(dla_gerfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv,
                 colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy,
                 y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info)
end

function dla_gerpvgrw(n, ncols, a, lda, af, ldaf)
    return ccall((@blasfunc(dla_gerpvgrw_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}),
                 n, ncols, a, lda, af, ldaf)
end

function dla_lin_berr(n, nz, nrhs, res, ayb, berr)
    return ccall((@blasfunc(dla_lin_berr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}), n, nz, nrhs, res, ayb, berr)
end

function dla_porcond(uplo, n, a, lda, af, ldaf, cmode, c, info, work, iwork)
    return ccall((@blasfunc(dla_porcond_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Clong), uplo,
                 n, a, lda, af, ldaf, cmode, c, info, work, iwork, 1)
end

function dla_porfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb,
                             y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res,
                             ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise,
                             info)
    return ccall((@blasfunc(dla_porfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y,
                 ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail,
                 rcond, ithresh, rthresh, dz_ub, ignore_cwise, info, 1)
end

function dla_porpvgrw(uplo, ncols, a, lda, af, ldaf, work)
    return ccall((@blasfunc(dla_porpvgrw_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Clong), uplo, ncols, a, lda, af, ldaf, work, 1)
end

function dla_syamv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(dla_syamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}), uplo, n, alpha, a, lda,
                 x, incx, beta, y, incy)
end

function dla_syrcond(uplo, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork)
    return ccall((@blasfunc(dla_syrcond_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Clong), uplo, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork, 1)
end

function dla_syrfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                             ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp,
                             res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                             ignore_cwise, info)
    return ccall((@blasfunc(dla_syrfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                 ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy,
                 y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info, 1)
end

function dla_syrpvgrw(uplo, n, info, a, lda, af, ldaf, ipiv, work)
    return ccall((@blasfunc(dla_syrpvgrw_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Clong), uplo, n, info, a, lda, af,
                 ldaf, ipiv, work, 1)
end

function dla_wwaddw(n, x, y, w)
    return ccall((@blasfunc(dla_wwaddw_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), n, x, y, w)
end

function dlabad(small, large)
    return ccall((@blasfunc(dlabad_), libblastrampoline), Cvoid, (Ref{Float64}, Ref{Float64}),
                 small, large)
end

function dlabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy)
    return ccall((@blasfunc(dlabrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}), m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y,
                 ldy)
end

function dlacn2(n, v, x, isgn, est, kase, isave)
    return ccall((@blasfunc(dlacn2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}), n, v, x, isgn, est, kase, isave)
end

function dlacon(n, v, x, isgn, est, kase)
    return ccall((@blasfunc(dlacon_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{Float64},
                  Ref{BlasInt}), n, v, x, isgn, est, kase)
end

function dlacpy(uplo, m, n, a, lda, b, ldb)
    return ccall((@blasfunc(dlacpy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Clong), uplo, m, n, a, lda, b, ldb, 1)
end

function dladiv(a, b, c, d, p, q)
    return ccall((@blasfunc(dladiv_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}), a, b, c, d, p, q)
end

function dladiv1(a, b, c, d, p, q)
    return ccall((@blasfunc(dladiv1_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}), a, b, c, d, p, q)
end

function dladiv2(a, b, c, d, r, t)
    return ccall((@blasfunc(dladiv2_), libblastrampoline), Float64,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}), a, b, c, d, r, t)
end

function dlae2(a, b, c, rt1, rt2)
    return ccall((@blasfunc(dlae2_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), a,
                 b, c, rt1, rt2)
end

function dlaebz(ijob, nitmax, n, mmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval,
                ab, c, mout, nab, work, iwork, info)
    return ccall((@blasfunc(dlaebz_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), ijob, nitmax, n, mmax, minp, nbmin,
                 abstol, reltol, pivmin, d, e, e2, nval, ab, c, mout, nab, work, iwork,
                 info)
end

function dlaed0(icompq, qsiz, n, d, e, q, ldq, qstore, ldqs, work, iwork, info)
    return ccall((@blasfunc(dlaed0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}),
                 icompq, qsiz, n, d, e, q, ldq, qstore, ldqs, work, iwork, info)
end

function dlaed1(n, d, q, ldq, indxq, rho, cutpnt, work, iwork, info)
    return ccall((@blasfunc(dlaed1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), n, d, q, ldq, indxq, rho,
                 cutpnt, work, iwork, info)
end

function dlaed2(k, n, n1, d, q, ldq, indxq, rho, z, dlambda, w, q2, indx, indxc, indxp,
                coltyp, info)
    return ccall((@blasfunc(dlaed2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), k,
                 n, n1, d, q, ldq, indxq, rho, z, dlambda, w, q2, indx, indxc, indxp,
                 coltyp, info)
end

function dlaed3(k, n, n1, d, q, ldq, rho, dlambda, q2, indx, ctot, w, s, info)
    return ccall((@blasfunc(dlaed3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), k, n, n1, d, q, ldq, rho, dlambda,
                 q2, indx, ctot, w, s, info)
end

function dlaed4(n, i, d, z, delta, rho, dlam, info)
    return ccall((@blasfunc(dlaed4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}), n, i, d, z, delta, rho, dlam,
                 info)
end

function dlaed5(i, d, z, delta, rho, dlam)
    return ccall((@blasfunc(dlaed5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}), i, d, z, delta, rho, dlam)
end

function dlaed6(kniter, orgati, rho, d, z, finit, tau, info)
    return ccall((@blasfunc(dlaed6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}), kniter, orgati, rho, d, z, finit,
                 tau, info)
end

function dlaed7(icompq, n, qsiz, tlvls, curlvl, curpbm, d, q, ldq, indxq, rho, cutpnt,
                qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, iwork, info)
    return ccall((@blasfunc(dlaed7_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), icompq, n, qsiz, tlvls,
                 curlvl, curpbm, d, q, ldq, indxq, rho, cutpnt, qstore, qptr, prmptr, perm,
                 givptr, givcol, givnum, work, iwork, info)
end

function dlaed8(icompq, k, n, qsiz, d, q, ldq, indxq, rho, cutpnt, z, dlambda, q2, ldq2, w,
                perm, givptr, givcol, givnum, indxp, indx, info)
    return ccall((@blasfunc(dlaed8_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), icompq, k, n, qsiz, d, q,
                 ldq, indxq, rho, cutpnt, z, dlambda, q2, ldq2, w, perm, givptr, givcol,
                 givnum, indxp, indx, info)
end

function dlaed9(k, kstart, kstop, n, d, q, ldq, rho, dlambda, w, s, lds, info)
    return ccall((@blasfunc(dlaed9_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), k, kstart, kstop, n, d, q, ldq, rho, dlambda, w, s,
                 lds, info)
end

function dlaeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z,
                ztemp, info)
    return ccall((@blasfunc(dlaeda_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}), n, tlvls, curlvl, curpbm, prmptr, perm, givptr,
                 givcol, givnum, q, qptr, z, ztemp, info)
end

function dlaein(rightv, noinit, n, h, ldh, wr, wi, vr, vi, b, ldb, work, eps3, smlnum,
                bignum, info)
    return ccall((@blasfunc(dlaein_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}),
                 rightv, noinit, n, h, ldh, wr, wi, vr, vi, b, ldb, work, eps3, smlnum,
                 bignum, info)
end

function dlaev2(a, b, c, rt1, rt2, cs1, sn1)
    return ccall((@blasfunc(dlaev2_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}), a, b, c, rt1, rt2, cs1, sn1)
end

function dlaexc(wantq, n, t, ldt, q, ldq, j1, n1, n2, work, info)
    return ccall((@blasfunc(dlaexc_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), wantq, n, t,
                 ldt, q, ldq, j1, n1, n2, work, info)
end

function dlag2(a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi)
    return ccall((@blasfunc(dlag2_), libblastrampoline), Cvoid,
                 (Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), a,
                 lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi)
end

function dlag2s(m, n, a, lda, sa, ldsa, info)
    return ccall((@blasfunc(dlag2s_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, a, lda, sa, ldsa, info)
end

function dlags2(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq)
    return ccall((@blasfunc(dlags2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}), upper, a1, a2, a3, b1, b2, b3,
                 csu, snu, csv, snv, csq, snq)
end

function dlagtf(n, a, lambda, b, c, tol, d, in, info)
    return ccall((@blasfunc(dlagtf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), n, a, lambda, b, c,
                 tol, d, in, info)
end

function dlagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb)
    return ccall((@blasfunc(dlagtm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Clong), trans, n, nrhs, alpha, dl, d, du, x, ldx,
                 beta, b, ldb, 1)
end

function dlagts(job, n, a, b, c, d, in, y, tol, info)
    return ccall((@blasfunc(dlagts_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{Float64}, Ref{BlasInt}), job, n,
                 a, b, c, d, in, y, tol, info)
end

function dlagv2(a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr)
    return ccall((@blasfunc(dlagv2_), libblastrampoline), Cvoid,
                 (Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}), a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr)
end

function dlahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info)
    return ccall((@blasfunc(dlahqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz,
                 ihiz, z, ldz, info)
end

function dlahr2(n, k, nb, a, lda, tau, t, ldt, y, ldy)
    return ccall((@blasfunc(dlahr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, k, nb, a, lda, tau,
                 t, ldt, y, ldy)
end

function dlaic1(job, j, x, sest, w, gamma, sestpr, s, c)
    return ccall((@blasfunc(dlaic1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), job, j, x, sest,
                 w, gamma, sestpr, s, c)
end

function dlaln2(ltrans, na, nw, smin, ca, a, lda, d1, d2, b, ldb, wr, wi, x, ldx, scale,
                xnorm, info)
    return ccall((@blasfunc(dlaln2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}), ltrans, na, nw, smin, ca, a, lda, d1, d2, b,
                 ldb, wr, wi, x, ldx, scale, xnorm, info)
end

function dlals0(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol,
                givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info)
    return ccall((@blasfunc(dlals0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}), icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx,
                 perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c,
                 s, work, info)
end

function dlalsa(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z,
                poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info)
    return ccall((@blasfunc(dlalsa_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), icompq, smlsiz, n, nrhs, b, ldb, bx,
                 ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm,
                 givnum, c, s, work, iwork, info)
end

function dlalsd(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, iwork, info)
    return ccall((@blasfunc(dlalsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong), uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work,
                 iwork, info, 1)
end

function dlamrg(n1, n2, a, dtrd1, dtrd2, index)
    return ccall((@blasfunc(dlamrg_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}), n1,
                 n2, a, dtrd1, dtrd2, index)
end

function dlamswlq(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dlamswlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans,
                 m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info, 1, 1)
end

function dlamtsqr(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dlamtsqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans,
                 m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info, 1, 1)
end

function dlaneg(n, d, lld, sigma, pivmin, r)
    return ccall((@blasfunc(dlaneg_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{BlasInt}), n, d, lld, sigma, pivmin, r)
end

function dlangb(norm, n, kl, ku, ab, ldab, work)
    return ccall((@blasfunc(dlangb_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Clong), norm, n, kl, ku, ab, ldab, work, 1)
end

function dlange(norm, m, n, a, lda, work)
    return ccall((@blasfunc(dlange_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Clong), norm, m, n, a, lda, work, 1)
end

function dlangt(norm, n, dl, d, du)
    return ccall((@blasfunc(dlangt_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Clong),
                 norm, n, dl, d, du, 1)
end

function dlanhs(norm, n, a, lda, work)
    return ccall((@blasfunc(dlanhs_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong),
                 norm, n, a, lda, work, 1)
end

function dlansb(norm, uplo, n, k, ab, ldab, work)
    return ccall((@blasfunc(dlansb_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Clong, Clong), norm, uplo, n, k, ab, ldab, work, 1, 1)
end

function dlansf(norm, transr, uplo, n, a, work)
    return ccall((@blasfunc(dlansf_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Clong, Clong, Clong), norm, transr, uplo, n, a, work, 1, 1, 1)
end

function dlansp(norm, uplo, n, ap, work)
    return ccall((@blasfunc(dlansp_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Clong,
                  Clong), norm, uplo, n, ap, work, 1, 1)
end

function dlanst(norm, n, d, e)
    return ccall((@blasfunc(dlanst_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Clong), norm, n, d, e,
                 1)
end

function dlansy(norm, uplo, n, a, lda, work)
    return ccall((@blasfunc(dlansy_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Clong, Clong), norm, uplo, n, a, lda, work, 1, 1)
end

function dlantb(norm, uplo, diag, n, k, ab, ldab, work)
    return ccall((@blasfunc(dlantb_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Clong, Clong, Clong), norm, uplo, diag, n, k, ab,
                 ldab, work, 1, 1, 1)
end

function dlantp(norm, uplo, diag, n, ap, work)
    return ccall((@blasfunc(dlantp_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Clong, Clong, Clong), norm, uplo, diag, n, ap, work, 1, 1, 1)
end

function dlantr(norm, uplo, diag, m, n, a, lda, work)
    return ccall((@blasfunc(dlantr_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Clong, Clong, Clong), norm, uplo, diag, m, n, a,
                 lda, work, 1, 1, 1)
end

function dlanv2(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn)
    return ccall((@blasfunc(dlanv2_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), a,
                 b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn)
end

function dlaorhr_col_getrfnp(m, n, a, lda, d, info)
    return ccall((@blasfunc(dlaorhr_col_getrfnp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}),
                 m, n, a, lda, d, info)
end

function dlaorhr_col_getrfnp2(m, n, a, lda, d, info)
    return ccall((@blasfunc(dlaorhr_col_getrfnp2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}),
                 m, n, a, lda, d, info)
end

function dlapll(n, x, incx, y, incy, ssmin)
    return ccall((@blasfunc(dlapll_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}), n, x, incx, y, incy, ssmin)
end

function dlapmr(forwrd, m, n, x, ldx, k)
    return ccall((@blasfunc(dlapmr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}),
                 forwrd, m, n, x, ldx, k)
end

function dlapmt(forwrd, m, n, x, ldx, k)
    return ccall((@blasfunc(dlapmt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}),
                 forwrd, m, n, x, ldx, k)
end

function dlapy2(x, y)
    return ccall((@blasfunc(dlapy2_), libblastrampoline), Float64, (Ref{Float64}, Ref{Float64}),
                 x, y)
end

function dlapy3(x, y, z)
    return ccall((@blasfunc(dlapy3_), libblastrampoline), Float64,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}), x, y, z)
end

function dlaqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed)
    return ccall((@blasfunc(dlaqgb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{UInt8}, Clong), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax,
                 equed, 1)
end

function dlaqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed)
    return ccall((@blasfunc(dlaqge_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{UInt8}, Clong), m, n, a,
                 lda, r, c, rowcnd, colcnd, amax, equed, 1)
end

function dlaqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work)
    return ccall((@blasfunc(dlaqp2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), m, n, offset, a,
                 lda, jpvt, tau, vn1, vn2, work)
end

function dlaqp2rk(m, n, nrhs, ioffset, kmax, abstol, reltol, kp1, maxc2nrm, a, lda, k,
                  maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2, work, info)
    return ccall((@blasfunc(dlaqp2rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), m, n, nrhs, ioffset, kmax, abstol,
                 reltol, kp1, maxc2nrm, a, lda, k, maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1,
                 vn2, work, info)
end

function dlaqp3rk(m, n, nrhs, ioffset, nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb,
                  maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info)
    return ccall((@blasfunc(dlaqp3rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs, ioffset, nb, abstol, reltol, kp1,
                 maxc2nrm, a, lda, done, kb, maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2,
                 auxv, f, ldf, iwork, info)
end

function dlaqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf)
    return ccall((@blasfunc(dlaqps_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), m, n, offset, nb, kb, a, lda,
                 jpvt, tau, vn1, vn2, auxv, f, ldf)
end

function dlaqr0(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork,
                info)
    return ccall((@blasfunc(dlaqr0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo, ihi,
                 h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info)
end

function dlaqr1(n, h, ldh, sr1, si1, sr2, si2, v)
    return ccall((@blasfunc(dlaqr1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}), n, h, ldh, sr1, si1, sr2, si2,
                 v)
end

function dlaqr2(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si,
                v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork)
    return ccall((@blasfunc(dlaqr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}), wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
                 ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work,
                 lwork)
end

function dlaqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si,
                v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork)
    return ccall((@blasfunc(dlaqr3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}), wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
                 ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work,
                 lwork)
end

function dlaqr4(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork,
                info)
    return ccall((@blasfunc(dlaqr4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo, ihi,
                 h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info)
end

function dlaqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z,
                ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh)
    return ccall((@blasfunc(dlaqr5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}), wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh,
                 iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh)
end

function dlaqsb(uplo, n, kd, ab, ldab, s, scond, amax, equed)
    return ccall((@blasfunc(dlaqsb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{UInt8}, Clong, Clong), uplo, n, kd, ab,
                 ldab, s, scond, amax, equed, 1, 1)
end

function dlaqsp(uplo, n, ap, s, scond, amax, equed)
    return ccall((@blasfunc(dlaqsp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{UInt8}, Clong, Clong), uplo, n, ap, s, scond, amax,
                 equed, 1, 1)
end

function dlaqsy(uplo, n, a, lda, s, scond, amax, equed)
    return ccall((@blasfunc(dlaqsy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{UInt8}, Clong, Clong), uplo, n, a, lda, s,
                 scond, amax, equed, 1, 1)
end

function dlaqtr(ltran, lreal, n, t, ldt, b, w, scale, x, work, info)
    return ccall((@blasfunc(dlaqtr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), ltran,
                 lreal, n, t, ldt, b, w, scale, x, work, info)
end

function dlar1v(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz,
                mingma, r, isuppz, nrminv, resid, rqcorr, work)
    return ccall((@blasfunc(dlar1v_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64}), n, b1, bn,
                 lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r,
                 isuppz, nrminv, resid, rqcorr, work)
end

function dlar2v(n, x, y, z, incx, c, s, incc)
    return ccall((@blasfunc(dlar2v_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), n, x, y, z, incx, c, s, incc)
end

function dlarf(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(dlarf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong), side, m, n, v, incv, tau,
                 c, ldc, work, 1)
end

function dlarf1f(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(dlarf1f_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong), side, m, n, v, incv, tau,
                 c, ldc, work, 1)
end

function dlarf1l(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(dlarf1l_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong), side, m, n, v, incv, tau,
                 c, ldc, work, 1)
end

function dlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork)
    return ccall((@blasfunc(dlarfb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong, Clong), side,
                 trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork, 1, 1,
                 1, 1)
end

function dlarfb_gett(ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork)
    return ccall((@blasfunc(dlarfb_gett_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Clong), ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork, 1)
end

function dlarfg(n, alpha, x, incx, tau)
    return ccall((@blasfunc(dlarfg_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}), n, alpha,
                 x, incx, tau)
end

function dlarfgp(n, alpha, x, incx, tau)
    return ccall((@blasfunc(dlarfgp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}), n, alpha,
                 x, incx, tau)
end

function dlarft(direct, storev, n, k, v, ldv, tau, t, ldt)
    return ccall((@blasfunc(dlarft_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), direct, storev, n,
                 k, v, ldv, tau, t, ldt, 1, 1)
end

function dlarfx(side, m, n, v, tau, c, ldc, work)
    return ccall((@blasfunc(dlarfx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong), side, m, n, v, tau, c, ldc,
                 work, 1)
end

function dlarfy(uplo, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(dlarfy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong), uplo, n, v, incv, tau, c,
                 ldc, work, 1)
end

function dlargv(n, x, incx, y, incy, c, incc)
    return ccall((@blasfunc(dlargv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}), n, x, incx, y, incy, c, incc)
end

function dlarmm(anorm, bnorm, cnorm)
    return ccall((@blasfunc(dlarmm_), libblastrampoline), Float64,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}), anorm, bnorm, cnorm)
end

function dlarnv(idist, iseed, n, x)
    return ccall((@blasfunc(dlarnv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}), idist, iseed, n, x)
end

function dlarra(n, d, e, e2, spltol, tnrm, nsplit, isplit, info)
    return ccall((@blasfunc(dlarra_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), n, d, e, e2, spltol, tnrm,
                 nsplit, isplit, info)
end

function dlarrb(n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, work, iwork,
                pivmin, spdiam, twist, info)
    return ccall((@blasfunc(dlarrb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr,
                 work, iwork, pivmin, spdiam, twist, info)
end

function dlarrc(jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info)
    return ccall((@blasfunc(dlarrc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info, 1)
end

function dlarrd(range, order, n, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit,
                isplit, m, w, werr, wl, wu, iblock, indexw, work, iwork, info)
    return ccall((@blasfunc(dlarrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), range, order, n, vl,
                 vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w, werr, wl,
                 wu, iblock, indexw, work, iwork, info, 1, 1)
end

function dlarre(range, n, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m,
                w, werr, wgap, iblock, indexw, gers, pivmin, work, iwork, info)
    return ccall((@blasfunc(dlarre_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), range, n, vl, vu, il, iu, d,
                 e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr, wgap, iblock,
                 indexw, gers, pivmin, work, iwork, info, 1)
end

function dlarrf(n, d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin,
                sigma, dplus, lplus, work, info)
    return ccall((@blasfunc(dlarrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}), n, d, l, ld, clstrt, clend, w, wgap, werr,
                 spdiam, clgapl, clgapr, pivmin, sigma, dplus, lplus, work, info)
end

function dlarrj(n, d, e2, ifirst, ilast, rtol, offset, w, werr, work, iwork, pivmin, spdiam,
                info)
    return ccall((@blasfunc(dlarrj_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}), n, d, e2, ifirst, ilast, rtol,
                 offset, w, werr, work, iwork, pivmin, spdiam, info)
end

function dlarrk(n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info)
    return ccall((@blasfunc(dlarrk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{BlasInt}), n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info)
end

function dlarrr(n, d, e, info)
    return ccall((@blasfunc(dlarrr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), n, d, e, info)
end

function dlarrv(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr,
                wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info)
    return ccall((@blasfunc(dlarrv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), n, vl, vu, d, l, pivmin, isplit, m,
                 dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z,
                 ldz, isuppz, work, iwork, info)
end

function dlarscl2(m, n, d, x, ldx)
    return ccall((@blasfunc(dlarscl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), m, n, d, x,
                 ldx)
end

function dlartgp(f, g, cs, sn, r)
    return ccall((@blasfunc(dlartgp_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), f,
                 g, cs, sn, r)
end

function dlartgs(x, y, sigma, cs, sn)
    return ccall((@blasfunc(dlartgs_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), x,
                 y, sigma, cs, sn)
end

function dlartv(n, x, incx, y, incy, c, s, incc)
    return ccall((@blasfunc(dlartv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}), n, x, incx, y, incy, c, s, incc)
end

function dlaruv(iseed, n, x)
    return ccall((@blasfunc(dlaruv_), libblastrampoline), Cvoid,
                 (Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}), iseed, n, x)
end

function dlarz(side, m, n, l, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(dlarz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong), side, m, n,
                 l, v, incv, tau, c, ldc, work, 1)
end

function dlarzb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work,
                ldwork)
    return ccall((@blasfunc(dlarzb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc,
                 work, ldwork, 1, 1, 1, 1)
end

function dlarzt(direct, storev, n, k, v, ldv, tau, t, ldt)
    return ccall((@blasfunc(dlarzt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), direct, storev, n,
                 k, v, ldv, tau, t, ldt, 1, 1)
end

function dlas2(f, g, h, ssmin, ssmax)
    return ccall((@blasfunc(dlas2_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), f,
                 g, h, ssmin, ssmax)
end

function dlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info)
    return ccall((@blasfunc(dlascl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), type, kl, ku,
                 cfrom, cto, m, n, a, lda, info, 1)
end

function dlascl2(m, n, d, x, ldx)
    return ccall((@blasfunc(dlascl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), m, n, d, x,
                 ldx)
end

function dlasd0(n, sqre, d, e, u, ldu, vt, ldvt, smlsiz, iwork, work, info)
    return ccall((@blasfunc(dlasd0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}),
                 n, sqre, d, e, u, ldu, vt, ldvt, smlsiz, iwork, work, info)
end

function dlasd1(nl, nr, sqre, d, alpha, beta, u, ldu, vt, ldvt, idxq, iwork, work, info)
    return ccall((@blasfunc(dlasd1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}), nl, nr, sqre, d, alpha, beta, u, ldu, vt, ldvt,
                 idxq, iwork, work, info)
end

function dlasd2(nl, nr, sqre, k, d, z, alpha, beta, u, ldu, vt, ldvt, dsigma, u2, ldu2, vt2,
                ldvt2, idxp, idx, idxc, idxq, coltyp, info)
    return ccall((@blasfunc(dlasd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), nl, nr,
                 sqre, k, d, z, alpha, beta, u, ldu, vt, ldvt, dsigma, u2, ldu2, vt2, ldvt2,
                 idxp, idx, idxc, idxq, coltyp, info)
end

function dlasd3(nl, nr, sqre, k, d, q, ldq, dsigma, u, ldu, u2, ldu2, vt, ldvt, vt2, ldvt2,
                idxc, ctot, z, info)
    return ccall((@blasfunc(dlasd3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}), nl, nr, sqre, k, d, q, ldq, dsigma, u, ldu, u2,
                 ldu2, vt, ldvt, vt2, ldvt2, idxc, ctot, z, info)
end

function dlasd4(n, i, d, z, delta, rho, sigma, work, info)
    return ccall((@blasfunc(dlasd4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}), n, i, d, z, delta,
                 rho, sigma, work, info)
end

function dlasd5(i, d, z, delta, rho, dsigma, work)
    return ccall((@blasfunc(dlasd5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}), i, d, z, delta, rho, dsigma, work)
end

function dlasd6(icompq, nl, nr, sqre, d, vf, vl, alpha, beta, idxq, perm, givptr, givcol,
                ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, iwork, info)
    return ccall((@blasfunc(dlasd6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), icompq, nl, nr, sqre, d, vf, vl,
                 alpha, beta, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles,
                 difl, difr, z, k, c, s, work, iwork, info)
end

function dlasd7(icompq, nl, nr, sqre, k, d, z, zw, vf, vfw, vl, vlw, alpha, beta, dsigma,
                idx, idxp, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info)
    return ccall((@blasfunc(dlasd7_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}), icompq,
                 nl, nr, sqre, k, d, z, zw, vf, vfw, vl, vlw, alpha, beta, dsigma, idx,
                 idxp, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info)
end

function dlasd8(icompq, k, d, z, vf, vl, difl, difr, lddifr, dsigma, work, info)
    return ccall((@blasfunc(dlasd8_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}), icompq, k, d, z, vf, vl, difl, difr, lddifr,
                 dsigma, work, info)
end

function dlasda(icompq, smlsiz, n, sqre, d, e, u, ldu, vt, k, difl, difr, z, poles, givptr,
                givcol, ldgcol, perm, givnum, c, s, work, iwork, info)
    return ccall((@blasfunc(dlasda_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}), icompq, smlsiz, n, sqre, d, e, u, ldu, vt, k, difl,
                 difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork,
                 info)
end

function dlasdq(uplo, sqre, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info)
    return ccall((@blasfunc(dlasdq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo,
                 sqre, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info, 1)
end

function dlasdt(n, lvl, nd, inode, ndiml, ndimr, msub)
    return ccall((@blasfunc(dlasdt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, lvl, nd, inode, ndiml, ndimr, msub)
end

function dlaset(uplo, m, n, alpha, beta, a, lda)
    return ccall((@blasfunc(dlaset_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Clong), uplo, m, n, alpha, beta, a, lda, 1)
end

function dlasq1(n, d, e, work, info)
    return ccall((@blasfunc(dlasq1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), n, d, e,
                 work, info)
end

function dlasq2(n, z, info)
    return ccall((@blasfunc(dlasq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), n, z, info)
end

function dlasq3(i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype,
                dmin1, dmin2, dn, dn1, dn2, g, tau)
    return ccall((@blasfunc(dlasq3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}), i0, n0, z, pp, dmin, sigma,
                 desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g,
                 tau)
end

function dlasq4(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g)
    return ccall((@blasfunc(dlasq4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}), i0, n0, z, pp, n0in, dmin, dmin1,
                 dmin2, dn, dn1, dn2, tau, ttype, g)
end

function dlasq5(i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2, ieee, eps)
    return ccall((@blasfunc(dlasq5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}), i0, n0, z, pp, tau, sigma, dmin,
                 dmin1, dmin2, dn, dnm1, dnm2, ieee, eps)
end

function dlasq6(i0, n0, z, pp, dmin, dmin1, dmin2, dn, dnm1, dnm2)
    return ccall((@blasfunc(dlasq6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), i0, n0, z, pp,
                 dmin, dmin1, dmin2, dn, dnm1, dnm2)
end

function dlasr(side, pivot, direct, m, n, c, s, a, lda)
    return ccall((@blasfunc(dlasr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong), side, pivot,
                 direct, m, n, c, s, a, lda, 1, 1, 1)
end

function dlasrt(id, n, d, info)
    return ccall((@blasfunc(dlasrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), id, n, d, info, 1)
end

function dlasv2(f, g, h, ssmin, ssmax, snr, csr, snl, csl)
    return ccall((@blasfunc(dlasv2_), libblastrampoline), Cvoid,
                 (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}), f, g, h, ssmin,
                 ssmax, snr, csr, snl, csl)
end

function dlaswlq(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(dlaswlq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, mb,
                 nb, a, lda, t, ldt, work, lwork, info)
end

function dlaswp(n, a, lda, k1, k2, ipiv, incx)
    return ccall((@blasfunc(dlaswp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, a, lda, k1, k2, ipiv, incx)
end

function dlasy2(ltranl, ltranr, isgn, n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx,
                xnorm, info)
    return ccall((@blasfunc(dlasy2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}), ltranl, ltranr, isgn,
                 n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx, xnorm, info)
end

function dlasyf(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(dlasyf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb, a,
                 lda, ipiv, w, ldw, info, 1)
end

function dlasyf_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work)
    return ccall((@blasfunc(dlasyf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong), uplo, j1, m, nb,
                 a, lda, ipiv, h, ldh, work, 1)
end

function dlasyf_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info)
    return ccall((@blasfunc(dlasyf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nb, kb, a, lda, e, ipiv, w, ldw, info, 1)
end

function dlasyf_rook(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(dlasyf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb, a,
                 lda, ipiv, w, ldw, info, 1)
end

function dlat2s(uplo, n, a, lda, sa, ldsa, info)
    return ccall((@blasfunc(dlat2s_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, a, lda, sa, ldsa, info, 1)
end

function dlatbs(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info)
    return ccall((@blasfunc(dlatbs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong), uplo, trans, diag, normin, n, kd,
                 ab, ldab, x, scale, cnorm, info, 1, 1, 1, 1)
end

function dlatdf(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv)
    return ccall((@blasfunc(dlatdf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ptr{BlasInt}, Ptr{BlasInt}), ijob, n, z, ldz, rhs, rdsum, rdscal,
                 ipiv, jpiv)
end

function dlatps(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info)
    return ccall((@blasfunc(dlatps_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), uplo, trans, diag, normin, n, ap, x, scale, cnorm, info, 1, 1, 1,
                 1)
end

function dlatrd(uplo, n, nb, a, lda, e, tau, w, ldw)
    return ccall((@blasfunc(dlatrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, nb, a, lda, e,
                 tau, w, ldw, 1)
end

function dlatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info)
    return ccall((@blasfunc(dlatrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong), uplo, trans, diag, normin, n, a, lda, x, scale,
                 cnorm, info, 1, 1, 1, 1)
end

function dlatrs3(uplo, trans, diag, normin, n, nrhs, a, lda, x, ldx, scale, cnorm, work,
                 lwork, info)
    return ccall((@blasfunc(dlatrs3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), uplo, trans, diag, normin, n, nrhs, a, lda, x, ldx, scale, cnorm,
                 work, lwork, info, 1, 1, 1, 1)
end

function dlatrz(m, n, l, a, lda, tau, work)
    return ccall((@blasfunc(dlatrz_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}), m, n, l, a, lda, tau, work)
end

function dlatsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(dlatsqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, mb,
                 nb, a, lda, t, ldt, work, lwork, info)
end

function dlauu2(uplo, n, a, lda, info)
    return ccall((@blasfunc(dlauu2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function dlauum(uplo, n, a, lda, info)
    return ccall((@blasfunc(dlauum_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function dopgtr(uplo, n, ap, tau, q, ldq, work, info)
    return ccall((@blasfunc(dopgtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, ap, tau, q, ldq,
                 work, info, 1)
end

function dopmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info)
    return ccall((@blasfunc(dopmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong, Clong), side, uplo, trans, m, n, ap, tau, c, ldc, work, info, 1, 1,
                 1)
end

function dorbdb(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22,
                theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info)
    return ccall((@blasfunc(dorbdb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22,
                 ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info, 1, 1)
end

function dorbdb1(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(dorbdb1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, p, q, x11, ldx11,
                 x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info)
end

function dorbdb2(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(dorbdb2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, p, q, x11, ldx11,
                 x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info)
end

function dorbdb3(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(dorbdb3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, p, q, x11, ldx11,
                 x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info)
end

function dorbdb4(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom,
                 work, lwork, info)
    return ccall((@blasfunc(dorbdb4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, p, q,
                 x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work,
                 lwork, info)
end

function dorbdb5(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
    return ccall((@blasfunc(dorbdb5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2,
                 ldq2, work, lwork, info)
end

function dorbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
    return ccall((@blasfunc(dorbdb6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2,
                 ldq2, work, lwork, info)
end

function dorcsd(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12,
                x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t,
                work, lwork, iwork, info)
    return ccall((@blasfunc(dorcsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong, Clong, Clong), jobu1, jobu2, jobv1t, jobv2t,
                 trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22,
                 theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, iwork,
                 info, 1, 1, 1, 1, 1, 1)
end

function dorcsd2by1(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1,
                    u2, ldu2, v1t, ldv1t, work, lwork, iwork, info)
    return ccall((@blasfunc(dorcsd2by1_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2,
                 ldu2, v1t, ldv1t, work, lwork, iwork, info, 1, 1, 1)
end

function dorg2l(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(dorg2l_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}), m, n, k, a, lda, tau, work, info)
end

function dorg2r(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(dorg2r_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}), m, n, k, a, lda, tau, work, info)
end

function dorgbr(vect, m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dorgbr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), vect, m, n, k,
                 a, lda, tau, work, lwork, info, 1)
end

function dorghr(n, ilo, ihi, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dorghr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), n, ilo, ihi, a, lda, tau, work,
                 lwork, info)
end

function dorgl2(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(dorgl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}), m, n, k, a, lda, tau, work, info)
end

function dorglq(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dorglq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda, tau, work, lwork,
                 info)
end

function dorgql(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dorgql_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda, tau, work, lwork,
                 info)
end

function dorgqr(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dorgqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda, tau, work, lwork,
                 info)
end

function dorgr2(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(dorgr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}), m, n, k, a, lda, tau, work, info)
end

function dorgrq(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dorgrq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda, tau, work, lwork,
                 info)
end

function dorgtr(uplo, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dorgtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, tau, work,
                 lwork, info, 1)
end

function dorgtsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(dorgtsqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, mb,
                 nb, a, lda, t, ldt, work, lwork, info)
end

function dorgtsqr_row(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(dorgtsqr_row_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, mb,
                 nb, a, lda, t, ldt, work, lwork, info)
end

function dorhr_col(m, n, nb, a, lda, t, ldt, d, info)
    return ccall((@blasfunc(dorhr_col_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), m, n, nb, a, lda, t, ldt, d, info)
end

function dorm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dorm22_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, n1, n2, q, ldq, c, ldc, work,
                 lwork, info, 1, 1)
end

function dorm2l(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(dorm2l_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work, info, 1,
                 1)
end

function dorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(dorm2r_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work, info, 1,
                 1)
end

function dormbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dormbr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), vect, side,
                 trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info, 1, 1, 1)
end

function dormhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dormhr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, ilo,
                 ihi, a, lda, tau, c, ldc, work, lwork, info, 1, 1)
end

function dorml2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(dorml2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work, info, 1,
                 1)
end

function dormlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dormlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 lwork, info, 1, 1)
end

function dormql(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dormql_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 lwork, info, 1, 1)
end

function dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dormqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 lwork, info, 1, 1)
end

function dormr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(dormr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work, info, 1,
                 1)
end

function dormr3(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(dormr3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, l, a, lda,
                 tau, c, ldc, work, info, 1, 1)
end

function dormrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dormrq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 lwork, info, 1, 1)
end

function dormrz(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dormrz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k,
                 l, a, lda, tau, c, ldc, work, lwork, info, 1, 1)
end

function dormtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(dormtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), side, uplo, trans, m, n, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1, 1)
end

function dpbcon(uplo, n, kd, ab, ldab, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dpbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, ab,
                 ldab, anorm, rcond, work, iwork, info, 1)
end

function dpbequ(uplo, n, kd, ab, ldab, s, scond, amax, info)
    return ccall((@blasfunc(dpbequ_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Clong), uplo, n, kd, ab, ldab, s,
                 scond, amax, info, 1)
end

function dpbrfs(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(dpbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1)
end

function dpbstf(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(dpbstf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function dpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(dpbsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab, ldab,
                 b, ldb, info, 1)
end

function dpbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx,
                rcond, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(dpbsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b,
                 ldb, x, ldx, rcond, ferr, berr, work, iwork, info, 1, 1, 1)
end

function dpbtf2(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(dpbtf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function dpbtrf(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(dpbtrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function dpbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(dpbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab, ldab,
                 b, ldb, info, 1)
end

function dpftrf(transr, uplo, n, a, info)
    return ccall((@blasfunc(dpftrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong),
                 transr, uplo, n, a, info, 1, 1)
end

function dpftri(transr, uplo, n, a, info)
    return ccall((@blasfunc(dpftri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong),
                 transr, uplo, n, a, info, 1, 1)
end

function dpftrs(transr, uplo, n, nrhs, a, b, ldb, info)
    return ccall((@blasfunc(dpftrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), transr, uplo, n, nrhs, a, b, ldb,
                 info, 1, 1)
end

function dpocon(uplo, n, a, lda, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dpocon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 anorm, rcond, work, iwork, info, 1)
end

function dpoequ(n, a, lda, s, scond, amax, info)
    return ccall((@blasfunc(dpoequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}), n, a, lda, s, scond, amax, info)
end

function dpoequb(n, a, lda, s, scond, amax, info)
    return ccall((@blasfunc(dpoequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}), n, a, lda, s, scond, amax, info)
end

function dporfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, iwork,
                info)
    return ccall((@blasfunc(dporfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, iwork, info, 1)
end

function dporfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr,
                 n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork,
                 info)
    return ccall((@blasfunc(dporfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, equed, n, nrhs, a, lda, af,
                 ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, iwork, info, 1, 1)
end

function dposv(uplo, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(dposv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, b, ldb, info, 1)
end

function dposvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                ferr, berr, work, iwork, info)
    return ccall((@blasfunc(dposvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                 ferr, berr, work, iwork, info, 1, 1, 1)
end

function dposvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                 rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params,
                 work, iwork, info)
    return ccall((@blasfunc(dposvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                 rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params,
                 work, iwork, info, 1, 1, 1)
end

function dpotf2(uplo, n, a, lda, info)
    return ccall((@blasfunc(dpotf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function dpotrf(uplo, n, a, lda, info)
    return ccall((@blasfunc(dpotrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function dpotrf2(uplo, n, a, lda, info)
    return ccall((@blasfunc(dpotrf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function dpotri(uplo, n, a, lda, info)
    return ccall((@blasfunc(dpotri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(dpotrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, b, ldb, info, 1)
end

function dppcon(uplo, n, ap, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dppcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, ap, anorm, rcond,
                 work, iwork, info, 1)
end

function dppequ(uplo, n, ap, s, scond, amax, info)
    return ccall((@blasfunc(dppequ_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Clong), uplo, n, ap, s, scond, amax, info, 1)
end

function dpprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(dpprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 ap, afp, b, ldb, x, ldx, ferr, berr, work, iwork, info, 1)
end

function dppsv(uplo, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(dppsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, nrhs, ap, b, ldb, info, 1)
end

function dppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr,
                work, iwork, info)
    return ccall((@blasfunc(dppsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{UInt8}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), fact, uplo, n, nrhs, ap, afp,
                 equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, 1, 1, 1)
end

function dpptrf(uplo, n, ap, info)
    return ccall((@blasfunc(dpptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, ap, info,
                 1)
end

function dpptri(uplo, n, ap, info)
    return ccall((@blasfunc(dpptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, ap, info,
                 1)
end

function dpptrs(uplo, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(dpptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, nrhs, ap, b, ldb, info, 1)
end

function dpstf2(uplo, n, a, lda, piv, rank, tol, work, info)
    return ccall((@blasfunc(dpstf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, a, lda, piv, rank,
                 tol, work, info, 1)
end

function dpstrf(uplo, n, a, lda, piv, rank, tol, work, info)
    return ccall((@blasfunc(dpstrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, a, lda, piv, rank,
                 tol, work, info, 1)
end

function dptcon(n, d, e, anorm, rcond, work, info)
    return ccall((@blasfunc(dptcon_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}), n, d, e, anorm, rcond, work, info)
end

function dpteqr(compz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(dpteqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work,
                 info, 1)
end

function dptrfs(n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, info)
    return ccall((@blasfunc(dptrfs_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), n, nrhs, d, e, df,
                 ef, b, ldb, x, ldx, ferr, berr, work, info)
end

function dptsv(n, nrhs, d, e, b, ldb, info)
    return ccall((@blasfunc(dptsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), n, nrhs, d, e, b, ldb, info)
end

function dptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, info)
    return ccall((@blasfunc(dptsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Clong), fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond,
                 ferr, berr, work, info, 1)
end

function dpttrf(n, d, e, info)
    return ccall((@blasfunc(dpttrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), n, d, e, info)
end

function dpttrs(n, nrhs, d, e, b, ldb, info)
    return ccall((@blasfunc(dpttrs_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), n, nrhs, d, e, b, ldb, info)
end

function dptts2(n, nrhs, d, e, b, ldb)
    return ccall((@blasfunc(dptts2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}), n, nrhs, d, e, b, ldb)
end

function drscl(n, sa, sx, incx)
    return ccall((@blasfunc(drscl_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}), n, sa, sx, incx)
end

function dsb2st_kernels(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt,
                        work)
    return ccall((@blasfunc(dsb2st_kernels_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Clong), uplo, wantz, ttype, st, ed,
                 sweep, n, nb, ib, a, lda, v, tau, ldvt, work, 1)
end

function dsbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, info)
    return ccall((@blasfunc(dsbev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong), jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, info, 1, 1)
end

function dsbev_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, info)
    return ccall((@blasfunc(dsbev_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, info,
                 1, 1)
end

function dsbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dsbevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z,
                 ldz, work, lwork, iwork, liwork, info, 1, 1)
end

function dsbevd_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork,
                       info)
    return ccall((@blasfunc(dsbevd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z,
                 ldz, work, lwork, iwork, liwork, info, 1, 1)
end

function dsbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z,
                ldz, work, iwork, ifail, info)
    return ccall((@blasfunc(dsbevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz,
                 range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz,
                 work, iwork, ifail, info, 1, 1, 1)
end

function dsbevx_2stage(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol,
                       m, w, z, ldz, work, lwork, iwork, ifail, info)
    return ccall((@blasfunc(dsbevx_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu,
                 abstol, m, w, z, ldz, work, lwork, iwork, ifail, info, 1, 1, 1)
end

function dsbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, info)
    return ccall((@blasfunc(dsbgst_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong), vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x,
                 ldx, work, info, 1, 1)
end

function dsbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, info)
    return ccall((@blasfunc(dsbgv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, ka, kb, ab, ldab,
                 bb, ldbb, w, z, ldz, work, info, 1, 1)
end

function dsbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, iwork,
                liwork, info)
    return ccall((@blasfunc(dsbgvd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, iwork,
                 liwork, info, 1, 1)
end

function dsbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu,
                abstol, m, w, z, ldz, work, iwork, ifail, info)
    return ccall((@blasfunc(dsbgvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo, n, ka, kb, ab, ldab,
                 bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail,
                 info, 1, 1, 1)
end

function dsbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info)
    return ccall((@blasfunc(dsbtrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong), vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work,
                 info, 1, 1)
end

function dsfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c)
    return ccall((@blasfunc(dsfrk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Clong, Clong, Clong),
                 transr, uplo, trans, n, k, alpha, a, lda, beta, c, 1, 1, 1)
end

function dsgesv(n, nrhs, a, lda, ipiv, b, ldb, x, ldx, work, swork, iter, info)
    return ccall((@blasfunc(dsgesv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), n, nrhs, a, lda, ipiv, b, ldb, x, ldx, work, swork, iter,
                 info)
end

function dspcon(uplo, n, ap, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dspcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, ap,
                 ipiv, anorm, rcond, work, iwork, info, 1)
end

function dspev(jobz, uplo, n, ap, w, z, ldz, work, info)
    return ccall((@blasfunc(dspev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), jobz,
                 uplo, n, ap, w, z, ldz, work, info, 1, 1)
end

function dspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dspevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), jobz, uplo, n, ap, w, z, ldz, work, lwork,
                 iwork, liwork, info, 1, 1)
end

function dspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork,
                ifail, info)
    return ccall((@blasfunc(dspevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m,
                 w, z, ldz, work, iwork, ifail, info, 1, 1, 1)
end

function dspgst(itype, uplo, n, ap, bp, info)
    return ccall((@blasfunc(dspgst_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Clong), itype, uplo, n, ap, bp, info, 1)
end

function dspgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info)
    return ccall((@blasfunc(dspgv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong), itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info, 1, 1)
end

function dspgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dspgvd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, ap, bp, w, z,
                 ldz, work, lwork, iwork, liwork, info, 1, 1)
end

function dspgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz,
                work, iwork, ifail, info)
    return ccall((@blasfunc(dspgvx_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz,
                 work, iwork, ifail, info, 1, 1, 1)
end

function dsposv(uplo, n, nrhs, a, lda, b, ldb, x, ldx, work, swork, iter, info)
    return ccall((@blasfunc(dsposv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, b, ldb, x, ldx, work, swork,
                 iter, info, 1)
end

function dsprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(dsprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, 1)
end

function dspsv(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(dspsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b, ldb, info, 1)
end

function dspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(dspsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr,
                 berr, work, iwork, info, 1, 1)
end

function dsptrd(uplo, n, ap, d, e, tau, info)
    return ccall((@blasfunc(dsptrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, ap, d, e, tau, info, 1)
end

function dsptrf(uplo, n, ap, ipiv, info)
    return ccall((@blasfunc(dsptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, ap, ipiv, info, 1)
end

function dsptri(uplo, n, ap, ipiv, work, info)
    return ccall((@blasfunc(dsptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Clong), uplo, n, ap, ipiv, work, info, 1)
end

function dsptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(dsptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b, ldb, info, 1)
end

function dstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit,
                work, iwork, info)
    return ccall((@blasfunc(dstebz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong), range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit,
                 w, iblock, isplit, work, iwork, info, 1, 1)
end

function dstedc(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dstedc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info, 1)
end

function dstegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
                lwork, iwork, liwork, info)
    return ccall((@blasfunc(dstegr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n, d, e, vl,
                 vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info,
                 1, 1)
end

function dstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info)
    return ccall((@blasfunc(dstein_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail,
                 info)
end

function dstemr(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac,
                work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dstemr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n,
                 d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork,
                 iwork, liwork, info, 1, 1)
end

function dsteqr(compz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(dsteqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work,
                 info, 1)
end

function dsterf(n, d, e, info)
    return ccall((@blasfunc(dsterf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}), n, d, e, info)
end

function dstev(jobz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(dstev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), jobz, n, d, e, z, ldz, work,
                 info, 1)
end

function dstevd(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dstevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info, 1)
end

function dstevr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
                lwork, iwork, liwork, info)
    return ccall((@blasfunc(dstevr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n, d, e, vl,
                 vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info,
                 1, 1)
end

function dstevx(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork,
                ifail, info)
    return ccall((@blasfunc(dstevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong), jobz, range, n, d, e, vl, vu, il, iu, abstol, m,
                 w, z, ldz, work, iwork, ifail, info, 1, 1)
end

function dsycon(uplo, n, a, lda, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dsycon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 ipiv, anorm, rcond, work, iwork, info, 1)
end

function dsycon_3(uplo, n, a, lda, e, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dsycon_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, e, ipiv, anorm, rcond, work, iwork, info, 1)
end

function dsycon_rook(uplo, n, a, lda, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(dsycon_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 ipiv, anorm, rcond, work, iwork, info, 1)
end

function dsyconv(uplo, way, n, a, lda, ipiv, e, info)
    return ccall((@blasfunc(dsyconv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a, lda, ipiv, e,
                 info, 1, 1)
end

function dsyconvf(uplo, way, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(dsyconvf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a, lda, e, ipiv, info,
                 1, 1)
end

function dsyconvf_rook(uplo, way, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(dsyconvf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a, lda, e, ipiv, info,
                 1, 1)
end

function dsyequb(uplo, n, a, lda, s, scond, amax, work, info)
    return ccall((@blasfunc(dsyequb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, a,
                 lda, s, scond, amax, work, info, 1)
end

function dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
    return ccall((@blasfunc(dsyev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, a, lda,
                 w, work, lwork, info, 1, 1)
end

function dsyev_2stage(jobz, uplo, n, a, lda, w, work, lwork, info)
    return ccall((@blasfunc(dsyev_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, a, lda,
                 w, work, lwork, info, 1, 1)
end

function dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dsyevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info, 1, 1)
end

function dsyevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dsyevd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info, 1, 1)
end

function dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz,
                work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dsyevr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo,
                 n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,
                 iwork, liwork, info, 1, 1, 1)
end

function dsyevr_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                       isuppz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dsyevr_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo,
                 n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,
                 iwork, liwork, info, 1, 1, 1)
end

function dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work,
                lwork, iwork, ifail, info)
    return ccall((@blasfunc(dsyevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo, n, a, lda,
                 vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info, 1,
                 1, 1)
end

function dsyevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                       work, lwork, iwork, ifail, info)
    return ccall((@blasfunc(dsyevx_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo, n, a, lda,
                 vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info, 1,
                 1, 1)
end

function dsygs2(itype, uplo, n, a, lda, b, ldb, info)
    return ccall((@blasfunc(dsygs2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), itype, uplo, n, a, lda, b, ldb, info, 1)
end

function dsygst(itype, uplo, n, a, lda, b, ldb, info)
    return ccall((@blasfunc(dsygst_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), itype, uplo, n, a, lda, b, ldb, info, 1)
end

function dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
    return ccall((@blasfunc(dsygv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info,
                 1, 1)
end

function dsygv_2stage(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
    return ccall((@blasfunc(dsygv_2stage_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info,
                 1, 1)
end

function dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dsygvd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, a, lda, b, ldb,
                 w, work, lwork, iwork, liwork, info, 1, 1)
end

function dsygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w,
                z, ldz, work, lwork, iwork, ifail, info)
    return ccall((@blasfunc(dsygvx_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu,
                 abstol, m, w, z, ldz, work, lwork, iwork, ifail, info, 1, 1, 1)
end

function dsyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(dsyrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1)
end

function dsyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 iwork, info)
    return ccall((@blasfunc(dsyrfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, equed, n, nrhs, a, lda, af,
                 ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, iwork, info, 1, 1)
end

function dsysv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(dsysv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function dsysv_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(dsysv_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function dsysv_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork,
                         info)
    return ccall((@blasfunc(dsysv_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2,
                 b, ldb, work, lwork, info, 1)
end

function dsysv_rk(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(dsysv_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info, 1)
end

function dsysv_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(dsysv_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function dsysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr,
                berr, work, lwork, iwork, info)
    return ccall((@blasfunc(dsysvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), fact, uplo, n, nrhs, a,
                 lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, lwork, iwork,
                 info, 1, 1)
end

function dsysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, iwork, info)
    return ccall((@blasfunc(dsysvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x,
                 ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, iwork, info, 1, 1, 1)
end

function dsyswapr(uplo, n, a, lda, i1, i2)
    return ccall((@blasfunc(dsyswapr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, i1, i2, 1)
end

function dsytd2(uplo, n, a, lda, d, e, tau, info)
    return ccall((@blasfunc(dsytd2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, a, lda, d, e, tau,
                 info, 1)
end

function dsytf2(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(dsytf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function dsytf2_rk(uplo, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(dsytf2_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv, info, 1)
end

function dsytf2_rook(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(dsytf2_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function dsytrd(uplo, n, a, lda, d, e, tau, work, lwork, info)
    return ccall((@blasfunc(dsytrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, d, e, tau, work, lwork, info, 1)
end

function dsytrd_2stage(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info)
    return ccall((@blasfunc(dsytrd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), vect, uplo, n, a, lda, d, e, tau,
                 hous2, lhous2, work, lwork, info, 1, 1)
end

function dsytrd_sb2st(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork,
                      info)
    return ccall((@blasfunc(dsytrd_sb2st_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), stage1, vect,
                 uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info, 1, 1, 1)
end

function dsytrd_sy2sb(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info)
    return ccall((@blasfunc(dsytrd_sy2sb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, kd, a, lda, ab, ldab, tau, work, lwork, info, 1)
end

function dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(dsytrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, lwork, info, 1)
end

function dsytrf_aa(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(dsytrf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, lwork, info, 1)
end

function dsytrf_aa_2stage(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info)
    return ccall((@blasfunc(dsytrf_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info, 1)
end

function dsytrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(dsytrf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv,
                 work, lwork, info, 1)
end

function dsytrf_rook(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(dsytrf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, lwork, info, 1)
end

function dsytri(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(dsytri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function dsytri2(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(dsytri2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, lwork, info, 1)
end

function dsytri2x(uplo, n, a, lda, ipiv, work, nb, info)
    return ccall((@blasfunc(dsytri2x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, nb, info, 1)
end

function dsytri_3(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(dsytri_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv,
                 work, lwork, info, 1)
end

function dsytri_3x(uplo, n, a, lda, e, ipiv, work, nb, info)
    return ccall((@blasfunc(dsytri_3x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv,
                 work, nb, info, 1)
end

function dsytri_rook(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(dsytri_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(dsytrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, ipiv,
                 b, ldb, info, 1)
end

function dsytrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info)
    return ccall((@blasfunc(dsytrs2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 a, lda, ipiv, b, ldb, work, info, 1)
end

function dsytrs_3(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info)
    return ccall((@blasfunc(dsytrs_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a,
                 lda, e, ipiv, b, ldb, info, 1)
end

function dsytrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(dsytrs_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function dsytrs_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info)
    return ccall((@blasfunc(dsytrs_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info, 1)
end

function dsytrs_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(dsytrs_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, ipiv,
                 b, ldb, info, 1)
end

function dtbcon(norm, uplo, diag, n, kd, ab, ldab, rcond, work, iwork, info)
    return ccall((@blasfunc(dtbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), norm, uplo, diag, n, kd, ab, ldab, rcond, work, iwork, info, 1, 1,
                 1)
end

function dtbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(dtbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx,
                 ferr, berr, work, iwork, info, 1, 1, 1)
end

function dtbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(dtbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info, 1, 1, 1)
end

function dtfsm(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb)
    return ccall((@blasfunc(dtfsm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong, Clong), transr, side, uplo, trans, diag, m, n, alpha,
                 a, b, ldb, 1, 1, 1, 1, 1)
end

function dtftri(transr, uplo, diag, n, a, info)
    return ccall((@blasfunc(dtftri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong, Clong), transr, uplo, diag, n, a, info, 1, 1, 1)
end

function dtfttp(transr, uplo, n, arf, ap, info)
    return ccall((@blasfunc(dtfttp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), transr, uplo, n, arf, ap, info, 1, 1)
end

function dtfttr(transr, uplo, n, arf, a, lda, info)
    return ccall((@blasfunc(dtfttr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), transr, uplo, n, arf, a, lda, info, 1, 1)
end

function dtgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, m, work,
                info)
    return ccall((@blasfunc(dtgevc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), side,
                 howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, m, work, info,
                 1, 1)
end

function dtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork,
                info)
    return ccall((@blasfunc(dtgex2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), wantq, wantz,
                 n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork, info)
end

function dtgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork,
                info)
    return ccall((@blasfunc(dtgexc_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), wantq, wantz, n, a, lda,
                 b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork, info)
end

function dtgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alphar, alphai, beta, q, ldq,
                z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(dtgsen_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}), ijob, wantq, wantz, select, n, a, lda,
                 b, ldb, alphar, alphai, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork,
                 iwork, liwork, info)
end

function dtgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u,
                ldu, v, ldv, q, ldq, work, ncycle, info)
    return ccall((@blasfunc(dtgsja_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobu, jobv, jobq, m, p, n, k,
                 l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work,
                 ncycle, info, 1, 1, 1)
end

function dtgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m,
                work, lwork, iwork, info)
    return ccall((@blasfunc(dtgsna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), job, howmny, select, n, a, lda, b,
                 ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info, 1, 1)
end

function dtgsy2(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                rdsum, rdscal, iwork, pq, info)
    return ccall((@blasfunc(dtgsy2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                 rdsum, rdscal, iwork, pq, info, 1)
end

function dtgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                dif, work, lwork, iwork, info)
    return ccall((@blasfunc(dtgsyl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                 dif, work, lwork, iwork, info, 1)
end

function dtpcon(norm, uplo, diag, n, ap, rcond, work, iwork, info)
    return ccall((@blasfunc(dtpcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), norm, uplo,
                 diag, n, ap, rcond, work, iwork, info, 1, 1, 1)
end

function dtplqt(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info)
    return ccall((@blasfunc(dtplqt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}), m, n, l, mb, a, lda, b, ldb, t, ldt, work, info)
end

function dtplqt2(m, n, l, a, lda, b, ldb, t, ldt, info)
    return ccall((@blasfunc(dtplqt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, l, a, lda, b, ldb,
                 t, ldt, info)
end

function dtpmlqt(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info)
    return ccall((@blasfunc(dtpmlqt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong), side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work,
                 info, 1, 1)
end

function dtpmqrt(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info)
    return ccall((@blasfunc(dtpmqrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong), side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work,
                 info, 1, 1)
end

function dtpqrt(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info)
    return ccall((@blasfunc(dtpqrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}), m, n, l, nb, a, lda, b, ldb, t, ldt, work, info)
end

function dtpqrt2(m, n, l, a, lda, b, ldb, t, ldt, info)
    return ccall((@blasfunc(dtpqrt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), m, n, l, a, lda, b, ldb,
                 t, ldt, info)
end

function dtprfb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb,
                work, ldwork)
    return ccall((@blasfunc(dtprfb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), side, trans, direct, storev, m, n, k, l, v,
                 ldv, t, ldt, a, lda, b, ldb, work, ldwork, 1, 1, 1, 1)
end

function dtprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work, iwork,
                info)
    return ccall((@blasfunc(dtprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work, iwork,
                 info, 1, 1, 1)
end

function dtptri(uplo, diag, n, ap, info)
    return ccall((@blasfunc(dtptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong),
                 uplo, diag, n, ap, info, 1, 1)
end

function dtptrs(uplo, trans, diag, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(dtptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, nrhs, ap, b, ldb, info, 1, 1, 1)
end

function dtpttf(transr, uplo, n, ap, arf, info)
    return ccall((@blasfunc(dtpttf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), transr, uplo, n, ap, arf, info, 1, 1)
end

function dtpttr(uplo, n, ap, a, lda, info)
    return ccall((@blasfunc(dtpttr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, ap, a, lda, info, 1)
end

function dtrcon(norm, uplo, diag, n, a, lda, rcond, work, iwork, info)
    return ccall((@blasfunc(dtrcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 norm, uplo, diag, n, a, lda, rcond, work, iwork, info, 1, 1, 1)
end

function dtrevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, info)
    return ccall((@blasfunc(dtrevc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), side, howmny, select, n, t, ldt,
                 vl, ldvl, vr, ldvr, mm, m, work, info, 1, 1)
end

function dtrevc3(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork,
                 info)
    return ccall((@blasfunc(dtrevc3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, howmny, select,
                 n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info, 1, 1)
end

function dtrexc(compq, n, t, ldt, q, ldq, ifst, ilst, work, info)
    return ccall((@blasfunc(dtrexc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), compq, n, t, ldt,
                 q, ldq, ifst, ilst, work, info, 1)
end

function dtrrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr, work, iwork,
                info)
    return ccall((@blasfunc(dtrrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1, 1, 1)
end

function dtrsen(job, compq, select, n, t, ldt, q, ldq, wr, wi, m, s, sep, work, lwork,
                iwork, liwork, info)
    return ccall((@blasfunc(dtrsen_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), job, compq, select, n, t, ldt, q, ldq, wr, wi,
                 m, s, sep, work, lwork, iwork, liwork, info, 1, 1)
end

function dtrsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, m, work,
                ldwork, iwork, info)
    return ccall((@blasfunc(dtrsna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong), job, howmny, select, n, t, ldt, vl, ldvl, vr,
                 ldvr, s, sep, mm, m, work, ldwork, iwork, info, 1, 1)
end

function dtrsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info)
    return ccall((@blasfunc(dtrsyl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{BlasInt}, Clong, Clong), trana, tranb, isgn, m, n, a, lda, b, ldb, c,
                 ldc, scale, info, 1, 1)
end

function dtrsyl3(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, iwork, liwork,
                 swork, ldswork, info)
    return ccall((@blasfunc(dtrsyl3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, iwork, liwork,
                 swork, ldswork, info, 1, 1)
end

function dtrti2(uplo, diag, n, a, lda, info)
    return ccall((@blasfunc(dtrti2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), uplo, diag, n, a, lda, info, 1, 1)
end

function dtrtri(uplo, diag, n, a, lda, info)
    return ccall((@blasfunc(dtrtri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), uplo, diag, n, a, lda, info, 1, 1)
end

function dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(dtrtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), uplo,
                 trans, diag, n, nrhs, a, lda, b, ldb, info, 1, 1, 1)
end

function dtrttf(transr, uplo, n, a, lda, arf, info)
    return ccall((@blasfunc(dtrttf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong), transr, uplo, n, a, lda, arf, info, 1, 1)
end

function dtrttp(uplo, n, a, lda, ap, info)
    return ccall((@blasfunc(dtrttp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ap, info, 1)
end

function dtzrzf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(dtzrzf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function dzsum1(n, cx, incx)
    return ccall((@blasfunc(dzsum1_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, cx, incx)
end

function icmax1(n, cx, incx)
    return ccall((@blasfunc(icmax1_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, cx, incx)
end

function ieeeck(ispec, zero, one)
    return ccall((@blasfunc(ieeeck_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ref{Float32}, Ref{Float32}), ispec, zero, one)
end

function ilaclc(m, n, a, lda)
    return ccall((@blasfunc(ilaclc_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, a, lda)
end

function ilaclr(m, n, a, lda)
    return ccall((@blasfunc(ilaclr_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), m, n, a, lda)
end

function iladiag(diag)
    return ccall((@blasfunc(iladiag_), libblastrampoline), BlasInt, (Ref{UInt8}, Clong), diag, 1)
end

function iladlc(m, n, a, lda)
    return ccall((@blasfunc(iladlc_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), m, n, a, lda)
end

function iladlr(m, n, a, lda)
    return ccall((@blasfunc(iladlr_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), m, n, a, lda)
end

function ilaprec(prec)
    return ccall((@blasfunc(ilaprec_), libblastrampoline), BlasInt, (Ref{UInt8}, Clong), prec, 1)
end

function ilaslc(m, n, a, lda)
    return ccall((@blasfunc(ilaslc_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), m, n, a, lda)
end

function ilaslr(m, n, a, lda)
    return ccall((@blasfunc(ilaslr_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), m, n, a, lda)
end

function ilatrans(trans)
    return ccall((@blasfunc(ilatrans_), libblastrampoline), BlasInt, (Ref{UInt8}, Clong), trans, 1)
end

function ilauplo(uplo)
    return ccall((@blasfunc(ilauplo_), libblastrampoline), BlasInt, (Ref{UInt8}, Clong), uplo, 1)
end

function ilazlc(m, n, a, lda)
    return ccall((@blasfunc(ilazlc_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, a, lda)
end

function ilazlr(m, n, a, lda)
    return ccall((@blasfunc(ilazlr_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, a, lda)
end

function izmax1(n, zx, incx)
    return ccall((@blasfunc(izmax1_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, zx, incx)
end

function sbbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2,
                ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d,
                b22e, work, lwork, info)
    return ccall((@blasfunc(sbbcsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong), jobu1, jobu2, jobv1t,
                 jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t,
                 ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d, b22e, work, lwork, info,
                 1, 1, 1, 1, 1)
end

function sbdsdc(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info)
    return ccall((@blasfunc(sbdsdc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, compq, n, d, e,
                 u, ldu, vt, ldvt, q, iq, work, iwork, info, 1, 1)
end

function sbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info)
    return ccall((@blasfunc(sbdsqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, ncvt,
                 nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info, 1)
end

function sbdsvdx(uplo, jobz, range, n, d, e, vl, vu, il, iu, ns, s, z, ldz, work, iwork,
                 info)
    return ccall((@blasfunc(sbdsvdx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), uplo, jobz, range, n, d, e, vl, vu, il, iu, ns, s, z, ldz, work,
                 iwork, info, 1, 1, 1)
end

function scsum1(n, cx, incx)
    return ccall((@blasfunc(scsum1_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ptr{ComplexF32}, Ref{BlasInt}), n, cx, incx)
end

function sdisna(job, m, n, d, sep, info)
    return ccall((@blasfunc(sdisna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Clong), job, m, n, d, sep, info, 1)
end

function sgbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work,
                info)
    return ccall((@blasfunc(sgbbrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Clong), vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt,
                 ldpt, c, ldc, work, info, 1)
end

function sgbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(sgbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, iwork, info,
                 1)
end

function sgbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(sgbequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{BlasInt}), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
end

function sgbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(sgbequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{BlasInt}), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
end

function sgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr,
                berr, work, iwork, info)
    return ccall((@blasfunc(sgbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), trans, n, kl, ku, nrhs, ab, ldab, afb,
                 ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, 1)
end

function sgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x,
                 ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, iwork, info)
    return ccall((@blasfunc(sgbrfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c,
                 b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, iwork, info, 1, 1)
end

function sgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    return ccall((@blasfunc(sgbsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), n, kl, ku, nrhs, ab, ldab,
                 ipiv, b, ldb, info)
end

function sgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
                ldb, x, ldx, rcond, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(sgbsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, kl, ku, nrhs,
                 ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr,
                 work, iwork, info, 1, 1, 1)
end

function sgbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
                 ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, iwork, info)
    return ccall((@blasfunc(sgbsvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, kl, ku, nrhs,
                 ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 iwork, info, 1, 1, 1)
end

function sgbtf2(m, n, kl, ku, ab, ldab, ipiv, info)
    return ccall((@blasfunc(sgbtf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, kl, ku, ab, ldab, ipiv, info)
end

function sgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
    return ccall((@blasfunc(sgbtrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, kl, ku, ab, ldab, ipiv, info)
end

function sgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    return ccall((@blasfunc(sgbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), trans,
                 n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info, 1)
end

function sgebak(job, side, n, ilo, ihi, scale, m, v, ldv, info)
    return ccall((@blasfunc(sgebak_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job, side,
                 n, ilo, ihi, scale, m, v, ldv, info, 1, 1)
end

function sgebal(job, n, a, lda, ilo, ihi, scale, info)
    return ccall((@blasfunc(sgebal_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong), job, n, a, lda, ilo, ihi, scale, info, 1)
end

function sgebd2(m, n, a, lda, d, e, tauq, taup, work, info)
    return ccall((@blasfunc(sgebd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), m, n, a, lda, d, e,
                 tauq, taup, work, info)
end

function sgebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info)
    return ccall((@blasfunc(sgebrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a,
                 lda, d, e, tauq, taup, work, lwork, info)
end

function sgecon(norm, n, a, lda, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(sgecon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), norm, n, a, lda,
                 anorm, rcond, work, iwork, info, 1)
end

function sgeequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(sgeequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}), m, n, a, lda, r, c,
                 rowcnd, colcnd, amax, info)
end

function sgeequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(sgeequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}), m, n, a, lda, r, c,
                 rowcnd, colcnd, amax, info)
end

function sgees(jobvs, sort, select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork,
               info)
    return ccall((@blasfunc(sgees_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), jobvs, sort,
                 select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info, 1, 1)
end

function sgeesx(jobvs, sort, select, sense, n, a, lda, sdim, wr, wi, vs, ldvs, rconde,
                rcondv, work, lwork, iwork, liwork, bwork, info)
    return ccall((@blasfunc(sgeesx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobvs, sort, select, sense, n,
                 a, lda, sdim, wr, wi, vs, ldvs, rconde, rcondv, work, lwork, iwork, liwork,
                 bwork, info, 1, 1, 1)
end

function sgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    return ccall((@blasfunc(sgeev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobvl, jobvr, n, a,
                 lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info, 1, 1)
end

function sgeevx(balanc, jobvl, jobvr, sense, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, ilo,
                ihi, scale, abnrm, rconde, rcondv, work, lwork, iwork, info)
    return ccall((@blasfunc(sgeevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), balanc, jobvl, jobvr, sense, n, a, lda, wr,
                 wi, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm, rconde, rcondv, work,
                 lwork, iwork, info, 1, 1, 1, 1)
end

function sgehd2(n, ilo, ihi, a, lda, tau, work, info)
    return ccall((@blasfunc(sgehd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}), n, ilo, ihi, a, lda, tau, work, info)
end

function sgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sgehrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), n, ilo, ihi, a, lda, tau, work,
                 lwork, info)
end

function sgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work,
                lwork, iwork, info)
    return ccall((@blasfunc(sgejsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong, Clong), joba, jobu, jobv,
                 jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv, work, lwork, iwork,
                 info, 1, 1, 1, 1, 1, 1)
end

function sgelq(m, n, a, lda, t, tsize, work, lwork, info)
    return ccall((@blasfunc(sgelq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, tsize, work, lwork,
                 info)
end

function sgelq2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(sgelq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function sgelqf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sgelqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function sgelqt(m, n, mb, a, lda, t, ldt, work, info)
    return ccall((@blasfunc(sgelqt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), m, n, mb, a, lda, t, ldt, work, info)
end

function sgelqt3(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(sgelqt3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function sgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(sgels_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function sgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
    return ccall((@blasfunc(sgelsd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work,
                 lwork, iwork, info)
end

function sgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
    return ccall((@blasfunc(sgelss_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
end

function sgelst(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(sgelst_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function sgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
    return ccall((@blasfunc(sgelsy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork,
                 info)
end

function sgemlq(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sgemlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, t,
                 tsize, c, ldc, work, lwork, info, 1, 1)
end

function sgemlqt(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info)
    return ccall((@blasfunc(sgemlqt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, mb, v, ldv,
                 t, ldt, c, ldc, work, info, 1, 1)
end

function sgemqr(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sgemqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, t,
                 tsize, c, ldc, work, lwork, info, 1, 1)
end

function sgemqrt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info)
    return ccall((@blasfunc(sgemqrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, nb, v, ldv,
                 t, ldt, c, ldc, work, info, 1, 1)
end

function sgeql2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(sgeql2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function sgeqlf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sgeqlf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function sgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info)
    return ccall((@blasfunc(sgeqp3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, jpvt, tau, work, lwork,
                 info)
end

function sgeqp3rk(m, n, nrhs, kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk,
                  jpiv, tau, work, lwork, iwork, info)
    return ccall((@blasfunc(sgeqp3rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs,
                 kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk, jpiv, tau, work,
                 lwork, iwork, info)
end

function sgeqr(m, n, a, lda, t, tsize, work, lwork, info)
    return ccall((@blasfunc(sgeqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, tsize, work, lwork,
                 info)
end

function sgeqr2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(sgeqr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function sgeqr2p(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(sgeqr2p_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function sgeqrf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sgeqrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function sgeqrfp(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sgeqrfp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function sgeqrt(m, n, nb, a, lda, t, ldt, work, info)
    return ccall((@blasfunc(sgeqrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), m, n, nb, a, lda, t, ldt, work, info)
end

function sgeqrt2(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(sgeqrt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function sgeqrt3(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(sgeqrt3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function sgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(sgerfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1)
end

function sgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 iwork, info)
    return ccall((@blasfunc(sgerfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), trans,
                 equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond, berr,
                 n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork,
                 info, 1, 1)
end

function sgerq2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(sgerq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function sgerqf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sgerqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function sgesc2(n, a, lda, rhs, ipiv, jpiv, scale)
    return ccall((@blasfunc(sgesc2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{Float32}), n, a, lda, rhs, ipiv, jpiv, scale)
end

function sgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
    return ccall((@blasfunc(sgesdd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), jobz, m, n, a, lda, s, u, ldu, vt, ldvt,
                 work, lwork, iwork, info, 1)
end

function sgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(sgesv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), n, nrhs, a, lda, ipiv, b, ldb, info)
end

function sgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    return ccall((@blasfunc(sgesvd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobu, jobvt, m, n, a,
                 lda, s, u, ldu, vt, ldvt, work, lwork, info, 1, 1)
end

function sgesvdq(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank,
                 iwork, liwork, work, lwork, rwork, lrwork, info)
    return ccall((@blasfunc(sgesvdq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong,
                  Clong), joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv,
                 numrank, iwork, liwork, work, lwork, rwork, lrwork, info, 1, 1, 1, 1, 1)
end

function sgesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt,
                 work, lwork, iwork, info)
    return ccall((@blasfunc(sgesvdx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobu,
                 jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt, work,
                 lwork, iwork, info, 1, 1, 1)
end

function sgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info)
    return ccall((@blasfunc(sgesvj_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), joba, jobu, jobv, m, n, a,
                 lda, sva, mv, v, ldv, work, lwork, info, 1, 1, 1)
end

function sgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                rcond, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(sgesvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, nrhs, a, lda, af, ldaf,
                 ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, 1,
                 1, 1)
end

function sgesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, iwork, info)
    return ccall((@blasfunc(sgesvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed,
                 r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, iwork, info, 1, 1, 1)
end

function sgetc2(n, a, lda, ipiv, jpiv, info)
    return ccall((@blasfunc(sgetc2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), n,
                 a, lda, ipiv, jpiv, info)
end

function sgetf2(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(sgetf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m,
                 n, a, lda, ipiv, info)
end

function sgetrf(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(sgetrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m,
                 n, a, lda, ipiv, info)
end

function sgetrf2(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(sgetrf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m,
                 n, a, lda, ipiv, info)
end

function sgetri(n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(sgetri_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), n, a, lda, ipiv, work, lwork, info)
end

function sgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(sgetrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), trans, n, nrhs, a, lda, ipiv,
                 b, ldb, info, 1)
end

function sgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(sgetsls_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function sgetsqrhrt(m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(sgetsqrhrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}),
                 m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info)
end

function sggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info)
    return ccall((@blasfunc(sggbak_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info, 1, 1)
end

function sggbal(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info)
    return ccall((@blasfunc(sggbal_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Clong), job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info, 1)
end

function sgges(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta,
               vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info)
    return ccall((@blasfunc(sgges_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai,
                 beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info, 1, 1, 1)
end

function sgges3(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai, beta,
                vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info)
    return ccall((@blasfunc(sgges3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alphar, alphai,
                 beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, info, 1, 1, 1)
end

function sggesx(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alphar,
                alphai, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork,
                liwork, bwork, info)
    return ccall((@blasfunc(sggesx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong), jobvsl,
                 jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alphar, alphai, beta,
                 vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, iwork, liwork, bwork,
                 info, 1, 1, 1, 1)
end

function sggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr,
               work, lwork, info)
    return ccall((@blasfunc(sggev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                 ldvl, vr, ldvr, work, lwork, info, 1, 1)
end

function sggev3(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr,
                work, lwork, info)
    return ccall((@blasfunc(sggev3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                 ldvl, vr, ldvr, work, lwork, info, 1, 1)
end

function sggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alphar, alphai, beta, vl,
                ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv,
                work, lwork, iwork, bwork, info)
    return ccall((@blasfunc(sggevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong), balanc, jobvl, jobvr, sense, n, a,
                 lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, ilo, ihi, lscale,
                 rscale, abnrm, bbnrm, rconde, rcondv, work, lwork, iwork, bwork, info, 1,
                 1, 1, 1)
end

function sggglm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info)
    return ccall((@blasfunc(sggglm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), n, m, p, a, lda, b, ldb, d, x, y, work, lwork,
                 info)
end

function sgghd3(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork,
                info)
    return ccall((@blasfunc(sgghd3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), compq,
                 compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork, info, 1,
                 1)
end

function sgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info)
    return ccall((@blasfunc(sgghrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), compq, compz, n, ilo, ihi, a, lda, b,
                 ldb, q, ldq, z, ldz, info, 1, 1)
end

function sgglse(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    return ccall((@blasfunc(sgglse_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, p, a, lda, b, ldb, c, d, x, work, lwork,
                 info)
end

function sggqrf(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info)
    return ccall((@blasfunc(sggqrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info)
end

function sggrqf(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info)
    return ccall((@blasfunc(sggrqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info)
end

function sggsvd3(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v,
                 ldv, q, ldq, work, lwork, iwork, info)
    return ccall((@blasfunc(sggsvd3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), jobu, jobv, jobq, m, n, p, k, l, a, lda,
                 b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info, 1,
                 1, 1)
end

function sggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v,
                 ldv, q, ldq, iwork, tau, work, lwork, info)
    return ccall((@blasfunc(sggsvp3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobu, jobv, jobq, m, p, n, a,
                 lda, b, ldb, tola, tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, tau, work,
                 lwork, info, 1, 1, 1)
end

function sgsvj0(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work,
                lwork, info)
    return ccall((@blasfunc(sgsvj0_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep,
                 work, lwork, info, 1)
end

function sgsvj1(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work,
                lwork, info)
    return ccall((@blasfunc(sgsvj1_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), jobv, m, n, n1, a, lda, d, sva, mv, v, ldv,
                 eps, sfmin, tol, nsweep, work, lwork, info, 1)
end

function sgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(sgtcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), norm, n, dl, d, du, du2, ipiv, anorm, rcond,
                 work, iwork, info, 1)
end

function sgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr,
                berr, work, iwork, info)
    return ccall((@blasfunc(sgtrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), trans, n, nrhs,
                 dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1)
end

function sgtsv(n, nrhs, dl, d, du, b, ldb, info)
    return ccall((@blasfunc(sgtsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), n, nrhs, dl, d, du, b, ldb, info)
end

function sgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx,
                rcond, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(sgtsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb,
                 x, ldx, rcond, ferr, berr, work, iwork, info, 1, 1)
end

function sgttrf(n, dl, d, du, du2, ipiv, info)
    return ccall((@blasfunc(sgttrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}), n, dl, d, du, du2, ipiv, info)
end

function sgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info)
    return ccall((@blasfunc(sgttrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info, 1)
end

function sgtts2(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb)
    return ccall((@blasfunc(sgtts2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}), itrans, n, nrhs, dl, d,
                 du, du2, ipiv, b, ldb)
end

function shgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q,
                ldq, z, ldz, work, lwork, info)
    return ccall((@blasfunc(shgeqz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), job,
                 compq, compz, n, ilo, ihi, h, ldh, t, ldt, alphar, alphai, beta, q, ldq, z,
                 ldz, work, lwork, info, 1, 1, 1)
end

function shsein(side, eigsrc, initv, select, n, h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m,
                work, ifaill, ifailr, info)
    return ccall((@blasfunc(shsein_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), side, eigsrc, initv, select,
                 n, h, ldh, wr, wi, vl, ldvl, vr, ldvr, mm, m, work, ifaill, ifailr, info,
                 1, 1, 1)
end

function shseqr(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info)
    return ccall((@blasfunc(shseqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job, compz, n, ilo,
                 ihi, h, ldh, wr, wi, z, ldz, work, lwork, info, 1, 1)
end

function sla_gbamv(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy)
    return ccall((@blasfunc(sla_gbamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}), trans, m, n, kl, ku, alpha, ab, ldab, x, incx,
                 beta, y, incy)
end

function sla_gbrcond(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, cmode, c, info, work,
                     iwork)
    return ccall((@blasfunc(sla_gbrcond_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Clong), trans, n, kl, ku, ab, ldab, afb, ldafb,
                 ipiv, cmode, c, info, work, iwork, 1)
end

function sla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb,
                             ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms,
                             err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond,
                             ithresh, rthresh, dz_ub, ignore_cwise, info)
    return ccall((@blasfunc(sla_gbrfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}), prec_type, trans_type,
                 n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, colequ, c, b, ldb, y, ldy,
                 berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail,
                 rcond, ithresh, rthresh, dz_ub, ignore_cwise, info)
end

function sla_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb)
    return ccall((@blasfunc(sla_gbrpvgrw_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}), n, kl, ku, ncols, ab, ldab, afb, ldafb)
end

function sla_geamv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(sla_geamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}), trans, m,
                 n, alpha, a, lda, x, incx, beta, y, incy)
end

function sla_gercond(trans, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork)
    return ccall((@blasfunc(sla_gercond_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Clong), trans, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork, 1)
end

function sla_gerfsx_extended(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ,
                             c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb,
                             dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info)
    return ccall((@blasfunc(sla_gerfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv,
                 colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb, dy,
                 y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info)
end

function sla_gerpvgrw(n, ncols, a, lda, af, ldaf)
    return ccall((@blasfunc(sla_gerpvgrw_), libblastrampoline), Float32,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}),
                 n, ncols, a, lda, af, ldaf)
end

function sla_lin_berr(n, nz, nrhs, res, ayb, berr)
    return ccall((@blasfunc(sla_lin_berr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}), n, nz, nrhs, res, ayb, berr)
end

function sla_porcond(uplo, n, a, lda, af, ldaf, cmode, c, info, work, iwork)
    return ccall((@blasfunc(sla_porcond_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Clong), uplo,
                 n, a, lda, af, ldaf, cmode, c, info, work, iwork, 1)
end

function sla_porfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb,
                             y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res,
                             ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise,
                             info)
    return ccall((@blasfunc(sla_porfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb, y,
                 ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail,
                 rcond, ithresh, rthresh, dz_ub, ignore_cwise, info, 1)
end

function sla_porpvgrw(uplo, ncols, a, lda, af, ldaf, work)
    return ccall((@blasfunc(sla_porpvgrw_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Clong), uplo, ncols, a, lda, af, ldaf, work, 1)
end

function sla_syamv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(sla_syamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}), uplo, n, alpha, a, lda,
                 x, incx, beta, y, incy)
end

function sla_syrcond(uplo, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork)
    return ccall((@blasfunc(sla_syrcond_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Clong), uplo, n, a, lda, af, ldaf, ipiv, cmode, c, info, work, iwork, 1)
end

function sla_syrfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                             ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp,
                             res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                             ignore_cwise, info)
    return ccall((@blasfunc(sla_syrfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                 ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res, ayb, dy,
                 y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info, 1)
end

function sla_syrpvgrw(uplo, n, info, a, lda, af, ldaf, ipiv, work)
    return ccall((@blasfunc(sla_syrpvgrw_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Clong), uplo, n, info, a, lda, af,
                 ldaf, ipiv, work, 1)
end

function sla_wwaddw(n, x, y, w)
    return ccall((@blasfunc(sla_wwaddw_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}), n, x, y, w)
end

function slabad(small, large)
    return ccall((@blasfunc(slabad_), libblastrampoline), Cvoid, (Ref{Float32}, Ref{Float32}),
                 small, large)
end

function slabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy)
    return ccall((@blasfunc(slabrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}), m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y,
                 ldy)
end

function slacn2(n, v, x, isgn, est, kase, isave)
    return ccall((@blasfunc(slacn2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}), n, v, x, isgn, est, kase, isave)
end

function slacon(n, v, x, isgn, est, kase)
    return ccall((@blasfunc(slacon_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{Float32},
                  Ref{BlasInt}), n, v, x, isgn, est, kase)
end

function slacpy(uplo, m, n, a, lda, b, ldb)
    return ccall((@blasfunc(slacpy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Clong), uplo, m, n, a, lda, b, ldb, 1)
end

function sladiv(a, b, c, d, p, q)
    return ccall((@blasfunc(sladiv_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}), a, b, c, d, p, q)
end

function sladiv1(a, b, c, d, p, q)
    return ccall((@blasfunc(sladiv1_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}), a, b, c, d, p, q)
end

function sladiv2(a, b, c, d, r, t)
    return ccall((@blasfunc(sladiv2_), libblastrampoline), Float32,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}), a, b, c, d, r, t)
end

function slae2(a, b, c, rt1, rt2)
    return ccall((@blasfunc(slae2_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), a,
                 b, c, rt1, rt2)
end

function slaebz(ijob, nitmax, n, mmax, minp, nbmin, abstol, reltol, pivmin, d, e, e2, nval,
                ab, c, mout, nab, work, iwork, info)
    return ccall((@blasfunc(slaebz_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), ijob, nitmax, n, mmax, minp, nbmin,
                 abstol, reltol, pivmin, d, e, e2, nval, ab, c, mout, nab, work, iwork,
                 info)
end

function slaed0(icompq, qsiz, n, d, e, q, ldq, qstore, ldqs, work, iwork, info)
    return ccall((@blasfunc(slaed0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}),
                 icompq, qsiz, n, d, e, q, ldq, qstore, ldqs, work, iwork, info)
end

function slaed1(n, d, q, ldq, indxq, rho, cutpnt, work, iwork, info)
    return ccall((@blasfunc(slaed1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), n, d, q, ldq, indxq, rho,
                 cutpnt, work, iwork, info)
end

function slaed2(k, n, n1, d, q, ldq, indxq, rho, z, dlambda, w, q2, indx, indxc, indxp,
                coltyp, info)
    return ccall((@blasfunc(slaed2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), k,
                 n, n1, d, q, ldq, indxq, rho, z, dlambda, w, q2, indx, indxc, indxp,
                 coltyp, info)
end

function slaed3(k, n, n1, d, q, ldq, rho, dlambda, q2, indx, ctot, w, s, info)
    return ccall((@blasfunc(slaed3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), k, n, n1, d, q, ldq, rho, dlambda,
                 q2, indx, ctot, w, s, info)
end

function slaed4(n, i, d, z, delta, rho, dlam, info)
    return ccall((@blasfunc(slaed4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}), n, i, d, z, delta, rho, dlam,
                 info)
end

function slaed5(i, d, z, delta, rho, dlam)
    return ccall((@blasfunc(slaed5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}), i, d, z, delta, rho, dlam)
end

function slaed6(kniter, orgati, rho, d, z, finit, tau, info)
    return ccall((@blasfunc(slaed6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}), kniter, orgati, rho, d, z, finit,
                 tau, info)
end

function slaed7(icompq, n, qsiz, tlvls, curlvl, curpbm, d, q, ldq, indxq, rho, cutpnt,
                qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, iwork, info)
    return ccall((@blasfunc(slaed7_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), icompq, n, qsiz, tlvls,
                 curlvl, curpbm, d, q, ldq, indxq, rho, cutpnt, qstore, qptr, prmptr, perm,
                 givptr, givcol, givnum, work, iwork, info)
end

function slaed8(icompq, k, n, qsiz, d, q, ldq, indxq, rho, cutpnt, z, dlambda, q2, ldq2, w,
                perm, givptr, givcol, givnum, indxp, indx, info)
    return ccall((@blasfunc(slaed8_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), icompq, k, n, qsiz, d, q,
                 ldq, indxq, rho, cutpnt, z, dlambda, q2, ldq2, w, perm, givptr, givcol,
                 givnum, indxp, indx, info)
end

function slaed9(k, kstart, kstop, n, d, q, ldq, rho, dlambda, w, s, lds, info)
    return ccall((@blasfunc(slaed9_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), k, kstart, kstop, n, d, q, ldq, rho, dlambda, w, s,
                 lds, info)
end

function slaeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, q, qptr, z,
                ztemp, info)
    return ccall((@blasfunc(slaeda_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}), n, tlvls, curlvl, curpbm, prmptr, perm, givptr,
                 givcol, givnum, q, qptr, z, ztemp, info)
end

function slaein(rightv, noinit, n, h, ldh, wr, wi, vr, vi, b, ldb, work, eps3, smlnum,
                bignum, info)
    return ccall((@blasfunc(slaein_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}),
                 rightv, noinit, n, h, ldh, wr, wi, vr, vi, b, ldb, work, eps3, smlnum,
                 bignum, info)
end

function slaev2(a, b, c, rt1, rt2, cs1, sn1)
    return ccall((@blasfunc(slaev2_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}), a, b, c, rt1, rt2, cs1, sn1)
end

function slaexc(wantq, n, t, ldt, q, ldq, j1, n1, n2, work, info)
    return ccall((@blasfunc(slaexc_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), wantq, n, t,
                 ldt, q, ldq, j1, n1, n2, work, info)
end

function slag2(a, lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi)
    return ccall((@blasfunc(slag2_), libblastrampoline), Cvoid,
                 (Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), a,
                 lda, b, ldb, safmin, scale1, scale2, wr1, wr2, wi)
end

function slag2d(m, n, sa, ldsa, a, lda, info)
    return ccall((@blasfunc(slag2d_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, sa, ldsa, a, lda, info)
end

function slags2(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq)
    return ccall((@blasfunc(slags2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}), upper, a1, a2, a3, b1, b2, b3,
                 csu, snu, csv, snv, csq, snq)
end

function slagtf(n, a, lambda, b, c, tol, d, in, info)
    return ccall((@blasfunc(slagtf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), n, a, lambda, b, c,
                 tol, d, in, info)
end

function slagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb)
    return ccall((@blasfunc(slagtm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Clong), trans, n, nrhs, alpha, dl, d, du, x, ldx,
                 beta, b, ldb, 1)
end

function slagts(job, n, a, b, c, d, in, y, tol, info)
    return ccall((@blasfunc(slagts_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32}, Ref{Float32}, Ref{BlasInt}), job, n,
                 a, b, c, d, in, y, tol, info)
end

function slagv2(a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr)
    return ccall((@blasfunc(slagv2_), libblastrampoline), Cvoid,
                 (Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}), a, lda, b, ldb, alphar, alphai, beta, csl, snl, csr, snr)
end

function slahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info)
    return ccall((@blasfunc(slahqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz,
                 ihiz, z, ldz, info)
end

function slahr2(n, k, nb, a, lda, tau, t, ldt, y, ldy)
    return ccall((@blasfunc(slahr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), n, k, nb, a, lda, tau,
                 t, ldt, y, ldy)
end

function slaic1(job, j, x, sest, w, gamma, sestpr, s, c)
    return ccall((@blasfunc(slaic1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), job, j, x, sest,
                 w, gamma, sestpr, s, c)
end

function slaln2(ltrans, na, nw, smin, ca, a, lda, d1, d2, b, ldb, wr, wi, x, ldx, scale,
                xnorm, info)
    return ccall((@blasfunc(slaln2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}), ltrans, na, nw, smin, ca, a, lda, d1, d2, b,
                 ldb, wr, wi, x, ldx, scale, xnorm, info)
end

function slals0(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol,
                givnum, ldgnum, poles, difl, difr, z, k, c, s, work, info)
    return ccall((@blasfunc(slals0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}), icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx,
                 perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c,
                 s, work, info)
end

function slalsa(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z,
                poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork, info)
    return ccall((@blasfunc(slalsa_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), icompq, smlsiz, n, nrhs, b, ldb, bx,
                 ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr, givcol, ldgcol, perm,
                 givnum, c, s, work, iwork, info)
end

function slalsd(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, iwork, info)
    return ccall((@blasfunc(slalsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong), uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work,
                 iwork, info, 1)
end

function slamrg(n1, n2, a, strd1, strd2, index)
    return ccall((@blasfunc(slamrg_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}), n1,
                 n2, a, strd1, strd2, index)
end

function slamswlq(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info)
    return ccall((@blasfunc(slamswlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans,
                 m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info, 1, 1)
end

function slamtsqr(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info)
    return ccall((@blasfunc(slamtsqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans,
                 m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info, 1, 1)
end

function slaneg(n, d, lld, sigma, pivmin, r)
    return ccall((@blasfunc(slaneg_), libblastrampoline), BlasInt,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{BlasInt}), n, d, lld, sigma, pivmin, r)
end

function slangb(norm, n, kl, ku, ab, ldab, work)
    return ccall((@blasfunc(slangb_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Clong), norm, n, kl, ku, ab, ldab, work, 1)
end

function slange(norm, m, n, a, lda, work)
    return ccall((@blasfunc(slange_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Clong), norm, m, n, a, lda, work, 1)
end

function slangt(norm, n, dl, d, du)
    return ccall((@blasfunc(slangt_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Clong),
                 norm, n, dl, d, du, 1)
end

function slanhs(norm, n, a, lda, work)
    return ccall((@blasfunc(slanhs_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong),
                 norm, n, a, lda, work, 1)
end

function slansb(norm, uplo, n, k, ab, ldab, work)
    return ccall((@blasfunc(slansb_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Clong, Clong), norm, uplo, n, k, ab, ldab, work, 1, 1)
end

function slansf(norm, transr, uplo, n, a, work)
    return ccall((@blasfunc(slansf_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Clong, Clong, Clong), norm, transr, uplo, n, a, work, 1, 1, 1)
end

function slansp(norm, uplo, n, ap, work)
    return ccall((@blasfunc(slansp_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Clong,
                  Clong), norm, uplo, n, ap, work, 1, 1)
end

function slanst(norm, n, d, e)
    return ccall((@blasfunc(slanst_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Clong), norm, n, d, e,
                 1)
end

function slansy(norm, uplo, n, a, lda, work)
    return ccall((@blasfunc(slansy_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Clong, Clong), norm, uplo, n, a, lda, work, 1, 1)
end

function slantb(norm, uplo, diag, n, k, ab, ldab, work)
    return ccall((@blasfunc(slantb_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Clong, Clong, Clong), norm, uplo, diag, n, k, ab,
                 ldab, work, 1, 1, 1)
end

function slantp(norm, uplo, diag, n, ap, work)
    return ccall((@blasfunc(slantp_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Clong, Clong, Clong), norm, uplo, diag, n, ap, work, 1, 1, 1)
end

function slantr(norm, uplo, diag, m, n, a, lda, work)
    return ccall((@blasfunc(slantr_), libblastrampoline), Float32,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Clong, Clong, Clong), norm, uplo, diag, m, n, a,
                 lda, work, 1, 1, 1)
end

function slanv2(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn)
    return ccall((@blasfunc(slanv2_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), a,
                 b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn)
end

function slaorhr_col_getrfnp(m, n, a, lda, d, info)
    return ccall((@blasfunc(slaorhr_col_getrfnp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}),
                 m, n, a, lda, d, info)
end

function slaorhr_col_getrfnp2(m, n, a, lda, d, info)
    return ccall((@blasfunc(slaorhr_col_getrfnp2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}),
                 m, n, a, lda, d, info)
end

function slapll(n, x, incx, y, incy, ssmin)
    return ccall((@blasfunc(slapll_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}), n, x, incx, y, incy, ssmin)
end

function slapmr(forwrd, m, n, x, ldx, k)
    return ccall((@blasfunc(slapmr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}),
                 forwrd, m, n, x, ldx, k)
end

function slapmt(forwrd, m, n, x, ldx, k)
    return ccall((@blasfunc(slapmt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}),
                 forwrd, m, n, x, ldx, k)
end

function slapy2(x, y)
    return ccall((@blasfunc(slapy2_), libblastrampoline), Float32, (Ref{Float32}, Ref{Float32}),
                 x, y)
end

function slapy3(x, y, z)
    return ccall((@blasfunc(slapy3_), libblastrampoline), Float32,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}), x, y, z)
end

function slaqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed)
    return ccall((@blasfunc(slaqgb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{UInt8}, Clong), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax,
                 equed, 1)
end

function slaqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed)
    return ccall((@blasfunc(slaqge_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{UInt8}, Clong), m, n, a,
                 lda, r, c, rowcnd, colcnd, amax, equed, 1)
end

function slaqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work)
    return ccall((@blasfunc(slaqp2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}), m, n, offset, a,
                 lda, jpvt, tau, vn1, vn2, work)
end

function slaqp2rk(m, n, nrhs, ioffset, kmax, abstol, reltol, kp1, maxc2nrm, a, lda, k,
                  maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2, work, info)
    return ccall((@blasfunc(slaqp2rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), m, n, nrhs, ioffset, kmax, abstol,
                 reltol, kp1, maxc2nrm, a, lda, k, maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1,
                 vn2, work, info)
end

function slaqp3rk(m, n, nrhs, ioffset, nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb,
                  maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info)
    return ccall((@blasfunc(slaqp3rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs, ioffset, nb, abstol, reltol, kp1,
                 maxc2nrm, a, lda, done, kb, maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2,
                 auxv, f, ldf, iwork, info)
end

function slaqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf)
    return ccall((@blasfunc(slaqps_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), m, n, offset, nb, kb, a, lda,
                 jpvt, tau, vn1, vn2, auxv, f, ldf)
end

function slaqr0(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork,
                info)
    return ccall((@blasfunc(slaqr0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo, ihi,
                 h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info)
end

function slaqr1(n, h, ldh, sr1, si1, sr2, si2, v)
    return ccall((@blasfunc(slaqr1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}), n, h, ldh, sr1, si1, sr2, si2,
                 v)
end

function slaqr2(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si,
                v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork)
    return ccall((@blasfunc(slaqr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}), wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
                 ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work,
                 lwork)
end

function slaqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sr, si,
                v, ldv, nh, t, ldt, nv, wv, ldwv, work, lwork)
    return ccall((@blasfunc(slaqr3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}), wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
                 ihiz, z, ldz, ns, nd, sr, si, v, ldv, nh, t, ldt, nv, wv, ldwv, work,
                 lwork)
end

function slaqr4(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork,
                info)
    return ccall((@blasfunc(slaqr4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo, ihi,
                 h, ldh, wr, wi, iloz, ihiz, z, ldz, work, lwork, info)
end

function slaqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh, iloz, ihiz, z,
                ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh)
    return ccall((@blasfunc(slaqr5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}), wantt, wantz, kacc22, n, ktop, kbot, nshfts, sr, si, h, ldh,
                 iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh)
end

function slaqsb(uplo, n, kd, ab, ldab, s, scond, amax, equed)
    return ccall((@blasfunc(slaqsb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{UInt8}, Clong, Clong), uplo, n, kd, ab,
                 ldab, s, scond, amax, equed, 1, 1)
end

function slaqsp(uplo, n, ap, s, scond, amax, equed)
    return ccall((@blasfunc(slaqsp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{UInt8}, Clong, Clong), uplo, n, ap, s, scond, amax,
                 equed, 1, 1)
end

function slaqsy(uplo, n, a, lda, s, scond, amax, equed)
    return ccall((@blasfunc(slaqsy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{UInt8}, Clong, Clong), uplo, n, a, lda, s,
                 scond, amax, equed, 1, 1)
end

function slaqtr(ltran, lreal, n, t, ldt, b, w, scale, x, work, info)
    return ccall((@blasfunc(slaqtr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), ltran,
                 lreal, n, t, ldt, b, w, scale, x, work, info)
end

function slar1v(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz,
                mingma, r, isuppz, nrminv, resid, rqcorr, work)
    return ccall((@blasfunc(slar1v_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ptr{Float32}), n, b1, bn,
                 lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r,
                 isuppz, nrminv, resid, rqcorr, work)
end

function slar2v(n, x, y, z, incx, c, s, incc)
    return ccall((@blasfunc(slar2v_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), n, x, y, z, incx, c, s, incc)
end

function slarf(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(slarf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong), side, m, n, v, incv, tau,
                 c, ldc, work, 1)
end

function slarf1f(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(slarf1f_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong), side, m, n, v, incv, tau,
                 c, ldc, work, 1)
end

function slarf1l(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(slarf1l_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong), side, m, n, v, incv, tau,
                 c, ldc, work, 1)
end

function slarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork)
    return ccall((@blasfunc(slarfb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong, Clong), side,
                 trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork, 1, 1,
                 1, 1)
end

function slarfb_gett(ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork)
    return ccall((@blasfunc(slarfb_gett_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Clong), ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork, 1)
end

function slarfg(n, alpha, x, incx, tau)
    return ccall((@blasfunc(slarfg_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}), n, alpha,
                 x, incx, tau)
end

function slarfgp(n, alpha, x, incx, tau)
    return ccall((@blasfunc(slarfgp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}), n, alpha,
                 x, incx, tau)
end

function slarft(direct, storev, n, k, v, ldv, tau, t, ldt)
    return ccall((@blasfunc(slarft_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), direct, storev, n,
                 k, v, ldv, tau, t, ldt, 1, 1)
end

function slarfx(side, m, n, v, tau, c, ldc, work)
    return ccall((@blasfunc(slarfx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong), side, m, n, v, tau, c, ldc,
                 work, 1)
end

function slarfy(uplo, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(slarfy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong), uplo, n, v, incv, tau, c,
                 ldc, work, 1)
end

function slargv(n, x, incx, y, incy, c, incc)
    return ccall((@blasfunc(slargv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}), n, x, incx, y, incy, c, incc)
end

function slarmm(anorm, bnorm, cnorm)
    return ccall((@blasfunc(slarmm_), libblastrampoline), Float32,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}), anorm, bnorm, cnorm)
end

function slarnv(idist, iseed, n, x)
    return ccall((@blasfunc(slarnv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}), idist, iseed, n, x)
end

function slarra(n, d, e, e2, spltol, tnrm, nsplit, isplit, info)
    return ccall((@blasfunc(slarra_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), n, d, e, e2, spltol, tnrm,
                 nsplit, isplit, info)
end

function slarrb(n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr, work, iwork,
                pivmin, spdiam, twist, info)
    return ccall((@blasfunc(slarrb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), n, d, lld, ifirst, ilast, rtol1, rtol2, offset, w, wgap, werr,
                 work, iwork, pivmin, spdiam, twist, info)
end

function slarrc(jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info)
    return ccall((@blasfunc(slarrc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), jobt, n, vl, vu, d, e, pivmin, eigcnt, lcnt, rcnt, info, 1)
end

function slarrd(range, order, n, vl, vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit,
                isplit, m, w, werr, wl, wu, iblock, indexw, work, iwork, info)
    return ccall((@blasfunc(slarrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), range, order, n, vl,
                 vu, il, iu, gers, reltol, d, e, e2, pivmin, nsplit, isplit, m, w, werr, wl,
                 wu, iblock, indexw, work, iwork, info, 1, 1)
end

function slarre(range, n, vl, vu, il, iu, d, e, e2, rtol1, rtol2, spltol, nsplit, isplit, m,
                w, werr, wgap, iblock, indexw, gers, pivmin, work, iwork, info)
    return ccall((@blasfunc(slarre_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), range, n, vl, vu, il, iu, d,
                 e, e2, rtol1, rtol2, spltol, nsplit, isplit, m, w, werr, wgap, iblock,
                 indexw, gers, pivmin, work, iwork, info, 1)
end

function slarrf(n, d, l, ld, clstrt, clend, w, wgap, werr, spdiam, clgapl, clgapr, pivmin,
                sigma, dplus, lplus, work, info)
    return ccall((@blasfunc(slarrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}), n, d, l, ld, clstrt, clend, w, wgap, werr,
                 spdiam, clgapl, clgapr, pivmin, sigma, dplus, lplus, work, info)
end

function slarrj(n, d, e2, ifirst, ilast, rtol, offset, w, werr, work, iwork, pivmin, spdiam,
                info)
    return ccall((@blasfunc(slarrj_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}), n, d, e2, ifirst, ilast, rtol,
                 offset, w, werr, work, iwork, pivmin, spdiam, info)
end

function slarrk(n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info)
    return ccall((@blasfunc(slarrk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{BlasInt}), n, iw, gl, gu, d, e2, pivmin, reltol, w, werr, info)
end

function slarrr(n, d, e, info)
    return ccall((@blasfunc(slarrr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), n, d, e, info)
end

function slarrv(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr,
                wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info)
    return ccall((@blasfunc(slarrv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), n, vl, vu, d, l, pivmin, isplit, m,
                 dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z,
                 ldz, isuppz, work, iwork, info)
end

function slarscl2(m, n, d, x, ldx)
    return ccall((@blasfunc(slarscl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), m, n, d, x,
                 ldx)
end

function slartgp(f, g, cs, sn, r)
    return ccall((@blasfunc(slartgp_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), f,
                 g, cs, sn, r)
end

function slartgs(x, y, sigma, cs, sn)
    return ccall((@blasfunc(slartgs_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), x,
                 y, sigma, cs, sn)
end

function slartv(n, x, incx, y, incy, c, s, incc)
    return ccall((@blasfunc(slartv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}), n, x, incx, y, incy, c, s, incc)
end

function slaruv(iseed, n, x)
    return ccall((@blasfunc(slaruv_), libblastrampoline), Cvoid,
                 (Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}), iseed, n, x)
end

function slarz(side, m, n, l, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(slarz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong), side, m, n,
                 l, v, incv, tau, c, ldc, work, 1)
end

function slarzb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work,
                ldwork)
    return ccall((@blasfunc(slarzb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc,
                 work, ldwork, 1, 1, 1, 1)
end

function slarzt(direct, storev, n, k, v, ldv, tau, t, ldt)
    return ccall((@blasfunc(slarzt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), direct, storev, n,
                 k, v, ldv, tau, t, ldt, 1, 1)
end

function slas2(f, g, h, ssmin, ssmax)
    return ccall((@blasfunc(slas2_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), f,
                 g, h, ssmin, ssmax)
end

function slascl(type, kl, ku, cfrom, cto, m, n, a, lda, info)
    return ccall((@blasfunc(slascl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), type, kl, ku,
                 cfrom, cto, m, n, a, lda, info, 1)
end

function slascl2(m, n, d, x, ldx)
    return ccall((@blasfunc(slascl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), m, n, d, x,
                 ldx)
end

function slasd0(n, sqre, d, e, u, ldu, vt, ldvt, smlsiz, iwork, work, info)
    return ccall((@blasfunc(slasd0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}),
                 n, sqre, d, e, u, ldu, vt, ldvt, smlsiz, iwork, work, info)
end

function slasd1(nl, nr, sqre, d, alpha, beta, u, ldu, vt, ldvt, idxq, iwork, work, info)
    return ccall((@blasfunc(slasd1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}), nl, nr, sqre, d, alpha, beta, u, ldu, vt, ldvt,
                 idxq, iwork, work, info)
end

function slasd2(nl, nr, sqre, k, d, z, alpha, beta, u, ldu, vt, ldvt, dsigma, u2, ldu2, vt2,
                ldvt2, idxp, idx, idxc, idxq, coltyp, info)
    return ccall((@blasfunc(slasd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), nl, nr,
                 sqre, k, d, z, alpha, beta, u, ldu, vt, ldvt, dsigma, u2, ldu2, vt2, ldvt2,
                 idxp, idx, idxc, idxq, coltyp, info)
end

function slasd3(nl, nr, sqre, k, d, q, ldq, dsigma, u, ldu, u2, ldu2, vt, ldvt, vt2, ldvt2,
                idxc, ctot, z, info)
    return ccall((@blasfunc(slasd3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}), nl, nr, sqre, k, d, q, ldq, dsigma, u, ldu, u2,
                 ldu2, vt, ldvt, vt2, ldvt2, idxc, ctot, z, info)
end

function slasd4(n, i, d, z, delta, rho, sigma, work, info)
    return ccall((@blasfunc(slasd4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}), n, i, d, z, delta,
                 rho, sigma, work, info)
end

function slasd5(i, d, z, delta, rho, dsigma, work)
    return ccall((@blasfunc(slasd5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}), i, d, z, delta, rho, dsigma, work)
end

function slasd6(icompq, nl, nr, sqre, d, vf, vl, alpha, beta, idxq, perm, givptr, givcol,
                ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c, s, work, iwork, info)
    return ccall((@blasfunc(slasd6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}), icompq, nl, nr, sqre, d, vf, vl,
                 alpha, beta, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, poles,
                 difl, difr, z, k, c, s, work, iwork, info)
end

function slasd7(icompq, nl, nr, sqre, k, d, z, zw, vf, vfw, vl, vlw, alpha, beta, dsigma,
                idx, idxp, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info)
    return ccall((@blasfunc(slasd7_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}), icompq,
                 nl, nr, sqre, k, d, z, zw, vf, vfw, vl, vlw, alpha, beta, dsigma, idx,
                 idxp, idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info)
end

function slasd8(icompq, k, d, z, vf, vl, difl, difr, lddifr, dsigma, work, info)
    return ccall((@blasfunc(slasd8_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}), icompq, k, d, z, vf, vl, difl, difr, lddifr,
                 dsigma, work, info)
end

function slasda(icompq, smlsiz, n, sqre, d, e, u, ldu, vt, k, difl, difr, z, poles, givptr,
                givcol, ldgcol, perm, givnum, c, s, work, iwork, info)
    return ccall((@blasfunc(slasda_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}), icompq, smlsiz, n, sqre, d, e, u, ldu, vt, k, difl,
                 difr, z, poles, givptr, givcol, ldgcol, perm, givnum, c, s, work, iwork,
                 info)
end

function slasdq(uplo, sqre, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info)
    return ccall((@blasfunc(slasdq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo,
                 sqre, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info, 1)
end

function slasdt(n, lvl, nd, inode, ndiml, ndimr, msub)
    return ccall((@blasfunc(slasdt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, lvl, nd, inode, ndiml, ndimr, msub)
end

function slaset(uplo, m, n, alpha, beta, a, lda)
    return ccall((@blasfunc(slaset_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Clong), uplo, m, n, alpha, beta, a, lda, 1)
end

function slasq1(n, d, e, work, info)
    return ccall((@blasfunc(slasq1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), n, d, e,
                 work, info)
end

function slasq2(n, z, info)
    return ccall((@blasfunc(slasq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), n, z, info)
end

function slasq3(i0, n0, z, pp, dmin, sigma, desig, qmax, nfail, iter, ndiv, ieee, ttype,
                dmin1, dmin2, dn, dn1, dn2, g, tau)
    return ccall((@blasfunc(slasq3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}), i0, n0, z, pp, dmin, sigma,
                 desig, qmax, nfail, iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g,
                 tau)
end

function slasq4(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g)
    return ccall((@blasfunc(slasq4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}), i0, n0, z, pp, n0in, dmin, dmin1,
                 dmin2, dn, dn1, dn2, tau, ttype, g)
end

function slasq5(i0, n0, z, pp, tau, sigma, dmin, dmin1, dmin2, dn, dnm1, dnm2, ieee, eps)
    return ccall((@blasfunc(slasq5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{Float32}), i0, n0, z, pp, tau, sigma, dmin,
                 dmin1, dmin2, dn, dnm1, dnm2, ieee, eps)
end

function slasq6(i0, n0, z, pp, dmin, dmin1, dmin2, dn, dnm1, dnm2)
    return ccall((@blasfunc(slasq6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), i0, n0, z, pp,
                 dmin, dmin1, dmin2, dn, dnm1, dnm2)
end

function slasr(side, pivot, direct, m, n, c, s, a, lda)
    return ccall((@blasfunc(slasr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong), side, pivot,
                 direct, m, n, c, s, a, lda, 1, 1, 1)
end

function slasrt(id, n, d, info)
    return ccall((@blasfunc(slasrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), id, n, d, info, 1)
end

function slasv2(f, g, h, ssmin, ssmax, snr, csr, snl, csl)
    return ccall((@blasfunc(slasv2_), libblastrampoline), Cvoid,
                 (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Float32}), f, g, h, ssmin,
                 ssmax, snr, csr, snl, csl)
end

function slaswlq(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(slaswlq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, mb,
                 nb, a, lda, t, ldt, work, lwork, info)
end

function slaswp(n, a, lda, k1, k2, ipiv, incx)
    return ccall((@blasfunc(slaswp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, a, lda, k1, k2, ipiv, incx)
end

function slasy2(ltranl, ltranr, isgn, n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx,
                xnorm, info)
    return ccall((@blasfunc(slasy2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}), ltranl, ltranr, isgn,
                 n1, n2, tl, ldtl, tr, ldtr, b, ldb, scale, x, ldx, xnorm, info)
end

function slasyf(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(slasyf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb, a,
                 lda, ipiv, w, ldw, info, 1)
end

function slasyf_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work)
    return ccall((@blasfunc(slasyf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong), uplo, j1, m, nb,
                 a, lda, ipiv, h, ldh, work, 1)
end

function slasyf_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info)
    return ccall((@blasfunc(slasyf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nb, kb, a, lda, e, ipiv, w, ldw, info, 1)
end

function slasyf_rook(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(slasyf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb, a,
                 lda, ipiv, w, ldw, info, 1)
end

function slatbs(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info)
    return ccall((@blasfunc(slatbs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong), uplo, trans, diag, normin, n, kd,
                 ab, ldab, x, scale, cnorm, info, 1, 1, 1, 1)
end

function slatdf(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv)
    return ccall((@blasfunc(slatdf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ptr{BlasInt}, Ptr{BlasInt}), ijob, n, z, ldz, rhs, rdsum, rdscal,
                 ipiv, jpiv)
end

function slatps(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info)
    return ccall((@blasfunc(slatps_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), uplo, trans, diag, normin, n, ap, x, scale, cnorm, info, 1, 1, 1,
                 1)
end

function slatrd(uplo, n, nb, a, lda, e, tau, w, ldw)
    return ccall((@blasfunc(slatrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, nb, a, lda, e,
                 tau, w, ldw, 1)
end

function slatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info)
    return ccall((@blasfunc(slatrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong), uplo, trans, diag, normin, n, a, lda, x, scale,
                 cnorm, info, 1, 1, 1, 1)
end

function slatrs3(uplo, trans, diag, normin, n, nrhs, a, lda, x, ldx, scale, cnorm, work,
                 lwork, info)
    return ccall((@blasfunc(slatrs3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), uplo, trans, diag, normin, n, nrhs, a, lda, x, ldx, scale, cnorm,
                 work, lwork, info, 1, 1, 1, 1)
end

function slatrz(m, n, l, a, lda, tau, work)
    return ccall((@blasfunc(slatrz_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}), m, n, l, a, lda, tau, work)
end

function slatsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(slatsqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, mb,
                 nb, a, lda, t, ldt, work, lwork, info)
end

function slauu2(uplo, n, a, lda, info)
    return ccall((@blasfunc(slauu2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function slauum(uplo, n, a, lda, info)
    return ccall((@blasfunc(slauum_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function sopgtr(uplo, n, ap, tau, q, ldq, work, info)
    return ccall((@blasfunc(sopgtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, ap, tau, q, ldq,
                 work, info, 1)
end

function sopmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info)
    return ccall((@blasfunc(sopmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong, Clong), side, uplo, trans, m, n, ap, tau, c, ldc, work, info, 1, 1,
                 1)
end

function sorbdb(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22,
                theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info)
    return ccall((@blasfunc(sorbdb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22,
                 ldx22, theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info, 1, 1)
end

function sorbdb1(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(sorbdb1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, p, q, x11, ldx11,
                 x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info)
end

function sorbdb2(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(sorbdb2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, p, q, x11, ldx11,
                 x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info)
end

function sorbdb3(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(sorbdb3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, p, q, x11, ldx11,
                 x21, ldx21, theta, phi, taup1, taup2, tauq1, work, lwork, info)
end

function sorbdb4(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom,
                 work, lwork, info)
    return ccall((@blasfunc(sorbdb4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, p, q,
                 x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom, work,
                 lwork, info)
end

function sorbdb5(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
    return ccall((@blasfunc(sorbdb5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2,
                 ldq2, work, lwork, info)
end

function sorbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
    return ccall((@blasfunc(sorbdb6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2,
                 ldq2, work, lwork, info)
end

function sorcsd(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12,
                x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t,
                work, lwork, iwork, info)
    return ccall((@blasfunc(sorcsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong, Clong, Clong), jobu1, jobu2, jobv1t, jobv2t,
                 trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22,
                 theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, iwork,
                 info, 1, 1, 1, 1, 1, 1)
end

function sorcsd2by1(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1,
                    u2, ldu2, v1t, ldv1t, work, lwork, iwork, info)
    return ccall((@blasfunc(sorcsd2by1_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1, u2,
                 ldu2, v1t, ldv1t, work, lwork, iwork, info, 1, 1, 1)
end

function sorg2l(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(sorg2l_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}), m, n, k, a, lda, tau, work, info)
end

function sorg2r(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(sorg2r_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}), m, n, k, a, lda, tau, work, info)
end

function sorgbr(vect, m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sorgbr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), vect, m, n, k,
                 a, lda, tau, work, lwork, info, 1)
end

function sorghr(n, ilo, ihi, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sorghr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), n, ilo, ihi, a, lda, tau, work,
                 lwork, info)
end

function sorgl2(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(sorgl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}), m, n, k, a, lda, tau, work, info)
end

function sorglq(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sorglq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda, tau, work, lwork,
                 info)
end

function sorgql(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sorgql_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda, tau, work, lwork,
                 info)
end

function sorgqr(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sorgqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda, tau, work, lwork,
                 info)
end

function sorgr2(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(sorgr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}), m, n, k, a, lda, tau, work, info)
end

function sorgrq(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sorgrq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda, tau, work, lwork,
                 info)
end

function sorgtr(uplo, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(sorgtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, tau, work,
                 lwork, info, 1)
end

function sorgtsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(sorgtsqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, mb,
                 nb, a, lda, t, ldt, work, lwork, info)
end

function sorgtsqr_row(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(sorgtsqr_row_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, mb,
                 nb, a, lda, t, ldt, work, lwork, info)
end

function sorhr_col(m, n, nb, a, lda, t, ldt, d, info)
    return ccall((@blasfunc(sorhr_col_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}), m, n, nb, a, lda, t, ldt, d, info)
end

function sorm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sorm22_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, n1, n2, q, ldq, c, ldc, work,
                 lwork, info, 1, 1)
end

function sorm2l(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(sorm2l_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work, info, 1,
                 1)
end

function sorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(sorm2r_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work, info, 1,
                 1)
end

function sormbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sormbr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), vect, side,
                 trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info, 1, 1, 1)
end

function sormhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sormhr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, ilo,
                 ihi, a, lda, tau, c, ldc, work, lwork, info, 1, 1)
end

function sorml2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(sorml2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work, info, 1,
                 1)
end

function sormlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sormlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 lwork, info, 1, 1)
end

function sormql(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sormql_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 lwork, info, 1, 1)
end

function sormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sormqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 lwork, info, 1, 1)
end

function sormr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(sormr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work, info, 1,
                 1)
end

function sormr3(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(sormr3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, l, a, lda,
                 tau, c, ldc, work, info, 1, 1)
end

function sormrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sormrq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 lwork, info, 1, 1)
end

function sormrz(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sormrz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k,
                 l, a, lda, tau, c, ldc, work, lwork, info, 1, 1)
end

function sormtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(sormtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), side, uplo, trans, m, n, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1, 1)
end

function spbcon(uplo, n, kd, ab, ldab, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(spbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, ab,
                 ldab, anorm, rcond, work, iwork, info, 1)
end

function spbequ(uplo, n, kd, ab, ldab, s, scond, amax, info)
    return ccall((@blasfunc(spbequ_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Clong), uplo, n, kd, ab, ldab, s,
                 scond, amax, info, 1)
end

function spbrfs(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(spbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1)
end

function spbstf(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(spbstf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function spbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(spbsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab, ldab,
                 b, ldb, info, 1)
end

function spbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx,
                rcond, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(spbsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b,
                 ldb, x, ldx, rcond, ferr, berr, work, iwork, info, 1, 1, 1)
end

function spbtf2(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(spbtf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function spbtrf(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(spbtrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function spbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(spbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab, ldab,
                 b, ldb, info, 1)
end

function spftrf(transr, uplo, n, a, info)
    return ccall((@blasfunc(spftrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong),
                 transr, uplo, n, a, info, 1, 1)
end

function spftri(transr, uplo, n, a, info)
    return ccall((@blasfunc(spftri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong),
                 transr, uplo, n, a, info, 1, 1)
end

function spftrs(transr, uplo, n, nrhs, a, b, ldb, info)
    return ccall((@blasfunc(spftrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), transr, uplo, n, nrhs, a, b, ldb,
                 info, 1, 1)
end

function spocon(uplo, n, a, lda, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(spocon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 anorm, rcond, work, iwork, info, 1)
end

function spoequ(n, a, lda, s, scond, amax, info)
    return ccall((@blasfunc(spoequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}), n, a, lda, s, scond, amax, info)
end

function spoequb(n, a, lda, s, scond, amax, info)
    return ccall((@blasfunc(spoequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}), n, a, lda, s, scond, amax, info)
end

function sporfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, iwork,
                info)
    return ccall((@blasfunc(sporfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, iwork, info, 1)
end

function sporfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr,
                 n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, iwork,
                 info)
    return ccall((@blasfunc(sporfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, equed, n, nrhs, a, lda, af,
                 ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, iwork, info, 1, 1)
end

function sposv(uplo, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(sposv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, b, ldb, info, 1)
end

function sposvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                ferr, berr, work, iwork, info)
    return ccall((@blasfunc(sposvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                 ferr, berr, work, iwork, info, 1, 1, 1)
end

function sposvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                 rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params,
                 work, iwork, info)
    return ccall((@blasfunc(sposvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                 rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params,
                 work, iwork, info, 1, 1, 1)
end

function spotf2(uplo, n, a, lda, info)
    return ccall((@blasfunc(spotf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function spotrf(uplo, n, a, lda, info)
    return ccall((@blasfunc(spotrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function spotrf2(uplo, n, a, lda, info)
    return ccall((@blasfunc(spotrf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function spotri(uplo, n, a, lda, info)
    return ccall((@blasfunc(spotri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, a, lda, info, 1)
end

function spotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(spotrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, b, ldb, info, 1)
end

function sppcon(uplo, n, ap, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(sppcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, ap, anorm, rcond,
                 work, iwork, info, 1)
end

function sppequ(uplo, n, ap, s, scond, amax, info)
    return ccall((@blasfunc(sppequ_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Clong), uplo, n, ap, s, scond, amax, info, 1)
end

function spprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(spprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 ap, afp, b, ldb, x, ldx, ferr, berr, work, iwork, info, 1)
end

function sppsv(uplo, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(sppsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, nrhs, ap, b, ldb, info, 1)
end

function sppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr,
                work, iwork, info)
    return ccall((@blasfunc(sppsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{UInt8}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), fact, uplo, n, nrhs, ap, afp,
                 equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, iwork, info, 1, 1, 1)
end

function spptrf(uplo, n, ap, info)
    return ccall((@blasfunc(spptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, ap, info,
                 1)
end

function spptri(uplo, n, ap, info)
    return ccall((@blasfunc(spptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, ap, info,
                 1)
end

function spptrs(uplo, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(spptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, nrhs, ap, b, ldb, info, 1)
end

function spstf2(uplo, n, a, lda, piv, rank, tol, work, info)
    return ccall((@blasfunc(spstf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, a, lda, piv, rank,
                 tol, work, info, 1)
end

function spstrf(uplo, n, a, lda, piv, rank, tol, work, info)
    return ccall((@blasfunc(spstrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, a, lda, piv, rank,
                 tol, work, info, 1)
end

function sptcon(n, d, e, anorm, rcond, work, info)
    return ccall((@blasfunc(sptcon_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{Float32}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}), n, d, e, anorm, rcond, work, info)
end

function spteqr(compz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(spteqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work,
                 info, 1)
end

function sptrfs(n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, info)
    return ccall((@blasfunc(sptrfs_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), n, nrhs, d, e, df,
                 ef, b, ldb, x, ldx, ferr, berr, work, info)
end

function sptsv(n, nrhs, d, e, b, ldb, info)
    return ccall((@blasfunc(sptsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), n, nrhs, d, e, b, ldb, info)
end

function sptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, info)
    return ccall((@blasfunc(sptsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Clong), fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond,
                 ferr, berr, work, info, 1)
end

function spttrf(n, d, e, info)
    return ccall((@blasfunc(spttrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), n, d, e, info)
end

function spttrs(n, nrhs, d, e, b, ldb, info)
    return ccall((@blasfunc(spttrs_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}), n, nrhs, d, e, b, ldb, info)
end

function sptts2(n, nrhs, d, e, b, ldb)
    return ccall((@blasfunc(sptts2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}), n, nrhs, d, e, b, ldb)
end

function srscl(n, sa, sx, incx)
    return ccall((@blasfunc(srscl_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}), n, sa, sx, incx)
end

function ssb2st_kernels(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt,
                        work)
    return ccall((@blasfunc(ssb2st_kernels_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Clong), uplo, wantz, ttype, st, ed,
                 sweep, n, nb, ib, a, lda, v, tau, ldvt, work, 1)
end

function ssbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, info)
    return ccall((@blasfunc(ssbev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong), jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, info, 1, 1)
end

function ssbev_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, info)
    return ccall((@blasfunc(ssbev_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, info,
                 1, 1)
end

function ssbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(ssbevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z,
                 ldz, work, lwork, iwork, liwork, info, 1, 1)
end

function ssbevd_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork,
                       info)
    return ccall((@blasfunc(ssbevd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z,
                 ldz, work, lwork, iwork, liwork, info, 1, 1)
end

function ssbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z,
                ldz, work, iwork, ifail, info)
    return ccall((@blasfunc(ssbevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz,
                 range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz,
                 work, iwork, ifail, info, 1, 1, 1)
end

function ssbevx_2stage(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol,
                       m, w, z, ldz, work, lwork, iwork, ifail, info)
    return ccall((@blasfunc(ssbevx_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu,
                 abstol, m, w, z, ldz, work, lwork, iwork, ifail, info, 1, 1, 1)
end

function ssbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, info)
    return ccall((@blasfunc(ssbgst_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong), vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x,
                 ldx, work, info, 1, 1)
end

function ssbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, info)
    return ccall((@blasfunc(ssbgv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, ka, kb, ab, ldab,
                 bb, ldbb, w, z, ldz, work, info, 1, 1)
end

function ssbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, iwork,
                liwork, info)
    return ccall((@blasfunc(ssbgvd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, iwork,
                 liwork, info, 1, 1)
end

function ssbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu,
                abstol, m, w, z, ldz, work, iwork, ifail, info)
    return ccall((@blasfunc(ssbgvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo, n, ka, kb, ab, ldab,
                 bb, ldbb, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail,
                 info, 1, 1, 1)
end

function ssbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info)
    return ccall((@blasfunc(ssbtrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong), vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work,
                 info, 1, 1)
end

function ssfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c)
    return ccall((@blasfunc(ssfrk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Clong, Clong, Clong),
                 transr, uplo, trans, n, k, alpha, a, lda, beta, c, 1, 1, 1)
end

function sspcon(uplo, n, ap, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(sspcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, ap,
                 ipiv, anorm, rcond, work, iwork, info, 1)
end

function sspev(jobz, uplo, n, ap, w, z, ldz, work, info)
    return ccall((@blasfunc(sspev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), jobz,
                 uplo, n, ap, w, z, ldz, work, info, 1, 1)
end

function sspevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(sspevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), jobz, uplo, n, ap, w, z, ldz, work, lwork,
                 iwork, liwork, info, 1, 1)
end

function sspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork,
                ifail, info)
    return ccall((@blasfunc(sspevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m,
                 w, z, ldz, work, iwork, ifail, info, 1, 1, 1)
end

function sspgst(itype, uplo, n, ap, bp, info)
    return ccall((@blasfunc(sspgst_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Clong), itype, uplo, n, ap, bp, info, 1)
end

function sspgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info)
    return ccall((@blasfunc(sspgv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong), itype, jobz, uplo, n, ap, bp, w, z, ldz, work, info, 1, 1)
end

function sspgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(sspgvd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, ap, bp, w, z,
                 ldz, work, lwork, iwork, liwork, info, 1, 1)
end

function sspgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz,
                work, iwork, ifail, info)
    return ccall((@blasfunc(sspgvx_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz,
                 work, iwork, ifail, info, 1, 1, 1)
end

function ssprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info)
    return ccall((@blasfunc(ssprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, info, 1)
end

function sspsv(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(sspsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b, ldb, info, 1)
end

function sspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(sspsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr,
                 berr, work, iwork, info, 1, 1)
end

function ssptrd(uplo, n, ap, d, e, tau, info)
    return ccall((@blasfunc(ssptrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, ap, d, e, tau, info, 1)
end

function ssptrf(uplo, n, ap, ipiv, info)
    return ccall((@blasfunc(ssptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, ap, ipiv, info, 1)
end

function ssptri(uplo, n, ap, ipiv, work, info)
    return ccall((@blasfunc(ssptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Clong), uplo, n, ap, ipiv, work, info, 1)
end

function ssptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(ssptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b, ldb, info, 1)
end

function sstebz(range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit, w, iblock, isplit,
                work, iwork, info)
    return ccall((@blasfunc(sstebz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong), range, order, n, vl, vu, il, iu, abstol, d, e, m, nsplit,
                 w, iblock, isplit, work, iwork, info, 1, 1)
end

function sstedc(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(sstedc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info, 1)
end

function sstegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
                lwork, iwork, liwork, info)
    return ccall((@blasfunc(sstegr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n, d, e, vl,
                 vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info,
                 1, 1)
end

function sstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info)
    return ccall((@blasfunc(sstein_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail,
                 info)
end

function sstemr(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac,
                work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(sstemr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n,
                 d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork,
                 iwork, liwork, info, 1, 1)
end

function ssteqr(compz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(ssteqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work,
                 info, 1)
end

function ssterf(n, d, e, info)
    return ccall((@blasfunc(ssterf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}), n, d, e, info)
end

function sstev(jobz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(sstev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), jobz, n, d, e, z, ldz, work,
                 info, 1)
end

function sstevd(jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(sstevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), jobz, n, d, e, z, ldz, work, lwork, iwork, liwork, info, 1)
end

function sstevr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
                lwork, iwork, liwork, info)
    return ccall((@blasfunc(sstevr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n, d, e, vl,
                 vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info,
                 1, 1)
end

function sstevx(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork,
                ifail, info)
    return ccall((@blasfunc(sstevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong), jobz, range, n, d, e, vl, vu, il, iu, abstol, m,
                 w, z, ldz, work, iwork, ifail, info, 1, 1)
end

function ssycon(uplo, n, a, lda, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(ssycon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 ipiv, anorm, rcond, work, iwork, info, 1)
end

function ssycon_3(uplo, n, a, lda, e, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(ssycon_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, e, ipiv, anorm, rcond, work, iwork, info, 1)
end

function ssycon_rook(uplo, n, a, lda, ipiv, anorm, rcond, work, iwork, info)
    return ccall((@blasfunc(ssycon_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 ipiv, anorm, rcond, work, iwork, info, 1)
end

function ssyconv(uplo, way, n, a, lda, ipiv, e, info)
    return ccall((@blasfunc(ssyconv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a, lda, ipiv, e,
                 info, 1, 1)
end

function ssyconvf(uplo, way, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(ssyconvf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a, lda, e, ipiv, info,
                 1, 1)
end

function ssyconvf_rook(uplo, way, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(ssyconvf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a, lda, e, ipiv, info,
                 1, 1)
end

function ssyequb(uplo, n, a, lda, s, scond, amax, work, info)
    return ccall((@blasfunc(ssyequb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, a,
                 lda, s, scond, amax, work, info, 1)
end

function ssyev(jobz, uplo, n, a, lda, w, work, lwork, info)
    return ccall((@blasfunc(ssyev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, a, lda,
                 w, work, lwork, info, 1, 1)
end

function ssyev_2stage(jobz, uplo, n, a, lda, w, work, lwork, info)
    return ccall((@blasfunc(ssyev_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, a, lda,
                 w, work, lwork, info, 1, 1)
end

function ssyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(ssyevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info, 1, 1)
end

function ssyevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(ssyevd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info, 1, 1)
end

function ssyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz,
                work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(ssyevr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo,
                 n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,
                 iwork, liwork, info, 1, 1, 1)
end

function ssyevr_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                       isuppz, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(ssyevr_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo,
                 n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork,
                 iwork, liwork, info, 1, 1, 1)
end

function ssyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work,
                lwork, iwork, ifail, info)
    return ccall((@blasfunc(ssyevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo, n, a, lda,
                 vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info, 1,
                 1, 1)
end

function ssyevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                       work, lwork, iwork, ifail, info)
    return ccall((@blasfunc(ssyevx_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo, n, a, lda,
                 vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info, 1,
                 1, 1)
end

function ssygs2(itype, uplo, n, a, lda, b, ldb, info)
    return ccall((@blasfunc(ssygs2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), itype, uplo, n, a, lda, b, ldb, info, 1)
end

function ssygst(itype, uplo, n, a, lda, b, ldb, info)
    return ccall((@blasfunc(ssygst_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), itype, uplo, n, a, lda, b, ldb, info, 1)
end

function ssygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
    return ccall((@blasfunc(ssygv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info,
                 1, 1)
end

function ssygv_2stage(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
    return ccall((@blasfunc(ssygv_2stage_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info,
                 1, 1)
end

function ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(ssygvd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, a, lda, b, ldb,
                 w, work, lwork, iwork, liwork, info, 1, 1)
end

function ssygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w,
                z, ldz, work, lwork, iwork, ifail, info)
    return ccall((@blasfunc(ssygvx_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ref{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu,
                 abstol, m, w, z, ldz, work, lwork, iwork, ifail, info, 1, 1, 1)
end

function ssyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(ssyrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1)
end

function ssyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 iwork, info)
    return ccall((@blasfunc(ssyrfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, equed, n, nrhs, a, lda, af,
                 ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, iwork, info, 1, 1)
end

function ssysv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(ssysv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function ssysv_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(ssysv_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function ssysv_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork,
                         info)
    return ccall((@blasfunc(ssysv_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2,
                 b, ldb, work, lwork, info, 1)
end

function ssysv_rk(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(ssysv_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info, 1)
end

function ssysv_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(ssysv_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function ssysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr,
                berr, work, lwork, iwork, info)
    return ccall((@blasfunc(ssysvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), fact, uplo, n, nrhs, a,
                 lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, lwork, iwork,
                 info, 1, 1)
end

function ssysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, iwork, info)
    return ccall((@blasfunc(ssysvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x,
                 ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, iwork, info, 1, 1, 1)
end

function ssyswapr(uplo, n, a, lda, i1, i2)
    return ccall((@blasfunc(ssyswapr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, i1, i2, 1)
end

function ssytd2(uplo, n, a, lda, d, e, tau, info)
    return ccall((@blasfunc(ssytd2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, a, lda, d, e, tau,
                 info, 1)
end

function ssytf2(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(ssytf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function ssytf2_rk(uplo, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(ssytf2_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv, info, 1)
end

function ssytf2_rook(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(ssytf2_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function ssytrd(uplo, n, a, lda, d, e, tau, work, lwork, info)
    return ccall((@blasfunc(ssytrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, d, e, tau, work, lwork, info, 1)
end

function ssytrd_2stage(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info)
    return ccall((@blasfunc(ssytrd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), vect, uplo, n, a, lda, d, e, tau,
                 hous2, lhous2, work, lwork, info, 1, 1)
end

function ssytrd_sb2st(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork,
                      info)
    return ccall((@blasfunc(ssytrd_sb2st_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), stage1, vect,
                 uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info, 1, 1, 1)
end

function ssytrd_sy2sb(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info)
    return ccall((@blasfunc(ssytrd_sy2sb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, kd, a, lda, ab, ldab, tau, work, lwork, info, 1)
end

function ssytrf(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(ssytrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, lwork, info, 1)
end

function ssytrf_aa(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(ssytrf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, lwork, info, 1)
end

function ssytrf_aa_2stage(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info)
    return ccall((@blasfunc(ssytrf_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info, 1)
end

function ssytrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(ssytrf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv,
                 work, lwork, info, 1)
end

function ssytrf_rook(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(ssytrf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, lwork, info, 1)
end

function ssytri(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(ssytri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function ssytri2(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(ssytri2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, lwork, info, 1)
end

function ssytri2x(uplo, n, a, lda, ipiv, work, nb, info)
    return ccall((@blasfunc(ssytri2x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, nb, info, 1)
end

function ssytri_3(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(ssytri_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv,
                 work, lwork, info, 1)
end

function ssytri_3x(uplo, n, a, lda, e, ipiv, work, nb, info)
    return ccall((@blasfunc(ssytri_3x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv,
                 work, nb, info, 1)
end

function ssytri_rook(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(ssytri_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function ssytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(ssytrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, ipiv,
                 b, ldb, info, 1)
end

function ssytrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info)
    return ccall((@blasfunc(ssytrs2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 a, lda, ipiv, b, ldb, work, info, 1)
end

function ssytrs_3(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info)
    return ccall((@blasfunc(ssytrs_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a,
                 lda, e, ipiv, b, ldb, info, 1)
end

function ssytrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(ssytrs_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo,
                 n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function ssytrs_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info)
    return ccall((@blasfunc(ssytrs_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info, 1)
end

function ssytrs_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(ssytrs_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, ipiv,
                 b, ldb, info, 1)
end

function stbcon(norm, uplo, diag, n, kd, ab, ldab, rcond, work, iwork, info)
    return ccall((@blasfunc(stbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), norm, uplo, diag, n, kd, ab, ldab, rcond, work, iwork, info, 1, 1,
                 1)
end

function stbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr, work,
                iwork, info)
    return ccall((@blasfunc(stbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx,
                 ferr, berr, work, iwork, info, 1, 1, 1)
end

function stbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(stbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong), uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info, 1, 1, 1)
end

function stfsm(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb)
    return ccall((@blasfunc(stfsm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong, Clong), transr, side, uplo, trans, diag, m, n, alpha,
                 a, b, ldb, 1, 1, 1, 1, 1)
end

function stftri(transr, uplo, diag, n, a, info)
    return ccall((@blasfunc(stftri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong, Clong), transr, uplo, diag, n, a, info, 1, 1, 1)
end

function stfttp(transr, uplo, n, arf, ap, info)
    return ccall((@blasfunc(stfttp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), transr, uplo, n, arf, ap, info, 1, 1)
end

function stfttr(transr, uplo, n, arf, a, lda, info)
    return ccall((@blasfunc(stfttr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), transr, uplo, n, arf, a, lda, info, 1, 1)
end

function stgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, m, work,
                info)
    return ccall((@blasfunc(stgevc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong), side,
                 howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, m, work, info,
                 1, 1)
end

function stgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork,
                info)
    return ccall((@blasfunc(stgex2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), wantq, wantz,
                 n, a, lda, b, ldb, q, ldq, z, ldz, j1, n1, n2, work, lwork, info)
end

function stgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork,
                info)
    return ccall((@blasfunc(stgexc_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), wantq, wantz, n, a, lda,
                 b, ldb, q, ldq, z, ldz, ifst, ilst, work, lwork, info)
end

function stgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alphar, alphai, beta, q, ldq,
                z, ldz, m, pl, pr, dif, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(stgsen_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}), ijob, wantq, wantz, select, n, a, lda,
                 b, ldb, alphar, alphai, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork,
                 iwork, liwork, info)
end

function stgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u,
                ldu, v, ldv, q, ldq, work, ncycle, info)
    return ccall((@blasfunc(stgsja_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobu, jobv, jobq, m, p, n, k,
                 l, a, lda, b, ldb, tola, tolb, alpha, beta, u, ldu, v, ldv, q, ldq, work,
                 ncycle, info, 1, 1, 1)
end

function stgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m,
                work, lwork, iwork, info)
    return ccall((@blasfunc(stgsna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), job, howmny, select, n, a, lda, b,
                 ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work, lwork, iwork, info, 1, 1)
end

function stgsy2(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                rdsum, rdscal, iwork, pq, info)
    return ccall((@blasfunc(stgsy2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ref{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                 rdsum, rdscal, iwork, pq, info, 1)
end

function stgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                dif, work, lwork, iwork, info)
    return ccall((@blasfunc(stgsyl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                 dif, work, lwork, iwork, info, 1)
end

function stpcon(norm, uplo, diag, n, ap, rcond, work, iwork, info)
    return ccall((@blasfunc(stpcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{Float32},
                  Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), norm, uplo,
                 diag, n, ap, rcond, work, iwork, info, 1, 1, 1)
end

function stplqt(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info)
    return ccall((@blasfunc(stplqt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}), m, n, l, mb, a, lda, b, ldb, t, ldt, work, info)
end

function stplqt2(m, n, l, a, lda, b, ldb, t, ldt, info)
    return ccall((@blasfunc(stplqt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, l, a, lda, b, ldb,
                 t, ldt, info)
end

function stpmlqt(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info)
    return ccall((@blasfunc(stpmlqt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong), side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work,
                 info, 1, 1)
end

function stpmqrt(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info)
    return ccall((@blasfunc(stpmqrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong,
                  Clong), side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work,
                 info, 1, 1)
end

function stpqrt(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info)
    return ccall((@blasfunc(stpqrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}), m, n, l, nb, a, lda, b, ldb, t, ldt, work, info)
end

function stpqrt2(m, n, l, a, lda, b, ldb, t, ldt, info)
    return ccall((@blasfunc(stpqrt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}), m, n, l, a, lda, b, ldb,
                 t, ldt, info)
end

function stprfb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb,
                work, ldwork)
    return ccall((@blasfunc(stprfb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), side, trans, direct, storev, m, n, k, l, v,
                 ldv, t, ldt, a, lda, b, ldb, work, ldwork, 1, 1, 1, 1)
end

function stprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work, iwork,
                info)
    return ccall((@blasfunc(stprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work, iwork,
                 info, 1, 1, 1)
end

function stptri(uplo, diag, n, ap, info)
    return ccall((@blasfunc(stptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong, Clong),
                 uplo, diag, n, ap, info, 1, 1)
end

function stptrs(uplo, trans, diag, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(stptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, nrhs, ap, b, ldb, info, 1, 1, 1)
end

function stpttf(transr, uplo, n, ap, arf, info)
    return ccall((@blasfunc(stpttf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Clong, Clong), transr, uplo, n, ap, arf, info, 1, 1)
end

function stpttr(uplo, n, ap, a, lda, info)
    return ccall((@blasfunc(stpttr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, ap, a, lda, info, 1)
end

function strcon(norm, uplo, diag, n, a, lda, rcond, work, iwork, info)
    return ccall((@blasfunc(strcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 norm, uplo, diag, n, a, lda, rcond, work, iwork, info, 1, 1, 1)
end

function strevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, info)
    return ccall((@blasfunc(strevc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Clong, Clong), side, howmny, select, n, t, ldt,
                 vl, ldvl, vr, ldvr, mm, m, work, info, 1, 1)
end

function strevc3(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork,
                 info)
    return ccall((@blasfunc(strevc3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, howmny, select,
                 n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork, info, 1, 1)
end

function strexc(compq, n, t, ldt, q, ldq, ifst, ilst, work, info)
    return ccall((@blasfunc(strexc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Clong), compq, n, t, ldt,
                 q, ldq, ifst, ilst, work, info, 1)
end

function strrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr, work, iwork,
                info)
    return ccall((@blasfunc(strrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ptr{Float32}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr, work,
                 iwork, info, 1, 1, 1)
end

function strsen(job, compq, select, n, t, ldt, q, ldq, wr, wi, m, s, sep, work, lwork,
                iwork, liwork, info)
    return ccall((@blasfunc(strsen_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32}, Ref{BlasInt},
                  Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), job, compq, select, n, t, ldt, q, ldq, wr, wi,
                 m, s, sep, work, lwork, iwork, liwork, info, 1, 1)
end

function strsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, m, work,
                ldwork, iwork, info)
    return ccall((@blasfunc(strsna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong), job, howmny, select, n, t, ldt, vl, ldvl, vr,
                 ldvr, s, sep, mm, m, work, ldwork, iwork, info, 1, 1)
end

function strsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info)
    return ccall((@blasfunc(strsyl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ref{BlasInt}, Clong, Clong), trana, tranb, isgn, m, n, a, lda, b, ldb, c,
                 ldc, scale, info, 1, 1)
end

function strsyl3(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, iwork, liwork,
                 swork, ldswork, info)
    return ccall((@blasfunc(strsyl3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{Float32},
                  Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, iwork, liwork,
                 swork, ldswork, info, 1, 1)
end

function strti2(uplo, diag, n, a, lda, info)
    return ccall((@blasfunc(strti2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), uplo, diag, n, a, lda, info, 1, 1)
end

function strtri(uplo, diag, n, a, lda, info)
    return ccall((@blasfunc(strtri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), uplo, diag, n, a, lda, info, 1, 1)
end

function strtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(strtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), uplo,
                 trans, diag, n, nrhs, a, lda, b, ldb, info, 1, 1, 1)
end

function strttf(transr, uplo, n, a, lda, arf, info)
    return ccall((@blasfunc(strttf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32},
                  Ref{BlasInt}, Clong, Clong), transr, uplo, n, a, lda, arf, info, 1, 1)
end

function strttp(uplo, n, a, lda, ap, info)
    return ccall((@blasfunc(strttp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ap, info, 1)
end

function stzrzf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(stzrzf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float32}, Ref{BlasInt}, Ptr{Float32}, Ptr{Float32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork, info)
end

function zbbcsd(jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2,
                ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d,
                b22e, rwork, lrwork, info)
    return ccall((@blasfunc(zbbcsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong, Clong, Clong),
                 jobu1, jobu2, jobv1t, jobv2t, trans, m, p, q, theta, phi, u1, ldu1, u2,
                 ldu2, v1t, ldv1t, v2t, ldv2t, b11d, b11e, b12d, b12e, b21d, b21e, b22d,
                 b22e, rwork, lrwork, info, 1, 1, 1, 1, 1)
end

function zbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info)
    return ccall((@blasfunc(zbdsqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n,
                 ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, rwork, info, 1)
end

function zcgesv(n, nrhs, a, lda, ipiv, b, ldb, x, ldx, work, swork, rwork, iter, info)
    return ccall((@blasfunc(zcgesv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF32}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}), n, nrhs, a, lda,
                 ipiv, b, ldb, x, ldx, work, swork, rwork, iter, info)
end

function zcposv(uplo, n, nrhs, a, lda, b, ldb, x, ldx, work, swork, rwork, iter, info)
    return ccall((@blasfunc(zcposv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF32}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n,
                 nrhs, a, lda, b, ldb, x, ldx, work, swork, rwork, iter, info, 1)
end

function zdrscl(n, sa, sx, incx)
    return ccall((@blasfunc(zdrscl_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}), n, sa, sx, incx)
end

function zgbbrd(vect, m, n, ncc, kl, ku, ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work,
                rwork, info)
    return ccall((@blasfunc(zgbbrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong), vect, m, n, ncc, kl, ku,
                 ab, ldab, d, e, q, ldq, pt, ldpt, c, ldc, work, rwork, info, 1)
end

function zgbcon(norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(zgbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ref{BlasInt}, Clong), norm, n, kl, ku, ab, ldab, ipiv, anorm, rcond, work,
                 rwork, info, 1)
end

function zgbequ(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(zgbequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{BlasInt}), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
end

function zgbequb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(zgbequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{BlasInt}), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, info)
end

function zgbrfs(trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr,
                berr, work, rwork, info)
    return ccall((@blasfunc(zgbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong), trans, n, kl, ku, nrhs,
                 ab, ldab, afb, ldafb, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info,
                 1)
end

function zgbrfsx(trans, equed, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, r, c, b, ldb, x,
                 ldx, rcond, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info)
    return ccall((@blasfunc(zgbrfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong), trans, equed, n, kl, ku, nrhs, ab, ldab, afb,
                 ldafb, ipiv, r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, rwork, info, 1, 1)
end

function zgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    return ccall((@blasfunc(zgbsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), n, kl, ku, nrhs, ab,
                 ldab, ipiv, b, ldb, info)
end

function zgbsvx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
                ldb, x, ldx, rcond, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(zgbsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{UInt8}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong), fact,
                 trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x,
                 ldx, rcond, ferr, berr, work, rwork, info, 1, 1, 1)
end

function zgbsvxx(fact, trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b,
                 ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, rwork, info)
    return ccall((@blasfunc(zgbsvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{UInt8}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong), fact,
                 trans, n, kl, ku, nrhs, ab, ldab, afb, ldafb, ipiv, equed, r, c, b, ldb, x,
                 ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, rwork, info, 1, 1, 1)
end

function zgbtf2(m, n, kl, ku, ab, ldab, ipiv, info)
    return ccall((@blasfunc(zgbtf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, kl, ku, ab, ldab, ipiv, info)
end

function zgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
    return ccall((@blasfunc(zgbtrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, kl, ku, ab, ldab, ipiv, info)
end

function zgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    return ccall((@blasfunc(zgbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info, 1)
end

function zgebak(job, side, n, ilo, ihi, scale, m, v, ldv, info)
    return ccall((@blasfunc(zgebak_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job,
                 side, n, ilo, ihi, scale, m, v, ldv, info, 1, 1)
end

function zgebal(job, n, a, lda, ilo, ihi, scale, info)
    return ccall((@blasfunc(zgebal_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong), job, n, a, lda, ilo, ihi, scale, info, 1)
end

function zgebd2(m, n, a, lda, d, e, tauq, taup, work, info)
    return ccall((@blasfunc(zgebd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}), m, n, a, lda, d, e, tauq, taup, work, info)
end

function zgebrd(m, n, a, lda, d, e, tauq, taup, work, lwork, info)
    return ccall((@blasfunc(zgebrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, d, e, tauq, taup, work, lwork, info)
end

function zgecon(norm, n, a, lda, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(zgecon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong), norm, n,
                 a, lda, anorm, rcond, work, rwork, info, 1)
end

function zgeequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(zgeequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}), m, n,
                 a, lda, r, c, rowcnd, colcnd, amax, info)
end

function zgeequb(m, n, a, lda, r, c, rowcnd, colcnd, amax, info)
    return ccall((@blasfunc(zgeequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}), m, n,
                 a, lda, r, c, rowcnd, colcnd, amax, info)
end

function zgees(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork,
               info)
    return ccall((@blasfunc(zgees_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), jobvs, sort,
                 select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info, 1,
                 1)
end

function zgeesx(jobvs, sort, select, sense, n, a, lda, sdim, w, vs, ldvs, rconde, rcondv,
                work, lwork, rwork, bwork, info)
    return ccall((@blasfunc(zgeesx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobvs, sort, select, sense, n,
                 a, lda, sdim, w, vs, ldvs, rconde, rcondv, work, lwork, rwork, bwork, info,
                 1, 1, 1)
end

function zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    return ccall((@blasfunc(zgeev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), jobvl,
                 jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info, 1, 1)
end

function zgeevx(balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi,
                scale, abnrm, rconde, rcondv, work, lwork, rwork, info)
    return ccall((@blasfunc(zgeevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong, Clong), balanc, jobvl,
                 jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, ilo, ihi, scale, abnrm,
                 rconde, rcondv, work, lwork, rwork, info, 1, 1, 1, 1)
end

function zgehd2(n, ilo, ihi, a, lda, tau, work, info)
    return ccall((@blasfunc(zgehd2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}), n, ilo, ihi, a, lda, tau,
                 work, info)
end

function zgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zgehrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), n, ilo, ihi, a,
                 lda, tau, work, lwork, info)
end

function zgejsv(joba, jobu, jobv, jobr, jobt, jobp, m, n, a, lda, sva, u, ldu, v, ldv,
                cwork, lwork, rwork, lrwork, iwork, info)
    return ccall((@blasfunc(zgejsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong,
                  Clong, Clong, Clong, Clong), joba, jobu, jobv, jobr, jobt, jobp, m, n, a,
                 lda, sva, u, ldu, v, ldv, cwork, lwork, rwork, lrwork, iwork, info, 1, 1,
                 1, 1, 1, 1)
end

function zgelq(m, n, a, lda, t, tsize, work, lwork, info)
    return ccall((@blasfunc(zgelq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, tsize,
                 work, lwork, info)
end

function zgelq2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(zgelq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function zgelqf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zgelqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function zgelqt(m, n, mb, a, lda, t, ldt, work, info)
    return ccall((@blasfunc(zgelqt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, mb, a, lda,
                 t, ldt, work, info)
end

function zgelqt3(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(zgelqt3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function zgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zgels_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function zgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(zgelsd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), m, n,
                 nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, iwork, info)
end

function zgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, rwork, info)
    return ccall((@blasfunc(zgelss_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), m, n, nrhs, a, lda,
                 b, ldb, s, rcond, rank, work, lwork, rwork, info)
end

function zgelst(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zgelst_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function zgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, rwork, info)
    return ccall((@blasfunc(zgelsy_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), m, n, nrhs, a, lda,
                 b, ldb, jpvt, rcond, rank, work, lwork, rwork, info)
end

function zgemlq(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zgemlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, a, lda, t, tsize, c, ldc, work, lwork, info, 1, 1)
end

function zgemlqt(side, trans, m, n, k, mb, v, ldv, t, ldt, c, ldc, work, info)
    return ccall((@blasfunc(zgemlqt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, mb, v, ldv, t, ldt, c, ldc, work, info, 1, 1)
end

function zgemqr(side, trans, m, n, k, a, lda, t, tsize, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zgemqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, a, lda, t, tsize, c, ldc, work, lwork, info, 1, 1)
end

function zgemqrt(side, trans, m, n, k, nb, v, ldv, t, ldt, c, ldc, work, info)
    return ccall((@blasfunc(zgemqrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, nb, v, ldv, t, ldt, c, ldc, work, info, 1, 1)
end

function zgeql2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(zgeql2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function zgeqlf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zgeqlf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function zgeqp3(m, n, a, lda, jpvt, tau, work, lwork, rwork, info)
    return ccall((@blasfunc(zgeqp3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}), m,
                 n, a, lda, jpvt, tau, work, lwork, rwork, info)
end

function zgeqp3rk(m, n, nrhs, kmax, abstol, reltol, a, lda, k, maxc2nrmk, relmaxc2nrmk,
                  jpiv, tau, work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(zgeqp3rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs, kmax, abstol, reltol, a, lda, k,
                 maxc2nrmk, relmaxc2nrmk, jpiv, tau, work, lwork, rwork, iwork, info)
end

function zgeqr(m, n, a, lda, t, tsize, work, lwork, info)
    return ccall((@blasfunc(zgeqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, tsize,
                 work, lwork, info)
end

function zgeqr2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(zgeqr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function zgeqr2p(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(zgeqr2p_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function zgeqrf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zgeqrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function zgeqrfp(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zgeqrfp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function zgeqrt(m, n, nb, a, lda, t, ldt, work, info)
    return ccall((@blasfunc(zgeqrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, nb, a, lda,
                 t, ldt, work, info)
end

function zgeqrt2(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(zgeqrt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function zgeqrt3(m, n, a, lda, t, ldt, info)
    return ccall((@blasfunc(zgeqrt3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, t, ldt, info)
end

function zgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(zgerfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{Float64}, Ref{BlasInt}, Clong), trans, n, nrhs, a, lda, af, ldaf, ipiv,
                 b, ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function zgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 rwork, info)
    return ccall((@blasfunc(zgerfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong),
                 trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, r, c, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 rwork, info, 1, 1)
end

function zgerq2(m, n, a, lda, tau, work, info)
    return ccall((@blasfunc(zgerq2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}), m, n, a, lda, tau, work, info)
end

function zgerqf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zgerqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function zgesc2(n, a, lda, rhs, ipiv, jpiv, scale)
    return ccall((@blasfunc(zgesc2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{Float64}), n, a, lda, rhs, ipiv, jpiv, scale)
end

function zgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(zgesdd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info,
                 1)
end

function zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(zgesv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), n, nrhs, a, lda, ipiv, b, ldb,
                 info)
end

function zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
    return ccall((@blasfunc(zgesvd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), jobu,
                 jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info, 1, 1)
end

function zgesvdq(joba, jobp, jobr, jobu, jobv, m, n, a, lda, s, u, ldu, v, ldv, numrank,
                 iwork, liwork, cwork, lcwork, rwork, lrwork, info)
    return ccall((@blasfunc(zgesvdq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong, Clong), joba, jobp, jobr, jobu, jobv, m, n, a, lda,
                 s, u, ldu, v, ldv, numrank, iwork, liwork, cwork, lcwork, rwork, lrwork,
                 info, 1, 1, 1, 1, 1)
end

function zgesvdx(jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u, ldu, vt, ldvt,
                 work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(zgesvdx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), jobu, jobvt, range, m, n, a, lda, vl, vu, il, iu, ns, s, u,
                 ldu, vt, ldvt, work, lwork, rwork, iwork, info, 1, 1, 1)
end

function zgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, cwork, lwork, rwork,
                lrwork, info)
    return ccall((@blasfunc(zgesvj_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, cwork,
                 lwork, rwork, lrwork, info, 1, 1, 1)
end

function zgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                rcond, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(zgesvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, nrhs, a, lda, af, ldaf,
                 ipiv, equed, r, c, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1,
                 1, 1)
end

function zgesvxx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info)
    return ccall((@blasfunc(zgesvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong, Clong), fact, trans, n, nrhs, a, lda, af, ldaf,
                 ipiv, equed, r, c, b, ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds,
                 err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, 1, 1, 1)
end

function zgetc2(n, a, lda, ipiv, jpiv, info)
    return ccall((@blasfunc(zgetc2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                 n, a, lda, ipiv, jpiv, info)
end

function zgetf2(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(zgetf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                 m, n, a, lda, ipiv, info)
end

function zgetrf(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(zgetrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                 m, n, a, lda, ipiv, info)
end

function zgetrf2(m, n, a, lda, ipiv, info)
    return ccall((@blasfunc(zgetrf2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}),
                 m, n, a, lda, ipiv, info)
end

function zgetri(n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(zgetri_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}), n, a, lda, ipiv, work, lwork, info)
end

function zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(zgetrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), trans, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function zgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zgetsls_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info, 1)
end

function zgetsqrhrt(m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(zgetsqrhrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{BlasInt}), m, n, mb1, nb1, nb2, a, lda, t, ldt, work, lwork, info)
end

function zggbak(job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info)
    return ccall((@blasfunc(zggbak_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), job, side, n, ilo, ihi, lscale, rscale, m, v, ldv, info, 1, 1)
end

function zggbal(job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work, info)
    return ccall((@blasfunc(zggbal_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Clong), job, n, a, lda, b, ldb, ilo, ihi, lscale, rscale, work,
                 info, 1)
end

function zgges(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl,
               ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info)
    return ccall((@blasfunc(zgges_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim,
                 alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info, 1, 1,
                 1)
end

function zgges3(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, alpha, beta, vsl,
                ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info)
    return ccall((@blasfunc(zgges3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim,
                 alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, bwork, info, 1, 1,
                 1)
end

function zggesx(jobvsl, jobvsr, sort, selctg, sense, n, a, lda, b, ldb, sdim, alpha, beta,
                vsl, ldvsl, vsr, ldvsr, rconde, rcondv, work, lwork, rwork, iwork, liwork,
                bwork, info)
    return ccall((@blasfunc(zggesx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), jobvsl, jobvsr, sort, selctg, sense, n, a,
                 lda, b, ldb, sdim, alpha, beta, vsl, ldvsl, vsr, ldvsr, rconde, rcondv,
                 work, lwork, rwork, iwork, liwork, bwork, info, 1, 1, 1, 1)
end

function zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work,
               lwork, rwork, info)
    return ccall((@blasfunc(zggev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), jobvl, jobvr, n, a,
                 lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info, 1,
                 1)
end

function zggev3(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work,
                lwork, rwork, info)
    return ccall((@blasfunc(zggev3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), jobvl, jobvr, n, a,
                 lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info, 1,
                 1)
end

function zggevx(balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr,
                ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, lwork,
                rwork, iwork, bwork, info)
    return ccall((@blasfunc(zggevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), balanc, jobvl, jobvr, sense, n, a, lda, b, ldb, alpha, beta, vl,
                 ldvl, vr, ldvr, ilo, ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv,
                 work, lwork, rwork, iwork, bwork, info, 1, 1, 1, 1)
end

function zggglm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info)
    return ccall((@blasfunc(zggglm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), n, m, p, a, lda,
                 b, ldb, d, x, y, work, lwork, info)
end

function zgghd3(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work, lwork,
                info)
    return ccall((@blasfunc(zgghd3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, work,
                 lwork, info, 1, 1)
end

function zgghrd(compq, compz, n, ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info)
    return ccall((@blasfunc(zgghrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), compq, compz, n,
                 ilo, ihi, a, lda, b, ldb, q, ldq, z, ldz, info, 1, 1)
end

function zgglse(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    return ccall((@blasfunc(zgglse_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, p, a, lda,
                 b, ldb, c, d, x, work, lwork, info)
end

function zggqrf(n, m, p, a, lda, taua, b, ldb, taub, work, lwork, info)
    return ccall((@blasfunc(zggqrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), n, m, p, a, lda, taua, b, ldb,
                 taub, work, lwork, info)
end

function zggrqf(m, p, n, a, lda, taua, b, ldb, taub, work, lwork, info)
    return ccall((@blasfunc(zggrqf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, p, n, a, lda, taua, b, ldb,
                 taub, work, lwork, info)
end

function zggsvd3(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v,
                 ldv, q, ldq, work, lwork, rwork, iwork, info)
    return ccall((@blasfunc(zggsvd3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobu,
                 jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q,
                 ldq, work, lwork, rwork, iwork, info, 1, 1, 1)
end

function zggsvp3(jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, v,
                 ldv, q, ldq, iwork, rwork, tau, work, lwork, info)
    return ccall((@blasfunc(zggsvp3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), jobu, jobv, jobq, m, p, n, a, lda, b, ldb, tola,
                 tolb, k, l, u, ldu, v, ldv, q, ldq, iwork, rwork, tau, work, lwork, info,
                 1, 1, 1)
end

function zgsvj0(jobv, m, n, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work,
                lwork, info)
    return ccall((@blasfunc(zgsvj0_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), jobv, m, n, a, lda, d, sva, mv, v, ldv, eps,
                 sfmin, tol, nsweep, work, lwork, info, 1)
end

function zgsvj1(jobv, m, n, n1, a, lda, d, sva, mv, v, ldv, eps, sfmin, tol, nsweep, work,
                lwork, info)
    return ccall((@blasfunc(zgsvj1_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), jobv, m, n, n1, a, lda, d, sva, mv, v, ldv,
                 eps, sfmin, tol, nsweep, work, lwork, info, 1)
end

function zgtcon(norm, n, dl, d, du, du2, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(zgtcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong), norm, n, dl, d, du, du2, ipiv, anorm, rcond, work,
                 info, 1)
end

function zgtrfs(trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, ferr,
                berr, work, rwork, info)
    return ccall((@blasfunc(zgtrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ref{BlasInt}, Clong), trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b,
                 ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function zgtsv(n, nrhs, dl, d, du, b, ldb, info)
    return ccall((@blasfunc(zgtsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), n, nrhs, dl, d, du, b, ldb, info)
end

function zgtsvx(fact, trans, n, nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx,
                rcond, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(zgtsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), fact, trans, n,
                 nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, ldb, x, ldx, rcond, ferr,
                 berr, work, rwork, info, 1, 1)
end

function zgttrf(n, dl, d, du, du2, ipiv, info)
    return ccall((@blasfunc(zgttrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ref{BlasInt}), n, dl, d, du, du2, ipiv, info)
end

function zgttrs(trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info)
    return ccall((@blasfunc(zgttrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), trans, n, nrhs, dl, d, du, du2, ipiv, b, ldb, info, 1)
end

function zgtts2(itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb)
    return ccall((@blasfunc(zgtts2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}),
                 itrans, n, nrhs, dl, d, du, du2, ipiv, b, ldb)
end

function zhb2st_kernels(uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt,
                        work)
    return ccall((@blasfunc(zhb2st_kernels_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong),
                 uplo, wantz, ttype, st, ed, sweep, n, nb, ib, a, lda, v, tau, ldvt, work,
                 1)
end

function zhbev(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, rwork, info)
    return ccall((@blasfunc(zhbev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z, ldz, work,
                 rwork, info, 1, 1)
end

function zhbev_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, info)
    return ccall((@blasfunc(zhbev_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, kd, ab, ldab, w, z,
                 ldz, work, lwork, rwork, info, 1, 1)
end

function zhbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork,
                liwork, info)
    return ccall((@blasfunc(zhbevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork,
                 liwork, info, 1, 1)
end

function zhbevd_2stage(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork,
                       iwork, liwork, info)
    return ccall((@blasfunc(zhbevd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, rwork, lrwork, iwork,
                 liwork, info, 1, 1)
end

function zhbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z,
                ldz, work, rwork, iwork, ifail, info)
    return ccall((@blasfunc(zhbevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo, n, kd, ab,
                 ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork,
                 ifail, info, 1, 1, 1)
end

function zhbevx_2stage(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol,
                       m, w, z, ldz, work, lwork, rwork, iwork, ifail, info)
    return ccall((@blasfunc(zhbevx_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo,
                 n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                 rwork, iwork, ifail, info, 1, 1, 1)
end

function zhbgst(vect, uplo, n, ka, kb, ab, ldab, bb, ldbb, x, ldx, work, rwork, info)
    return ccall((@blasfunc(zhbgst_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), vect, uplo, n,
                 ka, kb, ab, ldab, bb, ldbb, x, ldx, work, rwork, info, 1, 1)
end

function zhbgv(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, rwork, info)
    return ccall((@blasfunc(zhbgv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), jobz,
                 uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, rwork, info, 1, 1)
end

function zhbgvd(jobz, uplo, n, ka, kb, ab, ldab, bb, ldbb, w, z, ldz, work, lwork, rwork,
                lrwork, iwork, liwork, info)
    return ccall((@blasfunc(zhbgvd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, ka, kb, ab, ldab, bb,
                 ldbb, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function zhbgvx(jobz, range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu,
                abstol, m, w, z, ldz, work, rwork, iwork, ifail, info)
    return ccall((@blasfunc(zhbgvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz,
                 range, uplo, n, ka, kb, ab, ldab, bb, ldbb, q, ldq, vl, vu, il, iu, abstol,
                 m, w, z, ldz, work, rwork, iwork, ifail, info, 1, 1, 1)
end

function zhbtrd(vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work, info)
    return ccall((@blasfunc(zhbtrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong, Clong), vect, uplo, n, kd, ab, ldab, d, e, q, ldq, work,
                 info, 1, 1)
end

function zhecon(uplo, n, a, lda, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(zhecon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, ipiv, anorm, rcond, work, info, 1)
end

function zhecon_3(uplo, n, a, lda, e, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(zhecon_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, e, ipiv, anorm, rcond, work, info, 1)
end

function zhecon_rook(uplo, n, a, lda, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(zhecon_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, ipiv, anorm, rcond, work, info, 1)
end

function zheequb(uplo, n, a, lda, s, scond, amax, work, info)
    return ccall((@blasfunc(zheequb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, s, scond, amax, work, info, 1)
end

function zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
    return ccall((@blasfunc(zheev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong), jobz, uplo, n, a, lda, w, work, lwork, rwork, info, 1, 1)
end

function zheev_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
    return ccall((@blasfunc(zheev_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong), jobz, uplo, n, a, lda, w, work, lwork, rwork, info, 1, 1)
end

function zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    return ccall((@blasfunc(zheevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, a, lda, w,
                 work, lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function zheevd_2stage(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork,
                       info)
    return ccall((@blasfunc(zheevd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n, a, lda, w,
                 work, lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function zheevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz,
                work, lwork, rwork, lrwork, iwork, liwork, info)
    return ccall((@blasfunc(zheevr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), jobz, range, uplo, n, a, lda, vl, vu, il, iu,
                 abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork,
                 info, 1, 1, 1)
end

function zheevr_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                       isuppz, work, lwork, rwork, lrwork, iwork, liwork, info)
    return ccall((@blasfunc(zheevr_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong), jobz, range, uplo, n, a, lda, vl, vu, il, iu,
                 abstol, m, w, z, ldz, isuppz, work, lwork, rwork, lrwork, iwork, liwork,
                 info, 1, 1, 1)
end

function zheevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work,
                lwork, rwork, iwork, ifail, info)
    return ccall((@blasfunc(zheevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz,
                 range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                 rwork, iwork, ifail, info, 1, 1, 1)
end

function zheevx_2stage(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz,
                       work, lwork, rwork, iwork, ifail, info)
    return ccall((@blasfunc(zheevx_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz,
                 range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                 rwork, iwork, ifail, info, 1, 1, 1)
end

function zhegs2(itype, uplo, n, a, lda, b, ldb, info)
    return ccall((@blasfunc(zhegs2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), itype, uplo, n, a, lda, b,
                 ldb, info, 1)
end

function zhegst(itype, uplo, n, a, lda, b, ldb, info)
    return ccall((@blasfunc(zhegst_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), itype, uplo, n, a, lda, b,
                 ldb, info, 1)
end

function zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    return ccall((@blasfunc(zhegv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, a, lda, b,
                 ldb, w, work, lwork, rwork, info, 1, 1)
end

function zhegv_2stage(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    return ccall((@blasfunc(zhegv_2stage_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, a, lda, b,
                 ldb, w, work, lwork, rwork, info, 1, 1)
end

function zhegvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork,
                liwork, info)
    return ccall((@blasfunc(zhegvd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
                 itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork,
                 liwork, info, 1, 1)
end

function zhegvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w,
                z, ldz, work, lwork, rwork, iwork, ifail, info)
    return ccall((@blasfunc(zhegvx_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), itype, jobz, range,
                 uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork,
                 rwork, iwork, ifail, info, 1, 1, 1)
end

function zherfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(zherfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, af, ldaf, ipiv, b,
                 ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function zherfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 rwork, info)
    return ccall((@blasfunc(zherfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), uplo, equed, n,
                 nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds,
                 err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, 1, 1)
end

function zhesv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zhesv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function zhesv_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zhesv_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function zhesv_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork,
                         info)
    return ccall((@blasfunc(zhesv_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info, 1)
end

function zhesv_rk(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zhesv_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, e, ipiv, b, ldb,
                 work, lwork, info, 1)
end

function zhesv_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zhesv_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function zhesvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr,
                berr, work, lwork, rwork, info)
    return ccall((@blasfunc(zhesvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), fact,
                 uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr,
                 work, lwork, rwork, info, 1, 1)
end

function zhesvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info)
    return ccall((@blasfunc(zhesvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong, Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b,
                 ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, rwork, info, 1, 1, 1)
end

function zheswapr(uplo, n, a, lda, i1, i2)
    return ccall((@blasfunc(zheswapr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, i1, i2, 1)
end

function zhetd2(uplo, n, a, lda, d, e, tau, info)
    return ccall((@blasfunc(zhetd2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, a, lda, d, e,
                 tau, info, 1)
end

function zhetf2(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(zhetf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function zhetf2_rk(uplo, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(zhetf2_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv, info, 1)
end

function zhetf2_rook(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(zhetf2_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function zhetrd(uplo, n, a, lda, d, e, tau, work, lwork, info)
    return ccall((@blasfunc(zhetrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, d, e, tau, work, lwork, info, 1)
end

function zhetrd_2stage(vect, uplo, n, a, lda, d, e, tau, hous2, lhous2, work, lwork, info)
    return ccall((@blasfunc(zhetrd_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), vect, uplo, n, a,
                 lda, d, e, tau, hous2, lhous2, work, lwork, info, 1, 1)
end

function zhetrd_hb2st(stage1, vect, uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork,
                      info)
    return ccall((@blasfunc(zhetrd_hb2st_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), stage1, vect,
                 uplo, n, kd, ab, ldab, d, e, hous, lhous, work, lwork, info, 1, 1, 1)
end

function zhetrd_he2hb(uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info)
    return ccall((@blasfunc(zhetrd_he2hb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, kd, a, lda, ab, ldab, tau, work, lwork, info,
                 1)
end

function zhetrf(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(zhetrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function zhetrf_aa(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(zhetrf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function zhetrf_aa_2stage(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info)
    return ccall((@blasfunc(zhetrf_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info, 1)
end

function zhetrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(zhetrf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, lwork, info, 1)
end

function zhetrf_rook(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(zhetrf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function zhetri(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(zhetri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function zhetri2(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(zhetri2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function zhetri2x(uplo, n, a, lda, ipiv, work, nb, info)
    return ccall((@blasfunc(zhetri2x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, nb, info, 1)
end

function zhetri_3(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(zhetri_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, lwork, info, 1)
end

function zhetri_3x(uplo, n, a, lda, e, ipiv, work, nb, info)
    return ccall((@blasfunc(zhetri_3x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, nb, info, 1)
end

function zhetri_rook(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(zhetri_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function zhetrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(zhetrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function zhetrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info)
    return ccall((@blasfunc(zhetrs2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 nrhs, a, lda, ipiv, b, ldb, work, info, 1)
end

function zhetrs_3(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info)
    return ccall((@blasfunc(zhetrs_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info, 1)
end

function zhetrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zhetrs_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function zhetrs_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info)
    return ccall((@blasfunc(zhetrs_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2,
                 b, ldb, info, 1)
end

function zhetrs_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(zhetrs_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function zhfrk(transr, uplo, trans, n, k, alpha, a, lda, beta, c)
    return ccall((@blasfunc(zhfrk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Clong, Clong,
                  Clong), transr, uplo, trans, n, k, alpha, a, lda, beta, c, 1, 1, 1)
end

function zhgeqz(job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z, ldz,
                work, lwork, rwork, info)
    return ccall((@blasfunc(zhgeqz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong),
                 job, compq, compz, n, ilo, ihi, h, ldh, t, ldt, alpha, beta, q, ldq, z,
                 ldz, work, lwork, rwork, info, 1, 1, 1)
end

function zhpcon(uplo, n, ap, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(zhpcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, ap, ipiv,
                 anorm, rcond, work, info, 1)
end

function zhpev(jobz, uplo, n, ap, w, z, ldz, work, rwork, info)
    return ccall((@blasfunc(zhpev_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), jobz, uplo, n, ap, w, z, ldz, work, rwork, info, 1, 1)
end

function zhpevd(jobz, uplo, n, ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork,
                info)
    return ccall((@blasfunc(zhpevd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, uplo, n,
                 ap, w, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function zhpevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork,
                iwork, ifail, info)
    return ccall((@blasfunc(zhpevx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), jobz, range, uplo,
                 n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, rwork, iwork, ifail,
                 info, 1, 1, 1)
end

function zhpgst(itype, uplo, n, ap, bp, info)
    return ccall((@blasfunc(zhpgst_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong), itype, uplo, n, ap, bp, info, 1)
end

function zhpgv(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, rwork, info)
    return ccall((@blasfunc(zhpgv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), itype, jobz,
                 uplo, n, ap, bp, w, z, ldz, work, rwork, info, 1, 1)
end

function zhpgvd(itype, jobz, uplo, n, ap, bp, w, z, ldz, work, lwork, rwork, lrwork, iwork,
                liwork, info)
    return ccall((@blasfunc(zhpgvd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Clong, Clong), itype, jobz, uplo, n, ap, bp, w, z, ldz, work,
                 lwork, rwork, lrwork, iwork, liwork, info, 1, 1)
end

function zhpgvx(itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu, abstol, m, w, z, ldz,
                work, rwork, iwork, ifail, info)
    return ccall((@blasfunc(zhpgvx_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), itype, jobz, range, uplo, n, ap, bp, vl, vu, il, iu,
                 abstol, m, w, z, ldz, work, rwork, iwork, ifail, info, 1, 1, 1)
end

function zhprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(zhprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt},
                  Clong), uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 rwork, info, 1)
end

function zhpsv(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(zhpsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b,
                 ldb, info, 1)
end

function zhpsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(zhpsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), fact, uplo, n, nrhs, ap, afp,
                 ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1, 1)
end

function zhptrd(uplo, n, ap, d, e, tau, info)
    return ccall((@blasfunc(zhptrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, ap, d, e, tau, info, 1)
end

function zhptrf(uplo, n, ap, ipiv, info)
    return ccall((@blasfunc(zhptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, ap, ipiv, info, 1)
end

function zhptri(uplo, n, ap, ipiv, work, info)
    return ccall((@blasfunc(zhptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong), uplo, n, ap, ipiv, work, info, 1)
end

function zhptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(zhptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b,
                 ldb, info, 1)
end

function zhsein(side, eigsrc, initv, select, n, h, ldh, w, vl, ldvl, vr, ldvr, mm, m, work,
                rwork, ifaill, ifailr, info)
    return ccall((@blasfunc(zhsein_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), side, eigsrc, initv, select,
                 n, h, ldh, w, vl, ldvl, vr, ldvr, mm, m, work, rwork, ifaill, ifailr, info,
                 1, 1, 1)
end

function zhseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info)
    return ccall((@blasfunc(zhseqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job, compz, n, ilo, ihi, h, ldh, w,
                 z, ldz, work, lwork, info, 1, 1)
end

function zla_gbamv(trans, m, n, kl, ku, alpha, ab, ldab, x, incx, beta, y, incy)
    return ccall((@blasfunc(zla_gbamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}), trans, m, n, kl, ku, alpha, ab, ldab, x, incx,
                 beta, y, incy)
end

function zla_gbrcond_c(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, c, capply, info, work,
                       rwork)
    return ccall((@blasfunc(zla_gbrcond_c_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Clong), trans, n, kl, ku, ab, ldab, afb,
                 ldafb, ipiv, c, capply, info, work, rwork, 1)
end

function zla_gbrcond_x(trans, n, kl, ku, ab, ldab, afb, ldafb, ipiv, x, info, work, rwork)
    return ccall((@blasfunc(zla_gbrcond_x_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Clong), trans, n, kl, ku, ab, ldab, afb,
                 ldafb, ipiv, x, info, work, rwork, 1)
end

function zla_gbrfsx_extended(prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb,
                             ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms,
                             err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond,
                             ithresh, rthresh, dz_ub, ignore_cwise, info)
    return ccall((@blasfunc(zla_gbrfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}), prec_type, trans_type, n, kl, ku, nrhs, ab, ldab, afb, ldafb,
                 ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm,
                 err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                 ignore_cwise, info)
end

function zla_gbrpvgrw(n, kl, ku, ncols, ab, ldab, afb, ldafb)
    return ccall((@blasfunc(zla_gbrpvgrw_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}), n, kl, ku, ncols, ab, ldab, afb, ldafb)
end

function zla_geamv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(zla_geamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}), trans,
                 m, n, alpha, a, lda, x, incx, beta, y, incy)
end

function zla_gercond_c(trans, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork)
    return ccall((@blasfunc(zla_gercond_c_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Clong), trans, n, a, lda, af, ldaf, ipiv, c, capply, info,
                 work, rwork, 1)
end

function zla_gercond_x(trans, n, a, lda, af, ldaf, ipiv, x, info, work, rwork)
    return ccall((@blasfunc(zla_gercond_x_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Clong), trans, n, a, lda, af, ldaf, ipiv, x, info, work,
                 rwork, 1)
end

function zla_gerfsx_extended(prec_type, trans_type, n, nrhs, a, lda, af, ldaf, ipiv, colequ,
                             c, b, ldb, y, ldy, berr_out, n_norms, errs_n, errs_c, res, ayb,
                             dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise, info)
    return ccall((@blasfunc(zla_gerfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}), prec_type, trans_type, n, nrhs, a,
                 lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms, errs_n,
                 errs_c, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise,
                 info)
end

function zla_gerpvgrw(n, ncols, a, lda, af, ldaf)
    return ccall((@blasfunc(zla_gerpvgrw_), libblastrampoline), Float64,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}), n, ncols, a, lda, af, ldaf)
end

function zla_heamv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(zla_heamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}), uplo,
                 n, alpha, a, lda, x, incx, beta, y, incy)
end

function zla_hercond_c(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork)
    return ccall((@blasfunc(zla_hercond_c_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Clong), uplo, n, a, lda, af, ldaf, ipiv, c, capply, info,
                 work, rwork, 1)
end

function zla_hercond_x(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork)
    return ccall((@blasfunc(zla_hercond_x_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Clong), uplo, n, a, lda, af, ldaf, ipiv, x, info, work,
                 rwork, 1)
end

function zla_herfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                             ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp,
                             res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                             ignore_cwise, info)
    return ccall((@blasfunc(zla_herfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), prec_type, uplo, n, nrhs, a,
                 lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms,
                 err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh,
                 rthresh, dz_ub, ignore_cwise, info, 1)
end

function zla_herpvgrw(uplo, n, info, a, lda, af, ldaf, ipiv, work)
    return ccall((@blasfunc(zla_herpvgrw_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Clong), uplo, n,
                 info, a, lda, af, ldaf, ipiv, work, 1)
end

function zla_lin_berr(n, nz, nrhs, res, ayb, berr)
    return ccall((@blasfunc(zla_lin_berr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64},
                  Ptr{Float64}), n, nz, nrhs, res, ayb, berr)
end

function zla_porcond_c(uplo, n, a, lda, af, ldaf, c, capply, info, work, rwork)
    return ccall((@blasfunc(zla_porcond_c_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Clong), uplo, n, a, lda, af, ldaf, c, capply, info, work,
                 rwork, 1)
end

function zla_porcond_x(uplo, n, a, lda, af, ldaf, x, info, work, rwork)
    return ccall((@blasfunc(zla_porcond_x_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64},
                  Clong), uplo, n, a, lda, af, ldaf, x, info, work, rwork, 1)
end

function zla_porfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, colequ, c, b, ldb,
                             y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp, res,
                             ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub, ignore_cwise,
                             info)
    return ccall((@blasfunc(zla_porfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), prec_type, uplo, n, nrhs, a,
                 lda, af, ldaf, colequ, c, b, ldb, y, ldy, berr_out, n_norms, err_bnds_norm,
                 err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                 ignore_cwise, info, 1)
end

function zla_porpvgrw(uplo, ncols, a, lda, af, ldaf, work)
    return ccall((@blasfunc(zla_porpvgrw_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Clong), uplo, ncols, a, lda, af, ldaf, work, 1)
end

function zla_syamv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(zla_syamv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt}), uplo,
                 n, alpha, a, lda, x, incx, beta, y, incy)
end

function zla_syrcond_c(uplo, n, a, lda, af, ldaf, ipiv, c, capply, info, work, rwork)
    return ccall((@blasfunc(zla_syrcond_c_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Clong), uplo, n, a, lda, af, ldaf, ipiv, c, capply, info,
                 work, rwork, 1)
end

function zla_syrcond_x(uplo, n, a, lda, af, ldaf, ipiv, x, info, work, rwork)
    return ccall((@blasfunc(zla_syrcond_x_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Clong), uplo, n, a, lda, af, ldaf, ipiv, x, info, work,
                 rwork, 1)
end

function zla_syrfsx_extended(prec_type, uplo, n, nrhs, a, lda, af, ldaf, ipiv, colequ, c, b,
                             ldb, y, ldy, berr_out, n_norms, err_bnds_norm, err_bnds_comp,
                             res, ayb, dy, y_tail, rcond, ithresh, rthresh, dz_ub,
                             ignore_cwise, info)
    return ccall((@blasfunc(zla_syrfsx_extended_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{Float64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong), prec_type, uplo, n, nrhs, a,
                 lda, af, ldaf, ipiv, colequ, c, b, ldb, y, ldy, berr_out, n_norms,
                 err_bnds_norm, err_bnds_comp, res, ayb, dy, y_tail, rcond, ithresh,
                 rthresh, dz_ub, ignore_cwise, info, 1)
end

function zla_syrpvgrw(uplo, n, info, a, lda, af, ldaf, ipiv, work)
    return ccall((@blasfunc(zla_syrpvgrw_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Clong), uplo, n,
                 info, a, lda, af, ldaf, ipiv, work, 1)
end

function zla_wwaddw(n, x, y, w)
    return ccall((@blasfunc(zla_wwaddw_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}), n, x, y, w)
end

function zlabrd(m, n, nb, a, lda, d, e, tauq, taup, x, ldx, y, ldy)
    return ccall((@blasfunc(zlabrd_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, nb, a, lda, d, e, tauq,
                 taup, x, ldx, y, ldy)
end

function zlacgv(n, x, incx)
    return ccall((@blasfunc(zlacgv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), n, x, incx)
end

function zlacn2(n, v, x, est, kase, isave)
    return ccall((@blasfunc(zlacn2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{Float64}, Ref{BlasInt},
                  Ptr{BlasInt}), n, v, x, est, kase, isave)
end

function zlacon(n, v, x, est, kase)
    return ccall((@blasfunc(zlacon_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{Float64}, Ref{BlasInt}), n,
                 v, x, est, kase)
end

function zlacp2(uplo, m, n, a, lda, b, ldb)
    return ccall((@blasfunc(zlacp2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, m, n, a, lda, b, ldb, 1)
end

function zlacpy(uplo, m, n, a, lda, b, ldb)
    return ccall((@blasfunc(zlacpy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, m, n, a, lda, b, ldb, 1)
end

function zlacrm(m, n, a, lda, b, ldb, c, ldc, rwork)
    return ccall((@blasfunc(zlacrm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}), m, n, a, lda, b, ldb, c, ldc,
                 rwork)
end

function zlacrt(n, cx, incx, cy, incy, c, s)
    return ccall((@blasfunc(zlacrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{ComplexF64}, Ref{ComplexF64}), n, cx, incx, cy, incy, c, s)
end

function zladiv(x, y)
    return ccall((@blasfunc(zladiv_), libblastrampoline), ComplexF64,
                 (Ref{ComplexF64}, Ref{ComplexF64}), x, y)
end

function zlaed0(qsiz, n, d, e, q, ldq, qstore, ldqs, rwork, iwork, info)
    return ccall((@blasfunc(zlaed0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ref{BlasInt}), qsiz, n, d, e, q, ldq, qstore, ldqs, rwork, iwork, info)
end

function zlaed7(n, cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, ldq, rho, indxq, qstore, qptr,
                prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, info)
    return ccall((@blasfunc(zlaed7_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), n,
                 cutpnt, qsiz, tlvls, curlvl, curpbm, d, q, ldq, rho, indxq, qstore, qptr,
                 prmptr, perm, givptr, givcol, givnum, work, rwork, iwork, info)
end

function zlaed8(k, n, qsiz, q, ldq, d, rho, cutpnt, z, dlambda, q2, ldq2, w, indxp, indx,
                indxq, perm, givptr, givcol, givnum, info)
    return ccall((@blasfunc(zlaed8_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}), k, n, qsiz, q, ldq, d,
                 rho, cutpnt, z, dlambda, q2, ldq2, w, indxp, indx, indxq, perm, givptr,
                 givcol, givnum, info)
end

function zlaein(rightv, noinit, n, h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info)
    return ccall((@blasfunc(zlaein_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}), rightv, noinit, n,
                 h, ldh, w, v, b, ldb, rwork, eps3, smlnum, info)
end

function zlaesy(a, b, c, rt1, rt2, evscal, cs1, sn1)
    return ccall((@blasfunc(zlaesy_), libblastrampoline), Cvoid,
                 (Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64},
                  Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}), a, b,
                 c, rt1, rt2, evscal, cs1, sn1)
end

function zlaev2(a, b, c, rt1, rt2, cs1, sn1)
    return ccall((@blasfunc(zlaev2_), libblastrampoline), Cvoid,
                 (Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{ComplexF64}), a, b, c, rt1, rt2, cs1, sn1)
end

function zlag2c(m, n, a, lda, sa, ldsa, info)
    return ccall((@blasfunc(zlag2c_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, sa, ldsa, info)
end

function zlags2(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq)
    return ccall((@blasfunc(zlags2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ref{ComplexF64}, Ref{Float64}, Ref{Float64},
                  Ref{ComplexF64}, Ref{Float64}, Ref{Float64}, Ref{ComplexF64},
                  Ref{Float64}, Ref{ComplexF64}, Ref{Float64}, Ref{ComplexF64}), upper, a1,
                 a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq)
end

function zlagtm(trans, n, nrhs, alpha, dl, d, du, x, ldx, beta, b, ldb)
    return ccall((@blasfunc(zlagtm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), trans, n, nrhs, alpha,
                 dl, d, du, x, ldx, beta, b, ldb, 1)
end

function zlahef(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(zlahef_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb,
                 a, lda, ipiv, w, ldw, info, 1)
end

function zlahef_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work)
    return ccall((@blasfunc(zlahef_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong), uplo, j1,
                 m, nb, a, lda, ipiv, h, ldh, work, 1)
end

function zlahef_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info)
    return ccall((@blasfunc(zlahef_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info, 1)
end

function zlahef_rook(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(zlahef_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb,
                 a, lda, ipiv, w, ldw, info, 1)
end

function zlahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info)
    return ccall((@blasfunc(zlahqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz,
                 z, ldz, info)
end

function zlahr2(n, k, nb, a, lda, tau, t, ldt, y, ldy)
    return ccall((@blasfunc(zlahr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}),
                 n, k, nb, a, lda, tau, t, ldt, y, ldy)
end

function zlaic1(job, j, x, sest, w, gamma, sestpr, s, c)
    return ccall((@blasfunc(zlaic1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{Float64}, Ptr{ComplexF64},
                  Ref{ComplexF64}, Ref{Float64}, Ref{ComplexF64}, Ref{ComplexF64}), job, j,
                 x, sest, w, gamma, sestpr, s, c)
end

function zlals0(icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx, perm, givptr, givcol, ldgcol,
                givnum, ldgnum, poles, difl, difr, z, k, c, s, rwork, info)
    return ccall((@blasfunc(zlals0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}), icompq, nl, nr, sqre, nrhs, b, ldb, bx, ldbx,
                 perm, givptr, givcol, ldgcol, givnum, ldgnum, poles, difl, difr, z, k, c,
                 s, rwork, info)
end

function zlalsa(icompq, smlsiz, n, nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z,
                poles, givptr, givcol, ldgcol, perm, givnum, c, s, rwork, iwork, info)
    return ccall((@blasfunc(zlalsa_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), icompq, smlsiz, n,
                 nrhs, b, ldb, bx, ldbx, u, ldu, vt, k, difl, difr, z, poles, givptr,
                 givcol, ldgcol, perm, givnum, c, s, rwork, iwork, info)
end

function zlalsd(uplo, smlsiz, n, nrhs, d, e, b, ldb, rcond, rank, work, rwork, iwork, info)
    return ccall((@blasfunc(zlalsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, smlsiz, n, nrhs, d, e,
                 b, ldb, rcond, rank, work, rwork, iwork, info, 1)
end

function zlamswlq(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zlamswlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork,
                 info, 1, 1)
end

function zlamtsqr(side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zlamtsqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong), side, trans, m, n, k, mb, nb, a, lda, t, ldt, c, ldc, work, lwork,
                 info, 1, 1)
end

function zlangb(norm, n, kl, ku, ab, ldab, work)
    return ccall((@blasfunc(zlangb_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Clong), norm, n, kl, ku, ab, ldab, work, 1)
end

function zlange(norm, m, n, a, lda, work)
    return ccall((@blasfunc(zlange_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Clong), norm, m, n, a, lda, work, 1)
end

function zlangt(norm, n, dl, d, du)
    return ccall((@blasfunc(zlangt_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Clong), norm, n, dl, d, du, 1)
end

function zlanhb(norm, uplo, n, k, ab, ldab, work)
    return ccall((@blasfunc(zlanhb_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Clong, Clong), norm, uplo, n, k, ab, ldab, work, 1, 1)
end

function zlanhe(norm, uplo, n, a, lda, work)
    return ccall((@blasfunc(zlanhe_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Clong, Clong), norm, uplo, n, a, lda, work, 1, 1)
end

function zlanhf(norm, transr, uplo, n, a, work)
    return ccall((@blasfunc(zlanhf_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Clong, Clong, Clong), norm, transr, uplo, n, a, work, 1, 1,
                 1)
end

function zlanhp(norm, uplo, n, ap, work)
    return ccall((@blasfunc(zlanhp_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Clong,
                  Clong), norm, uplo, n, ap, work, 1, 1)
end

function zlanhs(norm, n, a, lda, work)
    return ccall((@blasfunc(zlanhs_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Clong),
                 norm, n, a, lda, work, 1)
end

function zlanht(norm, n, d, e)
    return ccall((@blasfunc(zlanht_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Clong), norm, n, d,
                 e, 1)
end

function zlansb(norm, uplo, n, k, ab, ldab, work)
    return ccall((@blasfunc(zlansb_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Clong, Clong), norm, uplo, n, k, ab, ldab, work, 1, 1)
end

function zlansp(norm, uplo, n, ap, work)
    return ccall((@blasfunc(zlansp_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Clong,
                  Clong), norm, uplo, n, ap, work, 1, 1)
end

function zlansy(norm, uplo, n, a, lda, work)
    return ccall((@blasfunc(zlansy_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Clong, Clong), norm, uplo, n, a, lda, work, 1, 1)
end

function zlantb(norm, uplo, diag, n, k, ab, ldab, work)
    return ccall((@blasfunc(zlantb_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Clong, Clong, Clong), norm, uplo, diag, n, k, ab,
                 ldab, work, 1, 1, 1)
end

function zlantp(norm, uplo, diag, n, ap, work)
    return ccall((@blasfunc(zlantp_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Clong, Clong, Clong), norm, uplo, diag, n, ap, work, 1, 1,
                 1)
end

function zlantr(norm, uplo, diag, m, n, a, lda, work)
    return ccall((@blasfunc(zlantr_), libblastrampoline), Float64,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Clong, Clong, Clong), norm, uplo, diag, m, n, a,
                 lda, work, 1, 1, 1)
end

function zlapll(n, x, incx, y, incy, ssmin)
    return ccall((@blasfunc(zlapll_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}), n, x, incx, y, incy, ssmin)
end

function zlapmr(forwrd, m, n, x, ldx, k)
    return ccall((@blasfunc(zlapmr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}),
                 forwrd, m, n, x, ldx, k)
end

function zlapmt(forwrd, m, n, x, ldx, k)
    return ccall((@blasfunc(zlapmt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}),
                 forwrd, m, n, x, ldx, k)
end

function zlaqgb(m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax, equed)
    return ccall((@blasfunc(zlaqgb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
                  Ref{UInt8}, Clong), m, n, kl, ku, ab, ldab, r, c, rowcnd, colcnd, amax,
                 equed, 1)
end

function zlaqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed)
    return ccall((@blasfunc(zlaqge_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{UInt8},
                  Clong), m, n, a, lda, r, c, rowcnd, colcnd, amax, equed, 1)
end

function zlaqhb(uplo, n, kd, ab, ldab, s, scond, amax, equed)
    return ccall((@blasfunc(zlaqhb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{UInt8}, Clong, Clong), uplo,
                 n, kd, ab, ldab, s, scond, amax, equed, 1, 1)
end

function zlaqhe(uplo, n, a, lda, s, scond, amax, equed)
    return ccall((@blasfunc(zlaqhe_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{UInt8}, Clong, Clong), uplo, n, a, lda, s,
                 scond, amax, equed, 1, 1)
end

function zlaqhp(uplo, n, ap, s, scond, amax, equed)
    return ccall((@blasfunc(zlaqhp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{UInt8}, Clong, Clong), uplo, n, ap, s, scond, amax,
                 equed, 1, 1)
end

function zlaqp2(m, n, offset, a, lda, jpvt, tau, vn1, vn2, work)
    return ccall((@blasfunc(zlaqp2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}), m, n,
                 offset, a, lda, jpvt, tau, vn1, vn2, work)
end

function zlaqp2rk(m, n, nrhs, ioffset, kmax, abstol, reltol, kp1, maxc2nrm, a, lda, k,
                  maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2, work, info)
    return ccall((@blasfunc(zlaqp2rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, nrhs,
                 ioffset, kmax, abstol, reltol, kp1, maxc2nrm, a, lda, k, maxc2nrmk,
                 relmaxc2nrmk, jpiv, tau, vn1, vn2, work, info)
end

function zlaqp3rk(m, n, nrhs, ioffset, nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb,
                  maxc2nrmk, relmaxc2nrmk, jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info)
    return ccall((@blasfunc(zlaqp3rk_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}), m, n, nrhs, ioffset,
                 nb, abstol, reltol, kp1, maxc2nrm, a, lda, done, kb, maxc2nrmk,
                 relmaxc2nrmk, jpiv, tau, vn1, vn2, auxv, f, ldf, iwork, info)
end

function zlaqps(m, n, offset, nb, kb, a, lda, jpvt, tau, vn1, vn2, auxv, f, ldf)
    return ccall((@blasfunc(zlaqps_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, offset, nb, kb, a,
                 lda, jpvt, tau, vn1, vn2, auxv, f, ldf)
end

function zlaqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info)
    return ccall((@blasfunc(zlaqr0_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo,
                 ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info)
end

function zlaqr1(n, h, ldh, s1, s2, v)
    return ccall((@blasfunc(zlaqr1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}, Ref{ComplexF64},
                  Ptr{ComplexF64}), n, h, ldh, s1, s2, v)
end

function zlaqr2(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v,
                ldv, nh, t, ldt, nv, wv, ldwv, work, lwork)
    return ccall((@blasfunc(zlaqr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), wantt, wantz, n,
                 ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt,
                 nv, wv, ldwv, work, lwork)
end

function zlaqr3(wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v,
                ldv, nh, t, ldt, nv, wv, ldwv, work, lwork)
    return ccall((@blasfunc(zlaqr3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), wantt, wantz, n,
                 ktop, kbot, nw, h, ldh, iloz, ihiz, z, ldz, ns, nd, sh, v, ldv, nh, t, ldt,
                 nv, wv, ldwv, work, lwork)
end

function zlaqr4(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info)
    return ccall((@blasfunc(zlaqr4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), wantt, wantz, n, ilo,
                 ihi, h, ldh, w, iloz, ihiz, z, ldz, work, lwork, info)
end

function zlaqr5(wantt, wantz, kacc22, n, ktop, kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz,
                v, ldv, u, ldu, nv, wv, ldwv, nh, wh, ldwh)
    return ccall((@blasfunc(zlaqr5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), wantt, wantz, kacc22, n, ktop,
                 kbot, nshfts, s, h, ldh, iloz, ihiz, z, ldz, v, ldv, u, ldu, nv, wv, ldwv,
                 nh, wh, ldwh)
end

function zlaqsb(uplo, n, kd, ab, ldab, s, scond, amax, equed)
    return ccall((@blasfunc(zlaqsb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{UInt8}, Clong, Clong), uplo,
                 n, kd, ab, ldab, s, scond, amax, equed, 1, 1)
end

function zlaqsp(uplo, n, ap, s, scond, amax, equed)
    return ccall((@blasfunc(zlaqsp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{UInt8}, Clong, Clong), uplo, n, ap, s, scond, amax,
                 equed, 1, 1)
end

function zlaqsy(uplo, n, a, lda, s, scond, amax, equed)
    return ccall((@blasfunc(zlaqsy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{UInt8}, Clong, Clong), uplo, n, a, lda, s,
                 scond, amax, equed, 1, 1)
end

function zlar1v(n, b1, bn, lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz,
                mingma, r, isuppz, nrminv, resid, rqcorr, work)
    return ccall((@blasfunc(zlar1v_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64}), n, b1, bn,
                 lambda, d, l, ld, lld, pivmin, gaptol, z, wantnc, negcnt, ztz, mingma, r,
                 isuppz, nrminv, resid, rqcorr, work)
end

function zlar2v(n, x, y, z, incx, c, s, incc)
    return ccall((@blasfunc(zlar2v_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}), n, x, y, z, incx, c, s, incc)
end

function zlarcm(m, n, a, lda, b, ldb, c, ldc, rwork)
    return ccall((@blasfunc(zlarcm_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}), m, n, a, lda, b, ldb, c, ldc,
                 rwork)
end

function zlarf(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(zlarf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong),
                 side, m, n, v, incv, tau, c, ldc, work, 1)
end

function zlarf1f(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(zlarf1f_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong),
                 side, m, n, v, incv, tau, c, ldc, work, 1)
end

function zlarf1l(side, m, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(zlarf1l_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong),
                 side, m, n, v, incv, tau, c, ldc, work, 1)
end

function zlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork)
    return ccall((@blasfunc(zlarfb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong,
                  Clong, Clong), side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c,
                 ldc, work, ldwork, 1, 1, 1, 1)
end

function zlarfb_gett(ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork)
    return ccall((@blasfunc(zlarfb_gett_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong), ident, m, n, k, t, ldt, a, lda, b, ldb, work, ldwork,
                 1)
end

function zlarfg(n, alpha, x, incx, tau)
    return ccall((@blasfunc(zlarfg_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}),
                 n, alpha, x, incx, tau)
end

function zlarfgp(n, alpha, x, incx, tau)
    return ccall((@blasfunc(zlarfgp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}),
                 n, alpha, x, incx, tau)
end

function zlarft(direct, storev, n, k, v, ldv, tau, t, ldt)
    return ccall((@blasfunc(zlarft_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), direct,
                 storev, n, k, v, ldv, tau, t, ldt, 1, 1)
end

function zlarfx(side, m, n, v, tau, c, ldc, work)
    return ccall((@blasfunc(zlarfx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong), side, m, n, v, tau,
                 c, ldc, work, 1)
end

function zlarfy(uplo, n, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(zlarfy_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong), uplo, n, v, incv,
                 tau, c, ldc, work, 1)
end

function zlargv(n, x, incx, y, incy, c, incc)
    return ccall((@blasfunc(zlargv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}), n, x, incx, y, incy, c, incc)
end

function zlarnv(idist, iseed, n, x)
    return ccall((@blasfunc(zlarnv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}), idist, iseed, n, x)
end

function zlarrv(n, vl, vu, d, l, pivmin, isplit, m, dol, dou, minrgp, rtol1, rtol2, w, werr,
                wgap, iblock, indexw, gers, z, ldz, isuppz, work, iwork, info)
    return ccall((@blasfunc(zlarrv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{Float64}, Ptr{BlasInt}, Ref{BlasInt}), n, vl, vu, d, l, pivmin, isplit, m,
                 dol, dou, minrgp, rtol1, rtol2, w, werr, wgap, iblock, indexw, gers, z,
                 ldz, isuppz, work, iwork, info)
end

function zlarscl2(m, n, d, x, ldx)
    return ccall((@blasfunc(zlarscl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, d,
                 x, ldx)
end

function zlartv(n, x, incx, y, incy, c, s, incc)
    return ccall((@blasfunc(zlartv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}), n, x, incx, y, incy, c, s,
                 incc)
end

function zlarz(side, m, n, l, v, incv, tau, c, ldc, work)
    return ccall((@blasfunc(zlarz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong),
                 side, m, n, l, v, incv, tau, c, ldc, work, 1)
end

function zlarzb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, c, ldc, work,
                ldwork)
    return ccall((@blasfunc(zlarzb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong,
                  Clong, Clong, Clong), side, trans, direct, storev, m, n, k, l, v, ldv, t,
                 ldt, c, ldc, work, ldwork, 1, 1, 1, 1)
end

function zlarzt(direct, storev, n, k, v, ldv, tau, t, ldt)
    return ccall((@blasfunc(zlarzt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), direct,
                 storev, n, k, v, ldv, tau, t, ldt, 1, 1)
end

function zlascl(type, kl, ku, cfrom, cto, m, n, a, lda, info)
    return ccall((@blasfunc(zlascl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), type, kl, ku,
                 cfrom, cto, m, n, a, lda, info, 1)
end

function zlascl2(m, n, d, x, ldx)
    return ccall((@blasfunc(zlascl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, d,
                 x, ldx)
end

function zlaset(uplo, m, n, alpha, beta, a, lda)
    return ccall((@blasfunc(zlaset_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{ComplexF64}, Ref{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, m, n, alpha, beta, a, lda, 1)
end

function zlasr(side, pivot, direct, m, n, c, s, a, lda)
    return ccall((@blasfunc(zlasr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong), side,
                 pivot, direct, m, n, c, s, a, lda, 1, 1, 1)
end

function zlaswlq(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(zlaswlq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 mb, nb, a, lda, t, ldt, work, lwork, info)
end

function zlaswp(n, a, lda, k1, k2, ipiv, incx)
    return ccall((@blasfunc(zlaswp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, a, lda, k1, k2, ipiv, incx)
end

function zlasyf(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(zlasyf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb,
                 a, lda, ipiv, w, ldw, info, 1)
end

function zlasyf_aa(uplo, j1, m, nb, a, lda, ipiv, h, ldh, work)
    return ccall((@blasfunc(zlasyf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Clong), uplo, j1,
                 m, nb, a, lda, ipiv, h, ldh, work, 1)
end

function zlasyf_rk(uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info)
    return ccall((@blasfunc(zlasyf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nb, kb, a, lda, e, ipiv, w, ldw, info, 1)
end

function zlasyf_rook(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
    return ccall((@blasfunc(zlasyf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nb, kb,
                 a, lda, ipiv, w, ldw, info, 1)
end

function zlat2c(uplo, n, a, lda, sa, ldsa, info)
    return ccall((@blasfunc(zlat2c_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF32},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, sa, ldsa, info, 1)
end

function zlatbs(uplo, trans, diag, normin, n, kd, ab, ldab, x, scale, cnorm, info)
    return ccall((@blasfunc(zlatbs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong), uplo, trans, diag, normin, n, kd,
                 ab, ldab, x, scale, cnorm, info, 1, 1, 1, 1)
end

function zlatdf(ijob, n, z, ldz, rhs, rdsum, rdscal, ipiv, jpiv)
    return ccall((@blasfunc(zlatdf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{Float64}, Ref{Float64}, Ptr{BlasInt}, Ptr{BlasInt}), ijob, n, z, ldz, rhs,
                 rdsum, rdscal, ipiv, jpiv)
end

function zlatps(uplo, trans, diag, normin, n, ap, x, scale, cnorm, info)
    return ccall((@blasfunc(zlatps_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong), uplo, trans, diag, normin, n, ap, x, scale,
                 cnorm, info, 1, 1, 1, 1)
end

function zlatrd(uplo, n, nb, a, lda, e, tau, w, ldw)
    return ccall((@blasfunc(zlatrd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo,
                 n, nb, a, lda, e, tau, w, ldw, 1)
end

function zlatrs(uplo, trans, diag, normin, n, a, lda, x, scale, cnorm, info)
    return ccall((@blasfunc(zlatrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong, Clong, Clong), uplo, trans, diag, normin, n, a,
                 lda, x, scale, cnorm, info, 1, 1, 1, 1)
end

function zlatrs3(uplo, trans, diag, normin, n, nrhs, a, lda, x, ldx, scale, cnorm, work,
                 lwork, info)
    return ccall((@blasfunc(zlatrs3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong,
                  Clong), uplo, trans, diag, normin, n, nrhs, a, lda, x, ldx, scale, cnorm,
                 work, lwork, info, 1, 1, 1, 1)
end

function zlatrz(m, n, l, a, lda, tau, work)
    return ccall((@blasfunc(zlatrz_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}), m, n, l, a, lda, tau, work)
end

function zlatsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(zlatsqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 mb, nb, a, lda, t, ldt, work, lwork, info)
end

function zlaunhr_col_getrfnp(m, n, a, lda, d, info)
    return ccall((@blasfunc(zlaunhr_col_getrfnp_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}), m, n, a, lda, d, info)
end

function zlaunhr_col_getrfnp2(m, n, a, lda, d, info)
    return ccall((@blasfunc(zlaunhr_col_getrfnp2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}), m, n, a, lda, d, info)
end

function zlauu2(uplo, n, a, lda, info)
    return ccall((@blasfunc(zlauu2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function zlauum(uplo, n, a, lda, info)
    return ccall((@blasfunc(zlauum_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function zpbcon(uplo, n, kd, ab, ldab, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(zpbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, anorm, rcond, work, rwork, info, 1)
end

function zpbequ(uplo, n, kd, ab, ldab, s, scond, amax, info)
    return ccall((@blasfunc(zpbequ_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Clong), uplo, n, kd,
                 ab, ldab, s, scond, amax, info, 1)
end

function zpbrfs(uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(zpbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab, ldab, afb, ldafb, b, ldb, x,
                 ldx, ferr, berr, work, rwork, info, 1)
end

function zpbstf(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(zpbstf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function zpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(zpbsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab,
                 ldab, b, ldb, info, 1)
end

function zpbsvx(fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb, equed, s, b, ldb, x, ldx,
                rcond, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(zpbsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong, Clong), fact, uplo, n, kd, nrhs, ab, ldab, afb, ldafb,
                 equed, s, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1, 1, 1)
end

function zpbtf2(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(zpbtf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function zpbtrf(uplo, n, kd, ab, ldab, info)
    return ccall((@blasfunc(zpbtrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, kd, ab, ldab, info, 1)
end

function zpbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(zpbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, kd, nrhs, ab,
                 ldab, b, ldb, info, 1)
end

function zpftrf(transr, uplo, n, a, info)
    return ccall((@blasfunc(zpftrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong,
                  Clong), transr, uplo, n, a, info, 1, 1)
end

function zpftri(transr, uplo, n, a, info)
    return ccall((@blasfunc(zpftri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong,
                  Clong), transr, uplo, n, a, info, 1, 1)
end

function zpftrs(transr, uplo, n, nrhs, a, b, ldb, info)
    return ccall((@blasfunc(zpftrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), transr, uplo, n,
                 nrhs, a, b, ldb, info, 1, 1)
end

function zpocon(uplo, n, a, lda, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(zpocon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, anorm, rcond, work, rwork, info, 1)
end

function zpoequ(n, a, lda, s, scond, amax, info)
    return ccall((@blasfunc(zpoequ_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}), n, a, lda, s, scond, amax, info)
end

function zpoequb(n, a, lda, s, scond, amax, info)
    return ccall((@blasfunc(zpoequb_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}), n, a, lda, s, scond, amax, info)
end

function zporfs(uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr, berr, work, rwork,
                info)
    return ccall((@blasfunc(zporfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, af, ldaf, b, ldb, x, ldx, ferr,
                 berr, work, rwork, info, 1)
end

function zporfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, s, b, ldb, x, ldx, rcond, berr,
                 n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work, rwork,
                 info)
    return ccall((@blasfunc(zporfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), uplo, equed, n, nrhs, a, lda, af,
                 ldaf, s, b, ldb, x, ldx, rcond, berr, n_err_bnds, err_bnds_norm,
                 err_bnds_comp, nparams, params, work, rwork, info, 1, 1)
end

function zposv(uplo, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(zposv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, b,
                 ldb, info, 1)
end

function zposvx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                ferr, berr, work, rwork, info)
    return ccall((@blasfunc(zposvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong,
                  Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx,
                 rcond, ferr, berr, work, rwork, info, 1, 1, 1)
end

function zposvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx, rcond,
                 rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params,
                 work, rwork, info)
    return ccall((@blasfunc(zposvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{UInt8}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong,
                  Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, equed, s, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info, 1, 1, 1)
end

function zpotf2(uplo, n, a, lda, info)
    return ccall((@blasfunc(zpotf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function zpotrf(uplo, n, a, lda, info)
    return ccall((@blasfunc(zpotrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function zpotrf2(uplo, n, a, lda, info)
    return ccall((@blasfunc(zpotrf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function zpotri(uplo, n, a, lda, info)
    return ccall((@blasfunc(zpotri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, info, 1)
end

function zpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(zpotrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, b,
                 ldb, info, 1)
end

function zppcon(uplo, n, ap, anorm, rcond, work, rwork, info)
    return ccall((@blasfunc(zppcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{Float64}, Ref{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, ap, anorm,
                 rcond, work, rwork, info, 1)
end

function zppequ(uplo, n, ap, s, scond, amax, info)
    return ccall((@blasfunc(zppequ_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ref{Float64},
                  Ref{Float64}, Ref{BlasInt}, Clong), uplo, n, ap, s, scond, amax, info, 1)
end

function zpprfs(uplo, n, nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(zpprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n,
                 nrhs, ap, afp, b, ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function zppsv(uplo, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(zppsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, b, ldb, info, 1)
end

function zppsvx(fact, uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr,
                work, rwork, info)
    return ccall((@blasfunc(zppsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{UInt8}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong, Clong), fact,
                 uplo, n, nrhs, ap, afp, equed, s, b, ldb, x, ldx, rcond, ferr, berr, work,
                 rwork, info, 1, 1, 1)
end

function zpptrf(uplo, n, ap, info)
    return ccall((@blasfunc(zpptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, ap,
                 info, 1)
end

function zpptri(uplo, n, ap, info)
    return ccall((@blasfunc(zpptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, ap,
                 info, 1)
end

function zpptrs(uplo, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(zpptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, b, ldb, info, 1)
end

function zpstf2(uplo, n, a, lda, piv, rank, tol, work, info)
    return ccall((@blasfunc(zpstf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, a, lda, piv, rank,
                 tol, work, info, 1)
end

function zpstrf(uplo, n, a, lda, piv, rank, tol, work, info)
    return ccall((@blasfunc(zpstrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, a, lda, piv, rank,
                 tol, work, info, 1)
end

function zptcon(n, d, e, anorm, rcond, rwork, info)
    return ccall((@blasfunc(zptcon_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{Float64}, Ref{Float64},
                  Ptr{Float64}, Ref{BlasInt}), n, d, e, anorm, rcond, rwork, info)
end

function zpteqr(compz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(zpteqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work,
                 info, 1)
end

function zptrfs(uplo, n, nrhs, d, e, df, ef, b, ldb, x, ldx, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(zptrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, nrhs, d, e, df, ef, b, ldb, x,
                 ldx, ferr, berr, work, rwork, info, 1)
end

function zptsv(n, nrhs, d, e, b, ldb, info)
    return ccall((@blasfunc(zptsv_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}), n, nrhs, d, e, b, ldb, info)
end

function zptsvx(fact, n, nrhs, d, e, df, ef, b, ldb, x, ldx, rcond, ferr, berr, work, rwork,
                info)
    return ccall((@blasfunc(zptsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong), fact, n, nrhs, d, e, df,
                 ef, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1)
end

function zpttrf(n, d, e, info)
    return ccall((@blasfunc(zpttrf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}), n, d, e, info)
end

function zpttrs(uplo, n, nrhs, d, e, b, ldb, info)
    return ccall((@blasfunc(zpttrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, d, e, b,
                 ldb, info, 1)
end

function zptts2(iuplo, n, nrhs, d, e, b, ldb)
    return ccall((@blasfunc(zptts2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}), iuplo, n, nrhs, d, e, b, ldb)
end

function zrot(n, cx, incx, cy, incy, c, s)
    return ccall((@blasfunc(zrot_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{ComplexF64}), n, cx, incx, cy, incy, c, s)
end

function zrscl(n, a, x, incx)
    return ccall((@blasfunc(zrscl_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}), n, a, x, incx)
end

function zspcon(uplo, n, ap, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(zspcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, ap, ipiv,
                 anorm, rcond, work, info, 1)
end

function zspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    return ccall((@blasfunc(zspmv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 alpha, ap, x, incx, beta, y, incy, 1)
end

function zspr(uplo, n, alpha, x, incx, ap)
    return ccall((@blasfunc(zspr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Clong), uplo, n, alpha, x, incx, ap, 1)
end

function zsprfs(uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work, rwork, info)
    return ccall((@blasfunc(zsprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt},
                  Clong), uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, ferr, berr, work,
                 rwork, info, 1)
end

function zspsv(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(zspsv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b,
                 ldb, info, 1)
end

function zspsvx(fact, uplo, n, nrhs, ap, afp, ipiv, b, ldb, x, ldx, rcond, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(zspsvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), fact, uplo, n, nrhs, ap, afp,
                 ipiv, b, ldb, x, ldx, rcond, ferr, berr, work, rwork, info, 1, 1)
end

function zsptrf(uplo, n, ap, ipiv, info)
    return ccall((@blasfunc(zsptrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, ap, ipiv, info, 1)
end

function zsptri(uplo, n, ap, ipiv, work, info)
    return ccall((@blasfunc(zsptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong), uplo, n, ap, ipiv, work, info, 1)
end

function zsptrs(uplo, n, nrhs, ap, ipiv, b, ldb, info)
    return ccall((@blasfunc(zsptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, ap, ipiv, b,
                 ldb, info, 1)
end

function zstedc(compz, n, d, e, z, ldz, work, lwork, rwork, lrwork, iwork, liwork, info)
    return ccall((@blasfunc(zstedc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work, lwork, rwork,
                 lrwork, iwork, liwork, info, 1)
end

function zstegr(jobz, range, n, d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work,
                lwork, iwork, liwork, info)
    return ccall((@blasfunc(zstegr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n,
                 d, e, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork,
                 liwork, info, 1, 1)
end

function zstein(n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail, info)
    return ccall((@blasfunc(zstein_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}), n, d, e, m, w, iblock, isplit, z, ldz, work, iwork, ifail,
                 info)
end

function zstemr(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac,
                work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(zstemr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), jobz, range, n,
                 d, e, vl, vu, il, iu, m, w, z, ldz, nzc, isuppz, tryrac, work, lwork,
                 iwork, liwork, info, 1, 1)
end

function zsteqr(compz, n, d, e, z, ldz, work, info)
    return ccall((@blasfunc(zsteqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong), compz, n, d, e, z, ldz, work,
                 info, 1)
end

function zsycon(uplo, n, a, lda, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(zsycon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, ipiv, anorm, rcond, work, info, 1)
end

function zsycon_3(uplo, n, a, lda, e, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(zsycon_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong),
                 uplo, n, a, lda, e, ipiv, anorm, rcond, work, info, 1)
end

function zsycon_rook(uplo, n, a, lda, ipiv, anorm, rcond, work, info)
    return ccall((@blasfunc(zsycon_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, ipiv, anorm, rcond, work, info, 1)
end

function zsyconv(uplo, way, n, a, lda, ipiv, e, info)
    return ccall((@blasfunc(zsyconv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a, lda, ipiv, e,
                 info, 1, 1)
end

function zsyconvf(uplo, way, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(zsyconvf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a,
                 lda, e, ipiv, info, 1, 1)
end

function zsyconvf_rook(uplo, way, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(zsyconvf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), uplo, way, n, a,
                 lda, e, ipiv, info, 1, 1)
end

function zsyequb(uplo, n, a, lda, s, scond, amax, work, info)
    return ccall((@blasfunc(zsyequb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 a, lda, s, scond, amax, work, info, 1)
end

function zsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    return ccall((@blasfunc(zsymv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Clong), uplo, n, alpha, a, lda, x, incx, beta, y, incy, 1)
end

function zsyr(uplo, n, alpha, x, incx, a, lda)
    return ccall((@blasfunc(zsyr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, alpha, x, incx, a, lda, 1)
end

function zsyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(zsyrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{Float64}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, af, ldaf, ipiv, b,
                 ldb, x, ldx, ferr, berr, work, rwork, info, 1)
end

function zsyrfsx(uplo, equed, n, nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond,
                 berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams, params, work,
                 rwork, info)
    return ccall((@blasfunc(zsyrfsx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), uplo, equed, n,
                 nrhs, a, lda, af, ldaf, ipiv, s, b, ldb, x, ldx, rcond, berr, n_err_bnds,
                 err_bnds_norm, err_bnds_comp, nparams, params, work, rwork, info, 1, 1)
end

function zsysv(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zsysv_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function zsysv_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zsysv_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function zsysv_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork,
                         info)
    return ccall((@blasfunc(zsysv_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs,
                 a, lda, tb, ltb, ipiv, ipiv2, b, ldb, work, lwork, info, 1)
end

function zsysv_rk(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zsysv_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, e, ipiv, b, ldb,
                 work, lwork, info, 1)
end

function zsysv_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zsysv_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function zsysvx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr,
                berr, work, lwork, rwork, info)
    return ccall((@blasfunc(zsysvx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), fact,
                 uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, rcond, ferr, berr,
                 work, lwork, rwork, info, 1, 1)
end

function zsysvxx(fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b, ldb, x, ldx,
                 rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp, nparams,
                 params, work, rwork, info)
    return ccall((@blasfunc(zsysvxx_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{UInt8}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong, Clong), fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, equed, s, b,
                 ldb, x, ldx, rcond, rpvgrw, berr, n_err_bnds, err_bnds_norm, err_bnds_comp,
                 nparams, params, work, rwork, info, 1, 1, 1)
end

function zsyswapr(uplo, n, a, lda, i1, i2)
    return ccall((@blasfunc(zsyswapr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, i1, i2, 1)
end

function zsytf2(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(zsytf2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function zsytf2_rk(uplo, n, a, lda, e, ipiv, info)
    return ccall((@blasfunc(zsytf2_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, e, ipiv, info, 1)
end

function zsytf2_rook(uplo, n, a, lda, ipiv, info)
    return ccall((@blasfunc(zsytf2_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, ipiv, info, 1)
end

function zsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(zsytrf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function zsytrf_aa(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(zsytrf_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function zsytrf_aa_2stage(uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info)
    return ccall((@blasfunc(zsytrf_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong), uplo, n, a, lda, tb, ltb, ipiv, ipiv2, work, lwork, info, 1)
end

function zsytrf_rk(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(zsytrf_rk_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, lwork, info, 1)
end

function zsytrf_rook(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(zsytrf_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function zsytri(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(zsytri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function zsytri2(uplo, n, a, lda, ipiv, work, lwork, info)
    return ccall((@blasfunc(zsytri2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, lwork, info, 1)
end

function zsytri2x(uplo, n, a, lda, ipiv, work, nb, info)
    return ccall((@blasfunc(zsytri2x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv,
                 work, nb, info, 1)
end

function zsytri_3(uplo, n, a, lda, e, ipiv, work, lwork, info)
    return ccall((@blasfunc(zsytri_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, lwork, info, 1)
end

function zsytri_3x(uplo, n, a, lda, e, ipiv, work, nb, info)
    return ccall((@blasfunc(zsytri_3x_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda,
                 e, ipiv, work, nb, info, 1)
end

function zsytri_rook(uplo, n, a, lda, ipiv, work, info)
    return ccall((@blasfunc(zsytri_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, a, lda, ipiv, work, info, 1)
end

function zsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(zsytrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function zsytrs2(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, info)
    return ccall((@blasfunc(zsytrs2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n,
                 nrhs, a, lda, ipiv, b, ldb, work, info, 1)
end

function zsytrs_3(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info)
    return ccall((@blasfunc(zsytrs_3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info, 1)
end

function zsytrs_aa(uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info)
    return ccall((@blasfunc(zsytrs_aa_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong),
                 uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info, 1)
end

function zsytrs_aa_2stage(uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2, b, ldb, info)
    return ccall((@blasfunc(zsytrs_aa_2stage_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda, tb, ltb, ipiv, ipiv2,
                 b, ldb, info, 1)
end

function zsytrs_rook(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
    return ccall((@blasfunc(zsytrs_rook_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, nrhs, a, lda,
                 ipiv, b, ldb, info, 1)
end

function ztbcon(norm, uplo, diag, n, kd, ab, ldab, rcond, work, rwork, info)
    return ccall((@blasfunc(ztbcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong,
                  Clong, Clong), norm, uplo, diag, n, kd, ab, ldab, rcond, work, rwork,
                 info, 1, 1, 1)
end

function ztbrfs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, x, ldx, ferr, berr, work,
                rwork, info)
    return ccall((@blasfunc(ztbrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64},
                  Ref{BlasInt}, Clong, Clong, Clong), uplo, trans, diag, n, kd, nrhs, ab, ldab,
                 b, ldb, x, ldx, ferr, berr, work, rwork, info, 1, 1, 1)
end

function ztbtrs(uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info)
    return ccall((@blasfunc(ztbtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong,
                  Clong, Clong), uplo, trans, diag, n, kd, nrhs, ab, ldab, b, ldb, info, 1,
                 1, 1)
end

function ztfsm(transr, side, uplo, trans, diag, m, n, alpha, a, b, ldb)
    return ccall((@blasfunc(ztfsm_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong, Clong), transr, side, uplo, trans, diag, m, n,
                 alpha, a, b, ldb, 1, 1, 1, 1, 1)
end

function ztftri(transr, uplo, diag, n, a, info)
    return ccall((@blasfunc(ztftri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Clong, Clong, Clong), transr, uplo, diag, n, a, info, 1, 1, 1)
end

function ztfttp(transr, uplo, n, arf, ap, info)
    return ccall((@blasfunc(ztfttp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong, Clong), transr, uplo, n, arf, ap, info, 1, 1)
end

function ztfttr(transr, uplo, n, arf, a, lda, info)
    return ccall((@blasfunc(ztfttr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), transr, uplo, n, arf, a, lda, info,
                 1, 1)
end

function ztgevc(side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr, ldvr, mm, m, work,
                rwork, info)
    return ccall((@blasfunc(ztgevc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong), side, howmny, select, n, s, lds, p, ldp, vl, ldvl, vr,
                 ldvr, mm, m, work, rwork, info, 1, 1)
end

function ztgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, j1, info)
    return ccall((@blasfunc(ztgex2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}), wantq, wantz, n, a, lda, b, ldb, q, ldq,
                 z, ldz, j1, info)
end

function ztgexc(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, ifst, ilst, info)
    return ccall((@blasfunc(ztgexc_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}), wantq, wantz, n, a, lda, b,
                 ldb, q, ldq, z, ldz, ifst, ilst, info)
end

function ztgsen(ijob, wantq, wantz, select, n, a, lda, b, ldb, alpha, beta, q, ldq, z, ldz,
                m, pl, pr, dif, work, lwork, iwork, liwork, info)
    return ccall((@blasfunc(ztgsen_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}), ijob, wantq, wantz, select, n, a, lda,
                 b, ldb, alpha, beta, q, ldq, z, ldz, m, pl, pr, dif, work, lwork, iwork,
                 liwork, info)
end

function ztgsja(jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta, u,
                ldu, v, ldv, q, ldq, work, ncycle, info)
    return ccall((@blasfunc(ztgsja_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 jobu, jobv, jobq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, beta,
                 u, ldu, v, ldv, q, ldq, work, ncycle, info, 1, 1, 1)
end

function ztgsna(job, howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m,
                work, lwork, iwork, info)
    return ccall((@blasfunc(ztgsna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Clong, Clong), job,
                 howmny, select, n, a, lda, b, ldb, vl, ldvl, vr, ldvr, s, dif, mm, m, work,
                 lwork, iwork, info, 1, 1)
end

function ztgsy2(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                rdsum, rdscal, info)
    return ccall((@blasfunc(ztgsy2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Clong), trans, ijob,
                 m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, rdsum, rdscal,
                 info, 1)
end

function ztgsyl(trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale,
                dif, work, lwork, iwork, info)
    return ccall((@blasfunc(ztgsyl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong), trans, ijob, m, n, a, lda, b, ldb, c, ldc, d, ldd, e,
                 lde, f, ldf, scale, dif, work, lwork, iwork, info, 1)
end

function ztpcon(norm, uplo, diag, n, ap, rcond, work, rwork, info)
    return ccall((@blasfunc(ztpcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong,
                  Clong), norm, uplo, diag, n, ap, rcond, work, rwork, info, 1, 1, 1)
end

function ztplqt(m, n, l, mb, a, lda, b, ldb, t, ldt, work, info)
    return ccall((@blasfunc(ztplqt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}), m, n, l, mb, a, lda, b, ldb, t, ldt, work, info)
end

function ztplqt2(m, n, l, a, lda, b, ldb, t, ldt, info)
    return ccall((@blasfunc(ztplqt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 l, a, lda, b, ldb, t, ldt, info)
end

function ztpmlqt(side, trans, m, n, k, l, mb, v, ldv, t, ldt, a, lda, b, ldb, work, info)
    return ccall((@blasfunc(ztpmlqt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, l, mb, v, ldv, t, ldt, a,
                 lda, b, ldb, work, info, 1, 1)
end

function ztpmqrt(side, trans, m, n, k, l, nb, v, ldv, t, ldt, a, lda, b, ldb, work, info)
    return ccall((@blasfunc(ztpmqrt_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, l, nb, v, ldv, t, ldt, a,
                 lda, b, ldb, work, info, 1, 1)
end

function ztpqrt(m, n, l, nb, a, lda, b, ldb, t, ldt, work, info)
    return ccall((@blasfunc(ztpqrt_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}), m, n, l, nb, a, lda, b, ldb, t, ldt, work, info)
end

function ztpqrt2(m, n, l, a, lda, b, ldb, t, ldt, info)
    return ccall((@blasfunc(ztpqrt2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 l, a, lda, b, ldb, t, ldt, info)
end

function ztprfb(side, trans, direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb,
                work, ldwork)
    return ccall((@blasfunc(ztprfb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong, Clong, Clong), side, trans,
                 direct, storev, m, n, k, l, v, ldv, t, ldt, a, lda, b, ldb, work, ldwork,
                 1, 1, 1, 1)
end

function ztprfs(uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work, rwork,
                info)
    return ccall((@blasfunc(ztprfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong,
                  Clong), uplo, trans, diag, n, nrhs, ap, b, ldb, x, ldx, ferr, berr, work,
                 rwork, info, 1, 1, 1)
end

function ztptri(uplo, diag, n, ap, info)
    return ccall((@blasfunc(ztptri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong,
                  Clong), uplo, diag, n, ap, info, 1, 1)
end

function ztptrs(uplo, trans, diag, n, nrhs, ap, b, ldb, info)
    return ccall((@blasfunc(ztptrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), uplo, trans,
                 diag, n, nrhs, ap, b, ldb, info, 1, 1, 1)
end

function ztpttf(transr, uplo, n, ap, arf, info)
    return ccall((@blasfunc(ztpttf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong, Clong), transr, uplo, n, ap, arf, info, 1, 1)
end

function ztpttr(uplo, n, ap, a, lda, info)
    return ccall((@blasfunc(ztpttr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{BlasInt}, Clong), uplo, n, ap, a, lda, info, 1)
end

function ztrcon(norm, uplo, diag, n, a, lda, rcond, work, rwork, info)
    return ccall((@blasfunc(ztrcon_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong,
                  Clong), norm, uplo, diag, n, a, lda, rcond, work, rwork, info, 1, 1, 1)
end

function ztrevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, rwork,
                info)
    return ccall((@blasfunc(ztrevc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt}, Clong, Clong), side,
                 howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, rwork, info, 1,
                 1)
end

function ztrevc3(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m, work, lwork,
                 rwork, lrwork, info)
    return ccall((@blasfunc(ztrevc3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m,
                 work, lwork, rwork, lrwork, info, 1, 1)
end

function ztrexc(compq, n, t, ldt, q, ldq, ifst, ilst, info)
    return ccall((@blasfunc(ztrexc_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong), compq, n, t, ldt, q,
                 ldq, ifst, ilst, info, 1)
end

function ztrrfs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx, ferr, berr, work, rwork,
                info)
    return ccall((@blasfunc(ztrrfs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64}, Ptr{Float64}, Ref{BlasInt},
                  Clong, Clong, Clong), uplo, trans, diag, n, nrhs, a, lda, b, ldb, x, ldx,
                 ferr, berr, work, rwork, info, 1, 1, 1)
end

function ztrsen(job, compq, select, n, t, ldt, q, ldq, w, m, s, sep, work, lwork, info)
    return ccall((@blasfunc(ztrsen_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64},
                  Ref{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), job,
                 compq, select, n, t, ldt, q, ldq, w, m, s, sep, work, lwork, info, 1, 1)
end

function ztrsna(job, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, s, sep, mm, m, work,
                ldwork, rwork, info)
    return ccall((@blasfunc(ztrsna_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{Float64}, Ref{BlasInt}, Clong, Clong), job, howmny, select, n, t, ldt,
                 vl, ldvl, vr, ldvr, s, sep, mm, m, work, ldwork, rwork, info, 1, 1)
end

function ztrsyl(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, info)
    return ccall((@blasfunc(ztrsyl_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ref{BlasInt}, Clong, Clong), trana, tranb, isgn, m, n, a, lda,
                 b, ldb, c, ldc, scale, info, 1, 1)
end

function ztrsyl3(trana, tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, swork, ldswork,
                 info)
    return ccall((@blasfunc(ztrsyl3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ref{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), trana,
                 tranb, isgn, m, n, a, lda, b, ldb, c, ldc, scale, swork, ldswork, info, 1,
                 1)
end

function ztrti2(uplo, diag, n, a, lda, info)
    return ccall((@blasfunc(ztrti2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), uplo, diag, n, a, lda, info, 1, 1)
end

function ztrtri(uplo, diag, n, a, lda, info)
    return ccall((@blasfunc(ztrtri_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt},
                  Clong, Clong), uplo, diag, n, a, lda, info, 1, 1)
end

function ztrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    return ccall((@blasfunc(ztrtrs_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
                 uplo, trans, diag, n, nrhs, a, lda, b, ldb, info, 1, 1, 1)
end

function ztrttf(transr, uplo, n, a, lda, arf, info)
    return ccall((@blasfunc(ztrttf_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), transr, uplo, n, a, lda, arf,
                 info, 1, 1)
end

function ztrttp(uplo, n, a, lda, ap, info)
    return ccall((@blasfunc(ztrttp_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong), uplo, n, a, lda, ap, info, 1)
end

function ztzrzf(m, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(ztzrzf_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, a, lda, tau, work, lwork,
                 info)
end

function zunbdb(trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22,
                theta, phi, taup1, taup2, tauq1, tauq2, work, lwork, info)
    return ccall((@blasfunc(zunbdb_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), trans, signs, m, p, q, x11, ldx11,
                 x12, ldx12, x21, ldx21, x22, ldx22, theta, phi, taup1, taup2, tauq1, tauq2,
                 work, lwork, info, 1, 1)
end

function zunbdb1(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(zunbdb1_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}),
                 m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
end

function zunbdb2(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(zunbdb2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}),
                 m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
end

function zunbdb3(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
    return ccall((@blasfunc(zunbdb3_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}),
                 m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, work,
                 lwork, info)
end

function zunbdb4(m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1, taup2, tauq1, phantom,
                 work, lwork, info)
    return ccall((@blasfunc(zunbdb4_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}), m, p, q, x11, ldx11, x21, ldx21, theta, phi, taup1,
                 taup2, tauq1, phantom, work, lwork, info)
end

function zunbdb5(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
    return ccall((@blasfunc(zunbdb5_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m1, m2, n, x1, incx1,
                 x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
end

function zunbdb6(m1, m2, n, x1, incx1, x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
    return ccall((@blasfunc(zunbdb6_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m1, m2, n, x1, incx1,
                 x2, incx2, q1, ldq1, q2, ldq2, work, lwork, info)
end

function zuncsd(jobu1, jobu2, jobv1t, jobv2t, trans, signs, m, p, q, x11, ldx11, x12, ldx12,
                x21, ldx21, x22, ldx22, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t,
                work, lwork, rwork, lrwork, iwork, info)
    return ccall((@blasfunc(zuncsd_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{UInt8},
                  Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
                  Clong, Clong, Clong, Clong, Clong, Clong), jobu1, jobu2, jobv1t, jobv2t,
                 trans, signs, m, p, q, x11, ldx11, x12, ldx12, x21, ldx21, x22, ldx22,
                 theta, u1, ldu1, u2, ldu2, v1t, ldv1t, v2t, ldv2t, work, lwork, rwork,
                 lrwork, iwork, info, 1, 1, 1, 1, 1, 1)
end

function zuncsd2by1(jobu1, jobu2, jobv1t, m, p, q, x11, ldx11, x21, ldx21, theta, u1, ldu1,
                    u2, ldu2, v1t, ldv1t, work, lwork, rwork, lrwork, iwork, info)
    return ccall((@blasfunc(zuncsd2by1_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
                  Ref{BlasInt}, Clong, Clong, Clong), jobu1, jobu2, jobv1t, m, p, q, x11,
                 ldx11, x21, ldx21, theta, u1, ldu1, u2, ldu2, v1t, ldv1t, work, lwork,
                 rwork, lrwork, iwork, info, 1, 1, 1)
end

function zung2l(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(zung2l_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, k, a, lda, tau, work,
                 info)
end

function zung2r(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(zung2r_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, k, a, lda, tau, work,
                 info)
end

function zungbr(vect, m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zungbr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), vect, m,
                 n, k, a, lda, tau, work, lwork, info, 1)
end

function zunghr(n, ilo, ihi, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zunghr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), n, ilo, ihi, a,
                 lda, tau, work, lwork, info)
end

function zungl2(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(zungl2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, k, a, lda, tau, work,
                 info)
end

function zunglq(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zunglq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda,
                 tau, work, lwork, info)
end

function zungql(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zungql_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda,
                 tau, work, lwork, info)
end

function zungqr(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zungqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda,
                 tau, work, lwork, info)
end

function zungr2(m, n, k, a, lda, tau, work, info)
    return ccall((@blasfunc(zungr2_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, k, a, lda, tau, work,
                 info)
end

function zungrq(m, n, k, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zungrq_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n, k, a, lda,
                 tau, work, lwork, info)
end

function zungtr(uplo, n, a, lda, tau, work, lwork, info)
    return ccall((@blasfunc(zungtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong), uplo, n, a, lda, tau, work,
                 lwork, info, 1)
end

function zungtsqr(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(zungtsqr_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 mb, nb, a, lda, t, ldt, work, lwork, info)
end

function zungtsqr_row(m, n, mb, nb, a, lda, t, ldt, work, lwork, info)
    return ccall((@blasfunc(zungtsqr_row_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}), m, n,
                 mb, nb, a, lda, t, ldt, work, lwork, info)
end

function zunhr_col(m, n, nb, a, lda, t, ldt, d, info)
    return ccall((@blasfunc(zunhr_col_), libblastrampoline), Cvoid,
                 (Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}), m, n, nb, a, lda,
                 t, ldt, d, info)
end

function zunm22(side, trans, m, n, n1, n2, q, ldq, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zunm22_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, n1, n2, q, ldq, c,
                 ldc, work, lwork, info, 1, 1)
end

function zunm2l(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(zunm2l_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 info, 1, 1)
end

function zunm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(zunm2r_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 info, 1, 1)
end

function zunmbr(vect, side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zunmbr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), vect, side,
                 trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info, 1, 1, 1)
end

function zunmhr(side, trans, m, n, ilo, ihi, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zunmhr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 ilo, ihi, a, lda, tau, c, ldc, work, lwork, info, 1, 1)
end

function zunml2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(zunml2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 info, 1, 1)
end

function zunmlq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zunmlq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1)
end

function zunmql(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zunmql_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1)
end

function zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zunmqr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1)
end

function zunmr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(zunmr2_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c, ldc, work,
                 info, 1, 1)
end

function zunmr3(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, info)
    return ccall((@blasfunc(zunmr3_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, l, a,
                 lda, tau, c, ldc, work, info, 1, 1)
end

function zunmrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zunmrq_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n, k, a, lda, tau, c,
                 ldc, work, lwork, info, 1, 1)
end

function zunmrz(side, trans, m, n, k, l, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zunmrz_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt},
                  Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong), side, trans, m, n,
                 k, l, a, lda, tau, c, ldc, work, lwork, info, 1, 1)
end

function zunmtr(side, uplo, trans, m, n, a, lda, tau, c, ldc, work, lwork, info)
    return ccall((@blasfunc(zunmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong), side, uplo, trans, m, n, a,
                 lda, tau, c, ldc, work, lwork, info, 1, 1, 1)
end

function zupgtr(uplo, n, ap, tau, q, ldq, work, info)
    return ccall((@blasfunc(zupgtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64}, Ptr{ComplexF64},
                  Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Clong), uplo, n, ap, tau, q, ldq,
                 work, info, 1)
end

function zupmtr(side, uplo, trans, m, n, ap, tau, c, ldc, work, info)
    return ccall((@blasfunc(zupmtr_), libblastrampoline), Cvoid,
                 (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ptr{ComplexF64},
                  Ptr{ComplexF64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
                  Clong, Clong, Clong), side, uplo, trans, m, n, ap, tau, c, ldc, work,
                 info, 1, 1, 1)
end
