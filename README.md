# LinearAlgebra

This package is part of the Julia standard library (stdlib).

`LinearAlgebra.jl` provides functionality for performing linear algebra operations in Julia.

| **Build Status**                    | **Coverage**                    |
|:-------------------------------:|:-------------------------------:|
[![Build status](https://badge.buildkite.com/4032f509b3a7354c8ce7b34b98c7496d66adc4504b0101cbcf.svg?branch=master)](https://buildkite.com/julialang/linearalgebra-dot-jl) | [![][codecov-img]][codecov-url] |

[codecov-img]: https://codecov.io/gh/JuliaLang/LinearAlgebra.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaLang/LinearAlgebra.jl

## Migrating a Pull Request (PR) from the Julia repository to this repository

Since this package was split out from the main Julia repository, you might have previously made a pull request (PR) to the Julia repo. You can easily migrate such PRs to this repository using the following steps:

1. Add the Julia repository (or your fork) as a new remote repository:
   ```bash
   git remote add juliarepo https://github.com/JuliaLang/julia
   ```

2. Fetch the commits from the Julia repository:
   ```bash
   git fetch juliarepo
   ```

3. Cherry-pick the relevant commits made in the Julia repository to this repository:
   ```bash
   git cherry-pick $JULIA_COMMIT
   ```

## Using development versions of this package

This package performs some type piracy and is also included in the sysimage, which makes using a development version slightly more complex than usual.

To use a development version of this package, you can choose one of the following methods:

1. **Change the UUID in the project file and load the package:**

   Change the UUID line in `Project.toml` as
   ```diff
   - uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
   + uuid = "27e2e46d-f89d-539d-b4ee-838fcccc9c8e"
   ```

   Start `julia` as
   ```console
   JULIA_PRUNE_OLD_LA=true julia +nightly --compiled-modules=existing --project
   ```
   where it is assumed that one is already within the `LinearAlgebra` directory (otherwise, adjust
   the project path accordingly). The `julia +nightly` command above assumes that `juliaup` is being used
   to launch `julia`, but one may substitute this with the path to the julia executable.

   Within the `julia` session, load the `LinearAlgebra` module after pruning the one in the sysimage. This may be done as
   ```julia
   include("test/prune_old_LA.jl") && using LinearAlgebra
   ```

   Note that loading the test files in the REPL will automatically carry out the pruning to ensure that the development version of the package is used in the tests.

   If you are contributing to the repo using this method, it may be convenient to ignore the local changes to `Project.toml` by running
   ```console
   git update-index --skip-worktree Project.toml
   ```

2. **Build Julia with the custom `LinearAlgebra` commit:**

   Modify the commit in [`stdlib/LinearAlgebra.version`](https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra.version) and build Julia. This requires one to push the development branch
   to `GitHub` or an equivalent platform.

3. **Build a custom sysimage with the new `LinearAlgebra`:**
   - Install `PackageCompiler`.
   - Load it and, with this project active, run:
     ```julia
     create_sysimage(["LinearAlgebra"]; sysimage_path="new_sysimage.so", incremental=false, filter_stdlibs=true)
     ```
   - Start Julia with the custom sysimage:
     ```bash
     julia -Jnew_sysimage.so
     ```
