module LinearAlgebraTestHelpers

export COMPREHENSIVE, DualNumbers, FillArrays, Furlongs, ImmutableArrays,
       InfiniteArrays, OffsetArrays, Quaternions, SizedArrays, StructArrays

# Comprehensive mode additionally runs the tests guarded with `@static if COMPREHENSIVE`
# (more element types, wrappers, sizes and corner cases). It is enabled by
# `LINALG_TEST_COMPREHENSIVE=true`, which `runtests.jl` sets for `--comprehensive`.
const COMPREHENSIVE = get(ENV, "LINALG_TEST_COMPREHENSIVE", "false") == "true"

include("DualNumbers.jl")
include("FillArrays.jl")
include("Furlongs.jl")
include("ImmutableArrays.jl")
include("InfiniteArrays.jl")
include("OffsetArrays.jl")
include("Quaternions.jl")
include("SizedArrays.jl")
include("StridedArrays.jl")
include("StructArrays.jl")

end
