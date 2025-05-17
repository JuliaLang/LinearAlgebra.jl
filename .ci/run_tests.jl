if haskey(ENV, "BUILDKITE")
    ncores = Sys.CPU_THREADS
else
    ncores = ceil(Int, Sys.CPU_THREADS / 2)
end

proj = abspath(joinpath(@__DIR__, ".."))
cmd = """Base.runtests(["LinearAlgebra"]; propagate_project=true, ncores=$ncores)"""
run(addenv(`$(Base.julia_cmd()) --project=$proj --compiled-modules=existing -e $cmd`,
    "JULIA_NUM_THREADS" => 1, "JULIA_PRUNE_OLD_LA" => true))
