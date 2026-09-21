# One-time setup for the Julia benchmarks.
#
# Project.toml already points BracketingNonlinearSolve / NonlinearSolveBase at the
# local NonlinearSolve.jl checkout (../../NonlinearSolve.jl) via [sources], so this
# just resolves and precompiles that environment.
#
#   julia --project=. setup_env.jl
#
# Afterwards, run the benchmarks with:
#   julia --project=. nonlinearsolve_modab_benchmark.jl
#   julia --project=. roots_modab_benchmark.jl
#   julia --project=. "NonlinearSolve.jl/nonlinearsolve_modab__speed_benchmark.jl"
using Pkg

Pkg.activate(@__DIR__)

const NLS_ROOT = normpath(joinpath(@__DIR__, "..", "..", "NonlinearSolve.jl"))
isdir(NLS_ROOT) || error("Local NonlinearSolve.jl checkout not found at $NLS_ROOT")

Pkg.resolve()
Pkg.instantiate()
Pkg.precompile()
Pkg.status()
