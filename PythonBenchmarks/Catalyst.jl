using Pkg

# Activate a local environment next to this script and only install packages
# if they are missing from that environment.
Pkg.activate(joinpath(@__DIR__, "catalyst_environment"))

function ensure_packages(pkgs)
	installed = Set(dep.name for dep in values(Pkg.dependencies()))
	missing = [pkg for pkg in pkgs if !(pkg in installed)]
	if !isempty(missing)
		Pkg.add(missing)
	end
end

ensure_packages([
	"Catalyst",
	"DifferentialEquations",
	"JumpProcesses",
])

using Catalyst
using DifferentialEquations
using JumpProcesses
using Random
using Statistics

# -----------------------------
# Basic Julia notes
# -----------------------------
# 1. Julia arrays are 1-based.
# 2. Dictionaries use =>, for example Dict(:X => 1).
# 3. Symbols like :X are commonly used for species/parameter names.
# 4. Functions are defined with the function ... end syntax.

# -----------------------------
# Helper utilities
# -----------------------------
function make_varmap(vars)
	varmap = Dict{Symbol, Any}()
	for var in vars
		var_str = string(var)
		varmap[Symbol(var_str)] = var
		base_str = replace(var_str, r"\(.*\)" => "")
		varmap[Symbol(base_str)] = var
	end
	return varmap
end

function normalize_pairs(varmap, pairs)
	for (name, _) in pairs
		if !haskey(varmap, name)
			error("Unknown Catalyst symbol: $(name)")
		end
	end
	return [varmap[name] => value for (name, value) in pairs]
end

function solve_ssa_ensemble(rn, u0, p, tspan; ntraj=100, seed=1)
	Random.seed!(seed)
	rn_complete = complete(rn)
	species_map = make_varmap(species(rn_complete))
	param_map = make_varmap(parameters(rn_complete))
	u0_norm = normalize_pairs(species_map, u0)
	p_norm = normalize_pairs(param_map, p)
	dprob = DiscreteProblem(rn_complete, u0_norm, tspan, p_norm)
	jprob = JumpProblem(rn_complete, dprob, Direct())
	eprob = EnsembleProblem(jprob)
	return solve(eprob, SSAStepper(), EnsembleSerial(); trajectories=ntraj)
end

function final_time_mean(ensemble_sol, species_names)
	means = Dict{Symbol, Float64}()
	for sp in species_names
		means[sp] = mean(sol[sp][end] for sol in ensemble_sol)
	end
	return means
end

function benchmark_ssa(rn, u0, p, tspan; ntraj=1000, seed=1)
	elapsed = @elapsed begin
		ensemble_sol = solve_ssa_ensemble(rn, u0, p, tspan; ntraj=ntraj, seed=seed)
		println("Final-time means: ", final_time_mean(ensemble_sol, [first(kv) for kv in u0]))
	end
	println("Elapsed time for ", ntraj, " trajectories: ", elapsed, " seconds")
end

# -----------------------------
# Example 1: Birth-death model
# -----------------------------
# SSIT version:
#   birth: bk
#   death: dk*X
birth_death_rn = @reaction_network begin
	bk, 0 --> X
	dk, X --> 0
end

u0_bd = [:X => 0]
p_bd = [:bk => 1000.0, :dk => 1.0]
tspan_bd = (0.0, 10.0)

bd_sol = solve_ssa_ensemble(birth_death_rn, u0_bd, p_bd, tspan_bd; ntraj=100)
println("Birth-death final mean: ", final_time_mean(bd_sol, [:X]))

# -----------------------------
# Example 2: Goutsias model
# -----------------------------
# This mirrors the SSIT benchmark model using mass-action reactions.
goutsias_rn = @reaction_network begin
	k1, RNA --> RNA + M
	k2, M --> 0
	k3, DNAD --> DNAD + RNA
	k4, RNA --> 0
	k5, DNA + D --> DNAD
	k6, DNAD --> DNA + D
	k7, DNAD + D --> DNA2D
	k8, DNA2D --> DNAD + D
	k9, 2M --> D
	k10, D --> 2M
end

u0_g = [
	:DNA => 2,
	:DNAD => 0,
	:DNA2D => 0,
	:RNA => 0,
	:M => 2,
	:D => 6,
]

p_g = [
	:k1 => 0.043,
	:k2 => 0.0007,
	:k3 => 0.0715,
	:k4 => 0.0039,
	:k5 => 0.012e9 / (6.0221415e23 * 1e-15),
	:k6 => 0.4791,
	:k7 => 0.00012e9 / (6.0221415e23 * 1e-15),
	:k8 => 0.8765e-11,
	# Catalyst mass-action for 2M -> D already uses the correct combinatorics,
	# so keep k9 as the underlying kinetic constant, not k9/2.
	:k9 => 0.05e9 / (6.0221415e23 * 1e-15),
	:k10 => 0.5,
]

tspan_g = (0.0, 100.0)

goutsias_sol = solve_ssa_ensemble(goutsias_rn, u0_g, p_g, tspan_g; ntraj=100)
println("Goutsias final means: ", final_time_mean(goutsias_sol, [:DNA, :DNAD, :DNA2D, :RNA, :M, :D]))

# -----------------------------
# Getting started workflow
# -----------------------------
# 1. Run this file in Julia.
# 2. Confirm the birth-death and Goutsias examples execute.
# 3. Change ntraj in solve_ssa_ensemble or benchmark_ssa.
# 4. Add a new @reaction_network for another SSIT benchmark model.
#
# Example benchmark call:
# benchmark_ssa(goutsias_rn, u0_g, p_g, tspan_g; ntraj=1000, seed=1)

