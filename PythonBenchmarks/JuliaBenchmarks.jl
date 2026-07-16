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

struct ModelSpec
	jump_prob::Any
	saveat::Vector{Float64}
	tspan::Tuple{Float64, Float64}
	output_species::Vector{Symbol}
	species_index::Dict{Symbol, Int}
end

function prepare_model_spec(rn, u0, p, tspan, output_species; save_points)
	rn_complete = complete(rn)
	species_vars = species(rn_complete)
	species_map = make_varmap(species_vars)
	param_map = make_varmap(parameters(rn_complete))
	u0_norm = Dict(normalize_pairs(species_map, u0))
	p_norm = Dict(normalize_pairs(param_map, p))
	species_order = [Symbol(replace(string(sp), r"\(.*\)" => "")) for sp in species_vars]
	species_index = Dict(sp => idx for (idx, sp) in enumerate(species_order))
	jump_prob = JumpProblem(rn_complete, u0_norm, tspan, p_norm; save_positions=(false, false))
	saveat = collect(range(tspan[1], tspan[2]; length=save_points))
	return ModelSpec(jump_prob, saveat, tspan, output_species, species_index)
end

normalize_model_name(model_name) = lowercase(replace(String(model_name), r"[\s_-]+" => ""))

const MODEL_ALIASES = Dict(
	"birthdeath" => :birth_death,
	"bd" => :birth_death,
	"goutsias" => :goutsias,
	"goutsias100" => :goutsias100,
	"goutsias300" => :goutsias300,
	"goutsias1000" => :goutsias1000,
	"dnaexpression" => :goutsias,
	"mm" => :michaelis_menten,
	"michaelismenten" => :michaelis_menten,
	"mapk" => :mapk,
	"toggle" => :toggle,
	"6speciestoggle" => :toggle,
	"phage" => :phage,
)

function build_birth_death_model()
	rn = @reaction_network begin
		bk, 0 --> X
		dk, X --> 0
	end

	return prepare_model_spec(
		rn,
		[:X => 0],
		[:bk => 1000.0, :dk => 1.0],
		(0.0, 10.0),
		[:X],
		save_points=11,
	)
end

function build_goutsias_model(tfinal::Float64=100.0)
	rn = @reaction_network begin
		k1, RNA --> RNA + M
		k2, M --> 0
		k3, DNAD --> DNAD + RNA
		k4, RNA --> 0
		k5, DNA + D --> DNAD
		k6, DNAD --> DNA + D
		k7, DNAD + D --> DNA2D
		k8, DNA2D --> DNAD + D
		k9, M + M --> D
		k10, D --> M + M
	end

	return prepare_model_spec(
		rn,
		[
			:DNA => 2,
			:DNAD => 0,
			:DNA2D => 0,
			:RNA => 0,
			:M => 2,
			:D => 6,
		],
		[
			:k1 => 0.043,
			:k2 => 0.0007,
			:k3 => 0.0715,
			:k4 => 0.0039,
			:k5 => 0.012e9 / (6.0221415e23 * 1e-15),
			:k6 => 0.4791,
			:k7 => 0.00012e9 / (6.0221415e23 * 1e-15),
			:k8 => 0.8765e-11,
			:k9 => 0.05e9 / (6.0221415e23 * 1e-15),
			:k10 => 0.5,
		],
			(0.0, tfinal),
		[:DNA, :DNAD, :DNA2D, :RNA, :M, :D],
		save_points=11,
	)
end

	const Goutsias100 = build_goutsias_model(100.0)
	const Goutsias300 = build_goutsias_model(300.0)
	const Goutsias1000 = build_goutsias_model(1000.0)

function build_michaelis_menten_model()
	rn = @reaction_network begin
		k1, S + E --> ES
		k2, ES --> S + E
		k3, ES --> E + P
	end

	return prepare_model_spec(
		rn,
		[:S => 100, :E => 1000, :ES => 0, :P => 0],
		[:k1 => 1.0, :k2 => 1.0, :k3 => 0.1],
		(0.0, 70.0),
		[:S, :E, :ES, :P],
		save_points=71,
	)
end

function build_mapk_model()
	rn = @reaction_network begin
		s2, 0 --> MEK
		d2, MEK --> 0
		s3, Kpp --> Kpp + MEK
		s1, 0 --> K
		d1, K --> 0
		d1, KpY --> 0
		d1, KpT --> 0
		d1, Kpp --> 0
		k1, K + MEK --> K_MEK_Y
		k2, K_MEK_Y --> KpY + MEK
		k3, KpY + MEK --> KpY_MEK
		k_3, KpY_MEK --> KpY + MEK
		k4, KpY_MEK --> Kpp + MEK
		k5, K + MEK --> K_MEK_T
		k_5, K_MEK_T --> K + MEK
		k6, K_MEK_T --> KpT + MEK
		k7, KpT + MEK --> KpT_MEK
		k_7, KpT_MEK --> KpT + MEK
		k8, KpT_MEK --> Kpp + MEK
		h1, Kpp + MKP3 --> Kpp_MKP3
		h_1, Kpp_MKP3 --> Kpp + MKP3
		h2, Kpp_MKP3 --> KpT_MKP3_Y
		h3, KpT_MKP3_Y --> KpT + MKP3
		h_3, KpT + MKP3 --> KpT_MKP3_Y
		h4, KpT + MKP3 --> KpT_MKP3_T
		h_4, KpT_MKP3_T --> KpT + MKP3
		h5, KpT_MKP3_T --> K_MKP3_T
		h6, K_MKP3_T --> K + MKP3
		h_6, K + MKP3 --> K_MKP3_T
		h7, KpY + MKP3 --> KpY_MKP3
		h_7, KpY_MKP3 --> KpY + MKP3
		h8, KpY_MKP3 --> K_MKP3_Y
		h9, K_MKP3_Y --> K + MKP3
		h_9, K + MKP3 --> K_MKP3_Y
	end

	return prepare_model_spec(
		rn,
		[
			:MKP3 => 1,
			:Kpp_MKP3 => 0,
			:KpT_MKP3_Y => 0,
			:KpT_MKP3_T => 0,
			:K_MKP3_T => 0,
			:KpY_MKP3 => 0,
			:K_MKP3_Y => 0,
			:KpT_MEK => 0,
			:KpT => 0,
			:KpY_MEK => 0,
			:K_MEK_T => 0,
			:KpY => 0,
			:K_MEK_Y => 0,
			:K => 3,
			:MEK => 0,
			:Kpp => 0,
		],
		[
			:s1 => 0.00024,
			:d1 => 0.0001,
			:s2 => 0.001,
			:d2 => 0.15,
			:s3 => 0.005,
			:k1 => 0.375,
			:k2 => 0.06,
			:k3 => 0.375,
			:k_3 => 1.0,
			:k4 => 4.5,
			:k5 => 0.375,
			:k_5 => 1.0,
			:k6 => 0.06,
			:k7 => 0.375,
			:k_7 => 1.0,
			:k8 => 4.5,
			:h1 => 0.015,
			:h_1 => 1.0,
			:h2 => 0.032,
			:h3 => 0.31,
			:h_3 => 0.01,
			:h4 => 0.01,
			:h_4 => 1.0,
			:h5 => 0.5,
			:h6 => 0.086,
			:h_6 => 0.0011,
			:h7 => 0.01,
			:h_7 => 1.0,
			:h8 => 0.47,
			:h9 => 0.14,
			:h_9 => 0.0018,
		],
		(0.0, 10.0),
		[:MKP3, :Kpp_MKP3, :KpT_MKP3_Y, :KpT_MKP3_T, :K_MKP3_T, :KpY_MKP3, :K_MKP3_Y, :KpT_MEK, :KpT, :KpY_MEK, :K_MEK_T, :KpY, :K_MEK_Y, :K, :MEK, :Kpp],
		save_points=101,
	)
end

function build_toggle_model()
	rn = @reaction_network begin
		k1, GeneA --> GeneA + A
		k2, GeneB --> GeneB + B
		k3, A --> 0
		k4, B --> 0
		k5, A + A + GeneB --> bGeneB
		k7, bGeneB --> A + A + GeneB
		k6, B + B + GeneA --> bGeneA
		k8, bGeneA --> B + B + GeneA
	end

	return prepare_model_spec(
		rn,
		[
			:GeneA => 1,
			:GeneB => 1,
			:A => 0,
			:B => 0,
			:bGeneB => 0,
			:bGeneA => 0,
		],
		[
			:k1 => 40.0,
			:k2 => 20.0,
			:k3 => 1.0,
			:k4 => 1.0,
			:k5 => 1e-5 / 2,
			:k6 => 3.5e-5 / 2,
			:k7 => 1.0,
			:k8 => 1.0,
		],
		(0.0, 30.0),
		[:GeneA, :GeneB, :A, :B, :bGeneB, :bGeneA],
		save_points=11,
	)
end

function build_phage_model()
	rn = @reaction_network begin
		k1, OR3 + OR2 --> OR3 + OR2 + CI2
		k2, OR3 + COR2 --> OR3 + COR2 + CI2
		k3, OR3 + ROR2 --> OR3 + ROR2 + CI2
		k4, OR1 + OR2 --> OR1 + OR2 + Cro2
		k5, CI2 --> 0
		k6, Cro2 --> 0
		k7, CI2 + OR1 --> ROR1
		k8, CI2 + OR2 --> ROR2
		k9, CI2 + OR3 --> ROR3
		k10, Cro2 + OR1 --> COR1
		k11, Cro2 + OR2 --> COR2
		k12, Cro2 + OR3 --> COR3
		k13, ROR1 + OR2 --> CI2 + OR1 + OR2
		k14, ROR1 + ROR2 + OR3 --> CI2 + OR1 + ROR2 + OR3
		k15, ROR1 + ROR2 + ROR3 --> CI2 + OR1 + ROR2 + ROR3
		k16, ROR1 + ROR2 + COR3 --> CI2 + OR1 + ROR2 + COR3
		k17, ROR1 + COR2 --> CI2 + OR1 + COR2
		k18, ROR2 + OR1 + OR3 --> CI2 + OR2 + OR1 + OR3
		k19, ROR2 + ROR1 + OR3 --> CI2 + OR2 + ROR1 + OR3
		k20, ROR2 + OR1 + ROR3 --> CI2 + OR2 + OR1 + ROR3
		k21, ROR2 + ROR1 + ROR3 --> CI2 + OR2 + ROR1 + ROR3
		k22, ROR2 + COR1 + OR3 --> CI2 + OR2 + COR1 + OR3
		k23, ROR2 + OR1 + COR3 --> CI2 + OR2 + OR1 + COR3
		k24, ROR2 + COR1 + COR3 --> CI2 + OR2 + COR1 + COR3
		k25, ROR2 + ROR1 + COR3 --> CI2 + OR2 + ROR1 + COR3
		k26, ROR2 + COR1 + ROR3 --> CI2 + OR2 + COR1 + ROR3
		k27, ROR3 + OR2 --> CI2 + OR3 + OR2
		k28, ROR3 + ROR2 + OR1 --> CI2 + OR3 + ROR2 + OR1
		k29, ROR3 + ROR2 + ROR1 --> CI2 + OR3 + ROR2 + ROR1
		k30, ROR3 + ROR2 + COR1 --> CI2 + OR3 + ROR2 + COR1
		k31, ROR3 + COR2 --> CI2 + OR3 + COR2
		k32, COR1 + OR2 --> Cro2 + OR1 + OR2
		k33, COR1 + ROR2 --> Cro2 + OR1 + ROR2
		k34, COR1 + COR2 + OR3 --> Cro2 + OR1 + COR2 + OR3
		k35, COR1 + COR2 + ROR3 --> Cro2 + OR1 + COR2 + ROR3
		k36, COR1 + COR2 + COR3 --> Cro2 + OR1 + COR2 + COR3
		k37, COR2 + OR1 + OR3 --> Cro2 + OR2 + OR1 + OR3
		k38, COR2 + ROR1 + OR3 --> Cro2 + OR2 + ROR1 + OR3
		k39, COR2 + OR1 + ROR3 --> Cro2 + OR2 + OR1 + ROR3
		k40, COR2 + ROR1 + ROR3 --> Cro2 + OR2 + ROR1 + ROR3
		k41, COR2 + COR1 + OR3 --> Cro2 + OR2 + COR1 + OR3
		k42, COR2 + OR1 + COR3 --> Cro2 + OR2 + OR1 + COR3
		k43, COR2 + COR1 + COR3 --> Cro2 + OR2 + COR1 + COR3
		k44, COR2 + ROR1 + COR3 --> Cro2 + OR2 + ROR1 + COR3
		k45, COR2 + COR1 + ROR3 --> Cro2 + OR2 + COR1 + ROR3
		k46, COR3 + OR2 --> Cro2 + OR3 + OR2
		k47, COR3 + ROR2 --> Cro2 + OR3 + ROR2
		k48, COR3 + COR2 + OR1 --> Cro2 + OR3 + COR2 + OR1
		k49, COR3 + COR2 + ROR1 --> Cro2 + OR3 + COR2 + ROR1
		k50, COR3 + COR2 + COR1 --> Cro2 + OR3 + COR2 + COR1
	end

	return prepare_model_spec(
		rn,
		[
			:OR1 => 1,
			:OR2 => 1,
			:OR3 => 1,
			:COR1 => 0,
			:COR2 => 0,
			:COR3 => 0,
			:ROR1 => 0,
			:ROR2 => 0,
			:ROR3 => 0,
			:CI2 => 0,
			:Cro2 => 0,
		],
		[
			:k1 => 0.0069,
			:k2 => 0.0069,
			:k3 => 0.069,
			:k4 => 0.0929,
			:k5 => 0.0026,
			:k6 => 0.0025,
			:k7 => 0.021,
			:k8 => 0.021,
			:k9 => 0.021,
			:k10 => 0.021,
			:k11 => 0.021,
			:k12 => 0.021,
			:k13 => 0.00898,
			:k14 => 0.00011,
			:k15 => 0.01242,
			:k16 => 0.00011,
			:k17 => 0.00898,
			:k18 => 0.2297,
			:k19 => 0.0029,
			:k20 => 0.0021,
			:k21 => 0.0029,
			:k22 => 0.2297,
			:k23 => 0.2297,
			:k24 => 0.2297,
			:k25 => 0.0029,
			:k26 => 0.0021,
			:k27 => 1.13,
			:k28 => 0.0106,
			:k29 => 0.0106,
			:k30 => 0.0106,
			:k31 => 1.13,
			:k32 => 0.0202,
			:k33 => 0.0202,
			:k34 => 0.0040,
			:k35 => 0.0040,
			:k36 => 0.0040,
			:k37 => 0.1413,
			:k38 => 0.1413,
			:k39 => 0.1413,
			:k40 => 0.1413,
			:k41 => 0.0279,
			:k42 => 0.053,
			:k43 => 0.0328,
			:k44 => 0.053,
			:k45 => 0.0279,
			:k46 => 0.0022,
			:k47 => 0.0022,
			:k48 => 0.0008,
			:k49 => 0.0008,
			:k50 => 0.003,
		],
		(0.0, 30.0),
		[:OR1, :OR2, :OR3, :COR1, :COR2, :COR3, :ROR1, :ROR2, :ROR3, :CI2, :Cro2],
		save_points=11,
	)
end

const MODEL_SPECS = Dict{Symbol, ModelSpec}(
	:birth_death => build_birth_death_model(),
	:goutsias => Goutsias100,
	:goutsias100 => Goutsias100,
	:goutsias300 => Goutsias300,
	:goutsias1000 => Goutsias1000,
	:michaelis_menten => build_michaelis_menten_model(),
	:mapk => build_mapk_model(),
	:toggle => build_toggle_model(),
	:phage => build_phage_model(),
)

function get_model_key(model_name)
	key = normalize_model_name(model_name)
	haskey(MODEL_ALIASES, key) || error("Unknown model '$(model_name)'. Available models: birth_death, goutsias, goutsias100, goutsias300, goutsias1000, michaelis_menten, mapk, toggle, phage.")
	return MODEL_ALIASES[key]
end

function get_model_spec(model_name)
	return MODEL_SPECS[get_model_key(model_name)]
end

default_ensemble_alg() = EnsembleSerial()

function solve_ssa_ensemble(spec::ModelSpec; ntraj=100, seed=1, ensemble_alg=default_ensemble_alg())
	Random.seed!(seed)
	eprob = EnsembleProblem(spec.jump_prob)
	return solve(
		eprob,
		SSAStepper(),
		ensemble_alg;
		trajectories=ntraj,
		saveat=spec.saveat,
		save_start=true,
	)
end

function solve_ssa_ensemble(rn, u0, p, tspan; ntraj=100, seed=1, ensemble_alg=default_ensemble_alg(), save_points=11)
	spec = prepare_model_spec(rn, u0, p, tspan, [first(kv) for kv in u0]; save_points=save_points)
	return solve_ssa_ensemble(spec; ntraj=ntraj, seed=seed, ensemble_alg=ensemble_alg)
end

function final_time_mean(ensemble_sol, species_names, species_index)
	means = Dict{Symbol, Float64}()
	for sp in species_names
		means[sp] = mean(sol.u[end][species_index[sp]] for sol in ensemble_sol.u)
	end
	return means
end

function benchmark_ssa(rn, u0, p, tspan; ntraj=1000, seed=1, save_points=11)
	elapsed = @elapsed begin
		spec = prepare_model_spec(rn, u0, p, tspan, [first(kv) for kv in u0]; save_points=save_points)
		ensemble_sol = solve_ssa_ensemble(spec; ntraj=ntraj, seed=seed)
		println("Final-time means: ", final_time_mean(ensemble_sol, spec.output_species, spec.species_index))
	end
	println("Elapsed time for ", ntraj, " trajectories: ", elapsed, " seconds")
end

function run_model(model_name; ntraj=100, seed=1)
	spec = get_model_spec(model_name)
	return solve_ssa_ensemble(spec; ntraj=ntraj, seed=seed)
end

function benchmark_ssa(model_name::Union{AbstractString, Symbol}; ntraj=1000, seed=1)
	spec = get_model_spec(model_name)
	elapsed = @elapsed begin
		ensemble_sol = solve_ssa_ensemble(spec; ntraj=ntraj, seed=seed)
		println("Final-time means: ", final_time_mean(ensemble_sol, spec.output_species, spec.species_index))
	end
	println("Elapsed time for ", ntraj, " trajectories: ", elapsed, " seconds")
end

const BENCHMARK_MODEL_ORDER = [
	:birth_death,
	:goutsias100,
	:goutsias300,
	:goutsias1000,
	:michaelis_menten,
	:mapk,
	:toggle,
	:phage,
]

function benchmark_all(; ntraj=10000, seed=1, models=BENCHMARK_MODEL_ORDER)
	results = Dict{Symbol, Float64}()
	println("Running ", length(models), " benchmarks sequentially.")
	for model in models
		println("\n=== Benchmark: ", model, " ===")
		spec = get_model_spec(model)
		elapsed = @elapsed begin
			ensemble_sol = solve_ssa_ensemble(spec; ntraj=ntraj, seed=seed)
			println("Final-time means: ", final_time_mean(ensemble_sol, spec.output_species, spec.species_index))
		end
		results[model] = elapsed
		println("Elapsed time for ", ntraj, " trajectories: ", elapsed, " seconds")
	end

	total_elapsed = sum(values(results))
	println("\n=== Benchmark Summary ===")
	for model in models
		println(model, " => ", results[model], " seconds")
	end
	println("Total elapsed: ", total_elapsed, " seconds")
	return results
end

function available_models()
	return collect(keys(MODEL_SPECS))
end

# Example usage:
# julia --project=catalyst_environment -e 'include("Catalyst.jl"); benchmark_ssa(:goutsias100; ntraj=10000, seed=1)'
# julia --project=catalyst_environment -e 'include("Catalyst.jl"); benchmark_ssa(:goutsias300; ntraj=10000, seed=1)'
# julia --project=catalyst_environment -e 'include("Catalyst.jl"); benchmark_ssa(:goutsias1000; ntraj=10000, seed=1)'
# julia --project=catalyst_environment -e 'include("Catalyst.jl"); benchmark_ssa(:birth_death; ntraj=10000, seed=1)'
# julia --project=catalyst_environment -e 'include("Catalyst.jl"); benchmark_ssa(:michaelis_menten; ntraj=10000, seed=1)'
# julia --project=catalyst_environment -e 'include("Catalyst.jl"); benchmark_ssa(:mapk; ntraj=10000, seed=1)'
# julia --project=catalyst_environment -e 'include("Catalyst.jl"); benchmark_ssa(:toggle; ntraj=10000, seed=1)'
# julia --project=catalyst_environment -e 'include("Catalyst.jl"); benchmark_ssa(:phage; ntraj=10000, seed=1)'
# julia --project=catalyst_environment -e 'include("Catalyst.jl"); benchmark_all(ntraj=10000, seed=1)'