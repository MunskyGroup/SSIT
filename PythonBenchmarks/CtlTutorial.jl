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
    "OrdinaryDiffEq",
    "Plots",
])

# Fetch required packages.
using Catalyst, OrdinaryDiffEq, Plots

# Create model.
model = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end

# Create an ODE that can be simulated.
u0 = [:S => 50.0, :E => 10.0, :SE => 0.0, :P => 0.0]
tspan = (0., 200.)
ps = [:kB => 0.01, :kD => 0.1, :kP => 0.1]
ode = ODEProblem(model, u0, tspan, ps)

# Simulate ODE and plot results.
sol = solve(ode)
plot(sol; lw = 5)

# SSA simulation.
# The initial conditions are now integers as we track exact populations for each species.
using JumpProcesses
u0_integers = [:S => 50, :E => 10, :SE => 0, :P => 0]
jprob = JumpProblem(model, u0_integers, tspan, ps)
jump_sol = solve(jprob)
plot(jump_sol; lw = 2)

# Run 10000 SSA, time the process, compute the mean
using Statistics
