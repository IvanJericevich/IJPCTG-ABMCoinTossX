#----- Dependencies -----#
using Printf
import StatsBase: var
import LinearAlgebra: rmul!
import Random: rand
import NLSolversBase: value, value!, value!!, NonDifferentiable
#---------------------------------------------------------------------------------------------------

#----- Structures -----#
struct OptimizationState{Tf <: Real}
    iteration::Int64
    value::Tf
    g_norm::Tf
    metadata::Dict
end
const OptimizationTrace{Tf} = Vector{OptimizationState{Tf}}
mutable struct OptimizationResults{Tx <: AbstractArray, Tf <: Real}
    initial_x::Tx
    minimizer::Tx
    minimum::Tf
    iterations::Int64
    iteration_converged::Bool
    f_abstol::Float64
    g_converged::Bool
    g_abstol::Float64
    f_increased::Bool
    trace::OptimizationTrace{Tf}
    time_limit::Float64
    time_run::Float64
    stopped_by_time_limit::Bool
end
struct Options{T <: Real}
    f_abstol::T
    g_abstol::T
    iterations::Int64
    store_trace::Bool
    trace_simplex::Bool
    show_trace::Bool
    extended_trace::Bool
    show_every::Int64
    time_limit::Float64
    ξ::Float64
    ta_rounds::Int64
end
function Options(; f_abstol::Real = 0.0, ta_rounds::Int64 = 0, g_abstol::Real = 1e-8, iterations::Int64 = 1_000, store_trace::Bool = false, trace_simplex::Bool = false, show_trace::Bool = false, extended_trace::Bool = false, show_every::Int64 = 1, time_limit = NaN, ξ = 0.0)
    Options(f_abstol, g_abstol, iterations, store_trace, trace_simplex, show_trace, extended_trace, show_every, Float64(time_limit), ξ, ta_rounds)
end
abstract type Simplexer end
struct AffineSimplexer <: Simplexer
    a::Float64
    b::Float64
    AffineSimplexer(; a = 0.025, b = 0.5) = new(a, b)
end
abstract type NMParameters end
struct AdaptiveParameters <: NMParameters
    α::Float64
    β::Float64
    γ::Float64
    δ::Float64
    AdaptiveParameters(; α = 1.0, β = 1.0, γ = 0.75 , δ = 1.0) = new(α, β, γ, δ)
end
struct FixedParameters <: NMParameters
    α::Float64
    β::Float64
    γ::Float64
    δ::Float64
    FixedParameters(; α = 1.0, β = 2.0, γ = 0.5, δ = 0.5) = new(α, β, γ, δ)
end
parameters(P::AdaptiveParameters, n::Integer) = (P.α, P.β + 2/n, P.γ - 1/2n, P.δ - 1/n)
parameters(P::FixedParameters, n::Integer) = (P.α, P.β, P.γ, P.δ)
struct NelderMead{Ts <: Simplexer, Tp <: NMParameters}
    initial_simplex::Ts
    parameters::Tp
end
function NelderMead(; kwargs...)
    KW = Dict(kwargs)
    if haskey(KW, :initial_simplex) || haskey(KW, :parameters)
        initial_simplex, parameters = AffineSimplexer(), AdaptiveParameters()
        haskey(KW, :initial_simplex) && (initial_simplex = KW[:initial_simplex])
        haskey(KW, :parameters) && (parameters = KW[:parameters])
        return NelderMead(initial_simplex, parameters)
    else
        return NelderMead(AffineSimplexer(), AdaptiveParameters())
    end
end
mutable struct NelderMeadState{Tx <: AbstractArray, T <: Real, Tfs <: AbstractArray}
    x::Tx
    m::Int
    simplex::Vector{Tx}
    x_centroid::Tx
    x_lowest::Tx
    x_second_highest::Tx
    x_highest::Tx
    x_reflect::Tx
    x_cache::Tx
    f_simplex::Tfs
    nm_x::T
    f_lowest::T
    i_order::Vector{Int}
    α::T
    β::T
    γ::T
    δ::T
    step_type::String
end
#---------------------------------------------------------------------------------------------------

#----- Nelder-Mead -----#
NelderMeadObjective(y::Vector, m::Integer, n::Integer) = sqrt(var(y) * (m / n))
function simplexer(S::AffineSimplexer, initial_x::Tx) where Tx <: AbstractArray
    n = length(initial_x)
    initial_simplex = Tx[copy(initial_x) for i = 1:(n + 1)]
    for j ∈ eachindex(initial_x)
        initial_simplex[j + 1][j] = (1 + S.b) * initial_simplex[j + 1][j] + S.a
    end
    initial_simplex
end
function centroid!(c::AbstractArray{T}, simplex, h::Integer = 0) where T # centroid except h-th vertex
    n = length(c)
    fill!(c, zero(T))
    for i in eachindex(simplex)
        if i != h
            xi = simplex[i]
            c .+= xi
        end
    end
    rmul!(c, T(1)/n)
end
centroid(simplex, h::Integer) = centroid!(similar(simplex[1]), simplex, h)
function InitialState(nelder_mead::NelderMead, f::NonDifferentiable, initial_x::Tx) where Tx <: AbstractArray
    T = eltype(initial_x)
    n = length(initial_x)
    m = n + 1
    simplex = simplexer(nelder_mead.initial_simplex, initial_x)
    f_simplex = zeros(T, m)
    value!!(f, first(simplex))
    f_simplex[1] = value(f)
    for i in 2:length(simplex)
        f_simplex[i] = value(f, simplex[i])
    end
    i_order = sortperm(f_simplex) # Get the indices that correspond to the ordering of the f values at the vertices. i_order[1] is the index in the simplex of the vertex with the lowest function value, and i_order[end] is the index in the simplex of the vertex with the highest function value
    α, β, γ, δ = parameters(nelder_mead.parameters, n)
    NelderMeadState(copy(initial_x), # Variable to hold final minimizer value for MultivariateOptimizationResults
        m, # Number of vertices in the simplex
        simplex, # Maintain simplex in state.simplex
        centroid(simplex,  i_order[m]), # Maintain centroid in state.centroid
        copy(initial_x), # Store cache in state.x_lowest
        copy(initial_x), # Store cache in state.x_second_highest
        copy(initial_x), # Store cache in state.x_highest
        copy(initial_x), # Store cache in state.x_reflect
        copy(initial_x), # Store cache in state.x_cache
        f_simplex, # Store objective values at the vertices in state.f_simplex
        T(NelderMeadObjective(f_simplex, n, m)), # Store NelderMeadObjective in state.nm_x
        f_simplex[i_order[1]], # Store lowest f in state.f_lowest
        i_order, # Store a vector of rankings of objective values
        T(α), T(β), T(γ), T(δ), "initial")
end
function ThresholdAccepting!(f::NonDifferentiable, state::NelderMeadState, thresholds)
    println("Hello")
    n, m = length(state.x), state.m
    for τ in thresholds # Should the solution not be generated from the best
        for param_index in 1:n # Perturb each parameter
            perturbation = (rand() - 0.5) * (sum(getindex.(state.simplex, param_index)) / n) # (sum(state.simplex[i][param_index] for i in 1:n) / n)
            for i in 2:m # Never perturb the best simplex?????????? Iterate over each solution
                @. state.x_cache = state.simplex[state.i_order[i]] + perturbation
                f_perturb = value(f, state.x_cache)
                if f_perturb < state.f_simplex[state.i_order[1]] + τ # Only accept if perturbation is now worse than the best + some threshold
                    # Update state
                    copyto!(state.simplex[state.i_order[i]], state.x_cache) # Replace solution with new solution
                    @inbounds state.f_simplex[state.i_order[i]] = f_perturb # Replace objective value with new value
                end
            end
        end
    end
    state.step_type = "thresholding"
    sortperm!(state.i_order, state.f_simplex) # Sort indeces of simplexes in ascending order of objective value
end
function UpdateState!(f::NonDifferentiable, state::NelderMeadState)
    shrink = false; n, m = length(state.x), state.m # Augment the iteration counter
    centroid!(state.x_centroid, state.simplex, state.i_order[m])
    copyto!(state.x_lowest, state.simplex[state.i_order[1]])
    copyto!(state.x_second_highest, state.simplex[state.i_order[n]])
    copyto!(state.x_highest, state.simplex[state.i_order[m]])
    state.f_lowest = state.f_simplex[state.i_order[1]]
    f_second_highest = state.f_simplex[state.i_order[n]]
    f_highest = state.f_simplex[state.i_order[m]]
    # Compute a reflection
    @. state.x_reflect = state.x_centroid + state.α * (state.x_centroid - state.x_highest)
    f_reflect = value(f, state.x_reflect)
    if f_reflect < state.f_lowest # Reflection has improved the objective
        # Compute an expansion
        @. state.x_cache = state.x_centroid + state.β * (state.x_reflect - state.x_centroid)
        f_expand = value(f, state.x_cache)
        if f_expand < f_reflect # Expansion has improved the objective
            # Update state
            copyto!(state.simplex[state.i_order[m]], state.x_cache)
            @inbounds state.f_simplex[state.i_order[m]] = f_expand
            state.step_type = "expansion"
        else # Expansion did not improve the objective
            # Update state
            copyto!(state.simplex[state.i_order[m]], state.x_reflect)
            @inbounds state.f_simplex[state.i_order[m]] = f_reflect
            state.step_type = "reflection"
        end
        # Shift all order indeces, and wrap the last one around to the first (update best objective)
        i_highest = state.i_order[m]
        @inbounds for i = m:-1:2
            state.i_order[i] = state.i_order[i-1]
        end
        state.i_order[1] = i_highest
    elseif f_reflect < f_second_highest # Reflection is better than the second worst
        # Update state
        copyto!(state.simplex[state.i_order[m]], state.x_reflect)
        @inbounds state.f_simplex[state.i_order[m]] = f_reflect
        state.step_type = "reflection"
        sortperm!(state.i_order, state.f_simplex)
    else
        if f_reflect < f_highest # Reflection is better than the worst but mot better than the second worst
            # Outside contraction
            @. state.x_cache = state.x_centroid + state.γ * (state.x_reflect - state.x_centroid)
            f_outside_contraction = value(f, state.x_cache)
            if f_outside_contraction < f_reflect
                # Update state
                copyto!(state.simplex[state.i_order[m]], state.x_cache)
                @inbounds state.f_simplex[state.i_order[m]] = f_outside_contraction
                state.step_type = "outside contraction"
                sortperm!(state.i_order, state.f_simplex)
            else
                shrink = true
            end
        else # f_reflect > f_highest - new worst
            # Inside constraction
            @. state.x_cache = state.x_centroid - state.γ *(state.x_reflect - state.x_centroid)
            f_inside_contraction = value(f, state.x_cache)
            if f_inside_contraction < f_highest
                # Update state
                copyto!(state.simplex[state.i_order[m]], state.x_cache)
                @inbounds state.f_simplex[state.i_order[m]] = f_inside_contraction
                state.step_type = "inside contraction"
                sortperm!(state.i_order, state.f_simplex)
            else
                shrink = true
            end
        end
    end
    # Apply shrinkage if the worst could not be improved
    if shrink
        for i = 2:m
            ord = state.i_order[i]
            # Update state
            copyto!(state.simplex[ord], state.x_lowest + state.δ*(state.simplex[ord]-state.x_lowest))
            state.f_simplex[ord] = value(f, state.simplex[ord])
        end
        step_type = "shrink"
        sortperm!(state.i_order, state.f_simplex)
    end
    state.nm_x = NelderMeadObjective(state.f_simplex, n, m)
end
function PostProcess!(f::NonDifferentiable, state::NelderMeadState)
    sortperm!(state.i_order, state.f_simplex)
    x_centroid_min = centroid(state.simplex, state.i_order[state.m])
    f_centroid_min = value(f, x_centroid_min)
    f_min, i_f_min = findmin(state.f_simplex)
    x_min = state.simplex[i_f_min]
    if f_centroid_min < f_min
        x_min = x_centroid_min
        f_min = f_centroid_min
    end
    f.F = f_min
    state.x .= x_min
end
# Convergence
function AssessConvergence(state::NelderMeadState, options::Options)
    g_converged = state.nm_x <= options.g_abstol # Stopping criterior
    return g_converged, false
end
function InitialConvergence(f::NonDifferentiable, state::NelderMeadState, initial_x::Tx, options::Options) where Tx <: AbstractArray
    nmo = NelderMeadObjective(state.f_simplex, state.m, length(initial_x))
    !isfinite(value(f)) || nmo <= options.g_abstol#, !isfinite(nmo)
end
#---------------------------------------------------------------------------------------------------

#----- Trace -----#
function Trace!(tr::OptimizationTrace, state::NelderMeadState, iteration::Int64, options::Options, curr_time = time())
    dt = Dict()
    dt["time"] = curr_time
    if options.extended_trace
        dt["centroid"] = copy(state.x_centroid)
        dt["step_type"] = state.step_type
    end
    if options.trace_simplex
        dt["simplex"] = state.simplex
        dt["simplex_values"] = state.f_simplex
    end
    os = OptimizationState(iteration, state.f_lowest, state.nm_x, dt)
    if options.store_trace
        push!(tr, os)
    end
    if options.show_trace
        if iteration % options.show_every == 0
            show(os)
            flush(stdout)
        end
    end
end
#---------------------------------------------------------------------------------------------------

#----- Optimization -----#
function Optimize(f::NonDifferentiable{Tf, Tx}, initial_x::Tx, options::Options{T} = Options(), state::NelderMeadState = InitialState(NelderMead(), f, initial_x)) where {Tx <: AbstractArray, Tf <: Real, T <: Real}
    t₀ = time() # Initial time stamp used to control early stopping by options.time_limit
    tr = OptimizationTrace{Tf}() # Store optimization trace
    thresholds = options.ta_rounds != 0 ? reverse(range(0, stop = options.f_abstol, length = options.ta_rounds)[2:end]) : nothing
    tracing = options.store_trace || options.show_trace || options.extended_trace
    stopped, stopped_by_time_limit, f_increased = false, false, false
    g_converged = InitialConvergence(f, state, initial_x, options) # Converged if criterion is met or f is infinite
    iteration = 0 # Counter
    if options.show_trace # Print header
        @printf "Iter     Function value    √(Σ(yᵢ-ȳ)²)/n \n"
        @printf "------   --------------    --------------\n"
    end
    t = time()
    Trace!(tr, state, iteration, options, t - t₀)
    while !g_converged && !stopped_by_time_limit && iteration < options.iterations
        iteration += 1
        if rand() < options.ξ
            ThresholdAccepting!(f, state, thresholds)
        else
            UpdateState!(f, state)
        end
        g_converged, f_increased = AssessConvergence(state, options)
        if tracing
            Trace!(tr, state, iteration, options, time() - t₀)
        end
        t = time()
        stopped_by_time_limit = t - t₀ > options.time_limit
    end
    PostProcess!(f, state)
    return OptimizationResults{Tx, Tf}(initial_x, state.x, value(f), iteration, iteration == options.iterations, Float64(options.f_abstol), g_converged, Float64(options.g_abstol), f_increased, tr, options.time_limit, t - t₀, stopped_by_time_limit)
end
#---------------------------------------------------------------------------------------------------

#----- API -----#
minimizer(r::OptimizationResults) = r.minimizer
minimum(r::OptimizationResults) = r.minimum
iterations(r::OptimizationResults) = r.iterations
iteration_limit_reached(r::OptimizationResults) = r.iteration_converged
trace(r::OptimizationResults) = length(r.trace) > 0 ? r.trace : error("No trace in optimization results. To get a trace, run optimize() with store_trace = true.")
converged(r::OptimizationResults) = r.g_converged
f_abstol(r::OptimizationResults) = r.f_abstol
g_abstol(r::OptimizationResults) = r.g_abstol
initial_state(r::OptimizationResults) = r.initial_x
time_limit(r::OptimizationResults) = r.time_limit
time_run(r::OptimizationResults) = r.time_run
g_norm_trace(r::OptimizationResults) = [ state.g_norm for state in trace(r) ]
f_trace(r::OptimizationResults) = [ state.value for state in trace(r) ]
function centroid_trace(r::OptimizationResults)
    tr = trace(r)
    !haskey(tr[1].metadata, "centroid") && error("Trace does not contain centroid. To get a trace of the centroid, run optimize() with extended_trace = true")
    [ state.metadata["centroid"] for state in tr ]
end
function simplex_trace(r::OptimizationResults)
    tr = trace(r)
    !haskey(tr[1].metadata, "simplex") && error("Trace does not contain simplex. To get a trace of the simplex, run optimize() with trace_simplex = true")
    [ state.metadata["simplex"] for state in tr ]
end
function simplex_value_trace(r::OptimizationResults)
    tr = trace(r)
    !haskey(tr[1].metadata, "simplex_values") && error("Trace does not contain objective values at the simplex. To get a trace of the simplex values, run optimize() with trace_simplex = true")
    [ state.metadata["simplex_values"] for state in tr ]
end
function Base.show(io::IO, r::OptimizationResults)
    take = Iterators.take
    failure_string = "failure"
    if iteration_limit_reached(r)
        failure_string *= " (reached maximum number of iterations)"
    end
    if time_run(r) > time_limit(r)
        failure_string *= " (exceeded time limit of $(time_limit(r)))"
    end
    @printf io " * Status: %s\n\n" converged(r) ? "success" : failure_string
    @printf io " * Final objective value:     %e\n" minimum(r)
    @printf io "\n"
    @printf io " * Convergence measures\n"
    @printf io "    √(Σ(yᵢ-ȳ)²)/n %s %.1e\n" converged(r) ? "≤" : "≰" g_abstol(r)
    @printf io "\n"
    @printf io " * Work counters\n"
    @printf io "    Seconds run:   %d  (vs limit %d)\n" time_run(r) isnan(time_limit(r)) ? Inf : time_limit(r)
    @printf io "    Iterations:    %d\n" iterations(r)
    return
end
function Base.show(io::IO, trace::OptimizationTrace{<:Real})
    @printf io "Iter     Function value    √(Σ(yᵢ-ȳ)²)/n \n"
    @printf io "------   --------------    --------------\n"
    for state in trace.states
        show(io, state)
    end
    return
end
function Base.show(io::IO, t::OptimizationState{<:Real})
    @printf io "%6d   %14e    %14e\n" t.iteration t.value t.g_norm
    if !isempty(t.metadata)
        for (key, value) in t.metadata
            @printf io " * %s: %s\n" key value
        end
    end
    return
end
#---------------------------------------------------------------------------------------------------

#=
z = NonDifferentiable(x -> (1.0 / 2.0) * (x[1]^2 + 0.9 * x[2]^2), [127.0, 921.0])
a = Optimize(z, [127.0, 921.0])
minimizer(a)
time_run(a)
trace(a)

function reset!(state::NelderMeadState, obj, x)
    state.simplex = simplexer(method.initial_simplex, x)
    value!(obj, first(state.simplex))
    state.f_simplex[1] = value(obj)
    for i in 2:length(state.simplex)
        state.f_simplex[i] = value(obj, state.simplex[i])
    end
    # Get the indices that correspond to the ordering of the f values at the vertices. i_order[1] is the index in the simplex of the vertex with the lowest function value, and i_order[end] is the index in the simplex of the vertex with the highest function value
    state.i_order = sortperm(state.f_simplex)
end

value(z, [2, 2])
a = optimize(x -> (1.0 / 2.0) * (x[1]^2 + 0.9 * x[2]^2), [127.0, 921.0], NelderMead())
=#
