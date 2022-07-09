
using CairoMakie, ...plotting_style, Jedi

plotting_style.default_makie!()

abstract type inference_model  end


struct exponential <: inference_model
    λ_params::AbstractVector
    y0_params::AbstractVector
    σ_params::AbstractVector

end
exponential(;λ_params=[0, 0.005], y0_params=[-5, 0.5], σ_params=[-4, 2]) = exponential(λ_params, y0_params, σ_params)

function evaluate(
        x::T, 
        y::T,
        model::inference_model;
        chains=4,
        procs=1
    ) where {T <: AbstractVector} 
    #=
    if procs > 1
        # Get up to 4 workers, comment out if only one worker is wanted
        addprocs(procs - nprocs())
    end

    if typeof(model) == exponential
        @everywhere begin
            dir = @__DIR__
            exponential() = exponential()
            include("/$dir/inference/exponential_model.jl")
        end
    end
    =#
    # Import model
    dir = @__DIR__
    include("/$dir/inference/exponential_model.jl")
    
    # Run model
    mod_func = inference_model(x, y, model)
    chain = sample(mod_func, NUTS(0.65), MCMCThreads(), 1000, procs)
    chains_params = Turing.MCMCChains.get_sections(chain, :parameters)
    gen = generated_quantities(mod_func, chains_params)

    return chain, gen
end


"""
    function prior_check(
        model::inference_model
        n_samples::Int=1000
    )

Perform prior predictive check on a model.
"""
function prior_check(
        model::inference_model,
        x::AbstractVector=[],
        n_samples::Int=1000
    )
    dir = @__DIR__
    if typeof(model) == exponential
        include("/$dir/inference/exponential_model.jl")
    else

    end

    # Run Prior Predictive Check
    if length(x) == 0
        mod_func = inference_model([0], [0], model)
    else
        mod_func = inference_model(x, zeros(length(x)), model)
    end


    chain = sample(mod_func, Prior(), n_samples)

    chains_params = Turing.MCMCChains.get_sections(chain, :parameters)
    gen = generated_quantities(mod_func, chains_params) |> vec
    
    # Plot Results
    parameters = fieldnames(typeof(model))
    num_params = length(parameters)

    # Check for prior predictive checks
    num_ppc = 0
    ppc_inds = []
    if length(x) > 0
        ppc_inds = filter(x -> occursin("_ppc", x), keys(gen[1]))
        num_ppc = length(ppc_inds)
    end 

    fig = Figure(resolution=(300, 300 * (num_params + num_ppc)))
    
    for (i, p) in enumerate(parameters)
        param = split(string(p), "_")[1]
        ax = Axis(fig[i, 1])
        ax.xlabel = param
        ax.ylabel = "ECDF"

        # Check if parameter is in chain or generated quantities
        if Symbol(param) in keys(chain)
            _x = chain[Symbol(param)].data |> vec
        elseif param in keys(gen[1])
            _x = [g[param] for g in gen]
        else
            throw(ErrorException("Parameter $(param) could not be found."))
        end

        lines!(ax, sort(_x), 1 / length(_x) : 1 / length(_x) : 1)
    end

    for (i, y_ppc) in enumerate(ppc_inds)
        ax = Axis(fig[i + num_params, 1])
        ax.xlabel = y_ppc
        ax.ylabel = "ECDF"

        # Check if parameter is in chain or generated quantities
        y = [g[y_ppc] |> vec for g in gen]
        for _y in y[1:10:end]
            lines!(ax, x, _y, alpha=0.2, linewidth=1, color="gray")
        end
    end



    return fig, chain
end
