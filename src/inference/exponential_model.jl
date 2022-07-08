using Turing, ...plotting_style, CairoMakie

plotting_style.default_makie!()

"""
    @model function exponential(t, y) 

"""
@model function model(
        t::Union{AbstractVector, Missing}, 
        y::Union{AbstractVector, Missing},
        τ_params,
        y_0_params,
        σ_params      
    )

    ## Priors
    # Doubling time
    log_10_τ ~ Normal(τ_params...)
    # Initial OD
    y_0 ~ LogNormal(y_0_params...)
    # Standard deviation of likelihood
    σ ~ LogNormal(σ_params...)


    # Growth rate
    λ = log(2) / 10^log_10_τ

    # Normal Likelihood
    y ~ MvNormal(y_0 .* exp.(λ .* t), σ^2)
end


"""
    function prior_check(
        τ_params=[2, 0.5],
        y_0_params=[-3, 0.5],
        σ_params=[-3, 2],
        n_samples::Int=1000
    )

Perform prior predictive check on exponential growth model.
"""
function prior_check(
        τ_params=[2, 0.5],
        y_0_params=[-5, 0.5],
        σ_params=[-3, 2],
        n_samples::Int=1000
    )

    # Run Prior Predictive Check
    mod_func = model([0], [0], τ_params, y_0_params, σ_params)
    chain = sample(mod_func, Prior(), n_samples)
    
    # Plot Results
    fig = Figure(resolution=(900, 300))
    
    ax1 = Axis(fig[1, 1])
    ax1.xlabel = "λ"
    ax1.ylabel = "ECDF"
    x1 = chain[:log_10_τ].data |> vec
    x1 = log(2) ./ 10 .^ x1

    lines!(ax1, sort(x1), 1/length(x1):1/length(x1):1)

    ax2 = Axis(fig[1, 2])
    ax2.xlabel = "y_0"
    ax2.ylabel = "ECDF"
    x2 = chain[:y_0].data |> vec
    lines!(ax2, sort(x2), 1/length(x2):1/length(x2):1)

    ax3 = Axis(fig[1, 3])
    ax3.xlabel = "log10 σ"
    ax3.ylabel = "ECDF"
    x3 = chain[:σ].data |> vec 
    lines!(ax3, log10.(sort(x3)), 1/length(x3):1/length(x3):1)

    return fig, chain
end


function evaluate(
        t::AbstractVector, 
        y::AbstractVector;
        τ_params=[1.75, 0.25],
        y_0_params=[-2, 0.5],
        σ_params=[-3, 2]
    )

    mod_func = model(t, y, τ_params, y_0_params, σ_params)
    chain = sample(mod_func, NUTS(0.65), 1000)
    return chain
end