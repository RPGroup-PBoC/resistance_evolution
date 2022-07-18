using Turing, LinearAlgebra, Distributions, jlArchetype

@model function inference_model(
        t::Vector{<:Real}, 
        y::Vector{<:Real},
        model::exponential     
    )

    ## Priors
    # Growth Rate
    λ ~ truncated(Normal(model.λ_params...), lower=model.λ_params[1])
    # Initial OD
    y0 ~ truncated(Normal(model.y0_params...), lower=model.y0_params[1])
    # Standard deviation of likelihood
    σ ~ truncated(Normal(model.σ_params...), lower=model.σ_params[1])




    # Normal Likelihood
    y ~ Turing.MvNormal(log(y0)  .+ λ .* t, σ)
    if length(t) > 1
        y_ppc = rand(Distributions.MvNormal(log(y0)  .+ λ .* t, σ), 1) |> collect
    else 
        y_ppc = rand(Distributions.Normal(log(y0)  + λ * t[1], σ), 1)
    end
    # Return transformed parameters
    return Dict("y_ppc" => exp.(y_ppc))
end


# This function was heavily influenced by Manuel Razo-Mejia
@model function inference_model(
        x::Vector{<:Real}, 
        y::Vector{<:Real},
        model::gaussian_process
    )

    # Priors
    α ~ Turing.LogNormal(model.α_params...)
    ρ ~ Turing.LogNormal(model.ρ_params...)
    σ ~ Turing.LogNormal(model.σ_params...)

    # Compute Exponentiated quadratic covariance function
    cov_exp = jlArchetype.bayes.cov_exp_quad(x, α, ρ; offset=1E-10) + 
              (LinearAlgebra.I(length(x)) * σ^2)

    # Check if covariance matrix is positive definite. If not, set probability
    # to -∞. NOTE: I tried everything to make sure the matrix was positive 
    # definite, but I couldn't get it to work. This is a way to get around the
    # problem
    if !LinearAlgebra.isposdef(cov_exp)
        Turing.@addlogprob! -Inf
        return 
    end

    # Likelihood
    y ~ Turing.MvNormal(zeros(length(x)), cov_exp)
    return Dict()
end