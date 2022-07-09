using Turing, Distributions


@model function inference_model(
        t, 
        y,
        model::exponential     
    )

    ## Priors
    # Doubling time
    λ ~ truncated(Normal(model.λ_params...), lower=0)
    # Initial OD
    y0 ~ LogNormal(model.y0_params...)
    # Standard deviation of likelihood
    σ ~ LogNormal(model.σ_params...)




    # Normal Likelihood
    y ~ Turing.MvNormal(y0 .* exp.(λ .* t), σ^2)
    if length(t) > 1
        y_ppc = rand(Distributions.MvNormal(y0 .* exp.(λ .* t), σ^2), 1)
    else 
        y_ppc = rand(Distributions.Normal.(y0 .* exp.(λ * t), σ^2), 1)
    end
    # Return transformed parameters
    return Dict("y_ppc" => y_ppc)
end