using Turing, jlArchetype, LinearAlgebra


@model function growth_rate_gp(
        y, t, model::inference_model
    )
    
    # Defien parameter types
    α = Float64[]
    ρ = Float64[]
    σ = Float64[]

    # Priors
    # α ~ Turing.truncated(Turing.Normal(α_prior...), 0, Inf)
    α ~ Turing.LogNormal(inference_model.α_params...)
    ρ ~ Turing.LogNormal(inference_model.ρ_params...)
    σ ~ Turing.LogNormal(inference_model.σ_params...)

    # Compute Exponentiated quadratic covariance function
    cov_exp = cov_exp_quad(t, α, ρ; offset=1E-10) + 
              (LinearAlgebra.I(length(t)) * σ^2)

    # Check if covariance matrix is positive definite. If not, set probability
    # to -∞. NOTE: I tried everything to make sure the matrix was positive 
    # definite, but I couldn't get it to work. This is a way to get around the
    # problem
    if !LinearAlgebra.isposdef(cov_exp)
        Turing.@addlogprob! -Inf
        return 
    end

    # Likelihood
    return y ~ Turing.MvNormal(zeros(length(t)), cov_exp)
end # @model function