using Turing, jlArchetype, LinearAlgebra, Distributions

@model function inference_model(
        t, 
        y,
        model::exponential     
    )

    ## Priors
    # Doubling time
    λ ~ truncated(Normal(model.λ_params...), lower=model.λ_params[1])
    # Initial OD
    y0 ~ truncated(Normal(model.y0_params...), lower=y0_params[1])
    # Standard deviation of likelihood
    σ ~ LogNormal(model.σ_params...)




    # Normal Likelihood
    y ~ Turing.MvNormal(log(y0)  .+ λ .* t, σ^2)
    if length(t) > 1
        y_ppc = rand(Distributions.MvNormal(log(y0)  .+ λ .* t, σ^2), 1)
    else 
        y_ppc = rand(Distributions.Normal.(log(y0)  .+ λ .* t, σ^2), 1)
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
    offset=10^-10

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

    ## Posterior predictive check

    # Define variables for evaluation

    # number of data points
    N_data = length(y)
    # number of time points on ppc
    x_ppc = model.x_ppc
    N_ppc = length(x_ppc)


    # Build necessary covariance matrices
    ## 1. Build bottom left covariance matrix K₂₂
    ### 1.1 Build Kx*x*
    K_xs_xs = jlArchetype.bayes.cov_exp_quad(x_ppc, α, ρ; offset)

    ### 1.2 Initialize d1x_Kx*x* and d2x_Kx*x*
    d1x_K_xs_xs = Matrix{Float64}(undef, N_ppc, N_ppc)
    d2x_K_xs_xs = Matrix{Float64}(undef, N_ppc, N_ppc)

    ### 1.3 Inititalize dxx_Kx*x*
    d2xx_K_xs_xs = Matrix{Float64}(undef, N_ppc, N_ppc)

    ### 1.4 Compute derivatives of the matrices by multiplying by corresponding
    ### prefactors
    for i=1:N_ppc
        for j=1:N_ppc
            d1x_K_xs_xs[i, j] = -1 / ρ^2 * (x_ppc[i] - x_ppc[j]) * K_xs_xs[i, j]
            d2x_K_xs_xs[i, j] = 1 / ρ^2 * (x_ppc[i] - x_ppc[j]) * K_xs_xs[i, j]
            d2xx_K_xs_xs[i, j] = 1 / ρ^2 * 
                            (1 - (x_ppc[i] - x_ppc[j])^2 / ρ^2) * K_xs_xs[i, j]
        end # for
    end # for

    ### 1.5 Concatenate matrices
    K_22_top = hcat(K_xs_xs, d1x_K_xs_xs');
    K_22_bottom = hcat(d1x_K_xs_xs, d2xx_K_xs_xs);
    K_22 = vcat(K_22_top, K_22_bottom);


    ## 2. Compute top right and bottom left matrices K₁₂, K₂₁jj
    ### 2.1 Build Kxx*
    K_x_xs = jlArchetype.bayes.cov_exp_quad(x, x_ppc, α, ρ);
    K_xs_x = jlArchetype.bayes.cov_exp_quad(x_ppc, x, α, ρ);

    ### 2.2 Initialize d1x_Kx*x and d2x_Kxx*
    d2x_K_x_xs = Matrix{Float64}(undef, N_data, N_ppc)
    d1x_K_xs_x = Matrix{Float64}(undef, N_ppc, N_data)

    ### 2.3 Compute derivative of matrices by multiplying by corresonding
    ### prefactors
    for i = 1:N_data
        for j = 1:N_ppc
            d2x_K_x_xs[i, j] = 1 / ρ^2 * (x[i] - x_ppc[j]) * K_x_xs[i, j]
            d1x_K_xs_x[j, i] = - 1 / ρ^2 * (x_ppc[j] - x[i]) * K_xs_x[j, i]
        end # for
    end # for

    ### 2.5 Concatenate matrices
    K_12 = hcat(K_x_xs, d2x_K_x_xs);
    K_21 = vcat(K_xs_x, d1x_K_xs_x);

    ## 3. Solve equation Kxx * a = y
    ### 3.1 Generate covariance matrix for the data Kxx
    K_x_x = jlArchetype.bayes.cov_exp_quad(x, α, ρ) .+ LinearAlgebra.I(N_data) * σ^2

    ### 3.2 Perform Cholesky decomposition Kxx = Lxx * Lxx'
    L_x_x = LinearAlgebra.cholesky(K_x_x).L

    ### 3.3 Solve for b = inv(Lxx) y taking advantage that Lxx is a triangular
    ### matrix
    b = L_x_x\y

    ### 3.4 Solve a = inv(Lxx') b taking advantage that Lxx is a triangular
    ### matrix. Recall that a = inv(Kxx) y
    a = L_x_x'\b

    ## 4. Compute conditional mean ⟨[f(x*), dx*f(x*)] | f(x)⟩
    mean_conditional = K_21 * a

    ## 5. Evaluate v = inv(Lxx) * Kxx*
    v = L_x_x \ K_12

    ## 6. Evaluate v' = inv(Lxx) * Kx*x
    v_prime = (L_x_x \ K_21')'
    
    ## 7. Compute conditional covariance
    cov_conditional = K_22 - v_prime * v .+ LinearAlgebra.I(2 * N_ppc) * offset

    # Generate random samples given the conditional mean and covariance
    f_predict = rand(Distributions.MvNormal(
        mean_conditional, LinearAlgebra.Symmetric(cov_conditional)
        ))
    y_predict = rand(
        Distributions.Normal(f_predict[1:N_ppc], σ_)
    ) 
    dy_predict = rand(
        Distributions.Normal(f_predict[N_predict + 1:end], 1E-10)
    )
    return Dict("y_predict" => y_predict, "dy_predict" => dy_predict)
end # @model function