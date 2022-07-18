using CairoMakie, CSV, DataFrames, Jedi, Statistics, Printf, ColorSchemes, ...plotting_style, ...inference, LinearAlgebra, Distributions

CairoMakie.activate!()
plotting_style.default_makie!()


"""
    function plot_samples_exp(df, data, fig, i, title=nothing)

Add a predictive regression plot to an existing figure.

# Parameters
------------
- `df`: Dataframe containing posterior predictive check. Has to contain a column `y_ppc`.
- `data`: Array. First dimension is x-coordinate, second dimension is y-coordinate.
- `fig`: Makie Figure
- `i`: Row to add figure in.
- `title`: Title of the panel, default no title

# Returns
---------
- `fig`: Updated figure
"""
function plot_samples_exp(df, data, fig, i, title=nothing)
    # Extract variable names
    y_ppc = filter(x -> occursin("y_ppc", String(x)), names(df))

    # Prepare axes
    ax1 = Axis(fig[i, 1])
    ax1.xlabel = "time [min]"
    ax1.ylabel = "OD 600"

    if ~isnothing(title)
        ax1.title = title
    end
    ax1 = Jedi.viz.predictive_regression(
        [df[!, x] for x in y_ppc],
        data[1],
        ax1,
        data=data,
        data_kwargs=Dict(:markersize => 6)
        )
        
    return fig
end


"""
    function run_exponential_model(
            file; 
            lower_bound=exp(-4), 
            upper_bound=exp(-2),
            λ_params::AbstractVector=[0, 0.005],
            y0_params::AbstractVector=[0, 0.001],
            σ_params::AbstractVector=[-3, 2]
        )

Fit an exponential model to a dataset from the plate reader. Can give lower and upper bound.

# Parameters
------------
- `file`: Filename of data. Has to be string.
- `lower_bound`: lower bound used for inference of exponential growth model, default exp(-4).
- `upper_bound`: lower bound used for inference of exponential growth model, default exp(-2).
- `λ_params` : parameters for HalfNormal prior for growth rate, default `μ=0`, `σ=0.05`.
- `y0_params` : parameters for HalfNormal prior for initial OD, default `μ=0`, `σ=0.001`.
- `σ_params` : parameters for LogNormal prior for likelihood standard deviation, default `μ=-3`, `σ=2`.
"""
function run_exponential_model(
        file; 
        lower_bound=exp(-4), 
        upper_bound=exp(-2),
        λ_params::AbstractVector=[0, 0.005],
        y0_params::AbstractVector=[0, 0.001],
        σ_params::AbstractVector=[0, 0.01]
    )
    
    dir = @__DIR__
    home_dir = joinpath(split(dir, '/')[1:end-2]...)

    # Find date
    date = split(file, "_")[1]
    run = split(file, "_")[2]

    # Create file path
    file_path ="/$home_dir/processing/plate_reader/$file/growth_plate.csv"

    # Read output
    data_df = CSV.read(file_path, DataFrame)
    # Get wells
    wells = data_df.well |> unique
    wells = wells

    # Define plotting canvas
    fig_exp = Figure(resolution=(600, 350*length(wells)))

    # Define DataFrames
    return_sum_df = DataFrames.DataFrame()

    println("Analyzing data...")
    for (i, well) in enumerate(wells)
        println(" Running well $well...")
        # Choosing data for well
        sub_df = data_df[data_df.well .== well, :]
        x = sub_df[!, "time_min"]
        y = sub_df[!, "OD600_norm"]
        
        # Create dataframe for maximum growth rate
        #strain, rep = split(unique(sub_df.strain, "_")[1])[1, 2]
        _df = DataFrames.DataFrame(
            strain=sub_df.strain |> unique,
            pos_selection=sub_df.pos_selection |> unique, 
            well=well,
        )

        # Run Exponential model
        println("  Running Exponential Growth Model...")

        ind1 = findfirst(t -> t > lower_bound, y)
        ind2 = findfirst(t -> t > upper_bound, y)
        if ~isnothing(ind1) && ~isnothing(ind2)
            x_exp = x[ind1:ind2]
            y_exp = y[ind1:ind2]

            model = inference.exponential(
                λ_params = λ_params, 
                y0_params = y0_params, 
                σ_params = σ_params
                )
            chn, gen = inference.evaluate(x_exp, log.(y_exp), model)

            insertcols!(
                _df, 
                :exp_growth_rate=>mean(chn[:λ]),
            )
        else
            continue
        end
        append!(return_sum_df, _df)
        println(mean(_df.exp_growth_rate))

        # Prepare axes
        ax = Axis(fig_exp[i, 1])
        ax.xlabel = "time [min]"
        ax.ylabel = "OD 600"

        ax_log = Axis(fig_exp[i, 2])
        ax_log.xlabel = "time [min]"
        ax_log.ylabel = "log OD 600"

        y_ppc = [g["y_ppc"] |> vec for g in gen]

        ax = Jedi.viz.predictive_regression(
            hcat(y_ppc...),
            x_exp,
            ax,
            data=[x_exp, y_exp],
            data_kwargs=Dict(:markersize => 6)
            )

        ax_log = Jedi.viz.predictive_regression(
            log.(hcat(y_ppc...)),
            x_exp,
            ax_log,
            data=[x_exp, log.(y_exp)],
            data_kwargs=Dict(:markersize => 6)
            )
        mean_λ = mean(chn[:λ])
        mean_y0 = mean(chn[:y0]) 
        mean_sigma = mean(chn[:σ])
        ax.title = "$well, $(unique(sub_df.strain)[1]), $(unique(sub_df.pos_selection)[1])"
        ax_log.title = @sprintf "%.5f, %.5f, %.5f" mean_λ mean_y0 mean_sigma
            
    end
    insertcols!(return_sum_df, 1, :run=>run)
    CSV.write("/$home_dir/processing/plate_reader/$file/exp_analysis_summary.csv", return_sum_df)

    save("/$home_dir/processing/plate_reader/$file/exp_model_analysis.pdf", fig_exp) 

    return fig_exp
end

    
    #=
    fig = Figure(resolution=(800, 400 * (data_df.strain |> unique |> length)))


    # iterater through strains
    for (i, strain) in enumerate(data_df.strain |> unique)
        # Get data from this strain
        sub_df = data_df[data_df.strain .== strain, :]
        # Make subplots
        ax1 = Axis(fig[i, 1])
        ax1.xlabel = "time [min]"
        ax1.ylabel = "normalized OD600"
        ax1.title = "$strain"

        ax2 = Axis(fig[i, 2])
        ax2.xlabel = "time [min]"
        ax2.ylabel = "log normalized OD600"
        # Iterate through drug concentrations
        for (pos_selection, color) in zip(sub_df.pos_selection |> unique, ColorSchemes.seaborn_colorblind)
            # Get data
            sub_sub_df = sub_df[sub_df.pos_selection .== pos_selection, :]
            # Get exponential growth rate
            λ = return_sum_df[(return_sum_df.strain .== strain) .& (return_sum_df.pos_selection .== pos_selection), "exp_growth_rate"]
            mean_λ = mean(λ)
            std_λ = std(λ)
            # Compute mean OD per time step
            gdf = groupby(sub_sub_df, :time_min)
            mean_df = combine(gdf, :OD600_norm => (x -> (OD_mean=mean(x), OD_std=std(x))) => AsTable)
            insertcols!(mean_df, "strain"=>strain)

            # Find time points where exponential growth rate was computed
            ind1 = findfirst(x->x > lower_bound, mean_df.OD_mean)
            ind2 = findfirst(x->x > upper_bound, mean_df.OD_mean)
            # Compute offset
            b = log(mean_df.OD_mean[ind1]) - mean_λ * mean_df.time_min[ind1]
            log_y_exp = mean_λ .* mean_df.time_min[ind1:ind2] .+ b

            # Plot growth curve on linear scale
            #lines!(ax1, mean_df.time_min, mean_df.OD_mean, label=pos_selection)
            
            scatter!(ax1, mean_df.time_min[ind1:ind2], mean_df.OD_mean[ind1:ind2], markersize=8, label=pos_selection)
            lines!(ax1, mean_df.time_min[ind1:ind2], exp.(log_y_exp), linestyle="--", color="orange", linewidth=4)
            # Plot growth curve on log scale with exponential growth prediction
            #lines!(ax2, mean_df.time_min, log.(mean_df.OD_mean))
            
            scatter!(ax2, mean_df.time_min[ind1:ind2], log.(mean_df.OD_mean)[ind1:ind2], markersize=8)
            lines!(ax2, mean_df.time_min[ind1:ind2], log_y_exp, linestyle="--", color="orange", linewidth=4)
            errorbars!(ax1, mean_df.time_min[ind1:ind2], mean_df.OD_mean[ind1:ind2], mean_df.OD_std[ind1:ind2])
            ax2.title = @sprintf "%.5f  ± %.5f 1/min" mean_λ std_λ

        end
        axislegend(ax1, position=:lt)#, merge = merge, unique = unique)
    end


    save("/$home_dir/processing/plate_reader/$file/all_curves_with_th.pdf", fig)
    
    println("Done!")
end


"""
    exponential_model()

Use Turing.jl to fit an exponential model to growth rate data.
"""
function exponential_model()


end

=#

#=
"""
    function run_gp_model(file;)

Fit an gaussian process model to a dataset from the plate reader. Can give lower and upper bound.

# Parameters
------------
- `file`: Filename of data. Has to be string.
"""
function run_gp_model(
        file::AbstractString;
        α_params::AbstractVector=[0., 1.], 
        ρ_params::AbstractVector=[2., 1.], 
        σ_params::AbstractVector=[0., 1.]
    )
    
    dir = @__DIR__
    home_dir = joinpath(split(dir, '/')[1:end-2]...)

    # Find date
    date = split(file, "_")[1]
    run = split(file, "_")[2]

    # Create file path
    file_path ="/$home_dir/processing/plate_reader/$file/growth_plate.csv"

    # Read output
    data_df = CSV.read(file_path, DataFrame)
    # Get wells
    wells = data_df.well |> unique
    wells = wells

    # Define plotting canvas
    fig_exp = Figure(resolution=(900, 350*length(wells)))

    # Define DataFrames
    return_sum_df = DataFrames.DataFrame()

    println("Analyzing data...")
    for (i, well) in enumerate(wells)
        println(" Running well $well...")
        # Choosing data for well
        sub_df = data_df[data_df.well .== well, :]
        _x = sub_df[!, "time_min"]
        _y = sub_df[!, "OD600_norm"]

        x = _x[1:2:end]
        y = _y[1:2:end]
    
        
        # Create dataframe for maximum growth rate
        _df = DataFrames.DataFrame(
            strain=sub_df.strain |> unique, 
            pos_selection=sub_df.pos_selection |> unique, 
            well=well,
        )

        # Run Exponential model
        println("  Running Gaussian Process Growth Model...")

        # Center data
        x_cen = (x .- mean(x)) ./ std(x)
        y_cen = (y .- mean(y)) ./ std(y)

        model = inference.gaussian_process(
            x_cen,
            α_params = α_params,
            ρ_params = ρ_params,
            σ_params = σ_params,
            )

        chn, gen = inference.evaluate(x_cen, y_cen, model)
        α_arr = chn[:α].data |> vec
        ρ_arr = chn[:ρ].data |> vec
        σ_arr = chn[:σ].data |> vec

        #=
        fig = Figure(resolution=(300, 900))
        ax1 = Axis(fig[1, 1])
        ax1.xlabel = "α"
        ax1.ylabel = "ECDF"
        lines!(ax1, sort(α_arr), 1/length(α_arr):1/length(α_arr):1) 

        ax2 = Axis(fig[2, 1])
        ax2.xlabel = "ρ"
        ax2.ylabel = "ECDF"
        lines!(ax2, sort(ρ_arr), 1/length(ρ_arr):1/length(ρ_arr):1)

        ax3 = Axis(fig[3, 1])
        ax3.xlabel = "σ"
        ax3.ylabel = "ECDF"
        lines!(ax3, sort(σ_arr), 1/length(σ_arr):1/length(σ_arr):1)
        return fig
        =#

        y_ppc = []
        dy_ppc = []
        f_predict = []

        for (α, ρ, σ) in zip(α_arr, ρ_arr, σ_arr)
            _y_ppc, _dy_ppc, f_pred = gp_posterior_predictive_check(
                x, y, model.x_ppc, α, ρ, σ
            )
            # Rescale y and dy
            _y_ppc = ( _y_ppc .+ mean(y)) .* std(y)
            _dy_ppc = std(y) ./ std(x) .* _dy_ppc
            push!(y_ppc, _y_ppc)
            push!(dy_ppc, _dy_ppc)
            push!(f_predict, f_pred[1:length(x)])
        end
        
        

        λ = [maximum(_y_ppc ./ _dy_ppc) for (_y_ppc, _dy_ppc) in zip(y_ppc, dy_ppc)]

        insertcols!(
            _df, 
            :growth_rate=>mean(λ),
        )

        append!(return_sum_df, _df)

        # Prepare axes
        ax = Axis(fig_exp[i, 1])
        ax.xlabel = "time [min]"
        ax.ylabel = "OD 600"

        ax = Jedi.viz.predictive_regression(
            hcat(y_ppc...),
            x,
            ax,
            data=[x, y],
            data_kwargs=Dict(:markersize => 6)
            )

        ax2 = Axis(fig_exp[i, 2])
        ax2.xlabel = "time [min]"
        ax2.ylabel = "dOD 600/dt"
        
        ax2 = Jedi.viz.predictive_regression(
            hcat(dy_ppc...),
            x,
            ax2,
            data_kwargs=Dict(:markersize => 6)
            )

        ax3 = Axis(fig_exp[i, 3])
        ax3.xlabel = "time [min]"
        ax3.ylabel = "f predict"
        
        ax3 = Jedi.viz.predictive_regression(
            hcat(f_predict...),
            x,
            ax3,
            data_kwargs=Dict(:markersize => 6)
            )

        mean_λ = mean(λ)
        ax.title = @sprintf "%.5f" mean_λ 
            
    end
    insertcols!(return_sum_df, 1, :run=>run)
    CSV.write("/$home_dir/processing/plate_reader/$file/gp_analysis_summary.csv", return_sum_df)

    save("/$home_dir/processing/plate_reader/$file/gp_model_analysis.pdf", fig_exp) 

    return fig_exp
end


function gp_posterior_predictive_check(
    x::Vector{<:Real},
    y::Vector{<:Real},
    x_ppc::Vector{<:Real},
    α::Real,
    ρ::Real,
    σ::Real;
    offset::Real=10^-10
)
    # Define variables for evaluation

    # number of data points
    N_data = length(y)
    # number of time points on ppc
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
        Distributions.MvNormal(f_predict[1:N_ppc], σ)
    ) 
    dy_predict = rand(
        Distributions.MvNormal(f_predict[N_ppc + 1:end], 1E-10)
    )
    return y_predict, dy_predict, f_predict
end # @model function
    =#