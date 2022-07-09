using CairoMakie, CSV, DataFrames, Jedi, StanSample, Statistics, Printf, ColorSchemes, ...plotting_style, ...inference

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
    function run_exponential_model(file; lower_bound=exp(-4), upper_bound=exp(-2))

Fit an exponential model to a dataset from the plate reader. Can give lower and upper bound.

# Parameters
------------
- `file`: Filename of data. Has to be string.
- `lower_bound`: lower bound used for inference of exponential growth model, default exp(-4).
- `upper_bound`: lower bound used for inference of exponential growth model, default exp(-4).
"""
function run_exponential_model(file; lower_bound=exp(-4), upper_bound=exp(-2))
    
    dir = @__DIR__
    home_dir = joinpath(split(dir, '/')[1:end-2]...)

    # Find date
    date = split(file, "_")[1]
    run = split(file, "_")[2]

    # Create file path
    file_path ="/$home_dir/processing/plate_reader/$file/growth_plate.csv"

    exp_stan = open("/$home_dir/stan_code/exponential_model.stan") do file
        read(file, String)
    end

    # Compile stan code
    println("Compiling Stan files...")
    if ~isdir("/$home_dir/processing/plate_reader/$file/stan")
        mkdir("/$home_dir/processing/plate_reader/$file/stan")
    end
    sm = SampleModel("exponential_growth", exp_stan, "/$home_dir/processing/plate_reader/$file/stan")

    # Read output
    data_df = CSV.read(file_path, DataFrame)
    # Get wells
    wells = data_df.well |> unique
    wells = wells

    # Define plotting canvas
    fig_exp = Figure(resolution=(300, 350*length(wells)))

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
        _df = DataFrames.DataFrame(
            strain=sub_df.strain |> unique, 
            pos_selection=sub_df.pos_selection |> unique, 
            well=well,
        )

        # Run Exponential model
        println("  Running Exponential Growth Model...")

        ind1 = findfirst(t -> t > lower_bound, y)
        ind2 = findfirst(t -> t > upper_bound, y)
        x_exp = x[ind1:ind2]
        y_exp = y[ind1:ind2]

        chn = inference.exponential_model.evaluate(x_exp, y_exp)


        df_exp = vcat(df_exp_list...)

        if summary_exp[summary_exp.parameters .== :divergent__, "mean"][1] != 0
            println("There were divergences! $(summary_exp[summary_exp.parameters .== :divergent__, "mean"][1])")
        end

        insertcols!(
            _df, 
            :exp_growth_rate=>mean(df_exp[!, "lambda"]),
        )
        append!(return_sum_df, _df)
        println(mean(_df.exp_growth_rate))

        fig_exp = plot_samples_exp(df_exp, [x_exp, log.(y_exp)], fig_exp, i, "Exponential Growth: $well\n $(sub_df.strain[1])\n $(sub_df.pos_selection[1])")
    end

    insertcols!(return_sum_df, 1, :run=>run)
    CSV.write("/$home_dir/processing/plate_reader/$file/gp_analysis_summary.csv", return_sum_df)

    save("/$home_dir/processing/plate_reader/$file/exp_model_analysis.pdf", fig_exp) 
    
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

