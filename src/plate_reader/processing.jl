using XLSX, CairoMakie, DataFrames, ColorSchemes, CSV, Statistics, ...plotting_style, Colors



function plot_layout(file)

    # Locate script
    dir = @__DIR__
    home_dir = joinpath(split(dir, '/')[1:end-2]...)

    # Find date
    date = split(file, "_")[1]
    

    # Create file path
    file_path = "/" * home_dir * "/data/plate_reader/" * file * "/" * date * "_plate_layout.xlsx"

    # Check if file exists
    if ~ispath(file_path)
        throw(ArgumentError("Plate layout file does not exist. Check spelling: $file_path"))
    end

    # Import excel sheet
    xf = XLSX.readxlsx(file_path)

    # See all sheet names
    layout = XLSX.sheetnames(xf)

    # Read layout information
    layout_info = DataFrames.DataFrame()

    # Loop through layout info
    for table in layout
        flattened = vec(hcat(XLSX.getdata(xf[table])...))
        DataFrames.insertcols!(layout_info, table => flattened)
    end

    # Define rows and columns read on plate
    DataFrames.insertcols!(layout_info, "row" => [x[1] for x in layout_info[!, "well"]])
    DataFrames.insertcols!(layout_info, "col" => [x[2:end] for x in layout_info[!, "well"]])

    ##
    # Initialize plot

    fig = Figure(resolution=(900, 300 * ((length(layout)) ÷ 2)))

    # Loop through features
    global i=0
    for feature in layout
        
        # If well, don't add information
        if (feature == "well")
            continue
        end
        global i += 1

        ax = Axis(fig[(i-1) ÷ 2 + 1, ((i-1) % 2) + 1 ])
        # Extract strain data and pivot dataframe to obtain proper dimensions
        data = layout_info[!, ["row", "col", feature]]
        # Add code for categorical data
        cat_dict = Dict(unique(data[!, feature]) .=> collect(1:length(unique(data[!, feature]))))
        data[!, "feature_code"] = [cat_dict[x] for x in data[!, feature]]
        row_dict = Dict(unique(data[!, "row"]) .=> collect(1:length(unique(data[!, "row"]))))
        data[!, "row_code"] = [row_dict[x] for x in data[!, "row"]]
        col_dict = Dict(unique(data[!, "col"]) .=> collect(1:length(unique(data[!, "col"]))))
        data[!, "col_code"] = [col_dict[x] for x in data[!, "col"]]
        # Pivot dataframe for both code and label to generate heatmap

        # Hacky way of making nicer colormaps
        if length(unique(data.feature_code)) > 1
            cmap = ColorScheme(Colors.sequential_palette(rand(1:360), length(unique(data.feature_code)) * 2, w=0, c=1)[1:length(unique(data.feature_code))])
        else
            cmap = "white"
        end
        
        # Generate colormap
        if (data.feature_code |> unique |> length) > 1
            heatmap!(ax, data.col_code, data.row_code, data.feature_code, strokewidth = 1, strokecolor="black", colormap=cmap)
        else
            ax.backgroundcolor = "white"
        end

        text!(ax, string.(data[!, feature]), position = Point.(data.col_code, data.row_code), align = (:center, :center),
        offset = (0, 0), color = :black, fontsize=3)
        # Set plot title
        ax.title = feature
        ax.xgridcolor="white"
        ax.ygridcolor="white"
        ax.xgridwidth=5
        ax.xticks=(collect(1:12), string.(collect(1:12)))
        ax.yticks=(collect(1:8), ["A", "B", "C", "D", "E", "F", "G", "H"])
        ax.yreversed = true
    end

    # Save it to the output file
    if ~ispath("/$home_dir/processing/plate_reader/$file")  # Check if directory exists
        mkdir("/$home_dir/processing/plate_reader/$file")  # Generate directory if required
    end
    save("/$home_dir/processing/plate_reader/$file/plate_layout.pdf", fig)
    return fig
end


function pre_processing(file)
    # Locate script
    dir = @__DIR__
    home_dir = joinpath(split(dir, '/')[1:end-2]...)

    # Find date
    date = split(file, "_")[1]
    run = split(file, "_")[2]

    # Create file path
    file_path_layout = "/" * home_dir * "/data/plate_reader/" * file * "/" * date * "_plate_layout.xlsx"
    file_path_data = "/" * home_dir * "/data/plate_reader/" * file * "/" * file * ".txt"

    # Check if file exists
    if ~ispath(file_path_layout)
        throw(ArgumentError("Plate layout file does not exist. Check spelling: $file_path_layout"))
    end

    if ~ispath(file_path_data)
        throw(ArgumentError("Data file does not exist. Check spelling: $file_path_data"))
    end

    xf = XLSX.readxlsx(file_path_layout)

    meta_df = DataFrames.DataFrame()

    for table in XLSX.sheetnames(xf)
        flattened = vec(hcat(XLSX.getdata(xf[table])...))
        DataFrames.insertcols!(meta_df, table => flattened)
    end
    meta_df = sort(meta_df, :well)
    ##

    measurements = CSV.read(
        file_path_data, 
        DataFrames.DataFrame, 
        header=false,
        )[1:end-2,:]

    time_strings = measurements[!, "Column1"]
    time_points = Float64[]
    for _t in time_strings
        _t = string(_t)
        push!(time_points, parse(Float64, split(_t, ':')[1]) * 60 + parse(Float64, split(_t, ':')[2]) + parse(Float64, split(_t, ':')[3]) / 60)
    end
    

    temperatures = measurements[!, "Column2"]

    df = DataFrames.DataFrame()

    for i in 3:size(measurements)[2]
        OD = measurements[!, i]
        _df = DataFrames.DataFrame(time_min=time_points, temp=temperatures, OD600=OD)
        
        for column in names(meta_df)
            DataFrames.insertcols!(_df, column => meta_df[i-2, column])
        end
        if typeof(_df.strain) == Vector{Int64}
            _df.strain = string.(_df.strain)
        elseif typeof(_df.strain) == Vector{Float64}
            _df.strain = string.(_df.strain)
        end
        append!(df, DataFrames.dropmissing(_df))
    end

    fig = Figure(resolution=(1200, 800))
    g = fig[1, 1] = GridLayout()
    row_dict = Dict(
        'A' => 1, 
        'B' => 2,
        'C' => 3,
        'D' => 4,
        'E' => 5,
        'F' => 6,
        'G' => 7,
        'H' => 8,
    )
    ax_list = []
    for well in unique(df.well)
        x = df[df.well .== well, :time_min]
        y = df[df.well .== well, :OD600]
        row = row_dict[well[1]]
        col = parse(Int64, well[2:end])
        ax = Axis(
            g[row, col+1],
            xticklabelsize=8,
            yticklabelsize=8
        )
        if col != 1
            hideydecorations!(ax, grid=false, label=false, minorgrid=false)
        end
        if row != 8
            hidexdecorations!(ax, grid=false, label=false, minorgrid=false)
        end
        lines!(ax, x, y)
        push!(ax_list, ax)
    end

 

    # Link axes
    linkyaxes!(ax_list...)
    linkxaxes!(ax_list...)

    # Add labels for all x and y axes
    Label(g[5, 1], "OD600", valign = :center, halign=:right, rotation=pi/2, padding=(30, 30, 30, 30), tellwidth=false, tellheight=false)
    Label(g[9, 1:12, Top()], "time [min]", valign = :top, halign=:center)#, rotation=pi/2)
    for i in 1:12
        Label(g[1, i+1, Top()], "$i", valign = :top, padding=(0, 0, 5, 0))#, rotation=pi/2))
    end
    for (i, x) in enumerate(['A', 'B', 'C', 'D' , 'E', 'F', 'G', 'H'])
        Label(g[i, 14], "$x", tellheight=false, padding=(10, 0, 0, 0))
    end

    colsize!(g, 1, 50)
    rowsize!(g, 9, 20)

    # Remove any gap between plots
    colgap!(g, 0)
    rowgap!(g, 0)

    CairoMakie.save("/$home_dir/processing/plate_reader/$file/well_plate_growth_curves.pdf", fig)

    # Blank normalization
    df_blank = df[df.strain .== "blank", :]
    OD_mean = df_blank[!, :OD600] |> mean
    insertcols!(df, 4, :OD600_norm => max.(df.OD600 .- OD_mean,0))
   #df.OD600_norm = [maximum(x, 0) for x in df.OD600_norm]
    df = df[df.strain .!= "blank", :]
    df = df[df.strain .!= "n", :]


    # Check if output directory exists
    if ~ispath("/$home_dir/processing/plate_reader/$file")  # Check if directory exists
        mkdir("/$home_dir/processing/plate_reader/$file")  # Generate directory if required
    end
    insertcols!(df, 1, :run=>run)
    # Save output
    CSV.write("/$home_dir/processing/plate_reader/$file/growth_plate.csv", df)

    # Make plot
    insertcols!(df, 1, :tc => [parse(Float64, split(x, "_")[1]) for x in df.pos_selection])
    l = df.strain |> unique |> length
    fig = Figure(resolution=(800, 400 * l))
    df_means = DataFrames.DataFrame()

    for (i, (strain, color)) in enumerate(zip(df.strain |> unique, ColorSchemes.seaborn_colorblind))
        ax1 = Axis(fig[i, 1])
        ax1.xlabel = "time [min]"
        ax1.ylabel = "normalized OD600"
        ax1.title = "$strain"

        ax2 = Axis(fig[i, 2])
        ax2.xlabel = "time [min]"
        ax2.ylabel = "log(normalized OD600)"
        sub_df = df[df.strain .== strain, :]

        for tc in sub_df.tc |>  unique |> sort
            sub_sub_df = sub_df[sub_df.tc .== tc, :]
            _gdf = DataFrames.groupby(sub_sub_df, :time_min)
            mean_df = combine(_gdf, :OD600_norm => (x -> (OD_mean=mean(x), OD_std=std(x))) => AsTable)
            insertcols!(mean_df, "strain"=>strain, "pos_selection"=>tc)
            append!(df_means, mean_df)
            lines!(ax1, mean_df.time_min, mean_df.OD_mean, label="$tc")
            scatter!(ax1, mean_df.time_min, mean_df.OD_mean, markersize=5)

            # plot error bars more sparse
            errorbars!(ax1, mean_df.time_min[1:5:end], mean_df.OD_mean[1:5:end], mean_df.OD_std[1:5:end])
            
            lines!(ax2, mean_df.time_min, log.(mean_df.OD_mean), label="$tc")
            scatter!(ax2, mean_df.time_min, log.(mean_df.OD_mean), markersize=5)
            #errorbars!(ax2, mean_df.time_min, log.(mean_df.OD_mean), log.(mean_df.OD_std))
        end
        fig[i, 3] = axislegend("tc [µg/ml]")#, merge = merge, unique = unique)
    end

    CSV.write("/$home_dir/processing/plate_reader/$file/growth_mean.csv", df_means)

    CairoMakie.save("/$home_dir/processing/plate_reader/$file/all_curves_by_strains.pdf", fig)


    # Make plot by concentrations
    l = df.tc |> unique |> length
    fig = Figure(resolution=(800, 400 * l))
    df_means = DataFrames.DataFrame()

    for (i, (tc, color)) in enumerate(zip(df.tc |> unique |> sort, ColorSchemes.seaborn_colorblind))
        ax1 = Axis(fig[i, 1])
        ax1.xlabel = "time [min]"
        ax1.ylabel = "normalized OD600"
        ax1.title = "tetracycline: $tc [µg/ml]"

        ax2 = Axis(fig[i, 2])
        ax2.xlabel = "time [min]"
        ax2.ylabel = "log(normalized OD600)"
        sub_df = df[df.tc .== tc, :]

        for strain in sub_df.strain |>  unique |> sort
            sub_sub_df = sub_df[sub_df.strain .== strain, :]
            _gdf = DataFrames.groupby(sub_sub_df, :time_min)
            mean_df = combine(_gdf, :OD600_norm => (x -> (OD_mean=mean(x), OD_std=std(x))) => AsTable)
            insertcols!(mean_df, "strain"=>strain, "pos_selection"=>tc)
            append!(df_means, mean_df)
            lines!(ax1, mean_df.time_min, mean_df.OD_mean, label="$strain")
            scatter!(ax1, mean_df.time_min, mean_df.OD_mean, markersize=5)

            # plot error bars more sparse
            errorbars!(ax1, mean_df.time_min[1:5:end], mean_df.OD_mean[1:5:end], mean_df.OD_std[1:5:end])
            
            lines!(ax2, mean_df.time_min, log.(mean_df.OD_mean), label="$strain")
            scatter!(ax2, mean_df.time_min, log.(mean_df.OD_mean), markersize=5)
            #errorbars!(ax2, mean_df.time_min, log.(mean_df.OD_mean), log.(mean_df.OD_std))
        end
        fig[i, 3] = axislegend("strain")#, merge = merge, unique = unique)
    end

    CSV.write("/$home_dir/processing/plate_reader/$file/growth_mean.csv", df_means)

    CairoMakie.save("/$home_dir/processing/plate_reader/$file/all_curves_by_tc.pdf", fig)

    return fig
end