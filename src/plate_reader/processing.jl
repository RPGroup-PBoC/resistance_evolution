using XLSX, CairoMakie, DataFrames, ColorSchemes, Dates, CSV, Statistics, ...plotting_style

plotting_style.default_makie!()

sequential_cmaps = [
    "Blues",
    "Reds",
    "Greens",
    "Oranges",
    "Purples"
]


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
        cmap = rand(sequential_cmaps)
        
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
        cmap = colorschemes[Symbol("$(sequential_cmaps[i])_$(max(3, length(unique(data.feature_code))))")]
        # Generate colormap
        if (data.feature_code |> unique |> length) > 1
            heatmap!(ax, data.col_code, data.row_code, data.feature_code, strokewidth = 1, strokecolor="black", colormap=cmap)
        else
            ax.backgroundcolor = "white"
        end

        text!(ax, string.(data[!, feature]), position = Point.(data.col_code, data.row_code), align = (:center, :center),
        offset = (0, 0), color = :black, textsize=3)
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
        types=Dict("Column1" => Dates.Time)
        )[1:end-2,:]

    time_points = measurements[!, "Column1"]
    time_points = Dates.second.(time_points) ./60 .+ Dates.minute.(time_points) .+ 60 .* Dates.hour.(time_points)

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
        end
        append!(df, DataFrames.dropmissing(_df))
    end

    # Blank normalization
    df_blank = df[df.strain .== "blank", :]
    gdf = DataFrames.groupby(df_blank, [:pos_selection])
    mean_blank_df = combine(gdf, :OD600 => (x -> OD_mean=mean(x)))
    df = df[df.strain .!= "blank", :]
    insertcols!(df, 4, :OD600_norm => df.OD600)

    for pos in mean_blank_df.pos_selection |> unique
        df[df.pos_selection .== pos, "OD600_norm"] = df[df.pos_selection .== pos, "OD600"] .- mean_blank_df[mean_blank_df.pos_selection .== pos, "OD600_function"][1]
    end

    ##

    # Check if output directory exists
    if ~ispath("/$home_dir/processing/plate_reader/$file")  # Check if directory exists
        mkdir("/$home_dir/processing/plate_reader/$file")  # Generate directory if required
    end
    insertcols!(df, 1, :run=>run)
    # Save output
    CSV.write("/$home_dir/processing/plate_reader/$file/growth_plate.csv", df)

    # Make plot

    l = df.strain |> unique |> length
    fig = Figure(resolution=(800, 400 * l))
    df_means = DataFrames.DataFrame()

    for (i, (strain, color)) in enumerate(zip(df.strain |> unique, ColorSchemes.seaborn_colorblind))
        ax1 = Axis(fig[i, 1])
        ax1.xlabel = "time [min]"
        ax1.ylabel = "normalized OD600"

        ax2 = Axis(fig[i, 2])
        ax2.xlabel = "time [min]"
        ax2.ylabel = "log(normalized OD600)"
        sub_df = df[df.strain .== strain, :]
        for pos_selection in sub_df.pos_selection |> unique
            sub_sub_df = sub_df[sub_df.pos_selection .== pos_selection, :]
            _gdf = DataFrames.groupby(sub_sub_df, :time_min)
            mean_df = combine(_gdf, :OD600_norm => (x -> (OD_mean=mean(x), OD_std=std(x))) => AsTable)
            insertcols!(mean_df, "strain"=>strain, "pos_selection"=>pos_selection)
            append!(df_means, mean_df)
            lines!(ax1, mean_df.time_min, mean_df.OD_mean, label="$strain, $pos_selection")
            scatter!(ax1, mean_df.time_min, mean_df.OD_mean, markersize=5)
            errorbars!(ax1, mean_df.time_min, mean_df.OD_mean, mean_df.OD_std)
            
            lines!(ax2, mean_df.time_min, log.(mean_df.OD_mean), label="$strain, $pos_selection")
            scatter!(ax2, mean_df.time_min, log.(mean_df.OD_mean), markersize=5)
            #errorbars!(ax2, mean_df.time_min, log.(mean_df.OD_mean), log.(mean_df.OD_std))
        end
        axislegend(ax1, position=:lt)#, merge = merge, unique = unique)
    end

    CSV.write("/$home_dir/processing/plate_reader/$file/growth_mean.csv", df_means)

    CairoMakie.save("/$home_dir/processing/plate_reader/$file/all_curves.pdf", fig)
    return fig
end