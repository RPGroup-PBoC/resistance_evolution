using CairoMakie

my_color_dict = Dict(
    "orange1" => "#f47c20",
    "orange2" => "#fecc96",
    "orange3" => "#ffe4c6",
    "yellow1" => "#fce317",
    "yellow2" => "#fff182",
    "yellow3" => "#fff8c1",
    "green1" => "#a8cf38",
    "green2" => "#d1e39b",
    "green3" => "#e6f0cb",
    "blue1" => "#324fa2",
    "blue2" => "#8d92c8",
    "blue3" => "#dbddef",
    "purple1" => "#9f2260",
    "purple2" => "#cca6b6",
    "purple3" => "#e9d1da",
    "red1" => "#D14241",
    "red2" => "#E59C8C",
    "red3" => "#F0CABF"
    )


function default_makie!()
    if ~isfile(assetpath("fonts", "Lucida-sans-Unicode-Regular.ttf"))
        #@warn "Lucida sans Unicode Regular font not added to Makie Fonts. Add to `~/.julia/packages/Makie/gQOQF/assets/fonts/`. Defaulting to NotoSans."
        Font = assetpath("fonts", "Lato-Regular.tff")
    else
        Font = assetpath("fonts", "Lucida-Sans-Unicode-Regular.ttf")
    end
    Font = assetpath("fonts", "Lato-Regular.tff")
    # Seaborn colorblind
    colors = ["#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc", "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9"]
    colors_new = ["#607794", "#946091", "#947d60", "#609463",
    "#A7B5C9", "#C9A7C7", "#C9B9A7", "#A7C9A9"]
    


    theme = Theme(
        Axis = (
            backgroundcolor = "#E3E7E9",
 
            # Font sizes
            titlesize=13,
            xlabelsize=13,
            ylabelsize=13,
            xticklabelsize=10,
            yticklabelsize=10,

            # Font styles
            titlefont=Font,
            xticklabelfont=Font,
            yticklabelfont=Font,
            xlabelfont=Font,
            ylabelfont=Font,

            # Grid
            xgridwidth=1.25,
            ygridwidth=1.25,
            xgridcolor="white",
            ygridcolor="white",
            xminorgridcolor="white",
            yminorgridcolor="white",
            xminorgridvisible=true,
            xminorgridwidth=1.,
            yminorgridvisible=true,
            yminorgridwidth=1,

            # Box
            rightspinevisible=false,
            topspinevisible=false,

            # Colorscheme
            palette = (color = colors_new,)

        ),
        Legend = (
            titlesize=10,
            labelsize=10,
            bgcolor="#E3E7E9",
            rowgap=-5,
            labelfont=Font

        ),
        backgroundcolor="white",
        linewidth=1.25,

    )
    set_theme!(theme)
end


sns_colorblind = ["#607794", "#946091", "#947d60", "#609463",
    "#A7B5C9", "#C9A7C7", "#C9B9A7", "#A7C9A9"]