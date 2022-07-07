using CairoMakie

fit_seq_colors = ["#607794", "#946091", "#947d60", "#609463",
    "#A7B5C9", "#C9A7C7", "#C9B9A7", "#A7C9A9"]


"""
    `default_makie!()`

Set plotting default to that used in Physical Biology of the Cell, 2nd edition
for the `makie` plotting library. This can be for either the GLMakie or the
CairoMakie backends.
"""
function default_makie!()
    if ~isfile(assetpath("fonts", "Lucida-sans-Unicode-Regular.ttf"))
        #@warn "Lucida sans Unicode Regular font not added to Makie Fonts. Add to `~/.julia/packages/Makie/gQOQF/assets/fonts/`. Defaulting to NotoSans."
        Font = assetpath("fonts", "NotoSans-Regular.tff")
    else
        Font = assetpath("fonts", "Lucida-Sans-Unicode-Regular.ttf")
    end
    
    # Seaborn colorblind
    colors = ["#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc", "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9"]
    fit_seq_colors = ["#607794", "#946091", "#947d60", "#609463",
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
            palette = (color = fit_seq_colors,)

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