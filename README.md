# resistance_evolution

The code in this project is written in Julia, which can be downloaded [here](https://julialang.org/downloads/).

## Setting up Computational Environment

### Julia

To run the code in this project, you need to activate the custom environment, which can be done by starting Julia in this project folder. This can either be done by adding the Julia path as an Environment variable, or by starting the executable Julia file and navigating into the project folder (using the shell mode by typing `;`, shell mode can be left by pressing `esc`). Instructions on how Julia can be started from the command line can be found [here](https://julialang.org/downloads/platform/). Once the correct folder is selected, the working environment can be set by first entering package mode (pressing `]`) typing in the Julia REPL

```julia
julia> activate .
```

Once the environment is activated, all necessary packages can be installed with

```julia
julia> instantiate
```

Package mode can be left again by typing pressing `esc`.

To run a script use

```julia
julia> include("path/to/script.jl")
```