#!/usr/bin/env julia

# Available colors: http://juliagraphics.github.io/Luxor.jl/stable/howto/colors-styles/

colors = Dict(
    "N" => "royalblue4",       "Ru" => "steelblue",
    "O" => "red3",             "Re" => "pink",
    "H" => "gray90",
    "C" => "darkgrey",
    "Cl" => "yellow",
    "S" => "orange"
)
radii = Dict(
    "N" => 1.1,         "Ru" => 1.6, 
    "O" => 1.1,         "Re" => 1.6,
    "H" => 0.8,
    "C" => 1.0, 
    "Cl" => 1.5, 
    "S" => 1.2    
)

