#! /usr/bin/env julia 

#==============================================================================
            Load procedures:
==============================================================================#

# Load file readers:
include("./procedures/ADF_reader.jl")
include("./procedures/xyz_reader.jl")
include("./procedures/ORCA-hess_reader.jl")
include("./procedures/TURBOMOLE_reader.jl")
include("./procedures/MOLDEN_reader.jl")

# Load essential modules:
using Interact, Colors, Luxor, Blink
using LinearAlgebra, Printf
using PeriodicTable
using ArgParse

#==============================================================================
            Argument parsing:
==============================================================================#

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "xyz"
            help = "Path to the .xyz file"
            arg_type = String
            #required = true
        "--orca"
            help = "Path to orca.hess file."
            arg_type = String
        "--adf"
            help = "ADF outfile with freq. calc. (no .xyz file is needed)"
            arg_type = String
        "--turbomole"
            help = "Path to vib_normal_modes, vibspectrum is expected to be in the same directory"
            arg_type = String
        "--molden"
            help = "normal modes in Molden format"
            arg_type = String
        "--nomw"
            help = "Do not use mass-weighting of normal modes in molden format"
            action = :store_false
    end
    return parse_args(s)
end

args = parse_commandline()

#===============================================================
                Define `Molecule` constructor:
===============================================================#

mutable struct Molecule
    #= Molecule constructor with all necessary properties =#
    xyzs::Array{Float64}
    atoms::Array{String}
    masses
    colors
    radii
    norm_modes
    freqs
    q_color
    q_num::Number
    mw::Bool
end

function Molecule(;
    #= Set default molecule properties =#
    xyzs = [],
    atoms = [],
    masses = map( x -> elements[Symbol(x)].atomic_mass.val, atoms),
    colors = map( x -> x == "H" ? "gray90" : elements[Symbol(x)].cpk_hex, atoms),
    radii = map( x -> define_radius(Symbol(x)), atoms),
    norm_modes = false,
    freqs = false,
    q_color = "#000000",
    q_num = 1,
    mw = false)
    return Molecule(xyzs, atoms, masses, colors, radii, norm_modes, freqs, q_color, q_num, mw)
end


#==============================================================================
            Function section:
==============================================================================#

function create_point(xyz, ϕ::Float64, θ::Float64)
    #= Return 2D coordinates of 3D [x,y,z] point using the orthographic projection.
        Point of View in polar coordinates has to be specified.
    =#
    pov = [ cosd(ϕ)*sind(θ), sind(ϕ)*sind(θ), cosd(θ) ]
    u = [ cosd(ϕ)*sind(θ + 90), sind(ϕ)*sind(θ + 90), cosd(θ + 90)]
    v = cross(u, pov)
    mat = Matrix{Float64}(undef,3,3)
    mat[:,1] = u
    mat[:,2] = v
    mat[:,3] = -pov
    sol = mat\xyz
    return sol[1:2]
end

function disp_coors(xyzs::Array, Cmat::Matrix{Float64}, q, q_scale::Number)
    #= Displace molecular coordinates in the direction of qth normal mode.
        Normal modes collected as rows of 3N×3N matrix C has to be specified.
        Coordinates are read as N×3 array.
    =#
    n = length(eachrow(xyzs))
    disp_vecs = reshape(Cmat[q,:], (3,n))'
    disp_points = xyzs + disp_vecs*2.2*q_scale
    return disp_points
end

function define_radius(atom::Symbol)
    #= Define relative radius of an atom for visualization =#
    mass = elements[atom].atomic_mass.val
    sign = mass < 12 ? -1 : 1
    radius = 1 + sign * log(mass)/12
    if atom == :H
        radius = 0.8
    end
    return radius
end

function check_heavy_atoms(atoms::Array{String})
    #= Function return a list of heavy atoms with larger bond distances =#
    heavy_atoms = String[]
    for atom in atoms
        atom_S = Symbol(atom)
        mass = elements[atom_S].atomic_mass.val
        if mass > 28
            push!(heavy_atoms, atom)
        end
    end
    return heavy_atoms
end

function draw_legend(atoms::Array)
    #= Draw legend in the top-left corner.
        Depends on size of Drawing in fucntions make_plot and make_plot2 =#
    types = unique(atoms)
    Ntypes = length(types)
    for (n, type) in enumerate(types)
        scale = define_radius(Symbol(type))
        setcolor( type == "H" ? "gray90" : elements[Symbol(type)].cpk_hex )
        circle(Point(350, -280 + 20*(n-1)), 4*scale, :fillpreserve)
        sethue("black")
        fontface("Sans")
        text(type, Point(370, -280 + 20*(n-1)), halign=:center, valign=:middle)
    end
    setline(2)
    box_center = Point(360, -280 + (Ntypes-1)*10)
    box( box_center, 50, (Ntypes+0.5)*20, 10, action=:stroke )
end

function make_plot(molecule::Molecule, r::Int, ϕ::Float64, θ::Float64, rotate::Float64,
    mode::String, hydrogens::Symbol; q_scale::Number=1.0)
    #= Creates a drawing with Luxor package, that is saved as "one_mol_vis.svg".
        Requires the Molecule object and scene meta-data.
        Point of View in polar coordinates is also required.
    =#
    hydrogens == :on ? (noHs = false) : (noHs = true)
    # Define Point-Of-View:
    pov = [ r*cosd(ϕ)*sind(θ), r*sind(ϕ)*sind(θ), r*cosd(θ) ]
    rotM = [[ cosd(rotate), -sind(rotate)];; [sind(rotate), cosd(rotate)]]
    # Initiate drawing:
    drawing = Drawing(800, 600, "one_mol_vis.svg")
    origin()
    
    # Prepare 2D atom coordinates:
    point_coors = map( p -> create_point(p*30, ϕ, θ), eachrow(molecule.xyzs))
    point_coors = map( p -> rotM'*p, point_coors)
    dists = map( x -> norm(x-pov), eachrow(molecule.xyzs))
    points = map(p -> Point(p...), point_coors*r/25) #.*(dists*1.5))
    
    # Draw a skelet from bonds:
    sethue("black")
    setline(2.5)
    heavy_atoms = check_heavy_atoms(molecule.atoms)
    for i in eachindex(points)
        for j in eachindex(points)
            if j >= i
                continue
            end
            d = norm(xyzs[i,:]-xyzs[j,:])
            if noHs && (atoms[i] == "H" || atoms[j] == "H")
                continue
            end
            if atoms[i] == "H" && atoms[j] == "H"
                continue
            end
            if d < 1.6
                line(points[i], points[j], :stroke)
            elseif (atoms[i] in heavy_atoms || atoms[j] in heavy_atoms) && d < 2.4
                line(points[i], points[j], :stroke)
            end
        end
    end
    # Visualization of normal modes:
    if typeof(molecule.norm_modes) != Bool
        q = molecule.q_num
        C = molecule.norm_modes
        if molecule.mw # Add mass-weighting?
            M = vcat( [[m,m,m] for m in molecule.masses ]... )
            M = sqrt.(M)
            M = diagm(M)
            C = C * M
            map!( x -> x/norm(x), eachrow(C), eachrow(C))
        end
        # Displace molecular geometry to get the arrow tips
        disps = disp_coors(xyzs, C, q, q_scale)*30
        norms = map( x -> norm(x), eachrow(reshape(C[q,:],(3,length(atoms)))') )
        arr_heads = map( p -> create_point(p, ϕ, θ), eachrow(disps))
        arr_heads = map( p -> rotM'*p, arr_heads)
        arr_heads = map( p -> Point(p...), arr_heads*r/25)
        for (i ,(p, f, cnorm)) in enumerate(zip(points, arr_heads, norms))
            if norm(p-f) < 1 || (atoms[i] == "H" && noHs)
                continue
            end
            setcolor(molecule.q_color)
            arrow(p, f, arrowheadlength=28*cnorm, linewidth=2.8)
        end
        setcolor("azure4")
        freq = @sprintf "%.4f" freqs[q]
        fontsize(14)
        text(freq, Point(0,-190), halign=:center, valign=:bottom)
    end

    # Define ordering and draw the molecule:
    vis_order = sortperm(dists)
    for i in vis_order
        if molecule.atoms[i] == "H" && noHs
            continue
        end
        scale = molecule.radii[i]
        setcolor( molecule.colors[i] )
        if mode == "Legended"
            circle(points[i],  4*scale, :fillpreserve)
        elseif typeof(molecule.norm_modes) != Bool
            setcolor("black")
            circle(points[i],  4*scale, :fill)
            fontsize(12)
            fontface("Sans")
            label(molecule.atoms[i], :NE, points[i], offset=10)
        else
            circle(points[i],  (12*dists[i])/r*scale, :fillpreserve)
            setline(2)
            sethue("black")
            strokepath()
            fontsize(10 + (r*0.2))
            fontface("Sans")
            text(molecule.atoms[i], points[i], halign=:center, valign=:middle)
        end
    end
    if mode == "Legended"
        draw_legend(atoms)
    end
    finish()
    preview()
    return drawing
end

#==============================================================================
            Main body of the programme:
==============================================================================#

# Prepare molecule:
if args["xyz"] != nothing
    (xyzs, atoms) = xyz_reader(args["xyz"])
    molecule = Molecule(xyzs=xyzs, atoms=atoms)
elseif args["molden"] != nothing
    (xyzs, atoms, freqs, C) = molden_reader(args["molden"])
    molecule = Molecule(xyzs=xyzs, atoms=atoms, norm_modes=C, freqs=freqs)
end

# Create widgets:
H_widget = widget([:on, :off], label="Show hydrogens:")

# Create an interactive object:
if molecule.norm_modes != Bool # Visualization of many normal modes
    if length(molecule.freqs) > 1
        one_mol = @manipulate for r in 10:40, ϕ in 0:0.1:360, θ in 0:0.1:360, rotate in 0:0.1:360,
            mode in ["Legended", "Labeled"], hydrogens in H_widget, q in 1:length(molecule.freqs),
            color in colorant"darkorange1", q_scale in 0.8:0.05:1.5
            
            molecule.q_color = color
            molecule.q_num = q
            molecule.mw = args["nomw"]
            make_plot( molecule, r, ϕ, θ, rotate, mode, hydrogens; q_scale=q_scale )
        end
    else # Visualization of one vibrational motion (e.g. Marcus dimension => in-house feature)
        one_mol = @manipulate for r in 10:40, ϕ in 0:0.1:360, θ in 0:0.1:360, rotate in 0:0.1:360,
            mode in ["Legended", "Labeled"], hydrogens in H_widget, color in colorant"darkorange1"
            
            molecule.q_color = color
            make_plot( molecule, r, ϕ, θ, rotate, mode, hydrogens )
        end
    end
else # Visualization of a plain molecule
    one_mol = @manipulate for r in 10:40, ϕ in 0:0.1:360, θ in 0:0.1:360, rotate in 0:0.1:360,
        mode in ["Legended", "Labeled"], hydrogens in H_widget
 
            make_plot( molecule, r, ϕ, θ, rotate, mode, hydrogens )
    end
end


# Open a window with visualization:
w = Window()
title(w::Window, "One molecule visualiation")
body!(w, one_mol)

# Keep window active:
while active(w)
    sleep(1)
end
