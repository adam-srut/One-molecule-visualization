#! /usr/bin/env -S julia --project=one_mol_vis -J/home/adam/programs/one_mol_vis-dev/one_mol_vis/one_mol_vis.so



#==============================================================================
            Argument parsing:
==============================================================================#
using ArgParse
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "xyz"
            help = "Path to the .xyz file"
            arg_type = String
            default = "dummy"
            #required = true
        "--bond_thickness", "-b"
            help = "Specify bond thickness"
            arg_type = Int
            default = 4
        "--output_format", "-o"
            help = "Specify output format [svg::Default, png]"
            arg_type = String
            default = "svg"
        "--orca"
            help = "Path to orca.hess file."
            arg_type = String
        "--adf"
            help = "ADF outfile with freq. calc. (no .xyz file is needed)"
            arg_type = String
        "--turbomole"
            help = "Path to vib_normal_modes, vibspectrum is expected to be in the same directory"
            arg_type = String
        "--markusdim"
            help = "Path to markus-dimension.txt file"
            arg_type = String
        "--noHs"
            help = "Exlude hydrogens from visualization"
            action = :store_true
    end
    return parse_args(s)
end

# Prepare essential variables (path to .xyz file, path to NM file):
args = parse_commandline()
xyzfile = args["xyz"]
filename = split(xyzfile, '/')[end]
if filename[end-3:end] == ".xyz"
    basename = filename[1:end-4]
else
    println("File does not have an .xyz ending")
    println("Please check the file name.")
    basename = "one_mol_vis"
end
bond_thickness = args["bond_thickness"]
out_format = args["output_format"]
noHs = args["noHs"]

# Check what plotting routine is requested and set appropriete variables:
if args["adf"] != nothing       # Visualzation of norm. modes. with ADF (directly from .out file)
    adf_out = args["adf"]
    adf = true
    bond_thickness = 2
    norm_mode = true
elseif args["orca"] != nothing      # Visualzation of norm. modes. with ORCA
    bond_thickness = 2
    norm_mode = true
    NM_file = args["orca"]
    orca = true
    adf = false
elseif args["turbomole"] != nothing # Visualzation of norm. modes. with TURBOMOLE
    bond_thickness = 2
    norm_mode = true
    NM_file = args["turbomole"]
    dirpath = split(NM_file, '/')
    if length(dirpath) > 1
        freqfile = join(dirpath[1:end-1], '/') * '/' * "vibspectrum"
    else
        freqfile = "vibspectrum" # File vibspectrum has to be presented in working dir.
    end
    turbomole = true
    orca = false
    adf = false
elseif args["markusdim"] != nothing     # Visualization of Markus dimension (in house feature)
    norm_mode = true
    bond_thickness = 2
    adf = false
    orca = false
    turbomole = false
else    # Plain visualization of molecular geometry
    adf = false
    orca = false
    turbomole = false
    norm_mode = false
end


#==============================================================================
            Load procedures:
==============================================================================#

# Load file readers:
include("./procedures/ADF_reader.jl")
include("./procedures/xyz_reader.jl")
include("./procedures/ORCA-hess_reader.jl")
include("./procedures/TURBOMOLE_reader.jl")

# Load dictionaries with atomic types:
include("./atom_types.jl")

# Load essential modules:
using Interact, Colors, Luxor, Blink
using LinearAlgebra, Printf

#==============================================================================
            Function section:
==============================================================================#

function create_point(xyz, ??::Float64, ??::Float64)
    #= Return 2D coordinates of 3D [x,y,z] point using the orthographic projection.
        Point of View in polar coordinates has to be specified.
    =#
    pov = [ cosd(??)*sind(??), sind(??)*sind(??), cosd(??) ]
    u = [ cosd(??)*sind(?? + 90), sind(??)*sind(?? + 90), cosd(?? + 90)]
    v = cross(u, pov)
    mat = Matrix{Float64}(undef,3,3)
    mat[:,1] = u
    mat[:,2] = v
    mat[:,3] = -pov
    sol = mat\xyz
    return sol[1:2]
end

function disp_coors(xyzs::Array, Cmat::Matrix{Float64}, q::Int)
    #= Displace molecular coordinates in the direction of qth normal mode.
        Normal modes collected as rows of 3N??3N matrix C has to be specified.
        Coordinates are read as N??3 array.
    =#
    n = length(eachrow(xyzs))
    disp_vecs = reshape(Cmat[6+q,:], (3,n))'
    disp_points = xyzs + disp_vecs*2.2
    return disp_points
end


function make_plot(xyzs::Array, atoms::Array, r::Int, ??::Float64, ??::Float64, rotate::Float64)
    #= Creates a drawing with Luxor package, that is saved as "one_mol_vis.svg".
        Atom types as Array{String} and xyz coordinates Matrix{Float64} has to be supplied.
        Point of View in polar coordinates is also required.
        Change to .png is possible, .pdf exhibits certain issues.
    =#
    global bond_thickness
    global out_format
    pov = [ r*cosd(??)*sind(??), r*sind(??)*sind(??), r*cosd(??) ]
    rotM = [[ cosd(rotate), -sind(rotate)];; [sind(rotate), cosd(rotate)]]
    # Initiate drawing:
    drawing = Drawing(800, 600, "$basename.$out_format")
    origin()
    # Prepare points coordinates:
    point_coors = map( p -> create_point(p*50, ??, ??), eachrow(xyzs))
    point_coors = map( p -> rotM'*p, point_coors)
    dists = map( x -> norm(x-pov), eachrow(xyzs))
    points = map(p -> Point(p...), point_coors*r/25) #.*(dists*1.5))
    # Draw a skelet from bonds:
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
            if d < 1.5
                sethue("black")
                setline(bond_thickness)
                line(points[i], points[j], :stroke)
            end
        end
    end
    # Order atoms by their distatce to pov and plot as labeled circles
    to_plot = map( (atom, point, dist) -> (atom, point, dist), atoms, points, dists)
    sort!(to_plot, by = x -> x[3])
    for atom in to_plot
        name = atom[1]
        if name == "H" && noHs
            continue
        end
        scale = radii[name]
        setcolor(colors[name])
        circle(atom[2],  (0.8*atom[3])^2/r*scale, :fillpreserve)
        setline(2)
        sethue("black")
        strokepath()
        fontsize(10 + (r*0.2))
        fontface("Sans")
        text(name, atom[2], halign=:center, valign=:middle)
    end
    finish()
    preview()
    return drawing
end


function make_plot2(xyzs::Array, atoms::Array, ??::Float64, ??::Float64, rotate::Float64, q::Int)
    #= Creates a drawing with Luxor package, that is saved as "one_mol_vis.svg".
        Atom types as Array{String} and xyz coordinates Matrix{Float64} has to be supplied.
        Point of View in polar coordinates is also required.
        This function serves for visualization of normal modes!
    =#
    global noHs
    pov = [ cosd(??)*sind(??), sind(??)*sind(??), cosd(??) ]*20
    rotM = [[ cosd(rotate), -sind(rotate)];; [sind(rotate), cosd(rotate)]]
    # Initiate drawing:
    drawing = Drawing(800, 600, "$basename.$out_format")
    origin()
    # Prepare points coordinates:
    point_coors = map( p -> create_point(p*30, ??, ??), eachrow(xyzs))
    point_coors = map( p -> rotM'*p, point_coors)
    dists = map( x -> norm(x-pov), eachrow(xyzs))
    points = map(p -> Point(p...), point_coors)
    # Draw a skelet from bonds:
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
            if d < 1.65
                sethue("black")
                setline(bond_thickness)
                line(points[i], points[j], :stroke)
            end
        end
    end
    disps = disp_coors(xyzs, C, q)*30
    norms = map( x -> norm(x), eachrow(reshape(C[q+6,:],(3,length(atoms)))') )
    arr_heads = map( p -> create_point(p, ??, ??), eachrow(disps))
    arr_heads = map( p -> rotM'*p, arr_heads)
    arr_heads = map( p -> Point(p...), arr_heads)
    for (i ,(p, f, cnorm)) in enumerate(zip(points, arr_heads, norms))
        if norm(p-f) < 1 || (atoms[i] == "H" && noHs)
            continue
        end
        setcolor("cadetblue3")
        arrow(p, f, arrowheadlength=22*cnorm, linewidth=2.8)
    end
    # Order atoms by their distatce to pov and plot as labeled circles
    to_plot = map( (atom, point, dist) -> (atom, point, dist), atoms, points, dists)
    sort!(to_plot, by = x -> x[3])
    for atom in to_plot
        name = atom[1]
        if noHs && name == "H"
            continue
        end
        setcolor("black")
        circle(atom[2],  4, :fill)
        fontsize(12)
        fontface("Sans")
        label(name, :NE, atom[2], offset=10)
    end
    setcolor("azure4")
    freq = @sprintf "%.2f cm?????" freqs[q+6]
    fontsize(14)
    text(freq, Point(0,100), halign=:center, valign=:bottom)
    finish()
    preview()
    return drawing
end
#==============================================================================
            Main body of the programme:
==============================================================================#

# Read input parameters:
if adf
    xyzs, atoms, freqs, C = ADF_reader(adf_out)
else
    xyzs, atoms = xyz_reader(xyzfile)

    if orca
        (freqs, C) = orcaHess_reader(NM_file)
    elseif turbomole
        C = TURBOMOLE_reader(NM_file, length(atoms)*3)
        freqs = turbofreq_read(freqfile)
    end
end


# Create an interactive object:
if ! norm_mode # Plain visualization of molecular geometry
    one_mol = @manipulate for r in 10:40, ?? in 0:0.1:360, ?? in 0:0.1:360, rotate in 0:0.1:360

        make_plot( xyzs, atoms, r, ??, ??, rotate )

    end
elseif args["markusdim"] != nothing # Visualization of Marcus dimension
    include("./procedures/plot_markus.jl")
    mfile = args["markusdim"]
    basename = mfile[1:end-4]
    q_markus = read_markus_dimension(args["markusdim"], length(atoms))
    one_mol = @manipulate for ?? in 0:0.1:360, ?? in 0:0.1:360, rotate in 0:0.1:360, q_scale in 0.5:0.1:2.0, labels in ["atom", "index"]
        plot_markus(xyzs, atoms, ??, ??, rotate, q_markus, q_scale, basename, labels)
    end
else
    one_mol = @manipulate for # Visualization of normal modes
             ?? in 0:0.1:360, ?? in 0:0.1:360, rotate in 0:0.1:360, q in 1:(length(atoms)*3-6)
        make_plot2( xyzs, atoms, ??, ??, rotate, q )

    end
end

# Open a window with visualization:
w = Window()
title(w::Window, "One molecule visualiation")
body!(w, one_mol)

while active(w)
    sleep(1)
end
