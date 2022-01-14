#! /usr/bin/env julia

colors = Dict(
    "N" => "blue",
    "O" => "red",
    "H" => "gray",
    "C" => "cyan",
	"Cl" => "yellow",
	"S" => "orange"
)
radii = Dict(
    "N" => 1.1,	"O" => 1.1,	 "H" => 0.8,
	"C" => 1.0,	"Cl" => 1.5, "S" => 1.2 
)

using ArgParse
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "xyz"
            help = "Path to the .xyz file"
            arg_type = String
            required = true
        "--bond_thickness", "-b"
            help = "Specify bond thickness"
            arg_type = Int
            default = 4
        "--output_format", "-o"
			help = "Specify output format [svg::Default, png]"
            arg_type = String
            default = "svg"
		"--normal_modes", "-n"
			help = "Path to orca.hess file."
			arg_type = String
    end
    return parse_args(s)
end

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
if args["normal_modes"] == nothing
	norm_mode = false
else
	norm_mode = true
	NM_file = args["normal_modes"]
	bond_thickness = 2
end


using Interact, Colors, Luxor, Blink
using LinearAlgebra, Printf

#==============================================================================
            Function section:
==============================================================================#

function xyz_reader(filename::String)
	#= XYZ reader. Takes path to .xyz file as an input.
	Returns tuple ( coors, atom ):
		- coors is a N×3 matrix of atomic coordinates.
		- atoms is a N vector with atomic types.
	=#
    open(filename) do file
        n = parse(Int,readline(file))
        coors = Matrix{Float64}(undef,n,3)
        atoms = Array{String}(undef,n)
        readline(file)
        for i_atom in 1:n
            line = split(readline(file))
            atom = line[1]
            xyz = parse.(Float64, line[2:end])
            coors[i_atom,:] = xyz
            atoms[i_atom] = atom
        end
        cog = sum(eachrow(coors))/n
        coors = coors .- cog' # Shift to geometric centre
        return coors, atoms
    end
end

function orcaHess_reader(filename::String)
	#= Reads .hess file from ORCA calculation.
		Returns matrix C with normal modes as rows (including translation and rotation).
	=#
    open(filename) do file
        while true
            line = readline(file)
            if startswith(line , "\$vibrational_frequencies")
                break
            end
        end
        dof = parse(Int, readline(file))
        freqs = Array{Float64}(undef, dof)
        for i in 1:dof
            line =  split(readline(file))
            freq = parse(Float64, line[2])
            freqs[i] = freq
        end
		while true
            line = readline(file)
            if startswith(line, "\$normal_modes")
                break
            end
        end
		readline(file)
        #dof = parse(Int, split(readline(file))[1])
        C = Matrix{Float64}(undef, dof, dof)
        i_block = 1
        while i_block <= fld(dof, 5) + 1
            readline(file)
            i_dof = 1
            while i_dof <= dof
                c_ij = parse.(Float64, split(readline(file))[2:end])
                for (j, cij) in enumerate(c_ij)
                    C[ (i_block-1)*5 + j, i_dof ] = cij
                end
                i_dof += 1
            end
            i_block += 1
        end
		return (freqs, C)
    end
end

function create_point(xyz, ϕ::Float64, θ::Float64)
    #= Return point coordinates of [x,y,z] point in orthographic projection.
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

function disp_coors(xyzs::Array, Cmat::Matrix{Float64}, q::Int)
	n = length(eachrow(xyzs))
	disp_vecs = reshape(Cmat[6+q,:], (3,n))'
	disp_points = xyzs + disp_vecs*2.2
	return disp_points
end


function make_plot(xyzs::Array, atoms::Array, r::Int, ϕ::Float64, θ::Float64, rotate::Float64)
    #= Creates a drawing with Luxor package, that is saved as "one_mol_vis.svg".
        Atom types as Array{String} and xyz coordinates Matrix{Float64} has to be supplied.
        Point of View in polar coordinates is also required.
        Change to .png is possible, .pdf exhibits certain issues.
    =#
    global bond_thickness
    global out_format
    pov = [ r*cosd(ϕ)*sind(θ), r*sind(ϕ)*sind(θ), r*cosd(θ) ]
    rotM = [[ cosd(rotate), -sind(rotate)];; [sind(rotate), cosd(rotate)]]
    # Initiate drawing:
    drawing = Drawing(600, 600, "$basename.$out_format")
    origin()
    # Prepare points coordinates:
    point_coors = map( p -> create_point(p*50, ϕ, θ), eachrow(xyzs))
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
        scale = radii[name]
        setcolor(colors[name])
        circle(atom[2],  (0.8*atom[3])^2/r*scale, :fillpreserve)
        setline(2)
        sethue("black")
        strokepath()
        fontsize(12 + (r*0.2))
        fontface("Sans")
        text(name, atom[2], halign=:center, valign=:middle)
    end
    finish()
    preview()
    return drawing
end

function make_plot2(xyzs::Array, atoms::Array, ϕ::Float64, θ::Float64, rotate::Float64, q::Int)
    #= Creates a drawing with Luxor package, that is saved as "one_mol_vis.svg".
        Atom types as Array{String} and xyz coordinates Matrix{Float64} has to be supplied.
        Point of View in polar coordinates is also required.
    =#
    pov = [ cosd(ϕ)*sind(θ), sind(ϕ)*sind(θ), cosd(θ) ]*20
    rotM = [[ cosd(rotate), -sind(rotate)];; [sind(rotate), cosd(rotate)]]
    # Initiate drawing:
    drawing = Drawing(400, 400, "$basename.$out_format")
    origin()
    # Prepare points coordinates:
    point_coors = map( p -> create_point(p*40, ϕ, θ), eachrow(xyzs))
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
        setcolor("black")
        circle(atom[2],  2, :fill)
        fontsize(12)
        fontface("Sans")
        label(name, :NE, atom[2], offset=10)
    end
	disps = disp_coors(xyzs, C, q)*40
	dists2 = map( x -> norm(x), eachrow(reshape(C[q+6,:],(3,length(atoms)))') )
	arr_heads = map( p -> create_point(p, ϕ, θ), eachrow(disps))
	arr_heads = map( p -> rotM'*p, arr_heads)
	#arr_heads = map( (p, d) -> p*(40), arr_heads, dists2)
	arr_heads = map( p -> Point(p...), arr_heads)
	for (i, f, cnorm) in zip(points, arr_heads, dists2)
		setcolor("azure3")
		arrow(i, f, arrowheadlength=22*cnorm, linewidth=2)
	end
	setcolor("azure3")
	freq = @sprintf "%.2f cm⁻¹" freqs[q+6]
	fontsize(14)
	text(freq, Point(0,180), halign=:center, valign=:bottom)
    finish()
    preview()
    return drawing
end
#==============================================================================
            Main body of the programme:
==============================================================================#

xyzs, atoms = xyz_reader(xyzfile)

if norm_mode
	(freqs, C) = orcaHess_reader(NM_file)
end

# Create an interactive object:
if ! norm_mode
	one_mol = @manipulate for r in 10:40, ϕ in 0:0.1:360, θ in 0:0.1:360, rotate in 0:0.1:360

    	make_plot( xyzs, atoms, r, ϕ, θ, rotate )

	end
else
	one_mol = @manipulate for 
			 ϕ in 0:0.1:360, θ in 0:0.1:360, rotate in 0:0.1:360, q in 1:(length(atoms)*3-6)

        make_plot2( xyzs, atoms, ϕ, θ, rotate, q )

    end
end

# Open a window with visualization:
w = Window()
title(w::Window, "One molecule visualiation")
body!(w, one_mol)

while active(w)
    sleep(1)
end
