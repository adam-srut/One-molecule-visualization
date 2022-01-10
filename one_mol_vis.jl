#! /usr/bin/env julia

colors = Dict(
    "N" => "blue",
    "O" => "red",
    "H" => "gray",
    "C" => "cyan",
	"Cl" => "yellow"
)
radii = Dict(
    "N" => 1.1,
    "O" => 1.1,
    "H" => 0.8,
    "C" => 1.0,
	"Cl" => 1.5
)

using ArgParse
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--xyz", "-f"
            help = "Path to the .xyz file"
            arg_type = String
            required = true
        "--bond_thickness", "-b"
            help = "Specify bond thickness"
            arg_type = Int
            default = 8
        "--output_format", "-o"
            help = "Specify output format"
            arg_type = String
            default = "svg"
    end
    return parse_args(s)
end

args = parse_commandline()
xyzfile = args["xyz"]
bond_thickness = args["bond_thickness"]
out_format = args["output_format"]

using Interact, Colors, Luxor, Blink
using LinearAlgebra

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

function create_point(xyz, r::Int, ϕ::Float64, θ::Float64)
    #= Return point coordinates of [x,y,z] point in orthographic projection.
        Point of View in polar coordinates has to be specified.
    =#
    pov = [ cosd(ϕ)*sind(θ), sind(ϕ)*sind(θ), cosd(θ) ]
    u = [ cosd(ϕ)*sind(θ + 90), sind(ϕ)*sind(θ + 90), cosd(θ + 90)]
    v = cross(u, pov)
    mat = Matrix{Float64}(undef,3,3)
    #[mat[i,1] = u[i] for i in 1:3]
    mat[:,1] = u
    mat[:,2] = v
    mat[:,3] = -pov
    #[mat[i,2] = v[i] for i in 1:3]
    #[mat[i,3] = -pov[i] for i in 1:3]
    sol = mat\xyz
    return sol[1:2]
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
    drawing = Drawing(600, 600, "one_mol_vis." * out_format)
    origin()
    # Prepare points coordinates:
    point_coors = map( p -> create_point(p,r,ϕ,θ), eachrow(xyzs))
    point_coors = map( p -> rotM'*p, point_coors)
    dists = map( x -> norm(x-pov), eachrow(xyzs))
    points = map(p -> Point(p...), point_coors.*(dists*1.5))
    # Draw a skelet from bond:
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

#==============================================================================
            Main body of the programme:
==============================================================================#

xyzs, atoms = xyz_reader(xyzfile)

one_mol = @manipulate for r in 10:40, ϕ in 0:0.1:360, θ in 0:0.1:360, rotate in 0:0.1:360

    make_plot( xyzs, atoms, r, ϕ, θ, rotate )

end


w = Window()
body!(w, one_mol)

while active(w)
    sleep(1)
end
