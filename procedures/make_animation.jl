#!/usr/bin/env julia

using Interact, Colors, Luxor
using LinearAlgebra, Printf

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

function disp_coors(xyzs::Array, q::Array)
    #= Displace molecular coordinates in the direction of qth normal mode.
        Normal modes collected as rows of 3N×3N matrix C has to be specified.
        Coordinates are read as N×3 array.
    =#
    n = length(eachrow(xyzs))
    disp_vecs = reshape(q, (3,n))'
    disp_points = xyzs + disp_vecs*2.2
    return disp_points
end

function make_plot2(xyzs::Array, atoms::Array, ϕ::Float64, θ::Float64, rotate::Float64, q::Array)
    #= Creates a drawing with Luxor package, that is saved as "one_mol_vis.svg".
        Atom types as Array{String} and xyz coordinates Matrix{Float64} has to be supplied.
        Point of View in polar coordinates is also required.
    =#
    pov = [ cosd(ϕ)*sind(θ), sind(ϕ)*sind(θ), cosd(θ) ]*20
    rotM = [[ cosd(rotate), -sind(rotate)];; [sind(rotate), cosd(rotate)]]
    # Initiate drawing:
    origin()
    # Prepare points coordinates:
    point_coors = map( p -> create_point(p*30, ϕ, θ), eachrow(xyzs))
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
                setline(2)
                line(points[i], points[j], :stroke)
            end
        end
    end
    disps = disp_coors(xyzs, q)*30
    norms = map( x -> norm(x), eachrow(reshape(q,(3,length(atoms)))') )
    arr_heads = map( p -> create_point(p, ϕ, θ), eachrow(disps))
    arr_heads = map( p -> rotM'*p, arr_heads)
    arr_heads = map( p -> Point(p...), arr_heads)
    for (i, f, cnorm) in zip(points, arr_heads, norms)
        if norm(i-f) < 1
            continue
        end
        #println(norm(i-f))
        setcolor("cadetblue3")
        Luxor.arrow(i, f, arrowheadlength=22*cnorm, linewidth=2.8)
    end
    # Order atoms by their distatce to pov and plot as labeled circles
    to_plot = map( (atom, point, dist) -> (atom, point, dist), atoms, points, dists)
    sort!(to_plot, by = x -> x[3])
    for atom in to_plot
        name = atom[1]
        setcolor("black")
        circle(atom[2],  4, :fill)
        fontsize(12)
        fontface("Sans")
        label(name, :NE, atom[2], offset=10)
    end
    setcolor("azure4")
    fontsize(14)
    #finish()
    #preview()
end


function make_animation(xyzs, atoms, q, name)
    function frame(scene, framenumber)
        make_plot2(xyzs, atoms, 90.0, framenumber, 0.0, q)
    end
    anim = Movie(500, 500, "./markus_dimension-$name")
    Luxor.animate(anim, Scene( anim, frame, 1.0:360.0), creategif=true,pathname="./MDLCNM/markus_dimension-$name.gif")
end


