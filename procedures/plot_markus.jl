function read_markus_dimension(filepath::String, Natoms::Int)
    open(filepath) do file
        readline(file)
        readline(file)
        q_markus = Array{Float64}(undef, 0) 
        for line in eachline(file)
            dxyz = parse.(Float64, split(line)[2:end])
            append!(q_markus, dxyz)
        end
        return q_markus
    end
end



function plot_markus(xyzs::Array, atoms::Array, ϕ::Float64, θ::Float64, rotate::Float64, q::Array)
    #= Creates a drawing with Luxor package, that is saved as "one_mol_vis.svg".
        Atom types as Array{String} and xyz coordinates Matrix{Float64} has to be supplied.
        Point of View in polar coordinates is also required.
    =#
    pov = [ cosd(ϕ)*sind(θ), sind(ϕ)*sind(θ), cosd(θ) ]*20
    rotM = [[ cosd(rotate), -sind(rotate)];; [sind(rotate), cosd(rotate)]]
    # Initiate drawing:
    drawing = Drawing(800, 600, "markus_dimension.svg")
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
                setline(bond_thickness)
                line(points[i], points[j], :stroke)
            end
        end
    end
    n = length(atoms)
    disp_vecs = reshape(q, (3,n))'
    disps = (xyzs + disp_vecs*2.2)*30
    norms = map( x -> norm(x), eachrow(reshape(q, (3, n))') )
    arr_heads = map( p -> create_point(p, ϕ, θ), eachrow(disps))
    arr_heads = map( p -> rotM'*p, arr_heads)
    arr_heads = map( p -> Point(p...), arr_heads)
    for (i, f, cnorm) in zip(points, arr_heads, norms)
        if norm(i-f) < 1
            continue
        end
        #println(norm(i-f))
        setcolor("cadetblue3")
        arrow(i, f, arrowheadlength=22*cnorm, linewidth=2.8)
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
    finish()
    preview()
    return drawing
end

