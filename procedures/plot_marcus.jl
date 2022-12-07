function read_marcus_dimension(filepath::String, Natoms::Int)
    open(filepath) do file
        readline(file)
        readline(file)
        q_marcus = Array{Float64}(undef, 0) 
        for line in eachline(file)
            dxyz = parse.(Float64, split(line)[2:end])
            append!(q_marcus, dxyz)
        end
        return q_marcus
    end
end



function plot_marcus(xyzs::Array, atoms::Array, ϕ::Float64, θ::Float64, rotate::Float64, q::Array, scale::Float64, name::String, 
        labels::String, color, mode::String)
    #= Creates a drawing with Luxor package, that is saved as "one_mol_vis.svg".
        Atom types as Array{String} and xyz coordinates Matrix{Float64} has to be supplied.
        Point of View in polar coordinates is also required.
    =#
    pov = [ cosd(ϕ)*sind(θ), sind(ϕ)*sind(θ), cosd(θ) ]*20
    rotM = [[ cosd(rotate), -sind(rotate)];; [sind(rotate), cosd(rotate)]]
    # Initiate drawing:
    drawing = Drawing(800, 600, "$name.svg")
    origin()
    # Prepare points coordinates:
    point_coors = map( p -> create_point(p*30, ϕ, θ), eachrow(xyzs))
    point_coors = map( p -> rotM'*p, point_coors)
    dists = map( x -> norm(x-pov), eachrow(xyzs))
    points = map(p -> Point(p...), point_coors)
    # Different bonding distaces for heavy-atoms:
    heavy_atoms = ["Ru", "Re", "Fe"]
    sethue("black")
    setline(bond_thickness)
    # Draw a skelet from bonds:
    for i in eachindex(points)
        for j in eachindex(points)
            if j >= i
                continue
            end
            if noHs && (atoms[i] == "H" || atoms[j] == "H")
                continue
            end
            d = norm(xyzs[i,:]-xyzs[j,:])
            if atoms[i] == "H" && atoms[j] == "H"
                continue
            end
            if d < 1.5
                line(points[i], points[j], :stroke)
            elseif (atoms[i] in heavy_atoms || atoms[j] in heavy_atoms) && d < 2.5
                line(points[i], points[j], :stroke)               
            end
        end
    end
    n = length(atoms)
    indices = [i for (i,atom) in enumerate(atoms)]
    disp_vecs = reshape(q*scale, (3,n))'
    disps = (xyzs + disp_vecs*2.2)*30
    norms = map( x -> norm(x), eachrow(reshape(q*scale, (3, n))') )
    arr_heads = map( p -> create_point(p, ϕ, θ), eachrow(disps))
    arr_heads = map( p -> rotM'*p, arr_heads)
    arr_heads = map( p -> Point(p...), arr_heads)
    for (i ,(p, f, cnorm)) in enumerate(zip(points, arr_heads, norms))
        if norm(p-f) < 1 || (atoms[i] == "H" && noHs)
            continue
        end
        setcolor(color)
        arrow(p, f, arrowheadlength=22*cnorm, linewidth=2.8)
    end
    # Order atoms by their distatce to pov and plot as labeled circles
    to_plot = map( (atom, point, dist, index) -> (atom, point, dist, index), atoms, points, dists, indices)
    sort!(to_plot, by = x -> x[3])
    for (i, atom) in enumerate(to_plot)
        atom_scale = radii[atom[1]]
        if noHs && atom[1] == "H"
            continue
        end
        if mode == "Legended"
            setcolor(colors[atom[1]])
            circle(atom[2],  4*atom_scale, :fillpreserve)
        else
            if labels == "atom"
                name = atom[1]
            else
                name = string(atom[4])
            end
            setcolor("black")
            circle(atom[2],  4*atom_scale, :fill)
            fontsize(12)
            fontface("Sans")
            label(name, :NE, atom[2], offset=10)
        end
    end
    if mode == "Legended"
        draw_legend(atoms)
    end
    finish()
    preview()
    return drawing
end

