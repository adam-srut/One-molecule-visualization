#!/usr/bin/env julia


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

