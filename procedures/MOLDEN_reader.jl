#! /usr/bin/env julia


function molden_reader(filename::String)
    #= Reads molden input with arbitrary number of normal modes =#
    Ang_to_Bohr = 1.8897259886
    modes = []
    freqs = Float64[]
    atoms = String[]
    xyzs = []
    NAtoms = 0
    open(filename, "r") do file
        line = readline(file)
        while true
            if line == "[FR-COORD]"
                while true
                    line = readline(file)
                    xyzline = split(line)
                    if length(xyzline) < 4
                        break
                    end
                    atom = xyzline[1]
                    atom = uppercasefirst(atom)
                    push!(atoms, atom)
                    NAtoms += 1
                    xyz = parse.(Float64, xyzline[2:end])
                    xyz = xyz/Ang_to_Bohr
                    push!(xyzs, xyz)
                end
            elseif startswith(line, "[FREQ]")
                while true
                    line = readline(file)
                    if startswith(line, "[")
                        break
                    end
                    freq = parse(Float64, line)
                    append!(freqs, freq)
                end
            elseif startswith(strip(line), "vibration")
                mode = Array{Float64}(undef, NAtoms*3)
                for i_atom in 1:NAtoms
                    line = readline(file)
                    xyz = parse.(Float64, split(line))
                    xyz = xyz/Ang_to_Bohr
                    mode[ i_atom*3-2 : i_atom*3 ] = xyz
                end
                append!(modes, [mode])
            elseif line == ""
                break
            else
                line = readline(file)
            end
        end
    end
    xyzs = vcat(map(x -> x', xyzs)...)
    cog = sum(eachrow(xyzs))/NAtoms 
    xyzs = xyzs .- cog' # Shift coordinates to the geometric centre
    modes = vcat(map(x -> x', modes)...)'
    # Remove padding with zero frequencies if present:
    if length(freqs) > 6
        if freqs[1:6] == zeros(6)
            freqs = freqs[7:end]
            C_mat = zeros(length(freqs), NAtoms*3)
            for (i,mode) in enumerate(eachcol(modes))
                if i <= 6
                    continue
                end
                C_mat[i-6,:] = mode'
           end
        else
            C_mat = convert(Matrix, modes')
        end
    else
        C_mat = convert(Matrix, modes')
    end
    return (xyzs, atoms, freqs, C_mat)
end


#(xyzs, atoms, f, C_mat) = molden_reader("/work/Robin-Day/CTI/cti4/lh20t/d4disp/molden.input")
