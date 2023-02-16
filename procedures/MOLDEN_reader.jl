#! /usr/bin/env julia


function molden_reader(filename::String)
    #= Reads molden input with arbitrary number of normal modes =#
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
                    xyz = xyz/1.8897259886
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
                    xyz = xyz/1.8897259886
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
    modes = vcat(map(x -> x', modes)...)'
    if freqs[1:6] == zeros(6)
        C_mat = convert(Matrix, modes')
    else
        C_mat = zeros(6+length(freqs), NAtoms*3)
        freqs_pad = zeros(6+length(freqs))
        for (i,mode) in enumerate(eachcol(modes))
            C_mat[6+i,:] = mode'
        end
        for (i,freq) in enumerate(freqs)
            freqs_pad[i+6] = freq
        end
        freqs = freqs_pad
    end
    return (xyzs, atoms, freqs, C_mat, modes)
end


#(xyzs, atoms, f, C_mat, modes) = molden_reader("/work/Robin-Day/CTI/cti4/lh20t/d4disp/molden.input")
