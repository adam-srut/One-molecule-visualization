#! /usr/bin/env julia


function molden_reader(filename::String, NAtoms::Int)
    #= Reads molden input with arbitraru number of normal modes =#
    modes = []
    freqs = []
    open(filename, "r") do file
        while true
            line = readline(file)
            if startswith(line, "[FREQ]")
                while true
                    line = readline(file)
                    if line == "[FR-COORD]"
                        break
                    end
                    freq = parse(Float64, line)
                    append!(freqs, freq)
                end
            end
            if startswith(line, "vibration")
                mode = Array{Float64}(undef, NAtoms*3)
                for i_atom in 1:NAtoms
                    xyz = parse.(Float64, split(readline(file)))
                    xyz = xyz/1.8897259886
                    mode[ i_atom*3-2 : i_atom*3 ] = xyz
                end
                append!(modes, [mode])
            end
            if line == ""
                break
            end
        end
    end
    modes = vcat(map(x -> x', modes)...)'
    C_mat = zeros(6+length(freqs), NAtoms*3)
    freqs_pad = zeros(6+length(freqs))
    for (i,mode) in enumerate(eachcol(modes))
        C_mat[6+i,:] = mode'
    end
    for (i,freq) in enumerate(freqs)
        freqs_pad[i] = freq
    end
    return (freqs_pad, C_mat)
end


#(f,m,C_mat) = molden_reader("/work/Robin-Day/CTI/cti4/lvc/average_traj/pca-modes.molden", 52)
