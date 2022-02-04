#! /usr/bin/env julia


function TURBOMOLE_reader(filepath::String, dof::Int)
    i_dof = 1
    C = Matrix{Float64}(undef, dof, dof)
    i_mode = 1
    for line in eachline(filepath)
        if startswith(line, "\$")
            continue
        end
        line = split(line)
        modenum = parse(Int, line[1])
        if i_mode == modenum
            vals = parse.(Float64, line[3:end])
            for (i, x) in enumerate(vals)
                C[ i + (i_dof-1)*5, i_mode ] = x
            end
            i_dof += 1
        else
            i_mode += 1
            i_dof = 1
            vals = parse.(Float64, line[3:end])
            for (i, x) in enumerate(vals)
                C[ i + (i_dof-1)*5, i_mode ] = x
            end
            i_dof += 1
        end
    end
    return C 
end

function turbofreq_read(filepath::String)
    freqs = []
    for line in eachline(filepath)
        if startswith(line, '#') || startswith(line, '$')
            continue
        end
        line = split(line)
        freq = parse(Float64, line[3])
        append!(freqs, freq)
    end
    return freqs
end



if abspath(PROGRAM_FILE) == @__FILE__
    f = "/home/adam/Documents/turbo-test/vib_normal_modes"
    C = TURBOMOLE_reader(f,48)
    ff = "/home/adam/Documents/turbo-test/vibspectrum"
    freq = turbofreq_read(ff)
end

