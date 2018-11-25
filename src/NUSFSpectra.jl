module NUSFSpectra

using Dates, PyCall, JLD2, FileIO, DelimitedFiles

import Base: *, -, +, /

function __init__()
    push!(pyimport("sys")["path"], joinpath(@__DIR__, "../pythonmodules/pySfgProcess/"))
    global dfg = pyimport("dfg")
    global pscalib = pyimport("pscalib")
    global spectrum = pyimport("spectrum")
    global winspec = pyimport("winspec")
end


mutable struct NUSFSpectrum{T<:Number,N} <: AbstractArray{T,N}
    s::Array{T,N}
    wavenumbers::Array{T,1}
    pyobjs::Array{PyObject,1}
    name::String
    dates::Array{DateTime,1}

    """
    Inner constructor. Access the PyObjects data to put in wavelength and signal
    arrays directly into the object.
    """
    function NUSFSpectrum(pyobjs, name, dates)
        firstspectrum = pyobjs[1][:dfgTruncatedSum]
        spectrumlength = length(firstspectrum)
        if length(pyobjs) > 1
            s = Array{eltype(firstspectrum), 2}(undef, spectrumlength, length(pyobjs))
        elseif length(pyobjs) == 1
            s = Array{eltype(firstspectrum), 1}(undef, spectrumlength)
        else
            error("2D Spectra not supported.")
        end

        s[:,1] = firstspectrum
        for i = 2:length(pyobjs)
            s[:,i] = pyobjs[i][:dfgTruncatedSum]
        end

        wavenumbers = pyobjs[1][:fullwn]

        new{eltype(s), ndims(s)}(s, wavenumbers, pyobjs, name, dates)
    end
end

Base.size(s::NUSFSpectrum) = size(s.s)
Base.getindex(s::NUSFSpectrum, i::Int) = getindex(s.s, i)
Base.getindex(s::NUSFSpectrum{T,N}, I::Vararg{Int, N}) where {N,T} = getindex(s.s, I...)
Base.setindex!(s::NUSFSpectrum, v::Number, i::Int) = setindex!(s.s, v, i)
Base.setindex!(s::NUSFSpectrum{T,N}, v::Number, I::Vararg{Int, N}) where {N,T} = setindex!(s.s, v, I...)
Base.copy(s::NUSFSpectrum) = NUSFSpectrum(copy(s.id), copy(s.s))
Base.stride(s::NUSFSpectrum, k) = size(s, 1)
Base.strides(s::NUSFSpectrum) = (1, size(s, 1))
Base.convert(T::Array{Number}, s::NUSFSpectrum) = s.s |> T

+(s::NUSFSpectrum, t::NUSFSpectrum) = s.s + t.s
-(s::NUSFSpectrum, t::NUSFSpectrum) = s.s - t.s
*(s::NUSFSpectrum, a::Number)     = s.s * a
*(a::Number, s::NUSFSpectrum)     = a * s.s
/(s::NUSFSpectrum, a::Number)     = s.s / a
/(a::Number, s::NUSFSpectrum)     = a ./ s.s



function calibrate(filepath, peak1=(2825, 2860), peak2=(3045, 3070))
    calibfile = filepath

    ps = pscalib[:PScalib](calibfile)
    ps[:plot]()

    #starting and ending values of segment to fit
    val1 = peak1[1]
    val2 = peak1[2]

    #fit it with gaussian, calculate shift
    ps[:fitPeak](val1,val2,0)

    #starting and ending values of segment to fit
    val1 = peak2[1]
    val2 = peak2[2]

    #fit it with gaussian, calculate shift
    ps[:fitPeak](val1,val2,1)

    calib = ps[:evaluateShift]()
    println("Shift: $calib wavenumbers")

    path = splitdir(calibfile)[1]

    save(joinpath(path, "calibration.jld2"), "calib", calib)
    println("""Saved calibration to $(joinpath(path, "calibration.jld2"))""")

    return calib
end


function batchprocess(path, calibfile=joinpath(path, "calibration.jld2");
                        trunc=1e-6, doplot=false,
                        savefile=joinpath(path, "data.jld2"))
    (root, dirs, files) = first(walkdir(path))
    filter!(x -> !startswith(x, "."), dirs)
    paths = joinpath.(root, dirs)

    subpaths = Array{String}[]
    for i in 1:length(paths)
        root, dirs, files = first(walkdir(paths[i]))
        push!(subpaths, joinpath.(paths[i], dirs))
    end

    specs = NUSFSpectrum[]
    for p in subpaths
        numruns = length(p)
        dates = Array{DateTime}(undef, numruns)
        pyobjs = Array{PyObject}(undef, numruns)
        spectra = Array{Array{Float64,1},1}(undef, numruns)
        for (i,s) in enumerate(p)
            # p are the arrays of subpaths, so s are the individual runs
            specinfo = processgold(s, calibfile, doplot=doplot, trunc=trunc)
            dates[i] = specinfo[2][1]
            pyobjs[i] = specinfo[3]
        end

        name = split(p[1], "/")[end-1]

        ## Order spectra by date
        i = sortperm(dates)
        pyobjs = pyobjs[i]
        dates = dates[i]

        push!(specs, NUSFSpectrum(pyobjs, name, dates))
    end

    return specs

    splitext(savefile)[end] != ".jld2" && (savefile *= ".jld2")

    save(savefile, "spectra", specs)

    @info "Saved the following spectra to $savefile:"
    for n in getfield.(specs, :name) println(n) end

    return specs
end


function processgold(path, calibrationfile; trunc=1e-6, doplot=false)

    calib = load(calibrationfile, "calib")
    path
    # get all the .SPE files
    filenames = filter!(x -> splitext(x)[2] == ".SPE", readdir(path))
    filepaths = joinpath.(path, filenames)

    #get datetime of the files
    datetimes = Array{DateTime}(undef, length(filepaths))
    for i = 1:length(filepaths)
        f = winspec[:SpeFile](filepaths[i])
        datetimestr = f[:header][:date] * f[:header][:ExperimentTimeLocal]
        datetimes[i] = DateTime(datetimestr, dateformat"dduuuyyyyHHMMSS")
    end

    #create object, loads each sample and background DFG
    gold = spectrum[:Spectrum](path, shift=calib)

    #subtract the appropriate background DFG from each sample DFG
    gold[:subtractBGs]()

    #plot the imported sample DFGs
    doplot && gold[:plotDFGs]()

    #plot the imported background DFGs
    doplot && gold[:plotBGs]()

    #pad the dfgs with zeros so they align and can be summed up
    gold[:padDFGs]()

    #find the indices of where the reference signal falls off to guide truncation
    gold[:findTruncateIndices](trunc)

    #truncate
    gold[:truncateFullDFGs](gold)

    #plot the truncated DFGs
    doplot && gold[:plotTruncatedDFGs]()

    #sum the truncated DFGs
    gold[:sumTruncatedDFGs]()

    #plot the summed, truncated DFGS
    doplot && gold[:plotSumTruncatedDFG]()

    wavenumbers = gold[:fullwn]
    signal = gold[:dfgTruncatedSum]
    filename = "processed_spectrum.csv"
    savepath = joinpath(path, filename)

    writedlm(savepath, [wavenumbers signal])

    return ([wavenumbers signal], datetimes, gold)

end


export NUSFSpectrum, batchprocess, processgold, calibrate


end # module
