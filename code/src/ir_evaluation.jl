######################################################
###          MODULE FOR IR-DATA-EVALUATION         ###
###            author: Michi Haugeneder            ###
######################################################
module irev

using ReadWriteDlm2, MAT, DelimitedFiles, Dates, PyCall,
    ImageFiltering, ProgressMeter, LaTeXStrings, NCDatasets
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
#not on HYPERION
if (gethostname() == "LINUX24" || gethostname() == "Michi-T450s")
    using PyPlot
end
include("general.jl")
include("fourier.jl")
import .gen
import .ft

export loadfromNetCDF4, loadexcerptfromNetCDF4, loadcolsfromNetCDF4, findfiles, saveirasnetcdf,
    fromIRloadwinddata, fromdirloadsupptime, fromdirloadMATfile, fromdirloadmultMATfiles,
    fromdirloadMATwindfile, fromdirloadMATwindandmaxqual, loadmultwindandmaxqualfiles,
    exportseqtoMAT, exportwindandmaxqualtoMAT, exportsupplementary, filesinseq,
    subtracttimeaverage, subtractbacklookingtimeaverage, setsnowsurfacetonan,
    spatialgaussfilter, completekagapreprocessing, showheatmap, showheatmaptrad, removesnowsurface,
    averagewindfield, fetchdependfunc

######################################################
###                LOADING AND I/O                 ###
######################################################
"Load variable from NetCDF4 file"
function loadfromNetCDF4(file::String, varname::String)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = dat[:, :, :]
    close(netfile)
    return dattot
end

"Load excerpt of variable from NetCDF4 file"
function loadexcerptfromNetCDF4(file::String, varname::String,
    rowrange, colrange)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = dat[rowrange, colrange, :]
    close(netfile)
    return dattot
end

"Load excerpt of variable from NetCDF4 file"
function loadcolsfromNetCDF4(file::String, varname::String,
    colrange)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = dat[:, colrange, :]
    close(netfile)
    return dattot
end

"Return vector containing files in 'folder' that match
Regular Expression 'pat'"
function findfiles(folder::String, pat::Regex)::Vector{String}
    files = readdir(folder)
    return files[occursin.(pat, files)]
end

"Saving 3D IRdata-array as NetCDF4"
function saveirasnetcdf(data::Array, name::String, target::String, deflatelvl::Int64=5)
    @info("Saving IRdata output to NetCDF4-file")
    ds = NCDataset(target, "c")
    defDim(ds, "row", size(data, 1))
    defDim(ds, "col", size(data, 2))
    defDim(ds, "frame", size(data, 3))
    defVar(ds, name, data, ("row", "col", "frame"); shuffle=true, deflatelevel=deflatelvl)
    close(ds)
end

"Load and format ultrasonic-data from Logger and return DateTime-vector,
data array, and a columnname vector"
function fromIRloadwinddata(
    filename::String, startdate::DateTime, enddate::DateTime)
    @info("loading wind data from IR-ultrasonic...")
    df = readdlm(
        string(filename), ',', String, '\n';
        skipstart=4, use_mmap=true)
    #extract first column to be DateTime
    dateformat = DateFormat("yyyy-mm-dd HH:MM:SS")
    timeofmeasure = DateTime.(df[:, 1], dateformat)
    #convert to Float64
    df = parse.(Float64, df[startdate.<timeofmeasure.<enddate, 5:6])
    #specify a String vector with the column names
    columnnames = ["Winddirection", "Windspeed"]
    return timeofmeasure[startdate.<timeofmeasure.<enddate], df, columnnames
end

"Read starttime and endtime from idx.txt file"
function fromdirloadsupptime(pathtofile::String)
    file = open(string(pathtofile, "data\\idx.txt"), "r")
    temp = read(file, String)
    dateformat = DateFormat("yyy-mm-dd HH:MM:SS")
    starttime = DateTime.(temp[10:28], dateformat)
    endtime = DateTime.(temp[31:49], dateformat)
    close(file)
    return starttime, endtime
end

"Load the data from the .MAT file @ 'locationMATfile'"
function fromdirloadMATfile(filename::String)::Array
    @info("Reading data from .mat file...", filename)
    file = matopen(string(filename), "r")
    m = read(file, "IRdata")
    close(file)
    return m
end

"Load multiple single mat files (IRdata) and concentate them"
function fromdirloadmultMATfiles(source::String, rowrange::Vector, colrange::Vector, re=r"^irdata\d{2}.mat")::Array
    files = readdir(source)
    files = files[occursin.(re, files)]
    @info("Reading in multiple .mat files. May take a while...")
    outdata = irev.fromdirloadMATfile(joinpath(source, files[1]))[rowrange, colrange, :]
    for ifile in 2:length(files)
        tmp = irev.fromdirloadMATfile(joinpath(source, files[ifile]))[rowrange, colrange, :]
        outdata = cat(outdata, tmp, dims=3)
    end
    return outdata
end

"Load the winddata from the .MAT file @ 'locationMATfile'"
function fromdirloadMATwindfile(pathtofile::String, filename::String)::Array
    @info("Reading data from .mat file...")
    file = matopen(string(pathtofile, filename), "r")
    m = read(file, "windfield")
    close(file)
    return m
end

"Load the winddata and maxqual-factor from the .MAT file @ 'locationMATfile'"
function fromdirloadMATwindandmaxqual(pathtofile::String, filename::String)
    file = matopen(string(pathtofile, filename), "r")
    #@info("Reading winddata and maxqual from .mat file...", filename)
    wind = read(file, "windfield")
    maxqual = read(file, "maxqual")
    close(file)
    return wind, maxqual
end

"Loading all windfieldfiles in 'source' directory. Optionally specify
the RegEx-Expression 're' to match the windfield-files."
function loadmultwindandmaxqualfiles(source::String, re=r"^windfield_irdata\d{2}.mat")
    files = readdir(source)
    files = files[occursin.(re, files)]
    @info("Loading multiple windfield files and concentating...", files)
    (windout, maxqualout) = irev.fromdirloadMATwindandmaxqual(source, files[1])
    @showprogress for ifiles in 2:length(files)
        (windtmp, maxqualtmp) = irev.fromdirloadMATwindandmaxqual(source, files[ifiles])
        windout = cat(windout, windtmp, dims=3)
        maxqualout = cat(maxqualout, maxqualtmp, dims=3)
    end
    return windout, maxqualout
end

"Export the obtained 3-dim array (x,y: picture; z:sequence) to .MAT file"
function exportseqtoMAT(
    pathtofolder::String, a::Array, prefix::String, nrseq::Integer,
    idxstartpic::Integer, idxendpic::Integer)
    @info("Exporting sequence to .mat file...")
    filename =
        string(prefix, string(nrseq)) #, "_", string(idxstartpic), "to", string(idxendpic))
    file = matopen(string(pathtofolder, filename, ".mat"), "w")
    write(file, "IRdata", a)
    close(file)
end

"Export the windfield and the maximum quality factor (maxqual) to .MAT file"
function exportwindandmaxqualtoMAT(
    pathtofolder::String, filename::String, wind::Array, maxqual::Array)
    @info("Exporting windfield and maxqual to .mat file...")
    file = matopen(string(joinpath(pathtofolder, filename), ".mat"), "w")
    write(file, "windfield", wind)
    write(file, "maxqual", maxqual)
    close(file)
    println("running MATLAB script to compress exported .mat-file...")
    if gethostname() == "LINUX24"
        command = string("-nodisplay -r \"cd('/home/haugened/Documents/MATLAB/'); resave_file('", joinpath(pathtofolder, filename), "'); quit\"")
        matcommand = `matlab $command`
        run(matcommand)
    end

    println("Done")
end

"Export .dat file with supplementary information"
function exportsupplementary(
    pathtofolder::String, nrseq::Integer, idxstart::Integer, idxend::Integer, timestart,
    timeend, timestep::AbstractFloat)
    @info("Exporting supplementary file...")
    sup = ["          " "INDEX" "TIME"
        "START       " string(idxstart) string(timestart)
        "END       " string(idxend) string(timeend)
        "DURATION  " string(idxend - idxstart + 1) string(Second(timeend - timestart))
        "" "" ""
        "FRAMES/S  " string(timestep) ""]
    filename =
        string(string(nrseq), "_supplementary.dat")
    #string(string(nrseq), "_", string(idxstart), "to", string(idxend),"_supplementary.dat")
    open(string(pathtofolder, filename), "w") do io
        writedlm(io, sup)
    end
end

"Returning a vector containing all irdataxx.mat files
belonging to the sequence (stored in the seq folder)."
function filesinseq(folder::String)::Vector{String}
    files = readdir(folder)
    re = r"^irdata.*\.mat"
    return files[occursin.(re, files)]
end
######################################################
###                   PREPROCESSING                ###
######################################################
"Perform time-averaging for all points as due to 'Inagaki et al. (2013)'
'nrelements' refers to third dimension of 'IRdata' (time-dimension)"
function subtracttimeaverage(IRdatain::Array, nrelements::Integer)::Array
    IRdatareturn = similar(IRdatain)
    @info("Performing time averaging...")
    @sync for i in 1:size(IRdatain, 2)
        Threads.@spawn for j in 1:size(IRdatain, 1)
            IRdatareturn[j, i, :] = IRdatain[j, i, :] - gen.movingaverage(IRdatain[j, i, :], nrelements)
        end
    end
    return IRdatareturn
end

"Perform back-looking time-averaging for all points as due to 'Inagaki et al. (2013)'
'nrelements' refers to third dimension of 'IRdata' (time-dimension)"
function subtractbacklookingtimeaverage(IRdatain::Array, nrelements::Integer)::Array
    IRdatareturn = similar(IRdatain)
    @info("Performing back-looking time averaging...")
    @sync for i in 1:size(IRdatain, 2)
        Threads.@spawn for j in 1:size(IRdatain, 1)
            IRdatareturn[j, i, :] = IRdatain[j, i, :] - gen.backlookingmovavg(IRdatain[j, i, :], nrelements)
        end
    end
    println("Done")
    return IRdatareturn
end

"Set the pixels showing the snow surface (see condition) to NaN
and return row nr. for each column"
function setsnowsurfacetonan(irdata::Array, thresh::AbstractFloat=1.5; framewise=false)
    @info("Setting snow-surface to NaN")
    irout = copy(irdata)
    if !framewise
        snowsurf = zeros(Int64, size(irdata, 2))
    else
        snowsurf = zeros(Int64, size(irdata, 2), size(irdata, 3))
    end
    #only do it for the lowest 30% assuming top is row=1
    lim = round(Int, size(irdata, 1) * 0.7)
    if !framewise
        @sync for icol in 1:size(irdata, 2)
            Threads.@spawn for jrow in lim:size(irdata, 1)
                if any(x -> x < thresh, irdata[jrow, icol, :])
                    irout[jrow, icol, :] .= NaN
                end
            end
        end
        for jcol in 1:size(irdata, 2)
            for irow in lim:size(irdata, 1)
                if isnan(irdata[irow, jcol, round(Int, (1 + size(irdata, 3)) / 2)])
                    snowsurf[jcol] = irow
                    break
                end
            end
        end
    else #do it frame-by-frame
        irtmp = irdata[lim:end, :, :]
        irtmp[irtmp.<thresh] .= NaN
        irout[lim:end, :, :] = irtmp
        for fr in 1:size(irout, 3)
            for col in 1:size(irout, 2)
                tmp = maximum(filter(x->!isnothing(x), [0 findfirst(x -> !isnan(x), irout[end:-1:1, col, fr])]))
                snowsurf[col, fr] = size(irout, 1) - tmp + 1
            end
        end
    end
    return irout, snowsurf
end

"Apply picture-wise spatial filtering using Gauss-kernel with sigmarow and sigmacol"
function spatialgaussfilter(irdata::Array, sigmrow::Integer, sigmcol::Integer)
    :Array
    @info("Picture-wise spatial Gauss-filtering with sigma_row=", sigmrow,
        " and sigma_col=", sigmcol)
    out = similar(irdata)
    for i in 1:size(irdata, 3)
        out[:, :, i] = imfilter(irdata[:, :, i], Kernel.gaussian((sigmrow, sigmcol)), NA())
    end
    return out
end

"Apply complete preprocessing to (already cut) 'irin' array."
function completekagapreprocessing(irin::Array, gaussparam::Vector{Integer},
    framespersec, backlookingavgtime, fco)
    @info("Preprocessing IR-data for Kaga-windfield-estimation")
    #setting the snow to NaN and spatial Gauss filtering
    (irtmp, snowsurf) = irev.setsnowsurfacetonan(irin)
    irtmp = irev.spatialgaussfilter(irtmp, gaussparam[1], gaussparam[2])

    #high-pass filtering = subtract backlooking time average
    irtmp2 = irev.subtractbacklookingtimeaverage(irtmp, round(Int, backlookingavgtime * framespersec))

    #low-pass filtering using Fourier transformation
    (lpir,) = ft.fouriercutoffmatrix(irtmp2, round(Int, framespersec * size(irtmp2, 3)), fco)
    return lpir, snowsurf
end

######################################################
###                 PROCESSING IR                  ###
######################################################
"Show a heatmap of a 2D-array"
function showheatmap(data::Array, min, max, title::String="IRdata", cmperpxlu=0.606, cmperpxlw=0.549, xis0pxl=30, yis0pxl=360)
    fntsze = 18
    lblsze = fntsze - 5
    patches = pyimport("matplotlib.patches")
    PyPlot.pygui(true)
    fig = PyPlot.figure(figsize=(0.469366 * 16, 5.2))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title)
    xleft = -(xis0pxl - 1) * cmperpxlu / 100
    xright = size(data, 2) * cmperpxlu / 100 + xleft
    ybottom = (yis0pxl - size(data, 1)) * cmperpxlw / 100
    ytop = size(data, 1) * cmperpxlw / 100 + ybottom
    #rect = patches.Rectangle((xleft+302*cmperpxlu/100, (size(data,1)-289)*cmperpxlw/100), 0.103, 0.0933, edgecolor="black", facecolor="None")
    #rect.set_lw(2)
    #ax.add_patch(rect)
    hm = ax.imshow(data, extent=[xleft, xright, ybottom, ytop], cmap="turbo", vmin=min, vmax=max)
    ax.set_xlabel(L"x~\mathrm{[m]}", fontsize=fntsze)
    ax.set_ylabel(L"h~\mathrm{[m]}", fontsize=fntsze)
    ax.tick_params(labelsize=lblsze)
    cbar = fig.colorbar(hm, ax=ax)
    cbar.set_label(L"T~\mathrm{[^\circ C]}", fontsize=fntsze)
    cbar.ax.tick_params(labelsize=lblsze)
    PyPlot.show()
end

"Show a heatmap of a 2D-array traditional without meter stuff"
function showheatmaptrad(data::Array, min, max, title::String="IRdata")
    fntsze = 22
    lblsze = fntsze - 5
    patches = pyimport("matplotlib.patches")
    PyPlot.pygui(true)
    fig = PyPlot.figure(figsize=(16, 5.2))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title)
    cmap = PyPlot.get_cmap("tab10")
    #rect1 = patches.Rectangle((227+90,98+425), 17, 17, edgecolor=cmap(0), facecolor="None")
    #rect2 = patches.Rectangle((260,190), 10, 10, edgecolor=cmap(1), facecolor="None")
    #rect3 = patches.Rectangle((12,325), 10, 10, edgecolor=cmap(2), facecolor="None")
    #rect4 = patches.Rectangle((800,190), 10, 10, edgecolor=cmap(3), facecolor="None")
    #rect5 = patches.Rectangle((800,300), 10, 10, edgecolor=cmap(4), facecolor="None")
    #rect6 = patches.Rectangle((900,335), 10, 10, edgecolor=cmap(5), facecolor="None")
    #rect1.set_lw(1.5)
    #rect2.set_lw(1.5)
    #rect3.set_lw(1.5)
    #rect4.set_lw(1.5)
    #rect5.set_lw(1.5)
    #rect6.set_lw(1.5)
    #ax.add_patch(rect1)
    #ax.add_patch(rect2)
    #ax.add_patch(rect3)
    #ax.add_patch(rect4)
    #ax.add_patch(rect5)
    #ax.add_patch(rect6)
    #xleft = -(xis0pxl-1)*cmperpxl/100
    #xright = size(data, 2)*cmperpxl/100+xleft
    hm = ax.imshow(data, cmap="turbo", vmin=min, vmax=max)
    PyPlot.axis("off")
    #ax.tick_params(labelsize=lblsze)
    divider = mpl_axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cbar = fig.colorbar(hm, cax=cax)
    cbar.set_label(L"T~\mathrm{[^\circ C]}", fontsize=fntsze)
    cbar.ax.tick_params(labelsize=lblsze)
    PyPlot.show()
end

"Remove snow surface (set to NaN before)"
function removesnowsurface(irin::Array)::Array
    #only do it for the lowest 30% assuming top is row=1
    lim = round(Int, size(irin, 1) * 0.7)
    idxsnowsurf = zeros(Int64, size(irin, 2))
    picpos = round(Int, size(irin, 3) / 2)
    for icol in 1:size(irin, 2)
        idxsnowsurf[icol] = findfirst(x -> isnan(x), irin[lim:end, icol, picpos]) + lim - 1
        if isnothing(idxsnowsurf[icol])
            idxsnowsurf[icol] = size(irin, 1)
        end
    end
    highestsnowsurf = minimum(idxsnowsurf)
    irout = zeros(Float64, highestsnowsurf - 1, size(irin, 2), size(irin, 3))
    for icol in 1:size(irin, 2)
        for jpic in 1:size(irin, 3)
            irout[:, icol, jpic] = irin[(-highestsnowsurf+1:-1).+idxsnowsurf[icol], icol, jpic]
        end
    end
    return irout
end

######################################################
###                 POST-PROCESSING                ###
######################################################
"Gliding average for the windfield, 'nrelements' refers to 3rd dim of IRdata"
function averagewindfield(windfield::Array, nrelements::Integer)::Array
    @info("Gliding average for windfield")
    avgwind = similar(windfield)
    for icol in 1:size(windfield, 2)
        for jcol in 1:size(windfield, 1)
            avgwind[jcol, icol, :, 1] = gen.movingaverage(windfield[jcol, icol, :, 1], nrelements)
            avgwind[jcol, icol, :, 2] = gen.movingaverage(windfield[jcol, icol, :, 2], nrelements)
        end
    end
    return avgwind
end

######################################################
###                   EVALUATION                   ###
######################################################
"Create a value table for a fetch-dependend function 'f'. 'fetchis0' is the 0 point;
fetch increases to the right. Return array with 1st col x, 2nd col y."
function fetchdependfunc(f::Function, irsize, snowsurf::Vector, fetchis0::Integer, cmperpxl::AbstractFloat)::Array
    xleftmeter = (1 - fetchis0) * cmperpxl * 0.01
    xrightmeter = (irsize[2] - fetchis0) * cmperpxl * 0.01
    xinmeter = collect(xleftmeter:cmperpxl*0.01:xrightmeter)
    for i in 1:length(xinmeter)
        if xinmeter[i] < 0.0
            xinmeter[i] = NaN
        end
    end
    res = f.(xinmeter)
    for i in 1:length(res)
        if snowsurf[i] > 270
            res[i] = (irsize[1] - snowsurf[i] - 1) * cmperpxl / 100 + res[i]
        else
            res[i] = NaN
        end
    end
    return hcat(xinmeter, res)
end

end #module