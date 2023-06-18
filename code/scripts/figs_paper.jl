######################################################
###             FIGURES FOR PAPER 2023             ###
###            author: Michi Haugeneder            ###
######################################################
######################################################
###     FIG 1: CLIMATOLOGICAL FLUX FOOTPRINTS     ###
######################################################
#rest of Fig 1 is map. Check GIMP-file

using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter, Distributed
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
LogNorm = pyimport("matplotlib.colors")
mpimg = pyimport("matplotlib.image")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
    datapath = "/home/michi/Documents/slf/data/"
    tjkpath = "/home/michi/Documents/slf/data/tjk/tjk_data.csv"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
    datapath = "/home/haugened/Documents/data/"
    tjkpath = "/home/haugened/Documents/data/tjk/tjk_data.csv"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "kljun_ffp.jl"))
import .turb
import .gen
import .kljun
PyPlot.pygui(true)
@pyinclude(joinpath(importdir, "src", "kljun_ffp_climatology.py"))

#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
kaijo_period_file = joinpath(datapath, "kaijo", "21_kaijo_periods.txt")
kaijo_outfile_stam = joinpath(datapath, "kaijo", "kaijo")
outfile_stam = joinpath(datapath, "WindSonic_Duerrboden_2005xx/tjk_sonic_200504-200514")
tower_outfile_stam = joinpath(datapath, "tower", "preproc")

#select data and measurement period to be evaluated
evalstart = DateTime(2021, 05, 20, 23, 00, 00)
evalend   = DateTime(2021, 06, 10, 23, 00, 00)
#evalend = evalstart + Day(10)

#evaldf1 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t1irgdb.nc"), evalstart, evalend)
#evaldf2 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2irgdb.nc"), evalstart, evalend)
evaldf3 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2lcsatdb.nc"), evalstart, evalend)
#evaldf4 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2ucsatdb.nc"), evalstart, evalend)
#evaldf5 = turb.readturbasnetcdf(string(kaijo_outfile_stam, ".nc"), evalstart, evalend)
evaldf6 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "tjkdf.nc"), evalstart, evalend)

#apply NaN-mask to T1 & T2 when repositioned (for DR)
#turb.repositionnanmask!(evaldf1)
#turb.repositionnanmask!(evaldf2)
turb.repositionnanmask!(evaldf3)
#turb.repositionnanmask!(evaldf4)

#show statistics about missing data
#turb.printmissstats(evaldf1)
#turb.printmissstats(evaldf2)
turb.printmissstats(evaldf3)
#turb.printmissstats(evaldf4)
#turb.printmissstats(evaldf5)
turb.printmissstats(evaldf6)
#interpolate missing
#evaldf1 = turb.interpolatemissing(evaldf1)
#evaldf2 = turb.interpolatemissing(evaldf2)
evaldf3 = turb.interpolatemissing(evaldf3)
#evaldf4 = turb.interpolatemissing(evaldf4)
#=
try
    evaldf5 = turb.interpolatemissing(evaldf5)
catch LoadError
    @warn("Kaijo DataFrame empty")
end
=#
evaldf6 = turb.interpolatemissing(evaldf6)
#double rotation
#turb.drdf!(evaldf1)
#turb.drdf!(evaldf2)
turb.drdf!(evaldf3)
#turb.drdf!(evaldf4)
#turb.drdf!(evaldf5)
turb.drdf!(evaldf6, periodwise=false)

tjkmeteodata = turb.csvtodataframe(tjkpath)
tjkmeteodata = tjkmeteodata[evalstart.<=tjkmeteodata.time.<=evalend, :]

timestep = Millisecond(50)
ρ_air = 1.2 #kg m^{-3}
c_p = 1004 #J kg^{-1} K^{-1}
L_v = 2450e3 #J kg^{-1} (approx @0°C) 

#turb.missing2nan!(evaldf1)
#turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
#turb.missing2nan!(evaldf4)
#turb.missing2nan!(evaldf5)
turb.missing2nan!(evaldf6)
#=
#apply up/down valley wind condition
disallowmissing!(tjkmeteodata)
disallowmissing!(evaldf6)
winddir_fine = zeros(Float64, size(evaldf6, 1))
turb.extendtofinertimeseries!(winddir_fine, evaldf6.time, tjkmeteodata.wind_mean_vector_direction, tjkmeteodata.time)
windclass = fill(0, length(winddir_fine))
windclass[310 .<= winddir_fine .|| winddir_fine .<=20] .= 1
windclass[130 .<= winddir_fine .<= 200] .= 2
evaldf3 = evaldf3[windclass .== 1, :]
evaldf6 = evaldf6[windclass .== 1, :]
=#
#Reynolds averaging times
#ra1 = Second(100) #Second(330)
#ra2 = Second(100) #Second(330)
ra3 = Second(100) #Second(320)
#ra4 = Second(100) #Second(210)
#ra5 = Second(100) #Second(300)
ra6 = Second(50) #Second(55)

#fx1_raw = turb.turbflux(evaldf1, ra1)
#fx2_raw = turb.turbflux(evaldf2, ra2)
fx3_raw = turb.turbflux(evaldf3, ra3)
#fx4_raw = turb.turbflux(evaldf4, ra4)
#fx5_raw = turb.turbflux(evaldf5, ra5)
fx6_raw = turb.turbflux(evaldf6, ra6)

#averaging
#fx1 = turb.avgflux(fx1_raw, Second(600))
#fx2 = turb.avgflux(fx2_raw, Second(600))
fx3 = turb.avgflux(fx3_raw, Second(600))
#fx4 = turb.avgflux(fx4_raw, Second(600))
#fx5 = turb.avgflux(fx5_raw, Second(600))
fx6 = turb.avgflux(fx6_raw, Second(600))

#Obukhov-length
#L1 = turb.obukhov(fx1_raw, Minute(30))
#L2 = turb.obukhov(fx2_raw, Minute(10))
L3 = turb.obukhov(fx3_raw, Minute(10))
#L4 = turb.obukhov(fx4_raw, Minute(10))
#try
#    L5 = turb.obukhov(fx5_raw, Minute(30))
#catch e
#end
L6 = turb.obukhov(fx6_raw, Minute(10))

names = [:evaldf1, :evaldf2, :evaldf3, :evaldf4, :evaldf5, :evaldf6]
meas_heights = [1.2, 0.9, 1.9, 2.8, 0.3, 5]
pbl_height = 1000.0
Ls = [:L1, :L2, :L3, :L4, :L5, :L6]
fluxes = [:fx1, :fx2, :fx3, :fx4, :fx5, :fx6]
outnames = [:ffp1, :ffp2, :ffp3, :ffp4, :ffp5, :ffp6]
aggtime = Minute(30) #aggregation time

#optional input
domain = nothing
dx = nothing
dy = nothing
nx = nothing
ny = nothing
rs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] #levels for plotting
rslayer = 0 #measurement within roughness sublayer (theory not working properly)
smooth_data = 1
crop = false #crop output to maximum defined rs (max 0.9)
pulse = nothing
verbosity = 2
fig = false

output = nothing

for ix in [3, 6]#1:size(names, 1)
    println("Calculating footprint for ", String(names[ix]))

    ecdata = @eval $(names[ix])
    turb.missing2nan!(ecdata)
    fluxdata = @eval $(fluxes[ix])
    obukl = @eval $(Ls[ix])

    #aggregate data to half-hour intervals
    duration = ecdata.time[end]-ecdata.time[1]
    aggidcs = round(Int, aggtime/Millisecond(50))
    nrblocks = div(duration, aggtime)
    #nrblocks = div(size(ecdata, 1),aggidcs)
    println("Evaluating ", nrblocks, " blocks, last ", mod(duration, aggtime), " discarded.")
    umean = fill(NaN, nrblocks)
    h = fill(pbl_height, nrblocks)
    ol = fill(NaN, nrblocks)
    sigmav = fill(NaN, nrblocks)
    ustar = fill(NaN, nrblocks)
    wind_dir = fill(NaN, nrblocks)

    for j in 1:nrblocks
        six = 1 + (j - 1)*aggidcs
        eix = six + (aggidcs - 1)
        umean[j] = mean(filter(!isnan, sqrt.(ecdata.u[six:eix] .^2 .+ ecdata.v[six:eix] .^2 .+ ecdata.w[six:eix] .^2)))
        ol[j] = mean(filter(!isnan, obukl.L[six:eix]))
        sigmav[j] = std(filter(!isnan, ecdata.v[six:eix]))
        ustar[j] = mean(filter(!isnan, fluxdata.u_star[six:eix]))
        wind_dir[j] = mean(filter(!isnan, tjkmeteodata.wind_mean_vector_direction[ecdata.time[six] .<= tjkmeteodata.time .<= ecdata.time[eix]]))
    end

    output = py"FFP_climatology"(meas_heights[ix], nothing, PyVector(umean), PyVector(h), PyVector(ol),
    PyVector(sigmav), PyVector(ustar), PyVector(wind_dir), domain, dx, dy, nx, ny, PyVector(rs), rslayer, smooth_data, crop, pulse, verbosity, fig)
    if output["flag_err"] != 0
        @warn("Error flag set to true! There is an error!")
    end
    @eval $(outnames[ix]) = output
end

#plotting the footprint on the ortho-mosaic

fileorthomosaic = joinpath(datapath, "pics", "Drone_orthomosaic_20210603_10cm_excerpt2.png")
orthomosaic = mpimg.imread(fileorthomosaic)
#PyPlot.imshow(orthomosaic)
#location of flux measurements 1-6 in original image
#[row-location, col-location]

fluxloc = [1452 926; 1540 960; 1540 960; 836 1265; 1116 1387; 972 892]

#extend of background [row, col]
#bgextend_m = [280, 280] #in m from measuring in GIS: 279.9
bgextend_pxl = [size(orthomosaic, 1), size(orthomosaic, 2)] #in pxl

#calculate m/pxl from it
meterperpxl_row = 0.1 #bgextend_m[1] / bgextend_pxl[1]
meterperpxl_col = 0.1 #bgextend_m[2] / bgextend_pxl[2]

#origin of figure
figorigin = [1540 960] #tower 2

#calculate fluxloc in new coordinates [m]
fluxloc_final = Array{Float64}(undef, size(fluxloc, 1), size(fluxloc, 2))
fluxloc_final[:, 1] = (figorigin[1] .- fluxloc[:, 1]) .* meterperpxl_row
fluxloc_final[:, 2] = (fluxloc[:, 2] .- figorigin[2])  .* meterperpxl_col

#calculate extend in new coordinates
lft = (-figorigin[2]) * meterperpxl_col
rght = (bgextend_pxl[2]-1-figorigin[2]) *meterperpxl_col
tp = (bgextend_pxl[1]-(figorigin[1]- bgextend_pxl[1])-1) * meterperpxl_row
btm = (figorigin[1] - bgextend_pxl[1]) * meterperpxl_row

bgextend_final = (-figorigin[2], bgextend_pxl[2]-1-figorigin[2], -(bgextend_pxl[1]-figorigin[1]), bgextend_pxl[1]-(bgextend_pxl[1]-figorigin[1])-1).*meterperpxl_col
#(lft, rght, btm, tp)
##
ctab10 = PyPlot.cm.tab10
ffp_fig = PyPlot.figure(figsize=(5.5,7))
ax1 = ffp_fig.add_subplot(111)
#ax1.set_title("Flux footprints 70% - Kljun et al. (2015)")
bg = ax1.imshow(orthomosaic, extent=bgextend_final)
#bg = ax1.pcolormesh(orthomosaic)
ax1.set_xlabel("x [m]")
ax1.set_ylabel("y [m]")
locfx1 = ax1.plot(fluxloc_final[1, 2], fluxloc_final[1, 1], ".", color=ctab10(4), markersize=15)#, label="T1IRG")
#locfx2 = ax1.plot(fluxloc_final[2, 2], fluxloc_final[2, 1], ".", color=ctab10(1))#, label="T2IRG")
locfx3 = ax1.plot(fluxloc_final[3, 2], fluxloc_final[3, 1], ".", color=ctab10(1), markersize=15)#, label="T2LCSAT")
#locfx4 = ax1.plot(fluxloc_final[4, 2], fluxloc_final[4, 1], ".", color=ctab10(3))#, label="T2UCSAT")
#locfx5 = ax1.plot(fluxloc_final[5, 2], fluxloc_final[5, 1], ".", color=ctab10(4))#, label="Kaijo")
locfx6 = ax1.plot(fluxloc_final[6, 2], fluxloc_final[6, 1], ".", color=ctab10(2), markersize=15)#, label="TJK")
#fp1 = ax1.plot(ffp1["xr"][:, end-1] .+ fluxloc_final[1, 2], ffp1["yr"][:, end-1] .+ fluxloc_final[1, 1], color=ctab10(0), label = "T1IRG")
#fp2 = ax1.plot(ffp2["xr"][:, end-1] .+ fluxloc_final[2, 2], ffp2["yr"][:, end-1] .+ fluxloc_final[2, 1], color=ctab10(1), label = "T2IRG")
fp3 = ax1.plot(ffp3["xr"][end-2] .+ fluxloc_final[3, 2], ffp3["yr"][end-2] .+ fluxloc_final[3, 1], color=ctab10(1), label = "T2LCSAT", lw=3)
#fp4 = ax1.plot(ffp4["xr"][:, end-1] .+ fluxloc_final[4, 2], ffp4["yr"][:, end-1] .+ fluxloc_final[4, 1], color=ctab10(3), label = "T2UCSAT")
#fp5 = ax1.plot(ffp5["xr"][:, end-1] .+ fluxloc_final[5, 2], ffp5["yr"][:, end-1] .+ fluxloc_final[5, 1], color=ctab10(4), label = "Kaijo")
fp6 = ax1.plot(ffp6["xr"][end-2] .+ fluxloc_final[6, 2], ffp6["yr"][end-2] .+ fluxloc_final[6, 1], color=ctab10(2), label = "TJK", lw=3)
#ax1.legend()
##

######################################################
###              FIG 2: CMP MRDs                   ###
######################################################
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter, Distributed
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
    datapath = "/home/michi/Documents/slf/data/"
    tjkpath = "/home/michi/Documents/slf/data/tjk/tjk_data.csv"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
    datapath = "/home/haugened/Documents/data/"
    tjkpath = "/home/haugened/Documents/data/tjk/tjk_data.csv"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
import .turb
import .gen
import .MRD
PyPlot.pygui(true)

#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
kaijo_period_file = joinpath(datapath, "kaijo", "21_kaijo_periods.txt")
kaijo_outfile_stam = joinpath(datapath, "kaijo", "kaijo")
outfile_stam = joinpath(datapath, "WindSonic_Duerrboden_2005xx/tjk_sonic_200504-200514")
tower_outfile_stam = joinpath(datapath, "tower", "preproc")

#select data and measurement period to be evaluated
evalstart = DateTime(2021, 05, 29, 08, 00, 00)
evalend   = DateTime(2021, 05, 29, 20, 00, 00)
#evalend = evalstart + Day(10)

evaldf3 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2lcsatdb.nc"), evalstart, evalend)

#apply NaN-mask to T1 & T2 when repositioned (for DR)
turb.repositionnanmask!(evaldf3)

#show statistics about missing data
turb.printmissstats(evaldf3)
#interpolate missing
evaldf3 = turb.interpolatemissing(evaldf3)
#double rotation
turb.drdf!(evaldf3)

evaldf = evaldf3
turb.missing2nan!(evaldf)

#Traditional orthonogal Mulitresolution Flux Decomposition
##
M = 15
blocklength = ceil(2^M * timestep, Dates.Second)
(mrd_x, mrd_data, time_middle) = MRD.completemrd(evaldf, "w", "T", M, round(Int, 0.1 * 2^M), normed=false)
(mrd_time, mrd_D_median, mrd_D_min, mrd_D_max, mrd_D_quant1, mrd_D_quant3) =
    MRD.mrdpp(mrd_x, mrd_data)

#non-orthogonal MRD
nrmodes1 = 100
nrmodes2 = M

(nomrd_x_1, nomrd_data_1, nomrd_time_middle_1) = MRD.completenomrd(evaldf, "w", "T", M, round(Int, 0.1 * 2^M), nrmodes1, normed=true)
(nomrd_time_1, nomrd_D_median_1, nomrd_D_min_1, nomrd_D_max_1, nomrd_D_quant1_1, nomrd_D_quant3_1) =
    MRD.mrdpp(nomrd_x_1[:, 2], nomrd_data_1)
@info(length(nomrd_D_median_1))
(nomrd_x_2, nomrd_data_2, nomrd_time_middle_2) = MRD.completenomrd(evaldf, "w", "T", M, round(Int, 0.1 * 2^M), nrmodes2, normed=true)
(nomrd_time_2, nomrd_D_median_2, nomrd_D_min_2, nomrd_D_max_2, nomrd_D_quant1_2, nomrd_D_quant3_2) =
    MRD.mrdpp(nomrd_x_2[:, 2], nomrd_data_2)
##
cmap = PyPlot.get_cmap("tab10")
gs = gridspec.GridSpec(1,2)
fig = PyPlot.figure(figsize=(9.5,4.6))
ax = fig.add_subplot(gs[1,1])
ax.set_title("a)", fontfamily="serif", loc="left", fontsize="large")
ax.set_xlabel(L"$\tau~\mathrm{[s]}$")
ax.set_ylabel(L"$C_{wT_S}~\mathrm{[10^{-3}~K~m~s^{-1}]}$")
ax.grid(true)
PyPlot.xscale("log")
ax.plot(Dates.value.(mrd_time) ./ 1000, mrd_D_median .*1e3, ".-", label="traditional")
ax.fill_between(Dates.value.(mrd_time) ./ 1000, mrd_D_quant1 .* 1e3, mrd_D_quant3 .*1e3, alpha=0.3)
ax.plot(Dates.value.(nomrd_time_2) ./ 1000, nomrd_D_median_2 .* 1e3, ".--", label=string("continuous n = ", length(nomrd_D_median_2)))
ax.fill_between(Dates.value.(nomrd_time_2) ./ 1000, nomrd_D_quant1_2 .* 1e3, nomrd_D_quant3_2 .* 1e3, alpha=0.3)
ax.set_ylim(-7, 1)
#ax.legend()
ax2 = fig.add_subplot(gs[1,2], sharey=ax)
ax2.set_title("b)", fontfamily="serif", loc="left", fontsize="large")
ax2.set_xlabel(L"$\tau~\mathrm{[s]}$")
ax2.grid(true)
PyPlot.xscale("log")
ax2.plot(Dates.value.(nomrd_time_2) ./ 1000, nomrd_D_median_2 .* 1e3, ".--", label=string("continuous n = ", length(nomrd_D_median_2)), color=cmap(1))
ax2.fill_between(Dates.value.(nomrd_time_2) ./ 1000, nomrd_D_quant1_2 .* 1e3, nomrd_D_quant3_2 .* 1e3, alpha=0.4, color=cmap(1))
ax2.plot(Dates.value.(nomrd_time_1) ./ 1000, nomrd_D_median_1 .* 1e3, label=string("continuous n = ", length(nomrd_D_median_1)), color=cmap(2), alpha=0.7)
ax2.fill_between(Dates.value.(nomrd_time_1) ./ 1000, nomrd_D_quant1_1 .* 1e3, nomrd_D_quant3_1 .* 1e3, alpha=0.3, color=cmap(2))
##

######################################################
###              FIG 3: METEO-DATA                ###
######################################################

using Dates, PyCall, DataFrames, Statistics, LaTeXStrings
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
    datapath = "/home/michi/Documents/slf/data/"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
    datapath = "/home/haugened/Documents/data/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .gen
PyPlot.pygui(true)

tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))

#select data and measurement period to be evaluated
evaldf = tjkmeteo
#evaldf = t1irg

evalstart = DateTime(2021, 05, 20, 23, 00, 00)
evalend = DateTime(2021, 06, 10, 23, 00, 00)
#evalend = evalstart + Day(10)

#create winddirection classes
windclass = fill(NaN, length(tjkmeteo.wind_mean_vector_direction))
windclass[310 .<= tjkmeteo.wind_mean_vector_direction.||tjkmeteo.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkmeteo.wind_mean_vector_direction .<= 200] .= 7

evaldf = evaldf[evalstart.<=evaldf[:, 1].<=evalend, :]

evaldf.time .+= Hour(1) #shift from CET to CEST (LT)

SCFtime = [DateTime(2021, 05, 21, 00, 00, 00),
           DateTime(2021, 05, 26, 12, 50, 00), DateTime(2021, 05, 28, 12, 10, 00),
           DateTime(2021, 05, 30, 15, 30, 00), DateTime(2021, 06, 02, 12, 40, 00),
           DateTime(2021, 06, 03, 09, 10, 00), DateTime(2021, 06, 04, 13, 30, 00),
           DateTime(2021, 06, 07, 14, 15, 00), DateTime(2021, 06, 08, 10, 30, 00),
           DateTime(2021, 06, 09, 11, 30, 00), DateTime(2021, 06, 23, 10, 30, 00)]

SCFvalue = [1, 0.85, 0.82, 0.73, 0.47, 0.41, 0.21, 0.09, 0.08, 0.05, 0.01]
##

PyPlot.pygui(true)
nrrows = 6
cmap_1 = PyPlot.get_cmap("tab10")
majorloc = pydates.DayLocator()
minorloc = pydates.HourLocator(interval=6)
alphamajorgrid = 0.25
alphaminorgrid = 0.15
majorformatter = pydates.DateFormatter("%d.%m")
minorformatter = pydates.DateFormatter("%H")
fig = PyPlot.figure(constrained_layout=true, figsize=(11.2, 7.3))
#fig.suptitle("Meteodata Duerrboden TJK station")
gs = gridspec.GridSpec(nrrows, 1, height_ratios=[3, 1, 2, 3, 3, 2])
#last axes: snow height
axl = fig.add_subplot(gs[nrrows, 1])
axl.plot(SCFtime, SCFvalue, color=cmap_1(9))
axl.xaxis.set_major_locator(majorloc)
axl.xaxis.set_minor_locator(minorloc)
axl.xaxis.set_major_formatter(majorformatter)
axl.tick_params(axis="x", which="major", labelsize="large")
axl.xaxis.set_minor_formatter(minorformatter)
axl.tick_params(axis="x", which="minor", labelbottom=false)
axl.set_ylabel(L"f_{s}~\mathrm{[1]}")
axl.grid(b=true, which="major", color="black", alpha=alphamajorgrid)
axl.grid(b=true, which="minor", color="black", alpha=alphaminorgrid)
axl.set_xlim(DateTime(2021,05,21,00,00,00), DateTime(2021,06,11,00,00,00))
for label in axl.xaxis.get_ticklabels()[1:2:end]
    label.set_visible(false)
end
axl.set_ylim(0,1)
#tair_CS215_avg
ax1 = fig.add_subplot(gs[1, 1], sharex=axl)
ax1.plot(evaldf.time, evaldf.tair_CS215_avg, color=cmap_1(0))
ax1.set_ylabel(L"T_{air}~\mathrm{[^\circ C]}")
ax1.xaxis.set_major_locator(majorloc)
ax1.xaxis.set_minor_locator(minorloc)
ax1.tick_params(axis="x", labelbottom=false)
ax1.xaxis.set_major_formatter(majorformatter)
ax1.xaxis.set_minor_formatter(minorformatter)
ax1.tick_params(axis="x", which="minor", labelbottom=false)
ax1.grid(b=true, which="major", color="black", alpha=alphamajorgrid)
ax1.grid(b=true, which="minor", color="black", alpha=alphaminorgrid)
#rel. humidity
ax2 = fig.add_subplot(gs[2, 1], sharex=axl)
ax2.plot(evaldf.time, evaldf.humidity_rel_ClimaVUE50, color=cmap_1(1))
ax2.xaxis.set_major_locator(majorloc)
ax2.xaxis.set_minor_locator(minorloc)
ax2.xaxis.set_major_formatter(majorformatter)
ax2.xaxis.set_minor_formatter(minorformatter)
ax2.tick_params(axis="x", which="minor", labelbottom=false)
ax2.set_ylabel(L"RH~\mathrm{[\%]}")
ax2.tick_params(axis="x", labelbottom=false)
ax2.grid(b=true, which="major", color="black", alpha=alphamajorgrid)
ax2.grid(b=true, which="minor", color="black", alpha=alphaminorgrid)
#mean radiation in
ax3 = fig.add_subplot(gs[3, 1], sharex=axl)
ax3.plot(evaldf.time, evaldf.radiation_in_avg, color=cmap_1(2))
ax3.xaxis.set_major_locator(majorloc)
ax3.xaxis.set_minor_locator(minorloc)
ax3.xaxis.set_major_formatter(majorformatter)
ax3.xaxis.set_minor_formatter(minorformatter)
ax3.tick_params(axis="x", which="minor", labelbottom=false)
ax3.set_ylim(0,1000)
ax3.set_ylabel(L"P_{in}~\mathrm{[W~m^{-2}]}")
ax3.tick_params(axis="x", labelbottom=false)
ax3.grid(b=true, which="major", color="black", alpha=alphamajorgrid)
ax3.grid(b=true, which="minor", color="black", alpha=alphaminorgrid)
#wind speed
ax4 = fig.add_subplot(gs[4, 1], sharex=axl)
ax4.plot(evaldf.time, evaldf.wind_mean_scalar, color=cmap_1(3))
ax4.plot(evaldf.time, evaldf.wind_max, color=cmap_1(3), lw=1, alpha=0.6)
ax4.xaxis.set_major_locator(majorloc)
ax4.xaxis.set_minor_locator(minorloc)
ax4.xaxis.set_major_formatter(majorformatter)
ax4.xaxis.set_minor_formatter(minorformatter)
ax4.tick_params(axis="x", which="minor", labelbottom=false)
ax4.set_ylabel(L"u~\mathrm{[m~s^{-1}]}")
ax4.tick_params(axis="x", labelbottom=false)
ax4.grid(b=true, which="major", color="black", alpha=alphamajorgrid)
ax4.grid(b=true, which="minor", color="black", alpha=alphaminorgrid)
#wind direction
ax5 = fig.add_subplot(gs[5, 1], sharex=axl)
#ax5.plot(evaldf.time, evaldf.wind_std_direction .* 4, color="black", alpha=0.5, lw=1)
ax5.plot(evaldf.time, evaldf.wind_mean_vector_direction, color=cmap_1(4))
ax5.xaxis.set_major_locator(majorloc)
ax5.xaxis.set_minor_locator(minorloc)
ax5.xaxis.set_major_formatter(majorformatter)
ax5.xaxis.set_minor_formatter(minorformatter)
ax5.tick_params(axis="x", which="minor", labelbottom=false)
yl = ax5.get_ylim()
ax5.pcolormesh(pydates.date2num.(tjkmeteo.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax5.set_ylabel(L"dir~\mathrm{[^\circ]}")
ax5.tick_params(axis="x", labelbottom=false)
ax5.set_ylim(-20,380)
ax5.grid(b=true, which="major", color="black", alpha=alphamajorgrid)
ax5.grid(b=true, which="minor", color="black", alpha=alphaminorgrid)
##

######################################################
###       FIG 4: MRD FOR REYNOLDS-AVG TIME        ###
######################################################

using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

readfilename = Array{String}(undef, 2)
readfilename[2] = joinpath(datapath, "2dmrd", "db_tot", "tjk_wT_norm.nc")
readfilename[1] = joinpath(datapath, "2dmrd", "db_tot", "t2lcsat_wT_norm.nc")

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Millisecond(2^9*timestep),
    "KAIJO" => Second(400))

occurpos1 = 0.0
occurpos2 = 0.0
occurneg1 = 0.0
occurneg2 = 0.0

##
fig = PyPlot.figure(figsize=(12,4.8))
gs = gridspec.GridSpec(1,2)
for i in 1:length(readfilename)
    (time_middle, x, mrd_data, meanwind, evalstart, evalend,
        cols_for_mrd, nrpoints, blocklen, shift, sourcefile) =
        MRD.read2dmrdfromnetcdf(readfilename[i])

    figtitle = MRD.create2dmrdtitle(sourcefile, evalstart, evalend)

    fxavgperiod = blocklen * (evalend - evalstart)

    #load TJK data for wind direction info
    tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
    tjkdata = tjkmeteo[evalstart.<=tjkmeteo[:, 1].<=evalend, :]

    #create winddirection classes
    windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
    windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
    windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

    #load data and calc fluxes
    inst = MRD.instrumentnamefromsourcefile(sourcefile)
    fluxfile = inst2file[inst]
    ra = inst2ra[inst]

    #take subinterval from 2D-MRD
    timemiddle_start = DateTime(2021, 05, 01, 00, 00, 00)
    timemiddle_end = DateTime(2021, 07, 01, 00, 00, 00)

    mrd_data_sub = mrd_data[:, timemiddle_start.<=time_middle.<=timemiddle_end]
    time_middle_sub = time_middle[timemiddle_start.<=time_middle.<=timemiddle_end]
    x_sub = x[:, timemiddle_start.<=time_middle.<=timemiddle_end]

    #categorize wT-decomposition in stability classes
    mrd_data_todo = mrd_data_sub

    stabclass = zeros(Float64, size(mrd_data_todo, 2))
    sumupto = Second(50)
    sumuptoidx = findfirst(x -> x > sumupto, x[:, 1])
    for i in 1:size(mrd_data_todo, 2)
        stabclass[i] = sum(mrd_data_todo[1:sumuptoidx, i])
    end

    if i == 1
        occurpos1 = count(x -> x .>= 0, stabclass)/size(mrd_data_todo, 2)
        occurneg1 = count(x -> x .< 0, stabclass)/size(mrd_data_todo, 2)
    elseif i == 2
        occurpos2 = count(x -> x .>= 0, stabclass)/size(mrd_data_todo, 2)
        occurneg2 = count(x -> x .< 0, stabclass)/size(mrd_data_todo, 2)
    end

    posmrd_data = mrd_data_sub[:, stabclass.>=0]
    negmrd_data = mrd_data_sub[:, stabclass.<0]

    posmrd_med = zeros(Float64, size(posmrd_data, 1))
    posmrd_q1 = similar(posmrd_med)
    posmrd_q3 = similar(posmrd_med)
    negmrd_med = zeros(Float64, size(negmrd_data, 1))
    negmrd_q1 = similar(negmrd_med)
    negmrd_q3 = similar(negmrd_med)

    for i in 1:size(posmrd_data, 1)
        (posmrd_q1[i], posmrd_med[i], posmrd_q3[i]) = quantile(filter(!isnan, posmrd_data[i, :]), [0.25, 0.5, 0.75])
        (negmrd_q1[i], negmrd_med[i], negmrd_q3[i]) = quantile(filter(!isnan, negmrd_data[i, :]), [0.25, 0.5, 0.75])
    end

    title = "xxx"
    if inst == "TJK"
        title = "5 m (meteo station)"
    elseif inst == "T2LCSAT"
        title = "2 m (T2)"
    end

    if i == 1
        ax = fig.add_subplot(gs[1,1])
        ax.set_title(title)
        fig.text(0.12, 0.9, "a)", size="large", )
        PyPlot.axhline(0, color="grey", alpha=0.5)
        ax.plot(Dates.value.(x[:, 1]) ./ 1000, posmrd_med .* 1e3, color=cramericm.vik(0.9), label="unstable+neutral")
        ax.fill_between(Dates.value.(x[:, 1]) ./ 1000, posmrd_q1 .* 1e3, posmrd_q3 .* 1e3, color=cramericm.vik(0.9), alpha=0.5)
        ax.plot(Dates.value.(x[:, 1]) ./ 1000, negmrd_med .* 1e3, color=cramericm.vik(0.1), label="stable", ls="dotted")
        ax.fill_between(Dates.value.(x[:, 1]) ./ 1000, negmrd_q1 .* 1e3, negmrd_q3 .* 1e3, color=cramericm.vik(0.1), alpha=0.5)
        PyPlot.xscale("log")
        ax.grid(alpha=0.5)
        ax.set_xlabel(L"$\tau~\mathrm{[s]}$")
        #ax.legend()
        ax.set_ylabel(L"C_{wT_S} [10^{-3}~\mathrm{K~m~s^{-1}}]")
        ax.set_ylim(-5, 12)
        PyPlot.axvline(100, 0, 1, alpha=0.5, color="black", ls="--")
    elseif i==2
        ax2 = fig.add_subplot(gs[1,2])
        ax2.set_title(title)
        fig.text(0.544,0.9, "b)", size="large")
        PyPlot.axhline(0, color="grey", alpha=0.5)
        ax2.plot(Dates.value.(x[:, 1]) ./ 1000, posmrd_med .* 1e3, color=cramericm.vik(0.9), label="unstable+neutral")
        ax2.fill_between(Dates.value.(x[:, 1]) ./ 1000, posmrd_q1 .* 1e3, posmrd_q3 .* 1e3, color=cramericm.vik(0.9), alpha=0.5)
        ax2.plot(Dates.value.(x[:, 1]) ./ 1000, negmrd_med .* 1e3, color=cramericm.vik(0.1), label="stable", ls="dotted")
        ax2.fill_between(Dates.value.(x[:, 1]) ./ 1000, negmrd_q1 .* 1e3, negmrd_q3 .* 1e3, color=cramericm.vik(0.1), alpha=0.5)
        PyPlot.xscale("log")
        ax2.grid(alpha=0.5)
        ax2.set_xlabel(L"$\tau~\mathrm{[s]}$")
        #ax2.legend()
        ax2.set_ylim(-5, 12)
        PyPlot.axvline(30, 0, 1, alpha=0.5, color="black", ls="--")

    end
end

@show occurneg1
@show occurpos1
@show occurneg2
@show occurpos2
##
######################################################
###          FIG 5: WT MOVING-WINDOW MRD           ###
######################################################

using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
lines = pyimport("matplotlib.lines")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

readfilename = joinpath(datapath, "2dmrd", "db_tot", "t2lcsat_wT_norm.nc")

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Second(30),
    "KAIJO" => Second(400))

(time_middle, x, mrd_data, meanwind, evalstart, evalend,
    cols_for_mrd, nrpoints, blocklen, shift, sourcefile) =
    MRD.read2dmrdfromnetcdf(readfilename)

rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< time_middle .< DateTime(2021,06,05,19,35,00)

mrd_data[:, rain_for_2dmrd] .= NaN
    
mrd_data[:, DateTime(2021,06,05,18,35,00) .< time_middle .< DateTime(2021,06,05,19,35,00)] .= NaN

figtitle = MRD.create2dmrdtitle(sourcefile, evalstart, evalend)

fxavgperiod = blocklen * (evalend - evalstart)

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[evalstart.<=tjkmeteo[:, 1].<=evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

#load data and calc fluxes
inst = MRD.instrumentnamefromsourcefile(sourcefile)
fluxfile = inst2file[inst]
ra = inst2ra[inst]

fluxdf = turb.readturbasnetcdf(fluxfile, evalstart, evalend)
#show statistics about missing data
turb.printmissstats(fluxdf)
#interpolate missing
fluxdf = turb.interpolatemissing(fluxdf)
#double rotation
if inst == "TJK"
    turb.drdf!(fluxdf, periodwise=false)
else
    turb.drdf!(fluxdf)
end
turb.missing2nan!(fluxdf)
wndspd = gen.movingaverage(sqrt.(fluxdf.u .^ 2 .+ fluxdf.v .^ 2 .+ fluxdf.w .^ 2), 20 * 3600)
fx1_raw = turb.turbflux(fluxdf, ra)
fx1 = turb.avgflux(fx1_raw, fxavgperiod)
fx1_raw = nothing

fx1 = fx1[1:1200:end, :]

rain_for_fx =
DateTime(2021,05,22,00,08,00) .< fx1.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx1.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx1.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx1.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx1.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx1.time .< DateTime(2021,06,05,19,35,00)

fx1[rain_for_fx, 2:end] .= NaN

#-----------------------------------------------
#temporal y-scale

if "w" in cols_for_mrd
    if "u" in cols_for_mrd
        labeltext = L"$C_{uw} [10^{-3}~\mathrm{m^2~s^{-2}}]$"
        labelbottom = L"$\overline{u'w'} [\mathrm{m^2~s^{-2}}]$"
        flux = fx1.uw
    elseif "v" in cols_for_mrd
        labeltext = L"$C_{vw} [10^{-3}~\mathrm{m^2~s^{-2}}]$"
        labelbottom = L"$\overline{v'w'} [\mathrm{m^2~s^{-2}}]$"
        flux = fx1.vw
    elseif "T" in cols_for_mrd
        labeltext = L"$C_{wT_s} [10^{-3}~\mathrm{K~m~s^{-1}}]$"
        labelbottom = L"$\overline{w'T_s'} [\mathrm{K~m~s^{-1}}]$"
        flux = fx1.wT
    else
        @warn("neither 'u' nor 'v' nor 'T' in cols_for_mrd!")
    end
else
    @warn("no w in cols_for_mrd!")
end

firstbiggerra = findfirst(x -> x > ra, Millisecond.(Dates.value.(x[:, round(Int, size(x, 2)/2)])))

##
#plot temporal scale (y)
cmin = -5 #colorbar limits
cmax = -cmin
date_format = pydates.DateFormatter("%d.%m.")
majorloc = pydates.DayLocator(interval=2)
minorloc = pydates.DayLocator(interval=1)
fig = PyPlot.figure(figsize=(9, 9))#5.3))
fig.text(0.1, 0.895, "a)", size="large")
gs = gridspec.GridSpec(5, 4, width_ratios=[24, 3, 24, 1], height_ratios=[77, 615, 77, 40, 231], wspace=0.05)
ax4 = fig.add_subplot(py"$(gs)[0,0:3]")
#ax4.set_title(figtitle)
ax4.plot(fluxdf.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax4.get_ylim()
ax4.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax4.set_ylabel(L"| \mathbf{u}|~\mathrm{[m~s^{-1}]}")
ax4.tick_params(axis="x", labelbottom=false)
ax4.grid()
ax4.set_ylim(0,8)
ax = fig.add_subplot(py"$(gs)[1,0:3]", sharex=ax4)
mrdplot = ax.pcolormesh(time_middle, Dates.value.(x) ./ 10^3, mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax.set_ylabel(L"$\tau~\mathrm{[s]}$")
ax.xaxis_date()
ax.set_ylim([0.1, 2e3])
ax.axhline(Dates.value(Millisecond(ra))/1e3, ls="--", color="black", alpha=0.3)
ax.fill_betweenx((Dates.value(Millisecond(ra))/1e3, last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(py"$(gs)[1, 3]")
cbar = fig.colorbar(mrdplot, cax=ax2)#, orientation = "horizontal")
cbar.set_label(labeltext)
ax3 = fig.add_subplot(py"$(gs)[2,0:3]", sharex=ax4)
#ax3.axhline(color="black", alpha=0.6)
ax3.plot(fx1.time, flux, color="black")
ax3.set_ylim((-0.1, 0.2))
ax3.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax3.set_ylabel(labelbottom)
ax3.xaxis.set_major_locator(majorloc)
ax3.xaxis.set_minor_locator(minorloc)
ax3.xaxis.set_major_formatter(date_format)
ax3.grid()

fig.add_artist(lines.Line2D([0.1, 0.9], [0.288, 0.288], color="black"))
fig.add_artist(lines.Line2D([0.413, 0.484], [0.308, 0.308], color="black", alpha=0.6))
fig.text(0.444, 0.293, "b")
fig.add_artist(lines.Line2D([0.665, 0.734], [0.308, 0.308], color="black", alpha=0.6))
fig.text(0.699, 0.293, "c")

ax5 = fig.add_subplot(py"$(gs)[4, 0]")
fig.text(0.1, 0.265, "b)", size="large")
mrdplot = ax5.pcolormesh(time_middle, Dates.value.(x) ./ 10^3, mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax5.set_ylabel(L"$\tau~\mathrm{[s]}$")
ax5.xaxis_date()
ax5.set_xlim(DateTime(2021,05,29,00,00,00), DateTime(2021,05,31,00,00,00))
ax5.set_ylim([0.1, 2e3])
ax5.axhline(Dates.value(Millisecond(ra))/1e3, ls="--", color="black", alpha=0.3)
ax5.fill_betweenx((Dates.value(Millisecond(ra))/1e3, last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax5.xaxis.set_major_locator(pydates.DayLocator(interval=1))
ax5.xaxis.set_minor_locator(pydates.HourLocator(interval=6))
ax5.xaxis.set_major_formatter(date_format)
ax6 = fig.add_subplot(py"$(gs)[4, 2]", sharey=ax5)
fig.text(0.508, 0.265, "c)", size="large")
mrdplot = ax6.pcolormesh(time_middle, Dates.value.(x) ./ 10^3, mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax6.xaxis_date()
ax6.set_xlim(DateTime(2021,06,05,00,00,00), DateTime(2021,06,07,00,00,00))
ax6.set_ylim([0.1, 2e3])
ax6.axhline(Dates.value(Millisecond(ra))/1e3, ls="--", color="black", alpha=0.3)
ax6.fill_betweenx((Dates.value(Millisecond(ra))/1e3, last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax6.tick_params(axis="y", labelleft=false)
ax6.xaxis.set_major_locator(pydates.DayLocator(interval=1))
ax6.xaxis.set_minor_locator(pydates.HourLocator(interval=6))
ax6.xaxis.set_major_formatter(date_format)
#=ax7 = fig.add_subplot(gs[5, 4])
cbar = fig.colorbar(mrdplot, cax=ax7)#, orientation = "horizontal")
cbar.set_label(labeltext)=#
##

#py"$(gs)[2:-2, 1]"

######################################################
###       FIG 6: EVOLUTION MRD T VARIANCE          ###
######################################################
using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings, DataFrames
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
lines = pyimport("matplotlib.lines")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Second(30), #Millisecond(2^10*timestep),
    "KAIJO" => Millisecond(2^9*timestep))


fluxtype = "TT"
fluxfilename = fluxtype

if fluxtype == "wT" || fluxtype == "tke"
    fluxfilename = string(fluxfilename, "_norm.nc")
else
    fluxfilename = string(fluxfilename, "_norm.nc")
end

panelfile1 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxfilename))

(panel1_time_middle, panel1_x, panel1_mrd_data, panel1_meanwind,
    panel1_evalstart, panel1_evalend, panel1_cols_for_mrd,
    panel1_nrpoints, panel1_blocklen, panel1_shift, panel1_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile1)

rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel1_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel1_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel1_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel1_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel1_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel1_time_middle .< DateTime(2021,06,05,19,35,00)

panel1_mrd_data[:, rain_for_2dmrd] .= NaN

panelinst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)

panel1ra = inst2ra[panelinst1]

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[panel1_evalstart.<=tjkmeteo[:, 1].<=panel1_evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

#load data and calc fluxes
inst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)

fluxfile1 = inst2file[inst1]

fx1df = turb.readturbasnetcdf(fluxfile1, panel1_evalstart, panel1_evalend)
#show statistics about missing data
turb.printmissstats(fx1df)
#interpolate missing
fx1df = turb.interpolatemissing(fx1df)
#double rotation
if inst1 == "TJK"
    turb.drdf!(fx1df, periodwise=false)
else
    turb.drdf!(fx1df)
end
turb.missing2nan!(fx1df)
wndspd = gen.movingaverage(sqrt.(fx1df.u .^ 2 .+ fx1df.v .^ 2 .+ fx1df.w .^ 2), 20 * 3600)

fluxpath = "/home/haugened/Documents/data/fluxes/"
fx1file = joinpath(fluxpath, string(lowercase(inst1), "_fx.nc"))

fx1 = turb.readturbasnetcdf(fx1file, panel1_evalstart, panel1_evalend)

fx1 = fx1[1:60:end, :]

rain_for_fx1 =
DateTime(2021,05,22,00,08,00) .< fx1.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx1.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx1.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx1.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx1.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx1.time .< DateTime(2021,06,05,19,35,00)

fx1[rain_for_fx1, 2:end] .= NaN

multfactor = 1

if fluxtype == "wT"
    cmin = -5 #colorbar limits
    cmax = -cmin
    multfactor = 1000
elseif fluxtype == "tke"
    cmin = 0
    cmax = 15
    multfactor = 100
elseif fluxtype == "TT"
    cmin = 0
    cmax = 6
    multfactor = 100
else
    @error("fluxtype not correct")
end

##
#plot
date_format = pydates.DateFormatter("%d.%m.")
majorloc = pydates.DayLocator(interval=2)
minorloc = pydates.DayLocator(interval=1)
fig = PyPlot.figure(figsize=(9, 9))#5.3))
fig.text(0.1, 0.895, "a)", size="large")
gs = gridspec.GridSpec(5, 4, width_ratios=[24, 3, 24, 1], height_ratios=[77, 615, 77, 40, 231], wspace=0.05)
ax4 = fig.add_subplot(py"$(gs)[0,0:3]")
ax4.plot(fx1df.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax4.get_ylim()
ax4.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax4.set_ylabel(L"| \mathbf{u}|~\mathrm{[m~s^{-1}]}")
ax4.tick_params(axis="x", labelbottom=false)
ax4.grid()
ax4.set_ylim(0,8)
ax4.set_xlim(DateTime(2021,05,21,00,00,00), DateTime(2021,06,11,00,00,00))
#fig.text(0.825, 0.48, "1m", rotation=-90, fontweight="bold")
ax1 = fig.add_subplot(py"$(gs)[1,0:3]", sharex=ax4)
#ax1.set_title("5m meteo station")
mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* multfactor, cmap=cramericm.batlow, vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax1.xaxis_date()
ax1.set_ylim([0.1, 2e3])
ax1.tick_params(axis="x", labelbottom=false)
ax1.axhline(Dates.value(Millisecond(panel1ra))./1000, ls="--", color="black", alpha=0.3)
ax1.fill_betweenx((Dates.value(Millisecond(panel1ra))./1000, last(ax1.get_ylim())), first(ax4.get_xlim()), last(ax4.get_xlim()), color="white", alpha=0.4)
ax1.set_ylabel(L"$\tau~\mathrm{[s]}$")
ax2 = fig.add_subplot(py"$(gs)[1, 3]")
cbar = fig.colorbar(mrdplot, cax=ax2)#, orientation = "horizontal")
cbar.set_label(L"$C_{T_s^2} [\mathrm{10^{-2}~K^2}]$")
ax3 = fig.add_subplot(py"$(gs)[2,0:3]", sharex=ax4)
if fluxtype == "TT"
    ax3.plot(fx1.time, fx1.TT, label="1m", color="black")
    ax3.set_ylabel(L"$\overline{\left( T_s' \right)^2}~\mathrm{[K^2]}$")
    ax3.set_ylim(0,0.8)
end
#ax3.legend(loc="center left", bbox_to_anchor=(1,0.5))
ax3.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(0, 1.5 * maximum(filter(!isnan, fx1.TT)), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax3.grid()
ax3.xaxis.set_major_locator(majorloc)
ax3.xaxis.set_minor_locator(minorloc)
ax3.xaxis.set_major_formatter(date_format)
#=
sleep(Second(1))
labels1 = ax3.xaxis.get_majorticklabels()#[item.get_text() for item in ax6.get_xticklabels()]
newlabels = Vector{String}(undef, length(labels1))
newlabels[1] = labels1[1].get_text()
for idx in 1:length(labels1)
    tmp = labels1[idx].get_text()
    if mod(idx, 2)!=0
        newlabels[idx] = ""
    else
        newlabels[idx] = tmp
    end
end
ax3.xaxis.set_ticklabels(newlabels)
=#
fig.add_artist(lines.Line2D([0.1, 0.9], [0.288, 0.288], color="black"))
fig.add_artist(lines.Line2D([0.413, 0.484], [0.308, 0.308], color="black", alpha=0.6))
fig.text(0.444, 0.293, "b")
fig.add_artist(lines.Line2D([0.665, 0.734], [0.308, 0.308], color="black", alpha=0.6))
fig.text(0.699, 0.293, "c")

ax5 = fig.add_subplot(py"$(gs)[4, 0]")
fig.text(0.1, 0.265, "b)", size="large")
mrdplot = ax5.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 10^3, panel1_mrd_data .* multfactor, cmap=cramericm.batlow, vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax5.set_ylabel(L"$\tau~\mathrm{[s]}$")
ax5.xaxis_date()
ax5.set_xlim(DateTime(2021,05,29,00,00,00), DateTime(2021,05,31,00,00,00))
ax5.set_ylim([0.1, 2e3])
ax5.axhline(Dates.value(Millisecond(panel1ra))/1e3, ls="--", color="black", alpha=0.3)
ax5.fill_betweenx((Dates.value(Millisecond(panel1ra))/1e3, last(ax5.get_ylim())), first(ax5.get_xlim()), last(ax5.get_xlim()), color="white", alpha=0.4)
ax5.xaxis.set_major_locator(pydates.DayLocator(interval=1))
ax5.xaxis.set_minor_locator(pydates.HourLocator(interval=6))
ax5.xaxis.set_major_formatter(date_format)
ax6 = fig.add_subplot(py"$(gs)[4, 2]", sharey=ax5)
fig.text(0.508, 0.265, "c)", size="large")
mrdplot = ax6.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 10^3, panel1_mrd_data .* multfactor, cmap=cramericm.batlow, vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax6.xaxis_date()
ax6.set_xlim(DateTime(2021,06,05,00,00,00), DateTime(2021,06,07,00,00,00))
ax6.set_ylim([0.1, 2e3])
ax6.axhline(Dates.value(Millisecond(panel1ra))/1e3, ls="--", color="black", alpha=0.3)
ax6.fill_betweenx((Dates.value(Millisecond(panel1ra))/1e3, last(ax6.get_ylim())), first(ax6.get_xlim()), last(ax6.get_xlim()), color="white", alpha=0.4)
ax6.tick_params(axis="y", labelleft=false)
ax6.xaxis.set_major_locator(pydates.DayLocator(interval=1))
ax6.xaxis.set_minor_locator(pydates.HourLocator(interval=6))
ax6.xaxis.set_major_formatter(date_format)
#=ax7 = fig.add_subplot(gs[5, 4])
cbar = fig.colorbar(mrdplot, cax=ax7)#, orientation = "horizontal")
cbar.set_label(labeltext)=#
##

######################################################
###       FIG 7: MOVING WINDOW MRD FOR TKE        ###
######################################################

using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
lines = pyimport("matplotlib.lines")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

readfilename = joinpath(datapath, "2dmrd", "db_tot", "t2lcsat_tke_norm.nc")

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Second(30), #Millisecond(2^10*timestep),
    "KAIJO" => Second(400))

(time_middle, x, mrd_data, meanwind, evalstart, evalend,
    cols_for_mrd, nrpoints, blocklen, shift, sourcefile) =
    MRD.read2dmrdfromnetcdf(readfilename)

rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< time_middle .< DateTime(2021,06,05,19,35,00)
    
mrd_data[:, rain_for_2dmrd] .= NaN

fxavgperiod = blocklen * (evalend - evalstart)

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[evalstart.<=tjkmeteo[:, 1].<=evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

#load data and calc fluxes
inst = MRD.instrumentnamefromsourcefile(sourcefile)
fluxfile = inst2file[inst]
ra = inst2ra[inst]

fluxdf = turb.readturbasnetcdf(fluxfile, evalstart, evalend)
#show statistics about missing data
turb.printmissstats(fluxdf)
#interpolate missing
fluxdf = turb.interpolatemissing(fluxdf)
#double rotation
if inst == "TJK"
    turb.drdf!(fluxdf, periodwise=false)
else
    turb.drdf!(fluxdf)
end
turb.missing2nan!(fluxdf)
wndspd = gen.movingaverage(sqrt.(fluxdf.u .^ 2 .+ fluxdf.v .^ 2 .+ fluxdf.w .^ 2), 20 * 3600)
fx1_raw = turb.turbflux(fluxdf, ra)
fx1 = turb.avgflux(fx1_raw, fxavgperiod)
fx1_raw = nothing
#Obukhov length
#L1 = turb.obukhov(fx1, blocklen * (evalend - evalstart))

fx1 = fx1[1:1200:end, :]

rain_for_fx =
DateTime(2021,05,22,00,08,00) .< fx1.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx1.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx1.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx1.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx1.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx1.time .< DateTime(2021,06,05,19,35,00)

fx1[rain_for_fx, 2:end] .= NaN

firstbiggerra = findfirst(x -> x > ra, Millisecond.(Dates.value.(x[:, round(Int, size(x, 2)/2)])))

tke_total = zeros(Float64, size(mrd_data, 2))
for i in 1:length(tke_total)
    tke_total[i] = sum(mrd_data[1:firstbiggerra-1, i])
end

##
#plot temporal scale (y)
cmin = 0 #colorbar limits
cmax = 15
date_format = pydates.DateFormatter("%d.%m.")
majorloc = pydates.DayLocator(interval=2)
minorloc = pydates.DayLocator(interval=1)
fig = PyPlot.figure(figsize=(9, 9))
fig.text(0.1, 0.895, "a)", size="large")
gs = gridspec.GridSpec(5, 4, width_ratios=[24, 3, 24, 1], height_ratios=[77, 615, 77, 40, 231], wspace=0.05)
ax4 = fig.add_subplot(py"$(gs)[0,0:3]")
ax4.plot(fluxdf.time[1:20*600:end], wndspd[1:20*600:end], color="black")
yl = ax4.get_ylim()
ax4.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax4.set_ylabel(L"| \mathbf{u}|~\mathrm{[m~s^{-1}]}")
#ax4.set_ylabel("winddir [°]")
ax4.tick_params(axis="x", labelbottom=false)
ax4.grid()
ax4.set_ylim(0, 8)
ax = fig.add_subplot(py"$(gs)[1,0:3]", sharex=ax4)
mrdplot = ax.pcolormesh(time_middle, Dates.value.(x) ./ 10^3, mrd_data .*100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax.set_ylabel(L"$\tau~\mathrm{[s]}$")
ax.xaxis_date()
ax.set_ylim([0.1, 2e3])
ax.axhline(Dates.value(Millisecond(ra))./1e3, ls="--", color="black", alpha=0.3)
ax.fill_betweenx((Dates.value(Millisecond(ra))./1e3, last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax.tick_params(axis="x", labelbottom=false)
ax2 = fig.add_subplot(py"$(gs)[1, 3]")
cbar = fig.colorbar(mrdplot, cax=ax2)#, orientation = "horizontal")
cbar.set_label(L"$C_{e} [\mathrm{10^{-2}~J~kg^{-1}}]$")
ax3 = fig.add_subplot(py"$(gs)[2,0:3]", sharex=ax4)
#ax3.axhline(color="black", alpha=0.6)
#ax3.plot(pydates.date2num.(time_middle), tke_total, color="black")
ax3.plot(fx1.time, fx1.tke, color="black")
#ax3.set_ylim((0, 30))
ax3.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(0, 1.5 * maximum(filter(!isnan, tke_total)), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax3.set_ylabel(L"$e~\mathrm{[J~kg^{-1}]}$")
ax3.grid()
ax3.set_ylim(0,1.5)
ax3.xaxis.set_major_locator(majorloc)
ax3.xaxis.set_minor_locator(minorloc)
ax3.xaxis.set_major_formatter(date_format)

fig.add_artist(lines.Line2D([0.1, 0.9], [0.288, 0.288], color="black"))
fig.add_artist(lines.Line2D([0.413, 0.484], [0.308, 0.308], color="black", alpha=0.6))
fig.text(0.444, 0.293, "b")
fig.add_artist(lines.Line2D([0.665, 0.734], [0.308, 0.308], color="black", alpha=0.6))
fig.text(0.699, 0.293, "c")

ax5 = fig.add_subplot(py"$(gs)[4, 0]")
fig.text(0.1, 0.265, "b)", size="large")
mrdplot = ax5.pcolormesh(time_middle, Dates.value.(x) ./ 10^3, mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax5.set_ylabel(L"$\tau~\mathrm{[s]}$")
ax5.xaxis_date()
ax5.set_xlim(DateTime(2021,05,29,00,00,00), DateTime(2021,05,31,00,00,00))
ax5.set_ylim([0.1, 2e3])
ax5.axhline(Dates.value(Millisecond(ra))/1e3, ls="--", color="black", alpha=0.3)
ax5.fill_betweenx((Dates.value(Millisecond(ra))/1e3, last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax5.xaxis.set_major_locator(pydates.DayLocator(interval=1))
ax5.xaxis.set_minor_locator(pydates.HourLocator(interval=6))
ax5.xaxis.set_major_formatter(date_format)
ax6 = fig.add_subplot(py"$(gs)[4, 2]", sharey=ax5)
fig.text(0.508, 0.265, "c)", size="large")
mrdplot = ax6.pcolormesh(time_middle, Dates.value.(x) ./ 10^3, mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax6.xaxis_date()
ax6.set_xlim(DateTime(2021,06,05,00,00,00), DateTime(2021,06,07,00,00,00))
ax6.set_ylim([0.1, 2e3])
ax6.axhline(Dates.value(Millisecond(ra))/1e3, ls="--", color="black", alpha=0.3)
ax6.fill_betweenx((Dates.value(Millisecond(ra))/1e3, last(ax.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax6.tick_params(axis="y", labelleft=false)
ax6.xaxis.set_major_locator(pydates.DayLocator(interval=1))
ax6.xaxis.set_minor_locator(pydates.HourLocator(interval=6))
ax6.xaxis.set_major_formatter(date_format)
#=ax7 = fig.add_subplot(gs[5, 4])
cbar = fig.colorbar(mrdplot, cax=ax7)#, orientation = "horizontal")
cbar.set_label(labeltext)=#
##

######################################################
###    FIG 8: LINEAR MRD WIND DIR FOR T2LCSAT      ###
######################################################

using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings, DataFrames
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

readfilename1 = joinpath(datapath, "2dmrd", "db_tot", "t2lcsat_wT_norm.nc")
readfilename2 = joinpath(datapath, "2dmrd", "db_tot", "t2lcsat_TT_norm.nc")
readfilename3 = joinpath(datapath, "2dmrd", "db_tot", "t2lcsat_tke_norm.nc")

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Millisecond(2^10*timestep),
    "KAIJO" => Millisecond(2^8*timestep))


(time_middle1, x1, mrd_data1, meanwind1, evalstart1, evalend1,
    cols_for_mrd1, nrpoints1, blocklen1, shift1, sourcefile1) =
    MRD.read2dmrdfromnetcdf(readfilename1)

(time_middle2, x2, mrd_data2, meanwind2, evalstart2, evalend2,
    cols_for_mrd2, nrpoints2, blocklen2, shift2, sourcefile2) =
    MRD.read2dmrdfromnetcdf(readfilename2)

(time_middle3, x3, mrd_data3, meanwind3, evalstart3, evalend3,
    cols_for_mrd3, nrpoints3, blocklen3, shift3, sourcefile3) =
    MRD.read2dmrdfromnetcdf(readfilename3)

rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< time_middle1 .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< time_middle1 .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< time_middle1 .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< time_middle1 .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< time_middle1 .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< time_middle1 .< DateTime(2021,06,05,19,35,00)
    
mrd_data1[:, rain_for_2dmrd] .= NaN
mrd_data2[:, rain_for_2dmrd] .= NaN
mrd_data3[:, rain_for_2dmrd] .= NaN

figtitle1 = MRD.create2dmrdtitle(sourcefile1, evalstart1, evalend1)
figtitle2 = MRD.create2dmrdtitle(sourcefile2, evalstart2, evalend2)
figtitle3 = MRD.create2dmrdtitle(sourcefile3, evalstart3, evalend3)

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[evalstart1.<=tjkmeteo[:, 1].<=evalend1, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

timemiddle_start = [DateTime(2021, 05, 01, 00, 00, 00)]
timemiddle_end = DateTime(2021, 07, 01, 00, 00, 00)
daytime_start = Time(08,00,00)
daytime_end = Time(18,00,00)

mrd_data1[:, daytime_end .< Time.(time_middle1) .|| Time.(time_middle1) .< daytime_start] .= NaN
mrd_data2[:, daytime_end .< Time.(time_middle2) .|| Time.(time_middle2) .< daytime_start] .= NaN
mrd_data3[:, daytime_end .< Time.(time_middle3) .|| Time.(time_middle3) .< daytime_start] .= NaN

mrd_data_up1 = copy(mrd_data1[:, timemiddle_start.<=time_middle1.<=timemiddle_end])
mrd_data_down1 = copy(mrd_data1[:, timemiddle_start.<=time_middle1.<=timemiddle_end])
time_middle_sub1 = copy(time_middle1[timemiddle_start.<=time_middle1.<=timemiddle_end])
x_sub1 = copy(x1[:, timemiddle_start.<=time_middle1.<=timemiddle_end])

mrd_data_up2 = copy(mrd_data2[:, timemiddle_start.<=time_middle2.<=timemiddle_end])
mrd_data_down2 = copy(mrd_data2[:, timemiddle_start.<=time_middle2.<=timemiddle_end])
time_middle_sub2 = copy(time_middle2[timemiddle_start.<=time_middle2.<=timemiddle_end])
x_sub2 = copy(x2[:, timemiddle_start.<=time_middle2.<=timemiddle_end])

mrd_data_up3 = copy(mrd_data3[:, timemiddle_start.<=time_middle3.<=timemiddle_end])
mrd_data_down3 = copy(mrd_data3[:, timemiddle_start.<=time_middle3.<=timemiddle_end])
time_middle_sub3 = copy(time_middle3[timemiddle_start.<=time_middle3.<=timemiddle_end])
x_sub3 = copy(x3[:, timemiddle_start.<=time_middle3.<=timemiddle_end])

#----------------------------------------------------
#take only specific wind directions if wanted
disallowmissing!(tjkdata)
windclass_2dmrd1 = fill(NaN, size(time_middle_sub1, 1))
windclass_2dmrd2 = fill(NaN, size(time_middle_sub2, 1))
windclass_2dmrd3 = fill(NaN, size(time_middle_sub3, 1))
turb.extendtofinertimeseries!(windclass_2dmrd1, time_middle_sub1, windclass, tjkdata.time)
turb.extendtofinertimeseries!(windclass_2dmrd2, time_middle_sub2, windclass, tjkdata.time)
turb.extendtofinertimeseries!(windclass_2dmrd3, time_middle_sub3, windclass, tjkdata.time)

upwind1 = windclass_2dmrd1 .== 2   #1 = up valley wind; 2 = down valley wind
downwind1 = windclass_2dmrd1 .== 7   #1 = up valley wind; 2 = down valley wind
mrd_data_up1[:, .!upwind1] .= NaN
mrd_data_down1[:, .!downwind1] .= NaN

upwind2 = windclass_2dmrd2 .== 2   #1 = up valley wind; 2 = down valley wind
downwind2 = windclass_2dmrd2 .== 7   #1 = up valley wind; 2 = down valley wind
mrd_data_up2[:, .!upwind2] .= NaN
mrd_data_down2[:, .!downwind2] .= NaN

upwind3 = windclass_2dmrd3 .== 2   #1 = up valley wind; 2 = down valley wind
downwind3 = windclass_2dmrd3 .== 7   #1 = up valley wind; 2 = down valley wind
mrd_data_up3[:, .!upwind3] .= NaN
mrd_data_down3[:, .!downwind3] .= NaN

#----------------------------------------------------

per1_lin_mrd_up_med1 = fill(NaN, size(mrd_data_up1, 1))
per1_lin_mrd_up_q11 = fill(NaN, size(mrd_data_up1, 1))
per1_lin_mrd_up_q31 = fill(NaN, size(mrd_data_up1, 1))
per1_lin_mrd_down_med1 = fill(NaN, size(mrd_data_down1, 1))
per1_lin_mrd_down_q11 = fill(NaN, size(mrd_data_down1, 1))
per1_lin_mrd_down_q31 = fill(NaN, size(mrd_data_down1, 1))

per1_lin_mrd_up_med2 = fill(NaN, size(mrd_data_up2, 1))
per1_lin_mrd_up_q12 = fill(NaN, size(mrd_data_up2, 1))
per1_lin_mrd_up_q32 = fill(NaN, size(mrd_data_up2, 1))
per1_lin_mrd_down_med2 = fill(NaN, size(mrd_data_down2, 1))
per1_lin_mrd_down_q12 = fill(NaN, size(mrd_data_down2, 1))
per1_lin_mrd_down_q32 = fill(NaN, size(mrd_data_down2, 1))

per1_lin_mrd_up_med3 = fill(NaN, size(mrd_data_up3, 1))
per1_lin_mrd_up_q13 = fill(NaN, size(mrd_data_up3, 1))
per1_lin_mrd_up_q33 = fill(NaN, size(mrd_data_up3, 1))
per1_lin_mrd_down_med3 = fill(NaN, size(mrd_data_down3, 1))
per1_lin_mrd_down_q13 = fill(NaN, size(mrd_data_down3, 1))
per1_lin_mrd_down_q33 = fill(NaN, size(mrd_data_down3, 1))

per2_lin_mrd_up_med1 = fill(NaN, size(mrd_data_up1, 1))
per2_lin_mrd_up_q11 = fill(NaN, size(mrd_data_up1, 1))
per2_lin_mrd_up_q31 = fill(NaN, size(mrd_data_up1, 1))
per2_lin_mrd_down_med1 = fill(NaN, size(mrd_data_down1, 1))
per2_lin_mrd_down_q11 = fill(NaN, size(mrd_data_down1, 1))
per2_lin_mrd_down_q31 = fill(NaN, size(mrd_data_down1, 1))

per2_lin_mrd_up_med2 = fill(NaN, size(mrd_data_up2, 1))
per2_lin_mrd_up_q12 = fill(NaN, size(mrd_data_up2, 1))
per2_lin_mrd_up_q32 = fill(NaN, size(mrd_data_up2, 1))
per2_lin_mrd_down_med2 = fill(NaN, size(mrd_data_down2, 1))
per2_lin_mrd_down_q12 = fill(NaN, size(mrd_data_down2, 1))
per2_lin_mrd_down_q32 = fill(NaN, size(mrd_data_down2, 1))

per2_lin_mrd_up_med3 = fill(NaN, size(mrd_data_up3, 1))
per2_lin_mrd_up_q13 = fill(NaN, size(mrd_data_up3, 1))
per2_lin_mrd_up_q33 = fill(NaN, size(mrd_data_up3, 1))
per2_lin_mrd_down_med3 = fill(NaN, size(mrd_data_down3, 1))
per2_lin_mrd_down_q13 = fill(NaN, size(mrd_data_down3, 1))
per2_lin_mrd_down_q33 = fill(NaN, size(mrd_data_down3, 1))

#split in different times
times_start = [DateTime(2021,05,29,00,00,00), DateTime(2021,06,05,00,00,00)]
times_end = [DateTime(2021,05,31,00,00,00), DateTime(2021,06,07,00,00,00)]

for i in 1:size(mrd_data_up1, 1)
    try
        (per1_lin_mrd_up_q11[i], per1_lin_mrd_up_med1[i], per1_lin_mrd_up_q31[i]) = quantile(filter(!isnan, mrd_data_up1[i, times_start[1] .<= time_middle_sub1 .<= times_end[1]]), [0.25, 0.5, 0.75])
    catch e
    end
    try
        (per1_lin_mrd_down_q11[i], per1_lin_mrd_down_med1[i], per1_lin_mrd_down_q31[i]) = quantile(filter(!isnan, mrd_data_down1[i, times_start[1] .<= time_middle_sub1 .<= times_end[1]]), [0.25, 0.5, 0.75])
    catch e
    end
end

for i in 1:size(mrd_data_up2, 1)
    try
        (per1_lin_mrd_up_q12[i], per1_lin_mrd_up_med2[i], per1_lin_mrd_up_q32[i]) = quantile(filter(!isnan, mrd_data_up2[i, times_start[1] .<= time_middle_sub2 .<= times_end[1]]), [0.25, 0.5, 0.75])
    catch e
    end
    try
        (per1_lin_mrd_down_q12[i], per1_lin_mrd_down_med2[i], per1_lin_mrd_down_q32[i]) = quantile(filter(!isnan, mrd_data_down2[i, times_start[1] .<= time_middle_sub2 .<= times_end[1]]), [0.25, 0.5, 0.75])
    catch e
    end
end

for i in 1:size(mrd_data_up3, 1)
    try
        (per1_lin_mrd_up_q13[i], per1_lin_mrd_up_med3[i], per1_lin_mrd_up_q33[i]) = quantile(filter(!isnan, mrd_data_up3[i, times_start[1] .<= time_middle_sub3 .<= times_end[1]]), [0.25, 0.5, 0.75])
    catch e
    end
    try
        (per1_lin_mrd_down_q13[i], per1_lin_mrd_down_med3[i], per1_lin_mrd_down_q33[i]) = quantile(filter(!isnan, mrd_data_down3[i, times_start[1] .<= time_middle_sub3 .<= times_end[1]]), [0.25, 0.5, 0.75])
    catch e
    end
end


for i in 1:size(mrd_data_up1, 1)
    try
        (per2_lin_mrd_up_q11[i], per2_lin_mrd_up_med1[i], per2_lin_mrd_up_q31[i]) = quantile(filter(!isnan, mrd_data_up1[i, times_start[2] .<= time_middle_sub1 .<= times_end[2]]), [0.25, 0.5, 0.75])
    catch e
    end
    try
        (per2_lin_mrd_down_q11[i], per2_lin_mrd_down_med1[i], per2_lin_mrd_down_q31[i]) = quantile(filter(!isnan, mrd_data_down1[i, times_start[2] .<= time_middle_sub1 .<= times_end[2]]), [0.25, 0.5, 0.75])
    catch e
    end
end

for i in 1:size(mrd_data_up2, 1)
    try
        (per2_lin_mrd_up_q12[i], per2_lin_mrd_up_med2[i], per2_lin_mrd_up_q32[i]) = quantile(filter(!isnan, mrd_data_up2[i, times_start[2] .<= time_middle_sub2 .<= times_end[2]]), [0.25, 0.5, 0.75])
    catch e
    end
    try
        (per2_lin_mrd_down_q12[i], per2_lin_mrd_down_med2[i], per2_lin_mrd_down_q32[i]) = quantile(filter(!isnan, mrd_data_down2[i, times_start[2] .<= time_middle_sub2 .<= times_end[2]]), [0.25, 0.5, 0.75])
    catch e
    end
end

for i in 1:size(mrd_data_up3, 1)
    try
        (per2_lin_mrd_up_q13[i], per2_lin_mrd_up_med3[i], per2_lin_mrd_up_q33[i]) = quantile(filter(!isnan, mrd_data_up3[i, times_start[2] .<= time_middle_sub3 .<= times_end[2]]), [0.25, 0.5, 0.75])
    catch e
    end
    try
        (per2_lin_mrd_down_q13[i], per2_lin_mrd_down_med3[i], per2_lin_mrd_down_q33[i]) = quantile(filter(!isnan, mrd_data_down3[i, times_start[2] .<= time_middle_sub3 .<= times_end[2]]), [0.25, 0.5, 0.75])
    catch e
    end
end

##
fig = PyPlot.figure(figsize=(8,7))
gs = gridspec.GridSpec(3,2)
ax2 = fig.add_subplot(gs[1,1])
ax2.set_title(L"$29.5.-31.5.,~f_s=0.75$")
#ax2.set_title("5 m (meteo station)")
PyPlot.axhline(0, color="grey", alpha=0.5)
ax2.plot(Dates.value.(x1[:, 1]) ./ 1000, per1_lin_mrd_down_med1 .* 1e3, color=cm.tab10(7), label="down valley", ls="dotted")
ax2.fill_between(Dates.value.(x1[:, 1]) ./ 1000, per1_lin_mrd_down_q11 .* 1e3, per1_lin_mrd_down_q31 .* 1e3, color=cm.tab10(7), alpha=0.5)
ax2.plot(Dates.value.(x1[:, 1]) ./ 1000, per1_lin_mrd_up_med1 .* 1e3, color=cm.tab10(1), label="up valley")
ax2.fill_between(Dates.value.(x1[:, 1]) ./ 1000, per1_lin_mrd_up_q11 .* 1e3, per1_lin_mrd_up_q31 .* 1e3, color=cm.tab10(1), alpha=0.5)
PyPlot.xscale("log")
ax2.grid(alpha=0.5)
#ax2.set_xlabel(L"$\tau~\mathrm{[s]}$")
ax2.set_ylabel(L"$C_{wT_s}~[10^{-3}~\mathrm{K~m~s^{-1}}]$")
PyPlot.axvline(102, 0, 1, alpha=0.5, color="black", ls="--")
ax2.tick_params(axis="x", labelbottom=false)
#ax2.legend()
ax3 = fig.add_subplot(gs[2,1], sharex=ax2)
PyPlot.axhline(0, color="grey", alpha=0.5)
ax3.plot(Dates.value.(x2[:, 1]) ./ 1000, per1_lin_mrd_down_med2 .* 1e2, color=cm.tab10(7), label="down valley", ls="dotted")
ax3.fill_between(Dates.value.(x2[:, 1]) ./ 1000, per1_lin_mrd_down_q12 .* 1e2, per1_lin_mrd_down_q32 .* 1e2, color=cm.tab10(7), alpha=0.5)
ax3.plot(Dates.value.(x2[:, 1]) ./ 1000, per1_lin_mrd_up_med2 .* 1e2, color=cm.tab10(1), label="up valley")
ax3.fill_between(Dates.value.(x2[:, 1]) ./ 1000, per1_lin_mrd_up_q12 .* 1e2, per1_lin_mrd_up_q32 .* 1e2, color=cm.tab10(1), alpha=0.5)
PyPlot.xscale("log")
ax3.grid(alpha=0.5)
#ax3.set_xlabel(L"$\tau~\mathrm{[s]}$")
ax3.set_ylabel(L"$C_{T_s^2}~[10^{-2}~\mathrm{K^{2}}]$")
PyPlot.axvline(102, 0, 1, alpha=0.5, color="black", ls="--")
ax3.tick_params(axis="x", labelbottom=false)
ax4 = fig.add_subplot(gs[3,1], sharex=ax2)
PyPlot.axhline(0, color="grey", alpha=0.5)
ax4.plot(Dates.value.(x3[:, 1]) ./ 1000, per1_lin_mrd_down_med3 .* 1e2, color=cm.tab10(7), label="down valley", ls="dotted")
ax4.fill_between(Dates.value.(x3[:, 1]) ./ 1000, per1_lin_mrd_down_q13 .* 1e2, per1_lin_mrd_down_q33 .* 1e2, color=cm.tab10(7), alpha=0.5)
ax4.plot(Dates.value.(x3[:, 1]) ./ 1000, per1_lin_mrd_up_med3 .* 1e2, color=cm.tab10(1), label="up valley")
ax4.fill_between(Dates.value.(x3[:, 1]) ./ 1000, per1_lin_mrd_up_q13 .* 1e2, per1_lin_mrd_up_q33 .* 1e2, color=cm.tab10(1), alpha=0.5)
PyPlot.xscale("log")
ax4.grid(alpha=0.5)
ax4.set_xlabel(L"$\tau~\mathrm{[s]}$")
ax4.set_ylabel(L"$C_{e}~[10^{-2}~\mathrm{J~kg^{-1}}]$")
PyPlot.axvline(102, 0, 1, alpha=0.5, color="black", ls="--")

ax5 = fig.add_subplot(gs[1,2], sharey=ax2)
ax5.set_title(L"$5.6.-7.6.,~f_s=0.16$")
#ax2.set_title("5 m (meteo station)")
PyPlot.axhline(0, color="grey", alpha=0.5)
ax5.plot(Dates.value.(x1[:, 1]) ./ 1000, per2_lin_mrd_down_med1 .* 1e3, color=cm.tab10(7), label="down valley", ls="dotted")
ax5.fill_between(Dates.value.(x1[:, 1]) ./ 1000, per2_lin_mrd_down_q11 .* 1e3, per2_lin_mrd_down_q31 .* 1e3, color=cm.tab10(7), alpha=0.5)
ax5.plot(Dates.value.(x1[:, 1]) ./ 1000, per2_lin_mrd_up_med1 .* 1e3, color=cm.tab10(1), label="up valley")
ax5.fill_between(Dates.value.(x1[:, 1]) ./ 1000, per2_lin_mrd_up_q11 .* 1e3, per2_lin_mrd_up_q31 .* 1e3, color=cm.tab10(1), alpha=0.5)
PyPlot.xscale("log")
ax5.grid(alpha=0.5)
#ax2.set_xlabel(L"$\tau~\mathrm{[s]}$")
#ax5.set_ylabel(L"$C_{wT_s}~[10^{-3}~\mathrm{K~m~s^{-1}}]$")
PyPlot.axvline(102, 0, 1, alpha=0.5, color="black", ls="--")
ax5.set_ylim(-5,12)
ax5.tick_params(axis="x", labelbottom=false)
ax5.tick_params(axis="y", labelleft=false)
#ax2.legend()
ax6 = fig.add_subplot(gs[2,2], sharex=ax5, sharey=ax3)
ax6.plot(Dates.value.(x2[:, 1]) ./ 1000, per2_lin_mrd_down_med2 .* 1e2, color=cm.tab10(7), label="down valley", ls="dotted")
ax6.fill_between(Dates.value.(x2[:, 1]) ./ 1000, per2_lin_mrd_down_q12 .* 1e2, per2_lin_mrd_down_q32 .* 1e2, color=cm.tab10(7), alpha=0.5)
ax6.plot(Dates.value.(x2[:, 1]) ./ 1000, per2_lin_mrd_up_med2 .* 1e2, color=cm.tab10(1), label="up valley")
ax6.fill_between(Dates.value.(x2[:, 1]) ./ 1000, per2_lin_mrd_up_q12 .* 1e2, per2_lin_mrd_up_q32 .* 1e2, color=cm.tab10(1), alpha=0.5)
PyPlot.xscale("log")
ax6.grid(alpha=0.5)
#ax3.set_xlabel(L"$\tau~\mathrm{[s]}$")
#ax6.set_ylabel(L"$C_{T_s^2}~[10^{-2}~\mathrm{K^{2}}]$")
PyPlot.axvline(102, 0, 1, alpha=0.5, color="black", ls="--")
ax6.set_ylim(0,5.5)
ax6.tick_params(axis="x", labelbottom=false)
ax6.tick_params(axis="y", labelleft=false)
ax7 = fig.add_subplot(gs[3,2], sharex=ax5, sharey=ax4)
ax7.plot(Dates.value.(x3[:, 1]) ./ 1000, per2_lin_mrd_down_med3 .* 1e2, color=cm.tab10(7), label="down valley", ls="dotted")
ax7.fill_between(Dates.value.(x3[:, 1]) ./ 1000, per2_lin_mrd_down_q13 .* 1e2, per2_lin_mrd_down_q33 .* 1e2, color=cm.tab10(7), alpha=0.5)
ax7.plot(Dates.value.(x3[:, 1]) ./ 1000, per2_lin_mrd_up_med3 .* 1e2, color=cm.tab10(1), label="up valley")
ax7.fill_between(Dates.value.(x3[:, 1]) ./ 1000, per2_lin_mrd_up_q13 .* 1e2, per2_lin_mrd_up_q33 .* 1e2, color=cm.tab10(1), alpha=0.5)
PyPlot.xscale("log")
ax7.grid(alpha=0.5)
ax7.set_xlabel(L"$\tau~\mathrm{[s]}$")
ax7.tick_params(axis="y", labelleft=false)
ax7.set_ylim(0,11.5)
#ax7.set_ylabel(L"$C_{e}~[10^{-2}~\mathrm{J~kg^{-1}}]$")
PyPlot.axvline(102, 0, 1, alpha=0.5, color="black", ls="--")
##

######################################################
###       FIG 9: PROFILE DECOMPOSITION TKE         ###
######################################################
using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings, DataFrames
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

readfilename = joinpath(datapath, "2dmrd", "db_tot", "t2lcsat_wT_norm.nc")

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Second(30), #Millisecond(2^10*timestep),
    "KAIJO" => Millisecond(2^8*timestep))

fluxtype = "tke"

panelfile1 = joinpath(datapath, "2dmrd", "db_tot", string("tjk_", fluxtype, "_norm.nc"))
panelfile2 = joinpath(datapath, "2dmrd", "db_tot", string("t2ucsat_", fluxtype, "_norm.nc"))
panelfile3 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxtype, "_norm.nc"))
panelfile4 = joinpath(datapath, "2dmrd", "db_tot", string("t2irg_", fluxtype, "_norm.nc"))

(panel1_time_middle, panel1_x, panel1_mrd_data, panel1_meanwind,
    panel1_evalstart, panel1_evalend, panel1_cols_for_mrd,
    panel1_nrpoints, panel1_blocklen, panel1_shift, panel1_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile1)

(panel2_time_middle, panel2_x, panel2_mrd_data, panel2_meanwind,
    panel2_evalstart, panel2_evalend, panel2_cols_for_mrd,
    panel2_nrpoints, panel2_blocklen, panel2_shift, panel2_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile2)

(panel3_time_middle, panel3_x, panel3_mrd_data, panel3_meanwind,
    panel3_evalstart, panel3_evalend, panel3_cols_for_mrd,
    panel3_nrpoints, panel3_blocklen, panel3_shift, panel3_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile3)

(panel4_time_middle, panel4_x, panel4_mrd_data, panel4_meanwind,
    panel4_evalstart, panel4_evalend, panel4_cols_for_mrd,
    panel4_nrpoints, panel4_blocklen, panel4_shift, panel4_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile4)

panel1_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel1_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel1_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel1_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel1_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel1_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel1_time_middle .< DateTime(2021,06,05,19,35,00)

panel2_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel2_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel2_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel2_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel2_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel2_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel2_time_middle .< DateTime(2021,06,05,19,35,00)

panel3_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel3_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel3_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel3_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel3_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel3_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel3_time_middle .< DateTime(2021,06,05,19,35,00)

panel4_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel4_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel4_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel4_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel4_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel4_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel4_time_middle .< DateTime(2021,06,05,19,35,00)

panel1_mrd_data[:, panel1_rain_for_2dmrd] .= NaN
panel2_mrd_data[:, panel2_rain_for_2dmrd] .= NaN
panel3_mrd_data[:, panel3_rain_for_2dmrd] .= NaN
panel4_mrd_data[:, panel4_rain_for_2dmrd] .= NaN

panelinst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
panelinst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
panelinst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
panelinst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

panel1ra = inst2ra[panelinst1]
panel2ra = inst2ra[panelinst2]
panel3ra = inst2ra[panelinst3]
panel4ra = inst2ra[panelinst4]

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[panel1_evalstart.<=tjkmeteo[:, 1].<=panel2_evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

#load data and calc fluxes
inst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
inst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
inst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
inst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

fluxfile2 = inst2file[inst2]
ra2 = inst2ra[inst2]

fx2df = turb.readturbasnetcdf(fluxfile2, panel2_evalstart, panel2_evalend)
#show statistics about missing data
turb.printmissstats(fx2df)
#interpolate missing
fx2df = turb.interpolatemissing(fx2df)
#double rotation
if inst2 == "TJK"
    turb.drdf!(fx2, periodwise=false)
else
    turb.drdf!(fx2df)
end
turb.missing2nan!(fx2df)
wndspd = gen.movingaverage(sqrt.(fx2df.u .^ 2 .+ fx2df.v .^ 2 .+ fx2df.w .^ 2), 20 * 3600)

fluxpath = "/home/haugened/Documents/data/fluxes/"
fx1file = joinpath(fluxpath, string(lowercase(inst1), "_fx.nc"))
fx2file = joinpath(fluxpath, string(lowercase(inst2), "_fx.nc"))
fx3file = joinpath(fluxpath, string(lowercase(inst3), "_fx.nc"))
fx4file = joinpath(fluxpath, string(lowercase(inst4), "_fx.nc"))

fx1 = turb.readturbasnetcdf(fx1file, panel1_evalstart, panel1_evalend)
fx2 = turb.readturbasnetcdf(fx2file, panel2_evalstart, panel2_evalend)
fx3 = turb.readturbasnetcdf(fx3file, panel3_evalstart, panel3_evalend)
fx4 = turb.readturbasnetcdf(fx4file, panel4_evalstart, panel4_evalend)

rain_for_fx1 =
DateTime(2021,05,22,00,08,00) .< fx1.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx1.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx1.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx1.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx1.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx1.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx2 =
DateTime(2021,05,22,00,08,00) .< fx2.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx2.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx2.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx2.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx2.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx2.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx3 =
DateTime(2021,05,22,00,08,00) .< fx3.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx3.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx3.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx3.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx3.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx3.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx4 =
DateTime(2021,05,22,00,08,00) .< fx4.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx4.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx4.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx4.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx4.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx4.time .< DateTime(2021,06,05,19,35,00)

fx1[rain_for_fx1, 2:end] .= NaN
fx2[rain_for_fx2, 2:end] .= NaN
fx3[rain_for_fx3, 2:end] .= NaN
fx4[rain_for_fx4, 2:end] .= NaN

if fluxtype == "wT"
    cmin = -5 #colorbar limits
    cmax = -cmin
elseif fluxtype == "tke"
    cmin = 0
    cmax = 15
else
    @error("fluxtype not correct")
end

##
#plot
fig = PyPlot.figure(figsize=(11,8))
gs = gridspec.GridSpec(6, 2, height_ratios=[5, 15, 15, 15, 15, 8], width_ratios=[50,1])
date_format = pydates.DateFormatter("%d.%m.")
majorloc = pydates.DayLocator(interval=1)
minorloc = pydates.HourLocator([6,12,18])
ax = fig.add_subplot(gs[1, 1])
ax.plot(fx2df.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax.get_ylim()
ax.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax.set_ylabel(L"| \mathbf{u}|_{3\mathrm{m}}~\mathrm{[m~s^{-1}]}")
ax.tick_params(axis="x", labelbottom=false)
ax.grid()
ax.set_ylim(0,5.5)
ax.set_xlim(DateTime(2021,05,29,00,00,00), DateTime(2021,06,07,00,00,00))
fig.text(0.825, 0.73, "5m", rotation=-90, fontweight="bold")
ax1 = fig.add_subplot(gs[2, 1], sharex=ax)
#ax1.set_title("5m meteo station")
if fluxtype == "wT"
    mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
PyPlot.yscale("log")
ax1.xaxis_date()
ax1.set_ylim([0.1, 2e3])
ax1.tick_params(axis="x", labelbottom=false)
ax1.axhline(Dates.value(Millisecond(panel1ra))./1000, ls="--", color="black", alpha=0.3)
ax1.fill_betweenx((Dates.value(Millisecond(panel1ra))./1000, last(ax1.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.57, "3m", rotation=-90, fontweight="bold")
ax2 = fig.add_subplot(gs[3, 1], sharex=ax, sharey=ax1)
#ax2.set_title("3m T2")
if fluxtype == "wT"
    mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
ax2.tick_params(axis="x", labelbottom=false)
#ax2.tick_params(axis="y", labelleft=false)
ax2.axhline(Dates.value(Millisecond(panel2ra))./1000, ls="--", color="black", alpha=0.3)
ax2.fill_betweenx((Dates.value(Millisecond(panel2ra))./1000, last(ax2.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.42, "2m", rotation=-90, fontweight="bold")
ax3 = fig.add_subplot(gs[4, 1], sharex=ax, sharey=ax1)
#ax3.set_title("2m T2")
if fluxtype == "wT"
    mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
ax3.axhline(Dates.value(Millisecond(panel3ra))./1000, ls="--", color="black", alpha=0.3)
ax3.fill_betweenx((Dates.value(Millisecond(panel3ra))./1000, last(ax3.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax3.tick_params(axis="x", labelbottom=false)
ax4 = fig.add_subplot(gs[5, 1], sharex=ax, sharey=ax1)
#ax4.set_title("1m T2")
if fluxtype == "wT"
    mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
#ax4.tick_params(axis="y", labelleft=false)
ax4.axhline(Dates.value(Millisecond(panel4ra))./1000, ls="--", color="black", alpha=0.3)
ax4.fill_betweenx((Dates.value(Millisecond(panel4ra))./1000, last(ax4.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax4.tick_params(axis="x", labelbottom=false)
fig.text(0.825, 0.26, "1m", rotation=-90, fontweight="bold")
ax5 = fig.add_subplot(py"$(gs)[2:-2, 1]")
cbar = fig.colorbar(mrdplot, cax=ax5, orientation="vertical")
if fluxtype == "wT"
    cbar.set_label(L"$C_{wT_s} [10^{-3}~\mathrm{K~m~s^{-1}}]$")
elseif fluxtype == "tke"
    cbar.set_label(L"$C_{e}~[\mathrm{10^{-2}~m^2~s^{-2}}]$")
end
fig.text(0.06, 0.49, L"$\tau~\mathrm{[s]}$", rotation="vertical")
ax6 = fig.add_subplot(gs[6, 1], sharex = ax)
if fluxtype == "wT"
    ax6.plot(fx2.time, fx2.wT, label="3m")
    ax6.plot(fx4.time, fx4.wT, alpha=0.9, label="1m")
    ax6.set_ylabel(L"$\overline{w'T_s'}~[\mathrm{K~m~s^{-1}}]$")
    ax6.set_ylim(-0.05, 0.2)
elseif fluxtype == "tke"
    ax6.plot(fx2.time, fx2.tke, label="3m")
    ax6.plot(fx4.time, fx4.tke, alpha=0.9, label="1m")
    ax6.set_ylabel(L"$e~\mathrm{[m^2~s^{-2}]}$")
    ax6.set_ylim(0,1.2)
end
#yl = ax6.get_ylim()
#ax6.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax6.legend(loc="center left", bbox_to_anchor=(1,0.5))
ax6.grid()
ax6.xaxis.set_major_locator(majorloc)
ax6.xaxis.set_minor_locator(minorloc)
ax6.xaxis.set_major_formatter(date_format)

######################################################
###  FIG 10: PROFILE DECOMPOSITION TKE WITH KAIJO  ###
######################################################
using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings, DataFrames
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Second(30), #Millisecond(2^10*timestep),
    "KAIJO" => Millisecond(2^9*timestep))

fluxtype = "tke"

panelfile1 = joinpath(datapath, "2dmrd", "db_tot", string("t2ucsat_", fluxtype, "_norm.nc"))
panelfile2 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxtype, "_norm.nc"))
panelfile3 = joinpath(datapath, "2dmrd", "db_tot", string("t2irg_", fluxtype, "_norm.nc"))
panelfile4 = joinpath(datapath, "2dmrd", "db_tot", string("kaijo_", fluxtype, "_norm.nc"))

(panel1_time_middle, panel1_x, panel1_mrd_data, panel1_meanwind,
    panel1_evalstart, panel1_evalend, panel1_cols_for_mrd,
    panel1_nrpoints, panel1_blocklen, panel1_shift, panel1_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile1)

(panel2_time_middle, panel2_x, panel2_mrd_data, panel2_meanwind,
    panel2_evalstart, panel2_evalend, panel2_cols_for_mrd,
    panel2_nrpoints, panel2_blocklen, panel2_shift, panel2_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile2)

(panel3_time_middle, panel3_x, panel3_mrd_data, panel3_meanwind,
    panel3_evalstart, panel3_evalend, panel3_cols_for_mrd,
    panel3_nrpoints, panel3_blocklen, panel3_shift, panel3_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile3)

(panel4_time_middle, panel4_x, panel4_mrd_data, panel4_meanwind,
    panel4_evalstart, panel4_evalend, panel4_cols_for_mrd,
    panel4_nrpoints, panel4_blocklen, panel4_shift, panel4_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile4)

panel1_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel1_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel1_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel1_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel1_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel1_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel1_time_middle .< DateTime(2021,06,05,19,35,00)

panel2_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel2_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel2_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel2_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel2_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel2_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel2_time_middle .< DateTime(2021,06,05,19,35,00)

panel3_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel3_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel3_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel3_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel3_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel3_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel3_time_middle .< DateTime(2021,06,05,19,35,00)

panel4_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel4_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel4_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel4_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel4_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel4_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel4_time_middle .< DateTime(2021,06,05,19,35,00)

panel1_mrd_data[:, panel1_rain_for_2dmrd] .= NaN
panel2_mrd_data[:, panel2_rain_for_2dmrd] .= NaN
panel3_mrd_data[:, panel3_rain_for_2dmrd] .= NaN
panel4_mrd_data[:, panel4_rain_for_2dmrd] .= NaN
    
panelinst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
panelinst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
panelinst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
panelinst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

panel1ra = inst2ra[panelinst1]
panel2ra = inst2ra[panelinst2]
panel3ra = inst2ra[panelinst3]
panel4ra = inst2ra[panelinst4]

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[panel1_evalstart.<=tjkmeteo[:, 1].<=panel2_evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

#load data and calc fluxes
inst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
inst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
inst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
inst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

fluxfile1 = inst2file[inst1]
ra1 = inst2ra[inst1]

fx1df = turb.readturbasnetcdf(fluxfile1, panel1_evalstart, panel1_evalend)
#show statistics about missing data
turb.printmissstats(fx1df)
#interpolate missing
fx1df = turb.interpolatemissing(fx1df)
#double rotation
if inst1 == "TJK"
    turb.drdf!(fx1, periodwise=false)
else
    turb.drdf!(fx1df)
end
turb.missing2nan!(fx1df)
wndspd = gen.movingaverage(sqrt.(fx1df.u .^ 2 .+ fx1df.v .^ 2 .+ fx1df.w .^ 2), 20 * 3600)

fluxpath = "/home/haugened/Documents/data/fluxes/"
fx1file = joinpath(fluxpath, string(lowercase(inst1), "_fx.nc"))
fx2file = joinpath(fluxpath, string(lowercase(inst2), "_fx.nc"))
fx3file = joinpath(fluxpath, string(lowercase(inst3), "_fx.nc"))
fx4file = joinpath(fluxpath, string(lowercase(inst4), "_fx.nc"))

fx1 = turb.readturbasnetcdf(fx1file, panel1_evalstart, panel1_evalend)
fx2 = turb.readturbasnetcdf(fx2file, panel2_evalstart, panel2_evalend)
fx3 = turb.readturbasnetcdf(fx3file, panel3_evalstart, panel3_evalend)
fx4 = turb.readturbasnetcdf(fx4file, panel4_evalstart, panel4_evalend)

rain_for_fx1 =
DateTime(2021,05,22,00,08,00) .< fx1.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx1.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx1.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx1.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx1.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx1.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx2 =
DateTime(2021,05,22,00,08,00) .< fx2.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx2.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx2.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx2.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx2.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx2.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx3 =
DateTime(2021,05,22,00,08,00) .< fx3.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx3.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx3.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx3.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx3.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx3.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx4 =
DateTime(2021,05,22,00,08,00) .< fx4.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx4.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx4.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx4.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx4.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx4.time .< DateTime(2021,06,05,19,35,00)

fx1[rain_for_fx1, 2:end] .= NaN
fx2[rain_for_fx2, 2:end] .= NaN
fx3[rain_for_fx3, 2:end] .= NaN
fx4[rain_for_fx4, 2:end] .= NaN

if fluxtype == "wT"
    cmin = -5 #colorbar limits
    cmax = -cmin
elseif fluxtype == "tke"
    cmin = 0
    cmax = 10
else
    @error("fluxtype not correct")
end

##
#plot
fig = PyPlot.figure(figsize=(11,8))
gs = gridspec.GridSpec(6, 2, height_ratios=[5, 15, 15, 15, 15, 8], width_ratios=[50,1])
date_format = pydates.DateFormatter("%d.%m. %H:%M")
majorloc = pydates.HourLocator([0,6,12,18])
minorloc = pydates.HourLocator(interval=1)
ax = fig.add_subplot(gs[1, 1])
ax.plot(fx1df.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax.get_ylim()
ax.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax.set_ylabel(L"| \mathbf{u}|_{3\mathrm{m}}~\mathrm{[m~s^{-1}]}")
ax.tick_params(axis="x", labelbottom=false)
ax.grid()
ax.set_ylim(0,5.5)
ax.set_xlim(DateTime(2021,06,01,06,00,00), DateTime(2021,06,02,12,00,00))
fig.text(0.825, 0.73, "3m", rotation=-90, fontweight="bold")
ax1 = fig.add_subplot(gs[2, 1], sharex=ax)
#ax1.set_title("5m meteo station")
if fluxtype == "wT"
    mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
PyPlot.yscale("log")
ax1.xaxis_date()
ax1.set_ylim([0.1, 2e3])
ax1.tick_params(axis="x", labelbottom=false)
ax1.axhline(Dates.value(Millisecond(panel1ra))./1000, ls="--", color="black", alpha=0.3)
ax1.fill_betweenx((Dates.value(Millisecond(panel1ra))./1000, last(ax1.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.57, "2m", rotation=-90, fontweight="bold")
ax2 = fig.add_subplot(gs[3, 1], sharex=ax, sharey=ax1)
#ax2.set_title("3m T2")
if fluxtype == "wT"
    mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
ax2.tick_params(axis="x", labelbottom=false)
#ax2.tick_params(axis="y", labelleft=false)
ax2.axhline(Dates.value(Millisecond(panel2ra))./1000, ls="--", color="black", alpha=0.3)
ax2.fill_betweenx((Dates.value(Millisecond(panel2ra))./1000, last(ax2.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.42, "1m", rotation=-90, fontweight="bold")
ax3 = fig.add_subplot(gs[4, 1], sharex=ax, sharey=ax1)
#ax3.set_title("2m T2")
if fluxtype == "wT"
    mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
ax3.axhline(Dates.value(Millisecond(panel3ra))./1000, ls="--", color="black", alpha=0.3)
ax3.fill_betweenx((Dates.value(Millisecond(panel3ra))./1000, last(ax3.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax3.tick_params(axis="x", labelbottom=false)
ax4 = fig.add_subplot(gs[5, 1], sharex=ax, sharey=ax1)
#ax4.set_title("1m T2")
if fluxtype == "wT"
    mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
#ax4.tick_params(axis="y", labelleft=false)
ax4.axhline(Dates.value(Millisecond(panel4ra))./1000, ls="--", color="black", alpha=0.3)
ax4.fill_betweenx((Dates.value(Millisecond(panel4ra))./1000, last(ax4.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax4.tick_params(axis="x", labelbottom=false)
fig.text(0.825, 0.26, "0.3m", rotation=-90, fontweight="bold")
ax5 = fig.add_subplot(py"$(gs)[2:-2, 1]")
cbar = fig.colorbar(mrdplot, cax=ax5, orientation="vertical")
if fluxtype == "wT"
    cbar.set_label(L"$C_{wT_s} [10^{-3}~\mathrm{K~m~s^{-1}}]$")
elseif fluxtype == "tke"
    cbar.set_label(L"$C_{e}~[\mathrm{10^{-2}~m^2~s^{-2}}]$")
end
fig.text(0.06, 0.49, L"$\tau~\mathrm{[s]}$", rotation="vertical")
ax6 = fig.add_subplot(gs[6, 1], sharex = ax)
if fluxtype == "wT"
    ax6.plot(fx1.time, fx1.wT, label="3m")
    ax6.plot(fx4.time, fx4.wT, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$\overline{w'T_s'}~[\mathrm{K~m~s^{-1}}]$")
    ax6.set_ylim(-0.1, 0.1)
elseif fluxtype == "tke"
    ax6.plot(fx1.time, fx1.tke, label="3m")
    ax6.plot(fx4.time, fx4.tke, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$e~\mathrm{[m^2~s^{-2}]}$")
    ax6.set_ylim(0,1.0)
end
ax6.legend(loc="center left", bbox_to_anchor=(1,0.5))
ax6.grid()
ax6.xaxis.set_major_locator(majorloc)
ax6.xaxis.set_minor_locator(minorloc)
ax6.xaxis.set_major_formatter(date_format)
sleep(Second(1))
labels1 = ax6.xaxis.get_majorticklabels()#[item.get_text() for item in ax6.get_xticklabels()]
newlabels = Vector{String}(undef, length(labels1))
newlabels[1] = labels1[1].get_text()
for idx in 2:length(labels1)
    tmp = labels1[idx].get_text()
    if tmp[end-4:end-3] != "00" && tmp[3] == '.'
        newlabels[idx] = tmp[end-4:end]
    else
        newlabels[idx] = tmp
    end
end
ax6.xaxis.set_ticklabels(newlabels)
##
######################################################
###       FIG 11: PROFILE DECOMPOSITION wT          ###
######################################################

using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings, DataFrames
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

readfilename = joinpath(datapath, "2dmrd", "db_tot", "t2lcsat_wT_norm.nc")

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Second(30), #Millisecond(2^10*timestep),
    "KAIJO" => Millisecond(2^8*timestep))

fluxtype = "wT"

panelfile1 = joinpath(datapath, "2dmrd", "db_tot", string("tjk_", fluxtype, "_norm.nc"))
panelfile2 = joinpath(datapath, "2dmrd", "db_tot", string("t2ucsat_", fluxtype, "_norm.nc"))
panelfile3 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxtype, "_norm.nc"))
panelfile4 = joinpath(datapath, "2dmrd", "db_tot", string("t2irg_", fluxtype, "_norm.nc"))

(panel1_time_middle, panel1_x, panel1_mrd_data, panel1_meanwind,
    panel1_evalstart, panel1_evalend, panel1_cols_for_mrd,
    panel1_nrpoints, panel1_blocklen, panel1_shift, panel1_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile1)

(panel2_time_middle, panel2_x, panel2_mrd_data, panel2_meanwind,
    panel2_evalstart, panel2_evalend, panel2_cols_for_mrd,
    panel2_nrpoints, panel2_blocklen, panel2_shift, panel2_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile2)

(panel3_time_middle, panel3_x, panel3_mrd_data, panel3_meanwind,
    panel3_evalstart, panel3_evalend, panel3_cols_for_mrd,
    panel3_nrpoints, panel3_blocklen, panel3_shift, panel3_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile3)

(panel4_time_middle, panel4_x, panel4_mrd_data, panel4_meanwind,
    panel4_evalstart, panel4_evalend, panel4_cols_for_mrd,
    panel4_nrpoints, panel4_blocklen, panel4_shift, panel4_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile4)

panel1_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel1_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel1_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel1_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel1_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel1_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel1_time_middle .< DateTime(2021,06,05,19,35,00)

panel2_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel2_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel2_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel2_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel2_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel2_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel2_time_middle .< DateTime(2021,06,05,19,35,00)

panel3_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel3_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel3_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel3_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel3_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel3_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel3_time_middle .< DateTime(2021,06,05,19,35,00)

panel4_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel4_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel4_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel4_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel4_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel4_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel4_time_middle .< DateTime(2021,06,05,19,35,00)

panel1_mrd_data[:, panel1_rain_for_2dmrd] .= NaN
panel2_mrd_data[:, panel2_rain_for_2dmrd] .= NaN
panel3_mrd_data[:, panel3_rain_for_2dmrd] .= NaN
panel4_mrd_data[:, panel4_rain_for_2dmrd] .= NaN

panelinst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
panelinst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
panelinst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
panelinst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

panel1ra = inst2ra[panelinst1]
panel2ra = inst2ra[panelinst2]
panel3ra = inst2ra[panelinst3]
panel4ra = inst2ra[panelinst4]

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[panel1_evalstart.<=tjkmeteo[:, 1].<=panel2_evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

#load data and calc fluxes
inst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
inst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
inst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
inst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

fluxfile2 = inst2file[inst2]
ra2 = inst2ra[inst2]

fx2df = turb.readturbasnetcdf(fluxfile2, panel2_evalstart, panel2_evalend)
#show statistics about missing data
turb.printmissstats(fx2df)
#interpolate missing
fx2df = turb.interpolatemissing(fx2df)
#double rotation
if inst2 == "TJK"
    turb.drdf!(fx2, periodwise=false)
else
    turb.drdf!(fx2df)
end
turb.missing2nan!(fx2df)
wndspd = gen.movingaverage(sqrt.(fx2df.u .^ 2 .+ fx2df.v .^ 2 .+ fx2df.w .^ 2), 20 * 3600)

fluxpath = "/home/haugened/Documents/data/fluxes/"
fx1file = joinpath(fluxpath, string(lowercase(inst1), "_fx.nc"))
fx2file = joinpath(fluxpath, string(lowercase(inst2), "_fx.nc"))
fx3file = joinpath(fluxpath, string(lowercase(inst3), "_fx.nc"))
fx4file = joinpath(fluxpath, string(lowercase(inst4), "_fx.nc"))

fx1 = turb.readturbasnetcdf(fx1file, panel1_evalstart, panel1_evalend)
fx2 = turb.readturbasnetcdf(fx2file, panel2_evalstart, panel2_evalend)
fx3 = turb.readturbasnetcdf(fx3file, panel3_evalstart, panel3_evalend)
fx4 = turb.readturbasnetcdf(fx4file, panel4_evalstart, panel4_evalend)

rain_for_fx1 =
DateTime(2021,05,22,00,08,00) .< fx1.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx1.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx1.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx1.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx1.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx1.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx2 =
DateTime(2021,05,22,00,08,00) .< fx2.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx2.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx2.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx2.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx2.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx2.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx3 =
DateTime(2021,05,22,00,08,00) .< fx3.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx3.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx3.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx3.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx3.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx3.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx4 =
DateTime(2021,05,22,00,08,00) .< fx4.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx4.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx4.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx4.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx4.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx4.time .< DateTime(2021,06,05,19,35,00)

fx1[rain_for_fx1, 2:end] .= NaN
fx2[rain_for_fx2, 2:end] .= NaN
fx3[rain_for_fx3, 2:end] .= NaN
fx4[rain_for_fx4, 2:end] .= NaN

if fluxtype == "wT"
    cmin = -5 #colorbar limits
    cmax = -cmin
elseif fluxtype == "tke"
    cmin = 0
    cmax = 15
else
    @error("fluxtype not correct")
end

##
#plot
fig = PyPlot.figure(figsize=(11,8))
gs = gridspec.GridSpec(6, 2, height_ratios=[5, 15, 15, 15, 15, 8], width_ratios=[50,1])
date_format = pydates.DateFormatter("%d.%m.")
majorloc = pydates.DayLocator(interval=1)
minorloc = pydates.HourLocator([6,12,18])
ax = fig.add_subplot(gs[1, 1])
ax.plot(fx2df.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax.get_ylim()
ax.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax.set_ylabel(L"| \mathbf{u}|_{3\mathrm{m}}~\mathrm{[m~s^{-1}]}")
ax.tick_params(axis="x", labelbottom=false)
ax.grid()
ax.set_ylim(0,5.5)
ax.set_xlim(DateTime(2021,05,29,00,00,00), DateTime(2021,06,07,00,00,00))
fig.text(0.825, 0.73, "5m", rotation=-90, fontweight="bold")
ax1 = fig.add_subplot(gs[2, 1], sharex=ax)
#ax1.set_title("5m meteo station")
if fluxtype == "wT"
    mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
PyPlot.yscale("log")
ax1.xaxis_date()
ax1.set_ylim([0.1, 2e3])
ax1.tick_params(axis="x", labelbottom=false)
ax1.axhline(Dates.value(Millisecond(panel1ra))./1000, ls="--", color="black", alpha=0.3)
ax1.fill_betweenx((Dates.value(Millisecond(panel1ra))./1000, last(ax1.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.57, "3m", rotation=-90, fontweight="bold")
ax2 = fig.add_subplot(gs[3, 1], sharex=ax, sharey=ax1)
#ax2.set_title("3m T2")
if fluxtype == "wT"
    mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
ax2.tick_params(axis="x", labelbottom=false)
#ax2.tick_params(axis="y", labelleft=false)
ax2.axhline(Dates.value(Millisecond(panel2ra))./1000, ls="--", color="black", alpha=0.3)
ax2.fill_betweenx((Dates.value(Millisecond(panel2ra))./1000, last(ax2.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.42, "2m", rotation=-90, fontweight="bold")
ax3 = fig.add_subplot(gs[4, 1], sharex=ax, sharey=ax1)
#ax3.set_title("2m T2")
if fluxtype == "wT"
    mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
ax3.axhline(Dates.value(Millisecond(panel3ra))./1000, ls="--", color="black", alpha=0.3)
ax3.fill_betweenx((Dates.value(Millisecond(panel3ra))./1000, last(ax3.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax3.tick_params(axis="x", labelbottom=false)
ax4 = fig.add_subplot(gs[5, 1], sharex=ax, sharey=ax1)
#ax4.set_title("1m T2")
if fluxtype == "wT"
    mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* 100, cmap="plasma", vmin=cmin, vmax=cmax, shading="auto")
end
#ax4.tick_params(axis="y", labelleft=false)
ax4.axhline(Dates.value(Millisecond(panel4ra))./1000, ls="--", color="black", alpha=0.3)
ax4.fill_betweenx((Dates.value(Millisecond(panel4ra))./1000, last(ax4.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax4.tick_params(axis="x", labelbottom=false)
fig.text(0.825, 0.26, "1m", rotation=-90, fontweight="bold")
ax5 = fig.add_subplot(py"$(gs)[2:-2, 1]")
cbar = fig.colorbar(mrdplot, cax=ax5, orientation="vertical")
if fluxtype == "wT"
    cbar.set_label(L"$C_{wT_s} [10^{-3}~\mathrm{K~m~s^{-1}}]$")
elseif fluxtype == "tke"
    cbar.set_label(L"$C_{e}~[\mathrm{10^{-2}~m^2~s^{-2}}]$")
end
fig.text(0.06, 0.49, L"$\tau~\mathrm{[s]}$", rotation="vertical")
ax6 = fig.add_subplot(gs[6, 1], sharex = ax)
if fluxtype == "wT"
    ax6.plot(fx2.time, fx2.wT, label="3m")
    ax6.plot(fx4.time, fx4.wT, alpha=0.9, label="1m")
    ax6.set_ylabel(L"$\overline{w'T_s'}~[\mathrm{K~m~s^{-1}}]$")
    ax6.set_ylim(-0.05, 0.2)
elseif fluxtype == "tke"
    ax6.plot(fx2.time, fx2.tke, label="3m")
    ax6.plot(fx4.time, fx4.tke, alpha=0.9, label="1m")
    ax6.set_ylabel(L"$e~\mathrm{[m^2~s^{-2}]}$")
    ax6.set_ylim(0,1.2)
end
ax6.legend(loc="center left", bbox_to_anchor=(1,0.5))
ax6.grid()
ax6.xaxis.set_major_locator(majorloc)
ax6.xaxis.set_minor_locator(minorloc)
ax6.xaxis.set_major_formatter(date_format)

######################################################
###         FIG 12: T-VARIANCE WITH KAIJO          ###
######################################################
using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings, DataFrames
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

readfilename = joinpath(datapath, "2dmrd", "db_tot", "kaijo_TT.nc")

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Second(30), #Millisecond(2^10*timestep),
    "KAIJO" => Millisecond(2^9*timestep))

fluxtype = "TT"
fluxfilename = fluxtype

if fluxtype == "wT" || fluxtype == "tke"
    fluxfilename = string(fluxfilename, "_norm.nc")
else
    fluxfilename = string(fluxfilename, "_norm.nc")
end

panelfile1 = joinpath(datapath, "2dmrd", "db_tot", string("t2ucsat_", fluxfilename))
panelfile2 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxfilename))
panelfile3 = joinpath(datapath, "2dmrd", "db_tot", string("t2irg_", fluxfilename))
panelfile4 = joinpath(datapath, "2dmrd", "db_tot", string("kaijo_", fluxfilename))

(panel1_time_middle, panel1_x, panel1_mrd_data, panel1_meanwind,
    panel1_evalstart, panel1_evalend, panel1_cols_for_mrd,
    panel1_nrpoints, panel1_blocklen, panel1_shift, panel1_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile1)

(panel2_time_middle, panel2_x, panel2_mrd_data, panel2_meanwind,
    panel2_evalstart, panel2_evalend, panel2_cols_for_mrd,
    panel2_nrpoints, panel2_blocklen, panel2_shift, panel2_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile2)

(panel3_time_middle, panel3_x, panel3_mrd_data, panel3_meanwind,
    panel3_evalstart, panel3_evalend, panel3_cols_for_mrd,
    panel3_nrpoints, panel3_blocklen, panel3_shift, panel3_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile3)

(panel4_time_middle, panel4_x, panel4_mrd_data, panel4_meanwind,
    panel4_evalstart, panel4_evalend, panel4_cols_for_mrd,
    panel4_nrpoints, panel4_blocklen, panel4_shift, panel4_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile4)

panel1_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel1_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel1_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel1_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel1_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel1_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel1_time_middle .< DateTime(2021,06,05,19,35,00)

panel2_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel2_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel2_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel2_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel2_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel2_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel2_time_middle .< DateTime(2021,06,05,19,35,00)

panel3_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel3_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel3_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel3_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel3_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel3_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel3_time_middle .< DateTime(2021,06,05,19,35,00)

panel4_rain_for_2dmrd =
DateTime(2021,05,22,00,08,00) .< panel4_time_middle .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< panel4_time_middle .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< panel4_time_middle .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< panel4_time_middle .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< panel4_time_middle .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< panel4_time_middle .< DateTime(2021,06,05,19,35,00)

panel1_mrd_data[:, panel1_rain_for_2dmrd] .= NaN
panel2_mrd_data[:, panel2_rain_for_2dmrd] .= NaN
panel3_mrd_data[:, panel3_rain_for_2dmrd] .= NaN
panel4_mrd_data[:, panel4_rain_for_2dmrd] .= NaN

panelinst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
panelinst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
panelinst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
panelinst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

panel1ra = inst2ra[panelinst1]
panel2ra = inst2ra[panelinst2]
panel3ra = inst2ra[panelinst3]
panel4ra = inst2ra[panelinst4]

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[panel1_evalstart.<=tjkmeteo[:, 1].<=panel2_evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

#load data and calc fluxes
inst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
inst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
inst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
inst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

fluxfile2 = inst2file[inst2]
ra2 = inst2ra[inst2]

fx2df = turb.readturbasnetcdf(fluxfile2, panel2_evalstart, panel2_evalend)
#show statistics about missing data
turb.printmissstats(fx2df)
#interpolate missing
fx2df = turb.interpolatemissing(fx2df)
#double rotation
if inst2 == "TJK"
    turb.drdf!(fx2, periodwise=false)
else
    turb.drdf!(fx2df)
end
turb.missing2nan!(fx2df)
wndspd = gen.movingaverage(sqrt.(fx2df.u .^ 2 .+ fx2df.v .^ 2 .+ fx2df.w .^ 2), 20 * 3600)

fluxpath = "/home/haugened/Documents/data/fluxes/"
fx1file = joinpath(fluxpath, string(lowercase(inst1), "_fx.nc"))
fx2file = joinpath(fluxpath, string(lowercase(inst2), "_fx.nc"))
fx3file = joinpath(fluxpath, string(lowercase(inst3), "_fx.nc"))
fx4file = joinpath(fluxpath, string(lowercase(inst4), "_fx.nc"))

fx1 = turb.readturbasnetcdf(fx1file, panel1_evalstart, panel1_evalend)
fx2 = turb.readturbasnetcdf(fx2file, panel2_evalstart, panel2_evalend)
fx3 = turb.readturbasnetcdf(fx3file, panel3_evalstart, panel3_evalend)
fx4 = turb.readturbasnetcdf(fx4file, panel4_evalstart, panel4_evalend)

rain_for_fx1 =
DateTime(2021,05,22,00,08,00) .< fx1.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx1.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx1.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx1.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx1.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx1.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx2 =
DateTime(2021,05,22,00,08,00) .< fx2.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx2.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx2.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx2.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx2.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx2.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx3 =
DateTime(2021,05,22,00,08,00) .< fx3.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx3.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx3.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx3.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx3.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx3.time .< DateTime(2021,06,05,19,35,00)

rain_for_fx4 =
DateTime(2021,05,22,00,08,00) .< fx4.time .< DateTime(2021,05,22,00,38,00) .||
DateTime(2021,05,22,01,30,00) .< fx4.time .< DateTime(2021,05,22,02,10,00) .||
DateTime(2021,05,22,03,55,00) .< fx4.time .< DateTime(2021,05,22,04,55,00) .||
DateTime(2021,05,22,20,48,00) .< fx4.time .< DateTime(2021,05,22,21,48,00) .||
DateTime(2021,05,23,02,19,00) .< fx4.time .< DateTime(2021,05,23,03,28,00) .||
DateTime(2021,06,05,18,35,00) .< fx4.time .< DateTime(2021,06,05,19,35,00)

fx1[rain_for_fx1, 2:end] .= NaN
fx2[rain_for_fx2, 2:end] .= NaN
fx3[rain_for_fx3, 2:end] .= NaN
fx4[rain_for_fx4, 2:end] .= NaN

multfactor = 1

if fluxtype == "wT"
    cmin = -5 #colorbar limits
    cmax = -cmin
    multfactor = 1000
elseif fluxtype == "tke"
    cmin = 0
    cmax = 15
    multfactor = 100
elseif fluxtype == "TT"
    cmin = 0
    cmax = 8
    multfactor = 100
else
    @error("fluxtype not correct")
end

##
#plot
fig = PyPlot.figure(figsize=(11,8))
gs = gridspec.GridSpec(6, 2, height_ratios=[5, 15, 15, 15, 15, 8], width_ratios=[50,1])
date_format = pydates.DateFormatter("%d.%m. %H:%M")
majorloc = pydates.HourLocator([0,6,12,18])
minorloc = pydates.HourLocator(interval=1)
ax = fig.add_subplot(gs[1, 1])
ax.plot(fx2df.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax.get_ylim()
ax.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax.set_ylabel(L"| \mathbf{u}|_{3\mathrm{m}}~\mathrm{[m~s^{-1}]}")
ax.tick_params(axis="x", labelbottom=false)
ax.grid()
ax.set_ylim(0,5.5)
ax.set_xlim(DateTime(2021,06,01,06,00,00), DateTime(2021,06,02,12,00,00))
fig.text(0.825, 0.73, "3m", rotation=-90, fontweight="bold")
ax1 = fig.add_subplot(gs[2, 1], sharex=ax)
#ax1.set_title("5m meteo station")
mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* multfactor, cmap=cramericm.batlow, vmin=cmin, vmax=cmax, shading="auto")
PyPlot.yscale("log")
ax1.xaxis_date()
ax1.set_ylim([0.1, 2e3])
ax1.tick_params(axis="x", labelbottom=false)
ax1.axhline(Dates.value(Millisecond(panel1ra))./1000, ls="--", color="black", alpha=0.3)
ax1.fill_betweenx((Dates.value(Millisecond(panel1ra))./1000, last(ax1.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.57, "2m", rotation=-90, fontweight="bold")
ax2 = fig.add_subplot(gs[3, 1], sharex=ax, sharey=ax1)
#ax2.set_title("3m T2")
mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* multfactor, cmap=cramericm.batlow, vmin=cmin, vmax=cmax, shading="auto")
ax2.tick_params(axis="x", labelbottom=false)
#ax2.tick_params(axis="y", labelleft=false)
ax2.axhline(Dates.value(Millisecond(panel2ra))./1000, ls="--", color="black", alpha=0.3)
ax2.fill_betweenx((Dates.value(Millisecond(panel2ra))./1000, last(ax2.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.42, "1m", rotation=-90, fontweight="bold")
ax3 = fig.add_subplot(gs[4, 1], sharex=ax, sharey=ax1)
#ax3.set_title("2m T2")
mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* multfactor, cmap=cramericm.batlow, vmin=cmin, vmax=cmax, shading="auto")
ax3.axhline(Dates.value(Millisecond(panel3ra))./1000, ls="--", color="black", alpha=0.3)
ax3.fill_betweenx((Dates.value(Millisecond(panel3ra))./1000, last(ax3.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax3.tick_params(axis="x", labelbottom=false)
ax4 = fig.add_subplot(gs[5, 1], sharex=ax, sharey=ax1)
#ax4.set_title("1m T2")
mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* multfactor, cmap=cramericm.batlow, vmin=cmin, vmax=cmax, shading="auto")
#ax4.tick_params(axis="y", labelleft=false)
ax4.axhline(Dates.value(Millisecond(panel4ra))./1000, ls="--", color="black", alpha=0.3)
ax4.fill_betweenx((Dates.value(Millisecond(panel4ra))./1000, last(ax4.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax4.tick_params(axis="x", labelbottom=false)
fig.text(0.825, 0.26, "0.3m", rotation=-90, fontweight="bold")
ax5 = fig.add_subplot(py"$(gs)[2:-2, 1]")
cbar = fig.colorbar(mrdplot, cax=ax5, orientation="vertical")
if fluxtype == "wT"
    cbar.set_label(L"$C_{wT_s} [10^{-3}~\mathrm{K~m~s^{-1}}]$")
elseif fluxtype == "tke"
    cbar.set_label(L"$C_{e}~[\mathrm{10^{-2}~m^2~s^{-2}}]$")
elseif fluxtype == "TT"
    cbar.set_label(L"$C_{T_s^2}~[\mathrm{10^{-2}~K^2}]$")
end
fig.text(0.06, 0.49, L"$\tau~\mathrm{[s]}$", rotation="vertical")
ax6 = fig.add_subplot(gs[6, 1], sharex = ax)
if fluxtype == "wT"
    ax6.plot(fx2.time, fx2.wT, label="2m")
    ax6.plot(fx4.time, fx4.wT, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$\overline{w'T_s'}~[\mathrm{K~m~s^{-1}}]$")
    ax6.set_ylim(-0.05, 0.2)
elseif fluxtype == "tke"
    ax6.plot(fx2.time, fx2.tke, label="2m")
    ax6.plot(fx4.time, fx4.tke, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$e~\mathrm{[m^2~s^{-2}]}$")
    ax6.set_ylim(0,1.2)
elseif fluxtype == "TT"
    ax6.plot(fx2.time, fx2.TT, label="2m")
    ax6.plot(fx4.time, fx4.TT, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$\overline{\left( T_s' \right)^2}~\mathrm{[K^2]}$")
    ax6.set_ylim(0,0.8)
end
ax6.legend(loc="center left", bbox_to_anchor=(1,0.5))
ax6.grid()
ax6.xaxis.set_major_locator(majorloc)
ax6.xaxis.set_minor_locator(minorloc)
ax6.xaxis.set_major_formatter(date_format)
##
######################################################
###             FIG 13: SCREEN STUFF               ###
######################################################

using PyCall, Statistics, LaTeXStrings, Dates
import PyPlot
gridspec = pyimport("matplotlib.gridspec")
mpwidgets = pyimport("matplotlib.widgets")
animation = pyimport("matplotlib.animation")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
    pathtofile = "/home/michi/Documents/slf/data/"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
    pathtofile = "/home/haugened/Documents/data/ir/"
end
include(joinpath(importdir, "src", "ir_evaluation.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .irev
import .gen
PyPlot.pygui(true)

#completing the filepath
pathtofile = joinpath(pathtofile, "210531_143457/")
ncfile = "irdata_09_to_11.nc"
#for plotting
cmperpxl = 0.6
cmperpxlu = 0.606
cmperpxlw = 0.549
framespersec = 29.9
profilexmean = [92, 669] #[pixel]
profilexspread = [16, 16] #mean+-spread [pixel]
period = 3241:3270

IRdata1 = irev.loadfromNetCDF4(joinpath(pathtofile, ncfile), "irdata")
IRdata2 = IRdata1[:,:,period]
#cut
IRdata2[:,460:538,:] .= NaN
IRdata2 = IRdata2[26:end,23:1000,:]
(IRdata, snowsurf) = irev.setsnowsurfacetonan(IRdata2, 0.0)

irmedian = fill(NaN, size(IRdata, 1), size(IRdata, 2))
irvar = fill(NaN, size(IRdata, 1), size(IRdata, 2))
for i in 1:size(IRdata, 2)
    for j in 1:size(IRdata, 1)
        irmedian[j, i] = median(IRdata[j, i, :])
        irvar[j, i] = std(IRdata[j, i, :])
    end
end

#preallocate arrays
profilesraw = fill(NaN, size(IRdata, 1), maximum(profilexspread)*2+1, length(period), length(profilexmean))
medianprofile = fill(NaN, size(IRdata, 1), length(profilexmean))
medianquantile1 = similar(medianprofile)
medianquantile3 = similar(medianprofile)
lastnumber = zeros(Int64, length(profilexmean))

for j in 1:size(profilexmean, 1)
    profilesraw[:, :, :, j] = IRdata[:, (-profilexspread[j]:profilexspread[j]).+profilexmean[j], :]
    for idx in 1:size(profilesraw, 1)
        try
            (medianquantile1[idx, j], medianprofile[idx, j], medianquantile3[idx, j]) = quantile(vec(filter(!isnan, profilesraw[idx, :, :, j])), [0.25, 0.5, 0.75])
        catch e
        end
    end
    lastnumber[j] = findlast(x->!isnan(x), medianprofile[:, j])
end

#fntsze = 18
#lblsze = fntsze - 5
xis0pxl = 9
yis0pxl = 360
mini = 4
maxi = 8
ctab10 = PyPlot.cm.tab10
colors = [ctab10(1), ctab10(0)]
profilexpos = (profilexmean .- xis0pxl) .* cmperpxlu/100

fig = PyPlot.figure(figsize=(8.5,6.3))
fig.text(0.06, 0.85, "a)", size="large", )
gs = gridspec.GridSpec(2, 3, height_ratios=[5, 5], width_ratios=[6,24,16])
ax = fig.add_subplot(py"$(gs)[0,0:3]")
#ax.set_title(title)
xleft = -(xis0pxl - 1) * cmperpxlu / 100
xright = size(irmedian, 2) * cmperpxlu / 100 + xleft
ybottom = (yis0pxl - size(irmedian, 1)) * cmperpxlw / 100
ytop = size(irmedian, 1) * cmperpxlw / 100 + ybottom
hm = ax.imshow(irmedian, extent=[xleft, xright, ybottom, ytop], cmap="viridis", vmin=mini, vmax=maxi)
PyPlot.axvline(profilexpos[1], 0, 1, color=colors[1])
PyPlot.axvline(profilexpos[2], 0, 1, color=colors[2])
ax.set_xlabel(L"x~\mathrm{[m]}")#, fontsize=fntsze)
ax.set_ylabel(L"h~\mathrm{[cm]}")#, fontsize=fntsze)
IRwindarrow = ax.annotate("", xy=(0.7,2.1), xytext=(0,2.1), annotation_clip=false, arrowprops=Dict("lw" => 2, "arrowstyle" => "->", "facecolor" => "black"))
IRwindarrowtxt = ax.text(0.15, 2.16, "wind")
#ax.tick_params(labelsize=lblsze)
cbar = fig.colorbar(hm, ax=ax)
cbar.set_label(L"T~\mathrm{[^\circ C]}")#, fontsize=fntsze)
#cbar.ax.tick_params(labelsize=lblsze)
fig.text(0.06, 0.47, "b)", size="large", )
medfigax = fig.add_subplot(gs[2,2])
#medfigax.set_title("1")
labels = ["x = 0.5m", "x = 4m"]
for k in 1:size(medianprofile, 2)
    medfigax.plot(medianprofile[1:lastnumber[k],k], (lastnumber[k]:-1:1).*cmperpxlw/100, color=colors[k], label=labels[k])
    medfigax.fill_betweenx((lastnumber[k]:-1:1).*cmperpxlw/100, medianquantile1[1:lastnumber[k], k], medianquantile3[1:lastnumber[k], k], alpha=0.2, color=colors[k])
end
#hor1, = medfigax.plot(medianprofile[:,1], (size(medianprofile, 1):-1:1).*cmperpxlw, color=ctab10(0))
#horfill1 = medfigax.fill_betweenx((size(medianprofile, 1):-1:1).*cmperpxlw, medianquantile1[:, 1], medianquantile3[:, 1], alpha=0.2, color=ctab10(0))
#medfigax.set_title("Profile for 31.05.2021 @ 1m south of bare-snow")
medfigax.set_xlim(4.5,9.5)
medfigax.set_xlabel(L"T~\mathrm{[^\circ C]}")
medfigax.set_ylabel(L"\tilde{h}~\mathrm{[cm]}")
medfigax.grid(true)
medfigax.set_xlim(4,9.5)
medfigax.set_ylim(0,1.80)
medfigax.legend()

######################################################
###     FIG SUP1: MRD FOR REYNOLDS-AVG TIME       ###
######################################################

using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

readfilename = Array{String}(undef, 5)
readfilename[1] = joinpath(datapath, "2dmrd", "db_tot", "tjk_wT_norm.nc")
readfilename[2] = joinpath(datapath, "2dmrd", "db_tot", "t2ucsat_wT_norm.nc")
readfilename[3] = joinpath(datapath, "2dmrd", "db_tot", "t2lcsat_wT_norm.nc")
readfilename[4] = joinpath(datapath, "2dmrd", "db_tot", "t2irg_wT_norm.nc")
readfilename[5] = joinpath(datapath, "2dmrd", "db_tot", "t1irg_wT_norm.nc")

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Millisecond(2^10*timestep),
    "KAIJO" => Second(400))

occurpos1 = 0.0
occurpos2 = 0.0
occurpos3 = 0.0
occurpos4 = 0.0
occurpos5 = 0.0
occurneg1 = 0.0
occurneg2 = 0.0
occurneg3 = 0.0
occurneg4 = 0.0
occurneg5 = 0.0

fig = PyPlot.figure(figsize=(10,8.5))
gs = gridspec.GridSpec(3,2, hspace=0.3)
for i in 1:length(readfilename)
    (time_middle, x, mrd_data, meanwind, evalstart, evalend,
        cols_for_mrd, nrpoints, blocklen, shift, sourcefile) =
        MRD.read2dmrdfromnetcdf(readfilename[i])

    figtitle = MRD.create2dmrdtitle(sourcefile, evalstart, evalend)

    fxavgperiod = blocklen * (evalend - evalstart)

    #load TJK data for wind direction info
    tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
    tjkdata = tjkmeteo[evalstart.<=tjkmeteo[:, 1].<=evalend, :]

    #create winddirection classes
    windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
    windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
    windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

    #load data and calc fluxes
    inst = MRD.instrumentnamefromsourcefile(sourcefile)
    fluxfile = inst2file[inst]
    ra = inst2ra[inst]

    #take subinterval from 2D-MRD
    timemiddle_start = DateTime(2021, 04, 01, 00, 00, 00)
    timemiddle_end = DateTime(2021, 07, 01, 00, 00, 00)

    mrd_data_sub = mrd_data[:, timemiddle_start.<=time_middle.<=timemiddle_end]
    time_middle_sub = time_middle[timemiddle_start.<=time_middle.<=timemiddle_end]
    x_sub = x[:, timemiddle_start.<=time_middle.<=timemiddle_end]

    #categorize wT-decomposition in stability classes
    mrd_data_todo = mrd_data_sub

    stabclass = zeros(Float64, size(mrd_data_todo, 2))
    sumupto = Second(50)
    sumuptoidx = findfirst(x -> x > sumupto, x[:, 1])
    for i in 1:size(mrd_data_todo, 2)
        stabclass[i] = sum(mrd_data_todo[1:sumuptoidx, i])
    end

    if i == 1
        occurpos1 = count(x -> x .>= 0, stabclass)/size(mrd_data_todo, 2)
        occurneg1 = count(x -> x .< 0, stabclass)/size(mrd_data_todo, 2)
    elseif i == 2
        occurpos2 = count(x -> x .>= 0, stabclass)/size(mrd_data_todo, 2)
        occurneg2 = count(x -> x .< 0, stabclass)/size(mrd_data_todo, 2)
    elseif i == 3
        occurpos3 = count(x -> x .>= 0, stabclass)/size(mrd_data_todo, 2)
        occurneg3 = count(x -> x .< 0, stabclass)/size(mrd_data_todo, 2)
    elseif i == 4
        occurpos4 = count(x -> x .>= 0, stabclass)/size(mrd_data_todo, 2)
        occurneg4 = count(x -> x .< 0, stabclass)/size(mrd_data_todo, 2)
    elseif i == 5
        occurpos5 = count(x -> x .>= 0, stabclass)/size(mrd_data_todo, 2)
        occurneg5 = count(x -> x .< 0, stabclass)/size(mrd_data_todo, 2)
    end

    posmrd_data = mrd_data_sub[:, stabclass.>=0]
    negmrd_data = mrd_data_sub[:, stabclass.<0]

    posmrd_med = zeros(Float64, size(posmrd_data, 1))
    posmrd_q1 = similar(posmrd_med)
    posmrd_q3 = similar(posmrd_med)
    negmrd_med = zeros(Float64, size(negmrd_data, 1))
    negmrd_q1 = similar(negmrd_med)
    negmrd_q3 = similar(negmrd_med)

    for i in 1:size(posmrd_data, 1)
        (posmrd_q1[i], posmrd_med[i], posmrd_q3[i]) = quantile(filter(!isnan, posmrd_data[i, :]), [0.25, 0.5, 0.75])
        (negmrd_q1[i], negmrd_med[i], negmrd_q3[i]) = quantile(filter(!isnan, negmrd_data[i, :]), [0.25, 0.5, 0.75])
    end

    title = "xxx"
    if inst == "TJK"
        title = "5 m (meteo station)"
    elseif inst == "T2LCSAT"
        title = "2 m (T2)"
    elseif inst == "T2UCSAT"
        title = "3 m (T2)"
    elseif inst == "T2IRG"
        title = "1 m (T2)"
    elseif inst == "T1IRG"
        title = "1 m (T1)"
    end

    ylow = -5
    yhigh = 12

    if i == 1
        ax = fig.add_subplot(gs[1,1])
        ax.set_title(title)
        ax.plot(Dates.value.(x[:, 1]) ./ 1000, posmrd_med .* 1e3, color=cramericm.vik(0.9), label="unstable+neutral")
        ax.fill_between(Dates.value.(x[:, 1]) ./ 1000, posmrd_q1 .* 1e3, posmrd_q3 .* 1e3, color=cramericm.vik(0.9), alpha=0.5)
        ax.plot(Dates.value.(x[:, 1]) ./ 1000, negmrd_med .* 1e3, color=cramericm.vik(0.1), label="stable")
        ax.fill_between(Dates.value.(x[:, 1]) ./ 1000, negmrd_q1 .* 1e3, negmrd_q3 .* 1e3, color=cramericm.vik(0.1), alpha=0.5)
        PyPlot.xscale("log")
        ax.grid()
        #ax.set_xlabel(L"$\tau~\mathrm{[s]}$")
        ax.set_ylabel(L"C_{wT_S} [10^{-3}~\mathrm{K~m~s^{-1}}]")
        ax.set_ylim(ylow, yhigh)
        ax.tick_params(axis="x", labelbottom=false)
    elseif i==2
        ax2 = fig.add_subplot(gs[2,1])
        ax2.set_title(title)
        ax2.plot(Dates.value.(x[:, 1]) ./ 1000, posmrd_med .* 1e3, color=cramericm.vik(0.9), label="unstable+neutral")
        ax2.fill_between(Dates.value.(x[:, 1]) ./ 1000, posmrd_q1 .* 1e3, posmrd_q3 .* 1e3, color=cramericm.vik(0.9), alpha=0.5)
        ax2.plot(Dates.value.(x[:, 1]) ./ 1000, negmrd_med .* 1e3, color=cramericm.vik(0.1), label="stable")
        ax2.fill_between(Dates.value.(x[:, 1]) ./ 1000, negmrd_q1 .* 1e3, negmrd_q3 .* 1e3, color=cramericm.vik(0.1), alpha=0.5)
        PyPlot.xscale("log")
        ax2.grid()
        ax2.set_ylabel(L"C_{wT_S} [10^{-3}~\mathrm{K~m~s^{-1}}]")
        #ax2.set_xlabel(L"$\tau~\mathrm{[s]}$")
        ax2.set_ylim(ylow, yhigh)
        ax2.tick_params(axis="x", labelbottom=false)
    elseif i==3
        ax3 = fig.add_subplot(gs[2,2])
        ax3.set_title(title)
        ax3.plot(Dates.value.(x[:, 1]) ./ 1000, posmrd_med .* 1e3, color=cramericm.vik(0.9), label="unstable+neutral")
        ax3.fill_between(Dates.value.(x[:, 1]) ./ 1000, posmrd_q1 .* 1e3, posmrd_q3 .* 1e3, color=cramericm.vik(0.9), alpha=0.5)
        ax3.plot(Dates.value.(x[:, 1]) ./ 1000, negmrd_med .* 1e3, color=cramericm.vik(0.1), label="stable")
        ax3.fill_between(Dates.value.(x[:, 1]) ./ 1000, negmrd_q1 .* 1e3, negmrd_q3 .* 1e3, color=cramericm.vik(0.1), alpha=0.5)
        PyPlot.xscale("log")
        ax3.grid()
        #ax3.set_xlabel(L"$\tau~\mathrm{[s]}$")
        ax3.set_ylim(ylow, yhigh)
        ax3.tick_params(axis="x", labelbottom=false)
        ax3.tick_params(axis="y", labelleft=false)
    elseif i==4
        ax4 = fig.add_subplot(gs[3,1])
        ax4.set_title(title)
        ax4.plot(Dates.value.(x[:, 1]) ./ 1000, posmrd_med .* 1e3, color=cramericm.vik(0.9), label="unstable+neutral")
        ax4.fill_between(Dates.value.(x[:, 1]) ./ 1000, posmrd_q1 .* 1e3, posmrd_q3 .* 1e3, color=cramericm.vik(0.9), alpha=0.5)
        ax4.plot(Dates.value.(x[:, 1]) ./ 1000, negmrd_med .* 1e3, color=cramericm.vik(0.1), label="stable")
        ax4.fill_between(Dates.value.(x[:, 1]) ./ 1000, negmrd_q1 .* 1e3, negmrd_q3 .* 1e3, color=cramericm.vik(0.1), alpha=0.5)
        PyPlot.xscale("log")
        ax4.grid()
        ax4.set_ylabel(L"C_{wT_S} [10^{-3}~\mathrm{K~m~s^{-1}}]")
        ax4.set_xlabel(L"$\tau~\mathrm{[s]}$")
        ax4.set_ylim(ylow, yhigh)
    elseif i==5
        ax5 = fig.add_subplot(gs[3,2])
        ax5.set_title(title)
        ax5.plot(Dates.value.(x[:, 1]) ./ 1000, posmrd_med .* 1e3, color=cramericm.vik(0.9), label="unstable+neutral")
        ax5.fill_between(Dates.value.(x[:, 1]) ./ 1000, posmrd_q1 .* 1e3, posmrd_q3 .* 1e3, color=cramericm.vik(0.9), alpha=0.5)
        ax5.plot(Dates.value.(x[:, 1]) ./ 1000, negmrd_med .* 1e3, color=cramericm.vik(0.1), label="stable")
        ax5.fill_between(Dates.value.(x[:, 1]) ./ 1000, negmrd_q1 .* 1e3, negmrd_q3 .* 1e3, color=cramericm.vik(0.1), alpha=0.5)
        PyPlot.xscale("log")
        ax5.grid()
        ax5.set_xlabel(L"$\tau~\mathrm{[s]}$")
        ax5.set_ylim(ylow, yhigh)
        ax5.tick_params(axis="y", labelleft=false)
    end
end

@show occurneg1
@show occurpos1
@show occurneg2
@show occurpos2
@show occurneg3
@show occurpos3
@show occurneg4
@show occurpos4
@show occurneg5
@show occurpos5

######################################################
###      FIG SUP2: SPECTRAL DECOMPOSITION         ###
######################################################

using Dates, PyCall, DataFrames, Statistics, LaTeXStrings
using ProgressMeter, Distributed, LsqFit, FFTW
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
animation = pyimport("matplotlib.animation")
gridspec = pyimport("matplotlib.gridspec")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
    datapath = "/home/michi/Documents/slf/data/"
    tjkpath = "/home/michi/Documents/slf/data/tjk/tjk_data.csv"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
    datapath = "/home/haugened/Documents/data/"
    tjkpath = "/home/haugened/Documents/data/tjk/tjk_data.csv"
end

include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "fourier.jl"))
import .ft
import .turb
import .gen
PyPlot.pygui(true)

#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
kaijo_period_file = joinpath(datapath, "kaijo", "21_kaijo_periods.txt")
kaijo_outfile_stam = joinpath(datapath, "kaijo", "kaijo")
outfile_stam = joinpath(datapath, "WindSonic_Duerrboden_2005xx/tjk_sonic_200504-200514")
tower_outfile_stam = joinpath(datapath, "tower", "preproc")
ventair_stam = joinpath(datapath, "tower", "vent_air")

#select data and measurement period to be evaluated
evalstart = DateTime(2021, 05, 20, 23, 00, 00)
evalend   = DateTime(2021, 06, 10, 23, 00, 00)
#evalend = evalstart + Day(10)

evaldf1 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t1irgdb.nc"), evalstart, evalend)
evaldf2 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2irgdb.nc"), evalstart, evalend)
evaldf3 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2lcsatdb.nc"), evalstart, evalend)
evaldf4 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2ucsatdb.nc"), evalstart, evalend)
#evaldf5 = turb.readturbasnetcdf(string(kaijo_outfile_stam, ".nc"), evalstart, evalend)
evaldf6 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "tjkdf.nc"), evalstart, evalend)

#apply NaN-mask to T1 & T2 when repositioned (for DR)
turb.repositionnanmask!(evaldf1)
turb.repositionnanmask!(evaldf2)
turb.repositionnanmask!(evaldf3)
turb.repositionnanmask!(evaldf4)

#show statistics about missing data
turb.printmissstats(evaldf1)
turb.printmissstats(evaldf2)
turb.printmissstats(evaldf3)
turb.printmissstats(evaldf4)
#turb.printmissstats(evaldf5)
turb.printmissstats(evaldf6)
#interpolate missing
evaldf1 = turb.interpolatemissing(evaldf1)
evaldf2 = turb.interpolatemissing(evaldf2)
evaldf3 = turb.interpolatemissing(evaldf3)
evaldf4 = turb.interpolatemissing(evaldf4)
#=try
    evaldf5 = turb.interpolatemissing(evaldf5)
catch LoadError
    @warn("Kaijo DataFrame empty")
end=#
evaldf6 = turb.interpolatemissing(evaldf6)
#double rotation
turb.drdf!(evaldf1)
turb.drdf!(evaldf2)
turb.drdf!(evaldf3)
turb.drdf!(evaldf4)
#turb.drdf!(evaldf5)
turb.drdf!(evaldf6, periodwise=false)

#data to decompose
turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)
turb.missing2nan!(evaldf6)

disallowmissing!(evaldf1)
disallowmissing!(evaldf2)
disallowmissing!(evaldf3)
disallowmissing!(evaldf4)
disallowmissing!(evaldf6)

ftdata1 = evaldf1.T
ftdata2 = evaldf2.T
ftdata3 = evaldf3.T
ftdata4 = evaldf4.T
ftdata6 = evaldf6.T

#remove NaN-values
ftdata1 = ftdata1[.!isnan.(ftdata1)]
ftdata2 = ftdata2[.!isnan.(ftdata2)]
ftdata3 = ftdata3[.!isnan.(ftdata3)]
ftdata4 = ftdata4[.!isnan.(ftdata4)]
ftdata6 = ftdata6[.!isnan.(ftdata6)]

#duration in seconds
dur1 = length(ftdata1)*0.05
dur2 = length(ftdata2)*0.05
dur3 = length(ftdata3)*0.05
dur4 = length(ftdata4)*0.05
dur6 = length(ftdata6)*0.05

#skip for long datasets
#=
#detrend the vector (subtract linear fit)
vec1 = ft.detrend(ftdata1)
vec2 = ft.detrend(ftdata2)
vec3 = ft.detrend(ftdata3)
vec4 = ft.detrend(ftdata4)
vec6 = ft.detrend(ftdata6)
#apply bell taper
vect1 = ft.belltaper(vec1)
vect2 = ft.belltaper(vec2)
vect3 = ft.belltaper(vec3)
vect4 = ft.belltaper(vec4)
vect6 = ft.belltaper(vec6)
=#

println("Applying Fourier-Trafo...")
(FTraw1, FTvec1) = ft.fourierplus(ftdata1)
(FTraw2, FTvec2) = ft.fourierplus(ftdata2)
(FTraw3, FTvec3) = ft.fourierplus(ftdata3)
(FTraw4, FTvec4) = ft.fourierplus(ftdata4)
(FTraw6, FTvec6) = ft.fourierplus(ftdata6)
(Soff1, freq1) = ft.spectralenergydensity(FTvec1, dur1)
(Soff2, freq2) = ft.spectralenergydensity(FTvec2, dur2)
(Soff3, freq3) = ft.spectralenergydensity(FTvec3, dur3)
(Soff4, freq4) = ft.spectralenergydensity(FTvec4, dur4)
(Soff6, freq6) = ft.spectralenergydensity(FTvec6, dur6)

freqtimesSoff1 = freq1.*Soff1
freqtimesSoff2 = freq2.*Soff2
freqtimesSoff3 = freq3.*Soff3
freqtimesSoff4 = freq4.*Soff4
freqtimesSoff6 = freq6.*Soff6
(freq1, freqtimesSoff1) = ft.logavg(freq1, freqtimesSoff1, 0.05)
(freq2, freqtimesSoff2) = ft.logavg(freq2, freqtimesSoff2, 0.05)
(freq3, freqtimesSoff3) = ft.logavg(freq3, freqtimesSoff3, 0.05)
(freq4, freqtimesSoff4) = ft.logavg(freq4, freqtimesSoff4, 0.05)
(freq6, freqtimesSoff6) = ft.logavg(freq6, freqtimesSoff6, 0.05)

#theoretical curves
x1 = collect(2:0.05:8)
x2 = collect(2:0.05:8)
x3 = collect(2:0.05:8)
x4 = collect(2:0.05:8)
x6 = collect(2:0.05:8)
y1=exp.(-2*log.(x1)/3)./35
y2=exp.(-2*log.(x2)/3)./30
y3=exp.(-2*log.(x3)/3)./40
y4=exp.(-2*log.(x4)/3)./55
y6=exp.(-2*log.(x6)/3)./110

xlimit = (4e-4, 2e1)
ylimit = (1e-3, 1e-1)

fig2 = PyPlot.figure(figsize=(10,8.5))
gs = gridspec.GridSpec(3,2, hspace=0.3)
ax1 = fig2.add_subplot(gs[1,1])
ax1.set_title("5 m (meteostation)")
ax1.set_ylabel("fS(f)")
tmp = gen.movingaverage(freqtimesSoff6, 8)
ax1.set_xlim(xlimit)
ax1.set_ylim(ylimit)
ax1.plot(freq6, tmp)
ax1.plot(x6,y6, color="black")
PyPlot.xscale("log")
PyPlot.yscale("log")
ax1.tick_params(axis="x", labelbottom=false)
ax1.grid()

ax2 = fig2.add_subplot(gs[2,1], sharex=ax1)
ax2.set_title("3 m (T2)")
ax2.set_ylabel("fS(f)")
ax2.set_ylim(ylimit)
tmp = gen.movingaverage(freqtimesSoff4, 8)
ax2.plot(freq4, tmp)
ax2.plot(x4,y4, color="black")
PyPlot.xscale("log")
PyPlot.yscale("log")
ax2.tick_params(axis="x", labelbottom=false)
ax2.grid()

ax3 = fig2.add_subplot(gs[2,2], sharey=ax2)
ax3.set_title("2 m (T2)")
ax3.set_xlim(xlimit)
tmp = gen.movingaverage(freqtimesSoff3, 8)
ax3.plot(freq3, tmp)
ax3.plot(x3,y3, color="black")
PyPlot.xscale("log")
PyPlot.yscale("log")
ax3.tick_params(axis="x", labelbottom=false)
ax3.tick_params(axis="y", labelleft=false)
ax3.grid()

ax4 = fig2.add_subplot(gs[3,1], sharex=ax2)
ax4.set_title("1 m (T2)")
ax4.set_xlabel("f [Hz]")
ax4.set_ylabel("fS(f)")
tmp = gen.movingaverage(freqtimesSoff2, 8)
ax4.set_ylim(ylimit)
ax4.plot(freq2, tmp)
ax4.plot(x2,y2, color="black")
PyPlot.xscale("log")
PyPlot.yscale("log")
ax4.grid()

ax6 = fig2.add_subplot(gs[3,2], sharex=ax3, sharey=ax4)
ax6.set_title("1 m (T1)")
ax6.set_xlabel("f [Hz]")
tmp = gen.movingaverage(freqtimesSoff1, 8)
ax6.plot(freq1, tmp)
ax6.plot(x1,y1, color="black")
PyPlot.xscale("log")
PyPlot.yscale("log")
ax6.tick_params(axis="y", labelleft=false)
ax6.grid()
##

######################################################
###          VIDEO SUP3: SCREEN VIDEO              ###
######################################################
#ESM_3.mp4 supplementary: video

##
using PyCall, Statistics, LaTeXStrings, NCDatasets
import PyPlot
animation = pyimport("matplotlib.animation")
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

#functions
function loadfromNetCDF4(file::String, varname::String)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = dat[:,:,:]
    close(netfile)
    return dattot
end


#definitions
irfile = "/home/haugened/Documents/data/ir/210531_143457/210531_143457_irdata09to11_for_pub.nc"
cmperpxl = 0.6
cmperpxlu = 0.606
cmperpxlw = 0.549
framespersec = 29.9

#=
IRdata2 = copy(IRdata)
IRdata = IRdata[:,:,1800:2815]
IRdata[:,460:538,:] .= NaN
IRdata = IRdata[26:end, 23:1000,:]
(IRdata, snowsurf) = irev.setsnowsurfacetonan(IRdata, 0.0, framewise=true)

IRdata = copy(IRdata2)
=#

IRdata = loadfromNetCDF4(irfile, "irdata")

#video for publication supplementary
"Update-function for the animation"
function pubanimupdate(frame::Int64)
    #pubfig.suptitle("")
    hm.set_data(IRdata[:,:,frame])
    return hm
end

fntsze = 16
lblsze = fntsze-2
xis0pxl = 9
PyPlot.pygui(true)
pubfig = PyPlot.figure(figsize=(16,5.2))
pubax = pubfig.add_subplot(1,1,1)
PyPlot.suptitle("Coupling of Submeso-Scale Motion with Near-Surface Turbulence over Patchy Snow Cover", fontsize=fntsze)
PyPlot.title("Haugeneder et al. (2023)", fontsize=lblsze)
xleft = -(xis0pxl-1)*cmperpxlu/100
xright = size(IRdata, 2)*cmperpxlu/100+xleft
ybottom = 0
ytop = size(IRdata, 1) * cmperpxlw/100
hm = pubax.imshow(IRdata[:,:,1], extent = [xleft,xright,ybottom,ytop],cmap="turbo", vmin=4, vmax=8)
pubax.set_xlabel(L"x~\mathrm{[m]}", fontsize=fntsze)
pubax.set_ylabel(L"h~\mathrm{[m]}", fontsize=fntsze)
pubax.tick_params(labelsize=lblsze)
divider = mpl_axes_grid1.make_axes_locatable(pubax)
cpubax = divider.append_axes("right", size="3%", pad=0.1)
cbar = pubfig.colorbar(hm, cax=cpubax)
cbar.set_label(L"T~\mathrm{[^\circ C]}", fontsize=fntsze)
cbar.ax.tick_params(labelsize=lblsze)
PyPlot.show()
pubani = animation.FuncAnimation(pubfig, pubanimupdate, frames=collect(1:1800), interval = 5, repeat_delay=500)
pubani.save("/home/haugened/Documents/process_paper_23/figures/ESM_3_tmp.mp4", fps=framespersec)
##

######################################################
###  FIG SUP4: PROFILE DECOMPOSITION WT WITH KAIJO ###
######################################################
using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings, DataFrames
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
widgets = pyimport("matplotlib.widgets")
cramericm = pyimport("cmcrameri.cm")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    datapath = "/home/michi/Documents/slf/data/"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
elseif gethostname() == "LINUX24"
    datapath = "/home/haugened/Documents/data/"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "mrd.jl"))
include(joinpath(importdir, "src", "general.jl"))
import .turb
import .MRD
import .gen

PyPlot.pygui(true)
timestep = Millisecond(50)

tmp = joinpath(datapath, "tower", "preproc")
tmp_kaijo = joinpath(datapath, "kaijo", "kaijo.nc")
inst2file = Dict("T1IRG" => joinpath(tmp, "t1irgdb.nc"),
    "T2IRG" => joinpath(tmp, "t2irgdb.nc"),
    "T2LCSAT" => joinpath(tmp, "t2lcsatdb.nc"),
    "T2UCSAT" => joinpath(tmp, "t2ucsatdb.nc"),
    "TJK" => joinpath(tmp, "tjkdf.nc"),
    "KAIJO" => tmp_kaijo)

inst2ra = Dict("T1IRG" => Millisecond(2^11*timestep),
    "T2IRG" => Millisecond(2^11*timestep),
    "T2LCSAT" => Millisecond(2^11*timestep),
    "T2UCSAT" => Millisecond(2^11*timestep),
    "TJK" => Millisecond(2^10*timestep),
    "KAIJO" => Millisecond(2^9*timestep))

fluxtype = "wT"

panelfile1 = joinpath(datapath, "2dmrd", "db_tot", string("t2ucsat_", fluxtype, "_norm.nc"))
panelfile2 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxtype, "_norm.nc"))
panelfile3 = joinpath(datapath, "2dmrd", "db_tot", string("t2irg_", fluxtype, "_norm.nc"))
panelfile4 = joinpath(datapath, "2dmrd", "db_tot", string("kaijo_", fluxtype, "_norm.nc"))

(panel1_time_middle, panel1_x, panel1_mrd_data, panel1_meanwind,
    panel1_evalstart, panel1_evalend, panel1_cols_for_mrd,
    panel1_nrpoints, panel1_blocklen, panel1_shift, panel1_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile1)

(panel2_time_middle, panel2_x, panel2_mrd_data, panel2_meanwind,
    panel2_evalstart, panel2_evalend, panel2_cols_for_mrd,
    panel2_nrpoints, panel2_blocklen, panel2_shift, panel2_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile2)

(panel3_time_middle, panel3_x, panel3_mrd_data, panel3_meanwind,
    panel3_evalstart, panel3_evalend, panel3_cols_for_mrd,
    panel3_nrpoints, panel3_blocklen, panel3_shift, panel3_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile3)

(panel4_time_middle, panel4_x, panel4_mrd_data, panel4_meanwind,
    panel4_evalstart, panel4_evalend, panel4_cols_for_mrd,
    panel4_nrpoints, panel4_blocklen, panel4_shift, panel4_sourcefile) =
    MRD.read2dmrdfromnetcdf(panelfile4)

panelinst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
panelinst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
panelinst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
panelinst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

panel1ra = inst2ra[panelinst1]
panel2ra = inst2ra[panelinst2]
panel3ra = inst2ra[panelinst3]
panel4ra = inst2ra[panelinst4]

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[panel1_evalstart.<=tjkmeteo[:, 1].<=panel2_evalend, :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 2
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 7

#load data and calc fluxes
inst1 = MRD.instrumentnamefromsourcefile(panel1_sourcefile)
inst2 = MRD.instrumentnamefromsourcefile(panel2_sourcefile)
inst3 = MRD.instrumentnamefromsourcefile(panel3_sourcefile)
inst4 = MRD.instrumentnamefromsourcefile(panel4_sourcefile)

fluxfile1 = inst2file[inst1]
ra1 = inst2ra[inst1]

fx1df = turb.readturbasnetcdf(fluxfile1, panel1_evalstart, panel1_evalend)
#show statistics about missing data
turb.printmissstats(fx1df)
#interpolate missing
fx1df = turb.interpolatemissing(fx1df)
#double rotation
if inst1 == "TJK"
    turb.drdf!(fx1, periodwise=false)
else
    turb.drdf!(fx1df)
end
turb.missing2nan!(fx1df)
wndspd = gen.movingaverage(sqrt.(fx1df.u .^ 2 .+ fx1df.v .^ 2 .+ fx1df.w .^ 2), 20 * 3600)

fluxpath = "/home/haugened/Documents/data/fluxes/"
fx1file = joinpath(fluxpath, string(lowercase(inst1), "_fx.nc"))
fx2file = joinpath(fluxpath, string(lowercase(inst2), "_fx.nc"))
fx3file = joinpath(fluxpath, string(lowercase(inst3), "_fx.nc"))
fx4file = joinpath(fluxpath, string(lowercase(inst4), "_fx.nc"))

fx1 = turb.readturbasnetcdf(fx1file, panel1_evalstart, panel1_evalend)
fx2 = turb.readturbasnetcdf(fx2file, panel2_evalstart, panel2_evalend)
fx3 = turb.readturbasnetcdf(fx3file, panel3_evalstart, panel3_evalend)
fx4 = turb.readturbasnetcdf(fx4file, panel4_evalstart, panel4_evalend)

if fluxtype == "wT"
    cmin = -5 #colorbar limits
    cmax = -cmin
elseif fluxtype == "tke"
    cmin = 0
    cmax = 10
else
    @error("fluxtype not correct")
end

##
#plot
fig = PyPlot.figure(figsize=(11,8))
gs = gridspec.GridSpec(6, 2, height_ratios=[5, 15, 15, 15, 15, 8], width_ratios=[50,1])
date_format = pydates.DateFormatter("%d.%m. %H:%M")
majorloc = pydates.HourLocator([0,6,12,18])
minorloc = pydates.HourLocator(interval=1)
ax = fig.add_subplot(gs[1, 1])
ax.plot(fx1df.time[1:20*600:end], wndspd[1:20*600:end], color="black")
#ax4.plot(tjkdf.time, tjkdf.wind_max, color="black", lw=1, alpha=0.6)
yl = ax.get_ylim()
ax.pcolormesh(pydates.date2num.(tjkdata.time), collect(range(first(yl), last(yl), 2)), repeat(transpose(windclass), 2), shading="auto", cmap="tab10", alpha=0.5, vmin=1, vmax=9)
ax.set_ylabel(L"| \mathbf{u}|_{3\mathrm{m}}~\mathrm{[m~s^{-1}]}")
ax.tick_params(axis="x", labelbottom=false)
ax.grid()
ax.set_ylim(0,5.5)
ax.set_xlim(DateTime(2021,06,01,06,00,00), DateTime(2021,06,02,12,00,00))
fig.text(0.825, 0.73, "3m", rotation=-90, fontweight="bold")
ax1 = fig.add_subplot(gs[2, 1], sharex=ax)
#ax1.set_title("5m meteo station")
if fluxtype == "wT"
    mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot = ax1.pcolormesh(panel1_time_middle, Dates.value.(panel1_x) ./ 1000, panel1_mrd_data .* 100, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
end
PyPlot.yscale("log")
ax1.xaxis_date()
ax1.set_ylim([0.1, 2e3])
ax1.tick_params(axis="x", labelbottom=false)
ax1.axhline(Dates.value(Millisecond(panel1ra))./1000, ls="--", color="black", alpha=0.3)
ax1.fill_betweenx((Dates.value(Millisecond(panel1ra))./1000, last(ax1.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.57, "2m", rotation=-90, fontweight="bold")
ax2 = fig.add_subplot(gs[3, 1], sharex=ax, sharey=ax1)
#ax2.set_title("3m T2")
if fluxtype == "wT"
    mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot2 = ax2.pcolormesh(panel2_time_middle, Dates.value.(panel2_x) ./ 1000, panel2_mrd_data .* 100, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
end
ax2.tick_params(axis="x", labelbottom=false)
#ax2.tick_params(axis="y", labelleft=false)
ax2.axhline(Dates.value(Millisecond(panel2ra))./1000, ls="--", color="black", alpha=0.3)
ax2.fill_betweenx((Dates.value(Millisecond(panel2ra))./1000, last(ax2.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
fig.text(0.825, 0.42, "1m", rotation=-90, fontweight="bold")
ax3 = fig.add_subplot(gs[4, 1], sharex=ax, sharey=ax1)
#ax3.set_title("2m T2")
if fluxtype == "wT"
    mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot3 = ax3.pcolormesh(panel3_time_middle, Dates.value.(panel3_x) ./ 1000, panel3_mrd_data .* 100, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
end
ax3.axhline(Dates.value(Millisecond(panel3ra))./1000, ls="--", color="black", alpha=0.3)
ax3.fill_betweenx((Dates.value(Millisecond(panel3ra))./1000, last(ax3.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax3.tick_params(axis="x", labelbottom=false)
ax4 = fig.add_subplot(gs[5, 1], sharex=ax, sharey=ax1)
#ax4.set_title("1m T2")
if fluxtype == "wT"
    mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* 1000, cmap=cramericm.vik, vmin=cmin, vmax=cmax, shading="auto")
elseif fluxtype == "tke"
    mrdplot4 = ax4.pcolormesh(panel4_time_middle, Dates.value.(panel4_x) ./ 1000, panel4_mrd_data .* 100, cmap="turbo", vmin=cmin, vmax=cmax, shading="auto")
end
#ax4.tick_params(axis="y", labelleft=false)
ax4.axhline(Dates.value(Millisecond(panel4ra))./1000, ls="--", color="black", alpha=0.3)
ax4.fill_betweenx((Dates.value(Millisecond(panel4ra))./1000, last(ax4.get_ylim())), first(ax.get_xlim()), last(ax.get_xlim()), color="white", alpha=0.4)
ax4.tick_params(axis="x", labelbottom=false)
fig.text(0.825, 0.26, "0.3m", rotation=-90, fontweight="bold")
ax5 = fig.add_subplot(py"$(gs)[2:-2, 1]")
cbar = fig.colorbar(mrdplot, cax=ax5, orientation="vertical")
if fluxtype == "wT"
    cbar.set_label(L"$C_{wT_s} [10^{-3}~\mathrm{K~m~s^{-1}}]$")
elseif fluxtype == "tke"
    cbar.set_label(L"$C_{e}~[\mathrm{10^{-2}~m^2~s^{-2}}]$")
end
fig.text(0.06, 0.49, L"$\tau~\mathrm{[s]}$", rotation="vertical")
ax6 = fig.add_subplot(gs[6, 1], sharex = ax)
if fluxtype == "wT"
    ax6.plot(fx1.time, fx1.wT, label="3m")
    ax6.plot(fx4.time, fx4.wT, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$\overline{w'T_s'}~[\mathrm{K~m~s^{-1}}]$")
    ax6.set_ylim(-0.1, 0.1)
elseif fluxtype == "tke"
    ax6.plot(fx1.time, fx1.tke, label="3m")
    ax6.plot(fx4.time, fx4.tke, alpha=0.9, label="0.3m")
    ax6.set_ylabel(L"$e~\mathrm{[m^2~s^{-2}]}$")
    ax6.set_ylim(0,1.0)
end
ax6.legend(loc="center left", bbox_to_anchor=(1,0.5))
ax6.grid()
ax6.xaxis.set_major_locator(majorloc)
ax6.xaxis.set_minor_locator(minorloc)
ax6.xaxis.set_major_formatter(date_format)
sleep(Second(1))
labels1 = ax6.xaxis.get_majorticklabels()#[item.get_text() for item in ax6.get_xticklabels()]
newlabels = Vector{String}(undef, length(labels1))
newlabels[1] = labels1[1].get_text()
for idx in 2:length(labels1)
    tmp = labels1[idx].get_text()
    if tmp[end-4:end-3] != "00" && tmp[3] == '.'
        newlabels[idx] = tmp[end-4:end]
    else
        newlabels[idx] = tmp
    end
end
ax6.xaxis.set_ticklabels(newlabels)