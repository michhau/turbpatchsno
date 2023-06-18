######################################################
###   ADDITIONAL VARIABLES TO CLASSIFY TKE-MRDS    ###
######################################################
#=
Find additional dependencies (additional to wind direction)
for 'linear' MRDs to collaps better (Fig. 7)
(Iva's feedback for process_paper_23)
=#

using Dates, Statistics, NCDatasets, PyCall, LaTeXStrings, DataFrames
using ProgressMeter
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

fluxtype = "tke"

if !(fluxtype in(["tke", "wT", "TT"]))
    @error("fluxtype not recognized")
end

multfactor = 1
if fluxtype == "tke" || fluxtype == "TT"
    multfactor = 100
elseif fluxtype == "wT"
    multfactor = 1000
end

ylabelstring = L"not assigned"
if fluxtype == "tke"
    ylabelstring = L"$C_{e}~[10^{-2}~\mathrm{J~kg^{-1}}]$"
elseif fluxtype == "wT"
    ylabelstring = L"$C_{wT_s}~[10^{-3}~\mathrm{K~m~s^{-1}}]$"
elseif fluxtype == "TT"
    ylabelstring = L"$C_{T_sT_s}~[10^{-2}~\mathrm{K^2}]$"
end

readfilename1 = joinpath(datapath, "2dmrd", "db_tot", string("tjk_", fluxtype, "_norm.nc"))
readfilename2 = joinpath(datapath, "2dmrd", "db_tot", string("t2lcsat_", fluxtype, "_norm.nc"))

figsavefoler = string("/home/haugened/Documents/plots/mrd_tke_classify/days/", fluxtype)

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

figtitle1 = MRD.create2dmrdtitle(sourcefile1, evalstart1, evalend1)
figtitle2 = MRD.create2dmrdtitle(sourcefile2, evalstart2, evalend2)

#load TJK data for wind direction info
tjkmeteo = turb.csvtodataframe(joinpath(datapath, "tjk", "tjk_data.csv"))
tjkdata = tjkmeteo[minimum([evalstart1, evalstart2]).<=tjkmeteo[:, 1].<=maximum([evalend1, evalend2]), :]

#create winddirection classes
windclass = fill(NaN, length(tjkdata.wind_mean_vector_direction))
windclass[310 .<= tjkdata.wind_mean_vector_direction.||tjkdata.wind_mean_vector_direction.<=20] .= 1
windclass[130 .<= tjkdata.wind_mean_vector_direction .<= 200] .= 2

starttimefirstfig = DateTime(2021, 05, 21, 08, 00, 00)
starttimelastfig  = DateTime(2021, 06, 10, 08, 00, 00)
stepbtwfigs = Day(1)

timestostart = collect(starttimefirstfig:stepbtwfigs:starttimelastfig)

PyPlot.pygui(true)
@showprogress "Creating figures..." for i in axes(timestostart, 1)
    timemiddle_start = timestostart[i]
    timemiddle_end = timestostart[i]+Hour(10)

    filename = string(Dates.value(Year(timemiddle_start)),
    lpad(Dates.value(Month(timemiddle_start)), 2, "0"), lpad(Dates.value(Day(timemiddle_start)), 2, "0"))
    filename = joinpath(figsavefoler, filename)

    mrd_data_up1 = copy(mrd_data1[:, timemiddle_start.<=time_middle1.<=timemiddle_end])
    mrd_data_down1 = copy(mrd_data1[:, timemiddle_start.<=time_middle1.<=timemiddle_end])
    time_middle_sub1 = copy(time_middle1[timemiddle_start.<=time_middle1.<=timemiddle_end])
    x_sub1 = copy(x1[:, timemiddle_start.<=time_middle1.<=timemiddle_end])

    mrd_data_up2 = copy(mrd_data2[:, timemiddle_start.<=time_middle2.<=timemiddle_end])
    mrd_data_down2 = copy(mrd_data2[:, timemiddle_start.<=time_middle2.<=timemiddle_end])
    time_middle_sub2 = copy(time_middle2[timemiddle_start.<=time_middle2.<=timemiddle_end])
    x_sub2 = copy(x2[:, timemiddle_start.<=time_middle2.<=timemiddle_end])

    #----------------------------------------------------
    #take only specific wind directions if wanted
    disallowmissing!(tjkdata)
    windclass_2dmrd1 = fill(NaN, size(time_middle_sub1, 1))
    windclass_2dmrd2 = fill(NaN, size(time_middle_sub2, 1))
    turb.extendtofinertimeseries!(windclass_2dmrd1, time_middle_sub1, windclass, tjkdata.time)
    turb.extendtofinertimeseries!(windclass_2dmrd2, time_middle_sub2, windclass, tjkdata.time)

    upwind1 = windclass_2dmrd1 .== 1   #1 = up valley wind; 2 = down valley wind
    downwind1 = windclass_2dmrd1 .== 2   #1 = up valley wind; 2 = down valley wind
    upwind2 = windclass_2dmrd2 .== 1   #1 = up valley wind; 2 = down valley wind
    downwind2 = windclass_2dmrd2 .== 2   #1 = up valley wind; 2 = down valley wind
    mrd_data_up1[:, .!upwind1] .= NaN
    mrd_data_down1[:, .!downwind1] .= NaN
    mrd_data_up2[:, .!upwind2] .= NaN
    mrd_data_down2[:, .!downwind2] .= NaN

    #----------------------------------------------------

    lin_mrd_up_med1 = zeros(Float64, size(mrd_data_up1, 1))
    lin_mrd_up_q11 = zeros(Float64, size(mrd_data_up1, 1))
    lin_mrd_up_q31 = zeros(Float64, size(mrd_data_up1, 1))
    lin_mrd_down_med1 = zeros(Float64, size(mrd_data_down1, 1))
    lin_mrd_down_q11 = zeros(Float64, size(mrd_data_down1, 1))
    lin_mrd_down_q31 = zeros(Float64, size(mrd_data_down1, 1))

    lin_mrd_up_med2 = zeros(Float64, size(mrd_data_up2, 1))
    lin_mrd_up_q12 = zeros(Float64, size(mrd_data_up2, 1))
    lin_mrd_up_q32 = zeros(Float64, size(mrd_data_up2, 1))
    lin_mrd_down_med2 = zeros(Float64, size(mrd_data_down2, 1))
    lin_mrd_down_q12 = zeros(Float64, size(mrd_data_down2, 1))
    lin_mrd_down_q32 = zeros(Float64, size(mrd_data_down2, 1))

    for i in 1:size(mrd_data_up1, 1)
        try
            (lin_mrd_up_q11[i], lin_mrd_up_med1[i], lin_mrd_up_q31[i]) = quantile(filter(!isnan, mrd_data_up1[i, :]), [0.25, 0.5, 0.75])
        catch e
        end
        try
            (lin_mrd_down_q11[i], lin_mrd_down_med1[i], lin_mrd_down_q31[i]) = quantile(filter(!isnan, mrd_data_down1[i, :]), [0.25, 0.5, 0.75])
        catch e    
        end
    end
    for i in 1:size(mrd_data_up2, 1)
        try
           (lin_mrd_up_q12[i], lin_mrd_up_med2[i], lin_mrd_up_q32[i]) = quantile(filter(!isnan, mrd_data_up2[i, :]), [0.25, 0.5, 0.75]) 
        catch e
        end
        try
           (lin_mrd_down_q12[i], lin_mrd_down_med2[i], lin_mrd_down_q32[i]) = quantile(filter(!isnan, mrd_data_down2[i, :]), [0.25, 0.5, 0.75]) 
        catch e
        end
    end

    fig = PyPlot.figure(figsize=(12,4.8))
    gs = gridspec.GridSpec(1,2)
    ax = fig.add_subplot(gs[1,1])
    ax.set_title("2 m (T2)")
    fig.text(0.45, 0.93, Date(timemiddle_start), size="large")
    fig.text(0.12, 0.9, "a)", size="large", )
    ax.plot(Dates.value.(x2[:, 1]) ./ 1000, lin_mrd_down_med2 .* multfactor, color="blue")
    ax.fill_between(Dates.value.(x2[:, 1]) ./ 1000, lin_mrd_down_q12 .* multfactor, lin_mrd_down_q32 .* multfactor, color="blue", alpha=0.5)
    ax.plot(Dates.value.(x2[:, 1]) ./ 1000, lin_mrd_up_med2 .* multfactor, color="red")
    ax.fill_between(Dates.value.(x2[:, 1]) ./ 1000, lin_mrd_up_q12 .* multfactor, lin_mrd_up_q32 .* multfactor, color="red", alpha=0.5)
    PyPlot.xscale("log")
    ax.grid()
    ax.set_xlabel("time scale [s]")
    ax.set_ylabel(ylabelstring)
    ax2 = fig.add_subplot(gs[1,2], sharey=ax)
    ax2.set_title("5 m (meteo station)")
    fig.text(0.544,0.9, "b)", size="large")
    ax2.plot(Dates.value.(x1[:, 1]) ./ 1000, lin_mrd_down_med1 .* multfactor, color="blue")
    ax2.fill_between(Dates.value.(x1[:, 1]) ./ 1000, lin_mrd_down_q11 .* multfactor, lin_mrd_down_q31 .* multfactor, color="blue", alpha=0.5)
    ax2.plot(Dates.value.(x1[:, 1]) ./ 1000, lin_mrd_up_med1 .* multfactor, color="red")
    ax2.fill_between(Dates.value.(x1[:, 1]) ./ 1000, lin_mrd_up_q11 .* multfactor, lin_mrd_up_q31 .* multfactor, color="red", alpha=0.5)
    PyPlot.xscale("log")
    ax2.grid()
    ax2.set_xlabel("time scale [s]")
    if fluxtype !== "wT"
        upylimmax = last(ax2.get_ylim())
        ax.set_ylim(0, minimum([upylimmax, 50]))
    end
    #ax2.tick_params(axis="y", labelleft=false)
    #ax.legend()
    PyPlot.savefig(filename)
    PyPlot.close()
end


############################################
#plot for single period
timemiddle_start = DateTime(2021,06,04,08,00,00)
timemiddle_end = DateTime(2021,06,04,16,00,00)

filename = string(Dates.value(Year(timemiddle_start)),
lpad(Dates.value(Month(timemiddle_start)), 2, "0"), lpad(Dates.value(Day(timemiddle_start)), 2, "0"))
filename = joinpath(figsavefoler, filename)

mrd_data_up1 = copy(mrd_data1[:, timemiddle_start.<=time_middle1.<=timemiddle_end])
mrd_data_down1 = copy(mrd_data1[:, timemiddle_start.<=time_middle1.<=timemiddle_end])
time_middle_sub1 = copy(time_middle1[timemiddle_start.<=time_middle1.<=timemiddle_end])
x_sub1 = copy(x1[:, timemiddle_start.<=time_middle1.<=timemiddle_end])

mrd_data_up2 = copy(mrd_data2[:, timemiddle_start.<=time_middle2.<=timemiddle_end])
mrd_data_down2 = copy(mrd_data2[:, timemiddle_start.<=time_middle2.<=timemiddle_end])
time_middle_sub2 = copy(time_middle2[timemiddle_start.<=time_middle2.<=timemiddle_end])
x_sub2 = copy(x2[:, timemiddle_start.<=time_middle2.<=timemiddle_end])

#----------------------------------------------------
#take only specific wind directions if wanted
disallowmissing!(tjkdata)
windclass_2dmrd1 = fill(NaN, size(time_middle_sub1, 1))
windclass_2dmrd2 = fill(NaN, size(time_middle_sub2, 1))
turb.extendtofinertimeseries!(windclass_2dmrd1, time_middle_sub1, windclass, tjkdata.time)
turb.extendtofinertimeseries!(windclass_2dmrd2, time_middle_sub2, windclass, tjkdata.time)

upwind1 = windclass_2dmrd1 .== 1   #1 = up valley wind; 2 = down valley wind
downwind1 = windclass_2dmrd1 .== 2   #1 = up valley wind; 2 = down valley wind
upwind2 = windclass_2dmrd2 .== 1   #1 = up valley wind; 2 = down valley wind
downwind2 = windclass_2dmrd2 .== 2   #1 = up valley wind; 2 = down valley wind
mrd_data_up1[:, .!upwind1] .= NaN
mrd_data_down1[:, .!downwind1] .= NaN
mrd_data_up2[:, .!upwind2] .= NaN
mrd_data_down2[:, .!downwind2] .= NaN

#----------------------------------------------------

lin_mrd_up_med1 = zeros(Float64, size(mrd_data_up1, 1))
lin_mrd_up_q11 = zeros(Float64, size(mrd_data_up1, 1))
lin_mrd_up_q31 = zeros(Float64, size(mrd_data_up1, 1))
lin_mrd_down_med1 = zeros(Float64, size(mrd_data_down1, 1))
lin_mrd_down_q11 = zeros(Float64, size(mrd_data_down1, 1))
lin_mrd_down_q31 = zeros(Float64, size(mrd_data_down1, 1))

lin_mrd_up_med2 = zeros(Float64, size(mrd_data_up2, 1))
lin_mrd_up_q12 = zeros(Float64, size(mrd_data_up2, 1))
lin_mrd_up_q32 = zeros(Float64, size(mrd_data_up2, 1))
lin_mrd_down_med2 = zeros(Float64, size(mrd_data_down2, 1))
lin_mrd_down_q12 = zeros(Float64, size(mrd_data_down2, 1))
lin_mrd_down_q32 = zeros(Float64, size(mrd_data_down2, 1))

for i in 1:size(mrd_data_up1, 1)
    try
        (lin_mrd_up_q11[i], lin_mrd_up_med1[i], lin_mrd_up_q31[i]) = quantile(filter(!isnan, mrd_data_up1[i, :]), [0.25, 0.5, 0.75])
    catch e
    end
    try
        (lin_mrd_down_q11[i], lin_mrd_down_med1[i], lin_mrd_down_q31[i]) = quantile(filter(!isnan, mrd_data_down1[i, :]), [0.25, 0.5, 0.75])
    catch e    
    end
end
for i in 1:size(mrd_data_up2, 1)
    try
       (lin_mrd_up_q12[i], lin_mrd_up_med2[i], lin_mrd_up_q32[i]) = quantile(filter(!isnan, mrd_data_up2[i, :]), [0.25, 0.5, 0.75]) 
    catch e
    end
    try
       (lin_mrd_down_q12[i], lin_mrd_down_med2[i], lin_mrd_down_q32[i]) = quantile(filter(!isnan, mrd_data_down2[i, :]), [0.25, 0.5, 0.75]) 
    catch e
    end
end

fig = PyPlot.figure(figsize=(12,4.8))
gs = gridspec.GridSpec(1,2)
ax = fig.add_subplot(gs[1,1])
ax.set_title("2 m (T2)")
fig.text(0.45, 0.94, string(Date(timemiddle_start), " 09-17LT"), size="large")
fig.text(0.12, 0.9, "a)", size="large", )
ax.plot(Dates.value.(x2[:, 1]) ./ 1000, lin_mrd_down_med2 .* multfactor, color="blue")
ax.fill_between(Dates.value.(x2[:, 1]) ./ 1000, lin_mrd_down_q12 .* multfactor, lin_mrd_down_q32 .* multfactor, color="blue", alpha=0.5)
ax.plot(Dates.value.(x2[:, 1]) ./ 1000, lin_mrd_up_med2 .* multfactor, color="red")
ax.fill_between(Dates.value.(x2[:, 1]) ./ 1000, lin_mrd_up_q12 .* multfactor, lin_mrd_up_q32 .* multfactor, color="red", alpha=0.5)
PyPlot.xscale("log")
ax.grid()
ax.set_xlabel("time scale [s]")
ax.set_ylabel(ylabelstring)
ax2 = fig.add_subplot(gs[1,2], sharey=ax)
ax2.set_title("5 m (meteo station)")
fig.text(0.544,0.9, "b)", size="large")
ax2.plot(Dates.value.(x1[:, 1]) ./ 1000, lin_mrd_down_med1 .* multfactor, color="blue")
ax2.fill_between(Dates.value.(x1[:, 1]) ./ 1000, lin_mrd_down_q11 .* multfactor, lin_mrd_down_q31 .* multfactor, color="blue", alpha=0.5)
ax2.plot(Dates.value.(x1[:, 1]) ./ 1000, lin_mrd_up_med1 .* multfactor, color="red")
ax2.fill_between(Dates.value.(x1[:, 1]) ./ 1000, lin_mrd_up_q11 .* multfactor, lin_mrd_up_q31 .* multfactor, color="red", alpha=0.5)
PyPlot.xscale("log")
ax2.grid()
ax2.set_xlabel("time scale [s]")
if fluxtype !== "wT"
    upylimmax = last(ax2.get_ylim())
    ax.set_ylim(0, minimum([upylimmax, 50]))
end
#ax2.tick_params(axis="y", labelleft=false)
#ax.legend()