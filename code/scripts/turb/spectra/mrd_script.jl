######################################################
###       MULTI RESOLUTION FLUX DECOMPOSITION      ###
###            author: Michi Haugeneder            ###
######################################################
#=
Calculate 'normal' and quasi-continuos MRD. Load data with
load_data.jl
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
animation = pyimport("matplotlib.animation")
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
include(joinpath(importdir, "src", "mrd.jl"))
import .turb
import .gen
import .MRD
PyPlot.pygui(true)

timestep = Millisecond(50)

turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)
turb.missing2nan!(evaldf5)
turb.missing2nan!(evaldf6)

######################################################
###                  MRD                           ###
######################################################
#uncomment to process the data
#Mulitresolution Flux Decomposition
##
M = 17
blocklength = ceil(2^M * timestep, Dates.Second)
(mrd_x_1, mrd_data_1, time_middle_1) = MRD.completemrd(evaldf1, "w", "T", M, round(Int, 0.1 * 2^M))
(mrd_time_1, mrd_D_median_1, mrd_D_min_1, mrd_D_max_1, mrd_D_quant1_1, mrd_D_quant3_1) =
    MRD.mrdpp(mrd_x_1, mrd_data_1)
(mrd_x_2, mrd_data_2, time_middle_2) = MRD.completemrd(evaldf2, "w", "T", M, round(Int, 0.1 * 2^M))
(mrd_time_2, mrd_D_median_2, mrd_D_min_2, mrd_D_max_2, mrd_D_quant1_2, mrd_D_quant3_2) =
    MRD.mrdpp(mrd_x_2, mrd_data_2)
(mrd_x_3, mrd_data_3, time_middle_3) = MRD.completemrd(evaldf3, "w", "T", M, round(Int, 0.1 * 2^M))
(mrd_time_3, mrd_D_median_3, mrd_D_min_3, mrd_D_max_3, mrd_D_quant1_3, mrd_D_quant3_3) =
    MRD.mrdpp(mrd_x_3, mrd_data_3)
(mrd_x_4, mrd_data_4, time_middle_4) = MRD.completemrd(evaldf4, "w", "T", M, round(Int, 0.1 * 2^M))
(mrd_time_4, mrd_D_median_4, mrd_D_min_4, mrd_D_max_4, mrd_D_quant1_4, mrd_D_quant3_4) =
    MRD.mrdpp(mrd_x_4, mrd_data_4)
(mrd_x_5, mrd_data_5, time_middle_5) = MRD.completemrd(evaldf5, "w", "T", M, round(Int, 0.1 * 2^M))
(mrd_time_5, mrd_D_median_5, mrd_D_min_5, mrd_D_max_5, mrd_D_quant1_5, mrd_D_quant3_5) =
    MRD.mrdpp(mrd_x_5, mrd_data_5)
(mrd_x_6, mrd_data_6, time_middle_6) = MRD.completemrd(evaldf6, "w", "T", M, round(Int, 0.1 * 2^M))
(mrd_time_6, mrd_D_median_6, mrd_D_min_6, mrd_D_max_6, mrd_D_quant1_6, mrd_D_quant3_6) =
    MRD.mrdpp(mrd_x_6, mrd_data_6)
#------------------------------------------------------------------------
#Plot result of MRD
PyPlot.pygui(true)
fig = PyPlot.figure()#figsize=(12,7))
ax = fig.add_subplot(111)
ax.set_title(string("MRD", evaldf1.time[1], " - ", evaldf1.time[end]))
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$C_{w\theta} [\cdot 10^{-3} \mathrm{Kms^{-1}}]$")
ax.grid(true)
PyPlot.xscale("log")
ax.plot(Dates.value.(mrd_time_1) ./ 1000, mrd_D_median_1 .* 1000, label="T1IRG")
ax.fill_between(Dates.value.(mrd_time_1) ./ 1000, mrd_D_quant1_1 .* 1000, mrd_D_quant3_1 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_2) ./ 1000, mrd_D_median_2 .* 1000, label="T2IRG")
ax.fill_between(Dates.value.(mrd_time_2) ./ 1000, mrd_D_quant1_2 .* 1000, mrd_D_quant3_2 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_3) ./ 1000, mrd_D_median_3 .* 1000, label="T2LCSAT")
ax.fill_between(Dates.value.(mrd_time_3) ./ 1000, mrd_D_quant1_3 .* 1000, mrd_D_quant3_3 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_4) ./ 1000, mrd_D_median_4 .* 1000, label="T2UCSAT")
ax.fill_between(Dates.value.(mrd_time_4) ./ 1000, mrd_D_quant1_4 .* 1000, mrd_D_quant3_4 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_5) ./ 1000, mrd_D_median_5 .* 1000, label="Kaijo")
ax.fill_between(Dates.value.(mrd_time_5) ./ 1000, mrd_D_quant1_5 .* 1000, mrd_D_quant3_5 .* 1000, alpha=0.4)
ax.plot(Dates.value.(mrd_time_6) ./ 1000, mrd_D_median_6 .* 1000, label="TJK")
ax.fill_between(Dates.value.(mrd_time_6) ./ 1000, mrd_D_quant1_6 .* 1000, mrd_D_quant3_6 .* 1000, alpha=0.4)
ax.legend()
#------------------------------------------------------------------------
######################################################
#MRD time series
timestart = Time(13, 00, 00)
timeend = Time(17, 00, 00)
daystart = Date(2021, 05, 20)
dayend = Date(2021, 06, 10)

tjkdata = turb.readturbasnetcdf(joinpath(datapath, "tower", "preproc", "tjkdf.nc"), daystart + timestart - Hour(1), dayend + timeend + Hour(1))
tjkdata = turb.interpolatemissing(tjkdata)
turb.missing2nan!(tjkdata)

starts = (daystart:Day(1):dayend) .+ timestart
ends = (daystart:Day(1):dayend) .+ timeend

M = 16
mrd_time = zeros(Millisecond, M, length(starts))
mrd_D_median = zeros(Float64, M, length(starts))
mrd_D_min = similar(mrd_D_median)
mrd_D_max = similar(mrd_D_median)
mrd_D_quant1 = similar(mrd_D_median)
mrd_D_quant3 = similar(mrd_D_median)

for idx in 1:length(starts)
    datatouse = tjkdata[starts[idx].<=tjkdata.time.<=ends[idx], :]
    (mrd_x, mrd_data, time_middle) = MRD.completemrd(datatouse, "w", "T", M, round(Int, 0.1 * 2^M))
    (mrd_time[:, idx], mrd_D_median[:, idx], mrd_D_min[:, idx], mrd_D_max[:, idx], mrd_D_quant1[:, idx], mrd_D_quant3[:, idx]) =
        MRD.mrdpp(mrd_x, mrd_data)
end

#=
function animupdate(i::Int64)
    ax.set_title(string("MRD ", starts[i], " - ", timeend))
    pm.set_ydata(mrd_D_median[:, i] .* 1000)
    plq.set_ydata(mrd_D_quant1[:, i] .* 1000)
    puq.set_ydata(mrd_D_quant3[:, i] .* 1000)
    return pm, plq, puq
end
=#

#Plot
for i in 1:length(starts)
    fig = PyPlot.figure()#figsize=(12,7))
    ax = fig.add_subplot(111)
    ax.set_title(string("TJK MRD ", starts[i], " - ", timeend))
    ax.set_xlabel("avg. time [s]")
    ax.set_ylabel(L"$C_{w\theta} [\cdot 10^{-3} \mathrm{Kms^{-1}}]$")
    #ax.set_ylim(-10,10)
    ax.grid(true)
    PyPlot.xscale("log")
    pm, = ax.plot(Dates.value.(mrd_time[:, i]) ./ 1000, mrd_D_median[:, i] .* 1000, label="TJK")
    #puq, = ax.plot(Dates.value.(mrd_time[:, i]) ./ 1000, mrd_D_quant3[:, i] .* 1000, alpha=0.5, color=pm.get_color())
    #plq, = ax.plot(Dates.value.(mrd_time[:, i]) ./ 1000, mrd_D_quant1[:, i] .* 1000, alpha=0.5, color=pm.get_color())
    fl = ax.fill_between(Dates.value.(mrd_time[:, i]) ./ 1000, mrd_D_quant1[:, i] .* 1000, mrd_D_quant3[:, i] .* 1000, alpha=0.4)
    #fig.savefig(string("/home/haugened/Documents/plots/tjk_mrd/", string(starts[i])[1:10], ".png"))
end
#ani = animation.FuncAnimation(fig, animupdate, frames=collect(1:length(starts)), interval=5000, repeat_delay=500)
#ani.save(string("/home/haugened/Desktop/", "short_video_", matfileraw, ".mp4"), fps=framespersec)

##

######################################################
###                SIMPLE MRD                      ###
######################################################
M=15
evaldf = evaldf3
evaldf = evaldf[(1:2^M),:]
data_a = evaldf.w
data_b = evaldf.T

#orthogonal MRD
D = MRD.mrd(data_a, data_b, M, 0)
x = 2 .^ collect(1:M) .* Millisecond(50)

#non-orthogonal MRD
(nox, noD) = MRD.fastnomrd(data_a, data_b, 2 .^ collect(M:-1:0))

#plot
fig = PyPlot.figure()#figsize=(12,7))
ax = fig.add_subplot(111)
ax.set_xlabel("time scale [s]")
ax.set_ylabel(L"$C_{wT_s} [10^{-3} \mathrm{Kms^{-1}}]$")
ax.grid(true)
PyPlot.xscale("log")
ax.plot(Dates.value.(x) ./ 1000, D .* 1000, label="orthogonal")
ax.plot(Dates.value.(nox) ./ 1000, noD .* 1000, label="non-orthogonal")
ax.legend()
######################################################
println("------------D-O-N-E---------------")
println()

##
#=
#comparison of MRDs
turb.missing2nan!(evaldf2)
b = evaldf2[1:2^17,:]
x1 = 2 .^ collect(1:17)*Millisecond(50)
d1 = turb.mrd(b.w, b.T, 17, 0)

(x2, d2, ds2) = turb.nomrd(b.w, b.T, 1)
(x3, d3, ds3) = turb.nomrd(b.w, b.T, 0.5)
(x4, d4, ds4) = turb.nomrd(b.w, b.T, 0.2)

#normalize
d1 ./= abs(sum(d1))
d2 ./= abs(sum(d2))
ds2 ./=abs(sum(d2))
d3 ./= abs(sum(d3))
ds3 ./=abs(sum(d3))
d4 ./= abs(sum(d4))
ds4 ./=abs(sum(d4))

fig = PyPlot.figure()
ax = fig.add_subplot(111)
fig.suptitle("Comparisons MRDs 2^17 datapts from T2IRG (31.05.21 12)")
PyPlot.xscale("log")
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$\frac{w'T'}{\int w'T'} [1]$")
ax.plot(Dates.value.(x1)./1000, d1, label="orthogonal->17pts")
ax.errorbar(Dates.value.(x2)./1000, d2, yerr=ds2, label="n-o spacing 1->17pts")
ax.errorbar(Dates.value.(x3)./1000, d3, yerr=ds3, label="n-o spacing 0.5->32pts", alpha=0.5)
ax.errorbar(Dates.value.(x4)./1000, d4, yerr=ds4, label="n-o spacing 0.2->73pts", alpha=0.5)
ax.grid()
ax.legend()
=#
##
#=
##
a = PyPlot.figure()
ax = a.add_subplot(111)
ax.set_title("T2UCSAT 02.06. 15:30-17:00 shift 0.5% - periodicity")
ax.plot(time_middle, mrd[11, :] .* 1000, label="102.4s")
#ax.plot(time_middle, mrd[10,:].*1000, label="51.2s")
#ax.plot(time_middle, mrd[7,:].*1000, label="6.4s")
#ax.plot(time_middle, mrd[14,:].*1000, label="819.2s")
#ax.plot(time_middle, mrd[13,:].*1000, label="409.6s")
ax.legend()
ax.set_ylabel(L"$C_{w\theta} [\cdot 10^{-3} \mathrm{Kms^{-1}}]$")
date_format = pydates.DateFormatter("%H:%M:%S")
ax.xaxis.set_major_formatter(date_format)
fig.autofmt_xdate()
ax.grid()
##
=#