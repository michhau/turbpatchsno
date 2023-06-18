######################################################
###              MRD VERCAUTEREN PLOT              ###
###            author: Michi Haugeneder            ###
######################################################
#=
Calculate MRD with std at different scales according to
Nilsson et al. 2014 QJRMS
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

evaldf = evaldf4

turb.missing2nan!(evaldf)

######################################################
###                  MRD                           ###
######################################################
#uncomment to process the data
M = 15
starttime = DateTime(2021,06,06,14,00,00)
endtime = DateTime(2021,06,06,18,00,00)

nans = isnan.(evaldf.u)
data_tmp = evaldf[.!nans, :]

data_tmp = data_tmp[starttime - Minute(1) .<= data_tmp.time .<= endtime + Minute(1), :]

D = MRD.mrd_std_block(data_tmp, starttime, endtime, 15, "u", "w", 2)

fig = PyPlot.figure()
ax = fig.add_subplot(111)
ax.set_title("T2UCSAT 06.06. 14-18, M=15")
mrd_std_plot = ax.pcolormesh(collect(0.5:14.5), collect(0.5:14.5), D)
ax.set_xlabel("w^2-scale [2^x data points]")
ax.set_ylabel("u-scale [2^x data points]")
cbar = fig.colorbar(mrd_std_plot, ax=ax)#, orientation = "horizontal")