######################################################
###            PROCESS TURBULENCE DATA             ###
###            author: Michi Haugeneder            ###
######################################################
#=
NOTE: Rawdata has to be read-in with read_raw_data.jl before!
      This script works with exported .csv data-files and converts
      them to DataFrame objects

Workflow to my current (210709) understanding:
- import .csv file to create DataFrame objects
- read-in the measurement period information using turb.readperiodfile(file)
- apply MRD to (for example 30min) blocks
- define cospectral gap T_g from MRD analysis and check eval_block
- actual flux calculation for blocks of T_g (with DR)
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
#animation = pyimport("matplotlib.animation")
#GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")
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

#########################################################################
# change variables here
#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)

println()
println("-----------S-T-A-R-T-------------")
#########################################################################
###                        DATA SELECTION                             ###
#########################################################################
#select data and measurement period to be evaluated
evaldf = kaijodf
#evaldf = t1irg

evalstart = DateTime(2021,04,28,10,57,07)
#evalend = DateTime(2021,04,23,12,30,00)
evalend = evalstart + Minute(10)

evaldf = evaldf[evalstart .<= evaldf[:,1] .<= evalend,:]
turb.nan2missing!(evaldf)
#replace bad data with missing
turb.qualcontrolflags!(evaldf)
#show statistics about missing data
turb.printmissstats(evaldf)
#interpolate missing
evaldf = turb.interpolatemissing(evaldf)
#detrend
evaldf.u = turb.detrend(evaldf.u)
evaldf.v = turb.detrend(evaldf.v)
evaldf.w = turb.detrend(evaldf.w)
#double rotation if needed
turb.drdf!(evaldf)
##
#=
#for Kaijo
#read-in information about the measurement periods
measureperiod = turb.readperiodfile(kaijo_period_file)
#select measurement period number for following evaluation (check DataFrame for times)
sel_per = 4

evaldf = evaldf[measureperiod[sel_per,1] .<= evaldf[:,1] .<= measureperiod[sel_per,2],:]
=#
#########################################################################
#continuous flux calculation

turb.missing2nan!(evaldf)
cflux = gen.movingaverage(turb.contflux(evaldf.w, evaldf.T, 600), 200)
cfluxmean = mean(cflux)
##
PyPlot.pygui(true)
cfig = PyPlot.figure()
cax = cfig.add_subplot(111)
cax.set_xlabel("time")
dateform = pydates.DateFormatter("%H:%M")
cax.xaxis.set_major_formatter(dateform)
cax.set_xticks([DateTime(2021,04,28,10,58,00), DateTime(2021,04,28,11,00,00), DateTime(2021,04,28,11,02,00), DateTime(2021,04,28,11,04,00), DateTime(2021,04,28,11,06,00), DateTime(2021,04,28,11,08,00)])
cax.set_ylabel(L"$w'\theta' [\cdot 10^{-3} \mathrm{Kms^{-1}}]$")
cax.grid(true)
cax.axhline(cfluxmean, ls="--", alpha = 0.5)
cax.axvline(DateTime(2021,04,28,11,00,09), color="black")
cax.text(DateTime(2021,04,28,11,00,02),0.13,L"\mathrm{\mathbf{t_1}}", fontweight="bold")
cax.axvline(DateTime(2021,04,28,11,05,58), color="black")
cax.text(DateTime(2021,04,28,11,05,51),0.13,L"\mathrm{\mathbf{t_2}}", fontweight="bold")
cax.plot(evaldf.time, cflux, label="")
PyPlot.show()
##
#########################################################################
###                       BLOCK EVALUATION                            ###
#########################################################################
#uncomment to process the data

#Mulitresolution Flux Decomposition
M=11
blocklength = ceil(2^M*timestep, Dates.Second)
(mrd_x, mrd_data, time_middle) = turb.completemrd(evaldf, M, 2^M)
(mrd_time,mrd_D_median,mrd_D_min,mrd_D_max,mrd_D_quant1, mrd_D_quant3) = 
turb.mrdpp(mrd_x, mrd_data)
#------------------------------------------------------------------------
#Plot result of MRD
PyPlot.pygui(true)
fig = PyPlot.figure(figsize=(12,7))
ax = fig.add_subplot(111)
ax.set_title(string("MRFD T2 ", evaldf[1,1], " - ", evaldf[end,1]))
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$C_{w\theta} [\cdot 10^{-3} \mathrm{Kms^{-1}}]$")
ax.grid(true)
PyPlot.xscale("log")
#ax.plot(mrd_kaijo_time, mrd_D_kaijo_median.*1000, label="KAIJO")
#ax.fill_between(mrd_kaijo_time, mrd_D_kaijo_quant1.*1000, mrd_D_kaijo_quant3.*1000, alpha=0.4)
ax.plot(mrd_time, mrd_D_median.*1000, label="")
ax.fill_between(mrd_time, mrd_D_quant1.*1000, mrd_D_quant3.*1000, alpha=0.4)
#=ax.plot(mrd_t2irg_time, mrd_D_t2irg_median.*1000, label="T2IRG")
ax.fill_between(mrd_t2irg_time, mrd_D_t2irg_quant1.*1000, mrd_D_t2irg_quant3.*1000, alpha=0.4)
ax.plot(mrd_t2lcsat_time, mrd_D_t2lcsat_median.*1000, label="T2LCSAT")
ax.fill_between(mrd_t2lcsat_time, mrd_D_t2lcsat_quant1.*1000, mrd_D_t2lcsat_quant3.*1000, alpha=0.4)
ax.plot(mrd_t2ucsat_time, mrd_D_t2ucsat_median.*1000, label="T2UCSAT")
ax.fill_between(mrd_t2ucsat_time, mrd_D_t2ucsat_quant1.*1000, mrd_D_t2ucsat_quant3.*1000, alpha=0.4)
=#ax.legend()
PyPlot.show()
#PyPlot.plot(mrd_time, (mrd_D_mean+mrd_D_std).*100)
#------------------------------------------------------------------------
#########################################################################
(kaijo_block, datarot) = turb.blockevaluation(evaldf, Second(5), timestep)

#write block DataFrame to file
#CSV.write(string(outfile_stam, "_block.csv"), eval_block)

#write MRD-block-DataFrame to file
#CSV.write(string(outfile_stam, "_block_mrd.csv"), eval_block_mrd)

#########################################################################
#=
#read in the already exported files
kaijodf = turb.csvtodataframe(string(kaijo_outfile_stam, ".csv"))
eval_block = turb.csvtodataframe(string(outfile_stam, "_block.csv"))
eval_block_mrd = turb.csvtodataframe(string(outfile_stam, "_block_mrd.csv"))
=#
#########################################################################
###                           Plotting                                ###
#########################################################################

#########################################################################
#simple plotting stuff
turb.missing2nan!(evaldf)
PyPlot.pygui(true)
fig2 = PyPlot.figure()
ax = fig2.add_subplot(111)
ax.set_title("Pardenn Kaijo")
ax.set_xlabel("time")
ax.set_ylabel(L"\left(\overline{w'\theta'}\right)_{30s} \mathrm{[\frac{Km}{s}]}")
ax.grid(true)
#ax.set_ylim((-0.15,0.25))
#ax.plot(kaijo_block.time_middle, kaijo_block.wT, label="kaijo")
#ax.plot(evaldf.time, evaldf.w)
#ax.plot(t2lcsat_block.time_middle, t2lcsat_block.wT, label="lower CSAT", alpha=0.5)
#ax.plot(t2ucsat_block.time_middle, t2ucsat_block.wT, label="upper CSAT", alpha=0.5)
ax.legend()
#########################################################################
#plot flux profile
turb.missing2nan!(kaijo_block)
#turb.missing2nan!(t1irg_block)
#turb.missing2nan!(t2irg_block)
#turb.missing2nan!(t2lcsat_block)
#turb.missing2nan!(t2ucsat_block)

PyPlot.pygui(true)
fig2 = PyPlot.figure()
ax = fig2.add_subplot(111)
ax.set_title("Flux profile Dürrboden 03.06.2021")
ax.set_xlabel("time")
ax.set_ylabel(L"\left(\overline{w'\theta'}\right)_{1min} \mathrm{[\frac{Km}{s}]}")
ax.grid(true)
#ax.set_ylim((-0.15,0.25))
ax.plot(kaijo_block.time_middle, kaijo_block.wT, label="Kaijo 5s")
#ax.plot(t1irg_block.time_middle, t1irg_block.wT, label="T1IRG")
#ax.plot(t2irg_block.time_middle, t2irg_block.wT, label="T2IRG")
#ax.plot(t2lcsat_block.time_middle, t2lcsat_block.wT, label="T2lCSAT", alpha=0.7)
#ax.plot(t2ucsat_block.time_middle, t2ucsat_block.wT, label="T2uCSAT", alpha=0.7)
ax.legend()

#########################################################################
PyPlot.pygui(true)
#plotting wind direction histograms
periodstart = DateTime(2021,06,08,10,00,00)
periodend = DateTime(2021,06,11,08,17,15)

ec1 = evaldf4
ec2tmp = t2irgdb

ec1 = ec1tmp[periodstart .<= ec1tmp.time .<= periodend,:]
ec2 = ec2tmp[periodstart .<= ec2tmp.time .<= periodend,:]
turb.drdf!(ec1)
turb.drdf!(ec2)
#ec3 = tjkdf[periodstart .<= tjkdf.time .<= periodend,:]
#tjkmeteo = tjkmeteodata[]

winddir1 = turb.winddir(ec1)
winddir2 = turb.winddir(ec2)
#winddir3 = turb.winddir(tjkdf)

winddirfig = PyPlot.figure()
dirax = winddirfig.add_subplot(111)
dirax.set_xlabel("wind direction [deg in instrument coord.]")
dirax.set_title("Wind direction histogram 8.6.-11.6. - before DR")
dirax.hist(winddir1.α, bins=collect(0:0.5:360), density=true, label="T2UCSAT", alpha=0.5)
#dirax.hist(winddir2.α, bins=collect(0:0.5:360), density=true, label="T2IRG", alpha=0.5)
#dirax.hist(tjkmeteo.wind_mean_vector_direction, bins=collect(0:2:360), density=true, label="TJK_prop", alpha=0.5)
#dirax.hist(winddir3.α, bins=collect(0:0.5:360), density=true, label="TJK_3DUS", alpha=0.5)
dirax.legend()
PyPlot.show()

PyPlot.plot(winddir1.time[1:10000], abs.(winddir1.α[1:10000]-winddir2.α[1:10000]))

#########################################################################
println("------------D-O-N-E---------------")
println()
