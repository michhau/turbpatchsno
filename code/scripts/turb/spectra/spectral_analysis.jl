######################################################
###           SPECTRAL ANALYSIS OF EC-DATA         ###
###            author: Michi Haugeneder            ###
######################################################
#=
Compare spectra for different classifications. Load data with
load_data.jl
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
animation = pyimport("matplotlib.animation")
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

timestep = Millisecond(50)

turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)
turb.missing2nan!(evaldf5)
turb.missing2nan!(evaldf6)

#=
######################################################
###              CATEGORIZING DATA                 ###
######################################################

#define category conditions
function condition1(ecdata::DataFrame)
    return Time(08, 00, 00) .<= Time.(ecdata.time) .<= Time(17, 00, 00)
end

#apply categories to data
df1cond1 = evaldf1[condition1(evaldf1), :]
df2cond1 = evaldf2[condition1(evaldf2), :]
df3cond1 = evaldf3[condition1(evaldf3), :]
df4cond1 = evaldf4[condition1(evaldf4), :]
df5cond1 = evaldf5[condition1(evaldf5), :]
df6cond1 = evaldf6[condition1(evaldf6), :]
df1!cond1 = evaldf1[.!condition1(evaldf1), :]
df2!cond1 = evaldf2[.!condition1(evaldf2), :]
df3!cond1 = evaldf3[.!condition1(evaldf3), :]
df4!cond1 = evaldf4[.!condition1(evaldf4), :]
df5!cond1 = evaldf5[.!condition1(evaldf5), :]
df6!cond1 = evaldf6[.!condition1(evaldf6), :]


######################################################
###                 CREATE SPECTRA                 ###
######################################################
nrmodes = 100
(df1cond1_t, df1cond1_d_median, df1cond1_d_q1, df1cond1_d_q3, a) = turb.catdatanomrd(df1cond1, "w", "T", nrmodes)
(df2cond1_t, df2cond1_d_median, df2cond1_d_q1, df2cond1_d_q3) = turb.catdatanomrd(df2cond1, "w", "T", nrmodes)
(df3cond1_t, df3cond1_d_median, df3cond1_d_q1, df3cond1_d_q3) = turb.catdatanomrd(df3cond1, "w", "T", nrmodes)
(df4cond1_t, df4cond1_d_median, df4cond1_d_q1, df4cond1_d_q3) = turb.catdatanomrd(df4cond1, "w", "T", nrmodes)
(df5cond1_t, df5cond1_d_median, df5cond1_d_q1, df5cond1_d_q3) = turb.catdatanomrd(df5cond1, "w", "T", nrmodes)
(df6cond1_t, df6cond1_d_median, df6cond1_d_q1, df6cond1_d_q3) = turb.catdatanomrd(df6cond1, "w", "T", nrmodes)
(df1!cond1_t, df1!cond1_d_median, df1!cond1_d_q1, df1!cond1_d_q3) = turb.catdatanomrd(df1!cond1, "w", "T", nrmodes)
(df2!cond1_t, df2!cond1_d_median, df2!cond1_d_q1, df2!cond1_d_q3) = turb.catdatanomrd(df2!cond1, "w", "T", nrmodes)
(df3!cond1_t, df3!cond1_d_median, df3!cond1_d_q1, df3!cond1_d_q3) = turb.catdatanomrd(df3!cond1, "w", "T", nrmodes)
(df4!cond1_t, df4!cond1_d_median, df4!cond1_d_q1, df4!cond1_d_q3) = turb.catdatanomrd(df4!cond1, "w", "T", nrmodes)
(df5!cond1_t, df5!cond1_d_median, df5!cond1_d_q1, df5!cond1_d_q3) = turb.catdatanomrd(df5!cond1, "w", "T", nrmodes)
(df6!cond1_t, df6!cond1_d_median, df6!cond1_d_q1, df6!cond1_d_q3) = turb.catdatanomrd(df6!cond1, "w", "T", nrmodes)


#----------------------------------------------------
toplot_x = df6cond1_t
toplot_!x = df6!cond1_t
toplot_dmed = df6cond1_d_median
toplot_!dmed = df6!cond1_d_median
toplot_dq1 = df6cond1_d_q1
toplot_!dq1 = df6!cond1_d_q1
toplot_dq3 = df6cond1_d_q3
toplot_!dq3 = df6!cond1_d_q3
label_!cond = "1700-0800LT"
label_cond = "0800-1700LT"
#plot
PyPlot.pygui(true)
fig = PyPlot.figure()#figsize=(12,7))
ax = fig.add_subplot(111)
ax.set_title(string("NO-MRD TJK 26.5. - 10.6."))
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$\mathrm{cospectral~density~\frac{w'T'}{\int w'T'}}$")
ax.grid(true)
PyPlot.xscale("log")
#ax.plot(Dates.value.(mrd_time_1) ./ 1000, mrd_D_median_1 .* 1000, label="ref")
#ax.fill_between(Dates.value.(mrd_time_1) ./ 1000, mrd_D_quant1_1 .* 1000, mrd_D_quant3_1 .* 1000, alpha=0.4)
ax.plot(Dates.value.(toplot_!x) ./ 1000, toplot_!dmed, label=label_!cond)
ax.fill_between(Dates.value.(toplot_!x) ./ 1000, toplot_!dq1, toplot_!dq3, alpha=0.4)
ax.plot(Dates.value.(toplot_x) ./ 1000, toplot_dmed, label=label_cond)
ax.fill_between(Dates.value.(toplot_x) ./ 1000, toplot_dq1, toplot_dq3, alpha=0.4)
ax.legend()
PyPlot.show()

#----------------------------------------------------
#plot without condition
nrmodes = 100
(df1_t, df1_d_med, df1_q1, df1_q3) = turb.catdatanomrd(evaldf1, "w", "T", nrmodes)
(df2_t, df2_d_med, df2_q1, df2_q3) = turb.catdatanomrd(evaldf2, "w", "T", nrmodes)
(df3_t, df3_d_med, df3_q1, df3_q3) = turb.catdatanomrd(evaldf3, "w", "T", nrmodes)
(df4_t, df4_d_med, df4_q1, df4_q3) = turb.catdatanomrd(evaldf4, "w", "T", nrmodes)
(df5_t, df5_d_med, df5_q1, df5_q3) = turb.catdatanomrd(evaldf5, "w", "T", nrmodes)
(df6_t, df6_d_med, df6_q1, df6_q3) = turb.catdatanomrd(evaldf6, "w", "T", nrmodes)

PyPlot.pygui(true)
fig = PyPlot.figure()#figsize=(8,5.5))
ax = fig.add_subplot(111)
ax.set_title(string("NO-MRD Duerrboden flux profile 01.06. 1100-1700 (no norm)"))
ax.set_xlabel("avg. time [s]")
ax.set_ylabel(L"$\mathrm{cospectral~density~w'T'}$")#\frac{w'T'}{\sum w'T'}}$")
ax.grid(true)
PyPlot.xscale("log")
#ax.plot(Dates.value.(df1_t) ./ 1000, df1_d_med , label="1.2m")
#ax.plot(Dates.value.(df5_t) ./ 1000, df5_d_med , label="0.3m")
ax.plot(Dates.value.(df2_t) ./ 1000, df2_d_med, label="1.2m, T2")
ax.plot(Dates.value.(df3_t) ./ 1000, df3_d_med, label="2.5m")
#ax.plot(Dates.value.(a2_t) ./ 1000, b2_t , label="31.05.")
ax.plot(Dates.value.(df4_t) ./ 1000, df4_d_med, label="3.5m")
ax.plot(Dates.value.(df6_t) ./ 1000, df6_d_med, label="5m")
ax.legend()
PyPlot.show()

#----------------------------------------------------