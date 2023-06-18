######################################################
### PROFILES OF WIND/TEMPERATURE FOR THE EC SENSORS###
###            author: Michi Haugeneder            ###
######################################################
#=
Load data with load_data.jl
=#
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter, LsqFit
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

timestep = Millisecond(50)
ρ_air = 1.2 #kg m^{-3}
c_p = 1004 #J kg^{-1} K^{-1}
L_v = 2450e3 #J kg^{-1} (approx @0°C) 

#=
######################################################
###              TURBULENT FLUXES                  ###
######################################################
#Reynolds averaging times
ra1 = Second(30)
ra2 = Second(20)

fx1 = turb.turbflux(evaldf1, ra1)
fx2 = turb.turbflux(evaldf2, ra1)
fx3 = turb.turbflux(evaldf3, ra1)
fx4 = turb.turbflux(evaldf4, ra1)
fx5 = turb.turbflux(evaldf5, ra2)

#averaging
turb.avgflux!(fx1, ra1)
turb.avgflux!(fx2, ra1)
turb.avgflux!(fx3, ra1)
turb.avgflux!(fx4, ra1)
turb.avgflux!(fx5, ra2)

#create additional data for day and night
daystart = Time(08,00,00)
dayend = Time(18,00,00)
(fx1day, fx1night) = turb.splitdaynight(fx1, daystart, dayend)
(fx2day, fx2night) = turb.splitdaynight(fx2, daystart, dayend)
(fx3day, fx3night) = turb.splitdaynight(fx3, daystart, dayend)
(fx4day, fx4night) = turb.splitdaynight(fx4, daystart, dayend)
(evaldf1day, evaldf1night) = turb.splitdaynight(evaldf1, daystart, dayend)
(evaldf2day, evaldf2night) = turb.splitdaynight(evaldf2, daystart, dayend)
(evaldf3day, evaldf3night) = turb.splitdaynight(evaldf3, daystart, dayend)
(evaldf4day, evaldf4night) = turb.splitdaynight(evaldf4, daystart, dayend)
=#
######################################################
###                WIND PROFILE                    ###
######################################################
turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)
turb.missing2nan!(evaldf5)
turb.missing2nan!(evaldf6)

#heigths of instruments
h = [0.3, 1, 2.0, 2.85]

evaldf1umedian = median(filter(!isnan, evaldf1.u))
evaldf2umedian = median(filter(!isnan, evaldf2.u))
evaldf3umedian = median(filter(!isnan, evaldf3.u))
evaldf4umedian = median(filter(!isnan, evaldf4.u))
evaldf5umedian = median(filter(!isnan, evaldf5.u))

evaldf1uquant = abs.(reshape(quantile(filter(!isnan, evaldf1.u), [0.25, 0.75]), 2, 1) .- evaldf1umedian)
evaldf2uquant = abs.(reshape(quantile(filter(!isnan, evaldf2.u), [0.25, 0.75]), 2, 1) .- evaldf2umedian)
evaldf3uquant = abs.(reshape(quantile(filter(!isnan, evaldf3.u), [0.25, 0.75]), 2, 1) .- evaldf3umedian)
evaldf4uquant = abs.(reshape(quantile(filter(!isnan, evaldf4.u), [0.25, 0.75]), 2, 1) .- evaldf4umedian)
evaldf5uquant = abs.(reshape(quantile(filter(!isnan, evaldf5.u), [0.25, 0.75]), 2, 1) .- evaldf5umedian)

##
#fit log wind profile
logwindprof(z, p) = p[1] / 0.4 * log.(z / p[2]) #p[1]=u⋆, p[2]=z₀
p0 = [0.2, 1e-3] #initial parameters
logwindfit = curve_fit(logwindprof, h, [evaldf5umedian, evaldf2umedian, evaldf3umedian, evaldf4umedian], p0)
logwindparam = logwindfit.param
zdata = collect(0.01:0.01:2.9)

#profile plot of mean wind speeds and fit
wpfig = PyPlot.figure()
wpax = wpfig.add_subplot(111)
wpax.set_title("median (and quartile) wind speeds")
wpax.set_ylabel(L"h~\mathrm{[m]}")
wpax.set_xlabel(L"\overline{u}~\mathrm{[m~s^{-1}]}")
wpax.errorbar(evaldf2umedian, h[2], xerr=evaldf2uquant, fmt="o", label="T2IRG")
wpax.errorbar(evaldf3umedian, h[3], xerr=evaldf3uquant, fmt="o", label="T2LCSAT")
wpax.errorbar(evaldf4umedian, h[4], xerr=evaldf4uquant, fmt="o", label="T2UCSAT")
wpax.errorbar(evaldf5umedian, h[1], xerr=evaldf5uquant, fmt="o", label="KAIJO")
wpax.errorbar(evaldf1umedian, 1.2, xerr=evaldf1uquant, fmt="o", label="T1IRG (no fit)", alpha=0.5)
wpax.plot(logwindprof(zdata, logwindparam), zdata, color="black", label="log-fit")
wpax.text(0.5, 0.5,
    string(L"u^\ast = ", round(logwindparam[1], digits=4), L"~\mathrm{[m~s^{-1}]}", "\n", L"z_0 = ", round(logwindparam[2] * 1e3, digits=2), L"~\mathrm{[mm]}"))
wpax.grid()
wpax.legend()
##