#from windfield_kaga.jl importing the CSAT-data for HEFEX

rawdatafile = "F:\\measurements\\HEFEX_IR_measurements\\rawdata_CSAT\\TOB1_3643.Wind_50ms_20180822_all.dat"
#times for the extraction
starttime = DateTime(2018,08,22,11,31,05)
endtime = DateTime(2018,08,22,11,31,17)

"Load 3D-ultrasonic (A) raw data -> timestamp, w_raw, T_raw"
function loadultrasonicrawdata(
    rawdatafile::String, starttime::DateTime, endtime::DateTime)
    println("Loading ultrasonic (A) raw data for given time period")
    df = readdlm(rawdatafile, ',', String, '\n'; skipstart = 4, use_mmap = true)
    #extract first column to be DateTime
    dateformat = DateFormat("yyyy-mm-dd HH:MM:SS.ss")
    timeofmeasure = DateTime.(df[:, 1], dateformat)
    #extract only used time period and w_raw and T_raw
    df = parse.(Float64, df[starttime .< timeofmeasure .< endtime, [2,3,4,5]])
    println("Done")
    return df, timeofmeasure[starttime .< timeofmeasure .< endtime]
end

#load ultrasonicdata
#CSATmeanhorwndspd = mean(loadultrasonicrawdata(rawdatafile, starttime, endtime)[1][:,2])

######################################################
