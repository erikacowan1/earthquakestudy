import argparse
import matplotlib
matplotlib.use('Agg')
import pylab as plt
from gwpy.timeseries import TimeSeries
from gwpy.timeseries import TimeSeriesDict
from gwpy.plotter import TimeSeriesPlot
import ConfigParser
matplotlib.rcParams['agg.path.chunksize'] = 10000
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# command line parsing
parser = argparse.ArgumentParser(
        description='For questions or concerns, contact Erika Cowan at \
        erika.cowan@ligo.org')
parser.add_argument('EQ_start_time', type=int, help='Please enter GPS start time')
parser.add_argument('EQ_end_time', type=int, help='Please enter GPS end time')
parser.add_argument('channel1', type=str, help='Please enter valid channel')
parser.add_argument('channel2', type=str, help='Please enter valid channel')
parser.add_argument('time_win', type=int, help='Please enter valid time window around EQ')
parser.add_argument('date_str', type=str, help='Please enter string in YYYY-MM-DD Format')
args = parser.parse_args()

EQ_before = args.EQ_start_time - args.time_win
EQ_after = args.EQ_end_time + args.time_win

time_segs= [
    [EQ_before, args.EQ_start_time],
    [args.EQ_start_time, args.EQ_end_time],
    [args.EQ_end_time, EQ_after]
]

def timesseriesdict(start, end):
    return TimeSeriesDict.get([args.channel1, args.channel2], start, end, verbose=True, frametype='H1_R')

def common_mode(E_ARM, I_ARM):
    return 0.5*(E_ARM + I_ARM)

def differential_mode(E_ARM, I_ARM):
    return 0.5*(E_ARM - I_ARM)
def plot_2(p1, p2, start, end, lgd1, lgd2, title, x_label, y_label, fig_label):
    f, ax = plt.subplots()
    ax.plot(p1 , label= lgd1) 
    ax.plot(p2 , label= lgd2) 
    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    ax.legend()
    if plt_type == "ASD":
        ax.set_xscale('log')
        ax.set_yscale('log')
        f.savefig(fig_label + '.png')
        plt.close(f)
    else:
        f.savefig(fig_label + '.png')
        plt.close(f)
'''
#plt_type=ASD, BLRMS, TimeSeries, AVG_FFT
data = TimeSeriesDict.get([args.channel1, args.channel2], args.EQ_start_time, args.EQ_end_time, verbose=True, frametype='H1_R')

    #grabbing individual data from TimeSeriesDict
c1_b = data[args.channel1]
c2_b = data[args.channel2]

plt_type = 'TimeSeries'
c1_b_name = c1.name
c2_b_name = c2.name
c = c1_b_name.replace("_", "\_")
d = c2_b_name.replace("_", "\_")
plot_2(c1_b, c2_b, EQ_before, EQ_after, c, d, plt_type, args.date_str + ", time(sec)" , "velocity (nm/s)", "total" + plt_type + str(EQ_before) + "\_" + str(EQ_after))
'''
dict_data = []
for start, end in time_segs:
        data = TimeSeriesDict.get([args.channel1, args.channel2], start, end, verbose=True, frametype='H1_R')
        #print data
        dict_data.append(data)
print dict_data[0][args.channel1]
print dict_data[1][args.channel1]
print time_segs[0][0]
#plotting before, during, after the EQ       
    
#plt_type=ASD, BLRMS, TimeSeries, AVG_FFT
#data = TimeSeriesDict.get([args.channel1, args.channel2], start, end, verbose=True, frametype='H1_R')
#grabbing individual data from TimeSeriesDict
c1_b = dict_data[0][args.channel1]
c2_b = dict_data[0][args.channel2]
c1_d = dict_data[1][args.channel1]
c2_d = dict_data[1][args.channel2]
c1_a = dict_data[2][args.channel1]
c2_a = dict_data[2][args.channel2]

#plotting TimeSeries of E and I
plt_type = 'TimeSeries'
c1_b_name = c1_b.name
c2_b_name = c2_b.name
c1_d_name = c1_d.name
c2_d_name = c2_d.name
c1_a_name = c1_a.name
c2_a_name = c2_a.name
c = c1_b_name.replace("_", "\_")
d = c2_b_name.replace("_", "\_")
c2 = c1_d_name.replace("_", "\_")
d2 = c2_d_name.replace("_", "\_")
c3 = c1_a_name.replace("_", "\_")
d3 = c2_a_name.replace("_", "\_")


#grabbing individual data from TimeSeriesDict
plt_type = 'TimeSeries'
plot_2(c1_b, c2_b, time_segs[0][0], time_segs[0][1], c, d, plt_type + " Before EQ", args.date_str + ", time(sec)" , "velocity (nm/s)", "before_" +  plt_type + str(time_segs[0][0]) + "_" + str(time_segs[0][1]))
plot_2(c1_d, c2_d, time_segs[1][0], time_segs[1][1], c2, d2, plt_type + " During EQ", args.date_str + ", time(sec)" , "velocity (nm/s)", "during_" +  plt_type + str(time_segs[1][0]) + "_" + str(time_segs[1][1]))
plot_2(c1_a, c2_a, time_segs[2][0], time_segs[2][1], c3, d3, plt_type + " After EQ", args.date_str + ", time(sec)" , "velocity (nm/s)", "after_" +  plt_type + str(time_segs[2][0]) + "_" + str(time_segs[2][1]))

#plotting TimeSeries CM, DM
plt_type = 'CM\_DM\_TimeSeries'
plot_2(common_mode(c1_b, c2_b),differential_mode(c1_b, c2_b), time_segs[0][0], time_segs[0][1], "CM", "DM", plt_type + " CM\_DM" + " Before EQ", args.date_str + ", time(sec)","velocity (nm/s)", "before_CM_DM_TimeSeries_" + str(time_segs[0][0]) + "_" + str(time_segs[0][1]))
plt_type = 'CM\_DM\_TimeSeries'
plot_2(common_mode(c1_d, c2_d),differential_mode(c1_d, c2_d), time_segs[1][0], time_segs[1][1], "CM", "DM", plt_type + " CM\_DM" + " During EQ", args.date_str + ", time(sec)","velocity (nm/s)", "during_CM_DM_TimeSeries_" + str(time_segs[1][0]) + "_" + str(time_segs[1][1]))
plt_type = 'CM\_DM\_TimeSeries'
plot_2(common_mode(c1_a, c2_a),differential_mode(c1_a, c2_a), time_segs[2][0], time_segs[2][1], "CM", "DM", plt_type + " CM\_DM", args.date_str + ", time(sec)" + " After EQ","velocity (nm/s)", "after_CM_DM_TimeSeries_" + str(time_segs[2][0]) + "_" + str(time_segs[2][1]))


#plotting BLRMS
plt_type = 'BLRMS'
freq = [[.03,.1],[.1,.3],[.3,1]]
for start_freq, stop_freq in freq:
    blrms1 = c1_b.bandpass(start_freq, stop_freq, 3, 30).rms()
    blrms2 = c2_b.bandpass(start_freq, stop_freq, 3, 30).rms()
    blrms3 = c1_d.bandpass(start_freq, stop_freq, 3, 30).rms()
    blrms4 = c2_d.bandpass(start_freq, stop_freq, 3, 30).rms()
    blrms5 = c1_a.bandpass(start_freq, stop_freq, 3, 30).rms()
    blrms6 = c2_a.bandpass(start_freq, stop_freq, 3, 30).rms()
    plot_2(blrms1, blrms2, time_segs[0][0], time_segs[0][1], "CM", "DM", plt_type + " CM\_DM " + str(start_freq) + "Hz\_" + str(stop_freq) + "Hz" + " Before EQ", args.date_str + ", time(sec)", "velocity (nm/s)", "before_" + plt_type + "_CM_DM" + str(time_segs[0][0]) + "_" + str(time_segs[0][1]) + "_" + str(start_freq) + "Hz_" + str(stop_freq) + "Hz")
    plot_2(blrms3, blrms4, time_segs[1][0], time_segs[1][1], "CM", "DM", plt_type + " CM\_DM " + str(start_freq) + "Hz\_" + str(stop_freq) + "Hz" + " During EQ", args.date_str + ", time(sec)", "velocity (nm/s)", "during_" + plt_type + "_CM_DM" + str(time_segs[1][0]) + "_" + str(time_segs[1][1]) + "_" + str(start_freq) + "Hz_" + str(stop_freq) + "Hz")
    plot_2(blrms5, blrms6, start, end, "CM", "DM", plt_type + " CM\_DM " + str(start_freq) + "Hz\_" + str(stop_freq) + "Hz" + " After EQ", args.date_str + ", time(sec)", "velocity (nm/s)", "after_" + plt_type + "_CM_DM" + str(time_segs[2][0]) + "_" + str(time_segs[2][1]) + "_" + str(start_freq) + "Hz_" + str(stop_freq) + "Hz")

#plotting ASD
plt_type = 'ASD'
asd1 = common_mode(c1_b, c2_b).asd(4,2)
asd2 = differential_mode(c1_b, c2_b).asd(4,2)
asd3 = common_mode(c1_d, c2_d).asd(4,2)
asd4 = differential_mode(c1_d, c2_d).asd(4,2)
asd5 = common_mode(c1_a, c2_a).asd(4,2)
asd6 = differential_mode(c1_a, c2_a).asd(4,2)
plot_2(asd1, asd2, time_segs[0][0], time_segs[0][1], "CM", "DM", plt_type + " CM\_DM" + " Before EQ", args.date_str + ", frequency(Hz)", "ASD velocity [1/$\sqrt{\mathrm{Hz}}$]", "before_" + plt_type + "_CM_DM" + str(time_segs[0][0]) + "_" + str(time_segs[0][1]))
plot_2(asd3, asd4, time_segs[1][0], time_segs[1][1], "CM", "DM", plt_type + " CM\_DM", args.date_str + ", frequency(Hz)" + " During EQ", "ASD velocity [1/$\sqrt{\mathrm{Hz}}$]", "during_" + plt_type + "_CM_DM" + str(time_segs[1][0]) + "_" + str(time_segs[1][1]))
plot_2(asd5, asd6, time_segs[2][0], time_segs[2][1], "CM", "DM", plt_type + " CM\_DM", args.date_str + ", frequency(Hz)" + " After EQ", "ASD velocity [1/$\sqrt{\mathrm{Hz}}$]", "after_" + plt_type + "_CM_DM" + str(time_segs[2][0]) + "_" + str(time_segs[2][1]))

'''
f, ax = plt.subplots()
ax.plot(c1_b)
ax.plot(c2_b)
ax.plot(c1_d) 
ax.plot(c2_d)
#plot_2(c1_b, c2_b, time_segs[0][0], time_segs[0][1], c, d, plt_type, args.date_str + ", time(sec)" , "velocity (nm/s)", plt_type + str(time_segs[0][0]) + "\_" + str(time_segs[0][1]))
#plot_2(c1_d, c2_d, time_segs[1][0], time_segs[1][1], c2, d2, plt_type, args.date_str + ", time(sec)" , "velocity (nm/s)", plt_type + str(time_segs[1][0]) + "\_" + str(time_segs[1][1]))
f.savefig('Test.png')
#plot_2(c1_b, c2_b, start, end, c, d, plt_type, args.date_str + ", time(sec)" , "velocity (nm/s)", plt_type + str(start) + "\_" + str(end))
'''
#plotting AVG_FFT CM, DM 
#plt_type = 'AVG\_FFT'
#fft_p1 = common_mode(c1_b, c2_b).average_fft(500, 50, window='hamming')
#fft_p2 = differential_mode(c1_b, c2_b).average_fft(500, 50, window='hamming')
#plot_2(fft_p1, fft_p2, start, end, "CM", "DM", plt_type + " CM\_DM", args.date_str + ", time(sec)", "velocity (nm/s)", "AVG_FFT_CM_DM" + str(start) + "_" + str(end))

