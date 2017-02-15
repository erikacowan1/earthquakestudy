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
    ax.plot(p1 ,label= lgd1) 
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
c1 = data[args.channel1]
c2 = data[args.channel2]

plt_type = 'TimeSeries'
c1_name = c1.name
c2_name = c2.name
c = c1_name.replace("_", "\_")
d = c2_name.replace("_", "\_")
plot_2(c1, c2, EQ_before, EQ_after, c, d, plt_type, "2017-01-04, time(sec)" , "velocity (nm/s)", "total" + plt_type + str(EQ_before) + "\_" + str(EQ_after))
'''

#plotting before, during, after the EQ       
for start, end in time_segs:
    
    #plt_type=ASD, BLRMS, TimeSeries, AVG_FFT
    data = TimeSeriesDict.get([args.channel1, args.channel2], start, end, verbose=True, frametype='H1_R')
    
    #grabbing individual data from TimeSeriesDict
    c1 = data[args.channel1]
    c2 = data[args.channel2]

    #plotting TimeSeries of E and I
    plt_type = 'TimeSeries'
    c1_name = c1.name
    c2_name = c2.name
    c = c1_name.replace("_", "\_")
    d = c2_name.replace("_", "\_")
    plot_2(c1, c2, start, end, c, d, plt_type, "2017-01-04, time(sec)" , "velocity (nm/s)", plt_type + str(start) + "\_" + str(end))
  
    #plotting TimeSeries CM, DM
    plt_type = 'CM\_DM\_TimeSeries'
    plot_2(common_mode(c1, c2),differential_mode(c1, c2), start, end, "CM", "DM", plt_type + " CM\_DM", "2017-01-04, time(sec)","velocity (nm/s)", "CM_DM_TimeSeries_" + str(start) + "_" + str(end))
    '''
    #plotting AVG_FFT CM, DM 
    plt_type = 'AVG\_FFT'
    fft_p1 = common_mode(c1, c2).average_fft(500, 50, window='hamming')
    fft_p2 = differential_mode(c1, c2).average_fft(500, 50, window='hamming')
    plot_2(fft_p1, fft_p2, start, end, "CM", "DM", plt_type + " CM\_DM", "2017-01-04, time(sec)", "velocity (nm/s)", "AVG_FFT_CM_DM" + str(start) + "_" + str(end))
    '''

    #plotting BLRMS
    plt_type = 'BLRMS'
    freq = [[.03,.1],[.1,.3],[.3,1]]
    for start_freq, stop_freq in freq:
        blrms1 = c1.bandpass(start_freq, stop_freq, 3, 30).rms()
        blrms2 = c2.bandpass(start_freq, stop_freq, 3, 30).rms()
        plot_2(blrms1, blrms2, start, end, "CM", "DM", plt_type + " CM\_DM " + str(start_freq) + "Hz\_" + str(stop_freq) + "Hz", "2017-01-04, time(sec)", "velocity (nm/s)", plt_type + "_CM_DM" + str(start) + "_" + str(end) + "_" + str(start_freq) + "Hz_" + str(stop_freq) + "Hz")
    #plotting ASD
    plt_type = 'ASD'
    asd1 = common_mode(c1, c2).asd(4,2)
    asd2 = differential_mode(c1, c2).asd(4,2)
    plot_2(asd1, asd2, start, end, "CM", "DM", plt_type + " CM\_DM", "2017-01-04, frequency(Hz)", "ASD velocity [1/$\sqrt{\mathrm{Hz}}$]", plt_type + "_CM_DM" + str(start) + "_" + str(end))
