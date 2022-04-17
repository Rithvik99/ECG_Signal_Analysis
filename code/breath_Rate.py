import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.misc import electrocardiogram
from scipy.interpolate import splrep, splev
from mne.filter import filter_data, resample
from scipy.signal import detrend, find_peaks

sns.set(context='talk')
from mat4py import loadmat

data = loadmat("s0064lrem.mat")

ecg = data["val"][0][0:115200]+data["val"][1][0:115200]+data["val"][2][0:69600]
for i in range(0,300000):
    ecg[i] = ecg[i]/2000
sf_ori = 1000
sf = 100
dsf = sf / sf_ori
ecg = resample(ecg, dsf)
ecg = filter_data(ecg, sf, 2, 30, verbose=0)

# Select only a 20 sec window
window = 100
start = 0
ecg = ecg[int(start*sf):int((start+window)*sf)]

# R-R peaks detection
rr, _ = find_peaks(ecg, distance=40, height=-0.1)
plt.figure(0)
plt.plot(ecg)
plt.figure(0)
plt.plot(rr, ecg[rr], 'o')
plt.title('ECG signal')
plt.xlabel('Samples')
_ =plt.ylabel('Voltage')

rr = (rr / sf) * 1000
rri = np.diff(rr)

# Interpolate and compute HR
def interp_cubic_spline(rri, sf_up=4):
    
    rri_time = np.cumsum(rri) / 1000.0
    time_rri = rri_time - rri_time[0]
    time_rri_interp = np.arange(0, time_rri[-1], 1 / float(sf_up))
    tck = splrep(time_rri, rri, s=0)
    rri_interp = splev(time_rri_interp, tck, der=0)
    return rri_interp

sf_up = 4
rri_interp = interp_cubic_spline(rri, sf_up) 
hr = 1000 * (60 / rri_interp)

edr = detrend(hr)
edr = (edr - edr.mean()) / edr.std()

# Find respiratory peaks
resp_peaks, _ = find_peaks(edr, height=0, distance=sf_up)

# Convert to seconds
resp_peaks = resp_peaks
resp_peaks_diff = np.diff(resp_peaks) / sf_up

# Plot the EDR waveform
plt.figure(1)
plt.plot(edr, '-')
plt.figure(1)
plt.plot(resp_peaks, edr[resp_peaks], 'o')
_ = plt.title('ECG derived respiration')

mresprate = resp_peaks.size / window
print('Mean respiratory rate: %.2f Hz' % mresprate)
print('Mean respiratory period: %.2f seconds' % (1 / mresprate))
print('Respiration RMS: %.2f seconds' % np.sqrt(np.mean(resp_peaks_diff**2)))
print('Respiration STD: %.2f seconds' % np.std(resp_peaks_diff))


from scipy.signal import periodogram
freqs, psd = periodogram(edr, sf_up)
plt.figure(2)
plt.plot(freqs, psd)
plt.ylabel('Power spectral density')
plt.xlabel('Frequencies')
plt.title('EDR spectrum')
print('Maximum frequency: %.2f Hz' % freqs[np.argmax(psd)])

plt.show()