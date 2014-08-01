radar
=====

A bunch of these files are one time use scripts, but I included them on this repository anyway. The main file that does the "scene" dection currently is [main_scene_angle.m](main_scene_angle.m). I'll list each important file/directory here with a brief description

### [different_subcarriers](different_subcarries)
  A directory with a bunch of OFDM signals I generated and tested. Bad form to put data files on Github, but oh well. Each file is a csv that contains the time domain representation of some OFDM signal at 1 GS/s.
  
### [fft_problem](fft_problem)
  A one time use directory, probably shouldn't ever need again. It was used when I determined that there was some sort of channel response. This code was used to debug what was actually causing that.
  
### [good_signals](good_signals)
  Storing more csvs on Github, oh well. Each of the tx csvs has a bunch of leading zeroes, this was used to "bypass" the effects of the channel response. 
  - [good_subcarriers](good_signals/good_subcarriers) - This directory contains experimental data from 3 OFDM signals that I found to have the best returns.
  - [tx_signal_42.txt](good_signals/tx_signal_42.txt) - My personal favorite and new go-to signal.
  - [tx_signal_46.txt](good_signals/tx_signal_46.txt) - Another signal that performed well.
  - [tx_signal_47.txt](good_signals/tx_signal_47.txt) - Another signal that performed well.
  - 

### [scripts](scripts)
  Commonly used scripts, mostly downsampling and reading in data

#### [channel_response.m](channel_response.m)
#### [check_downsampling.m](check_downsampling.m)
#### [mainDetect.m](mainDetect.m)
#### [mainDetect_once.m](mainDetect_once.m)
#### [mainLateral.m](mainLateral.m)
#### [mainLateralBob.m](mainLateralBob.m)
#### [mainLateral_all.m](mainLateral_all.m)
