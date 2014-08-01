radar
=====

A bunch of these files are one time use scripts, but I included them on this repository anyway. The main file that does the "scene" dection currently is [main_scene_angle.m](main_scene_angle.m). I'll list each important file/directory here with a brief description

### [different_subcarriers](different_subcarriers)
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
  similar to fft_problem directory
#### [check_downsampling.m](check_downsampling.m)
  One time use file
#### [mainDetect.m](mainDetect.m)
  Brian's main detection
#### [mainDetect_once.m](mainDetect_once.m)
  Brians main detection
#### [mainLateral.m](mainLateral.m)
  Brians main lateral
#### [mainLateralBob.m](mainLateralBob.m) 
  My main lateral, used to generate stuff for radarcon
#### [mainLateral_all.m](mainLateral_all.m)
  Brian's stuff, not sure, check his documentation
#### [main_angle.m](main_angle.m)
  My main angle detection file, this file used only a single h0 as opposed to the multiple h0s I've been using recently in the scene detection file
#### [main_scene_angle.m](main_scene_angle.m)
  Rework of main_angle.m to incorporate multiple h0s
#### [test_new_signal.m](test_new_signal.m)
  One time use
#### [test_tx_signals.m](test_tx_signals.m)
  One time use
