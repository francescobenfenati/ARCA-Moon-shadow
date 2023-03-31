# ARCA cosmic ray Moon and Sun shadows
## Introduction
Cosmic rays are blocked by nearby celestial bodies such as the Moon and the Sun. This induces a deficit in the atmospheric muon flux coming from the direction of these objects. Its observation can be used to verify the pointing accuracy and angular resolution of detectors which are able to measure the secondary muon flux from cosmic ray interactions.
This analysis is based on a comparison between the observed atmospheric muon flux coming from the Moon/Sun directions to the expected flux in case of no Moon/Sun. The search for the shadow of the Moon/Sun is done in the phase space of the angular distancers  between the reconstructed track coordinates and the celestial object. Both a 1-dimensional and a 2-dimensional analysis have been performed. The 1-dimensional analysis uses the space angle between the Moon/Sun and the track. The 2-dimensional analysis uses two Cartesian coordinates: 

$x=(\alpha_{\mu} - \alpha_{MS}) h_{\mu}$

$y=h_{\mu} - h_{MS}$

where $\alpha$ and $h$ are the azimuth and altitude of the tracks and of the Moon/Sun. For the analysis, I used the strategy adopted in[^1][^2].

The significance of the shadowing effect of Moon/Sun is determined with a likelihood ratio test, by comparing the likelihood of a no-shadow hypothesis model *H0*, with a shadow hypothesis model *H1* which includes a shadowing effect. The Poisson likelihood

$\chi^2(H) = 2 \sum_i[ N_{i,H}-n_i +n_i ln(n_i / N_{i,H})]  $

is used, where $n_i$ stands for the event count in the i-th bin to be compared with the expectations $N_{i,H}$ under the *H0* and *H1* hypotheses. A difference in

$\Delta \chi^2_{H1/H0} = \chi^2(H1)−\chi^2(H0)$

values has been used to determine the probability to reject the null hypothesis and to extract the significance of the observation from it.
The Monte Carlo sample is used to define a parametrisation of the track distribution in the absence of a shadowing effect i.e. no-shadow hypothesis
*H0*. 
For the 1D map a constant track density per space angle, $h_0$, describes the data, i.e.

$H0 = h_0$

while the 2D maps can be parametrized with the polynomial function of the form

$p_3(x,y,\mathbf{h},\mathbf{k})= h_0\cdot (1 + h_1x + h_2 x^2 + h_3 x^3)\cdot(1 + k_1 y +k_2 y^2 + k_3 y^3)$

For the 1-dimensional histogram, the event density under hypothesis *H1* is parametrised as follows:

$\frac{dn}{d\delta^2} = k(1 - \frac{A \pi R_{ms}^2}{2 \pi \sigma_{res}^2}e^{\frac{-\delta^2}{2\sigma_{res}^2}}\,)$

where $\delta$ is the angular distance between the track and the Moon/Sun and the second term describes the Gaussian-shaped deficit of events due to the shadowing effects of Moon/Sun.
In the 2D map the *H1* hypothesis is instead defined as

$G(x,y,\mathbf{\theta})=p_3(x,y,\mathbf{h},\mathbf{k})\left(1 - \frac{A \pi R_{MS}^2}{2 \pi \sigma_{res}^2}e^{-\frac{(x-x_{MS})^2+(y-y_{MS})^2}{2\sigma_{res}^2}}\right)$

Here, $A$ is the relative shadow amplitude, $\sigma_{res}$ is the width of the shadow, $R_{MS}$ is the angular radius of the Moon/Sun and $x_{MS}$, $y_{MS}$ are their relative angular position with respect to the nominal one.

The $\sigma_{res}$ is fixed to the value of the angular resolution found in fitting the *H1* hypothesis at the nominal (0,0) Moon/Sun position in the field of view. Free parameters are hence the amplitude
of the Moon/Sun deficit $A$ and the assumed position of the Moon $(x_{MS}, y_{MS})$ in the field of view.
For each bin in the map, the value of the test statistic $\lambda$ is minimised finding
the best estimation of the parameter $A$. The coordinates $(x_{MS}, y_{MS})$ corresponging to the smallest value of $\lambda$ are taken as the best estimation of the Moon/Sun position.



## ARCA8 analysis
The data used in this analysis were collected from ARCA8 (i.e. ARCA detector configuration with 8 Detection Units) between September 26, 2021 and June 1, 2022 for a total of about 8 months of livetime. Datasets of events near the Sun and near the Moon have been merged in order to increase the statistics, obtaining the combined shadow of both celestial objects. Quality cuts on track reconstruction were applied to remove badly reconstructed events. Furthermore, the reconstructed tracks are required to have an elevation of at least 6◦ above horizon. The position of the Moon/Sun in the sky is obtained using ASTROPY package that relies on the International Celestial Reference System (ICRS) coordinates. Monte Carlo simulations are used to optimise the set of quality criteria for track selection, determine a parametrization of the atmospheric muon flux and predict the angular resolution of the cosmic ray deficit induced by the Sun and the Moon
After the application of preliminary cuts required in order to select down-going muons and reduce background signals, events are selected which fall within 12° of angular distance with respect to the Moon or the Sun. An additional more stringent cut on the estimated error on the track angular resolution (&beta;<sub>0</sub>) is applied to improve the quality of reconstructed events.
The cut value is optimized to yield a good balance between the amount of selected events and the quality of the events.

The angular distributions of atmospheric muons from MonteCarlo simulations, before and after applying the &beta;<sub>0</sub> cut is shown:

<img width="481" alt="image" src="https://user-images.githubusercontent.com/48324006/228802053-dde3e928-9b84-46ea-bdc8-584dfcf44ec7.png">

Cuts on other quality reconstruction parameters have been shown to be not effective in improving the angular resolution. 

# Data processing
ARCA data root files are stored in Lyon CC. Each detector run data is saved in a root file, which is then processed by the Data Quality Data Processing team of the Collaboration. After this step, a root files containing run events is produced. 
In order to process those files, select down-going muon track events near the Moon and Sun positions and apply the preliminary cuts one has to launch: 

`./extractor.sh <list_of_root_files> <name_of_output_file> <name_of_logfile>
`

This calls the script `slurm_extractor.sh` which starts a separate job on the Slurm system in Lyon CC for each root file; they in turn call the python script `arca_data_extractor.py`.
(Note: in the python script, change the directory where csv files are written according to the det_id)
Slurm allows for the parallelization of the jobs and the corresponding huge reduction of computing time.
For the MonteCarlo root files instead, use:

`./MC_extractor.sh <list_of_root_files> <name_of_output_file> <name_of_logfile>`

(e.g. `./MC_extractor.sh ARCA8/MC/arca8rbr_mu_v8_list_subsample.txt /sps/km3net/users/fbenfe/Moon_shadow/CNN/arca8rbr_mu_v8.txt /sps/km3net/users/fbenfe/Moon_shadow/logs/arca8/MC/ `)

which calls the script `MC_extract_and_cut.py`
The data processing is largely based on the usage of Pandas dataframes.
The outputs are merged csv files containing a row for each surviving track event. Each row contains a number of columns corresponding to a subsample of parameters of interest for that event.

### 1-dimensional and 2-dimensional histograms
The steps to produce the 1-dimensional histograms for the Moon and the Sun shadows are the following:
* Prepare list of csv files and plot names:
e.g.`ARCA8/data/v8.0/1d_list.txt` :
`
/sps/km3net/users/fbenfe/Moon_shadow/csv/arca21/data/arca21_grbv8_pre+anoise_sun.txt arca21_grbv8_0_sun
/sps/km3net/users/fbenfe/Moon_shadow/csv/arca21/data/arca21_grbv8_pre+anoise_moon.txt arca21_grbv8_0_moon
`

* launch the script:
`./launcher_1D.sh <list of csv> <beta_cut> <logfile_path>`

where the second argument is the &beta;<sub>0</sub> cut value you want to apply to your data, e.g.:

`./launcher_1D.sh ARCA8/data/v8.0/1d_list.txt 0.5 logs/arca8 `

The bash script executes `data_density_histogram.exe` which is obtained compiling the root program `data_density_histogram.cc`.
This will produce the 1-dimensional histogram of the event densities in the region of angular distance < 4° from the Moon/Sun, for each of the csv data files in the list. Fit parameters are logged onto a file, or they can be printed as standard output.

For the 2-dimensional histogram instead. use 

`./launcher_2D.sh <list of csv> <sigma_fit_start> <beta_cut> <maxdeg> <bin_size> <logfile_path>`

where 

* `<sigma_fit_start>` is the starting value of the sigma parameter in the Gaussian fit (this is usually taken equal to the sigma fit parameter obtained from the 1-dimensional fit)
* `<beta_cut>` is the &beta;<sub>0</sub>; cut value you want to apply to your data
* `<maxdeg>` is the maximum degree of the 2D map
* `<bin_size>` is the bin size of the histogram

e.g.
`./launcher_2D.sh ARCA8/data/v8.0/2d_list.txt 0.5 0.5 6 0.1 logs/arca8`

The bash script executes `contour_data.exe` which is obtained compiling the root program `arca_2d_data_analysis.cc`.
This will produce the 2-dimensional map of $\Delta \chi^2_{H1/H0}$ in the [x,y] plane, for each of the csv data files in the list. Fit parameters are logged onto a file located in the specified directory.

The 1D and 2D plots obtained for ARCA8 dataset with a selected cut on &beta;<sub>0</sub> = 0.3 are

<img width="475" alt="image" src="https://user-images.githubusercontent.com/48324006/228802699-fa35a0c5-50ed-407f-8d60-bfc17463fe31.png">

![image](https://user-images.githubusercontent.com/48324006/228802771-81afe80e-ef2f-4b1c-a3b2-36c02cab34c5.png)



## Neural network classification
A classifier neural network has been developed in order to improve the quality of the selected tracks.
The code of the network, which is a 5-layers multiperceptron, is available on Google Colab at this link: 
https://colab.research.google.com/drive/1wO1DhXK8ba4KpiFo78aWYm-gvUZitjCe#scrollTo=qmW3rLz4plU_

The network classifies events into two classes: "good" reconstructed events with angular resolution < 1° and "bad" reconstructed events with angular resolution > 1°. 

![image](https://user-images.githubusercontent.com/48324006/229098088-62cc63d4-04f8-48a2-bafa-d95e771f897e.png)

After choosing an output score value that balances the false positive contamination and the amount of surviving events, the angular distribution looks like

![image](https://user-images.githubusercontent.com/48324006/229099478-71a552e7-1008-4f11-ab3d-18ca322dc8af.png)

The new dataset which contains tracks selected by applying the NN shows a significant improvement in the angular resolution distribution which now looks peaked and not flat anymore. The median angle is 0.49°.

NN predictions on test sample and on real data have been compared showing consistency:

<img width="725" alt="image" src="https://user-images.githubusercontent.com/48324006/229106333-eb21a009-191b-46f9-afea-47473249e773.png">


Also, in the 1D and 2D plots the presence of a signal is more evident, although the statistics is still too low to observe a clear signal (Event densities of ~1500 /deg^2 are estimated to be required for a 3sigma  signal)

<img width="475" alt="image" src="https://user-images.githubusercontent.com/48324006/228802637-6708a10e-3116-4691-9c7f-00b92cf4179a.png">


![image](https://user-images.githubusercontent.com/48324006/228802543-1c52a40f-1bf4-4a61-9517-3baad40898f1.png)

##Outlook
The implementation of all thse above-mentioned scripts in one single routine to be used also by non-experts users is foreseen. This routine could be exploited to perform regular analysis (for example every N months) on the detector pointing accuracy. 
The analysis on ARCA8 could be also integrated with future analysis with the growing detector, and also with the foreseen introduction of dynamic calibrations that take into account the DOMs horizontal fluctuations along the Detection Units due to the sea currents.


### Other useful scripts
If you need to separate runs from a merged dst file (too big to be analyzed as one) use:

`python dst_merged_separator.py`

To analyze then the single h5 dataframes:

`arca_h5_extractor.py <input_file> <output_file>`



Plot 1D hist. for ORCA Sun and Moon:
1) Prepare list of csv and plot names:

e.g. ORCA6/list_orca.txt :
`/sps/km3net/users/fbenfe/Moon_shadow/csv/ORCA/orca_data_pre+anoise+sun.txt orca_v7_0_sun
/sps/km3net/users/fbenfe/Moon_shadow/csv/ORCA/orca_data_pre+anoise+moon.txt orca_v7_0_moon
`
3) Launch:

`./orca_launcher_1D.sh <list of csv> <beta_cut>`
e.g.
`./orca_launcher_1D.sh ORCA6/list_orca.txt 0.27`


`scan_beta.cc`

creates csv file with lambda, sigma, amp scanning for beta cuts
otuput ex: "beta_scan_12s_shad_loglik.txt"

`scan_beta_plot.cc`

plotting variables as a function of beta cuts

COMPILE ROOT SCRIPT:

```
g++ data_density_histogram.cc -o data_density_histogram.exe `root-config --cflags --glibs`
g++ arca_2d_data_analysis.cc -o contour_data.exe `root-config --cflags --glibs`
```

### References
[^1]: KM3NeT Collaboration, [arXiv:2211.08977 [astro-
ph.IM]]
[^2]: Antares Collaboration, Eur.Phys.J.C 78 (2018) 12, 1006 
