These scripts cover the analysis of simulation output of all non-supplemental figures (and some supplemental figures) of "Intraspecific trait variation in personality‐related movement behavior promotes coexistence" (Milles et al., 2020) 

The scripts have the following purpose:

analysis_tools.R: contains functions to process movement data and toroidal landscapes and to run simulations in NetLogo using R
1_process_trade-off.R: analyses the amount of resources and individual forages per clump and how many clumps and individual utilizes (trade off) which is part of Figure 4
2_plot_trait-foragingeffciency.R: generates the plot of foraging efficiency-trait relationship as it is shown in Fig. 4 including the trade off between exploring few clumps thoroughly and many clumps superficially and exemplary movement paths
3_plot_resource_ratio.R: generates Fig. 5 which shows the change in interspecific differences in competitive abilities between scenarios with and without ITV
4_process_simulationinR.R: runs 1,000 parameter sets generated using latin hypercube sampling from R (files need manually to be adjusted!)
5_process_movement_metrics.R: while the respective simulations of "4_process_simulationinR.R" are run, the movement metrics used for the comparison with empirical data are calculated
6_plot_fit2data.R: generates Fig. 3 which shows the fit of simulated movement metrics to empirically observed variations in personality-related movement behaviour (empirical data can't be provided!)
7_plot_CVA.R: generates Fig. 6 which shows mean time to extinction without ITV and change in mean time to extinction due to ITV