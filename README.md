# Preliminary studies in the search of a Dark Matter signal with LIME data
When no significant **WIMP-induced nuclear recoils** are observed above the expected background, it is only possible to set **exclusion limits** in the **cross-section versus WIMP mass** parameter space.  
To compute an exclusion limit, one must determine, for a given WIMP mass, the cross section corresponding to a given number of expected events.  
The experiment consists of **counting the number of detected events**. In this preliminary study, neither the **energy spectra** nor the **directional information** of the observed events are exploited.  
By estimating the **upper limit on the number of WIMP-induced events**, denoted as $\mu_s$, from the **90% confidence interval (C.I.)**, it becomes possible to calculate the upper limit on the **WIMP–nucleon elastic cross section** for both the **Spin-Independent (SI**) and **Spin-Dependent (SD)** cases.  
The expected number of Dark Matter events is defined as: ![n_evt](https://latex.codecogs.com/svg.image?&space;N_{DM_{evt}}=\mu_{s,90\%}). 

## Limit evaluation

This repository contains code to evaluate limits for both **Spin-Dependent (SD)** and **Spin-Independent (SI)** interactions.  
The expected number of events is computed using the following code: https://github.com/antorita/cygno_sensitivity/tree/main

### Required Libraries 
The following Python libraries are required:
- numpy  
- math  
- scipy.integrate  
- sys

### How to run the code
Run the main script using Python3:

```
python3 lim.py <option>
```
where &lt;option&gt; can be:  
SD -> evaluate the Spin-Dependent case  
SI -> evaluate the Spin-Independent case  

After running the script, a new file named lim_output.txt will be created containing all the computed values.

### How to visualize the results
To generate plots from the output, run:
```
python3 plot.py
```

### Configuration Parameters (lim.py)
The following variables **must be set by the user in lim.py**:
- Energy threshold [ee]: [LowThre](https://github.com/antorita/limit_evaluation/blob/944f0930cbc719cb323307b12198d9d35f2bb056/lim.py#L41)
- correct exposuretime [days]:: [Daqtime](https://github.com/antorita/limit_evaluation/blob/944f0930cbc719cb323307b12198d9d35f2bb056/lim.py#L44)
- working pressure: [workingpressure](https://github.com/antorita/limit_evaluation/blob/944f0930cbc719cb323307b12198d9d35f2bb056/lim.py#L45)
- sensitive volume: [volume](https://github.com/antorita/limit_evaluation/blob/944f0930cbc719cb323307b12198d9d35f2bb056/lim.py#L46)
- Expected number of Dark Matter events: [Nev](https://github.com/antorita/limit_evaluation/blob/944f0930cbc719cb323307b12198d9d35f2bb056/lim.py#L75)

## References
For more details, see **Chapter 7 of my PhD thesis** https://arxiv.org/abs/2510.01646
