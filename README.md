# limit_evaluation

This repository contains code to evaluate limits for spin-dependent (SD) and spin-independent (SI) cases.

## How to run the code
Run the main script using Python3:

'''bash
python3 lim.py <option>

where option can be: 
SD -> evaluate the Spin-Dependent case
SI -> evaluate the Spin-Independent case

After running the script, a new file named lim_output.txt will be created containing all the computed values.

## How to visualize the results

To generate plots from the output, run:

python3 plot.py
