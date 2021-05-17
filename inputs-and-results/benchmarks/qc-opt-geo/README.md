# QC Optimized Geometry Benchmarks

This directory contains the scripts needed to compute the ddE, RMSD, and TDF metrics for a set
of force fields against a QCArchive optimization data collection such as the `OpenFF Full Optimization Benchmark 1` 
QCArchive data set as was used in the [Lim and Hahn study](https://doi.org/10.26434/chemrxiv.12551867.v2) (see the 
Acknowledgments section below for more details).

## Running the benchmarks

1) Run the `01-setup.py` script to retrieve the relevant QC data from the QCArchive and filter out
   and incomplete or undesirable (e.g. records whereby the molecule connectivity changed during 
   the QM minimization) records. The dataset to retrieve and benchmark against must be specified.
   For the Lim and Hahn study:
   
   ```shell
   python 01-setup.py "OpenFF Full Optimization Benchmark 1"
   ```
   
   or for the Industrial partner benchmarks:

   ```shell
   python 01-setup.py "OpenFF Industry Benchmark Season 1 v1.0"
   ```

2) 

    a. Run the `02-a-chunk-qm.py` script to split the `01-processed-qm.sdf` file produced by step 1) into
       many smaller chunks. This allows the next step to be parallelized across many different CPUs / workers.
       By default the split chuncks will be stored in a new `02-chunks` directory.

    b. Perform an MM energy minimization for each of the chunked files and for each force field of interest, e.g.:

   ```shell
   python 02-b-minimize.py -ff openff-1.3.0.offxml \ 
                           -i 02-chunks/01-processed-qm-1.sdf \
                           -o 02-outputs/openff-1-3-0-1.sdf
   ```

    c. Concatenate the files produced by step 2b) into single SDF files:

   ```shell
   python 02-c-join-outputs.py -i openff-1-3-0 \ 
                               -dir 02-chunks \
                               -o 02-outputs/openff-1-3-0.sdf
   ```

3) Run `03-compare-ffs.py` to compute the ddE, RMSD, and TDF for the minimized structures. These will
   be stored in a `03-metrics.pkl` file.
   
4) Plot the metrics by running the `04-plot-metrics.py` script.

## Acknowledgments

The scripts in this directory are based off of those found in the `benchmarkff` repository
at commit hash [`6351878`](https://github.com/MobleyLab/benchmarkff/tree/6351878) which was produced 
by Victoria Lim under the following license:

    MIT License
    
    Copyright (c) 2019 Victoria Lim
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

