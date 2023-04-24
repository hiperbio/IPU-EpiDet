# IPU-EpiDet

*IPU-EpiDet* is a tool for performing exhaustive high-order epistasis detection searches on Graphcore's IPU-POD systems.
It represents an implementation of algorithms devised around the bulk synchronous parallel model of execution used by the IPU, where tasks are split into supersteps composed of computation, communication and synchronization phases.
The tool exploits parallelism at multiple levels, making efficient use of the general purpose computation capabilities of multiple Intelligence Processing Unit (IPU) processors, each capable of executing thousands of hardware threads.
While the tool targets and has only been extensively tested on second generation IPUs (MK2/GC200 microarchitecture), it should also support systems with IPU devices implementing other microarchitectures.

## Dependencies and Requirements

* Poplar SDK (tested on v.2.6); and Python for dataset generation / Gnuplot for plotting charts
* Graphcore IPU-POD system with MK2/GC200 IPUs (tested on IPU-POD64)

## Installation and Deployment Process

The tool can be compiled targeting a range of different IPU-POD configurations with a single execution of the `make` command; or to a specific target, e.g. `make pod64` to compile a binary targeting an IPU-POD64 system.

On a system without IPU devices, the tool can still be used through execution of a binary (`ipu-epidet_cpu`) targeting a virtual (i.e. emulated) IPU.
Notice, however, that performance will be orders of magnitude lower than that achievable on real IPU-POD systems.

## Usage Example and Scalability Test

Epistasis detection searches are performed passing a dataset as input to the tool, as follows:

```bash
$ ./bin/ipu-epidet_pod4 datasets/6048snps_4000samples.csv
```

Execution of the `run_small.sh` script performs a scalability test targeting an IPU-POD64, making use of 1, 2, 4, 8, or all its 16 IPU-M2000 compute engines (each with 4 IPUs) to process datasets with 6048 SNPs and 4000 or 8000 samples.
Alternatively, executing `run_medium.sh` or `run_large.sh` results in runs being performed with more SNPs (medium: 12096, large: 24192).
In addition to encapsulating commands related to compilation and execution of the tool, these scripts also abstract the generation of synthetic datasets with the required dimensions and the generation of charts for visualization.

