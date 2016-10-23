# Parallel MultiNest with OpenMP
MultiNest is a maximum likelihood estimator. Here are several versions which 
make use of multiple kernels without changing the outcome regarding the vanilla 
sampling approach.

## What is inside:
Here you can find different versions of MultiNest. The original version has been
downloaded from http://ccpforge.cse.rl.ac.uk/gf/project/multinest/frs/?action=FrsReleaseBrowse&frs_package_id=83
(Login required) and has been modified slightly to make it more readable (e.g
corrected the indentation, added some comments, deleted unused toy models and
functions).
The other versions are described in the Bachelorthesis **Scalability of Sampling
in MultiNest** from January, 26 2016.
The parallelization approaches use MultiNest beginning with 2 threads up to
the available number of processors via `OMP_get_num_procs()`. If you want to
change that you have to alter `main.f90` within each toy model folder.
The evaluated time will be saved in a file `benchmark_*` at folder 
**benchmark**.


## Usage:
The parameters for each toy model are set to the parameters used in the
Bachelorthesis but you may modify them by editing the file `params.f90` in each 
toy model folder. 
In the Bachelorthesis there is a call to `sleep()` used. This has been commented
in all versions. If you want to use this call you have to edit the file 
`like.f90` in each toy model folder.
For compiling please use:
``make   <-- For building MultiNest 
make eggbox
make gauss_shell
make himmelblau``

For cleaning please use:
``make clean  <-- For cleaning MultiNest
make clean_eggbox
make clean_gauss_shell
make clean_himmelblau``

For creating evidence files you need to edit `params.f90` for each toy model
and set `nest_outfile` to true. The evidence files will be saved in `chains`.
For recreating the pictures from the Bachelorthesis you need to use the files
ending with `*ev.dat` together with either `visualize2D.py` or 
`visualize3D.py` by putting the files in one folder and by typing
``python visualize2D.py data_name picture_name [Title of graph]``
or
``python visualize3D.py data_name picture_name [Title of graph]``
into your terminal.


## Additional information:
You need to install the LAPACK (Linear Algebra PACKage) library before using
MultiNest.
Inside the folder with the original code is another Readme with more information 
from the original authors Farhan Feroz and Mike Hobson.


## Is it any good?
A speedup depends on the likelihood region and the used approach. Just don't use 
the second approach as it does not always terminate. 
