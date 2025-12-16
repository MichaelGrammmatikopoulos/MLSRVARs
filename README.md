------------------------------------------------------------------------
MLSRVARs : MACHINE LEARNING SHADOW RATE VARs   
------------------------------------------------------------------------
Grammatikopoulos, M. 2025. "Forecasting With Machine Learning 
Shadow-Rate VARs." Journal of Forecasting 1–17. 
https://doi.org/10.1002/for.70041. 
------------------------------------------------------------------------
This code comes without technical support of any kind.  It is expected 
to reproduce the results reported in the paper. Under no circumstances 
will the author be held responsible for any use (or misuse) of this code 
in any way.
------------------------------------------------------------------------
The views in this paper are solely of the author and do not represent 
the views of Moody's Analytics or the Moody’s Corporation
------------------------------------------------------------------------

------------------------------------------------------------------------
A.) **Toolbox dependencies** 
------------------------------------------------------------------------

  1. Optimization toolbox
  2. Econometrics toolbox
  3. Statistics and Machine Learning toolbox
  4. Parallel computing toolbox

------------------------------------------------------------------------
B.) **Set up** 
------------------------------------------------------------------------
You should only have to edit line 60 with your directory.
There might be additional steps needed in case you do not have the above 
toolboxes installed, or installed in different locations.

------------------------------------------------------------------------
C.) **Parallel computing settings**
------------------------------------------------------------------------

The code has by default **_numWorkers = no_of_models;_**.
If you get this error: 

**_Error using parpool (line 133)
Too many workers requested. The cluster "Processes"
has the NumWorkers property set to a maximum of 6
workers but 32 workers were requested. Either
request a number of workers less than NumWorkers, or
increase the value of the NumWorkers property for
the cluster (up to a maximum of 512 for the Local
cluster)._**

go to HOME>parallel>preferences and change the number of workers. 
Another solution is to comment out lines 221-222.

------------------------------------------------------------------------
D.) **Reproductibility / Computation time**
------------------------------------------------------------------------

Replication of the paper results can be done with setting **_ndraws = 15000_**.
Computational time in an AMD EPYC 7763 64-Core Processor, 2445 Mhz, 8 Core(s), 16 Logical Processor(s)
is approximately 50 hours. 
