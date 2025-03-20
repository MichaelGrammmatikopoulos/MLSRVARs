%======================================================================
%
%          MLSRVARs : MACHINE LEARNING SHADOW RATE VARs                       
%
%          Code Written By: Michael Grammatikopoulos    
%                    email: Michael.Grammatikopoulos@moodys.com  
%                           
%======================================================================

A.) **Toolbox dependencies:** 
  1. Optimization toolbox
  2. Econometrics toolbox
  3. Statistics and Machine Learning toolbox
  4. Parallel computing toolbox
  
B.) The code has by default **_numWorkers = no_of_models;_**.
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

C.) Replication of the paper results can be done with setting **_ndraws = 15000_**.
Computational time in an AMD EPYC 7763 64-Core Processor, 2445 Mhz, 8 Core(s), 16 Logical Processor(s)
is approximately 50 hours. 
 
