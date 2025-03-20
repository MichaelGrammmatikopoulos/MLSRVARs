The code has by default no_of_workers == no_of_models.
If you get this error: 

_Error using parpool (line 133)
Too many workers requested. The cluster "Processes"
has the NumWorkers property set to a maximum of 6
workers but 32 workers were requested. Either
request a number of workers less than NumWorkers, or
increase the value of the NumWorkers property for
the cluster (up to a maximum of 512 for the Local
cluster)._

go to HOME>parallel>preferences and change the number of workers.
