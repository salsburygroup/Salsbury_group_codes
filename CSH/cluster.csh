#!/bin/tcsh

set total_idle = 0
set total_cpus = 0

foreach line (`sinfo -N -o "%N %C" | tail -n +2`)
  set idle_cpus = `echo $line | awk -F'/' '{print $2}'`
  set total_node_cpus = `echo $line | awk -F'/' '{print $4}'`

  if ("$idle_cpus" =~ [0-9]* && "$total_node_cpus" =~ [0-9]*) then
    @ total_idle += $idle_cpus
    @ total_cpus += $total_node_cpus
  endif
end

echo "Total CPUs: $total_cpus"
echo "Total idle CPUs: $total_idle"

