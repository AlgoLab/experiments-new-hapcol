tables of runtimes, memory, and in the case of a haplotyping method,
switch error

here we analyze the resulting datasets obatined from 3 pipelines :

1. downsampling by whatshap's (deterministic) heuristic
2. downsampling by a (pseudo-random) shuffle and greedy selection
3. merging using the red-blue graph, followed by (2.)

. obtain_wif.tab is the time/memory needed for (1.)

. obtain_random_downsampling.tab is the time/memory needed for the
  downsampling part of (2.) and (3.)

. obtain_merging.tab is the time/memory neede for the first (the
  merging) part of (3.)

. hapchat_whatshap_downsample.tab, hapchat_random_downsample.tab and
  hapchat_merge_random_downsample.tab are then the time/memory/swerr
  incurred when running hapchat on the respective datasets (1.) (2.)
  and (3.)

. whatshap.tab only concerns (1.), and is the time/memory/swerr
  incurred on this dataset
