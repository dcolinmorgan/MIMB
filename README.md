# MIMB

Chapter submission for MIMB Summer 2021


## Outline:
- Brief literature review
- Data formats -- what is required, depending on availability
  - Chr start stop
  - Cg annotation (illumina)
  - Other → some method for relating genome location to binding
- Pipeline:
  0. FIMO scan & open-access data
  1. (py)bedtools
  2. integrating
    1. Integrating into netzoo bipartite
    1. Integrating into other GRN framework
      1. Smaller GRN models (TF-gene subsets) could just merge with methyl-motif and reduce their estimates where overlap / downweight (GRN estimate x 1- avg meth ratio)
  3. Benchmarking against ChIP-seq