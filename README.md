# KIPS
Kiel Particle Simulator

The Kiel Particle Simulator will be an implementation of
H²-matrix data structures and compression algorithms for evaluating
particle interactions (potentials, forces, tensors, etc.) efficiently.

In order to be able to handle a wide range of applications, the code
will be based on tensor interpolation that can be combined with adaptive
algebraic recompression if the configuration of particles is static or
at least mostly static (i.e., when solving integral equations).
