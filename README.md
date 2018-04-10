# Ember Communication Pattern Library

Multi-node communication patterns underpin the scalability and parallel
performance of the Department of Energy, and broader HPC workloads.
Modeling of these patterns is as important aspect of extreme scaled
supercomputing systems. To date, many vendors have relied on
communication traces which can be difficult to obtain at scale, and take
significant I/O storage. For interconnect simulators, the reading and
replay of traces requires high-performance I/O subsystems which are
often expensive and may be unavailable. To this end, the Ember suite
provides communication patterns in a simplified setting (simplified by
the removal of application calculations, control flow etc.). This
enables more efficient traces to be captured, or in the cases of the
Structural Simulation Toolkit (SST, http://sst-simulator.org),
these patterns can be easily
replicated without tracing using the Ember/SST motif library. The
intention of Ember is to enable much larger-scale modeling of
high-performance interconnects to achieve DOE's goal of scalable
Exascale computing systems. The motifs contained in the suite are
intentionally simplified, and by design, do not capture every
permutation of the basic patterns within the DOE workload. When used
collectively, our experience working with leading industry vendors has
been that the motifs capture pertinent aspects of the network
interconnect.

Communication Patterns information:

* README.MPI.halo3d (Structured nearest neighbor-like)
* README.MPI.halo3d-26 (Unstructured nearest neighbor-like)
* README.MPI.incast (Multiple inbound messages, I/O-like)
* README.MPI.sweep3d (Communication sweeping)
* README.SHMEM.randominc (Uniform random network access)
* README.SHMEM.hotspotinc (Hotspot random network access)

The Ember Communication Pattern Library is developed by the Scalable
Computer Architectures group at Sandia National Laboratories, NM.
Funding for the development is provided by the DOE NNSA/ASC Computing
Program and the DOE's Exascale Computing Project Hardware Evaluation
(HE) team. Information relating to the use of this code can be found in
the LICENSE file or individual source code files.

For more information please contact: Simon Hammond (Sandia National
Laboratories, sdhammo@sandia.gov).
