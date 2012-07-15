============================================================================

 Description   This archive contains code to replicate the results in
	       "Tapping the Supercomputing Under Your Desk: Solving
	       Dynamic Equilibrium Models with Graphics Processors" by
	       Aldrich, Fernandez-Villaverde, Gallant and
	       Rubio-Ramirez.

 Development   The CUDA code was developed using version 2.2 of the
	       NVIDIA CUDA Development Tools on a DELL Precision
	       Workstation R5400 with two 2.66 GHz quad core Intel
	       Xeon CPUs and one NVIDIA GeForce GTX 280 GPU, running a
	       Linux operating system. Compatibility with other
	       operating systems and other versions of the CUDA
	       development tools is not guaranteed.

 Details       The GPU timing results in the paper were not obtained
	       with a fixed block size - the optimal block size varies
	       with the number of total grid points allocated for
	       capital.

 Use           To run the code, simply expland the TAR file on a
	       Linux/Unix machine by typing 'tar -xzvf vfi.tar' at a
	       shell prompt and then type 'make; vfi'. Timing results
	       are printed to the display, as wells as comparisons of
	       the value functions and policy functions computed on
	       the CPU and GPU.

 =============================================================================

 Author:       Eric M. Aldrich

 Contact:      ealdrich@gmail.com

 Date:         28 July 2011

 ============================================================================
