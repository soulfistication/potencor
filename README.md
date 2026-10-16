# POTENCOR MATLAB

Original program by Miguel Angel Guevara P. written in Pascal then Delphi. Computes power and correlation parameters on EEG signals
using 16 bit data.

Change lists:

* Converted code to C++ and Python
* Converted code to Matlab by replacing the custom FFT implementation (custom impl was nice but made the code quite unreadable and
  2026 FFT is standard in many programming languages), sin & cos precomputed functions and the Fisher-z transform
* Converted code from Object Oriented (useless in Matlab) to procedural

Matlab built in functions should run much quicker and be fail proof to the original Pascal fft implementation and have index
array bounds checks (right?).
