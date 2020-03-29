# Inverse Electro-Hydrodynamic Lithography

## Problem statement

Consider spatiotemporal state <img src="/tex/7500c95c8fa7036167bd50777fdb3bfb.svg?invert_in_darkmode&sanitize=true" align=middle width=78.73854614999999pt height=24.65753399999998pt/> which statisfies 
the evolution PDE

<p align="center"><img src="/tex/c731cb4694f783196dafbe68c542f6af.svg?invert_in_darkmode&sanitize=true" align=middle width=242.68094564999998pt height=39.452455349999994pt/></p>

Given state initial condition at <img src="/tex/7102f42feef69c607ecb407d204abc19.svg?invert_in_darkmode&sanitize=true" align=middle width=42.02615174999998pt height=22.465723500000017pt/> and 
target profile at terimal time <img src="/tex/6dc238b06d1a5e1c37fcea4fd529f24a.svg?invert_in_darkmode&sanitize=true" align=middle width=45.696254999999994pt height=22.465723500000017pt/> such as

<p align="center"><img src="/tex/3b7e4dc7d0ecfa7fbed6198e7dc89d80.svg?invert_in_darkmode&sanitize=true" align=middle width=314.7540924pt height=16.438356pt/></p>

find the optimal design of spatial control <img src="/tex/737ffbca1c20b24c88b098ee6f893e74.svg?invert_in_darkmode&sanitize=true" align=middle width=61.34934959999998pt height=24.65753399999998pt/> such that

<p align="center"><img src="/tex/be8625cc9e16e4d1e4b0194c47557cf3.svg?invert_in_darkmode&sanitize=true" align=middle width=172.9848912pt height=16.438356pt/></p>

or at least as close as possible.

## Dependencies

* Eigen
* OpenCV 3.2
* Revised INIReader from https://github.com/benhoyt/inih
* C++11

