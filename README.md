## BDFFT: software for finite-frequency traveltime tomography
BDFFT contains a set of fortran and c-shell written programs used to conduct (multiscale) finite-frequency traveltime tomography beneath the region of interest using cross-correlation measured P or S traveltime residuals from teleseismic earthquakes.  

Installation and tutorials are described in `readme_BDFFT`.

Authors: [`Shu-Huei Hung`] at Dept of Geosciences, National Taiwan University (codes of finite-frequency kernel computation) & [`Lin-Yun Chiao`] at Inst. of Oceanography, National Taiwan University (codes of multiscale wavelet-based inversion)

### REFERECES:
#### For the theory behind calculating the ray-theoretical finite-frequency traveltime kernels:
- Dahlen, F. A., Hung, S.-H., & Nolet, G. (2000). Fréchet kernels for finite-frequency traveltimes-I. Theory. Geophysical Journal International, 141(1), 157–174. https://doi.org/10.1046/j.1365-246x.2000.00070.x
- Hung, S.-H., Dahlen, F. A., & Nolet, G. (2000). Fréchet kernels for finite-frequency traveltimes-II. Examples. Geophysical Journal International, 141(1), 175–203. https://doi.org/10.1046/j.1365-246x.2000.00072.x

#### For multiscale model parameterization/regularization behind the tomographic inversion algorithm:
Chiao, L.-Y., & Kuo, B.-Y. (2001). Multiscale seismic tomography. Geophysical Journal International, 145(2), 517–527. https://doi.org/10.1046/j.0956-540x.2001.01403.x
- Chiao, L.-Y., & Liang, W.-T. (2003). Multiresolution parameterization for geophysical inverse problems. Geophysics, 68(1), 199–209. https://doi.org/10.1190/1.1543207

##### For applications of (multiscale) finite-frequency tomography:
- Hung, S., Shen, Y., & Chiao, L. (2004). Imaging seismic velocity structure beneath the Iceland hot spot: A finite frequency approach. Journal of Geophysical Research: Solid Earth (1978–2012), 109(B8). https://doi.org/10.1029/2003jb002889
- Hung, S., Chen, W., Chiao, L., & Tseng, T. (2010). First multi-scale, finite-frequency tomography illuminates 3-D anatomy of the Tibetan Plateau. Geophysical Research Letters, 37(6), n/a-n/a. https://doi.org/10.1029/2009gl041875
- Hung, S.-H., Chen, W.-P., & Chiao, L.-Y. (2011). A data-adaptive, multiscale approach of finite-frequency, traveltime tomography with special reference toPandSwave data from central Tibet. Journal of Geophysical Research, 116(B6). https://doi.org/10.1029/2010jb008190
