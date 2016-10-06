# Spherical-Array-Processing
A collection of MATLAB routines for acoustical array processing on spherical harmonic signals, commonly captured with a spherical microphone array.

---
>
>   Archontis Politis, 2014 
>
>   Department of Signal Processing and Acoustics, Aalto University, Finland 
>   archontis.politis@aalto.fi 
>
---

## Description <a name="description"></a>

This is a collection MATLAB routines that perform array processing
techniques on spatially transformed signals, commonly captured with 
a spherical microphone array. The routines fall into four main
categories:

1. obtain spherical harmonic (SH) signals with broadband characteristics,
as much as possible,

2. generate beamforming weights in the spherical harmonic domain (SHD)
for common signal-independent beampatterns,

3. demonstrate some adaptive beamforming methods in the SHD,

4. demonstrate some direction-of-arrival (DoA) estimation methods in the
SHD,

5. demonstrate methods for analysis of diffuseness in the sound-field

6. demonstrate flexible diffuse-field coherence modeling of arrays

The latest version of the library can be found at

  https://github.com/polarch/Spherical-Array-Processing

Detailed demonstration of the routines is given in TEST_SCRIPTS.m and at

  http://research.spa.aalto.fi/projects/spharrayproc-lib/spharrayproc.html

The library relies in the other two libraries of the author related to
acoustical array processing found at:

  https://github.com/polarch/Array-Response-Simulator
  https://github.com/polarch/Spherical-Harmonic-Transform

They need to be added to the MATLAB path for most functions to
work.

For any questions, comments, corrections, or general feedback, please
contact archontis.politis@aalto.fi

---
## Table of Contents

  * [Description](#description)
  * [Microphone Signals to Spherical Harmonic signals [refs 1-4]](#mic2sh)
  * [Signal-independent Beamforming in the Spherical Harmonic Domain [refs 5-10]](#sibeam)
  * [Signal-Dependent and Parametric Beamforming [refs 11-12]](#sdbeam)
  * [Direction-of-Arrival (DoA) Estimation in the SHD](#doa)
  * [Diffuseness and Direct-to-diffuse Ratio (DDR) Estimation [refs 13-17]](#diffuseness)  
  * [Diffuse-field Coherence of Directional Sensors/beamformers [refs 8,18]](#dfc)
  * [References](#references)

---
## Microphone Signals to Spherical Harmonic signals [refs 1-4] <a name="mic2sh"></a>

The first operation is to obtain the SH signals from the microphone
signals. That corresponds to two operations: a matrixing of the signals
that performs a discrete spherical harmonic transform (SHT) on the
sound pressure over the spherical array, followed by an equalization step
of the SH signals that extrapolates them from the finite array radius to
array-independent sound-field coefficients. This operation is limited by
physical considerations, and the inversion should be limited to avoid
excessive noise amplification in the SH signals. A few approaches are
included:

* Theoretical-array-based filters 
* Measurement-based filters

The measurement-based approach is demonstrated using directional responses of an Eigenmike array [ref2], measured in the anechoic chamber of the Acoustics laboratory, Aalto University, Finland. The measurements were conducted by the author and Lauri Mela, summer 2013.

---
## Signal-independent Beamforming in the Spherical Harmonic Domain [refs 5-10] <a name="sibeam"></a>

After the SH signals have been obtained, it is possible to perform
beamforming on the SHD. In the frequency band that the SH signals are
close to the ideal ones, beamforming is frequency-independent and it
corresponds to a weight-and-sum operation on the SH signals. The
beamforming weights can be derived analytically for various common
beampatterns and for the available order of the SH signals. Beampatterns
maintain their directivity for all directions in the SHD, and if they are
axisymmetric their rotation to an arbitrary direction becomes very
simple.

The following axisymmetric patterns are included in the library:

* cardioid [single null at opposite of the look-direction]
* supercardioid (up to 4th-order) [ref5]
[maximum front-to-back rejection ratio]
* hypercardioid/superdirective/regular/plane-wave decomposition beamformer
[maximum directivity factor]
* max-energy vector (almost super-cardioid) [ref6]
[maximum intensity vector under isotropic diffuse conditions]
* Dolph-Chebyshev  [ref7]
[sidelobe level control]
* arbitrary patterns of differential form [ref8]
[weighted cosine power series]
* patterns corresponding to real- and symmetrically-weighted linear array [ref9]

and some non-axisymmetric cases:

* closest beamformer to a given directional function
[best least-squares approximation]
* acoustic velocity beamformers for a given spatial filter [ref10]
[capture the acoustic velocity of a directionally-weighted soundfield]

---
## Signal-Dependent and Parametric Beamforming [refs 11-12] <a name="sdbeam"></a>

Contrary to the fixed beamformers of the previous section, parametric and
signal-dependent beamformers use information about the signals of
interest, given either in terms of acoustical parameters, such as
DoAs of the sources, or extracted through the second-order statistics of
the array signals given through their spatial covariance matrix (or a
combination of the two)

The following examples are included in the library:

* plane-wave decomposition (PWD) beamformer at desired DoA, with
nulls at other specified DoAs
* null-steering beamformer at specified DoAs with constraint on 
omnidirectional response
* minimum-variance distortioneless response (MVDR) in the SHD
* linearly-constrained minimum variance (LCMV) beamformer in the SHD
* informed parametric multi-wave multi-channel Wiener spatial filter
(iPMMW) in the SHD [ref12]

---
## Direction-of-Arrival (DoA) Estimation in the SHD <a name="doa"></a>

Direction of arrival (DoA) estimation can be done by a steered-response 
power approach, steering a beamformer on a grid and checking for peaks on
the power output, or by a subspace approach such as MUSIC. Another
alternative is to utilize the acoustic intensity vector, obtained from
the first-order signals, which its temporal and spatial statistics reveal
information about presence and distribution of sound sources.

A few examples of DoA estimation in the SHD are included in the library:

* Steered-response power DoA estimation, based on plane-wave 
decomposition (regular) beamforming
* Steered-response power DoA estimation, based on MVDR beamforming
* Acoustic intensity vector DoA estimation
* MUSIC DoA estimation in the SHD

---
## Diffuseness and Direct-to-diffuse Ratio (DDR) Estimation [refs 13-17] <a name = "diffuseness"></a>

Diffuseness is a measure of how close a sound-field represents ideal
diffuse field conditions, meaning a sound-field of plane waves with random 
amplitudes and incident from random directions, but with constant mean
energy density, or equivalently constant power distribution from all
directions (isotropy). There are measures of diffuseness that consider
point to point quantities, here however we focus on measures that
consider the directional distribution, and relate to the SHD. In the case
of a single source of power $P_s$ and an ideal diffuse field of power $P_d$, 
diffuseness $\psi$ is directly related to the direct-to-diffuse ratio (DDR) 
$\Gamma = P_s/P_d$, with the relation $\psi = P_d/(Ps+Pd) = 1/(1+\Gamma)$. 
In the case of multiple sources $P_i$, an ideal diffuseness can be defined 
as $\psi = P_d/(P_d + \sum_i P_i)$ and DDR as $\Gamma = \sum_i P_i/P_d$.
Diffuseness is useful in a number of tasks - for acoustic analysis, for 
parametric spatial sound rendering of room impulse responses (e.g SIRR) 
or spatial sound recordings (e.g. DirAC), or for constructing filters for 
diffuse sound suppresion.

The following diffuseness measures are implemented

* intensity-energy density ratio (IE) [ref.13]
* temporal variation of intensity vectors (TV) [ref.14]
* spherical variance of intensity DoAs (SV) [ref.15]
* directional power variance (DPV) [ref.16]
* COMEDIE estimator (CMD) [ref.17]

---
## Diffuse-field Coherence of Directional Sensors/beamformers [refs 8,18] <a name = "dfc"></a>

The diffuse-field coherence (DFC) matrix, under isotropic diffuse conditions,
is a fundamental quantity in acoustical array processing, since it models
approximately the second-order statistics of late reverberant sound between array signals, and 
is useful in a wide variety of beamforming and dereverberation tasks. 
The DFC matrix expresses the PSD matrix between the array sensors, or 
beamformers, for a diffuse-sound field, normalized with the diffuse sound 
power, and it depends only on the properties of the microphones, 
directionality, orientation and position in space.

Analytic expressions for the DFC exist only for omnidirectional sensors
at arbitrary positions, the well-known sinc function of the
wavenumber-distance product, and for first-order directional microphones
with arbitrary orientations, see e.g. [ref18]. 
For more general directivities, [ref.8] shows that the DCM can be 
pre-computed through the expansion of the microphone/beamformer patterns 
into SHD coefficients.


---
## References <a name="references"></a>

References mentioned in the code and the examples:

  1.  Moreau, S., Daniel, J., Bertet, S., 2006, 
      3D sound field recording with higher order ambisonics-objective measurements and validation of spherical microphone. 
      In Audio Engineering Society Convention 120.

  2.  Mh Acoustics Eigenmike, https://mhacoustics.com/products#eigenmike1

  3.  Bernsch?tz, B., P?rschmann, C., Spors, S., Weinzierl, S., Verst?rkung, B., 2011. 
      Soft-limiting der modalen amplitudenverst?rkung bei sph?rischen mikrofonarrays im plane wave decomposition verfahren. 
      Proceedings of the 37. Deutsche Jahrestagung f?r Akustik (DAGA 2011)

  4.  Jin, C.T., Epain, N. and Parthy, A., 2014. 
      Design, optimization and evaluation of a dual-radius spherical microphone array. 
      IEEE/ACM Transactions on Audio, Speech, and Language Processing, 22(1), pp.193-204.

  5.  Elko, G.W., 2004. Differential microphone arrays. 
      In Audio signal processing for next-generation multimedia communication systems (pp. 11-65). Springer.

  6.  Zotter, F., Pomberger, H. and Noisternig, M., 2012. 
      Energy-preserving ambisonic decoding. Acta Acustica united with Acustica, 98(1), pp.37-47.

  7.  Koretz, A. and Rafaely, B., 2009. 
      Dolph?Chebyshev beampattern design for spherical arrays. 
      IEEE Transactions on Signal Processing, 57(6), pp.2417-2420.

  8.  Politis, A., 2016. 
      Diffuse-field coherence of sensors with arbitrary directional responses. arXiv preprint arXiv:1608.07713.

  9.  Hafizovic, I., Nilsen, C.I.C. and Holm, S., 2012. 
      Transformation between uniform linear and spherical microphone arrays with symmetric responses. 
      IEEE Transactions on Audio, Speech, and Language Processing, 20(4), pp.1189-1195.

  10. Politis, A. and Pulkki, V., 2016. 
      Acoustic intensity, energy-density and diffuseness estimation in a directionally-constrained region. 
      arXiv preprint arXiv:1609.03409.

  11. Thiergart, O. and Habets, E.A., 2014. 
      Extracting reverberant sound using a linearly constrained minimum variance spatial filter. 
      IEEE Signal Processing Letters, 21(5), pp.630-634.

  12. Thiergart, O., Taseska, M. and Habets, E.A., 2014. 
      An informed parametric spatial filter based on instantaneous direction-of-arrival estimates. 
      IEEE/ACM Transactions on Audio, Speech, and Language Processing, 22(12), pp.2182-2196.

  13. Merimaa, J. and Pulkki, V., 2005. 
      Spatial impulse response rendering I: Analysis and synthesis. 
      Journal of the Audio Engineering Society, 53(12), pp.1115-1127.

  14. Ahonen, J. and Pulkki, V., 2009. 
      Diffuseness estimation using temporal variation of intensity vectors. 
      In 2009 IEEE Workshop on Applications of Signal Processing to Audio and Acoustics (WASPAA).

  15. Politis, A., Delikaris-Manias, S. and Pulkki, V., 2015. 
      Direction-of-arrival and diffuseness estimation above spatial aliasing for symmetrical directional microphone arrays. 
      In 2015 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP).

  16. Gover, B.N., Ryan, J.G. and Stinson, M.R., 2002. 
      Microphone array measurement system for analysis of directional and spatial variations of sound fields. 
      The Journal of the Acoustical Society of America, 112(5), pp.1980-1991.

  17. Epain, N. and Jin, C.T., 2016. 
      Spherical Harmonic Signal Covariance and Sound Field Diffuseness. 
      IEEE/ACM Transactions on Audio, Speech, and Language Processing, 24(10), pp.1796-1807.

  18. Elko, G.W., 2001. 
      Spatial coherence functions for differential microphones in isotropic noise fields. 
      In Microphone Arrays (pp. 61-85). Springer Berlin Heidelberg.
