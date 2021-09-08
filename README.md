## Description

TBilateral is a spatial smoothing filter that uses the bilateral filtering algorithm. It does a nice job of smoothing while retaining picture structure.

This is [a port of the VapourSynth plugin TBilateral](https://github.com/dubhater/vapoursynth-tbilateral).

### Requirements:

- AviSynth 2.60 / AviSynth+ 3.4 or later

- Microsoft VisualC++ Redistributable Package 2022 (can be downloaded from [here](https://github.com/abbodi1406/vcredist/releases))

### Usage:

```
vsTBilateral (clip, clip "ppclip", int "diameterY", int "diameterU", int "diameterV", float "sdevY", float "sdevU", float "sdevV", float "idevY", float "idevU", float "idevV", float "csY", float "csU", float "csV", bool "d2", int "kerns", int "kerni", int "restype", int "y", int "u", int "v")
```

### Parameters:

- clip\
    A clip to process. It must be in 8..16-bit planar format.
    
- ppclip\
    Specifies a pre-filtered clip for TBilateral to take pixel values from when doing the luminance difference calculations.\
    The general recommendation for pre-processing is a gaussian blur with standard deviation equal to the sdev settings being used. Using a prefiltered clip should help in removing impulse noise (i.e. outliers with very large pixel differences) which standard bilateral filtering will not touch. It does tend to produce artifacts sometimes, especially around very fine details. Another recommendation for pre-processing is a center-weighted median or adaptive median.\
    Must have the same format, dimensions, and number of frames as clip.\
    Default: None.
    
- diameterY, diameterU, diameterV\
    This sets the size of the diameter of the filtering window. Larger values mean more pixels will be included in the average, but are also slower.\
    Must be an odd number greater than 1. This must be less than the width of the video and less than the height of the video.\
    Default: diameterY = 5; diameterU = 5; diameterV = 5.
    
- sdevY, sdevU, sdevV\
    This sets the spatial deviations. The larger sdev is, the less effect distance will have in the weighting of pixels in the average. That is, as you increase sdev, distant pixels will have more weight.\
    To get a better idea of what this setting does try setting idev to high values, and then gradually increase sdev from 0 on up while keeping idev constant.\
    Increasing this setting will increase the strength of the smoothing.\
    Must be greater than 0.\
    Default: sdevY = 1.4; sdevU = 1.4; sdevV = 1.4.
    
- idevY, idevU, idevV\
    This sets the pixel intensity deviations (or color deviations in the case of the chroma planes). The larger idev is, the less effect pixel difference will have in the weighting of pixels in the average. That is, as you increase idev, pixels that are very different from the current pixel will have more weight.\
    Try increasing this setting while keeping sdev constant to get a better idea of what this does. Increasing this setting will increase the strength of the smoothing.\
    Must be greater than 0.\
    Default: idevY = 7.0; idevU = 7.0; idevV = 7.0.
    
- csY, csU, csV\
    This value is multiplied to the center pixel's spatial weight value. A value of 1 does nothing, less than 1 means the center pixel will have less weight than normal, greater than 1 means the center pixel will have more weight than normal, 0 gives the center pixel no weight.\
    Must be at least 0. Setting cs to 0 will give you SUSAN denoising.\
    Default: csY = 1.0; csU = 1.0; csV = 1.0.
    
- d2\
    This setting makes TBilateral use the second derivative instead of the first when doing the intensity calculations. Using d2 should give better results on smooth gradients or anything that fits the piecewise linear model. Setting d2 to False will give better results on images that have uniformly colored areas with sharp edges (anything that fits the piecewise constant model). The actual difference between the two is usually not big for most sources. The effect is rather subtle.\
    Default: False.
    
- kerns\
    This specifies what kernel is used for the domain weights.\
    0: Andrews' wave\
    1: El Fallah Ford\
    2: Gaussian\
    3: Huber’s mini-max\
    4: Lorentzian\
    5: Tukey bi-weight\
    6: Linear descent\
    7: Cosine\
    8: Flat\
    9: Inverse\
    See the following paper for a description of all the kernels and their properties: http://www.prasa.org/proceedings/2003/prasa03-09.pdf \
    Default: 2.
    
- kerni\
    This specifies what kernel is used for the range weights. The possible choices are the same as for kerns.\
    Default: 2.
    
- restype\
    This specifies how the weights and pixel values are combined to obtain the final result.\
    0: Mean (weighted average)\
    1: Median (weighted median)\
    2: CW-Median (weighted median + extra center pixel weight)\
    3: Multiple linear regression (best fit plane)\
    Default: 0.
    
- y, u, v\
    Planes to process.\
    1: Return garbage.\
    2: Copy plane.\
    3: Process plane. Always process planes when the clip is RGB.\
    Default: y = 3; u = 3; v = 3.
