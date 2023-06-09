ice2021 repository

This is a repository for ice surface data as determined by Scanning electron microscopy & Gauss-Newton Bayesian Framework (GNBF)

Reference:
Butterfield et al, Quantitative three-dimensional ice roughness from scanning electron microscopy, J. Geophys. Res.: Atmospheres, 122, 3023-3041 (2017).


How to install dependencies:

pip install numpy-stl # see https://pypi.python.org/pypi/numpy-stl:
pip install f90nml
pip install reliability

This is work done by:  
Steven Neshyba (prior to 2020)  
Steven Neshyba and Kayden Lynch (2020)  
Steven Neshyba and Natalie Reszka (2021)  
Steven Neshyba, Tia Bottger, and Rohan Crossland (2023)

Data flow

GNBF_makextlvecs
Inputs: Images of a smooth ice crystal 
Output: Xtlvecs.nml (nameless with vectors)
Notes: The vectors define the orientation of a hexagonal prism:

Looking at the image as displayed by GNBF_makextlvecs, the coordinate system is

+x is to the right
+y is down
+z is into the plane
xyz form a right-handed coordinate system.

The three colored vectors are:
bvec (green) is x'
avec (red) is y'
cvec (blue) is z'
x'y'z' also form a right-handed coordinate system.

We want to vary xpos, ypos, alpha, beta, and gamma, such that
-The origin (xpos & ypos) put the center on a corner of the crystal
-cvec (blue) points along the c-axis of the crystal (perpendicular to basal plane)
-bvec (green) points along the boundary between the basal and one pyramidal
-avec (red) points along the boundary between the basal and another pyramidal

Strategy:
1. Move xpos & ypos bars into position
2. Adjust gamma until blue lines up with a c-axis
3. Adjust alpha and beta until red and green line up

GNBF_calibrate_from_Xtlvecs
Input: Xtlvecs.nml
Output: calibration.nml (slopes and intercepts)
Notes: The slopes and intercepts define the response of each of four detectors to the s-variable.

GNBF_retrieve 
Inputs: calibration.nml, segments.nml, SEM immages
Outputs: Segments_compr.npz, Segments_retrieved.jpg
Notes: The output .npz file contains the surface morphology of the segments.

GNBF_grid2stlwskirt
Input: Segments_compr.npz
Output: Segments_retrievedwskirt.stl
Notes: The stereolithography file (.stl) is a useful way to inspect surface morphology. There are various ways of visualizing the result. A simple one is, once the data are pushed to GitHub, to let Firefox display by clicking on the file.

GNBF_flattenandfilter
Input: Segments_compr.npz
Outputs: Segments_compr_flat.npz, Segments_compr_flat_filt.npz, Segments_compr_flat_filt_0_vx5.stl
Notes: Filtering is done by Fourier band-pass. Each segment is also flattened first. Stereolithography file (.stl) files are for each segment.

GNBF_analyzeroughness_of_flattened_npz
Input: Segments_compr_flat_filt.npz
Output: Segments_compr_flat_filt_roughness.jpg
Notes: The .jpg includes sigma and eta (Weibull) of either all the segments combined, or a subset of them defined in the code as “accumlist.” There is also a lot of numerical data generated during the run.
