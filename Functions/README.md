The image "breast_ct.jpg" was obtained from a contrast-enhanced cone-beam breast CT image in work from the Diagnostic Breast Center Göttingen (https://doi.org/10.1016/j.tranon.2017.08.010) under the Creative Commons Attribution-NonCommercial-No Derivatives License (CC BY NC ND) (https://creativecommons.org/licenses/by-nc-nd/4.0/).
[Reference: Uhlig, J., Fischer, U., von Fintel, E., Stahnke, V., Perske, C., Lotz, J., & Wienbeck, S. (2017). Contrast Enhancement on Cone-Beam Breast-CT for Discrimination of Breast Cancer Immunohistochemical Subtypes. Translational oncology, 10(6), 904-910.]

Prior to running any MATLAB code, you need to build the MEX function “msfm2d” by running the command “mex -compatibleArrayDims msfm2d.c” —- I have included the compiled MEX function for Mac OS X and Linux, but for other operating systems you have to run “mex -compatibleArrayDims msfm2d.c”
The MEX code for the fast marching method came from the Accurate Fast Marching Code by Dirk-Jan Kroon: https://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching. See license information below:

Copyright (c) 2009, Dirk-Jan Kroon
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



