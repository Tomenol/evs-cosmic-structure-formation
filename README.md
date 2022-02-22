# EVS of cosmic structures
Based on a Fortran code developed by Pier-Stefano Corasaniti, the goal of this Python/C code is to provide an efficient yet easy to use library to compute the mass of the most massive galaxy cluster in a given redshift bin. 

By using Extreme Value Statistics and assuming the validity of the LCDM cosmology, one can build a statistical model to predict the mass of the most massive galaxy cluster at a given redshift z. This code can be used to compute those predictions under the thin redshift bin hypothesis (dz <= 0.1). 

# Results
Figure 1.0 : Influence of the comsological parameters over the EVS predictions

![image](https://user-images.githubusercontent.com/54234406/155146710-767388c7-648a-4891-b182-bfa92fca3831.png)

Figure 2.0 : Influence of the mass bias over the measured masses of the most massive clusters (PSZ2 galaaxy survey) and comparison with the EVS predictions (Planck 2015)

![image](https://user-images.githubusercontent.com/54234406/155146511-e45218f3-01f1-41c1-a9c7-a31cedab74f7.png)

# References :
This library is based on the Fortran code developed by Pier-Stefano Corasaniti

The following papers were used for the development of this code :
- Harrison & Coles, Testing Cosmology with Extreme Galaxy Clusters : https://arxiv.org/pdf/1111.1184.pdf
- Tinker et al., Toward a halo mass function for precision cosmology: the limits of universality : https://arxiv.org/pdf/0803.2706.pdf
- Chongchitnan, On the Abundance of Extreme Voids : https://arxiv.org/pdf/1502.07705.pdf
- Planck Collaboration, Planck 2018 results : https://arxiv.org/pdf/1807.06209.pdf

# Author :
Thomas MAYNADIE
