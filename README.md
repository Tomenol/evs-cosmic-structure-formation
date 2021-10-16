# EVS of cosmic structures
Based on a Fortran code developed by Pier-Stefano Corasaniti (LUTH), the goal of this Python/C code is to provide an efficient yet easy to use library to compute the mass of the most massive galaxy cluster in a given redshift bin. 

By using Extreme Value Statistics and assuming the validity of the LCDM cosmology, one can build a statistical model to predict the mass of the most massive galaxy cluster at a given redshift z. This code can be used to compute those predictions under the thin redshift bin hypothesis (dz <= 0.1). 

# Sources :
This library is based on the Fortran code developed by Pier-Stefano Corasaniti (LUTH)

The following papers were also used for the development of this code :
- Harrison & Coles, Testing Cosmology with Extreme Galaxy Clusters : https://arxiv.org/pdf/1111.1184.pdf
- Tinker et al., Toward a halo mass function for precision cosmology: the limits of universality : https://arxiv.org/pdf/0803.2706.pdf
- Chongchitnan, On the Abundance of Extreme Voids : https://arxiv.org/pdf/1502.07705.pdf
- Planck Collaboration, Planck 2018 results : https://arxiv.org/pdf/1807.06209.pdf

# Author :
Thomas MAYNADIE
