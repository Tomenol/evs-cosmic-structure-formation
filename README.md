# EVS of cosmic structures
Based on a Fortran code developed by Pier-Stefano Corasaniti (LUTH), the goal of this Python/C code is to provide an efficient yet easy to use library to compute the mass of the most massive galaxy cluster in a given redshift bin. 

By using Extreme Value Statistics and assuming the validity of the LCDM cosmology, one can build a statistical model to predict the mass of the most massive galaxy cluster at a given redshift z. This code can be used to compute those predictions under the thin redshift bin hypothesis (dz <= 0.1). 

# Results
By using Planck 2018 cosmological parameters and setting survey characteristics to the default values, we can obtain the following plot (see examples/plotMmax.py) :

![image](https://user-images.githubusercontent.com/54234406/137623097-647ceb83-5c92-4351-b26a-697ac3be04bc.png)

where the red curve shows the most likely mass of the most massive galaxy cluster as a function of z and the blue curves show the two sigma confidence interval.

# Sources :
This library is based on the Fortran code developed by Pier-Stefano Corasaniti (LUTH)

The following papers were used for the development of this code :
- Harrison & Coles, Testing Cosmology with Extreme Galaxy Clusters : https://arxiv.org/pdf/1111.1184.pdf
- Tinker et al., Toward a halo mass function for precision cosmology: the limits of universality : https://arxiv.org/pdf/0803.2706.pdf
- Chongchitnan, On the Abundance of Extreme Voids : https://arxiv.org/pdf/1502.07705.pdf
- Planck Collaboration, Planck 2018 results : https://arxiv.org/pdf/1807.06209.pdf

# Author :
Thomas MAYNADIE
