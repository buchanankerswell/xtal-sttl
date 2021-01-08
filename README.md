# xtal-sttl

## Introduction
`xtal-sttl` is an open source web-based app for calculating the settling velocity of spherical crystals in a silicate melt. This involves three calculations:

1. Calculate melt density using a ten-component system: SiO2-TiO2-Al2O3-Fe2O3-FeO-MgO-CaO-Na2O-K2O-H2O

2. Calculate melt viscosity using the VFT model of Hess and Dingwell (1996):


<img src="http://www.sciweavers.org/tex2img.php?eq=ln%28%5Ceta%29%3D%5Cleft%5Ba_1%2Ba_2~ln%28water%29%5Cright%5D%2B%5Cfrac%7B%5Cleft%5Bb_1%2Bb_2~ln%28water%29%5Cright%5D%7D%7BT-%5Cleft%5Bc_1%2Bc_2~ln%28water%29%5Cright%5D%7D&bc=White&fc=Black&im=png&fs=18&ff=fourier&edit=0" align="center" border="0" alt="ln(\eta)=\left[a_1+a_2~ln(water)\right]+\frac{\left[b_1+b_2~ln(water)\right]}{T-\left[c_1+c_2~ln(water)\right]}" width="548" height="58" />


3. Approximate settling velocity using Stoke's equation:


$$ v = frac{ 2 g radius_{xtal}^2 \left( \rho_{xtal} - \rho_{liq} \right) }{ 9 \nu_{liq} } $$
