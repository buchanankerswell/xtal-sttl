# xtal-sttl

## Introduction
`xtal-sttl` is an open source web-based app for calculating the settling velocity of spherical crystals in a silicate melt. This involves three calculations:

1. Calculate melt density using a ten-component system: SiO2-TiO2-Al2O3-Fe2O3-FeO-MgO-CaO-Na2O-K2O-H2O

2. Calculate melt viscosity using the VFT model of Hess and Dingwell (1996):


ln(\eta)=\left[a_1+a_2~ln(water)\right]+\frac{\left[b_1+b_2~ln(water)\right]}{T-\left[c_1+c_2~ln(water)\right]}


3. Approximate settling velocity using Stoke's equation:


$$ v = frac{ 2 g radius_{xtal}^2 \left( \rho_{xtal} - \rho_{liq} \right) }{ 9 \nu_{liq} } $$
