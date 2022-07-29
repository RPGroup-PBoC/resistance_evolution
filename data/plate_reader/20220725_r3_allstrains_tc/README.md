# 2022-07-22 Plate Reader Growth Measurement

## Purpose


## Strain Information

| Plasmid | Genotype | Host Strain | Shorthand |
| :------ | :------- | ----------: | --------: |
| -| - | MG1655 | MG1655 |
| -| galK<>3.19_RiboJ_tetA_GFP| MG1655 | 3.19 |
| -| galK<>IW_RiboJ_tetA_GFP| MG1655 | IW |
| -| galK<>WTlac_RiboJ_tetA_GFP| MG1655 | WT |
| -| galK<>lacUV5_RiboJ_tetA_GFP| MG1655 | UV5 |




## Plate Layout

**96 plate layout**

![plate layout](output/plate_layout.png)

## Notes & Observations
For this experiment we used new 96-well plates from Eppendorf that are build with an extra space outside the outer wells that can be filled with water to prevent evaporation in the wells on the outside. In this run we used a blank in each corner.

## Analysis Files

**Whole Plate Growth Curves**
![plate layout](output/20220722_r4_all_curves.pdf)

**Whole Plate Growth Rate Inferences**
![plate layout](output/20220722_r4_all_curves_with_th.pdf)

## Experimental Protocol

1. Pick a colony from a fresh plate (< 2 weeks) and grow in 3ml of LB in the incubator at 37C without shaking for 15h overnight.
2. Dilute 1:100 into prewarmed M9 + 0.4% glucose + 0.1% thiamine + trace elements and grow at 37C in the shaker.
3. Prepare well plate with prewarmed media (300µl per well) and prewarm plate reader.
4. Grow cultures to OD ~ 0.4.
5. Add 15µl of cultures to wells to achieve 1:20 dilution.


## Plate reader settings
| Temperature | 37C |
| Shaking Mode | linear |
| Shaking speed | 1096 | 
| Time between reads | 7 min |
| Total time | 24 h |

## Conclusion