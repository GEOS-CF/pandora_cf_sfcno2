Python script to read Pandora L2_rnvh3p1-8 surface NO2 observations and calculate corresponding quantities from GEOS-CF. These quantities are:
- cf_no2_sfcmr: GEOS-CF surface mixing ratio of NO2, in ppbv
- cf_no2_sfcconc: GEOS-CF surface concentration of NO2, in mol/m3
- cf_no2_l1col: GEOS-CF partial column amount in layer 1 [moles per square meter]
- cf_no2_pbl: GEOS-CF boundary layer height, [meters]

The corresponding Pandora quantities are also provided. Conversion between surface mixing ratios and concentrations is done using environmental quantities obtained from GEOS-CF.


EXAMPLE:

`wget https://data.pandonia-global-network.org/WashingtonDC/Pandora140s1/L2/Pandora140s1_WashingtonDC_L2_rnvh3p1-8.txt`

`python pandora_cf_sfcno2.py -i Pandora140s1_WashingtonDC_L2_rnvh3p1-8.txt`
