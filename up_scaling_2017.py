# -*- coding: utf-8 -*-
"""
Created on Thu May 11 09:32:35 2023

@author: mmo990
"""
import ast
import math
import matplotlib 
from matplotlib.ticker import ScalarFormatter

def surface_areas(fname, wbtype):     
    with open(fname) as f:
        s = f.read()
        areas_list = ast.literal_eval(s)
    areas_km = areas_list[0]
    areas_perc = areas_list[1]
           
    area = areas_km[wbtype]
           
    return area

# Read in surface area coverage data ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Areas convert from km2 to m2
areas_file = r'.../2015_wb_areas.txt'
lake_area_17 = surface_areas(areas_file, 'Lake') * 1000000
pond_area_17 = surface_areas(areas_file, 'Pond') * 1000000
fluvial_area_17 = surface_areas(areas_file, 'Stream') * 1000000
tundra_area_17 = surface_areas(areas_file, 'Tundra') * 1000000


areas_file_15_cropped = r'.../2015_cropped_wb_areas.txt'
pond_area_15_cropped = surface_areas(areas_file_15_cropped, 'Pond') * 1000000
fluvial_area_15_cropped = surface_areas(areas_file_15_cropped, 'Stream') * 1000000
areas_file_17 = r'.../2017_wb_areas.txt'
flood_area_17 = (surface_areas(areas_file_17, 'Pond') * 1000000 - pond_area_15_cropped) + (surface_areas(areas_file_17, 'Stream') * 1000000 - fluvial_area_15_cropped)
tundra_area_17 =  tundra_area_16 - flood_area_17
flood_perc = flood_area_17 / (lake_area_17 + pond_area_17 + fluvial_area_17 + tundra_area_17 + flood_area_17 )
tundra_perc = tundra_area_17 / (lake_area_17 + pond_area_17 + fluvial_area_17 + tundra_area_17 + flood_area_17 )

print('An estimated {:2.2%} of our study area was flooded in 2017'.format(flood_perc))
print('Tundra covers an estimated {:2.2%} of our study area in 2017'.format(tundra_perc))

# Calculate CO2 fluxes 2017 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Lake CO2
lake_co2 = Fluxes('diss CO2.ppm', 'CO2 air', 'Lake', '2017')
print(lake_co2.gas)
print(lake_co2.wb_type)
lake_co2.subset_df(df)
lake_co2.CO2_conc()
print(lake_co2.df)
lake_co2.outlier_removal(0.25, 0.75, '.../lake_co2_17.csv')
lake_co2.CO2_eq()
lake_co2.Sc()
lake_co2.k600_cc()
lake_co2.k()

lake_co2_s_boot = lake_co2.bootstrapping('diss co2', 1000)
lake_co2_eq_boot = lake_co2.bootstrapping('CO2 eq', 1000)
lake_co2_k_boot = lake_co2.bootstrapping('k', 1000)

lake_co2.fluxes_calc(lake_co2_s_boot, lake_co2_eq_boot, lake_co2_k_boot, ['diss co2', 'co2 eq', 'k'], 'flux co2', 'flux co2-c mmol', 'flux co2-2 mg', '', '')


# Pond CO2
pond_co2 = Fluxes('diss CO2.ppm', 'CO2 air', 'Pond', '2017')
print(pond_co2.gas)
print(pond_co2.wb_type)
pond_co2.subset_df(df)
pond_co2.CO2_conc()
print(pond_co2.df)
pond_co2.outlier_removal(0.25, 0.75, '.../pond_co2_17.csv')
pond_co2.CO2_eq()
pond_co2.Sc()
pond_co2.uniform_k600_lit('C:/Users/mmo990/surfdrive/Paper1/Data/k600_literature.xlsx', 1000)
pond_co2.k()

pond_co2_s_boot = pond_co2.bootstrapping('diss co2', 1000)
pond_co2_eq_boot = pond_co2.bootstrapping('CO2 eq', 1000)
pond_co2_k_boot = pond_co2.bootstrapping('k', 1000)

pond_co2.fluxes_calc(pond_co2_s_boot, pond_co2_eq_boot, pond_co2_k_boot, ['diss co2', 'co2 eq', 'k'], 'flux co2', 'flux co2-c mmol', 'flux co2-2 mg', '', '')


# Fluvial CO2 
fluvial_co2 = Fluxes('diss CO2.ppm', 'CO2 air', 'Fluvial', '2017')
print(fluvial_co2.gas)
print(fluvial_co2.wb_type)
fluvial_co2.subset_df(df)
fluvial_co2.CO2_conc()
print(fluvial_co2.df)
fluvial_co2.outlier_removal(0.25, 0.75, '.../fluvial_co2_17.csv')
fluvial_co2.CO2_eq()
fluvial_co2.Sc()
fluvial_co2.uniform_k600_lit('.../k600_literature.xlsx', 1000)
fluvial_co2.k()

fluvial_co2_s_boot = fluvial_co2.bootstrapping('diss co2', 1000)
fluvial_co2_eq_boot = fluvial_co2.bootstrapping('CO2 eq', 1000)
fluvial_co2_k_boot = fluvial_co2.bootstrapping('k', 1000)

fluvial_co2.fluxes_calc(fluvial_co2_s_boot, fluvial_co2_eq_boot, fluvial_co2_k_boot, ['diss co2', 'co2 eq', 'k'], 'flux co2', 'flux co2-c mmol', 'flux co2-2 mg', '', '')


# Flood CO2
flood_co2 = Fluxes('diss CO2.ppm', 'CO2 air', 'Flood', '2017')
print(flood_co2.gas)
print(flood_co2.wb_type)
flood_co2.subset_df(df)
flood_co2.CO2_conc()
print(flood_co2.df)
flood_co2.outlier_removal(0.25, 0.75, '.../flood_co2_17.csv')
flood_co2.CO2_eq()
flood_co2.Sc()
flood_co2.k600_cc()
flood_co2.k()

flood_co2_s_boot = flood_co2.bootstrapping('diss co2', 1000)
flood_co2_eq_boot = flood_co2.bootstrapping('CO2 eq', 1000)
flood_co2_k_boot = flood_co2.bootstrapping('k', 1000)

flood_co2.fluxes_calc(flood_co2_s_boot, flood_co2_eq_boot, flood_co2_k_boot, ['diss co2', 'co2 eq', 'k'], 'flux co2', 'flux co2-c mmol', 'flux co2-2 mg', '', '')

# Calculate CH4 fluxes 2017 --------------------------------------------------------------------------------------
# Lake CH4 
lake_ch4 = Fluxes('diss CH4', 'CH4 air', 'Lake', '2017')
lake_ch4.subset_df(df)
lake_ch4.outlier_removal(0.25, 0.75, '.../lake_ch4_17.csv')
lake_ch4.CH4_eq()
lake_ch4.Sc()
lake_ch4.k600_cc()
lake_ch4.k()

lake_ch4_s_boot = lake_ch4.bootstrapping('diss CH4', 1000)
lake_ch4_eq_boot = lake_ch4.bootstrapping('CH4 eq', 1000)
lake_ch4_k_boot = lake_ch4.bootstrapping('k', 1000)

lake_ch4.fluxes_calc(lake_ch4_s_boot, lake_ch4_eq_boot, lake_ch4_k_boot, ['diss CH4', 'CH4 eq', 'k'], 'flux ch4', 'flux ch4-c mmol', 'flux ch4-2 mg', 'flux ch4 co2eq', '')


# Pond CH4 
pond_ch4 = Fluxes('diss CH4', 'CH4 air', 'Pond', '2017')
pond_ch4.subset_df(df)
pond_ch4.outlier_removal(0.25, 0.75, '.../pond_ch4_17.csv')
pond_ch4.CH4_eq()
pond_ch4.Sc()
pond_ch4.uniform_k600_lit('.../k600_literature.xlsx', 1000)
pond_ch4.k()

pond_ch4_s_boot = pond_ch4.bootstrapping('diss CH4', 1000)
pond_ch4_eq_boot = pond_ch4.bootstrapping('CH4 eq', 1000)
pond_ch4_k_boot = pond_ch4.bootstrapping('k', 1000)

pond_ch4.fluxes_calc(pond_ch4_s_boot, pond_ch4_eq_boot, pond_ch4_k_boot, ['diss CH4', 'CH4 eq', 'k'], 'flux ch4', 'flux ch4-c mmol', 'flux ch4-2 mg', 'flux ch4 co2eq', '')

# Fluvial CH4
fluvial_ch4 = Fluxes('diss CH4', 'CH4 air', 'Fluvial', '2017')
fluvial_ch4.subset_df(df)
fluvial_ch4.outlier_removal(0.25, 0.75, '.../fluvial_ch4_17.csv')
fluvial_ch4.CH4_eq()
fluvial_ch4.Sc()
fluvial_ch4.uniform_k600_lit('...k600_literature.xlsx', 1000)
fluvial_ch4.k()

fluvial_ch4_s_boot = fluvial_ch4.bootstrapping('diss CH4', 1000)
fluvial_ch4_eq_boot = fluvial_ch4.bootstrapping('CH4 eq', 1000)
fluvial_ch4_k_boot = fluvial_ch4.bootstrapping('k', 1000)

fluvial_ch4.fluxes_calc(fluvial_ch4_s_boot, fluvial_ch4_eq_boot, fluvial_ch4_k_boot, ['diss CH4', 'CH4 eq', 'k'], 'flux ch4', 'flux ch4-c mmol', 'flux ch4-2 mg', 'flux ch4 co2eq', '')


# Flood Ch4
flood_ch4 = Fluxes('diss CH4', 'CH4 air', 'Flood', '2017')
flood_ch4.subset_df(df)
flood_ch4.outlier_removal(0.25, 0.75, '.../flood_ch4_17.csv')
flood_ch4.CH4_eq()
flood_ch4.Sc()
flood_ch4.k600_cc()
flood_ch4.k()

flood_ch4_s_boot = flood_ch4.bootstrapping('diss CH4', 1000)
flood_ch4_eq_boot = flood_ch4.bootstrapping('CH4 eq', 1000)
flood_ch4_k_boot = flood_ch4.bootstrapping('k', 1000)

flood_ch4.fluxes_calc(flood_ch4_s_boot, flood_ch4_eq_boot, flood_ch4_k_boot, ['diss CH4', 'CH4 eq', 'k'], 'flux ch4', 'flux ch4-c mmol', 'flux ch4-2 mg', 'flux ch4 co2eq', '')



# Calculate N2O fluxes 2017  --------------------------------------------------------------------------------------
# Lake N2O
lake_n2o = Fluxes('diss N2O', 'N2O air', 'Lake', '2017')
lake_n2o.subset_df(df)
lake_n2o.outlier_removal(0.25, 0.75, '.../lake_n2o_17.csv')
lake_n2o.N2O_eq()
lake_n2o.Sc()
lake_n2o.k600_cc()
lake_n2o.k()

lake_n2o_s_boot = lake_n2o.bootstrapping('diss N2O', 1000)
lake_n2o_eq_boot = lake_n2o.bootstrapping('N2O eq', 1000)
lake_n2o_k_boot = lake_n2o.bootstrapping('k', 1000)

lake_n2o.fluxes_calc(lake_n2o_s_boot, lake_n2o_eq_boot, lake_n2o_k_boot, ['diss N2O', 'N2O eq', 'k'], 'flux n2o', '', '', 'flux n2o co2eq', 'cooling_flux n2o co2eq')

# Lake N2O
pond_n2o = Fluxes('diss N2O', 'N2O air', 'Pond', '2017')
pond_n2o.subset_df(df)
pond_n2o.outlier_removal(0.25, 0.75, '...g/pond_n2o_17.csv')
pond_n2o.N2O_eq()
pond_n2o.Sc()
pond_n2o.uniform_k600_lit('.../k600_literature.xlsx', 1000)
pond_n2o.k()

pond_n2o_s_boot = pond_n2o.bootstrapping('diss N2O', 1000)
pond_n2o_eq_boot = pond_n2o.bootstrapping('N2O eq', 1000)
pond_n2o_k_boot = pond_n2o.bootstrapping('k', 1000)

pond_n2o.fluxes_calc(pond_n2o_s_boot, pond_n2o_eq_boot, pond_n2o_k_boot, ['diss N2O', 'N2O eq', 'k'], 'flux n2o', '', '', 'flux n2o co2eq', 'cooling_flux n2o co2eq')

# Fluvial N2O
fluvial_n2o = Fluxes('diss N2O', 'N2O air', 'Fluvial', '2017')
fluvial_n2o.subset_df(df)
fluvial_n2o.outlier_removal(0.25, 0.75, '.../fluvial_n2o_17.csv')
fluvial_n2o.N2O_eq()
fluvial_n2o.Sc()
fluvial_n2o.uniform_k600_lit('.../k600_literature.xlsx', 1000)
fluvial_n2o.k()

fluvial_n2o_s_boot = fluvial_n2o.bootstrapping('diss N2O', 1000)
fluvial_n2o_eq_boot = fluvial_n2o.bootstrapping('N2O eq', 1000)
fluvial_n2o_k_boot = fluvial_n2o.bootstrapping('k', 1000)

fluvial_n2o.fluxes_calc(fluvial_n2o_s_boot, fluvial_n2o_eq_boot, fluvial_n2o_k_boot, ['diss N2O', 'N2O eq', 'k'], 'flux n2o', '', '', 'flux n2o co2eq', 'cooling_flux n2o co2eq')


# Flood N2O
flood_n2o = Fluxes('diss N2O', 'N2O air', 'Flood', '2017')
flood_n2o.subset_df(df)
flood_n2o.outlier_removal(0.25, 0.75, '.../flood_n2o_17.csv')
flood_n2o.N2O_eq()
flood_n2o.Sc()
flood_n2o.k600_cc()
flood_n2o.k()

flood_n2o_s_boot = flood_n2o.bootstrapping('diss N2O', 1000)
flood_n2o_eq_boot = flood_n2o.bootstrapping('N2O eq', 1000)
flood_n2o_k_boot = flood_n2o.bootstrapping('k', 1000)

flood_n2o.fluxes_calc(flood_n2o_s_boot, flood_n2o_eq_boot, flood_n2o_k_boot, ['diss N2O', 'N2O eq', 'k'], 'flux n2o', '', '', 'flux n2o co2eq', 'cooling_flux n2o co2eq')


# C emissions inland waters CO2-C + CH4-C 2017 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Pond CO2-C + CH4-C emissions
pond_co2_17 = (pond_co2.df_boot['flux co2-2 mg'] * pond_area_17) *10**-9
std_pond_co2_17  = np.std(pond_co2_17)
unc_pond_co2_17  = std_pond_co2_17 + math.sqrt((0.15*np.median(pond_co2_17))**2 + (0.15*(pond_area_17/1000000))**2)
pond_ch4_17 = (pond_ch4.df_boot['flux ch4-2 mg'] * pond_area_17) *10**-9
std_pond_ch4_17  = np.std(pond_ch4_17)
unc_pond_ch4_17  = std_pond_ch4_17 + math.sqrt((0.15*np.median(pond_ch4_17))**2 + (0.15*(pond_area_17/1000000))**2)
pond_c_17 = pond_co2_17 + pond_ch4_17
median_pond_c_17 = np.median(pond_c_17)
conf_int_pond_c_17 = np.percentile(pond_c_17, [5, 95])
std_pond_c_17 = np.std(pond_c_17)
unc_pond_c_17 = std_pond_c_17 + math.sqrt((0.15*median_pond_c_17)**2 + (0.15*(pond_area_17/1000000))**2)
print('Median pond CO2-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(np.median(pond_co2_17), unc_pond_co2_17))
print('Median pond CH4-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(np.median(pond_ch4_17), unc_pond_ch4_17))
print('Median pond C emissions (CO2-C + CH4-C) are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(median_pond_c_17, unc_pond_c_17))    
print('Pond C emissions (CO2-C + CH4-C) range from {:0.2f} to {:0.2f} (Mg C d-1) [95% CI]'.format(conf_int_pond_c_17[0], conf_int_pond_c_17[1]))  

# Lake CO2-C + CH4-C emissions
lake_co2_17 = (lake_co2.df_boot['flux co2-2 mg'] * lake_area_17) *10**-9
std_lake_co2_17  = np.std(lake_co2_17)
unc_lake_co2_17  = std_lake_co2_17 + math.sqrt((0.15*np.median(lake_co2_17))**2 + (0.15*(lake_area_17/1000000))**2)
lake_ch4_17 = (lake_ch4.df_boot['flux ch4-2 mg'] * lake_area_17) *10**-9
std_lake_ch4_17  = np.std(lake_ch4_17)
unc_lake_ch4_17  = std_lake_ch4_17 + math.sqrt((0.15*np.median(lake_ch4_17))**2 + (0.15*(lake_area_17/1000000))**2)
lake_c_17 = lake_co2_17 + lake_ch4_17
median_lake_c_17 = np.median(lake_c_17)
conf_int_lake_c_17 = np.percentile(lake_c_17, [5, 95])
std_lake_c_17 = np.std(lake_c_17)
unc_lake_c_17 = std_lake_c_17 + math.sqrt((0.15*np.median(lake_c_17))**2 + (0.15*(lake_area_17/1000000))**2)
print('Median lake CO2-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(np.median(lake_co2_17), unc_lake_co2_17))
print('Median lake CH4-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(np.median(lake_ch4_17), unc_lake_ch4_17))
print('Median lake C emissions (CO2-C + CH4-C) are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(median_lake_c_17, unc_lake_c_17))    
print('Lake C emissions (CO2-C + CH4-C) range from {:0.2f} to {:0.2f} (Mg C d-1) [95% CI]'.format(conf_int_lake_c_17[0], conf_int_lake_c_17[1]))  

# Fluvial CO2-C + CH4-C emissions
fluvial_co2_17 = (fluvial_co2.df_boot['flux co2-2 mg'] * fluvial_area_17) *10**-9
std_fluvial_co2_17  = np.std(fluvial_co2_17)
unc_fluvial_co2_17  = std_fluvial_co2_17 + math.sqrt((0.15*np.median(fluvial_co2_17))**2 + (0.15*(fluvial_area_17/1000000))**2)
fluvial_ch4_17 = (fluvial_ch4.df_boot['flux ch4-2 mg'] * fluvial_area_17) *10**-9
std_fluvial_ch4_17  = np.std(fluvial_ch4_17)
unc_fluvial_ch4_17  = std_fluvial_ch4_17 + math.sqrt((0.15*np.median(fluvial_ch4_17))**2 + (0.15*(fluvial_area_17/1000000))**2)
fluvial_c_17 = fluvial_co2_17 + fluvial_ch4_17
median_fluvial_c_17 = np.median(fluvial_c_17)
conf_int_fluvial_c_17 = np.percentile(fluvial_c_17, [5, 95])
std_fluvial_c_17 = np.std(fluvial_c_17)
unc_fluvial_c_17 = std_fluvial_c_17 + math.sqrt((0.15*np.median(fluvial_c_17))**2 + (0.15*(fluvial_area_17/1000000))**2)
print('Median fluvial CO2-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(np.median(fluvial_co2_17), unc_fluvial_co2_17))
print('Median fluvial CH4-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(np.median(fluvial_ch4_17), unc_fluvial_ch4_17))
print('Median fluvial C emissions (CO2-C + CH4-C) are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(median_fluvial_c_17, unc_fluvial_c_17))    
print('Fluvial C emissions (CO2-C + CH4-C) range from {:0.2f} to {:0.2f} (Mg C d-1) [95% CI]'.format(conf_int_fluvial_c_17[0], conf_int_fluvial_c_17[1])) 

# Flood CO2-C + CH4-C emissions
flood_co2_17 = (flood_co2.df_boot['flux co2-2 mg'] * flood_area_17) *10**-9
std_flood_co2_17  = np.std(flood_co2_17)
unc_flood_co2_17  = std_flood_co2_17 + math.sqrt((0.15*np.median(flood_co2_17))**2 + (0.15*(flood_area_17/1000000))**2)
flood_ch4_17 = (flood_ch4.df_boot['flux ch4-2 mg'] * flood_area_17) *10**-9
std_flood_ch4_17  = np.std(flood_ch4_17)
unc_flood_ch4_17  = std_flood_ch4_17 + math.sqrt((0.15*np.median(flood_ch4_17))**2 + (0.15*(flood_area_17/1000000))**2)
flood_c_17 = flood_co2_17 + flood_ch4_17
median_flood_c_17 = np.median(flood_c_17)
conf_int_flood_c_17 = np.percentile(flood_c_17, [5, 95])
std_flood_c_17 = np.std(flood_c_17)
unc_flood_c_17 = std_flood_c_17 + math.sqrt((0.15*median_flood_c_17)**2 + (0.15*(flood_area_17/1000000))**2)
print('Median flood CO2-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(np.median(flood_co2_17), unc_flood_co2_17))
print('Median flood CH4-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(np.median(flood_ch4_17), unc_flood_ch4_17))
print('Median flood C emissions (CO2-C + CH4-C) are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(median_flood_c_17, unc_flood_c_17))    
print('Flood C emissions (CO2-C + CH4-C) range from {:0.2f} to {:0.2f} (Mg C d-1) [95% CI]'.format(conf_int_flood_c_17[0], conf_int_flood_c_17[1]))  


# Total inland water emissions
co2_c_17 =  (lake_co2.df_boot['flux co2-2 mg'] * lake_area_17 + pond_co2.df_boot['flux co2-2 mg'] * pond_area_17 + \
             fluvial_co2.df_boot['flux co2-2 mg'] * fluvial_area_17 + flood_co2.df_boot['flux co2-2 mg'] * flood_area_17) * 10**-9
co2_c_17_median = np.median(co2_c_17)
ch4_c_17 = (lake_ch4.df_boot['flux ch4-2 mg'] * lake_area_17 + pond_ch4.df_boot['flux ch4-2 mg'] * pond_area_17 + \
            fluvial_ch4.df_boot['flux ch4-2 mg'] * fluvial_area_17 + flood_ch4.df_boot['flux ch4-2 mg'] * flood_area_17) * 10**-9
ch4_c_17_median = np.median(ch4_c_17)

emissions_c_17 = co2_c_17 + ch4_c_17

median_emissions_c_17 = np.median(emissions_c_17)
conf_int_emissions_c_17 = np.percentile(emissions_c_17, [5, 95])
unc_emissions_c_17 = math.sqrt(unc_pond_c_17**2 + unc_lake_c_17**2 + unc_fluvial_c_17**2 + unc_flood_c_17**2)
print('Median inland water C emissions (CO2-C + CH4-C) are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(median_emissions_c_17, unc_emissions_c_17))    
print('Inland water C emissions (CO2-C + CH4-C) range from {:0.2f} to {:0.2f} (Mg C d-1) [95% CI]'.format(conf_int_emissions_c_17[0], conf_int_emissions_c_17[1]))  

print('Inland water CO2 emissions accounted for {:2.2%} of the inland water CO2-C + CH4-C emissions in 2017'.format(co2_c_17_median / median_emissions_c_17))
print('Inland water CH4 emissions accounted for {:2.2%} of the inland water CO2-C + CH4-C emissions in 2017'.format(ch4_c_17_median / median_emissions_c_17))

print('Inland water lake emissions accounted for {:2.2%} of the inland water CO2-C + CH4-C emissions in 2017'.format(median_lake_c_17 / median_emissions_c_17))
print('Inland water pond emissions accounted for {:2.2%} of the inland water CO2-C + CH4-C emissions in 2017'.format(median_pond_c_17 / median_emissions_c_17))
print('Inland water fluvial emissions accounted for {:2.2%} of the inland water CO2-C + CH4-C emissions in 2017'.format(median_fluvial_c_17 / median_emissions_c_17))
print('Inland water food emissions accounted for {:2.2%} of the inland water CO2-C + CH4-C emissions in 2017'.format(median_flood_c_17 / median_emissions_c_17))

# CO2, CH4 and N2O (GHG) inland water emissions in CO2-eq in 2017 --------------------------------------------------------------------------------------------------------
# Lake GHG emissions in CO2-eq
lake_co2_eq_17 = (lake_co2.df_boot['flux co2'] + lake_ch4.df_boot['flux ch4 co2eq'] + lake_n2o.df_boot['flux n2o co2eq'])* lake_area_17 *10**-6
median_lake_co2_eq_17 = np.median(lake_co2_eq_17)

lake_co2_co2_eq_17 = lake_co2.df_boot['flux co2'] * lake_area_17 *10**-6
std_lake_co2_co2_eq_17 = np.std(lake_co2_eq_17)
unc_lake_co2_co2_eq_17 = std_lake_co2_co2_eq_17 + math.sqrt((0.15*np.median(lake_co2_eq_17))**2 + (0.15*(lake_area_17/1000000))**2)
lake_ch4_co2_eq_17 = lake_ch4.df_boot['flux ch4 co2eq']  * lake_area_17 *10**-6
std_lake_ch4_co2_eq_17 = np.std(lake_ch4_co2_eq_17)
unc_lake_ch4_co2_eq_17 = std_lake_ch4_co2_eq_17 + math.sqrt((0.15*np.median(lake_ch4_co2_eq_17))**2 + (0.15*(lake_area_17/1000000))**2)
lake_n2o_co2_eq_17 = lake_n2o.df_boot['flux n2o co2eq']* lake_area_17 *10**-6
std_lake_n2o_co2_eq_17 = np.std(lake_n2o_co2_eq_17)
unc_lake_n2o_co2_eq_17 = std_lake_n2o_co2_eq_17 + math.sqrt((0.15*np.median(lake_n2o_co2_eq_17))**2 + (0.15*(lake_area_17/1000000))**2)
median_lake_n2o_co2_eq_17 = np.median(lake_n2o_co2_eq_17)

conf_int_lake_co2_eq_17 = np.percentile(lake_co2_eq_17, [5, 95])
std_lake_co2_eq_17 = np.std(lake_co2_eq_17)
unc_lake_co2_eq_17 = std_lake_co2_eq_17 + math.sqrt((0.15*median_lake_co2_eq_17)**2 + (0.15*(lake_area_17/1000000))**2)
print('Median lake N2O emissions  are {:0.2f}  (kilo CO2-eq mole d-1) in 2017'.format(median_lake_n2o_co2_eq_17))
print('Median lake GHG emissions  are {:0.2f} ± {:0.2f} (kilo CO2-eq mole d-1)'.format(median_lake_co2_eq_17, unc_lake_co2_eq_17))    
print('Lake GHG emissions range from {:0.2f} to {:0.2f} (kilo CO2-eq mole d-1) [95% CI]'.format(conf_int_lake_co2_eq_17 [0], conf_int_lake_co2_eq_17 [1]))  

# Pond GHG emissions in CO2-eq
pond_co2_eq_17 = (pond_co2.df_boot['flux co2'] + pond_ch4.df_boot['flux ch4 co2eq'] + pond_n2o.df_boot['flux n2o co2eq'])* pond_area_17 *10**-6
median_pond_co2_eq_17 = np.median(pond_co2_eq_17)

pond_co2_co2_eq_17 = pond_co2.df_boot['flux co2'] * pond_area_17 *10**-6
std_pond_co2_co2_eq_17 = np.std(pond_co2_co2_eq_17)
unc_pond_co2_co2_eq_17 = std_pond_co2_co2_eq_17 + math.sqrt((0.15*np.median(pond_co2_co2_eq_17))**2 + (0.15*(pond_area_17/1000000))**2)
pond_ch4_co2_eq_17 = pond_ch4.df_boot['flux ch4 co2eq']  * pond_area_17 *10**-6
std_pond_ch4_co2_eq_17 = np.std(pond_ch4_co2_eq_17)
unc_pond_ch4_co2_eq_17 = std_pond_ch4_co2_eq_17 + math.sqrt((0.15*np.median(pond_ch4_co2_eq_17))**2 + (0.15*(pond_area_17/1000000))**2)
pond_n2o_co2_eq_17 = pond_n2o.df_boot['flux n2o co2eq']* pond_area_17 *10**-6
std_pond_n2o_co2_eq_17 = np.std(pond_n2o_co2_eq_17)
unc_pond_n2o_co2_eq_17 = std_pond_n2o_co2_eq_17 + math.sqrt((0.15*np.median(pond_n2o_co2_eq_17))**2 + (0.15*(pond_area_17/1000000))**2)
median_pond_n2o_co2_eq_17 = np.median(pond_n2o.df_boot['flux n2o co2eq']* pond_area_17 *10**-6)
conf_int_pond_co2_eq_17 = np.percentile(pond_co2_eq_17, [5, 95])

std_pond_co2_eq_17 = np.std(pond_co2_eq_17)
unc_pond_co2_eq_17 = std_pond_co2_eq_17 + math.sqrt((0.15*median_pond_co2_eq_17)**2 + (0.15*(pond_area_17/1000000))**2)
print('Median pond N2O emissions  are {:0.2f} (kilo CO2-eq mole d-1) in 2017'.format(median_pond_n2o_co2_eq_17))
print('Median pond GHG emissions  are {:0.2f} ± {:0.2f} (kilo CO2-eq mole d-1)'.format(median_pond_co2_eq_17, unc_pond_co2_eq_17))    
print('Pond GHG emissions range from {:0.2f} to {:0.2f} (kilo CO2-eq mole d-1) [95% CI]'.format(conf_int_pond_co2_eq_17 [0], conf_int_pond_co2_eq_17 [1]))  

# Fluvial GHG emissions in CO2-eq
fluvial_co2_eq_17 = (fluvial_co2.df_boot['flux co2'] + fluvial_ch4.df_boot['flux ch4 co2eq'] + fluvial_n2o.df_boot['flux n2o co2eq'])* fluvial_area_17 *10**-6
median_fluvial_co2_eq_17 = np.median(fluvial_co2_eq_17)

fluvial_co2_co2_eq_17 = fluvial_co2.df_boot['flux co2'] * fluvial_area_17 *10**-6
std_fluvial_co2_co2_eq_17 = np.std(fluvial_co2_co2_eq_17)
unc_fluvial_co2_co2_eq_17 = std_fluvial_co2_co2_eq_17 + math.sqrt((0.15*np.median(fluvial_co2_co2_eq_17))**2 + (0.15*(fluvial_area_17/1000000))**2)
fluvial_ch4_co2_eq_17 = fluvial_ch4.df_boot['flux ch4 co2eq']  * fluvial_area_17 *10**-6
std_fluvial_ch4_co2_eq_17 = np.std(fluvial_ch4_co2_eq_17)
unc_fluvial_ch4_co2_eq_17 = std_fluvial_ch4_co2_eq_17 + math.sqrt((0.15*np.median(fluvial_ch4_co2_eq_17))**2 + (0.15*(fluvial_area_17/1000000))**2)
fluvial_n2o_co2_eq_17 = fluvial_n2o.df_boot['flux n2o co2eq']* fluvial_area_17 *10**-6
std_fluvial_n2o_co2_eq_17 = np.std(fluvial_n2o_co2_eq_17)
unc_fluvial_n2o_co2_eq_17 = std_fluvial_n2o_co2_eq_17 + math.sqrt((0.15*np.median(fluvial_n2o_co2_eq_17))**2 + (0.15*(fluvial_area_17/1000000))**2)
median_fluvial_n2o_co2_eq_17 = np.median(fluvial_n2o.df_boot['flux n2o co2eq']* fluvial_area_17 *10**-6)
conf_int_fluvial_co2_eq_17 = np.percentile(fluvial_co2_eq_17, [5, 95])

std_fluvial_co2_eq_17 = np.std(fluvial_co2_eq_17)
unc_fluvial_co2_eq_17 = std_fluvial_co2_eq_17 + math.sqrt((0.15*median_fluvial_co2_eq_17)**2 + (0.15*(fluvial_area_17/1000000))**2)
print('Median fluvial N2O emissions  are {:0.2f}  (kilo CO2-eq mole d-1) in 2017'.format(median_fluvial_n2o_co2_eq_17))
print('Median fluvial GHG emissions  are {:0.2f} ± {:0.2f} (kilo CO2-eq mole d-1)'.format(median_fluvial_co2_eq_17, unc_fluvial_co2_eq_17))    
print('Fluvial GHG emissions range from {:0.2f} to {:0.2f} (kilo CO2-eq mole d-1) [95% CI]'.format(conf_int_fluvial_co2_eq_17 [0], conf_int_fluvial_co2_eq_17 [1]))  

# Flood GHG emissions in CO2-eq
flood_co2_eq_17 = (flood_co2.df_boot['flux co2'] + flood_ch4.df_boot['flux ch4 co2eq'] + flood_n2o.df_boot['flux n2o co2eq'])* flood_area_17 *10**-6
median_flood_co2_eq_17 = np.median(flood_co2_eq_17)

flood_co2_co2_eq_17 = flood_co2.df_boot['flux co2'] * flood_area_17 *10**-6
std_flood_co2_co2_eq_17 = np.std(flood_co2_co2_eq_17)
unc_flood_co2_co2_eq_17 = std_flood_co2_co2_eq_17 + math.sqrt((0.15*np.median(flood_co2_co2_eq_17))**2 + (0.15*(flood_area_17/1000000))**2)
flood_ch4_co2_eq_17 = flood_ch4.df_boot['flux ch4 co2eq']  * flood_area_17 *10**-6
std_flood_ch4_co2_eq_17 = np.std(flood_ch4_co2_eq_17)
unc_flood_ch4_co2_eq_17 = std_flood_ch4_co2_eq_17 + math.sqrt((0.15*np.median(flood_ch4_co2_eq_17))**2 + (0.15*(flood_area_17/1000000))**2)
flood_n2o_co2_eq_17 = flood_n2o.df_boot['flux n2o co2eq']* flood_area_17 *10**-6
std_flood_n2o_co2_eq_17 = np.std(flood_n2o_co2_eq_17)
unc_flood_n2o_co2_eq_17 = std_flood_n2o_co2_eq_17 + math.sqrt((0.15*np.median(flood_n2o_co2_eq_17))**2 + (0.15*(flood_area_17/1000000))**2)
median_flood_n2o_co2_eq_17 = np.median(flood_n2o.df_boot['flux n2o co2eq']* flood_area_17 *10**-6)
conf_int_flood_co2_eq_17 = np.percentile(flood_co2_eq_17, [5, 95])

std_flood_co2_eq_17 = np.std(flood_co2_eq_17)
unc_flood_co2_eq_17 = std_flood_co2_eq_17 + math.sqrt((0.15*median_flood_co2_eq_17)**2 + (0.15*(flood_area_17/1000000))**2)
print('Median flood N2O emissions  are {:0.2f}  (kilo CO2-eq mole d-1) in 2017'.format(median_flood_n2o_co2_eq_17))
print('Median flood GHG emissions  are {:0.2f} ± {:0.2f} (kilo CO2-eq mole d-1)'.format(median_flood_co2_eq_17, unc_flood_co2_eq_17))    
print('flood GHG emissions range from {:0.2f} to {:0.2f} (kilo CO2-eq mole d-1) [95% CI]'.format(conf_int_flood_co2_eq_17 [0], conf_int_flood_co2_eq_17 [1]))  



# Total inland water GHG emissions in CO2-eq
co2_eq_17 = (lake_co2.df_boot['flux co2'] * lake_area_17 + pond_co2.df_boot['flux co2'] * pond_area_17 + \
             fluvial_co2.df_boot['flux co2'] * fluvial_area_17 + flood_co2.df_boot['flux co2'] * flood_area_17) *10**-6
co2_eq_17_median = np.median(co2_eq_17)
co2_eq_17_unc = math.sqrt(unc_lake_co2_co2_eq_17**2 + unc_pond_co2_co2_eq_17**2 + unc_fluvial_co2_co2_eq_17**2 + unc_flood_co2_co2_eq_17**2)
ch4_eq_17 = (lake_ch4.df_boot['flux ch4 co2eq'] * lake_area_17 + pond_ch4.df_boot['flux ch4 co2eq'] * pond_area_17 + \
             fluvial_ch4.df_boot['flux ch4 co2eq'] * fluvial_area_17 + flood_ch4.df_boot['flux ch4 co2eq'] * flood_area_17) *10**-6
ch4_eq_17_median = np.median(ch4_eq_17)
ch4_eq_17_unc = math.sqrt(unc_lake_ch4_co2_eq_17**2 + unc_pond_ch4_co2_eq_17**2 + unc_fluvial_ch4_co2_eq_17**2 + unc_flood_ch4_co2_eq_17**2)
n2o_eq_17 = (lake_n2o.df_boot['flux n2o co2eq'] * lake_area_17 + pond_n2o.df_boot['flux n2o co2eq'] * pond_area_17 + \
             fluvial_n2o.df_boot['flux n2o co2eq'] * fluvial_area_17 + flood_n2o.df_boot['flux n2o co2eq'] * flood_area_17) *10**-6
n2o_eq_17_median = np.median(n2o_eq_17)
n2o_eq_17_unc = math.sqrt(unc_lake_n2o_co2_eq_17**2 + unc_pond_n2o_co2_eq_17**2 + unc_fluvial_n2o_co2_eq_17**2 + unc_flood_n2o_co2_eq_17**2)

print('Median CO2 inland water GHG emissions are {:0.2f} ± {:0.2f}  (kilo CO2-eq mole d-1) in 2017'.format(co2_eq_17_median, co2_eq_17_unc))
print('Median CH4 inland water GHG emissions are {:0.2f} ± {:0.2f}  (kilo CO2-eq mole d-1) in 2017'.format(ch4_eq_17_median, ch4_eq_17_unc))
print('Median N2O inland water GHG emissions are {:0.2f} ± {:0.2f}  (kilo CO2-eq mole d-1) in 2017'.format(n2o_eq_17_median, n2o_eq_17_unc))

emissions_ghg_17 = co2_eq_17 + ch4_eq_17 + n2o_eq_17
emissions_ghg_17_median = np.median(emissions_ghg_17)
conf_int_emissions_ghg_17 = np.percentile(emissions_ghg_17, [5, 95])

unc_emissions_ghg_17 = math.sqrt(unc_pond_co2_eq_17**2 + unc_lake_co2_eq_17**2 + unc_fluvial_co2_eq_17**2)
unc_emissions_ghg_17_perc = unc_emissions_ghg_17 / emissions_ghg_17_median

print('Inland water CO2 emissions accounted for {:2.2%} of the inland water GHG emissions in CO2-equivalents in 2017'.format(co2_eq_17_median / emissions_ghg_17_median))
print('Inland water CH4 emissions accounted for {:2.2%} of the inland water GHG emissions in CO2-equivalents in 2017'.format(ch4_eq_17_median / emissions_ghg_17_median))
print('Inland water N2O emissions accounted for {:2.2%} of the inland water GHG emissions in CO2-equivalents in 2017'.format(n2o_eq_17_median / emissions_ghg_17_median))
print('The uncertainty of the inland water GHG emissions in CO2-equivalents in 2017 is {:2.2%}'.format(unc_emissions_ghg_17_perc))

print('Median inland water GHG emissions (CO2-eq) are {:0.2f} ± {:0.2f} (kilo CO2-eq mole d-1)'.format(emissions_ghg_17_median, unc_emissions_ghg_17))    

# Read in tundra GHG emissions --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in tundra CH4 emssions 
fname = r'...\FLX_RU-Cok_FLUXNET-CH4_HH_2008-2016_1-1.csv'
startdate_17 = '2016-07-30'
enddate_17 = '2016-08-10'

# Tundra CH4 Fluxes in CO2-eq
t_ch4_eq = tundra_flux('CH4')
t_ch4_f_eq = t_ch4_eq.readin_fluxes_tch4(fname, startdate_17, enddate_17, 'GWP')
# Convert from mmol CO2-eq m-2 d-1 to kilo mol CO2-eq m-2 d-1
t_ch4_boot_eq = t_ch4_eq.bootstrapping(1000) *10**-6

# Tundra CH4 Fluxes in mg CH4-C m-2 d-1
t_ch4_c = tundra_flux('CH4')
t_ch4_f_c = t_ch4_c.readin_fluxes_tch4(fname, startdate_17, enddate_17, 'CH4-C')
t_ch4_boot_c = t_ch4_c.bootstrapping(1000)
t_ch4_boot_c_median = np.median(t_ch4_boot_c)

# Tundra CH4 Fluxes in Mg CH4-C m-2 d-1
t_ch4_boot_c = t_ch4_boot_c * 10**-9


# Read in tundra CO2 emssions 
fname = r'C:\Users\mmo990\surfdrive\Paper1\Data\EC_Tower_CO2.csv'
startdate_17 = '2016-07-30'
enddate_17 = '2016-08-10'

# Tundra CO2 Fluxes in CO2-eq
t_co2_eq = tundra_flux('CO2')
t_co2_f_eq = t_co2_eq.readin_fluxes_tco2(fname, startdate_17, enddate_17, 'GWP')
# Convert from mmol CO2 m-2 d-1 to kilo mol CO2 m-2 d-1
t_co2_boot_eq = t_co2_eq.bootstrapping(1000) *10**-6

# Tundra CO2 Fluxes in mg CO2-C m-2 d-1
t_co2_c = tundra_flux('CO2')
t_co2_f_c = t_co2_c.readin_fluxes_tco2(fname, startdate_17, enddate_17, 'CO2-C')
t_co2_boot_c = t_co2_c.bootstrapping(1000)
t_co2_boot_c_median = np.median(t_co2_boot_c)

# Tundra CO2 Fluxes in Mg CO2-C m-2 d-1
t_co2_boot_c = t_co2_boot_c * 10**-9


# Tundra N2O flux
fname = r'...\n2o_tundra_fluxes.xlsx'

t_n2o = tundra_flux('N2O')
t_n2o.readin_fluxes_tn2o(fname)
t_n2o_boot = t_n2o.bootstrapping(1000)
# Convert from mmol CO2-eq m-2 d-1 to kilo mol CO2-eq m-2 d-1
t_n2o_boot = t_n2o_boot * 10**-6

# Total tundra emissions in CO2-eq in kilo mol CO2-eq d-1
t_co2eq = (t_co2_boot_eq + t_ch4_boot_eq + t_n2o_boot) * tundra_area_17
t_co2eq_median = np.median(t_co2eq)
std_landscape_tundra_co2eq_17 = np.std(t_co2eq)
unc_tundra_co2eq_17 = std_landscape_tundra_co2eq_17 + math.sqrt((0.15*t_co2eq_median)**2 + (0.15*(tundra_area_17/1000000))**2)
print('Median tundra GHG emissions are {:0.2f} ± {:0.2f} (kilo CO2-eq d-1)'.format(t_co2eq_median, unc_tundra_co2eq_17))

# C (CO2-C + CH4-C) landscape exchange mg C m-2 d-1---------------------------------------------------------------------------------------------------------------------------------------------------------
tundra_c_ch4_17 = t_ch4_boot_c * tundra_area_17
tundra_c_ch4_17_median = np.median(tundra_c_ch4_17)
tundra_c_ch4_17_std = np.std(tundra_c_ch4_17)
tundra_c_ch4_17_unc = tundra_c_ch4_17_std + math.sqrt((0.15*tundra_c_ch4_17_median)**2 + (0.15*(tundra_area_17/1000000))**2)
print('Median tundra CH4-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(tundra_c_ch4_17_median, tundra_c_ch4_17_unc))
tundra_c_co2_17 = t_co2_boot_c * tundra_area_17
tundra_c_co2_17_median = np.median(tundra_c_co2_17)
tundra_c_co2_17_std = np.std(tundra_c_co2_17)
tundra_c_co2_17_unc = tundra_c_co2_17_std + math.sqrt((0.15*tundra_c_co2_17_median)**2 + (0.15*(tundra_area_17/1000000))**2)
print('Median tundra CO2-C emissions are {:0.2f} ± {:0.2f} (Mg C d-1)'.format(tundra_c_co2_17_median, tundra_c_co2_17_unc))

tundra_c_17 = tundra_c_ch4_17 + tundra_c_co2_17
tundra_c_17_median = np.median(tundra_c_ch4_17 + tundra_c_co2_17)
std_landscape_tundra_c_17 = np.std(tundra_c_17)
unc_tundra_c_17 = std_landscape_tundra_c_17 + math.sqrt((0.15*np.median(tundra_c_17))**2 + (0.15*(tundra_area_17/1000000))**2)

landscape_c_exchange_17 = emissions_c_17 + tundra_c_17
landscape_c_exchange_17_median = np.median(landscape_c_exchange_17)
conf_int_landscape_c_exchange_17 = np.percentile(landscape_c_exchange_17, [5, 95])
unc_landscape_c_exchange_17 = math.sqrt(unc_pond_c_17**2 + unc_lake_c_17**2 + unc_fluvial_c_17**2 + unc_tundra_c_17**2)
unc_landscape_c_exchange_17_perc = unc_landscape_c_exchange_17/landscape_c_exchange_17_median

print('Median landscape C exchange (CO2-C + CH4-C) is {:0.2f} ± {:0.2f} (Mg C d-1)'.format(landscape_c_exchange_17_median, unc_landscape_c_exchange_17))    
print('Median landscape C exchange (CO2-C + CH4-C) range from {:0.2f} to {:0.2f} (Mg C d-1) [95% CI]'.format(conf_int_landscape_c_exchange_17[0], conf_int_landscape_c_exchange_17[1]))
print('The uncertainty of landscape C exchange in 2017 is {:2.2%}'.format(unc_landscape_c_exchange_17_perc))


print('Inland water lake emissions accounted for {:2.2%} of the landscape CO2-C + CH4-C exchange in 2017'.format(median_lake_c_17/landscape_c_exchange_17_median))
print('Inland water pond emissions accounted for {:2.2%} of the landscape CO2-C + CH4-C exchange in 2017'.format(median_pond_c_17/landscape_c_exchange_17_median))
print('Inland water fluvial emissions accounted for {:2.2%} of the landscape CO2-C + CH4-C exchange in 2017'.format(median_fluvial_c_17/landscape_c_exchange_17_median))


print('Inland water C emissions offset the landscape C sink by {:2.2%} in 2017'\
      .format((tundra_c_17_median - (tundra_c_17_median + median_emissions_c_17)) / tundra_c_17_median))


# CO2 + CH4 + N20 landscape exchange CO2-eq ------------------------------------------------------------------------------------------------------------
tundra_co2_co2eq_17 = t_co2_boot_eq * tundra_area_17
median_tundra_co2_co2eq_17 = np.median(tundra_co2_co2eq_17)
std_tundra_co2_co2eq_17 = np.std(tundra_co2_co2eq_17)
unc_tundra_co2_co2eq_17 = std_tundra_co2_co2eq_17 + math.sqrt((0.15*median_tundra_co2_co2eq_17)**2 + (0.15*(tundra_area_17/1000000))**2)

tundra_ch4_co2eq_17 = t_ch4_boot_eq * tundra_area_17
median_tundra_ch4_co2eq_17 = np.median(tundra_ch4_co2eq_17)
std_tundra_ch4_co2eq_17 = np.std(tundra_ch4_co2eq_17)
unc_tundra_ch4_co2eq_17 = std_tundra_ch4_co2eq_17 + math.sqrt((0.15*median_tundra_ch4_co2eq_17)**2 + (0.15*(tundra_area_17/1000000))**2)

tundra_n2o_co2eq_17 = t_n2o_boot * tundra_area_17
median_tundra_n2o_co2eq_17 = np.median(tundra_n2o_co2eq_17)
std_tundra_n2o_co2eq_17 = np.std(tundra_n2o_co2eq_17)
unc_tundra_n2o_co2eq_17 = std_tundra_n2o_co2eq_17 + math.sqrt((0.15*median_tundra_n2o_co2eq_17)**2 + (0.15*(tundra_area_17/1000000))**2)


landscape_co2eq_exchange_17 = emissions_ghg_17 + t_co2eq
landscape_co2eq_exchange_17_median = np.median(landscape_co2eq_exchange_17)
conf_int_landscape_co2eq_exchange_17 = np.percentile(landscape_co2eq_exchange_17, [5, 95])
unc_landscape_co2eq_exchange_17 = math.sqrt(unc_pond_co2_eq_17**2 + unc_lake_co2_eq_17**2 + unc_fluvial_co2_eq_17**2 + unc_tundra_co2eq_17**2)

print('Median tundra CO2 GHG exchange in 2017 is {:0.2f} (kilo mole CO2-eq d-1)'.format(median_tundra_co2_co2eq_17))   
print('Median tundra CH4 GHG exchange in 2017 is {:0.2f} (kilo mole CO2-eq d-1)'.format(median_tundra_ch4_co2eq_17))     
print('Median tundra N2O GHG exchange in 2017 is {:0.2f} (kilo mole CO2-eq d-1)'.format(median_tundra_n2o_co2eq_17))

print('Median landscape GHG exchange in 2017 is {:0.2f} ± {:0.2f} (kilo mole CO2-eq d-1)'.format(landscape_co2eq_exchange_17_median, unc_landscape_co2eq_exchange_17))    
print('Median landscape GHG exchange in 2017 range from {:0.2f} to {:0.2f} (kilo mole CO2-eq d-1) [95% CI]'.format(conf_int_landscape_co2eq_exchange_17[0], conf_int_landscape_co2eq_exchange_17[1]))

print('Tundra CH4 emissions offset tundra CO2 sink in terms of CO2 equivalent by {:2.2f}'\
      .format((median_tundra_co2_co2eq_17 - (median_tundra_co2_co2eq_17 + median_tundra_ch4_co2eq_17 )) / median_tundra_co2_co2eq_17))
print('Tundra CH4 emissions offset tundra CO2 sink in terms of CO2 equivalent by {:2.2%}'\
      .format((median_tundra_co2_co2eq_17 - (median_tundra_co2_co2eq_17 + median_tundra_ch4_co2eq_17 )) / median_tundra_co2_co2eq_17))

print('Inland water emissions accounted for {:2.2%} of the landscape GHG exchange in 2017'.format(emissions_ghg_17_median/landscape_co2eq_exchange_17_median))

print('Inland water lake emissions accounted for {:2.2%} of the landscape GHG exchange in 2017'.format(median_lake_co2_eq_17/landscape_co2eq_exchange_17_median))
print('Inland water pond emissions accounted for {:2.2%} of the landscape GHG exchange in 2017'.format(median_pond_co2_eq_17/landscape_co2eq_exchange_17_median))
print('Inland water fluvial emissions accounted for {:2.2%} of the landscape GHG exchange in 2017'.format(median_fluvial_co2_eq_17/landscape_co2eq_exchange_17_median))

print('Inland water CO2 emissions accounted for {:2.2%} of the landscape GHG exchange in CO2-equivalents in 2017'.format(co2_eq_17_median/landscape_co2eq_exchange_17_median))
print('Inland water CH4 emissions accounted for {:2.2%} of the landscape GHG exchange in CO2-equivalents in 2017'.format(ch4_eq_17_median/landscape_co2eq_exchange_17_median))
print('Inland water N2O emissions accounted for {:2.2%} of the landscape GHG exchange in CO2-equivalents in 2017'.format(n2o_eq_17_median/landscape_co2eq_exchange_17_median))

print('CO2 emissions accounted for {:2.2%} of the landscape GHG exchange in CO2-equivalents in 2017'.format((median_tundra_co2_co2eq_17 + co2_eq_17_median) / landscape_co2eq_exchange_17_median))
print('CH4 emissions accounted for {:2.2%} of the landscape GHG exchange in CO2-equivalents in 2017'.format((median_tundra_ch4_co2eq_17 + ch4_eq_17_median) / landscape_co2eq_exchange_17_median))
print('N2O emissions accounted for {:2.2%} of the landscape GHG exchange in CO2-equivalents in 2017'.format((median_tundra_n2o_co2eq_17 + n2o_eq_17_median) /landscape_co2eq_exchange_17_median))


# Plots ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot results
barWidth = 0.3

wbs = [ 'Tundra','Inland waters']
co2 = [ median_tundra_co2_co2eq_17, co2_eq_17_median]
ch4 = [ median_tundra_ch4_co2eq_17, ch4_eq_17_median]
n2o= [ median_tundra_n2o_co2eq_17, n2o_eq_17_median]


error_CO2 = [(unc_tundra_co2_co2eq_17,0), (0, co2_eq_17_unc)] 
error_CH4 = [(0,0), (unc_tundra_ch4_co2eq_17, ch4_eq_17_unc)] 
error_N2O = [(0,0), (unc_tundra_n2o_co2eq_17, n2o_eq_17_unc)] 

# Set position of bar on X axis
r1 = np.arange(len(co2))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
 
# Make the plot
fig, ax = plt.subplots(figsize=(7,8))

plt.bar(r1, co2, color='#ad936b', width=barWidth, edgecolor='black', label='CO2', yerr=error_CO2, capsize=3)
plt.bar(r2, ch4, color='#935b73', width=barWidth, edgecolor='black', label='CH4', yerr=error_CH4, capsize=3)
plt.bar(r3, n2o, color='#859f62', width=barWidth, edgecolor='black', label='N2O', yerr=error_N2O, capsize=3)

 
# Add xticks on the middle of the group bars
plt.ylabel('$\\bf{up-scaled \ GHG \ emissions}$' + '\n' + '$\ (kilo \ mole \ CO_2-equivalent \ d^{-1})$', fontsize = 14) 
plt.xticks([r + barWidth*0.5 for r in range(len(co2))], ['Tundra', 'Inland waters'])
#plt.yscale('log')
color_name = "black"
ax.spines["top"].set_color(color_name)
ax.spines["bottom"].set_color(color_name)
ax.spines["left"].set_color(color_name)
ax.spines["right"].set_color(color_name)
ax.patch.set_facecolor('white')

#locmaj = matplotlib.ticker.LogLocator(base=10,numticks=6) 
#ax.yaxis.set_major_locator(locmaj)
#locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=6)
#ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())

# Create legend & Show graphic
plt.legend(loc='best', bbox_to_anchor=(0.5, 0.5, 0.5, 0.5), facecolor='white')
plt.show()

  
# Plot results
barWidth = 0.4

wbs = [ 'Pond','Lake', 'Fluvial', 'Flood']
co2 = [ np.mean(pond_co2_17), np.mean(lake_co2_17), np.mean(fluvial_co2_17), np.mean(flood_co2_17)]
ch4 = [ np.mean(pond_ch4_17), np.mean(lake_ch4_17), np.mean(fluvial_ch4_17), np.mean(flood_ch4_17)]


error_CO2 = [(0,0,0), (unc_pond_co2_17, unc_lake_co2_17, unc_fluvial_co2_17, unc_flood_co2_17)] 
error_CH4 = [(0,0,0), (unc_pond_ch4_17, unc_lake_ch4_17, unc_fluvial_ch4_17, unc_flood_ch4_17)] 

# Set position of bar on X axis
r1 = np.arange(len(co2))
r2 = [x + barWidth for x in r1]
 
# Make the plot
fig, ax = plt.subplots(figsize=(7,8))
#kwargs_co2 = {"hatch":'/'} 
#kwargs_ch4 = {"hatch":'.'} 
plt.bar(r1, co2, color='#ad936b', width=barWidth, edgecolor='black', label='CO2', yerr=error_CO2, capsize=3)
plt.bar(r2, ch4, color='#935b73', width=barWidth, edgecolor='black', label='CH4', yerr=error_CH4, capsize=3)

 
# Add xticks on the middle of the group bars
plt.ylabel('$\\bf{up-scaled \ C \ emissions\ (Mg\ C\ d^{-1})}$', fontsize = 14) 
plt.xticks([r + barWidth*0.5 for r in range(len(co2))], ['Pond \n (10%)', 'Lake \n (5%)' , 'Fluvial \n (2%)', 'Flood \n (10%)'])
plt.yscale('log')
color_name = "black"
ax.spines["top"].set_color(color_name)
ax.spines["bottom"].set_color(color_name)
ax.spines["left"].set_color(color_name)
ax.spines["right"].set_color(color_name)
ax.patch.set_facecolor('white')

locmaj = matplotlib.ticker.LogLocator(base=10,numticks=6) 
ax.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=6)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
plt.ylim([-10, 10])


# Create legend & Show graphic
plt.legend(loc='best', bbox_to_anchor=(0.5, 0.5, 0.5, 0.5), facecolor='white')
plt.show()

        
