# -*- coding: utf-8 -*-
"""
The following code is part of the manuscript: 
Written by M. Martyn Rosco

Functions for bootstrapping are from or based on Statistical Thinking in Python (Part 2) from https://www.datacamp.com
"""

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
 


class Fluxes:
    def __init__(self, gas, gas_air, wb_type, year, df=0, df_k600=0, k600=0, k=0, bs_replicates=0, df_boot=0, flux_std=0):
        self.gas = gas
        self.gas_air = gas_air
        self.wb_type = wb_type
        self.year = year
        
    def subset_df(self, my_df):
        self.df = my_df[my_df['waterbody'].str.contains(self.wb_type)]
        self.df = self.df[self.df['year'].str.contains(self.year)]
        self.df = self.df[[self.gas, self.gas_air, 'Date', 'Area', 'T', 'T - K', 'pH', 'EC', 'Salinity', 'year', 'wind speed', 'wind direction']]
        self.df['Date'] = pd.to_datetime(self.df['Date']).dt.date
        
    def outlier_removal(self, ll, ul, fpath_outliers):
        """
        outlier_removal removes outliers based on the 1.5IQR rule. 

        :param ll: lower limit of IQR
        :param ul: upper limit of IQR   
        :param fpath_outliers: filepath and name where csv file of removed outliers is saved
        :return: dataframe without outliers based on 1.5IQR
        """ 
    
        Q1 = self.df[self.gas].quantile(ll)
        Q3 = self.df[self.gas].quantile(ul)
        IQR = Q3 - Q1
        min_limit = (Q1 - 1.5*IQR)
        max_limit = (Q3 + 1.5*IQR)
        df_out_removed = self.df.loc[(self.df[self.gas] > min_limit) & (self.df[self.gas] < max_limit)]
        df_outliers = self.df[~self.df.astype(str).apply(tuple, 1).isin(df_out_removed.astype(str).apply(tuple, 1))]
        df_outliers.to_csv(fpath_outliers)
        self.df = df_out_removed

    def CO2_conc(self):  
        """ Converts dissolved CO2 concentration from ppm to micromol/L """ 

        self.df['diss co2'] = self.df['diss CO2.ppm'] * (np.exp(-58.0931+90.5069*(100/self.df['T - K'])+22.294*np.log(self.df['T - K']/100)+self.df['Salinity']*(0.027766-0.025888*(self.df['T - K']/100)+0.0050578*np.square(self.df['T - K']/100))))      
        return self.df
    
    def CO2_eq(self):
        """ Calculates equivalent CO2 dissolved concentration of atmospheric values using Weiss, 1974 formula """       
        
        self.df['CO2 eq'] = self.df['CO2 air']*(np.exp(-58.0931+90.5069*(100/self.df['T - K'])+22.294*np.log(self.df['T - K']/100)+self.df['Salinity']*(0.027766-0.025888*(self.df['T - K']/100)+0.0050578*np.square(self.df['T - K']/100))))
        return self.df
        
    def CH4_eq(self):
        """ Calculates equivalent CH4 dissolved concentration of atmospheric values using Wiesenburg & Guinasso Jr, 1979 """         
        A1 = -415.2807
        A2 = 596.8104
        A3 = 379.2599
        A4 = -62.0757
        B1 = -0.059160
        B2 = 0.032174
        B3 = -0.0048198
        self.df['CH4 eq'] = np.exp(np.log(self.df['CH4 air']*0.000001) + A1 + A2*(100/self.df['T - K']) + A3*np.log(self.df['T - K']/100) + A4*(self.df['T - K']/100) + self.df['Salinity']*(B1 + B2*(self.df['T - K']/100) + B3*(np.square(self.df['T - K']/100))))*0.001

    def N2O_eq(self):
        """  Calculates equivalent N2O dissolved concentration of atmospheric values using Weiss & Price, 1980 """   
        
        A1 = -165.8806
        A2 = 222.8743
        A3 = 92.0792
        A4 = -1.48425
        B1 = -0.056235
        B2 = 0.031619
        B3 = 0.0048472
        self.df['N2O eq'] = self.df['N2O air']*np.exp(A1 + A2*(100/self.df['T - K']) + A3*np.log(self.df['T - K']/100) + A4*np.square(self.df['T - K']/100) + self.df['Salinity']*(B1 + B2*(self.df['T - K']/100) + B3*(np.square(self.df['T - K']/100))))
                             
    def Sc(self):
        """ Calculates Schmidt Number (Sc) based on polynomial from Wanninkhof, 2014 """ 
        
        if self.gas == 'CO2':
            self.df['Sc freshwater'] = 1923.6 - 125.06 * self.df['T'] + 4.3773*np.power(self.df['T'], 2) - 0.085681*np.power(self.df['T'], 3) + 0.00070284*np.power(self.df['T'], 4)
            self.df['Sc saltwater'] = 2116.8 - 136.25 * self.df['T'] + 4.7353 * np.power(self.df['T'], 2) - 0.092307 * np.power(self.df['T'], 3) + 0.0007555 * np.power(self.df['T'], 4)
        elif self.gas == 'CH4':
            self.df['Sc freshwater'] = 1909.4 - 120.78 * self.df['T'] + 4.1555 * np.power(self.df['T'], 2) - 0.080578 * np.power(self.df['T'], 3) + 0.00065777 * np.power(self.df['T'], 4)
            self.df['Sc saltwater'] = 2101.2 - 131.54 * self.df['T'] + 4.4931 * np.power(self.df['T'], 2) - 0.08676 * np.power(self.df['T'], 3) + 0.00070663 * np.power(self.df['T'], 4)
        else:
            self.df['Sc freshwater'] = 2141.2 - 152.56 * self.df['T'] + 5.8963 * np.power(self.df['T'], 2) - 0.12411 * np.power(self.df['T'], 3) + 0.0010655 * np.power(self.df['T'], 4)
            self.df['Sc saltwater'] = 2356.2 - 166.38 * self.df['T'] + 6.3952 * np.power(self.df['T'], 2) - 0.13422 * np.power(self.df['T'], 3) + 0.0011506 * np.power(self.df['T'], 4)
        
        self.df['Sc'] = self.df['Sc freshwater'] + ((self.df['Sc saltwater'] - self.df['Sc freshwater']) / 35) * self.df['Salinity']
        
        return self.df 
    
    
    def k600_cc(self):
        """
        Calculates k600 [cm/h] using relationship with wind from Cole & Caraco, 1998
        Only used for lake and flood data
 
        """                 
      
        self.df['k600'] = 2.07 + 0.215 * np.power(self.df['wind speed'], 1.7) 
       
        return self.df


    def uniform_k600_lit(self, fname, dis_n):
        """
        Monte Carlo uncertainties with uniform distribution for k600 values from literature
        these are then randomly assigned to each data entry (sample)
        
        :fname: excel file with literature k600 values to be used
        :dis_n: uniform distribution size generated
        :return: dataframe with randomly assigned generated k600 values to fluvial sites
        
        """ 
        
        self.df_k600 = pd.read_excel(fname, index_col=0, engine='openpyxl', skiprows=1)
        self.df_k600.reset_index(inplace=True)
        
        
        np.random.seed(42)
        df_k = self.df_k600[self.df_k600['waterbody'].str.contains(self.wb_type)]
        k_600 = df_k['k600']
        
        # Draw out a uniform distribution with k values
        self.k600 = np.random.uniform(low=k_600.min(), high=k_600.max(), size=dis_n)
        
        # Plot the PDF and label axes
        plt.hist(self.k600, bins=50, density=True, histtype='step')
        plt.xlabel(self.wb_type + ' k600 values [cm/hr]')
        plt.ylabel('CDF')
        
        # Show plot
        plt.show()
        
        # Create an ECDF from real data: x, y
        data = k_600.to_numpy()
        n = len(data)
        
        x = np.sort(data)
        y = np.arange(1, n+1) / n
              
        # Create a CDF from theoretical samples: x_theor, y_theor
        n_theor = len(self.k600)
        
        x_theor = np.sort(self.k600)
        y_theor = np.arange(1, n_theor+1) / n_theor
               
        # Overlay the plots
        plt.plot(x_theor, y_theor)
        plt.plot(x, y, marker='.', linestyle='none')
        
        plt.xlabel(self.wb_type + ' k600 values [cm/hr]')
        plt.ylabel('CDF')
        
        plt.show()        
        
        # Randomly assign k600 values to each data entry / sample
        self.df = self.df.assign(k600=np.random.choice(self.k600, size=len(self.df.index)))

        return self.df   


    def k(self):
        """Calculates k using k600 and Sc""" 

        for i in range(len(self.df)):
            row = self.df.iloc[i]
            if row['wind speed']<3:
                self.df['k'] = (self.df['k600']*np.power((self.df['Sc']/600), 0.67)) / (100/24)
            else:
                self.df['k'] = (self.df['k600']*np.power((self.df['Sc']/600), 0.5)) / (100/24)
        
        return self.df

    
    def bootstrap_replicate_1d(self, data, func):
        """ Generate bootstrap replicate of 1D data. """
        
        bs_sample = np.random.choice(data, len(data))
        
        return func(bs_sample)
    
    
    def draw_bs_reps(self, data, func, size):
        """  Draw bootstrap replicates  """ 
        
        # Initialise array or replicates
        bs_replicates = np.empty(size)
        
        # Generate replicates
        for i in range(size):
            bs_replicates[i] = self.bootstrap_replicate_1d(data, func)
            
        return bs_replicates

    
    
    def bootstrapping(self, varble, boot_n): 
        """
        Generates bootstrapped median samples 

        :param varble: variable to select from dataframe and bootstrap
        :param boot_n: number of samples to generate
        :return: array of boostrapped median samples

        """ 
        
        gas_conc = self.df[varble].to_numpy()
        
        bs_replicates = self.draw_bs_reps(gas_conc, np.median, boot_n)
        
        # Make a histogram of the results
        plt.hist(bs_replicates, bins=50, density=True)
        plt.xlabel('median ' + varble)
        plt.ylabel('PDF')
        
        plt.show()
        
        # Compute the 95% confidence interval: conf_int
        conf_int = np.percentile(bs_replicates, [5, 95])
        print('95% confidence interval of ', self.wb_type, varble, conf_int) 
        
        return bs_replicates
    
    
    def fluxes_calc(self, cs, ceq, ks, cols, flux_name, flux_c_name, flux_c_name_mg, flux_gwp, flux_sgcp):
        """"
        Calculates fluxes using bootstrapped medians
        Flux = k x (Csurf - Ceq)
        Then also converts fluxes from (mmol/m2*d) to (mmol C/m2*d) and (mg C/m2*d)
        CO2 and CH4 fluxes are in (mmol/m2*d), N2O fluxes are in (nmol/m2*d)
        
        :param cs: array of bootstrapped median dissolved gas concentrations
        :param ceq: array of bootstrapped median atmosphere gas concentrations in equilibrium with dissolved concentration   
        :param ks: array of bootstrapped median gas exchange coefficients
        :param cols: list of strings of column names of new dataframe made wtih cs, ceq and k arrays
        :param flux_name: string name (and column name) of flux being calculatres
        :param flux_c_name: string name (and column name) of flux in (mmol C/m2*d), note only applicable to CO2 and CH4 fluxes
        :param flux_c_name_mg: string name (and column name) of flux in ((mg C/m2*d), note only applicable to CO2 and CH4 fluxes
        :return: dataframe without outliers based on 1.5IQR
        """
        
        co2_to_c = 0.27
        ch4_to_c = 0.75
        
        # in g/mol
        co2_atomic_mass = 44.01
        ch4_atomic_mass = 16.04
        
        boot_data = {cols[0]: cs, cols[1]: ceq, cols[2]: ks}
        self.df_boot = pd.DataFrame(boot_data)
        self.df_boot[flux_name] = (self.df_boot[cols[0]] - self.df_boot[cols[1]]) * self.df_boot[cols[2]]
        
        # Compute the 95% confidence interval: conf_int
        flux_median = self.df_boot[flux_name].median()
        conf_int = self.df_boot[flux_name].quantile([0.05, 0.95])
        std_err = self.df_boot[flux_name].std()
        if self.gas == 'diss N2O':
            print('The median of {} in {} is {:0.2f} (μmol/m2*d)'.format(self.gas, self.wb_type, flux_median)) 
            print('The standard error of {} in {} is {:0.2f} (μmol/m2*d)'.format(self.gas, self.wb_type, std_err)) 
            print('95% confidence interval of {} in {} is {:0.2f} to {:0.2f} (μmol/m2*d)'.format(self.gas, self.wb_type, conf_int.iloc[0], conf_int.iloc[1])) 
        else:    
            print('The median of {} in {} is {:0.2f} (mmol/m2*d)'.format(self.gas, self.wb_type, flux_median)) 
            print('The standard error of {} in {} is {:0.2f} (mmol/m2*d)'.format(self.gas, self.wb_type, std_err)) 
            print('95% confidence interval of {} in {} is {:0.2f} to {:0.2f} (mmol/m2*d)'.format(self.gas, self.wb_type, conf_int.iloc[0], conf_int.iloc[1])) 
        
        # Convert flux from (mmol/m2*d) to (mmol C/m2*d)
        # Note this is only applicable to CO2 and CH4 fluxes
        if self.gas == 'diss CO2.ppm':
            self.df_boot[flux_c_name] = self.df_boot[flux_name] * co2_to_c
            print('The median of {} in {} is {:0.2f} (mmol C/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_c_name].median())) 
            print('The standard error of {} in {} is {:0.2f} (mmol C/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_c_name].std()))
            ci_c =  self.df_boot[flux_c_name].quantile([0.05, 0.95])
            print('95% confidence interval of {} in {} is {:0.2f} to {:0.2f} (mmol C/m2*d)'.format(self.gas, self.wb_type, ci_c.iloc[0], ci_c.iloc[1])) 
        elif self.gas == 'diss CH4':
            self.df_boot[flux_c_name] = self.df_boot[flux_name] * ch4_to_c
            print('The median of {} in {} is {:0.2f} (mmol C/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_c_name].median())) 
            print('The standard error of {} in {} is {:0.2f} (mmol C/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_c_name].std())) 
            ci_c =  self.df_boot[flux_c_name].quantile([0.05, 0.95])
            print('95% confidence interval of {} in {} is {:0.2f} to {:0.2f} (mmol C/m2*d)'.format(self.gas, self.wb_type, ci_c.iloc[0], ci_c.iloc[1])) 
        
        # Convert flux from (mmol C/m2*d) to (mg C/m2*d)
        # Note this is only applicable to CO2 and CH4 fluxes
        
        if self.gas == 'diss CO2.ppm':
            self.df_boot[flux_c_name_mg] = self.df_boot[flux_c_name] * co2_atomic_mass 
            self.flux_std = self.df_boot[flux_c_name_mg].std()
            print('The median of {} in {} is {:0.2f} (mg C/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_c_name_mg].median())) 
            print('The standard error of {} in {} is {:0.2f} (mg C/m2*d)'.format(self.gas, self.wb_type, self.flux_std)) 
            ci_mg =  self.df_boot[flux_c_name_mg].quantile([0.05, 0.95])
            print('95% confidence interval of {} in {} is {:0.2f} to {:0.2f} (mg C/m2*d)'.format(self.gas, self.wb_type, ci_mg.iloc[0], ci_mg.iloc[1])) 
        elif self.gas == 'diss CH4':
            self.df_boot[flux_c_name_mg] = self.df_boot[flux_c_name] * ch4_atomic_mass
            self.flux_std = self.df_boot[flux_c_name_mg].std()
            print('The median of {} in {} is {:0.2f} (mg C/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_c_name_mg].median())) 
            print('The standard error of {} in {} is {:0.2f} (mg C/m2*d)'.format(self.gas, self.wb_type, self.flux_std)) 
            ci_mg =  self.df_boot[flux_c_name_mg].quantile([0.05, 0.95])
            print('95% confidence interval of {} in {} is {:0.2f} to {:0.2f} (mg C/m2*d)'.format(self.gas, self.wb_type, ci_mg.iloc[0], ci_mg.iloc[1])) 
            
       # Convert CH4 and N2O to CO2-equivalent using SGWP and SGCP over 100-year timeframe (Neubauer & Megonigal, 2015)
        ch4_sgwp = 45
        ch4_sgcp = 203
        n2o_sgwp = 270
        n2o_sgcp = 349
        if self.gas == 'diss CH4':
                self.df_boot[flux_gwp] = (self.df_boot[flux_name]*ch4_sgwp).where(self.df_boot[flux_name]>0)
                print('The median of {} in {} is {:0.2f} (mmol CO2-eq/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_gwp].median())) 
                print('The standard error of {} in {} is {:0.2f} (mmol CO2-eq/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_gwp].std())) 
                ci_co2eq =  self.df_boot[flux_gwp].quantile([0.05, 0.95])
                print('95% confidence interval of {} in {} is {:0.2f} to {:0.2f} (mmol CO2-eq/m2*d)'.format(self.gas, self.wb_type, ci_co2eq.iloc[0], ci_co2eq.iloc[1])) 
                #self.df_boot[flux_gwp] = (self.df_boot[flux_name]*ch4_sgcp).where(self.df_boot[flux_name]<0)
        elif self.gas == 'diss N2O':
            # Also convert N2O flux from micromol/L to mmol/L
                self.df_boot[flux_gwp] = (self.df_boot[flux_name]*10**-3*n2o_sgwp).where(self.df_boot[flux_name]>0)
                self.df_boot[flux_sgcp] = (self.df_boot[flux_name]*10**-3*n2o_sgcp).where(self.df_boot[flux_name]<0)      
                self.df_boot[flux_gwp].fillna(self.df_boot[flux_sgcp], inplace=True)
                print('The median of {} in {} is {:0.2f} (mmol CO2-eq/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_gwp].median())) 
                print('The standard error of {} in {} is {:0.2f} (mmol CO2-eq/m2*d)'.format(self.gas, self.wb_type, self.df_boot[flux_gwp].std())) 
                ci_co2eq =  self.df_boot[flux_gwp].quantile([0.05, 0.95])
                print('95% confidence interval of {} in {} is {:0.2f} to {:0.2f} (mmol CO2-eq/m2*d)'.format(self.gas, self.wb_type, ci_co2eq.iloc[0], ci_co2eq.iloc[1])) 
                 
        return self.df_boot
    
    
    



class tundra_flux:
    def __init__(self, gas, flux_array=0, bs_replicates=0, df_boot=0, tundra_std=0):
        self.gas = gas

    def readin_fluxes_tch4(self, fname, startdate_, enddate_, usecase):
        """
        Returns an array of tundra CH4 fluxes in units of mg CH4-C m-2 d-1 or mmol CO2-eq m-2 d-1
        :param fname: file path and name of file with EC data
        :param startdate_: date to start EC data selection on
        :param enddate_: date to end EC data selection on
        :param usecase: indicate CH4-C to convert fluxes to mg CH4-C m-2 d-1 any other string to convert to mmol CO2-eq m-2 d-1
        :return: dataframe without outliers based on 1.5IQR
        """
        
        fluxnet = pd.read_csv(fname, delimiter=(','))
        fluxnet = fluxnet[['TIMESTAMP_START', 'FCH4_F', 'WD']]
        # tundra units nmolCH4 m-2 s-1

        # make Date column
        # make string version of original timestamp column, call it 'col'
        fluxnet['col'] = fluxnet['TIMESTAMP_START'].astype(str)

        # make the new columns using string indexing
        fluxnet['year'] = fluxnet['col'].str[0:4]
        fluxnet['month'] = fluxnet['col'].str[4:6]
        fluxnet['day'] = fluxnet['col'].str[6:8]
        fluxnet['hour'] = fluxnet['col'].str[8:10]
        fluxnet['minute'] = fluxnet['col'].str[10:12]
        cols_date = ['year', 'month', 'day']
        cols_time = ['hour', 'minute']
        fluxnet['Date'] = fluxnet[cols_date].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
        fluxnet['time'] = fluxnet[cols_time].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)

        # get rid of the extra variable (if you want)
        fluxnet.drop(['TIMESTAMP_START', 'col', 'year', 'month', 'day', 'hour', 'minute'], axis=1, inplace=True)
        
        fluxnet['Date'] =  pd.to_datetime(fluxnet['Date'], format='%Y-%m-%d').dt.date
        fluxnet['time'] =  pd.to_datetime(fluxnet['time'], format='%H:%M').dt.time
        
        # Make a mask for time, only select fluxes during day time
        starttime = pd.to_datetime('08:30').time()
        endtime = pd.to_datetime('17:00').time()
        mask_time = (fluxnet['time'] > starttime) & (fluxnet['time'] <= endtime)
        
        startdate = pd.to_datetime(startdate_).date()
        enddate = pd.to_datetime(enddate_).date()
        mask_16 = (fluxnet['Date'] > startdate) & (fluxnet['Date'] <= enddate)
        t_flux_16 = fluxnet.loc[mask_16]
        t_flux_16 = t_flux_16.loc[mask_time]
        
        # Remove rows where WD = -9999
        t_flux_16 = t_flux_16 [t_flux_16.WD != -9999]
        # Remove rows where FCH4_F = -9999
        t_flux_16 = t_flux_16 [t_flux_16.FCH4_F != -9999]
        
        if usecase == 'CH4-C':
            # Convert flux from units nmol CH4 m-2 s-1 to mg CH4-C m-2 d-1
            # in g/mol
            ch4_atomic_mass = 16.04
            ch4_to_c = 0.75
        
            t_flux_16['ch4_c_mg'] = t_flux_16['FCH4_F'] * 86400 * 10**-6 * ch4_atomic_mass * ch4_to_c

            # Select fluxes in low surface water (< 0.5%) wind areas, from Parmentier et al., 2001
            t_flux_16_dry1 = t_flux_16.drop(t_flux_16[(t_flux_16.WD < 325) | (t_flux_16.WD > 360)].index)
            t_flux_16_dry2 = t_flux_16.drop(t_flux_16[(t_flux_16.WD > 95)].index)
            t_flux_16_mixed = t_flux_16.drop(t_flux_16[(t_flux_16.WD > 250) | (t_flux_16.WD < 190)].index)
        
            # Create new array with fluxes from wind areas
            self.flux_array = np.array(t_flux_16_dry1['ch4_c_mg'])
            a = t_flux_16_dry2['ch4_c_mg'].to_numpy()
            b = t_flux_16_mixed['ch4_c_mg'].to_numpy()
            self.flux_array = np.append(self.flux_array, a)
            self.flux_array = np.append(self.flux_array, b)
            
        else:
            # Convert flux from units nmol CH4 m-2 s-1 to mmol CH4 m-2 d-1
            t_flux_16['ch4_mmol'] = t_flux_16['FCH4_F'] * 86400 * 10**-6 
            
            ch4_sgwp = 45
            ch4_sgcp = 203
            
            t_flux_16['ch4 flux co2-eq'] = (t_flux_16['ch4_mmol']*ch4_sgwp).where(t_flux_16['ch4_mmol']>0)
            t_flux_16['ch4 sgcp'] = (t_flux_16['ch4_mmol']*ch4_sgcp).where(t_flux_16['ch4_mmol']<0)
            t_flux_16['ch4 flux co2-eq'].fillna(t_flux_16['ch4 sgcp'], inplace=True)
            
            # Select fluxes in low surface water (< 0.5%) wind areas, from Parmentier et al., 2001
            t_flux_16_dry1 = t_flux_16.drop(t_flux_16[(t_flux_16.WD < 325) | (t_flux_16.WD > 360)].index)
            t_flux_16_dry2 = t_flux_16.drop(t_flux_16[(t_flux_16.WD > 95)].index)
            t_flux_16_mixed = t_flux_16.drop(t_flux_16[(t_flux_16.WD > 250) | (t_flux_16.WD < 190)].index)
            
            # Create new array with fluxes from wind areas in unit mmol CO2-eq m-2 d-1
            self.flux_array = np.array(t_flux_16_dry1['ch4 flux co2-eq'])
            a = t_flux_16_dry2['ch4 flux co2-eq'].to_numpy()
            b = t_flux_16_mixed['ch4 flux co2-eq'].to_numpy()
            self.flux_array = np.append(self.flux_array, a)
            self.flux_array = np.append(self.flux_array, b)
            
            
            
        
    def readin_fluxes_tn2o(self, fname):
        """ Returns an array of tundra N2O fluxes in units of mmol CO2-eq m-2 d-1 """
        
        n2o_tundra = pd.read_excel(fname, engine='openpyxl', skiprows=1)
        # Convert to mmol m-2 d-1
        n2o_atomic_mass = 44.013
        n2o_tundra['n2o flux'] =  n2o_tundra['flux'] * 10**-6 / n2o_atomic_mass  * 10**3
        
        n2o_sgwp = 270
        n2o_sgcp = 349
        
        n2o_tundra['n2o flux co2-eq'] = (n2o_tundra['n2o flux']*n2o_sgwp).where(n2o_tundra['n2o flux']>0)
        n2o_tundra['n2o flux sgcp'] = (n2o_tundra['n2o flux']*n2o_sgcp).where(n2o_tundra['n2o flux']<0)
        n2o_tundra['n2o flux co2-eq'].fillna(n2o_tundra['n2o flux sgcp'], inplace=True)
        
        self.flux_array = n2o_tundra['n2o flux co2-eq'].to_numpy()
        
    def readin_fluxes_tco2(self, fname, startdate_, enddate_, usecase):
        """
        Returns an array of tundra CO2 fluxes in units of mg CO2-C m-2 d-1 or mmol m-2 d-1
        :param fname: file path and name of file with EC data
        :param startdate_: date to start EC data selection on
        :param enddate_: date to end EC data selection on
        :param usecase: indicate CO2-C to convert fluxes to mg CO2-C m-2 d-1 any other string to convert to mmol  m-2 d-1
        """
        
        flux = pd.read_csv(fname, delimiter=(','), dtype={'TIMESTAMP1': str, 'WD_10': float, 'NEE_filter_ustar': float})
        # tundra units micromol CO2 m-2 s-1

        # make Date column        
        # make the new columns using string indexing
        flux['year'] = flux['TIMESTAMP1'].str[0:4]
        flux['month'] = flux['TIMESTAMP1'].str[5:7]
        flux['day'] = flux['TIMESTAMP1'].str[8:10]
        flux['hour'] = flux['TIMESTAMP1'].str[11:13]
        flux['minute'] = flux['TIMESTAMP1'].str[14:16]
        cols_date = ['year', 'month', 'day']
        cols_time = ['hour', 'minute']
        flux['Date'] = flux[cols_date].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
        flux['time'] = flux[cols_time].apply(lambda row: ':'.join(row.values.astype(str)), axis=1)

        # get rid of the extra variable (if you want)
        flux.drop(['TIMESTAMP1', 'year', 'month', 'day', 'hour', 'minute'], axis=1, inplace=True)
        
        flux['Date'] =  pd.to_datetime(flux['Date'], format='%Y-%m-%d').dt.date
        flux['time'] =  pd.to_datetime(flux['time'], format='%H:%M').dt.time
        
        
        # Make a mask for time, only select fluxes during day time
        starttime = pd.to_datetime('08:30').time()
        endtime = pd.to_datetime('17:00').time()
        mask_time = (flux['time'] > starttime) & (flux['time'] <= endtime)
        
        startdate = pd.to_datetime(startdate_).date()
        enddate = pd.to_datetime(enddate_).date()
        mask_ = (flux['Date'] > startdate) & (flux['Date'] <= enddate)
        t_flux = flux.loc[mask_]
        t_flux = t_flux.loc[mask_time]
        
        # RDrop columns with nans
        t_flux.dropna(subset=['WD_10', 'NEE_filter_ustar'], how='any', inplace=True)
        
        # Convert flux from units micromol CO2 m-2 s-1 to mmol CO2 m-2 d-1
        t_flux['co2_mmol'] = t_flux['NEE_filter_ustar'] * 86400 * 10**-3 
            
        if usecase == 'CO2-C':
            # Convert flux from units micromol CO2 m-2 s-1 to mg CO2-C m-2 d-1
            # in g/mol
            co2_atomic_mass = 44.01
            co2_to_c = 0.27
        
            t_flux['co2_c_mg'] = t_flux['NEE_filter_ustar'] * 86400 * 10**-3 * co2_atomic_mass * co2_to_c

            # Select fluxes in low surface water (< 0.5%) wind areas, from Parmentier et al., 2001
            t_flux_dry1 = t_flux.drop(t_flux[(t_flux.WD_10 < 325) | (t_flux.WD_10 > 360)].index)
            t_flux_dry2 = t_flux.drop(t_flux[(t_flux.WD_10 > 95)].index)
            t_flux_mixed = t_flux.drop(t_flux[(t_flux.WD_10 > 250) | (t_flux.WD_10 < 190)].index)
        
            # Create new array with fluxes from wind areas
            self.flux_array = np.array(t_flux_dry1['co2_c_mg'])
            a = t_flux_dry2['co2_c_mg'].to_numpy()
            b = t_flux_mixed['co2_c_mg'].to_numpy()
            self.flux_array = np.append(self.flux_array, a)
            self.flux_array = np.append(self.flux_array, b)
            
        else:
            # Convert flux from units micromol CO2 m-2 s-1 to mmol CO2 m-2 d-1
            t_flux['co2_mmol'] = t_flux['NEE_filter_ustar'] * 86400 * 10**-3 
            
            # Select fluxes in low surface water (< 0.5%) wind areas, from Parmentier et al., 2001
            t_flux_dry1 = t_flux.drop(t_flux[(t_flux.WD_10 < 325) | (t_flux.WD_10 > 360)].index)
            t_flux_dry2 = t_flux.drop(t_flux[(t_flux.WD_10 > 95)].index)
            t_flux_mixed = t_flux.drop(t_flux[(t_flux.WD_10 > 250) | (t_flux.WD_10 < 190)].index)
            
            # Create new array with fluxes from wind areas in unit mmol CO2-eq m-2 d-1
            self.flux_array = np.array(t_flux_dry1['co2_mmol'])
            a = t_flux_dry2['co2_mmol'].to_numpy()
            b = t_flux_mixed['co2_mmol'].to_numpy()
            self.flux_array = np.append(self.flux_array, a)
            self.flux_array = np.append(self.flux_array, b)
        
        
    def bootstrap_replicate_1d(self, data, func):
        """ Generate bootstrap replicate of 1D data. """
        
        bs_sample = np.random.choice(data, len(data))
        
        return func(bs_sample)
    
    
    def draw_bs_reps(self, data, func, size):
        """ Draw bootstrap replicates """ 
        
        # Initialise array or replicates
        bs_replicates = np.empty(size)
        
        # Generate replicates
        for i in range(size):
            bs_replicates[i] = self.bootstrap_replicate_1d(data, func)
            
        return bs_replicates
    
    
    def bootstrapping(self, boot_n):      
        """
        Generates bootstrapped median samples 

        :param varble: variable to select from dataframe and bootstrap
        :param boot_n: number of samples to generate
        :return: array of boostrapped median samples

        """ 
        
        bs_replicates = self.draw_bs_reps(self.flux_array, np.median, boot_n)
        
        # Make a histogram of the results
        plt.hist(bs_replicates, bins=50, density=True)
        plt.xlabel('median ' + self.gas + ' tundra flux')
        plt.ylabel('PDF')
        
        plt.show()
        
        # Compute the 95% confidence interval: conf_int
        conf_int = np.percentile(bs_replicates, [5, 95])
        tundra_median = np.median(bs_replicates)
        self.tundra_std =  np.std(bs_replicates)
        
        print('The bootstrapped median tundra {} fluxes is {:0.2f} '.format(self.gas, tundra_median)) 
        print('The standard error of tundra {} fluxes is {:0.2f} '.format(self.gas, self.tundra_std)) 
        print('95% confidence interval of tundra {} fluxes is {:0.2f} to {:0.2f}'.format(self.gas, conf_int[0], conf_int[1]))
        
        return bs_replicates        

         


      















