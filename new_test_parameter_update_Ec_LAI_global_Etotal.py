# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 20:33:33 2022

@author: DELL
"""
# for resolution of 360*720 and lat is -90-90 lon is -180 t0 180

import math
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
# import interpolation_lithology as ilith
import pandas as pd
from scipy.interpolate import griddata
from scipy.optimize import fsolve
from scipy import optimize
import sys
import copy
import gc
import multiprocessing as mp
from mpi4py import MPI


# ---------------------------------------#
def convert_lith_sediment_pnas2011033117(lithology_low):
    # Use the define of pnas2011033117
    # https://www.pnas.org/content/pnas/suppl/2020/09/24/2011033117.DCSupplemental/pnas.2011033117.sapp.pdf
    lithology_carbonate = lithology_low[5, :, :] + lithology_low[13, :, :]
    lithology_metamorphic = lithology_low[7, :, :] + lithology_low[14, :, :]
    lithology_sediment = lithology_low[2, :, :] + lithology_low[4, :, :] + lithology_low[0, :, :]
    lithology_mafic = lithology_low[1, :, :] + lithology_low[3, :, :]
    lithology_felsic = lithology_low[6, :, :] + lithology_low[8, :, :] + lithology_low[11, :, :]
    lithology_ice = lithology_low[10, :, :] + lithology_low[15, :, :]
    lithology_intermediate = lithology_low[9, :, :] + lithology_low[12, :, :]
    return lithology_carbonate, lithology_metamorphic, lithology_sediment, lithology_mafic, lithology_felsic, lithology_ice, lithology_intermediate


end_yrs = 200
conv_yr_s = 31557600.0
const_rEarth = 6.37e6
CO2 = 340.12
Te = 20 * 1e15  # g/yr
rou = 2500 * 1e3  # g/m3
m = 0.5
n = 1

directory_name = '/work4/zhy/We-carb/erosion_before/'
directory_name_file = '/work4/zhy/We-carb/erosion_before/file_area_new/'
directory_name_csv = '/work4/zhy/We-carb/erosion_before/file_csv_new/'
filename_glim = 'glim_wgs84_0point5deg.txt'
filename_area = 'land_area_360_720_PD_match_clim_slope_lith.nc'
filename_climate1 = 'modern_runoff_v2_interp.nc'
filename_slope = '0Ma_slope.nc'

f_area = Dataset(directory_name + filename_area)
area_raw = f_area['area'][:]
f_lnd = Dataset(directory_name + filename_climate1)
R = f_lnd['rnf'][:]  # m/yr
R = np.maximum(R, 0)
f_slope = Dataset(directory_name + filename_slope)
s = f_slope['slope'][:]

filename_veg = 'Vegetation_from_surfdata_new.nc'
f_v = Dataset(directory_name + filename_veg)
Veg = f_v['sum_LAI'][:]
Veg_index = np.exp(-np.minimum(2, Veg))

E_par2 = pow(R, m) * pow(s, n) * Veg_index
E_par2_0 = np.nansum(E_par2 * area_raw)
ke2 = Te / rou / E_par2_0
E = ke2 * E_par2


E.dump(directory_name_file+'Erosion_calculated_v2_LAI_global_Etotal')

lithology_names_char = open(directory_name + filename_glim, 'r')
lithology_names = []
for line in lithology_names_char.readlines():
    list1 = line.strip('\n').split(sep=' ')  ##数字之间的空格消除
    list1 = [x.strip() for x in list1 if x.strip() != '']
    if len(list1) != 0:
        lithology_names.append(list1)
lithology_names_char.close()
lithology_names = np.array(lithology_names)
lithology_names = lithology_names.astype(float)
nm = lithology_names.shape
lithology_names_1 = np.zeros(nm)
for i in range(nm[0]):
    for j in range(nm[1]):
        lithology_names_1[(nm[0] - 1) - i, j] = lithology_names[i, j]
lithology_names_1 = lithology_names_1.astype(float)
lithology_names_1[lithology_names_1 == -9999.] = np.nan
# --------------------------------------------#
lithology_low = np.zeros((16, 360, 720))
for i in range(nm[0]):
    for j in range(nm[1]):
        counts = {}
        p = lithology_names_1[i, j]
        counts[p] = counts.get(p, 0) + 1
        for l in range(16):
            lithology_low[l, i, j] = counts.get(l + 1, 0)
carbonate, metamorphic, sediment, mafic, felsic, ice, intermediate = convert_lith_sediment_pnas2011033117(
    lithology_low)
dum_lithology = np.zeros((7, carbonate.shape[0], carbonate.shape[1]))
dum_lithology[0, :, :] = carbonate
dum_lithology[1, :, :] = metamorphic
dum_lithology[2, :, :] = sediment
dum_lithology[3, :, :] = mafic
dum_lithology[4, :, :] = felsic
dum_lithology[5, :, :] = ice
dum_lithology[6, :, :] = intermediate
dum_lithology.dump(directory_name_file+'filename_lithology_v2_global_new_Etotal')



def main(kd, kw, krp, sigma, lithology_test, countnn, test, T, R, E_u, area_raw, G1999_basin_new, f_G1999, obs_new, ama_name1, f_amazon, obs_ama1):
    # ---------------------------------------Parameter--------------------------------------#
    yr = 0
    num = '01'
    Ea = 42000  # KJ/mol
    R_s = 8.314472
    T0 = 286.0  # K
    d0 = 2.73  # m
    P0 = krp * R * np.exp(-Ea / R_s * (1.0 / T - 1.0 / T0))
    h = np.maximum(0, d0 * np.log(P0 / E_u))
    K_dis = kd * (1 - np.exp(-kw * R)) * np.exp(-Ea / R_s * (1.0 / T - 1.0 / T0))
    dum_total_calcium_flux = lithology_test * E_u * (
            1 - np.exp(-K_dis / (sigma + 1) * pow(h / E_u, sigma + 1)))
    weather_fCaSiO3 = np.nansum(dum_total_calcium_flux * area_raw)  # mol/yr
    # np.save(filename_ini, dum_total_calcium_flux)
    # f_w = Dataset(directory_name_file + filename_ini, 'w', format='NETCDF4')
    # f_w.description = "kd=" + str(kd) + " kw=" + str(kw) + " krp=" + str(krp) + " sigma=" + str(
    #     sigma) + " and sum weathering flux is " + str(weather_fCaSiO3) + "mol/yr"
    # f_w.createDimension('time', 1)
    # f_w.createDimension('lat', len(lat))
    # f_w.createDimension('lon', len(lon))
    # f_w.createVariable('time', np.int, ('time'))
    # f_w.createVariable('lat', np.float32, ('lat'))
    # f_w.createVariable('lon', np.float32, ('lon'))
    # f_w.variables['time'][:] = yr
    # # can add a time=yr to make this much better
    # f_w.variables['lon'][:] = lon
    # f_w.variables['lat'][:] = lat
    # f_w.createVariable('R', np.float32, ('time', 'lat', 'lon'))
    # f_w.createVariable('T', np.float32, ('time', 'lat', 'lon'))
    # f_w.createVariable('weather_fCaSiO3_2D', np.float32, ('time', 'lon', 'lat'))
    # f_w.variables['R'][:] = R
    # f_w.variables['T'][:] = T
    # f_w.variables['weather_fCaSiO3_2D'][:] = dum_total_calcium_flux.T
    # f_w.close()

    # --------------------------------Second part-------------------------------------------
    # ----------------------------------------Third part-------------------------------
    flux_obs = np.zeros((len(G1999_basin_new)))
    for i in range(len(G1999_basin_new)):
        ng = G1999_basin_new[i]
        flux_obs[i] = np.nansum(dum_total_calcium_flux * f_G1999[ng][:])
    R_2_global = 1 - np.sum(pow((flux_obs - obs_new), 2)) / np.sum(pow((obs_new - (np.mean(obs_new))), 2))
    R_2_global_log = 1 - np.sum(pow((np.log10(flux_obs) - np.log10(obs_new)), 2)) / np.sum(
        pow((np.log10(obs_new) - np.mean(np.log10(obs_new))), 2))
    model_ama = np.zeros(len(ama_name1))
    for i in range(len(ama_name1)):
        nn = ama_name1[i]
        model_ama[i] = np.nansum(dum_total_calcium_flux * f_amazon[nn][:])
    # flux_obs1 = np.delete(flux_obs, 21)
    # obs1 = np.delete(obs, 21)
    model_flux = np.append(model_ama, flux_obs)
    obs_flux = np.append(obs_ama1 * 1e9, obs_new)
    R_2_total_log = 1 - np.sum(pow((np.log10(model_flux) - np.log10(obs_flux)), 2)) / np.sum(
        pow((np.log10(obs_flux) - np.mean(np.log10(obs_flux))), 2))
    R_2_total = 1 - np.sum(pow((model_flux - obs_flux), 2)) / np.sum(pow((obs_flux - (np.mean(obs_flux))), 2))
    return R_2_global, R_2_global_log, R_2_total_log, R_2_total, weather_fCaSiO3, model_flux


metamorphic_lith = np.array([1500, 2000, 2500, 3000, 3500, 4000])
sediment_lith = np.array([500, 1000, 1500, 2000, 2500, 3000])
xi_lithology_t = np.array([[0, i, j, 10317, 1521, 0, 4759] for i in metamorphic_lith for j in sediment_lith])
krp_a = np.array([1.2 * 1e-3, 2 * 1e-3, 3 * 1e-3, 4 * 1e-3, 5 * 1e-3, 6 * 1e-3])
krp_b = np.array([7 * 1e-3, 8 * 1e-3, 9 * 1e-3, 1 * 1e-2, 1.5 * 1e-2, 5 * 1e-2])
krp_t = np.array([krp_a, krp_b])
result = [[i, j] for i in xi_lithology_t for j in krp_t]
a = np.linspace(1, 72, 72).astype("int")
result = [(x, y) for x, y in zip(list(a), result)]
# result_new = [result[i:i+24] for i in range(0, 72, 24)]


def calculate_weathering(krp_xi_lithology):
    directory_name = '/work4/zhy/We-carb/erosion_before/'
    directory_name_csv = '/work4/zhy/We-carb/erosion_before/file_csv_new/'
    directory_name_file = '/work4/zhy/We-carb/erosion_before/file_area_new/'
    filename_climate1 = 'modern_runoff_v2_interp.nc'
    f_lnd = Dataset(directory_name + filename_climate1)
    filename_climate2 = 'modern_temp_interp.nc'
    f_cam = Dataset(directory_name + filename_climate2)
    T = f_cam['tmp'][:] + 273.15  # K
    R = f_lnd['rnf'][:]  # m/yr
    R = np.maximum(R, 0)
    E_u = np.load(directory_name_file+'Erosion_calculated_v2_LAI_global_Etotal', allow_pickle=True)
    # filename_veg = 'Vegetation_from_surfdata.nc'
    # f_v = Dataset(directory_name + filename_veg)
    # Veg = f_v['tf_LAI_'][:]
    # Veg_index = np.exp(-np.minimum(2, Veg))
    # E_u = E_u * Veg_index
    filename_area = 'land_area_360_720_PD_match_clim_slope_lith.nc'
    f_area = Dataset(directory_name + filename_area)
    area_raw = f_area['area'][:]
    filename_G1999 = 'G1999_basin_update.nc'
    f_G1999 = Dataset(directory_name + filename_G1999)
    filename_G1999_basin_new = 'G1999_basin_park_new.csv'
    f_G1999_basin_new = pd.read_csv(directory_name + filename_G1999_basin_new)
    G1999_basin_new = f_G1999_basin_new['Basin_name'].values
    G1999_basin_new_g = f_G1999_basin_new['Basins'].values
    # area_new = f_G1999_basin_new['Area_1e6km2'].values * 1e6  # km2
    # continent_new = f_G1999_basin_new['Continent'].values
    # region_new = f_G1999_basin_new['Region'].values
    obs_new = f_G1999_basin_new['Ca_Mg'].values * 1e9
    filename_amazon = 'amazon_basins_360_720_global.nc'
    f_amazon = Dataset(directory_name + filename_amazon)
    filename_amazon_obs = 'amazon_basin_results.csv'
    f_obs_ama = pd.read_csv(directory_name + filename_amazon_obs)
    ama_name1 = f_obs_ama['basin'].values
    obs_ama1 = f_obs_ama['CaMg_sil'].values
    name = np.append(ama_name1, G1999_basin_new_g)
    kd_a = np.array(
        [5 * 1e-6, 1 * 1e-5, 2 * 1e-5, 5 * 1e-5, 1 * 1e-4, 2 * 1e-4, 5 * 1e-4, 1e-3, 2 * 1e-3, 5 * 1e-3, 1 * 1e-2])
    kw_a = np.array([1 * 1e-3, 2 * 1e-3, 5 * 1e-3, 1 * 1e-2, 2 * 1e-2, 5 * 1e-2, 1 * 1e-1, 2 * 1e-1, 5 * 1e-1, 1])
    sigma_a = np.array([-0.5, -0.4, -0.2, -0.1, 0, 0.1, 0.3])
    dum_lithology = np.load(directory_name_file+'filename_lithology', allow_pickle=True)
    print(krp_xi_lithology)
    xi_lithology = krp_xi_lithology[1][0]
    # print(xi_lithology.shape)
    krp_p = krp_xi_lithology[1][1]
    dum_lithology_sum = np.zeros((360, 720))
    for k in range(7):
        dum_lithology_sum = dum_lithology_sum + xi_lithology[int(k)] * dum_lithology[int(k), :, :]
    countnn = 0
    R2_global = []
    R2_log_global = []
    R2_amazon_global = []
    R2_log_amazon_global = []
    kd_n = []
    kw_n = []
    krp_n = []
    sigma_n = []
    flux_global = []
    metamorphic_n = []
    sediment_n = []
    model_flux_river = []
    for p1 in range(len(kd_a)):
        for p2 in range(len(kw_a)):
            for p3 in range(len(krp_p)):
                for p4 in range(len(sigma_a)):
                    countnn = countnn + 1
                    R_2_global, R_2_global_log, R_2_total_log, R_2_total, weather_fCaSiO3, model_flux = main(kd_a[p1], kw_a[p2], krp_p[p3], sigma_a[p4], dum_lithology_sum, countnn, krp_xi_lithology[0], T, R, E_u, area_raw, G1999_basin_new, f_G1999, obs_new, ama_name1, f_amazon, obs_ama1)
                    R2_global.append(R_2_global)
                    R2_log_global.append(R_2_global_log)
                    R2_amazon_global.append(R_2_total)
                    R2_log_amazon_global.append(R_2_total_log)
                    kd_n.append(kd_a[p1])
                    kw_n.append(kw_a[p2])
                    krp_n.append(krp_p[p3])
                    sigma_n.append(sigma_a[p4])
                    metamorphic_n.append(krp_xi_lithology[1][0][1])
                    sediment_n.append(krp_xi_lithology[1][0][2])
                    flux_global.append(weather_fCaSiO3)
                    model_flux_river.append(model_flux)
    Data_mon_mean = {'kd': np.array(kd_n), 'kw': np.array(kw_n), 'krp': np.array(krp_n), 'sigma': np.array(sigma_n),
                     'metamorphic': np.array(metamorphic_n), 'sediment': np.array(sediment_n),
                     'R2_global': np.array(R2_global), 'R2_log_global': np.array(R2_log_global), 'R2_amazon_global':
                         np.array(R2_amazon_global), 'R2_log_amazon_global': np.array(R2_log_amazon_global),
                     'flux_global': np.array(flux_global)}
    river_flux = pd.DataFrame(np.array(model_flux_river))
    river_flux.columns = name
    DataFrame_mon_mean = pd.DataFrame(Data_mon_mean)
    DataFrame_mon_mean = pd.concat([DataFrame_mon_mean, river_flux], axis=1)
    DataFrame_mon_mean.to_csv(directory_name_csv + "obs_mon_mean_parameter_area_test_full_v2_LAI_global_Etotal_"+str(krp_xi_lithology[0])+".csv",
                              index=False, sep=',')
    # return R2_global, R2_amazon_global, R2_log_global, R2_log_amazon_global, flux_global


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
# print(result)
print(rank)
print("result is :")
print(result[rank])
with mp.Pool(processes=72) as pool:
    # 在每个进程中并行计算子列表中的参数，并将结果保存到results列表中
    results = pool.map(calculate_weathering, (result[rank],))

# all_results = comm.gather(results, root=0)
if rank == 0:
    print("Running is over")


