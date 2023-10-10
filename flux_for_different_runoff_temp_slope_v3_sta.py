import pdb

from netCDF4 import Dataset
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import copy
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import random
import pdb
from mpl_toolkits.basemap import Basemap
from scipy.optimize import fsolve
# from sympy import *
from scipy import integrate
import mpmath


# x = symbols('x')


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


conv_yr_s = 31557600.0
directory_name = '/work4/zhy/We-carb/different_runoff/'
directory_fig_name = '/work4/zhy/We-carb/different_runoff/fig/'
index1 = np.array(['Southern American', 'Southern Asia', 'Northern American', 'Northern Asia', 'Australia', 'Arc',
                   'Europe', 'African'])
index_color = np.array(['crimson', 'dimgrey', 'y', 'c', 'cornflowerblue', 'orchid', 'palegreen', 'orange'])
mapping_dict = dict(zip(index1, index_color))


def plot_bar_comp_sort(model, obs, name, color, p, unit, ind_u):
    # The model and obs don't need to be sorted, we will sort in this function
    sort = np.argsort(obs)
    obs_sort = np.zeros(p)
    model_sort = np.zeros(p)
    name_sort = np.zeros(p)
    name_sort = name_sort.astype(str)
    color_sort = np.zeros(p)
    color_sort = name_sort.astype(str)
    for i in range(p):
        obs_sort[i] = obs[sort[i]]
        model_sort[i] = model[sort[i]]
        name_sort[i] = name[sort[i]]
        color_sort[i] = color[sort[i]]
    index = np.linspace(0, p - 1, p)
    figure = plt.figure()
    figure.set_size_inches(16, 9)
    a = plt.bar(index, obs_sort[::-1], width=0.4, label='Obs ' + unit, color='white', edgecolor=color_sort[::-1],
                tick_label=name_sort[::-1], zorder=1)
    b = plt.scatter(index, model_sort[::-1], s=60, marker='*', color=color_sort[::-1], zorder=2)
    plt.xticks(rotation=70)
    plt.subplots_adjust(bottom=0.180)
    plt.legend()
    # plt.title(title)
    plt.savefig(directory_fig_name + ind_u, dpi=1200, bbox_inches="tight")
    # plt.show()


def plot_bar_comp_sort2(model, obs, obs2, name, color, p, unit, ind_u):
    # The model and obs don't need to be sorted, we will sort in this function
    sort = np.argsort(obs)
    obs_sort = np.zeros(p)
    obs2_sort = np.zeros(p)
    model_sort = np.zeros(p)
    name_sort = np.zeros(p)
    color_sort = np.zeros(p)
    name_sort = name_sort.astype(str)
    color_sort = color_sort.astype(str)
    for i in range(p):
        obs_sort[i] = obs[sort[i]]
        model_sort[i] = model[sort[i]]
        obs2_sort[i] = obs2[sort[i]]
        name_sort[i] = name[sort[i]]
        color_sort[i] = color[sort[i]]
    print(color_sort)
    index = np.linspace(0, p - 1, p)
    figure = plt.figure()
    figure.set_size_inches(16, 9)
    a = plt.bar(index, obs2_sort[::-1], width=0.2, label='observation from Moon ' + unit, color='white',
                edgecolor=color_sort[::-1],
                hatch='//')
    index = index + 0.2
    c = plt.bar(index, obs_sort[::-1], width=0.2, label='observation from Gaillardet ' + unit, color='white',
                edgecolor=color_sort[::-1], tick_label=name_sort[::-1], zorder=1)
    # for i in range(len(index)):
    #     index[i] = index[i] + 0.2
    b = plt.scatter(index, model_sort[::-1], s=60, marker='*', color=color_sort[::-1], zorder=2)
    plt.xticks(rotation=70)
    plt.subplots_adjust(bottom=0.180)
    plt.legend(loc='upper right')
    # plt.title(title)
    plt.savefig(directory_fig_name + ind_u, dpi=1200, bbox_inches="tight")
    # plt.show()


def plot_residuals(obs, model, obs_per, model_per, ind_nt, ind_tropical, ind_u):
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axes[0, 0].plot(np.log10(obs[ind_tropical]), (model - obs)[ind_tropical], 'o')
    axes[0, 0].axhline(y=0, color='r', linestyle='-')
    axes[0, 0].set_xlim(8, 11.2)
    # axes[0, 0].set_ylim(-1e11, 3e11)
    axes[0, 0].set_ylim(-1.3e11, 1.3e11)

    axes[0, 1].plot(np.log10(obs[ind_nt]), (model - obs)[ind_nt], 'o')
    axes[0, 1].axhline(y=0, color='r', linestyle='-')
    axes[0, 1].set_xlim(8, 11.2)
    # axes[0, 1].set_ylim(-1e11, 3e11)
    axes[0, 1].set_ylim(-1.3e11, 1.3e11)

    axes[0, 2].plot(np.log10(obs), (model - obs), 'o')
    axes[0, 2].axhline(y=0, color='r', linestyle='-')
    # axes[0, 2].set_xlim(0, 1.7e11)
    axes[0, 2].set_xlim(8, 11.2)
    # axes[0, 2].set_ylim(-1e11, 3e11)
    axes[0, 2].set_ylim(-1.3e11, 1.3e11)

    axes[1, 0].plot(np.log10(obs[ind_tropical]), (np.log10(model) - np.log10(obs))[ind_tropical], 'o')
    axes[1, 0].axhline(y=0, color='r', linestyle='-')
    axes[1, 0].set_xlim(8, 11.2)
    axes[1, 0].set_ylim(-1.2, 1.2)

    axes[1, 1].plot(np.log10(obs[ind_nt]), (np.log10(model) - np.log10(obs))[ind_nt], 'o')
    axes[1, 1].axhline(y=0, color='r', linestyle='-')
    axes[1, 1].set_xlim(8, 11.2)
    axes[1, 1].set_ylim(-1.2, 1.2)

    axes[1, 2].plot(np.log10(obs), (np.log10(model) - np.log10(obs)), 'o')
    axes[1, 2].axhline(y=0, color='r', linestyle='-')
    axes[1, 2].set_xlim(8, 11.2)
    axes[1, 2].set_ylim(-1.2, 1.2)

    # axes[1, 0].plot(obs_per[ind_tropical], (model_per - obs_per)[ind_tropical], 'o')
    # axes[1, 0].axhline(y=0, color='r', linestyle='-')
    # axes[1, 0].set_xlim(0, 450000)
    # axes[1, 0].set_ylim(-200000, 250000)
    #
    # axes[1, 1].plot(obs_per[ind_nt], (model_per - obs_per)[ind_nt], 'o')
    # axes[1, 1].axhline(y=0, color='r', linestyle='-')
    # axes[1, 1].set_xlim(0, 450000)
    # axes[1, 1].set_ylim(-200000, 250000)
    #
    # axes[1, 2].plot(obs_per, (model_per - obs_per), 'o')
    # axes[1, 2].axhline(y=0, color='r', linestyle='-')
    # axes[1, 2].set_xlim(0, 450000)
    # axes[1, 2].set_ylim(-200000, 250000)
    #
    # fig.suptitle(name)
    # plt.show()
    plt.savefig(directory_fig_name + ind_u, dpi=1200, bbox_inches="tight")


def plot_residuals_percent(obs, model, obs_per, model_per, ind_nt, ind_tropical, ind_u):
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    axes[0, 0].plot(obs[ind_tropical], ((model - obs) / obs)[ind_tropical], 'o')
    axes[0, 0].axhline(y=0, color='r', linestyle='-')
    axes[0, 0].set_xlim(0, 1.7e11)
    axes[0, 0].set_ylim(-2, 12)

    axes[0, 1].plot(obs[ind_nt], ((model - obs) / obs)[ind_nt], 'o')
    axes[0, 1].axhline(y=0, color='r', linestyle='-')
    axes[0, 1].set_xlim(0, 1.7e11)
    axes[0, 1].set_ylim(-2, 12)

    axes[0, 2].plot(obs, ((model - obs) / obs), 'o')
    axes[0, 2].axhline(y=0, color='r', linestyle='-')
    axes[0, 2].set_xlim(0, 1.7e11)
    axes[0, 2].set_ylim(-2, 12)

    axes[1, 0].plot(obs_per[ind_tropical], ((model_per - obs_per) / obs_per)[ind_tropical], 'o')
    axes[1, 0].axhline(y=0, color='r', linestyle='-')
    axes[1, 0].set_xlim(0, 450000)
    axes[1, 0].set_ylim(-2, 12)

    axes[1, 1].plot(obs_per[ind_nt], ((model_per - obs_per) / obs_per)[ind_nt], 'o')
    axes[1, 1].axhline(y=0, color='r', linestyle='-')
    axes[1, 1].set_xlim(0, 450000)
    axes[1, 1].set_ylim(-2, 12)

    axes[1, 2].plot(obs_per, ((model_per - obs_per) / obs_per), 'o')
    axes[1, 2].axhline(y=0, color='r', linestyle='-')
    axes[1, 2].set_xlim(0, 450000)
    axes[1, 2].set_ylim(-2, 12)
    # plt.savefig(name, dpi=1200, bbox_inches="tight"
    # fig.suptitle(name)
    # plt.show()
    plt.savefig(directory_fig_name + ind_u, dpi=1200, bbox_inches="tight")


def t_main(T, R, s, Be, sta, B, dum_lithology, fig_name):
    kd = 5 * 1e-4
    kw = 1
    Ea = 42000  # KJ/mol
    R_s = 8.314472
    T0 = 286.0  # K
    Te = 20 * 1e15  # g/yr
    rou = 2500 * 1e3  # g/m3
    m = 0.5
    n = 1
    krp = 0.01
    d0 = 2.73  # m
    sigma = -0.4
    xi_lithology = np.array([0, 2500, 2000, 10317, 1521, 0, 4759])
    conv_yr_s = 31557600.0
    const_rEarth = 6.37e6
    CO2 = 340.12

    if sta == 1:
        filename_G1999 = 'G1999_basin_update.nc'
        f_G1999 = Dataset(directory_name + filename_G1999)
        filename_G1999_basin_new = 'G1999_basin_park_new.csv'
        f_G1999_basin_new = pd.read_csv(directory_name + filename_G1999_basin_new)
        G1999_basin_new = f_G1999_basin_new['Basin_name'].values
        G1999_basin_new_g = f_G1999_basin_new['Basins'].values
        area_new = f_G1999_basin_new['Area_1e6km2'].values * 1e6  # km2
        continent_new = f_G1999_basin_new['Continent'].map(mapping_dict).values
        region_new = f_G1999_basin_new['Region'].values
        obs_new = f_G1999_basin_new['Ca_Mg'].values * 1e9
    else:
        filename_G1999 = 'G1999_basins_360_720.nc'
        f_G1999 = Dataset(directory_name + filename_G1999)
        filename_G1999_basin_new = 'G1999_basin_park_from_park.csv'
        f_G1999_basin_new = pd.read_csv(directory_name + filename_G1999_basin_new)
        G1999_basin_new = f_G1999_basin_new['Basin_name'].values
        G1999_basin_new_g = f_G1999_basin_new['Basins'].values
        area_new = f_G1999_basin_new['Area_1e6km2'].values * 1e6  # km2
        continent_new = f_G1999_basin_new['Continent'].map(mapping_dict).values
        region_new = f_G1999_basin_new['Region'].values
        obs_new = f_G1999_basin_new['Ca_Mg'].values * 1e9
    filename_amazon = 'amazon_basins_360_720_global.nc'
    f_amazon = Dataset(directory_name + filename_amazon)
    filename_amazon_obs = 'amazon_basin_results.csv'
    f_obs_ama = pd.read_csv(directory_name + filename_amazon_obs)
    ama_name1 = f_obs_ama['basin'].values
    obs_ama1 = f_obs_ama['CaMg_sil'].values * 1e9
    area_ama_new = f_obs_ama['area'].values * 1e6  # km2
    continent_new_ama = pd.Series(['Southern American'] * len(ama_name1)).map(mapping_dict).values
    region_new_ama = ['Low'] * len(ama_name1)
    obs_new_per_value = obs_new / area_new
    obs_ama1_per_value = obs_ama1 / area_ama_new

    name = np.append(G1999_basin_new_g, ama_name1)
    obs_total = np.append(obs_new, obs_ama1)
    obs_per_value_total = np.append(obs_new_per_value, obs_ama1_per_value)
    region_total = np.append(region_new, region_new_ama)
    area_total = np.append(area_new, area_ama_new)

    Moon_obs_filename = 'Moon_table2.csv'
    Moon_obs = pd.read_csv(directory_name + Moon_obs_filename)['Silicate'].values * 1e9  # mol/yr
    Moon_obs_name = pd.read_csv(directory_name + Moon_obs_filename)['Name'].values
    Moon_obs_area = pd.read_csv(directory_name + Moon_obs_filename)['Area (106km2)'].values * 1e6  # km2
    Moon_obs_u = np.zeros((len(G1999_basin_new_g)))
    Moon_area_u = np.zeros((len(G1999_basin_new_g)))
    for i in range(len(G1999_basin_new_g)):
        ind2 = np.where(Moon_obs_name == G1999_basin_new_g[i])
        if len(ind2[0]) != 0:
            Moon_obs_u[i] = Moon_obs[ind2]
            Moon_area_u[i] = Moon_obs_area[ind2]
    Moon_obs_per_value = Moon_obs_u / Moon_area_u
    filename_obs_climate = 'climate_data_new2.csv'
    obs_name = pd.read_csv(directory_name + filename_obs_climate)['River'].values  # t/km2/yr
    obs_E = pd.read_csv(directory_name + filename_obs_climate)['Erosion'].values  # t/km2/yr
    filename_area = 'land_area_360_720_PD_match_clim_slope_lith.nc'
    f_area = Dataset(directory_name + filename_area)
    area_raw = f_area['area'][:]
    if Be == 1:
        obs_E[np.where(obs_name == 'Mackenzie')] = 222
        obs_E[np.where(obs_name == 'Columbia')] = 224
        obs_E[np.where(obs_name == 'Mississippi')] = 91
        obs_E[np.where(obs_name == 'Orinoco')] = 60
        obs_E[np.where(obs_name == 'Amazon')] = 81
        obs_E[np.where(obs_name == 'Rhine')] = 267
        obs_E[np.where(obs_name == 'Rhone')] = 183
        obs_E[np.where(obs_name == 'Po')] = 564
        obs_E[np.where(obs_name == 'Ganges')] = 491
        obs_E[np.where(obs_name == 'Irrawady')] = 379
        obs_E[np.where(obs_name == 'Mekong')] = 88
        obs_E[np.where(obs_name == 'Hong He')] = 144
        obs_E[np.where(obs_name == 'Xijiang')] = 68
        obs_E[np.where(obs_name == 'Changjiang')] = 139
        obs_E[np.where(obs_name == 'Huanghe')] = 222
        obs_E[np.where(obs_name == 'Zambese')] = 10.8
        obs_E[np.where(obs_name == 'Nile')] = 20
        obs_E[np.where(obs_name == 'Limpopo')] = 8

    def routing_sum(var):
        p = G1999_basin_new.shape[0]
        E_model_g = np.zeros(p)
        for i in range(p):
            ng = G1999_basin_new[i]
            E_model_g[i] = np.nansum(var * f_G1999[ng][:] / 1e6)
        return E_model_g

    def routing_average(var):
        p = G1999_basin_new.shape[0]
        T_model_g = np.zeros(p)
        area_model_g = routing_sum(1)
        for i in range(p):
            ng = G1999_basin_new[i]
            T_model_g[i] = np.nansum(var * f_G1999[ng][:] / 1e6) / area_model_g[i]
        return T_model_g

    # filename_E = 'E_R_v1_s_v2.nc'
    # filename_E_u = 'E_update_1_global.nc'
    # filename_E_Be_u = 'E_update_1_global_Be.nc'
    # f_E = Dataset(directory_name + filename_E)
    # E = f_E['E'][:]  # m/yr
    # f_E_u = Dataset(directory_name + filename_E_u)
    # E_u1 = f_E_u['E'][:]
    # E_Be_u = Dataset(directory_name + filename_E_Be_u)['E'][:]

    E_par2 = pow(R, m) * pow(s, n)  # m/yr
    E_par2_0 = np.sum(E_par2 * area_raw)
    ke2 = Te / rou / E_par2_0
    E2_c = ke2 * E_par2
    print(ke2)
    E2 = np.maximum(0, E2_c)  # m/yr t/km2/yr
    print(np.nansum(E2_c * area_raw) * rou)
    if B == 1:
        # pdb.set_trace()
        E_rout = routing_average(E2) * rou
        E_u = copy.deepcopy(E2)
        for j, name in enumerate(G1999_basin_new_g):
            ind_name = np.where(obs_name == name)[0]
            if not np.isnan(obs_E[ind_name]):
                ratio_1 = obs_E[ind_name] / E_rout[j]
                print(ratio_1)
                E_u[np.where(f_G1999[G1999_basin_new[j]][:] != 0)] = E_u[np.where(
                    f_G1999[G1999_basin_new[j]][:] != 0)] * ratio_1
                # pdb.set_trace()
            print(np.nansum(E_u * area_raw) * rou)
    else:
        E_u = copy.deepcopy(E2)
    # pdb.set_trace()
    print(np.nansum(E_u * area_raw) * rou)
    def mac_model(R, T):
        # ------------------------------------------------------------------------------
        # Some parameters
        # Lithology decides
        CO2 = 340.12  # ppm
        # --------------------------
        n = 0.316  # thermodynamic pCO2 dependence [-]
        Λ = 1.4e-3  # thermodynamic coefficient for Ceq [-]
        ρ = 12728.0  # mineral mass to fluid volume ratio [kg m⁻³]
        k0 = 8.7e-6  # reference rate constant [mol m⁻² yr⁻¹]
        m = 0.27  # mineral molar mass [kg/mol]
        X = 0.36  # reactive mineral conc. in fresh rock [-]
        # Soil thickness decides, slope
        # Thinking the soil thickness maybe decided by the eroison part
        # --------------------------
        L = 1.0  # flow path length [m]
        # Soil thickness and bedrock lithology and physical weathering, erosion, soil age(time)
        # --------------------------
        ϕ = 0.1  # porosity [-]
        A = 1e2  # specific surface area (not weathering surface area) [m²kg⁻¹]
        ts = 1e5  # soil age [yr]
        # --------------------------
        μ = np.exp(2)  # scaling constant [-]
        β = 0.2  # pCO2 scaling [-]

        # -----------------------------------------------------------------------------
        # Calculating
        r = R / conv_yr_s  # from m/yr into m/s
        pCO2 = CO2 / 1e6  # carbon dioxide concentration [units must be same as pCO2₀]
        pCO20 = 285e-6
        Te = 11.1  # scaling of temperature dependence [K]
        T0 = 288.15  # reference temperature [K]

        α = L * ϕ * ρ * A * X * μ
        Ceq = 1e3 * Λ * pow(pCO2, n)  # mol/L
        a = np.exp((T - T0) / Te)
        b = pow(pCO2 / pCO20, β)
        d = 1 / (k0 * a * b) + m * A * ts + α / (r * conv_yr_s * Ceq)
        Weath = (α / d)
        return Weath

    P0 = krp * R * np.exp(-Ea / R_s * (1.0 / T - 1.0 / T0))
    K_dis = kd * (1 - np.exp(-kw * R)) * np.exp(-Ea / R_s * (1.0 / T - 1.0 / T0))
    P0_new = np.zeros((360, 720))
    P0_new[P0.mask == False] = P0[P0.mask == False]
    K_dis_new = np.zeros((360, 720))
    K_dis_new[K_dis.mask == False] = K_dis[K_dis.mask == False]
    h = np.maximum(0, d0 * np.log(P0_new / E_u))
    h_new = np.zeros((360, 720))
    h_new[h.mask == False] = h[h.mask == False]
    E_new = np.zeros((360, 720))
    E_new[E_u.mask == False] = E_u[E_u.mask == False]

    Weath = mac_model(R, T)

    dum_lithology_sum = np.zeros((n_j, n_i))
    for k in range(7):
        dum_lithology_sum = dum_lithology_sum + xi_lithology[k] * dum_lithology[k, :, :]
    print(np.nansum(E_new * area_raw) * rou)
    dum_total_calcium_flux = dum_lithology_sum * E_new * (
            1 - np.exp(-K_dis_new / (sigma + 1) * pow(h_new / E_new, sigma + 1)))
    dum_total_calcium_flux[np.where((E_new == 0) & (R.mask == False))] = 0
    print(np.nansum(dum_total_calcium_flux * area_raw))
    weather_fCaSiO3 = np.nansum(dum_total_calcium_flux * area_raw)
    # pdb.set_trace()
    flux_obs = np.zeros((len(G1999_basin_new)))
    for i in range(len(G1999_basin_new)):
        ng = G1999_basin_new[i]
        flux_obs[i] = np.nansum(dum_total_calcium_flux * f_G1999[ng][:])
    R_2_global = 1 - np.sum(pow((flux_obs - obs_new), 2)) / np.sum(pow((obs_new - (np.mean(obs_new))), 2))
    R_2_global_log = 1 - np.sum(pow((np.log10(flux_obs) - np.log10(obs_new)), 2)) / np.sum(
        pow((np.log10(obs_new) - np.mean(np.log10(obs_new))), 2))
    # np.save(fig_name + '.npy', arr=flux_obs)
    # pdb.set_trace()
    model_ama = np.zeros(len(ama_name1))
    for i in range(len(ama_name1)):
        nn = ama_name1[i]
        model_ama[i] = np.nansum(dum_total_calcium_flux * f_amazon[nn][:])
    model_flux = np.append(flux_obs, model_ama)
    obs_flux = np.append(obs_new, obs_ama1)
    R_2_total_log = 1 - np.sum(pow((np.log10(model_flux) - np.log10(obs_flux)), 2)) / np.sum(
        pow((np.log10(obs_flux) - np.mean(np.log10(obs_flux))), 2))
    R_2_total = 1 - np.sum(pow((model_flux - obs_flux), 2)) / np.sum(pow((obs_flux - (np.mean(obs_flux))), 2))
    plot_bar_comp_sort2(flux_obs / area_new, obs_new / area_new, Moon_obs_per_value,
                        G1999_basin_new_g, continent_new, len(G1999_basin_new_g), 'Model(mol/km' + r'${^2}$' + '/yr)',
                        fig_name + '_01.jpg')
    plot_bar_comp_sort(model_ama / area_ama_new, obs_ama1 / area_ama_new, ama_name1,
                       continent_new_ama, len(obs_ama1), 'Model(mol/km' + r'${^2}$' + '/yr)', fig_name + '_02.jpg')
    plot_bar_comp_sort2(flux_obs, obs_new, Moon_obs_u,
                        G1999_basin_new_g, continent_new, len(G1999_basin_new_g), 'mol/yr', fig_name + '_03.jpg')
    plot_bar_comp_sort(model_ama, obs_ama1, ama_name1,
                       continent_new_ama, len(obs_ama1), 'mol/yr', fig_name + '_04.jpg')
    ind_tropical = np.where(region_new == 'Low')
    ind_nt = np.where(region_new == 'High')
    plot_residuals(obs_new, flux_obs, obs_new / area_new, flux_obs / area_new, ind_tropical, ind_nt,
                   fig_name + '_05.jpg')
    plot_residuals_percent(obs_new, flux_obs, obs_new / area_new, flux_obs / area_new, ind_tropical, ind_nt,
                           fig_name + '_06.jpg')
    ind_tropical_total = np.where(region_total == 'Low')
    ind_nt_total = np.where(region_total == 'High')
    plot_residuals(obs_total, model_flux, obs_total / area_total, model_flux / area_total, ind_tropical_total,
                   ind_nt_total, fig_name + '_07.jpg')
    plot_residuals_percent(obs_total, model_flux, obs_total / area_total, model_flux / area_total, ind_tropical_total,
                           ind_nt_total, fig_name + '_08.jpg')

    return R_2_global, R_2_global_log, R_2_total_log, R_2_total, weather_fCaSiO3


filename_slope = 'modern_slope.nc'
f_slope = Dataset(directory_name + filename_slope)
s = f_slope['slope'][:]
s = np.maximum(s, 0)  # m/m
filename_slope2 = '0Ma_slope.nc'
f_slope2 = Dataset(directory_name + filename_slope2)
s2 = f_slope2['slope'][:]
s2 = np.maximum(s2, 0)

filename_climate2 = 'obs_month_data_regrid05.nc'
f_clm = Dataset(directory_name + filename_climate2)
T_clm = f_clm['T'][:]
R_clm = f_clm['R'][:]
R_clm = R_clm * 1000 / (24 * 60 * 60)  # change m into mm/s
R_clm = R_clm * conv_yr_s * 1e-3  # from mm/s into m/yr
R_clm = np.maximum(R_clm, 0)
max_T = np.max(T_clm, axis=0)
min_T = np.min(T_clm, axis=0)
Tk = (max_T - min_T) / 2. - 7.5
# pdb.set_trace()
K_Tg = ((np.exp(Tk) - np.exp(-Tk)) / (np.exp(Tk) + np.exp(-Tk)) + 1.001) * 0.1  # Temperature difference

# filename_E_elm = 'erosion_elm_regrid05.nc'
# filename_E_rusle = 'erosion_rusle_2001_regrid05.nc'

filename_climate1 = 'modern_runoff_v1_interp.nc'
filename_climate3 = 'modern_temp_interp.nc'

filename_ERA1 = '50-21_year_mean_regrid05.nc'
filename_ERA2 = '50-97_year_mean_regrid05.nc'
filename_ERA3 = '50-79_year_mean_regrid05.nc'
filename_GRUN1 = 'GRUN_1902_2014.nc'
filename_GRUN2 = 'GRUN_1902_1996.nc'
filename_GRUN3 = 'GRUN_1902_1950.nc'
filename_climate4 = 'modern_runoff_v2_interp.nc'
filename_UNH = 'UNH_GRDC_runoff_regrid05.nc'

# filename_E_update_01 = 'E_update_rnf01_global.nc'
# filename_E_update_02 = 'E_update_rnf02_global.nc'
# filename_E_update_03 = 'E_update_rnf03_global.nc'
# filename_E_update_04 = 'E_update_rnf04_global.nc'
# filename_E_update_05 = 'E_update_rnf05_global.nc'
# filename_E_update_06 = 'E_update_rnf06_global.nc'
# filename_E_update_07 = 'E_update_rnf07_global.nc'
# filename_E_update_08 = 'E_update_rnf08_global.nc'
# filename_E_update_09 = 'E_update_rnf09_global.nc'

f_ERA1 = Dataset(directory_name + filename_ERA1)
f_ERA2 = Dataset(directory_name + filename_ERA2)
f_ERA3 = Dataset(directory_name + filename_ERA3)
R_ERA1 = f_ERA1['R'][:]  # m
R_ERA1 = R_ERA1 / 24 / 60 / 60 * conv_yr_s  # m/yr
R_ERA1 = np.maximum(R_ERA1, 0)
T_ERA1 = f_ERA1['T'][:]
R_ERA2 = f_ERA2['R'][:]  # m
R_ERA2 = R_ERA2 / 24 / 60 / 60 * conv_yr_s  # m/yr
R_ERA2 = np.maximum(R_ERA2, 0)
T_ERA2 = f_ERA2['T'][:]
R_ERA3 = f_ERA3['R'][:]  # m
R_ERA3 = R_ERA3 / 24 / 60 / 60 * conv_yr_s  # m/yr
R_ERA3 = np.maximum(R_ERA3, 0)
T_ERA3 = f_ERA3['T'][:]

f_GRUN1 = Dataset(directory_name + filename_GRUN1)
R_GRUN1 = f_GRUN1['Runoff'][:][0, :, :]  # mm/day
R_GRUN1 = R_GRUN1 / 1000 / 24 / 60 / 60 * conv_yr_s  # m/yr
R_GRUN1 = np.maximum(R_GRUN1, 0)
f_GRUN2 = Dataset(directory_name + filename_GRUN2)
R_GRUN2 = f_GRUN2['Runoff'][:][0, :, :]  # mm/day
R_GRUN2 = R_GRUN2 / 1000 / 24 / 60 / 60 * conv_yr_s  # m/yr
R_GRUN2 = np.maximum(R_GRUN2, 0)
f_GRUN3 = Dataset(directory_name + filename_GRUN3)
R_GRUN3 = f_GRUN3['Runoff'][:][0, :, :]  # mm/day
R_GRUN3 = R_GRUN3 / 1000 / 24 / 60 / 60 * conv_yr_s  # m/yr
R_GRUN3 = np.maximum(R_GRUN3, 0)
f_UNH = Dataset(directory_name + filename_UNH)
R_UNH = f_UNH['R'][:] / 1000  # m/yr
R_UNH = np.maximum(R_UNH, 0)

# E_u_01 = Dataset(directory_name + filename_E_update_01)['E'][:]
# E_u_02 = Dataset(directory_name + filename_E_update_02)['E'][:]
# E_u_03 = Dataset(directory_name + filename_E_update_03)['E'][:]
# E_u_04 = Dataset(directory_name + filename_E_update_04)['E'][:]
# E_u_05 = Dataset(directory_name + filename_E_update_05)['E'][:]
# E_u_06 = Dataset(directory_name + filename_E_update_06)['E'][:]
# E_u_07 = Dataset(directory_name + filename_E_update_07)['E'][:]
# E_u_08 = Dataset(directory_name + filename_E_update_08)['E'][:]
# E_u_09 = Dataset(directory_name + filename_E_update_09)['E'][:]


rou = 2500 * 1e3  # g/m3
# f_E_elm = Dataset(directory_name + filename_E_elm)
# E_elm_m = f_E_elm['erosion'][:][0:12, :, :] / rou
# E_elm = np.nanmean(E_elm_m, axis=0)
# f_E_rusle = Dataset(directory_name + filename_E_rusle)
# E_rusle = f_E_rusle['erosion_2001'][:] / rou

# filename_season = '50-21_month_mean.nc'


f_lnd = Dataset(directory_name + filename_climate1)
f_cam = Dataset(directory_name + filename_climate3)
T = f_cam['tmp'][:] + 273.15  # K
R = f_lnd['rnf'][:]  # m/yr
f_lnd4 = Dataset(directory_name + filename_climate4)
R_yves = f_lnd4['rnf'][:]  # m/yr
lat = f_lnd['lat'][:]
lon = f_lnd['lon'][:]

R = np.maximum(R, 0)
R_yves = np.maximum(R_yves, 0)
avg_runoff = np.mean(R)
n_i = len(lon)
n_j = len(lat)
avg_T = np.mean(T)

dst = np.zeros((n_j, n_i)).shape

filename_glim = 'glim_wgs84_0point5deg.txt'
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

lithology_low = np.zeros((16, 360, 720))
for i in range(nm[0]):
    for j in range(nm[1]):
        counts = {}
        p = lithology_names_1[i, j]
        counts[p] = counts.get(p, 0) + 1
        for l in range(16):
            lithology_low[l, i, j] = counts.get(l + 1, 0)

carbonate, metamorphic, sediment, mafic, felsic, ice, intermediate = convert_lith_sediment_pnas2011033117(lithology_low)
dum_lithology = np.zeros((7, carbonate.shape[0], carbonate.shape[1]))
dum_lithology[0, :, :] = carbonate
dum_lithology[1, :, :] = metamorphic
dum_lithology[2, :, :] = sediment
dum_lithology[3, :, :] = mafic
dum_lithology[4, :, :] = felsic
dum_lithology[5, :, :] = ice
dum_lithology[6, :, :] = intermediate
lithology_max = np.argmax(dum_lithology, axis=0)
lithology_max = lithology_max.astype(float)
for i in range(dst[0]):
    for j in range(dst[1]):
        if np.all(dum_lithology[:, i, j] == 0):
            lithology_max[i, j] = np.nan

# R = copy.deepcopy(R_GRUN3)
B = 1
Be = 1
sta = 1
T_list = [T, T_ERA1, T_ERA2, T_ERA3]
R_list = [R_yves, R_ERA1, R_ERA2, R_ERA3, R_GRUN1, R_GRUN2, R_GRUN3, R_UNH, R]
s_list = [s, s2]
# print(T_list)
# pdb.set_trace()
R_2_global_t = []
R_2_global_log_t = []
R_2_total_log_t = []
R_2_total_t = []
flux_t = []
# a = t_main(T_list[0], R_list[0], s_list[0], Be, sta, B, dum_lithology, 'control_case_T_' + str(0) + '_R_' + str(0) + '_s_'+ str(0))[4]
# print(a)
pdb.set_trace()
for ii in range(4):
    for jj in range(9):
        for kk in range(2):
            R_2_global, R_2_global_log, R_2_total_log, R_2_total, weather_fCaSiO3 = t_main(T_list[ii], R_list[jj], s_list[kk], Be, sta, B, dum_lithology, 'Be_case_T_' + str(ii) + '_R_' + str(jj) + '_s_'+ str(kk))
            R_2_global_t.append(R_2_global)
            R_2_global_log_t.append(R_2_global_log)
            R_2_total_log_t.append(R_2_total_log)
            R_2_total_t.append(R_2_total)
            flux_t.append(weather_fCaSiO3)
Data_mon_mean = {'R2_global': np.array(R_2_global_t), 'R2_log_global': np.array(R_2_global_log_t), 'R2_amazon_global':
    np.array(R_2_total_t), 'R2_log_amazon_global': np.array(R_2_total_log_t), 'flux_global': np.array(flux_t)}
DataFrame_mon_mean = pd.DataFrame(Data_mon_mean)
DataFrame_mon_mean.to_csv(directory_name + "Be_case.csv", index=False, sep=',')
