from pyXSteam.XSteam import XSteam
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Create an instance of the XSteam class
steam = XSteam(XSteam.UNIT_SYSTEM_MKS)  # m/kg/sec/°C/bar/W

#"Knobs"
guess_evap_tank_pressure = 0.05  # bar
guess_compressor_exit_pressure = 1  # bar
guess_mass_flow_rate_seawater_entrance = 1  # kg/s
guess_mass_flow_rate_fresh_water = 0.3 # kg/s
guess_steam_ratio = guess_mass_flow_rate_fresh_water / guess_mass_flow_rate_seawater_entrance
#not supposed to be knobs
# hexB_A = 3  # m²
# hexD_A = 3  # m²

# Define constants and input values
fresh_exit_temp = 35 #C
compressor_cost = 1000 #$/kW
HEX_cost = 2000 #$/m2
energy_cost = 0.11 #$/kWh
mdot1 = guess_mass_flow_rate_seawater_entrance
mdot5 = guess_mass_flow_rate_fresh_water
U = 1000  # W/m²·°C
t_seawater_in = 25  # °C
cp_water = 4180  # J/kg·°C



def calculations_galore(variables):
    tank_pressure, compressor_exit_pressure, steam_ratio = variables
    mdotSea = mdot1
    print(tank_pressure, compressor_exit_pressure)
    mdotFresh = steam_ratio * mdotSea
    t_cold = steam.tsat_p(tank_pressure)
    t_hot = steam.tsat_p(compressor_exit_pressure)
    s_compressor = steam.sV_p(tank_pressure)
    Qdot_cold = (mdotFresh * (steam.hV_t(t_cold) - steam.hL_t(t_cold)) +
                 mdotSea * (steam.hL_t(t_cold) - steam.h_pt(tank_pressure, t_seawater_in)))
    Qdot_hot_toV = (mdotFresh * (steam.hV_t(t_hot) - steam.hL_t(t_hot)) +
                mdotFresh * (steam.h_ps(compressor_exit_pressure, s_compressor)-steam.hV_t(t_hot)))

    hC = steam.hL_p(compressor_exit_pressure) - (Qdot_cold - Qdot_hot_toV)/mdotFresh
    if hC <= 0:
        return None
    Qdot_hot = (mdotFresh * (steam.hV_t(t_hot) - steam.hL_t(t_hot)) +
                mdotFresh * (steam.h_ps(compressor_exit_pressure, s_compressor)-steam.hV_t(t_hot)) +
                mdotFresh * (steam.hL_p(compressor_exit_pressure) - hC))
    A_tank = Qdot_hot*1000 / (U * (t_hot-t_cold))
    t_fresh_exit = steam.t_ph(compressor_exit_pressure, hC)
    compressor_power = mdotFresh * (steam.h_ps(compressor_exit_pressure, s_compressor) - steam.hV_p(tank_pressure))

    LCOW = ((compressor_power * compressor_power) +
            (A_tank * HEX_cost) +
            (energy_cost * (compressor_power * 8790 * 15))) / ((mdot5 / steam.rhoL_t(t_fresh_exit)) * 8760 * 3600 * 15)

    # print(Qdot_cold, Qdot_hot, A_tank, compressor_power, LCOW)
    return t_hot, t_cold, Qdot_hot, Qdot_cold, A_tank, compressor_power, LCOW

def plotter(tps,cps,frs):
    tank_pressure_list = []
    comp_pressure_list = []
    ratio_list = []
    result_list = []
    results = []

    lowest_LCOW = float('-inf')
    lowest_specs = []

    # Iterate through the range of values for each variable
    for tp in tps:
        for cp in cps:
            for fr in frs:
                # for mdotSea in mdotSs:
                LCOW_value = calculations_galore([tp, cp, fr])
                if LCOW_value is not None:
                    # print(tp, cp, fr, LCOW_value[5])
                    # variables = (tps, cps, fr, mdotSea)
                    result = LCOW_value
                    results.append((tp, cp, fr) + result)
                    tank_pressure_list.append(tp)
                    comp_pressure_list.append(cp)
                    ratio_list.append(fr)
                    result_list.append(LCOW_value[5])
                    # if LCOW_value < lowest_LCOW:
                    #     lowest_specs = [LCOW_value, tp, cp, ]


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot
    sc = ax.scatter(tank_pressure_list, comp_pressure_list, ratio_list, c=result_list,
                    cmap='viridis')

    # Add color bar
    plt.colorbar(sc, label='LCOW $/m3')

    # Set labels
    ax.set_xlabel('Tank Pressure')
    ax.set_ylabel('Compressor Pressure')
    ax.set_zlabel('Freshwater Ratio')
    ax.set_title('SimpleSystem')
    plt.show()

    results.sort(key=lambda x: x[-1])
    top_20_results = results[:20]

    # Create a DataFrame and save to Excel
    df = pd.DataFrame(top_20_results,
                      columns=['Tank Pressure(Bar)', 'Compressor Exit Pressure(Bar)', 'Freshwater Ratio',
                               'T Hot (Celsius)','T Cold(celsius)', 'Qdot Cold(KW)', 'Qdot Hot(KW)',
                               'A Tank(m2)', 'Compressor Power(KW)', 'LCOW($/m3)'])
    df.to_excel('top_20_results.xlsx', index=False)


# def generate_excel(tank_pressures, compressor_pressures, freshwater_ratios):
#     results = []
#
#     mdotSea = 1  # Placeholder value, please replace with actual value
#
#     for tank_pressure in tank_pressures:
#         for compressor_exit_pressure in compressor_pressures:
#             for freshwater_ratio in freshwater_ratios:
#                 variables = (tank_pressure, compressor_exit_pressure, freshwater_ratio, mdotSea)
#                 result = calculations_galore(variables)
#                 results.append((tank_pressure, compressor_exit_pressure, freshwater_ratio) + result)
#
#     # Sort results by LCOW and get the top 20
#     results.sort(key=lambda x: x[-1])
#     top_20_results = results[:20]
#
#     # Create a DataFrame and save to Excel
#     df = pd.DataFrame(top_20_results,
#                       columns=['Tank Pressure', 'Compressor Exit Pressure', 'Freshwater Ratio', 'Qdot Cold', 'Qdot Hot',
#                                'A Tank', 'Compressor Power', 'LCOW'])
#     df.to_excel('top_20_results.xlsx', index=False)




tank_pressures = np.arange(0.05, 1, .05)
compressor_pressures = np.arange(1, 20, 1)
freshwater_ratios = np.arange(0.30, 0.60, .05)




plotter(tank_pressures, compressor_pressures, freshwater_ratios)
# generate_excel(tank_pressures, compressor_pressures, freshwater_ratios)


# LCOW_value = calculations_galore([.1, 1, .55,1])
