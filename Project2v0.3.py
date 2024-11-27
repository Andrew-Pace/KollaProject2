from pyXSteam.XSteam import XSteam
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Create an instance of the XSteam class
steam = XSteam(XSteam.UNIT_SYSTEM_MKS)  # m/kg/sec/°C/bar/W

#"Knobs"
guess_evap_tank_pressure = 0.5  # bar
guess_compressor_exit_pressure = 10  # bar
guess_mass_flow_rate_seawater_entrance = 1  # kg/s
guess_mass_flow_rate_fresh_water = 0.3 # kg/s

#not supposed to be knobs
# hexB_A = 3  # m²
# hexD_A = 3  # m²

# Define constants and input values
fresh_exit_temp = 35 #C
compressor_cost = 1000 #$/kW
HEX_cost = 2000 #$/m2
total_energy_cost = 0.11 #$/kWh
# mdot1 = mass_flow_rate_seawater_entrance
# mdot5 = mass_flow_rate_fresh_water
U = 1000  # W/m²·°C
t_seawater_in = 25  # °C
cp_water = 4180  # J/kg·°C

# Function to calculate Cmin and Cmax
def calculate_Cmin_Cmax(mass_flow_rate1, mass_flow_rate2, Cp1=cp_water, Cp2=cp_water):
    C1 = mass_flow_rate1 * Cp1
    C2 = mass_flow_rate2 * Cp2
    Cmin = min(C1, C2)
    Cmax = max(C1, C2)
    return Cmin, Cmax

# Function to calculate effectiveness for a counterflow heat exchanger
def calculate_effectiveness(NTU, Cr):
    effectiveness = (1 - np.exp(-NTU * (1 - Cr))) / (1 - Cr * np.exp(-NTU * (1 - Cr)))
    return effectiveness


def calculations(variables):
    evap_tank_pressure, compressor_exit_pressure, mdot1, mdot5 = variables

    # State 1: Seawater inlet enthalpy
    h_seawater_in = steam.h_pt(1.013, t_seawater_in)  # kJ/kg

    # State F2: Spray enthalpy and temperature
    h_spray = steam.hL_p(evap_tank_pressure)  # kJ/kg
    t_spray = steam.t_ph(evap_tank_pressure, h_spray)  # °C

    # State 5: Compressor entrance enthalpy and entropy
    t_compressor_entrance = t_spray
    s_compressor_entrance = steam.sV_p(evap_tank_pressure)
    h_compressor_entrance = steam.h_ps(evap_tank_pressure, s_compressor_entrance)  # kJ/kg

    # State 6: Compressor exit temperature and enthalpy
    t_compressor_exit = steam.tsat_p(compressor_exit_pressure)  # °C
    h_compressor_exit = steam.hV_p(compressor_exit_pressure)  # kJ/kg
    W_compressor = mdot5 * (h_compressor_exit - h_compressor_entrance)

    # State 9: Brine temperature and enthalpy
    t_brine = t_spray
    h_brine = steam.hL_p(evap_tank_pressure)  # kJ/kg

    # HEXB calculations
    CminB, CmaxB = calculate_Cmin_Cmax(mdot1, mdot1 - mdot5)
    QmaxB = ((mdot1 - mdot5) * cp_water / 1000) * (t_brine - t_seawater_in)  # kW (converted from W to kW)

    NTUB = 3
    A_hexB = NTUB / (U / CminB)
    CrB = CminB / CmaxB
    eff_hexB = calculate_effectiveness(NTUB, CrB)
    QactB = eff_hexB * QmaxB  # kW

    deltat_HEXB_hot = QactB / ((mdot1 - mdot5) * cp_water / 1000)  # °C

    h_state_3 = QactB / (0.5 * mdot1) + h_seawater_in  # kJ/kg
    h_state_2 = 2 * h_spray - h_state_3  # kJ/kg
    t_HEXB_exit = t_brine - deltat_HEXB_hot  # °C

    CminD, CmaxD = calculate_Cmin_Cmax(mdot1 * 0.5, mdot5)
    QactD = 0.5 * mdot1 * (h_state_2 - h_seawater_in)  # kW (converted from W to kW)
    deltat_HEXD_hot = QactD / (mdot5 * cp_water / 1000)  # °C
    t_HEXD_entrance = deltat_HEXD_hot + fresh_exit_temp  # °C (35°C exit for fresh water)

    if t_HEXD_entrance > steam.tsat_p(compressor_exit_pressure):
        print(
            f"Warning: t_HEXD_entrance ({t_HEXD_entrance}°C) is higher than the saturation temperature at {compressor_exit_pressure} bar.")
        return None, None

    h_state7 = steam.h_pt(compressor_exit_pressure, t_HEXD_entrance)  # kJ/kg

    if h_state7 > h_compressor_exit:
        print(f"Warning: h_state7 ({h_state7} kJ/kg) is higher than h_compressor_exit ({h_compressor_exit} kJ/kg).")
        return None, None

    NTUD = 3
    A_hexD = NTUD / (U / CminD)

    A_evap_tank = (mdot5 * (h_compressor_exit - h_state7) * 1000) / (
                U * (t_compressor_exit - t_HEXD_entrance))  # *1000 for KJ to J conversion

    LCOW_num = compressor_cost * W_compressor + HEX_cost * (A_hexB + A_hexD + A_evap_tank) + total_energy_cost * (
                W_compressor * 8790 * 15)

    if LCOW_num < 0:
        print(f"Warning: LCOW_num is negative ({LCOW_num}) for variables {variables}")
        return None, None

    LCOW_dem = (mdot5 / steam.rhoL_t(fresh_exit_temp)) * 8760 * 3600 * 15

    if LCOW_dem <= 0:
        print(f"Warning: LCOW_dem is zero or negative ({LCOW_dem}) for variables {variables}")
        return None, None

    LCOW_value = LCOW_num / LCOW_dem

    if LCOW_value < 0:
        print(f"Warning: LCOW_value is negative ({LCOW_value}) for variables {variables}")
        return None, None

    state_info = {
        "t_seawater_in_celsius": t_seawater_in,
        "h_seawater_in_kJ/kg": h_seawater_in,
        "t_spray_celsius": t_spray,
        "h_spray_kJ/kg": h_spray,
        "t_compressor_entrance_celsius": t_compressor_entrance,
        "s_compressor_entrance_kJ/(kg·K)": s_compressor_entrance,
        "h_compressor_entrance_kJ/kg": h_compressor_entrance,
        "t_compressor_exit_celsius": t_compressor_exit,
        "h_compressor_exit_kJ/kg": h_compressor_exit,
        "W_compressor_kW": W_compressor,
        "t_brine_celsius": t_brine,
        "h_brine_kJ/kg": h_brine,
        "CminB": CminB,
        "CmaxB": CmaxB,
        "QmaxB_kW": QmaxB,
        "NTUB": NTUB,
        "CrB": CrB,
        "eff_hexB": eff_hexB,
        "QactB_kW": QactB,
        "deltat_HEXB_hot": deltat_HEXB_hot,
        "h_state_3_sea_hexB_exit_kJ/kg": h_state_3,
        "h_state_2_sea_hexD_exit_kJ/kg": h_state_2,
        "t_HEXB_brinewaste_exit_celsius": t_HEXB_exit,
        "CminD": CminD,
        "CmaxD": CmaxD,
        "QactD_kW": QactD,
        "deltat_HEXD_hot": deltat_HEXD_hot,
        "t_HEXD_fresh_hot_entrance_celsius": t_HEXD_entrance,
        "h_state7_hexD_hot_entrance_kJ/kg": h_state7,
        "NTUD": NTUD,
        "A_hexB_m2": A_hexB,
        "A_hexD_m2": A_hexD,
        "A_evap_tank_m2": A_evap_tank,
        "LCOW_numerator": LCOW_num,
        "LCOW_denominator": LCOW_dem,
        "LCOW_value": LCOW_value
    }
    return LCOW_value, state_info

# Define the range of values for each variable
evaptankpressure_values = np.arange(0.1, 1.5, .1)
compressorpressure_values = np.arange(5, 20, .5)
freshwater_values = np.arange(0.25, 0.55, .05)

def excel_print(evaptankpressure_values, compressorpressure_values, freshwater_values):
    # Initialize lists to store the results
    epressure_list = []
    cpressure_list = []
    rat_list = []
    result_list = []
    state_info_list = []

    # Iterate through the range of values for each variable
    for epressure in evaptankpressure_values:
        for cpressure in compressorpressure_values:
            for rat in freshwater_values:
                LCOW_value, state_info = calculations([epressure, cpressure, guess_mass_flow_rate_seawater_entrance, rat])
                if LCOW_value is not None:
                    epressure_list.append(epressure)
                    cpressure_list.append(cpressure)
                    rat_list.append(rat)
                    result_list.append(LCOW_value)
                    state_info_list.append(state_info)

    # Convert lists to numpy arrays for easier plotting
    epressure_array = np.array(epressure_list)
    cpressure_array = np.array(cpressure_list)
    rat_array = np.array(rat_list)
    result_array = np.array(result_list)
    state_info_array = np.array(state_info_list)

    # Find the threshold value for the lowest 10% of LCOW values
    threshold_value = np.percentile(result_array, 10)

    # Filter the arrays to include only the lowest 10% of LCOW values
    mask = result_array <= threshold_value
    epressure_array_filtered = epressure_array[mask]
    cpressure_array_filtered = cpressure_array[mask]
    rat_array_filtered = rat_array[mask]
    result_array_filtered = result_array[mask]
    state_info_filtered = state_info_array[mask]

    # Create a DataFrame with the filtered results
    df_filtered_results = pd.DataFrame({
        'Evap Tank Pressure': epressure_array_filtered,
        'Compressor Pressure': cpressure_array_filtered,
        'Freshwater Ratio': rat_array_filtered,
        'LCOW Value': result_array_filtered,
        'State Info': state_info_filtered
    })

    # Sort the DataFrame by LCOW Value to get the best 10 results
    df_best_10_results = df_filtered_results.nsmallest(10, 'LCOW Value')

    # Save the best 10 results to an Excel file
    df_best_10_results.to_excel('best_10_results.xlsx', index=False)

    print("The best 10 results have been saved to 'best_10_results.xlsx'.")

def plot_grid(evaptankpressure_values, compressorpressure_values, freshwater_values):
    # Define the range of values for each variable
    evaptankpressure_values = np.arange(0.1, 1.5, .1)
    compressorpressure_values = np.arange(5, 20, .5)
    freshwater_values = np.arange(0.25, 0.55, .05)

    # Initialize lists to store the results
    epressure_list = []
    cpressure_list = []
    rat_list = []
    result_list = []

    # Iterate through the range of values for each variable
    for epressure in evaptankpressure_values:
        for cpressure in compressorpressure_values:
            for rat in freshwater_values:
                LCOW_value, state_info = calculations(
                    [epressure, cpressure, guess_mass_flow_rate_seawater_entrance, rat])
                if LCOW_value is not None:
                    epressure_list.append(epressure)
                    cpressure_list.append(cpressure)
                    rat_list.append(rat)
                    result_list.append(LCOW_value)

    # Convert lists to numpy arrays for easier plotting
    epressure_array = np.array(epressure_list)
    cpressure_array = np.array(cpressure_list)
    rat_array = np.array(rat_list)
    result_array = np.array(result_list)

    # Find the threshold value for the lowest 10% of LCOW values
    threshold_value = np.percentile(result_array, 10)

    # Filter the arrays to include only the lowest 10% of LCOW values
    mask = result_array <= threshold_value
    epressure_array_filtered = epressure_array[mask]
    cpressure_array_filtered = cpressure_array[mask]
    rat_array_filtered = rat_array[mask]
    result_array_filtered = result_array[mask]

    # Plot the results
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot
    sc = ax.scatter(epressure_array_filtered, cpressure_array_filtered, rat_array_filtered, c=result_array_filtered,
                    cmap='viridis')

    # Add color bar
    plt.colorbar(sc, label='LCOW $/m3')

    # Set labels
    ax.set_xlabel('Evap Tank Pressure')
    ax.set_ylabel('Compressor Pressure')
    ax.set_zlabel('Freshwater Ratio')
    ax.set_title('Function Results (Lowest 10% LCOW Values)')

    plt.show()


excel_print(evaptankpressure_values, compressorpressure_values, freshwater_values)
plot_grid(evaptankpressure_values, compressorpressure_values, freshwater_values)