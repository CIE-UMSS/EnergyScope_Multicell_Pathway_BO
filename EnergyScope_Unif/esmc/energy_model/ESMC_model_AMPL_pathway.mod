# -------------------------------------------------------------------------------------------------------------------------													
#	EnergyScope Multi-Cell (ESMC) with Pathway Integration
#	Multi-regional energy system optimization with temporal pathway analysis
#	Based on ESMC (multi-regional single year) and PESTD (single region pathway) models
#													
#	Copyright (C) <2018-2024> <Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland and Université catholique de Louvain (UCLouvain), Belgium>
#													
#	Licensed under the Apache License, Version 2.0 (the "License");												
# -------------------------------------------------------------------------------------------------------------------------		

#########################
###  PATHWAY SETS     ###
#########################

## NEW SETS FOR PATHWAY:
set YEARS;                              # Years in the analysis (e.g., 2025, 2030, 2035...)
set PHASE;                              # Transition phases between years
set PHASE_START {PHASE} within YEARS;   # Starting year of each phase
set PHASE_STOP {PHASE} within YEARS;    # Ending year of each phase
set YEARS_WND within YEARS;             # Years in optimization window
set PHASE_WND within PHASE;             # Phases in optimization window
set YEARS_UP_TO within YEARS;           # Years up to optimization window
set PHASE_UP_TO within PHASE;           # Phases up to optimization window
set YEAR_ONE within YEARS;              # First year in optimization
set YEAR_ONE_NEXT within YEARS;         # Year after first year

## SETS FOR ANTERIOR PHASES
set ANT_INSTALL_PERIODS;                # Installation periods before optimization window
set ANT_TECHNOLOGIES;                   # Technologies installed in anterior periods

#########################
###  ORIGINAL SETS    ###
#########################

## MAIN SETS: Sets whose elements are input directly in the data file
set REGIONS;                            # regions
set RWITHOUTDAM within REGIONS;         # regions that don't have hydro dams

set PERIODS := 1 .. 8760;              # time periods (hours of the year)
set HOURS := 1 .. 24;                   # hours of the day
param nbr_tds > 0;                      # number of typical days
set TYPICAL_DAYS := 1 .. nbr_tds ordered; # typical days
set T_H_TD within {PERIODS, HOURS, TYPICAL_DAYS}; # set linking periods, hours, days, typical days
set SECTORS;                            # sectors of the energy system
set END_USES_INPUT;                     # Types of demand (end-uses). Input to the model
set END_USES_CATEGORIES;                # Categories of demand (end-uses): electricity, heat, mobility
set END_USES_TYPES_OF_CATEGORY {END_USES_CATEGORIES}; # Types of demand (end-uses).
set RESOURCES;                          # Resources: fuels (renewables and fossils) and electricity imports
set RE_FUELS within RESOURCES;          # imported biofuels.
set EXPORT within RESOURCES;            # exported resources
set NOT_LAYERS within RESOURCES;        # resources which are not a layer
set RES_IMPORT_CONSTANT within RESOURCES; # resources imported at constant rate
set END_USES_TYPES := setof {i in END_USES_CATEGORIES, j in END_USES_TYPES_OF_CATEGORY [i]} j; # secondary set
set TECHNOLOGIES_OF_END_USES_TYPE {END_USES_TYPES}; # set all energy conversion technologies
set STORAGE_TECH;                       # set of storage technologies 
set STORAGE_OF_END_USES_TYPE {END_USES_TYPES} within STORAGE_TECH; # storage tech per end-use type
set INFRASTRUCTURE;                     # Infrastructure: DHN, grid, etc.

## SECONDARY SETS
set LAYERS := (RESOURCES diff NOT_LAYERS) union END_USES_TYPES; # Layers for balancing
set TECHNOLOGIES := (setof {i in END_USES_TYPES, j in TECHNOLOGIES_OF_END_USES_TYPE [i]} j) union STORAGE_TECH union INFRASTRUCTURE; 
set TECHNOLOGIES_OF_END_USES_CATEGORY {i in END_USES_CATEGORIES} within TECHNOLOGIES := setof {j in END_USES_TYPES_OF_CATEGORY[i], k in TECHNOLOGIES_OF_END_USES_TYPE [j]} k;
set RE_RESOURCES within RESOURCES;      # List of RE resources
set NOEXCHANGES within RESOURCES;       # Resources that cannot be exchanged
set V2G within TECHNOLOGIES;            # EVs for vehicle-to-grid
set EVs_BATT within STORAGE_TECH;       # EV batteries
set EVs_BATT_OF_V2G {V2G};              # Link between EV batteries and V2G
set STORAGE_DAILY within STORAGE_TECH;  # Daily storage technologies
set TS_OF_DEC_TECH {TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}}; # Thermal storage links

## Additional SETS
set TYPICAL_DAY_OF_PERIOD {t in PERIODS} := setof {h in HOURS, td in TYPICAL_DAYS: (t,h,td) in T_H_TD} td;
set HOUR_OF_PERIOD {t in PERIODS} := setof {h in HOURS, td in TYPICAL_DAYS: (t,h,td) in T_H_TD} h;

set COGEN within TECHNOLOGIES;          # cogeneration technologies
set BOILERS within TECHNOLOGIES;        # boiler technologies

## Exchange-related sets
set EXCHANGE_FREIGHT_R within RESOURCES; # resources transported by freight
set EXCHANGE_NETWORK_R within RESOURCES; # resources exchanged through network
set EXCHANGE_R := EXCHANGE_FREIGHT_R union EXCHANGE_NETWORK_R; # all exchanged resources
set NETWORK_TYPE {EXCHANGE_NETWORK_R} within TECHNOLOGIES; # network types

# Grouped regions sets (for handling special cases like SIN in Bolivia)
set GROUPED_REGIONS := {'SIN-LL', 'SIN-VL', 'SIN-HL'};
set RELATED_PAIRS within {GROUPED_REGIONS, REGIONS} = {
    ('SIN-LL', 'SA-NO'), ('SIN-VL', 'SA-TJ'), ('SIN-LL', 'SA-SC')
};

## AGE TRACKING
set AGE {TECHNOLOGIES, PHASE} within PHASE union {"2015_2021"} union {"STILL_IN_USE"} default {"STILL_IN_USE"};

#################################
### PATHWAY PARAMETERS        ###
#################################

## Financial and temporal parameters
param max_inv_phase {PHASE} default Infinity;     # Maximum investment per phase
param t_phase {PHASE};                            # Duration of each phase [years]
param diff_2015_phase {PHASE};                    # Years from 2015 to phase
param gwp_limit_transition >= 0 default Infinity; # CO2 budget for transition
param decom_allowed {PHASE, PHASE union {"2015_2021"}, TECHNOLOGIES} default 0; # Decommissioning allowed
param remaining_years {TECHNOLOGIES, PHASE} >= 0 default 0; # Remaining years after 2050
param phase_to_year {PHASE};                      # Mapping phase to year

## Renovation/change limits
param limit_LT_renovation >= 0 default 0.33;     # Low-T heating renovation limit
param limit_pass_mob_changes >= 0 default 0.5;   # Passenger mobility change limit
param limit_freight_changes >= 0 default 0.5;    # Freight transport change limit
param limit_HT_renovation >= 0 default 0.33;     # High-T heating renovation limit
param limit_cooking_changes >= 0 default 0.33;   # Cooking change limit
param limit_mech_comm_changes >= 0 default 0.33; # Commercial mechanical change limit
param limit_mech_mov_agr_changes >= 0 default 0.33; # Mobile agricultural mechanical change limit
param limit_mech_fix_agr_changes >= 0 default 0.33; # Fixed agricultural mechanical change limit
param limit_mech_min_changes >= 0 default 0.33;  # Mining mechanical change limit
param limit_mech_fish_changes >= 0 default 0.33; # Fishing mechanical change limit

## Anterior phase parameters
param ant_install_year {ANT_INSTALL_PERIODS} default 0;
param ant_install_capacity {ANT_TECHNOLOGIES, ANT_INSTALL_PERIODS} default 0;
param ant_technology_lifespan {ANT_TECHNOLOGIES, ANT_INSTALL_PERIODS} default 25;

#################################
### ORIGINAL PARAMETERS       ###
#################################

### Parameters independent from REGIONS set (NOW WITH YEARS WHERE NEEDED)
param i_rate > 0 default 0.015;                       # discount rate
## Annualization factor
param annualised_factor {p in PHASE} := 1 / ((1 + i_rate)^diff_2015_phase[p]);
param gwp_limit_overall {YEARS} >= 0;   # maximum GWP emissions for global system
param t_op {HOURS, TYPICAL_DAYS} default 1; # operating time

# Technology and resource attributes (NOW WITH YEARS)
param c_op_exterior {YEARS, RESOURCES} >= 0;
param layers_in_out {YEARS, RESOURCES union TECHNOLOGIES diff STORAGE_TECH, LAYERS} default 0;
param gwp_op_exterior {YEARS, RESOURCES} >= 0;
param co2_net {YEARS, RESOURCES} >= 0;

# EV parameters
param batt_per_car {V2G} >= 0;
param state_of_charge_ev {V2G, HOURS} >= 0;

# Storage parameters (NOW WITH YEARS)
param storage_eff_in {YEARS, STORAGE_TECH, LAYERS} >= 0, <= 1;
param storage_eff_out {YEARS, STORAGE_TECH, LAYERS} >= 0, <= 1;
param storage_losses {YEARS, STORAGE_TECH} >= 0, <= 1;
param storage_availability {YEARS, STORAGE_TECH} >= 0, default 1;

# Solar parameters
param power_density_pv >= 0 default 0;
param power_density_solar_thermal >= 0 default 0;

# Network parameters
param loss_network {YEARS, END_USES_TYPES} >= 0 default 0;
param c_grid_extra >= 0, default 359;

# Exchange parameters (NOW WITH YEARS)
param exchange_losses {YEARS, EXCHANGE_R} >= 0 default 0;
param lhv {YEARS, EXCHANGE_FREIGHT_R} >= 0;

param total_time := sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} (t_op [h, td]);

### Parameters depending on REGIONS and YEARS
## Time series (now potentially year-dependent)
param electricity_time_series {REGIONS, HOURS, TYPICAL_DAYS} >= 0, <= 1;
param heating_time_series {REGIONS, HOURS, TYPICAL_DAYS} >= 0, <= 1;
param cooling_time_series {REGIONS, HOURS, TYPICAL_DAYS} >= 0, <= 1;
param mob_pass_time_series {REGIONS, HOURS, TYPICAL_DAYS} >= 0, <= 1;
param mob_freight_time_series {REGIONS, HOURS, TYPICAL_DAYS} >= 0, <= 1;
param c_p_t {TECHNOLOGIES, REGIONS, HOURS, TYPICAL_DAYS} >= 0 default 1;

## Demand and scenario parameters
param end_uses_demand_year {YEARS, REGIONS, END_USES_INPUT, SECTORS} >= 0 default 0;
param end_uses_input {y in YEARS, c in REGIONS, i in END_USES_INPUT} := 
    sum {s in SECTORS} (end_uses_demand_year [y, c, i, s]);
param re_share_primary {YEARS, REGIONS} >= 0;
param gwp_limit {YEARS, REGIONS} >= 0;

# Share parameters (NOW ALL WITH YEARS)
param share_mobility_public_min {YEARS, REGIONS} >= 0, <= 1;
param share_mobility_public_max {YEARS, REGIONS} >= 0, <= 1;
param share_short_haul_flights_min {YEARS, REGIONS} >= 0, <= 1;
param share_short_haul_flights_max {YEARS, REGIONS} >= 0, <= 1;
param share_freight_train_min {YEARS, REGIONS} >= 0, <= 1;
param share_freight_train_max {YEARS, REGIONS} >= 0, <= 1;
param share_freight_road_min {YEARS, REGIONS} >= 0, <= 1;
param share_freight_road_max {YEARS, REGIONS} >= 0, <= 1;
param share_freight_boat_min {YEARS, REGIONS} >= 0, <= 1;
param share_freight_boat_max {YEARS, REGIONS} >= 0, <= 1;
param share_heat_dhn_min {YEARS, REGIONS} >= 0, <= 1;
param share_heat_dhn_max {YEARS, REGIONS} >= 0, <= 1;
param share_ned {YEARS, REGIONS, END_USES_TYPES_OF_CATEGORY["NON_ENERGY"]} >= 0, <= 1;

# Technology parameters (ALL NOW WITH YEARS)
param f_max {YEARS, REGIONS, TECHNOLOGIES} >= 0 default 0;
param f_min {YEARS, REGIONS, TECHNOLOGIES} >= 0 default 0;
param fmax_perc {YEARS, REGIONS, TECHNOLOGIES} >= 0, <= 1 default 1;
param fmin_perc {YEARS, REGIONS, TECHNOLOGIES} >= 0, <= 1 default 0;
param avail_local {YEARS, REGIONS, RESOURCES} >= 0 default 0;
param avail_exterior {YEARS, REGIONS, RESOURCES} >= 0 default 0;
param c_op_local {YEARS, REGIONS, RESOURCES} >= 0 default 0;
param vehicule_capacity {YEARS, TECHNOLOGIES} >= 0, default 0;
param peak_sh_factor {REGIONS} >= 0 default 0;
param peak_sc_factor {REGIONS} >= 0 default 0;
param c_inv {YEARS, REGIONS, TECHNOLOGIES} >= 0 default 0;
param c_maint {YEARS, REGIONS, TECHNOLOGIES} >= 0 default 0;
param lifetime {YEARS, REGIONS, TECHNOLOGIES} >= 0 default 25;
param tau {y in YEARS, c in REGIONS, i in TECHNOLOGIES} := 
    i_rate * (1 + i_rate)^lifetime[y,c,i] / (((1 + i_rate)^lifetime[y,c,i]) - 1);
param gwp_constr {YEARS, REGIONS, TECHNOLOGIES} >= 0 default 0;
param gwp_op_local {YEARS, REGIONS, RESOURCES} >= 0 default 0;
param c_p {YEARS, REGIONS, TECHNOLOGIES} >= 0, <= 1 default 1;
param tc_min {YEARS, REGIONS, REGIONS, i in EXCHANGE_NETWORK_R, NETWORK_TYPE[i]} >= 0, default 0;
param tc_max {YEARS, REGIONS, REGIONS, i in EXCHANGE_NETWORK_R, NETWORK_TYPE[i]} >= 0, default 0;
param retro_gas_to_h2 {YEARS} >= 0 default 0;

# Storage regional parameters (NOW WITH YEARS)
param storage_charge_time {YEARS, REGIONS, STORAGE_TECH} >= 0;
param storage_discharge_time {YEARS, REGIONS, STORAGE_TECH} >= 0;

# Other regional parameters (NOW WITH YEARS)
param import_capacity {YEARS, REGIONS} >= 0;
param solar_area_rooftop {YEARS, REGIONS} >= 0;
param solar_area_ground {YEARS, REGIONS} >= 0;
param solar_area_ground_high_irr {YEARS, REGIONS} >= 0;
param sm_max >= 0 default 4;

# Distance and freight parameters (NOW WITH YEARS)
param dist {YEARS, REGIONS, REGIONS} >= 0 default 0;
param exch_freight {YEARS, REGIONS} >= 0 default 0;

#################################
###  VARIABLES                ###
#################################

## PATHWAY VARIABLES
var F_new {PHASE union {"2015_2021"}, REGIONS, TECHNOLOGIES} >= 0; # New capacity installed in phase
var F_decom {PHASE, PHASE union {"2015_2021"}, REGIONS, TECHNOLOGIES} >= 0; # Decommissioned capacity
var F_old {PHASE, REGIONS, TECHNOLOGIES} >= 0, default 0; # Retired capacity during phase
var F_ant_old {PHASE, REGIONS, TECHNOLOGIES} >= 0, default 0; # Retired capacity from anterior phases
var Total_F_ant_old {REGIONS, TECHNOLOGIES} >= 0, default 0; # Total retired anterior capacity

# Pathway cost variables
var C_inv_phase {PHASE} >= 0; # Total annualized investment cost per phase
var C_inv_phase_tech {PHASE, REGIONS, TECHNOLOGIES} >= 0; # Investment per tech/region/phase
var C_op_phase_tech {PHASE, REGIONS, TECHNOLOGIES} >= 0; # Operation cost per tech/region/phase
var C_op_phase_res {PHASE, REGIONS, RESOURCES} >= 0; # Resource cost per region/phase
var C_inv_return {REGIONS, TECHNOLOGIES} >= 0; # Investment return after 2050
var C_opex {YEARS} >= 0; # Annual OPEX
var C_tot_opex >= 0; # Total OPEX
var C_tot_capex >= 0; # Total CAPEX
var TotalGWPTransition >= 0; # Total transition GWP
var Delta_change {PHASE, REGIONS, TECHNOLOGIES} >= 0; # Technology changes per phase

## ORIGINAL VARIABLES (now with YEARS dimension where needed)
# Independent variables
var Share_mobility_public {y in YEARS, c in REGIONS} >= share_mobility_public_min[y,c], <= share_mobility_public_max[y,c];
var Share_short_haul_flights {y in YEARS, c in REGIONS} >= share_short_haul_flights_min[y,c], <= share_short_haul_flights_max[y,c];
var Share_freight_train {y in YEARS, c in REGIONS} >= share_freight_train_min[y,c], <= share_freight_train_max[y,c];
var Share_freight_road {y in YEARS, c in REGIONS} >= share_freight_road_min[y,c], <= share_freight_road_max[y,c];
var Share_freight_boat {y in YEARS, c in REGIONS} >= share_freight_boat_min[y,c], <= share_freight_boat_max[y,c];
var Share_heat_dhn {y in YEARS, c in REGIONS} >= share_heat_dhn_min[y,c], <= share_heat_dhn_max[y,c];

var F {YEARS, REGIONS, TECHNOLOGIES} >= 0; # Installed capacity
var F_t {YEARS, REGIONS, TECHNOLOGIES, HOURS, TYPICAL_DAYS} >= 0; # Operation in each period
var R_t_local {YEARS, REGIONS, RESOURCES, HOURS, TYPICAL_DAYS} >= 0; # Local resource use
var R_t_exterior {YEARS, REGIONS, RESOURCES, HOURS, TYPICAL_DAYS} >= 0; # Exterior resource import
var Storage_in {YEARS, REGIONS, i in STORAGE_TECH, LAYERS, HOURS, TYPICAL_DAYS} >= 0;
var Storage_out {YEARS, REGIONS, i in STORAGE_TECH, LAYERS, HOURS, TYPICAL_DAYS} >= 0;
var Transfer_capacity {YEARS, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_NETWORK_R, n in NETWORK_TYPE[i]} >= 0;
var Exch_imp {YEARS, REGIONS, REGIONS, EXCHANGE_R, HOURS, TYPICAL_DAYS} >= 0;
var Exch_exp {YEARS, REGIONS, REGIONS, EXCHANGE_R, HOURS, TYPICAL_DAYS} >= 0;

var Power_nuclear {YEARS, REGIONS} >= 0;
var Shares_mobility_passenger {YEARS, REGIONS, TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_PASSENGER"]} >= 0;
var Shares_mobility_freight {YEARS, REGIONS, TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_FREIGHT"]} >= 0;
var Shares_shipping {YEARS, REGIONS, TECHNOLOGIES_OF_END_USES_CATEGORY["SHIPPING"]} >= 0;
var Shares_lowT_dec {YEARS, REGIONS, TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}} >= 0;
var F_solar {YEARS, REGIONS, TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}} >= 0;
var F_t_solar {YEARS, REGIONS, TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}, h in HOURS, td in TYPICAL_DAYS} >= 0;

# Dependent variables
var End_uses {YEARS, REGIONS, LAYERS, HOURS, TYPICAL_DAYS} >= 0;
var TotalCost {YEARS, REGIONS} >= 0;
var C_inv {YEARS, REGIONS, TECHNOLOGIES} >= 0;
var C_maint {YEARS, REGIONS, TECHNOLOGIES} >= 0;
var C_op {YEARS, REGIONS, RESOURCES} >= 0;
var TotalGWP {YEARS, REGIONS} >= 0;
var GWP_constr {YEARS, REGIONS, TECHNOLOGIES} >= 0;
var GWP_op {YEARS, REGIONS, RESOURCES} >= 0;
var CO2_net {YEARS, REGIONS, RESOURCES} >= 0;
var Curt {YEARS, REGIONS, TECHNOLOGIES, HOURS, TYPICAL_DAYS} >= 0;
var Network_losses {YEARS, REGIONS, END_USES_TYPES, HOURS, TYPICAL_DAYS} >= 0;
var Storage_level {YEARS, REGIONS, STORAGE_TECH, PERIODS} >= 0;
var R_t_import {YEARS, REGIONS, RESOURCES, HOURS, TYPICAL_DAYS} >= 0;
var R_t_export {YEARS, REGIONS, RESOURCES, HOURS, TYPICAL_DAYS} >= 0;
var Exch_freight_border {YEARS, REGIONS, REGIONS} >= 0;
var Exch_freight {YEARS, REGIONS} >= 0;
var Import_constant {YEARS, REGIONS, RES_IMPORT_CONSTANT} >= 0;
var Exch_network_imp {YEARS, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_NETWORK_R, h in HOURS, td in TYPICAL_DAYS} >= 0 default 0;
var Exch_network_exp {YEARS, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_NETWORK_R, h in HOURS, td in TYPICAL_DAYS} >= 0 default 0;
var Exch_extra_freight_imp {YEARS, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_FREIGHT_R, h in HOURS, td in TYPICAL_DAYS} >= 0 default 0;
var Exch_extra_freight_exp {YEARS, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_FREIGHT_R, h in HOURS, td in TYPICAL_DAYS} >= 0 default 0;

#########################################
###      CONSTRAINTS                  ###
#########################################

## End-uses demand calculation constraints 
subject to end_uses_t {y in YEARS_WND diff YEAR_ONE, c in REGIONS, l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
	End_uses [y, c, l, h, td] = (if l == "ELECTRICITY"
		then
			(end_uses_input[y,c,l] * electricity_time_series [c, h, td] / t_op [h, td] ) + Network_losses [y,c,l,h,td]
		else (if l == "LIGHTING_R_C" then
			(end_uses_input[y,c,"LIGHTING_R_C"] * electricity_time_series [c, h, td] / t_op [h, td])
		else (if l == "LIGHTING_P" then
			(end_uses_input[y,c,"LIGHTING_P"] * electricity_time_series [c, h, td] / t_op [h, td])
		else (if l == "HEAT_LOW_T_DHN" then
			(end_uses_input[y,c,"HEAT_LOW_T_HW"] / total_time + end_uses_input[y,c,"HEAT_LOW_T_SH"] * heating_time_series [c, h, td] / t_op [h, td] ) * Share_heat_dhn[y,c] + Network_losses [y,c,l,h,td]
		else (if l == "HEAT_LOW_T_DECEN" then
			(end_uses_input[y,c,"HEAT_LOW_T_HW"] / total_time + end_uses_input[y,c,"HEAT_LOW_T_SH"] * heating_time_series [c, h, td] / t_op [h, td] ) * (1 - Share_heat_dhn[y,c])
		else (if l == "MOB_PUBLIC" then
			(end_uses_input[y,c,"MOBILITY_PASSENGER"] * mob_pass_time_series [c, h, td] / t_op [h, td]  ) * Share_mobility_public[y,c]
		else (if l == "AVIATION_SHORT_HAUL" then
			(end_uses_input[y,c,"MOBILITY_PASSENGER"] * mob_pass_time_series [c, h, td] / t_op [h, td]  ) * Share_short_haul_flights[y,c]
		else (if l == "AVIATION_LONG_HAUL" then
			(end_uses_input[y,c,"AVIATION_LONG_HAUL"] * mob_pass_time_series [c, h, td] / t_op [h, td]  ) 
		else (if l == "MOB_PRIVATE" then
			(end_uses_input[y,c,"MOBILITY_PASSENGER"] * mob_pass_time_series [c, h, td] / t_op [h, td]  ) * (1 - Share_mobility_public[y,c]-Share_short_haul_flights[y,c])
		else (if l == "MOB_FREIGHT_RAIL" then
			(end_uses_input[y,c,"MOBILITY_FREIGHT"]   * mob_freight_time_series [c, h, td] / t_op [h, td] ) *  Share_freight_train[y,c]
		else (if l == "MOB_FREIGHT_ROAD" then
			((end_uses_input[y,c,"MOBILITY_FREIGHT"]   * mob_freight_time_series [c, h, td] / t_op [h, td] ) *  Share_freight_road[y,c] + exch_freight [y,c] / total_time)
		else (if l == "MOB_FREIGHT_BOAT" then
			(end_uses_input[y,c,"MOBILITY_FREIGHT"]   * mob_freight_time_series [c, h, td] / t_op [h, td] ) *  Share_freight_boat[y,c]
		else (if l == "SHIPPING" then
			end_uses_input[y,c,l] / total_time
		else (if l == "HEAT_HIGH_T" then
			end_uses_input[y,c,l] / total_time
		else (if l == "FOOD_PRESERVATION" then
			end_uses_input[y,c,l] / total_time
		else (if l == "COOKING" then
			end_uses_input[y,c,l] / total_time
		else (if l == "SPACE_COOLING" then
			end_uses_input[y,c,l] * cooling_time_series [c, h, td] / t_op [h, td]
		else (if l == "PROCESS_COOLING" then
			end_uses_input[y,c,l] / total_time
		else (if l == "MECHANICAL_ENERGY_COMM" then
			end_uses_input[y,c,l] / total_time
		else (if l == "MECHANICAL_ENERGY_IND" then
			end_uses_input[y,c,l] / total_time
		else (if l == "MECHANICAL_ENERGY_MOV_AGR" then
			end_uses_input[y,c,l] / total_time
		else (if l == "MECHANICAL_ENERGY_FIX_AGR" then
			end_uses_input[y,c,l] / total_time
		else (if l == "MECHANICAL_ENERGY_MIN" then
			end_uses_input[y,c,l] / total_time
		else (if l == "MECHANICAL_ENERGY_FISH_OTHERS" then
			end_uses_input[y,c,l] / total_time			
		else (if l == "HVC" then
			end_uses_input[y,c,"NON_ENERGY"] * share_ned [y,c, "HVC"] / total_time
		else (if l == "AMMONIA" then
			end_uses_input[y,c, "NON_ENERGY"] * share_ned [y,c, "AMMONIA"] / total_time
		else (if l == "METHANOL" then
			end_uses_input[y,c, "NON_ENERGY"] * share_ned [y,c, "METHANOL"] / total_time
		else 
			0 ))))))))))))))))))))))))))); # For all layers which don't have an end-use demand

## Cost constraints
subject to totalcost_cal{y in YEARS_UP_TO union YEARS_WND, c in REGIONS}:
	TotalCost[y,c] = sum {j in TECHNOLOGIES} (tau [y,c,j]  * C_inv [y,c,j] + C_maint [y,c,j]) + sum {i in RESOURCES} C_op [y,c,i];
	
subject to investment_cost_calc {y in YEARS_UP_TO union YEARS_WND, c in REGIONS, j in TECHNOLOGIES}: 
	C_inv [y,c,j] = c_inv [y,c,j] * F [y,c,j];
		
subject to main_cost_calc {y in YEARS_UP_TO union YEARS_WND, c in REGIONS, j in TECHNOLOGIES}: 
	C_maint [y,c,j] = c_maint [y,c,j] * F [y,c,j];		

subject to op_cost_calc {y in YEARS_UP_TO union YEARS_WND, c in REGIONS, i in RESOURCES}:
	C_op [y,c,i] = sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} 
		(c_op_local [y,c, i] * R_t_local [y,c, i, h, td] * t_op [h, td] + c_op_exterior [y, i] * R_t_exterior [y,c, i, h, td] * t_op [h, td] ) ;

## Emissions constraints
subject to totalGWP_calc {y in YEARS_UP_TO union YEARS_WND, c in REGIONS}:
	TotalGWP[y,c] = sum {j in TECHNOLOGIES} (GWP_constr [y,c,j] / lifetime [y,c,j]) + sum {i in RESOURCES} GWP_op [y,c,i];
	
subject to gwp_constr_calc {y in YEARS_UP_TO union YEARS_WND, c in REGIONS, j in TECHNOLOGIES}:
	GWP_constr [y,c,j] = gwp_constr [y,c,j] * F [y,c,j];

subject to gwp_op_calc {y in YEARS_UP_TO union YEARS_WND, c in REGIONS, i in RESOURCES}:
	GWP_op [y,c,i] = sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} 
		(R_t_local [y,c, i, h, td] * gwp_op_local [y,c, i] * t_op [h, td] + R_t_exterior [y,c, i, h, td] * gwp_op_exterior [y, i] * t_op [h, td] );	

subject to co2_net_calc {y in YEARS_UP_TO union YEARS_WND, c in REGIONS, i in RESOURCES}:
	CO2_net [y,c,i] = sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} 
		(R_t_local [y,c, i, h, td] * co2_net [y, i] * t_op [h, td] + R_t_exterior [y,c, i, h, td] * co2_net [y, i] * t_op [h, td] );

## Capacity and operation constraints
subject to size_limit {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES}:
	f_min [y,c,j] <= F [y,c,j] <= f_max [y,c,j];

subject to capacity_factor_t {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES, h in HOURS, td in TYPICAL_DAYS}:
	F_t [y,c,j, h, td] + Curt[y,c, j, h, td] = F [y,c,j] * c_p_t [j, c, h, td];
	
subject to capacity_factor {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD [t], td in TYPICAL_DAY_OF_PERIOD [t]} (F_t [y,c, j, h, td] * t_op [h, td]) <= F [y,c, j] * c_p [y,c, j] * total_time;

## Resources constraints
subject to resource_availability_local {y in YEARS_WND diff YEAR_ONE, c in REGIONS, i in RESOURCES}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (R_t_local [y,c, i, h, td] * t_op [h, td]) <= avail_local [y,c, i];
	
subject to resource_availability_exterior {y in YEARS_WND diff YEAR_ONE, c in REGIONS, i in RESOURCES}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (R_t_exterior [y,c, i, h, td] * t_op [h, td]) <= avail_exterior [y,c, i];
	
subject to resource_constant_import {y in YEARS_WND diff YEAR_ONE, c in REGIONS, i in RES_IMPORT_CONSTANT, h in HOURS, td in TYPICAL_DAYS}:
	R_t_exterior [y,c, i, h, td] * t_op [h, td] = Import_constant [y,c, i];

## Layer balance
subject to layer_balance {y in YEARS_WND diff YEAR_ONE, c in REGIONS, l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
	sum {i in RESOURCES} (layers_in_out[y, i, l] * (R_t_local [y,c, i, h, td] + R_t_exterior [y,c, i, h, td] - R_t_export[y,c, i, h, td] + R_t_import[y,c, i, h, td])) 
	+ sum {k in TECHNOLOGIES diff STORAGE_TECH} (layers_in_out[y, k, l] * F_t [y,c, k, h, td]) 
	+ sum {j in STORAGE_TECH} ( Storage_out [y,c, j, l, h, td] - Storage_in [y,c, j, l, h, td] )
	- End_uses [y,c, l, h, td]
	= 0;
		
## Storage constraints	
subject to storage_level {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in STORAGE_TECH, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]}:
	Storage_level [y,c, j, t] = (if t == 1 then
		Storage_level [y,c, j, card(PERIODS)] * (1.0 -  storage_losses[y, j])
		+ t_op [h, td] * (   (sum {l in LAYERS: storage_eff_in [y, j,l] > 0}  (Storage_in [y,c, j, l, h, td]  * storage_eff_in  [y, j, l])) 
		                   - (sum {l in LAYERS: storage_eff_out [y, j,l] > 0} (Storage_out [y,c, j, l, h, td] / storage_eff_out [y, j, l])))
	else
		Storage_level [y,c, j, t-1] * (1.0 -  storage_losses[y, j])
		+ t_op [h, td] * (   (sum {l in LAYERS: storage_eff_in [y, j,l] > 0}  (Storage_in [y,c, j, l, h, td]  * storage_eff_in  [y, j, l])) 
		                   - (sum {l in LAYERS: storage_eff_out [y, j,l] > 0} (Storage_out [y,c, j, l, h, td] / storage_eff_out [y, j, l])))
		);

subject to impose_daily_storage {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in STORAGE_DAILY, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]}:
	Storage_level [y,c, j, t] = F_t [y,c, j, h, td];

subject to impose_non_daily_storage {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in STORAGE_TECH diff STORAGE_DAILY, h in HOURS, td in TYPICAL_DAYS}:
	F_t [y,c, j, h, td] = 0;
	
subject to limit_energy_stored_to_maximum {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in STORAGE_TECH diff STORAGE_DAILY , t in PERIODS}:
	Storage_level [y,c, j, t] <= F [y,c, j];
	
subject to storage_layer_in {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in STORAGE_TECH, l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
	Storage_in [y,c, j, l, h, td] * (ceil (storage_eff_in [y, j, l]) - 1) = 0;

subject to storage_layer_out {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in STORAGE_TECH, l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
	Storage_out [y,c, j, l, h, td] * (ceil (storage_eff_out [y, j, l]) - 1) = 0;

subject to limit_energy_to_power_ratio {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in STORAGE_TECH diff {"BEV_BATT","PHEV_BATT","SUV_ELEC_PRIVATE_BATT","MOTORCYCLE_ELEC_PRIVATE_BATT","BUS_ELEC_PUBLIC_BATT"}, l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
	Storage_in [y,c, j, l, h, td] * storage_charge_time[y,c,j] + Storage_out [y,c, j, l, h, td] * storage_discharge_time[y,c,j] <=  F [y,c, j] * storage_availability[y, j];

## Exchange constraints
subject to reciprocity_of_Exchanges {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_R, h in HOURS, td in TYPICAL_DAYS} :
	Exch_imp[y,c1,c2,i,h,td] * (1 + exchange_losses[y, i]* dist[y, c1, c2]/1000) - Exch_exp[y,c1,c2,i,h,td] = - Exch_imp[y,c2,c1,i,h,td] * (1 + exchange_losses[y, i] * dist[y, c2, c1]/1000) + Exch_exp[y,c2,c1,i,h,td];

subject to importation {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, i in EXCHANGE_R, h in HOURS, td in TYPICAL_DAYS}:
	R_t_import[y,c1, i, h, td]  = sum{c2 in REGIONS} Exch_imp[y,c1,c2,i,h,td];

subject to exportation {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, i in EXCHANGE_R, h in HOURS, td in TYPICAL_DAYS}:
	R_t_export[y,c1, i, h, td]  = sum{c2 in REGIONS} Exch_exp[y,c1,c2,i,h,td];
	
subject to exchanges_only_between_neighbours {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_R, h in HOURS, td in TYPICAL_DAYS : dist[y, c1, c2] == 0} :
	Exch_imp[y,c1,c2,i,h,td] = 0;

subject to exchanges_only_between_neighbours2 {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_R, h in HOURS, td in TYPICAL_DAYS : dist[y, c1, c2] == 0} :
	Exch_exp[y,c1,c2,i,h,td] = 0;

subject to resources_no_exchanges {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, n in NOEXCHANGES, h in HOURS, td in TYPICAL_DAYS} :
	R_t_import [y,c1, n, h, td] = 0;

subject to resources_no_exchanges2 {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, n in NOEXCHANGES, h in HOURS, td in TYPICAL_DAYS} :
	R_t_export [y,c1, n, h, td] = 0;

subject to exch_imp_def {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_R, h in HOURS, td in TYPICAL_DAYS} :
    Exch_imp[y,c1,c2,i,h,td] =
        (if i in EXCHANGE_NETWORK_R then Exch_network_imp[y,c1,c2,i,h,td] else 0)
      + (if i in EXCHANGE_FREIGHT_R then Exch_extra_freight_imp[y,c1,c2,i,h,td] else 0);

subject to exch_exp_def {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_R, h in HOURS, td in TYPICAL_DAYS} :
    Exch_exp[y,c1,c2,i,h,td] =
        (if i in EXCHANGE_NETWORK_R then Exch_network_exp[y,c1,c2,i,h,td] else 0)
      + (if i in EXCHANGE_FREIGHT_R then Exch_extra_freight_exp[y,c1,c2,i,h,td] else 0);

## Network exchanges
subject to network_capacity_limit_imp {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_NETWORK_R, h in HOURS, td in TYPICAL_DAYS} :
    Exch_network_imp[y,c1,c2,i,h,td] <= sum{n in NETWORK_TYPE[i]} Transfer_capacity[y,c2,c1,i,n];

subject to network_capacity_limit_exp {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_NETWORK_R, h in HOURS, td in TYPICAL_DAYS} :
    Exch_network_exp[y,c1,c2,i,h,td] <= sum{n in NETWORK_TYPE[i]} Transfer_capacity[y,c1,c2,i,n];

subject to transfer_capacity_bounds {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_NETWORK_R, n in NETWORK_TYPE[i]}:
	tc_min[y, c1, c2, i, n] <= Transfer_capacity[y,c1, c2, i, n] <= tc_max[y, c1, c2, i, n];

subject to bidirectonal_exchanges {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_NETWORK_R, n in NETWORK_TYPE[i]}:
	Transfer_capacity [y,c1,c2,i,n] = Transfer_capacity [y,c2,c1,i,n];

# Gas pipeline retrofitting to hydrogen
subject to transfer_capacity_bounds_gas_pipeline {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS}:
	tc_min[y, c1, c2, "GAS", "GAS_PIPELINE"] <=
	Transfer_capacity[y,c1, c2, "GAS", "GAS_PIPELINE"] + Transfer_capacity[y,c1, c2, "H2", "H2_RETROFITTED"] / retro_gas_to_h2[y]
	<= tc_max[y, c1, c2, "GAS", "GAS_PIPELINE"];

subject to transfer_capacity_bounds_gas_subsea {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS}:
	tc_min[y, c1, c2, "GAS", "GAS_SUBSEA"] <=
	Transfer_capacity[y,c1, c2, "GAS", "GAS_SUBSEA"] + Transfer_capacity[y,c1, c2, "H2", "H2_SUBSEA_RETRO"] / retro_gas_to_h2[y]
	<= tc_max[y, c1, c2, "GAS", "GAS_SUBSEA"];

subject to networks_infra {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, i in EXCHANGE_NETWORK_R, n in NETWORK_TYPE[i]}:
	F[y,c1, n] >= sum{c2 in REGIONS} (dist[y, c1, c2] * Transfer_capacity [y,c1,c2,i,n]/2);

# Freight exchanges
subject to exch_imp_freight_def {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_FREIGHT_R, h in HOURS, td in TYPICAL_DAYS} :
    Exch_extra_freight_imp[y,c1,c2,i,h,td] >= Exch_imp[y,c1,c2,i,h,td] -
        (if i in EXCHANGE_NETWORK_R then sum{n in NETWORK_TYPE[i]} Transfer_capacity[y,c2,c1,i,n] else 0);

subject to define_exch_extra_freight_imp_ub {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_FREIGHT_R, h in HOURS, td in TYPICAL_DAYS}:
    Exch_extra_freight_imp[y,c1,c2,i,h,td] <= Exch_imp[y,c1,c2,i,h,td];

subject to exch_imp_freight_def_non_negative {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_FREIGHT_R, h in HOURS, td in TYPICAL_DAYS} :
    Exch_extra_freight_imp[y,c1,c2,i,h,td] >= 0;

subject to exch_exp_freight_def {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_FREIGHT_R, h in HOURS, td in TYPICAL_DAYS} :
    Exch_extra_freight_exp[y,c1,c2,i,h,td] >= Exch_exp[y,c1,c2,i,h,td] -
        (if i in EXCHANGE_NETWORK_R then sum{n in NETWORK_TYPE[i]} Transfer_capacity[y,c1,c2,i,n] else 0);

subject to define_exch_extra_freight_exp_ub {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_FREIGHT_R, h in HOURS, td in TYPICAL_DAYS}:
    Exch_extra_freight_exp[y,c1,c2,i,h,td] <= Exch_exp[y,c1,c2,i,h,td];

subject to exch_exp_freight_def_non_negative {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS, i in EXCHANGE_FREIGHT_R, h in HOURS, td in TYPICAL_DAYS} :
    Exch_extra_freight_exp[y,c1,c2,i,h,td] >= 0;

subject to freight_of_exchanges_border {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS, c2 in REGIONS} :
    Exch_freight_border[y,c1,c2] = dist[y, c1,c2] * sum{r in EXCHANGE_FREIGHT_R, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]}
    	((Exch_extra_freight_imp[y,c1,c2,r,h,td] + Exch_extra_freight_exp[y,c1,c2,r,h,td]) / lhv[y, r]);

subject to limit_exch_freight_positive {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS: end_uses_input[y,c1, "MOBILITY_FREIGHT"] > 0} :
    Exch_freight[y,c1] <= end_uses_input[y,c1, "MOBILITY_FREIGHT"] * Share_freight_road[y,c1];

subject to freight_of_exchanges_grouped {y in YEARS_WND diff YEAR_ONE, c1 in GROUPED_REGIONS}:
    Exch_freight[y,c1] =
        sum {c2 in GROUPED_REGIONS} Exch_freight_border[y,c1, c2] / 2
      + sum {c2 in REGIONS diff GROUPED_REGIONS} (Exch_freight_border[y,c1, c2] + Exch_freight_border[y,c2, c1])
      + sum {c2 in REGIONS diff GROUPED_REGIONS, c3 in REGIONS diff GROUPED_REGIONS:
             (c1, c2) in RELATED_PAIRS} Exch_freight_border[y,c2, c3];

subject to freight_of_exchanges_not_grouped {y in YEARS_WND diff YEAR_ONE, c1 in REGIONS diff GROUPED_REGIONS}:
    Exch_freight[y,c1] = 0;

## Local networks constraints
subject to network_losses {y in YEARS_WND diff YEAR_ONE, c in REGIONS, eut in END_USES_TYPES, h in HOURS, td in TYPICAL_DAYS}:
	Network_losses [y,c,eut,h,td] = (sum {j in TECHNOLOGIES diff STORAGE_TECH: layers_in_out [y, j, eut] > 0} ((layers_in_out[y, j, eut]) * F_t [y,c, j, h, td])) * loss_network [y,eut] 
		+ (sum {i in RESOURCES} ((layers_in_out[y, i, eut]) * R_t_exterior[y,c, i, h, td])) * loss_network [y,eut];

subject to extra_grid{y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
    F [y,c,"GRID"] = 1 + (c_grid_extra / c_inv[y,c,"GRID"]) *( (F [y,c,"WIND_ONSHORE"] + F [y,c,"WIND_OFFSHORE"] + F [y,c,"PV_ROOFTOP"] + F [y,c,"PV_UTILITY"])
					                                     - (f_min [y,c,"WIND_ONSHORE"] + f_min [y,c,"WIND_OFFSHORE"] + f_min [y,c,"PV_ROOFTOP"] + f_min [y,c,"PV_UTILITY"]) );

subject to extra_dhn{y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
	F [y,c,"DHN"] = sum {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"]} (F [y,c,j]);
	
## Mobility shares
subject to operating_strategy_mob_passenger{y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_PASSENGER"], h in HOURS, td in TYPICAL_DAYS}:
	F_t [y,c, j, h, td]   = Shares_mobility_passenger [y,c,j] * (end_uses_input[y,c,"MOBILITY_PASSENGER"] * mob_pass_time_series [c, h, td] / t_op [h, td] );

subject to operating_strategy_mobility_freight{y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_FREIGHT"], h in HOURS, td in TYPICAL_DAYS}:
	F_t [y,c, j, h, td]   = Shares_mobility_freight [y,c,j] * (end_uses_input[y,c,"MOBILITY_FREIGHT"] * mob_freight_time_series [c, h, td] / t_op [h, td] );
	
subject to operating_strategy_shipping{y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES_OF_END_USES_CATEGORY["SHIPPING"], h in HOURS, td in TYPICAL_DAYS}:
	F_t [y,c, j, h, td]   = Shares_shipping [y,c,j] * (end_uses_input[y,c,"SHIPPING"] / total_time);

subject to Freight_shares {y in YEARS_WND diff YEAR_ONE, c in REGIONS} :
	Share_freight_train[y,c] + Share_freight_road[y,c] + Share_freight_boat[y,c] = 1;
	
## V2G constraints
subject to EV_storage_for_V2G_demand {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in V2G, i in EVs_BATT_OF_V2G[j], h in HOURS, td in TYPICAL_DAYS}:
	Storage_out [y,c,i,"ELECTRICITY",h,td] >=  - layers_in_out[y, j,"ELECTRICITY"]* F_t [y,c, j, h, td];

subject to EV_storage_size {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in V2G, i in EVs_BATT_OF_V2G[j]}:
	F [y,c,i] = F[y,c,j] / vehicule_capacity [y,j] * batt_per_car[j];

subject to limit_energy_to_power_ratio_bis {y in YEARS_WND diff YEAR_ONE, c in REGIONS, i in V2G, j in EVs_BATT_OF_V2G[i], l in LAYERS, h in HOURS, td in TYPICAL_DAYS}:
    Storage_in [y,c, j, l, h, td] * storage_charge_time[y,c,j] + (Storage_out [y,c, j, l, h, td] + layers_in_out[y, i,"ELECTRICITY"]* F_t [y,c, i, h, td]) * storage_discharge_time[y,c,j]  
    	<= ( F [y,c, j] - F_t [y,c,i,h,td] / vehicule_capacity [y,i] * batt_per_car[i] ) * storage_availability[y, j];
	
subject to EV_storage_min_SOC {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in V2G, i in EVs_BATT_OF_V2G[j], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]}:
	Storage_level [y,c, i, t] >= F [y,c,i] * state_of_charge_ev[j,h];

## Hydro dams
subject to storage_level_hydro_dams {y in YEARS_WND diff YEAR_ONE, c in REGIONS diff RWITHOUTDAM}:
	F [y,c,"DAM_STORAGE"] <= f_min [y,c,"DAM_STORAGE"] + (f_max [y,c,"DAM_STORAGE"]-f_min [y,c,"DAM_STORAGE"]) * (F [y,c,"HYDRO_DAM"] - f_min [y,c,"HYDRO_DAM"])/(f_max [y,c,"HYDRO_DAM"] - f_min [y,c,"HYDRO_DAM"]);

subject to impose_hydro_dams_inflow {y in YEARS_WND diff YEAR_ONE, c in REGIONS, h in HOURS, td in TYPICAL_DAYS}:
	Storage_in [y,c, "DAM_STORAGE", "ELECTRICITY", h, td] = F_t [y,c, "HYDRO_DAM", h, td];

subject to limit_hydro_dams_output {y in YEARS_WND diff YEAR_ONE, c in REGIONS, h in HOURS, td in TYPICAL_DAYS}:
	Storage_out [y,c, "DAM_STORAGE", "ELECTRICITY", h, td] <= F [y,c,"HYDRO_DAM"];

## CSP constraints
subject to sm_limit_solar_tower {y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
	-F[y,c,"ST_COLLECTOR"]/layers_in_out [y, "ST_POWER_BLOCK", "ST_HEAT"] <= sm_max * F[y,c,"ST_POWER_BLOCK"];

subject to sm_limit_parabolic_trough {y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
	-F[y,c,"PT_COLLECTOR"]/layers_in_out [y, "PT_POWER_BLOCK", "PT_HEAT"] <= sm_max * F[y,c,"PT_POWER_BLOCK"];
	
## Solar area constraints
subject to solar_area_rooftop_limited {y in YEARS_WND diff YEAR_ONE, c in REGIONS} :
	(F[y,c,"PV_ROOFTOP"])/power_density_pv +(F[y,c,"DEC_SOLAR"]+F[y,c,"DHN_SOLAR"])/power_density_solar_thermal <= solar_area_rooftop [y,c];

subject to solar_area_ground_limited {y in YEARS_WND diff YEAR_ONE, c in REGIONS} :
	(F[y,c,"PV_UTILITY"])/power_density_pv
		+ (-F[y,c,"PT_COLLECTOR"]/layers_in_out [y, "PT_POWER_BLOCK", "PT_HEAT"] - F[y,c,"ST_COLLECTOR"]/layers_in_out [y, "ST_POWER_BLOCK", "ST_HEAT"])/power_density_pv
	<= solar_area_ground [y,c];

subject to solar_area_ground_high_irr_limited {y in YEARS_WND diff YEAR_ONE, c in REGIONS} :
	(-F[y,c,"PT_COLLECTOR"]/layers_in_out [y, "PT_POWER_BLOCK", "PT_HEAT"]
		-F[y,c,"ST_COLLECTOR"]/layers_in_out [y, "ST_POWER_BLOCK", "ST_HEAT"])/power_density_pv
	<= solar_area_ground_high_irr [y,c];

## Decentralised heat production
subject to thermal_solar_capacity_factor {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}, h in HOURS, td in TYPICAL_DAYS}:
	F_t_solar [y,c, j, h, td] <= F_solar[y,c,j] * c_p_t["DEC_SOLAR", c, h, td];
	
subject to thermal_solar_total_capacity {y in YEARS_WND diff YEAR_ONE, c in REGIONS} :
	F [y,c,"DEC_SOLAR"] = sum {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}} F_solar[y,c,j];

subject to decentralised_heating_balance  {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}, i in TS_OF_DEC_TECH[j], h in HOURS, td in TYPICAL_DAYS}:
	F_t [y,c, j, h, td] + F_t_solar [y,c, j, h, td] + sum {l in LAYERS } ( Storage_out [y,c, i, l, h, td] - Storage_in [y,c, i, l, h, td])  
		= Shares_lowT_dec[y,c,j] * (end_uses_input[y,c,"HEAT_LOW_T_HW"] / total_time + end_uses_input[y,c,"HEAT_LOW_T_SH"] * heating_time_series [c, h, td] / t_op [h, td]);

## Peak demand constraints
subject to peak_lowT_dec {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"] diff {"DEC_SOLAR"}, h in HOURS, td in TYPICAL_DAYS}:
	F [y,c,j] >= peak_sh_factor[c] * F_t [y,c, j, h, td] ;

var Max_Heat_Demand{YEARS, REGIONS} >= 0;
subject to max_dhn_heat_demand {y in YEARS_WND diff YEAR_ONE, c in REGIONS, h in HOURS, td in TYPICAL_DAYS}:
	Max_Heat_Demand[y,c] >= End_uses [y,c,"HEAT_LOW_T_DHN", h, td];

subject to peak_lowT_dhn {y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
	sum {j in TECHNOLOGIES_OF_END_USES_TYPE ["HEAT_LOW_T_DHN"], i in STORAGE_OF_END_USES_TYPE["HEAT_LOW_T_DHN"]} (F [y,c,j] + F[y,c,i]/storage_discharge_time[y,c,i]) >= peak_sh_factor[c] * Max_Heat_Demand[y,c];
	
subject to peak_space_cooling {y in YEARS_WND diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES_OF_END_USES_TYPE["SPACE_COOLING"], h in HOURS, td in TYPICAL_DAYS}:
	F [y,c,j] >= peak_sc_factor[c] * F_t [y,c, j, h, td] ;
	
## Additional operational constraints
subject to constantNuc {y in YEARS_WND diff YEAR_ONE, c in REGIONS, h in HOURS, td in TYPICAL_DAYS}:
	F_t [y,c,"NUCLEAR", h, td] = Power_nuclear[y,c];
	
subject to extra_efficiency{y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
	F [y,c,"EFFICIENCY"] = 1 / (1 + i_rate);

## Policy constraints
subject to Minimum_GWP_reduction {y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
	sum{r in RESOURCES} (CO2_net [y,c,r]) <= gwp_limit[y,c];
	
subject to Minimum_GWP_reduction_global {y in YEARS_WND diff YEAR_ONE}:
	sum{c in REGIONS, r in RESOURCES} (CO2_net [y,c,r]) <= gwp_limit_overall [y];

subject to Minimum_RE_share {y in YEARS_WND diff YEAR_ONE, c in REGIONS} :
	sum {j in RE_RESOURCES, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} 
		(R_t_local [y,c, j, h, td]+ R_t_exterior [y,c, j, h, td] + R_t_import [y,c, j, h, td]) * t_op [h, td] 
	>=	re_share_primary[y,c] *
		sum {j in RESOURCES, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} 
		(R_t_local [y,c, j, h, td]+ R_t_exterior [y,c, j, h, td] + R_t_import [y,c, j, h, td]) * t_op [h, td];

subject to f_max_perc_train_pub {y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c,"TRAIN_PUB",h,td] * t_op[h,td]) 
	<= fmax_perc [y,c,"TRAIN_PUB"] * sum {j in TECHNOLOGIES_OF_END_USES_TYPE["MOB_PUBLIC"], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c, j, h, td] * t_op[h,td]);
	
subject to f_max_perc_tramway {y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c,"TRAMWAY_TROLLEY",h,td] * t_op[h,td]) 
	<= fmax_perc [y,c,"TRAMWAY_TROLLEY"] * sum {j in TECHNOLOGIES_OF_END_USES_TYPE["MOB_PUBLIC"], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c, j, h, td] * t_op[h,td]);

subject to f_max_perc_ind_direct_elec {y in YEARS_WND diff YEAR_ONE, c in REGIONS}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c,"IND_DIRECT_ELEC",h,td] * t_op[h,td]) 
	<= fmax_perc [y,c,"IND_DIRECT_ELEC"] * sum {j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_HIGH_T"], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c, j, h, td] * t_op[h,td]);

subject to f_max_perc {y in YEARS_WND diff YEAR_ONE, c in REGIONS, eut in END_USES_TYPES, j in TECHNOLOGIES_OF_END_USES_TYPE[eut]}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c,j,h,td] * t_op[h,td])
	<= fmax_perc [y,c,j] *(sum {j2 in TECHNOLOGIES_OF_END_USES_TYPE[eut], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c, j2, h, td] * t_op[h,td])
	+ sum {r in RESOURCES, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (layers_in_out [y, r, eut] * (R_t_local [y,c, r, h, td] + R_t_exterior [y,c, r, h, td] + R_t_import [y,c, r, h, td] - R_t_export [y,c, r, h, td])));

subject to f_min_perc {y in YEARS_WND diff YEAR_ONE, c in REGIONS, eut in END_USES_TYPES, j in TECHNOLOGIES_OF_END_USES_TYPE[eut]}:
	sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c,j,h,td] * t_op[h,td]) >= fmin_perc [y,c,j] * 
	(sum {j2 in TECHNOLOGIES_OF_END_USES_TYPE[eut], t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (F_t [y,c, j2, h, td] * t_op[h,td])
	+ sum {r in RESOURCES, t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} (layers_in_out [y, r, eut] * (R_t_local [y,c, r, h, td] + R_t_exterior [y,c, r, h, td] + R_t_import [y,c, r, h, td] - R_t_export [y,c, r, h, td])));
	
subject to max_elec_import {y in YEARS_WND diff YEAR_ONE, c in REGIONS, h in HOURS, td in TYPICAL_DAYS}:
	R_t_exterior [y,c,"ELECTRICITY", h, td] * t_op [h, td] <= import_capacity[y,c]; 

#########################################
###   PATHWAY CONSTRAINTS             ###
#########################################

## Phase transitions
subject to phase_new_build {p in PHASE_WND, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS, i in TECHNOLOGIES}:
	F[y_stop,c,i] = F[y_start,c,i] + F_new[p,c,i] - F_old[p,c,i] 
		- sum {p2 in {PHASE_WND union PHASE_UP_TO union {"2015_2021"}}} F_decom[p,p2,c,i];

subject to define_f_decom_properly {p_decom in PHASE, p_built in PHASE union {"2015_2021"}, c in REGIONS, i in TECHNOLOGIES}:
	if decom_allowed[p_decom,p_built,i] == 0 then F_decom[p_decom,p_built,c,i] = 0;

subject to F_new_initialisation {c in REGIONS, tech in TECHNOLOGIES}:
	F_new["2015_2021",c,tech] = F["YEAR_2021",c,tech];

subject to no_decom_if_no_built {c in REGIONS, i in TECHNOLOGIES, p in PHASE_WND union PHASE_UP_TO union {"2015_2021"}}:
	F_new[p,c,i] - sum {p2 in PHASE_WND union PHASE_UP_TO} F_decom[p2,p,c,i] >= 0;

## Technology change limits
var F_used_year_start{YEARS_WND, REGIONS, TECHNOLOGIES} >= 0;
subject to compute_F_used_year_start{p in PHASE_WND, y_start in PHASE_START[p] diff YEAR_ONE, c in REGIONS, j in TECHNOLOGIES} :
	F_used_year_start[y_start,c,j] = (sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} F_t[y_start,c,j,h,td] * t_op[h,td]);

subject to delta_change_definition {p in PHASE_WND, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS, j in TECHNOLOGIES} :
	Delta_change[p,c,j] >= F_used_year_start[y_start,c,j] - (sum {t in PERIODS, h in HOUR_OF_PERIOD[t], td in TYPICAL_DAY_OF_PERIOD[t]} F_t[y_stop,c,j,h,td] * t_op[h,td]);

subject to limit_changes_heat {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["HEAT_LOW_T"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_LT_renovation * (end_uses_input[y_start,c,"HEAT_LOW_T_HW"] + end_uses_input[y_start,c,"HEAT_LOW_T_SH"]);

subject to limit_changes_mob {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["MOBILITY_PASSENGER"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_pass_mob_changes * (end_uses_input[y_start,c,"MOBILITY_PASSENGER"]);

subject to limit_changes_freight {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["MOBILITY_FREIGHT"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_freight_changes * (end_uses_input[y_start,c,"MOBILITY_FREIGHT"] + exch_freight[y_start,c]);

subject to limit_changes_high_heat {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["HEAT_HIGH_T"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_HT_renovation * (end_uses_input[y_start,c,"HEAT_HIGH_T"]);

subject to limit_changes_cooking {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["COOKING"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_cooking_changes * (end_uses_input[y_start,c,"COOKING"]);

subject to limit_changes_mech_comm {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["MECHANICAL_ENERGY_COMM"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_mech_comm_changes * (end_uses_input[y_start,c,"MECHANICAL_ENERGY_COMM"]);

subject to limit_changes_mech_mov_agr {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["MECHANICAL_ENERGY_MOV_AGR"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_mech_mov_agr_changes * (end_uses_input[y_start,c,"MECHANICAL_ENERGY_MOV_AGR"]);

subject to limit_changes_mech_fix_agr {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["MECHANICAL_ENERGY_FIX_AGR"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_mech_fix_agr_changes * (end_uses_input[y_start,c,"MECHANICAL_ENERGY_FIX_AGR"]);

subject to limit_changes_mech_min {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["MECHANICAL_ENERGY_MIN"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_mech_min_changes * (end_uses_input[y_start,c,"MECHANICAL_ENERGY_MIN"]);

subject to limit_changes_mech_fish {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} :
	sum {euc in END_USES_TYPES_OF_CATEGORY["MECHANICAL_ENERGY_FISH_OTHERS"], j in TECHNOLOGIES_OF_END_USES_TYPE[euc]} Delta_change[p,c,j] 
		<= limit_mech_fish_changes * (end_uses_input[y_start,c,"MECHANICAL_ENERGY_FISH_OTHERS"]);

## Anterior phase calculations
subject to Calculate_F_ant_old_based_on_lifespan {p in PHASE, c in REGIONS, t in ANT_TECHNOLOGIES}:
	F_ant_old[p,c,t] = sum{inst in ANT_INSTALL_PERIODS: (phase_to_year[p] > (ant_install_year[inst] + ant_technology_lifespan[t,inst])) and (phase_to_year[p] <= (ant_install_year[inst] + ant_technology_lifespan[t,inst] + t_phase[p]))} ant_install_capacity[t,inst];

subject to Calculate_Total_F_ant_old {c in REGIONS, t in ANT_TECHNOLOGIES}:
	Total_F_ant_old[c,t] = sum {inst in ANT_INSTALL_PERIODS} ant_install_capacity[t,inst];

subject to Force_F_ant_old_to_Zero {p in PHASE, c in REGIONS, i in TECHNOLOGIES diff ANT_TECHNOLOGIES}:
	F_ant_old[p,c,i] = 0;

subject to Force_Total_F_ant_old_to_Zero {c in REGIONS, i in TECHNOLOGIES diff ANT_TECHNOLOGIES}:
	Total_F_ant_old[c,i] = 0;

subject to phase_out_assignment {c in REGIONS, i in TECHNOLOGIES, p in PHASE_WND, age in AGE[i,p]}:
	F_old[p,c,i] = F_ant_old[p,c,i] + (if age == "STILL_IN_USE" then 0
		else (F_new[age,c,i] - sum {p2 in PHASE_WND union PHASE_UP_TO} F_decom[p2,age,c,i] - (if age == "2015_2021" then Total_F_ant_old[c,i] else 0)));

## Cost calculations
subject to total_capex:
	C_tot_capex = sum{p in PHASE_WND union PHASE_UP_TO union {"2015_2021"}} C_inv_phase[p]
		- sum {c in REGIONS, i in TECHNOLOGIES} C_inv_return[c,i];

subject to investment_computation {p in PHASE_WND union PHASE_UP_TO union {"2015_2021"}, y_start in PHASE_START[p], y_stop in PHASE_STOP[p]}:
	C_inv_phase[p] = sum {c in REGIONS, i in TECHNOLOGIES} F_new[p,c,i] * annualised_factor[p] * (c_inv[y_start,c,i] + c_inv[y_stop,c,i]) / 2;

subject to investment_computation_tech {p in PHASE_WND union PHASE_UP_TO union {"2015_2021"}, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS, i in TECHNOLOGIES}:
	C_inv_phase_tech[p,c,i] = F_new[p,c,i] * annualised_factor[p] * (c_inv[y_start,c,i] + c_inv[y_stop,c,i]) / 2;

subject to investment_return {c in REGIONS, i in TECHNOLOGIES}:
	C_inv_return[c,i] = sum {p in PHASE_WND union PHASE_UP_TO union {"2015_2021"}, y_start in PHASE_START[p], y_stop in PHASE_STOP[p]} 
		(remaining_years[i,p] / lifetime[y_start,c,i] * (F_new[p,c,i] - sum {p2 in PHASE_WND union PHASE_UP_TO} F_decom[p2,p,c,i]) * annualised_factor[p] * (c_inv[y_start,c,i] + c_inv[y_stop,c,i]) / 2);

subject to Opex_tot_cost_calculation:
	C_tot_opex = C_opex["YEAR_2021"] 
		+ sum {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p]} 
			(t_phase[p] * (C_opex[y_start] + C_opex[y_stop])/2 * annualised_factor[p]);

subject to Opex_cost_calculation{y in YEARS_WND union YEARS_UP_TO}:
	C_opex[y] = sum {c in REGIONS, j in TECHNOLOGIES} C_maint[y,c,j] + sum {c in REGIONS, i in RESOURCES} C_op[y,c,i];

subject to operation_computation_tech {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS, i in TECHNOLOGIES}:
	C_op_phase_tech[p,c,i] = t_phase[p] * ((C_maint[y_start,c,i] + C_maint[y_stop,c,i])/2 * annualised_factor[p]);

subject to operation_computation_res {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS, i in RESOURCES}:
	C_op_phase_res[p,c,i] = t_phase[p] * ((C_op[y_start,c,i] + C_op[y_stop,c,i])/2 * annualised_factor[p]);

subject to maxInvestment {p in PHASE_WND}:
	C_inv_phase[p] <= max_inv_phase[p];

subject to totalGWPTransition_calculation:
	TotalGWPTransition = sum{c in REGIONS} TotalGWP["YEAR_2021",c] + sum {p in PHASE_UP_TO union PHASE_WND, y_start in PHASE_START[p], y_stop in PHASE_STOP[p]} 
		(t_phase[p] * (sum{c in REGIONS} TotalGWP[y_start,c] + sum{c in REGIONS} TotalGWP[y_stop,c])/2);

subject to minimum_GWP_transition:
	TotalGWPTransition <= gwp_limit_transition;

##########################
### OBJECTIVE FUNCTION ###
##########################

# Option 1: Minimize total transition cost (CAPEX + OPEX across all regions and phases)
minimize TotalTransitionCost: C_tot_capex + C_tot_opex;