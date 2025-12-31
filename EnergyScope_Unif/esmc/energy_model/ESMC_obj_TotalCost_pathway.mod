##########################
### OBJECTIVE FUNCTION ###
##########################

# Multi-regional Energy System with Pathway Optimization
# The objective function minimizes the total transition cost across all regions and phases

# Option 1: Minimize total transition cost (CAPEX + OPEX across all regions and phases)
#minimize TotalTransitionCost: C_tot_capex + C_tot_opex;





# Option 2: Minimize total transition cost with additional freight costs (if freight optimization is important)
# minimize TotalTransitionCostWithFreight: 
#     C_tot_capex + C_tot_opex + 
#     sum{y in YEARS_WND, c in REGIONS} (260 * Exch_freight[y,c] * annualised_factor[p]);

# Option 3: Minimize with GWP penalty (if carbon pricing is considered)
# param gwp_cost {YEARS} >= 0 default 0; # Cost per ton of CO2 [€/tCO2]
# var Gwp_tot_cost >= 0;
# 
# subject to Gwp_tot_cost_calculation:
#     Gwp_tot_cost = sum {p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p]} 
#         (t_phase[p] * (sum{c in REGIONS} (TotalGWP[y_start,c] * gwp_cost[y_start] + TotalGWP[y_stop,c] * gwp_cost[y_stop]))/2);
#
# minimize TotalTransitionCostWithGWP: C_tot_capex + C_tot_opex + Gwp_tot_cost;

# Option 4: Multi-objective with weighted freight transport costs by region
# param freight_cost_weight {REGIONS} >= 0 default 260; # Different freight costs per region
# minimize TotalTransitionCostWeighted: 
#     C_tot_capex + C_tot_opex + 
#     sum{p in PHASE_WND union PHASE_UP_TO, y_start in PHASE_START[p], y_stop in PHASE_STOP[p], c in REGIONS} 
#         (t_phase[p] * freight_cost_weight[c] * (Exch_freight[y_start,c] + Exch_freight[y_stop,c])/2 * annualised_factor[p]);

# Option 5: Minimize only GWP emissions across transition (for environmental scenarios)
# minimize TotalTransitionGWP: TotalGWPTransition;

# Additional options for debugging or specific analyses:

# Option 6: Minimize cost for a specific region (useful for regional studies)
# param focus_region symbolic in REGIONS;
# minimize RegionalCost: 
#     sum{p in PHASE_WND union PHASE_UP_TO union {"2015_2021"}} C_inv_phase_tech[p,focus_region,*] +
#     sum{p in PHASE_WND union PHASE_UP_TO, y in YEARS} C_maint[y,focus_region,*] * annualised_factor[p];

# Option 7: Minimize technology-specific costs (for technology assessment)
# param focus_tech symbolic in TECHNOLOGIES;
# minimize TechCost:
#     sum{p in PHASE_WND union PHASE_UP_TO union {"2015_2021"}, c in REGIONS} C_inv_phase_tech[p,c,focus_tech] +
#     sum{y in YEARS_WND union YEARS_UP_TO, c in REGIONS} C_maint[y,c,focus_tech];

# Note: The default option (Option 1) minimizes the total discounted system cost across all regions and phases,
# which is the most common objective for multi-regional energy system planning with pathway optimization.