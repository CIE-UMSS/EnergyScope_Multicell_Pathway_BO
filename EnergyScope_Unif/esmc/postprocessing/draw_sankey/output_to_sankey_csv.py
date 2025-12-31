"""This script post-processes the results to get data for the sankey per region

Author: Pablo Jimenez
"""

import pandas as pd
import numpy as np
from pathlib import Path


class Cell:

    def __init__(self, name, year_balance, storage):
        self.name = name
        self.year_balance = year_balance
        self.storage = storage

        self.year_balance = self.arrange_columns(self.year_balance.copy(), RegroupLayers)
        self.year_balance = self.arrange_rows(self.year_balance.copy(), RegroupElements)
        self.storage = self.arrange_columns(self.storage.copy(), StorageField)
        self.storage = self.arrange_rows(self.storage.copy(), EndUseStorage)

        self.update_year_balance()

    def arrange(self, dataframe, regroup_dict):

        for column_name, group_columns in regroup_dict.items():
            if len(group_columns) > 1:
                new_column = dataframe.pop(group_columns[0])
                for i in range(1, len(group_columns)):
                    new_column = new_column.add(dataframe.pop(group_columns[i]))
                new_column.name = column_name
                dataframe = pd.concat([dataframe, new_column], axis=1)

            else:
                dataframe = dataframe.rename(columns={group_columns[0]: column_name})
        for column in dataframe.columns:
            if column not in regroup_dict.keys():
                dataframe.pop(column)
        return dataframe

    def arrange_columns(self, data_frame, regroup_dict):
        return self.arrange(data_frame, regroup_dict)

    def arrange_rows(self, data_frame, regroup_dict):
        return self.arrange(data_frame.T, regroup_dict).T

    def update_year_balance(self):
        for sto in StorageLayer:
            self.year_balance.loc[sto, StorageLayer[sto]] = -self.storage.loc[sto, "Year energy flux"]


def read_regional_electricity_demand(proj_dir, space_id, case_study):
    """
    Lee las demandas eléctricas regionales desde el archivo reg_demands.dat
    
    Returns:
        dict: Diccionario con demandas por región y total
    """
    # Ruta al archivo de demandas regionales
    reg_demands_file = proj_dir / "case_studies" / space_id / case_study / "reg_demands.dat"
    
    # Diccionario para almacenar las demandas por región
    regional_demands = {}
    total_system_demand = 0
    
    try:
        with open(reg_demands_file, "r") as file:
            lines = file.readlines()
            
            # Buscar líneas que contengan demanda de ELECTRICITY
            for line in lines:
                # Saltar líneas de comentarios o encabezados
                if line.strip().startswith('param') or line.strip() == '' or line.strip().endswith(':='):
                    continue
                if line.strip() == ';':  # Final del archivo
                    break
                    
                # Dividir la línea por tabulaciones o espacios múltiples
                parts = line.strip().split()
                
                if len(parts) >= 3:  # Debe tener al menos región, tipo y un valor
                    region = parts[0]
                    demand_type = parts[1]
                    
                    # Solo procesar líneas de demanda eléctrica
                    if demand_type == "ELECTRICITY":
                        # Sumar todos los valores de los sectores (columnas 3 en adelante)
                        sector_values = [float(val) for val in parts[2:] if val != '0.0']
                        regional_demand = sum(float(val) for val in parts[2:])
                        
                        # Almacenar la demanda regional
                        regional_demands[region] = regional_demand
                        total_system_demand += regional_demand
                        
                        print(f"Región {region}: Demanda eléctrica = {regional_demand:.3f} TWh")
    
    except FileNotFoundError:
        print(f"Advertencia: No se encontró el archivo {reg_demands_file}")
        print("Se usará el método anterior de cálculo de demanda")
        return None
    
    # Agregar el total del sistema
    regional_demands["Total"] = total_system_demand
    print(f"Total del sistema: Demanda eléctrica = {total_system_demand:.3f} TWh")
    
    return regional_demands


def write_sankey_file(space_id, case_study):
    proj_dir = Path(__file__).parents[3]
    output_dir = proj_dir / "case_studies" / space_id / case_study / "outputs/regional_results/"

    # Leer las demandas eléctricas reales desde reg_demands.dat
    print("Leyendo demandas eléctricas regionales...")
    regional_elec_demands = read_regional_electricity_demand(proj_dir, space_id, case_study)

    with open(output_dir / "Year_balance.csv", "r") as year_balance_file:

        # Get the name of all the macrocells
        year_balance_file.readline()
        first_column = [line.split(";")[0] for line in year_balance_file.readlines()]
        cells_name = []
        for index_name in first_column:
            if index_name not in cells_name:
                cells_name.append(index_name)

        year_balance_file.seek(0)

        # Get year_balance file under the form of a panda dataframe
        all_data_balance = pd.read_csv(year_balance_file, index_col=[0, 1], sep=";")
        all_data_balance = all_data_balance.replace(to_replace=np.nan, value=0)

    with open(output_dir / "Sto_assets.csv", "r") as sto_assets_file:

        # Get Sto_assets file under the form of a panda dataframe
        all_data_sto = pd.read_csv(sto_assets_file, index_col=[0, 1], sep=";")
        all_data_sto = all_data_sto.replace(to_replace=np.nan, value=0)
    
    # Read Curt.csv into a pandas dataframe and calculate the total sum of the Curt column for each unique value of the first column
    curt_df = pd.read_csv(output_dir / "Curt.csv", delimiter=';')

    curt_sum_by_region = curt_df.groupby(curt_df.columns[0])["Curt"].sum()
    total_curt=curt_df["Curt"].sum()

    cells = {}
    for cell_name in cells_name:
        cells[cell_name] = Cell(cell_name, all_data_balance.loc[cell_name], all_data_sto.loc[cell_name])

    total_data_balance = all_data_balance.loc[cells_name[0]]
    total_data_sto = all_data_sto.loc[cells_name[0]]
    for i in range(1, len(cells_name)):
        total_data_balance = total_data_balance.add(all_data_balance.loc[cells_name[i]])
        total_data_sto = total_data_sto.add(all_data_sto.loc[cells_name[i]])

    cells["Total"] = Cell("Total", total_data_balance, total_data_sto)

    for cell in cells.values():
        file_name = "input2sankey_" + cell.name + ".csv"
        total_elec_demand = 0  # Generación eléctrica total de la región
        
        with open(output_dir / file_name, "w") as input2csv_file:
            print("source,target,realValue,layerID,layerColor,layerUnit", file=input2csv_file)

            for tech in cell.year_balance.index:
                tech_count = 0
                for layer in cell.year_balance.columns:
                    value = cell.year_balance.loc[tech, layer]
                    if value > 5 and tech != layer:
                        if tech == "End Use":
                            continue
                        else:
                            if tech not in TechLayer.keys():
                                print(
                                    "%s,%s,%f,%s,%s,%s" % (tech, layer, value / 1000, layer, LayerColor[layer], "TWh"),
                                    file=input2csv_file)
                            else:
                                print("%s,%s,%f,%s,%s,%s" %
                                      (tech, layer, value / 1000, TechLayer[tech], LayerColor[TechLayer[tech]], "TWh"),
                                      file=input2csv_file)

                    if value > 0:
                        tech_count = tech_count + 1

            for layer in cell.year_balance.T.index:
                for tech in cell.year_balance.T.columns:
                    value = cell.year_balance.T.loc[layer][tech]
                    if value < -5:
                        if layer in EndUseLayer:
                            continue
                        elif tech != "End Use":
                            print("%s,%s,%f,%s,%s,%s" % (layer, tech, -value / 1000, layer, LayerColor[layer], "TWh"),
                                  file=input2csv_file)

            for layer in EndUseName:
                value = cell.year_balance.loc["End Use"][layer]
                if value < -5:
                    print("%s,%s,%f,%s,%s,%s" %
                          (layer, EndUseName[layer], -value / 1000, layer, LayerColor[layer], "TWh"),
                          file=input2csv_file)

        # Leer el archivo generado para calcular la generación eléctrica total
        with open(output_dir / file_name, "r") as previous_file:
            for line in previous_file:
                parts = line.strip().split(',')
                if parts[0] == "Elec" and parts[1] == "Elec Demand":
                    total_elec_demand += float(parts[2]) * 1000  # Convertir de TWh a GWh

        # Obtener la demanda eléctrica real de la región
        if regional_elec_demands and cell.name in regional_elec_demands:
            # Usar la demanda real del archivo reg_demands.dat
            real_demand = regional_elec_demands[cell.name]  # Convertir TWh a GWh
        else:
            # Fallback: usar el método anterior si no se puede leer el archivo
            print(f"Advertencia: Usando método anterior para {cell.name}")
            if cell.name in curt_sum_by_region.index:
                real_demand = total_elec_demand - curt_sum_by_region[cell.name]
            elif cell.name == "Total":
                real_demand = total_elec_demand - total_curt
            else:
                real_demand = total_elec_demand

        # Calcular las pérdidas eléctricas
        elec_losses_value = total_elec_demand - real_demand
        
        print(f"\n--- Región {cell.name} ---")
        print(f"Generación eléctrica total: {total_elec_demand:.1f} GWh")
        print(f"Demanda eléctrica real: {real_demand:.1f} GWh")
        print(f"Pérdidas eléctricas: {elec_losses_value:.1f} GWh")

        # Escribir las entradas adicionales al archivo CSV
        with open(output_dir / file_name, "a") as input2csv_file:
            color = LayerColor["Elec"]
            
            # Escribir curtailment
            if cell.name in curt_sum_by_region.index:
                curt_value = curt_sum_by_region[cell.name]
                print("Elec,Curt.,%f,Curtailment,%s,TWh" % (curt_value / 1000, color), file=input2csv_file)
            elif cell.name == "Total":
                print("Elec,Curt.,%f,Curtailment,%s,TWh" % (total_curt / 1000, color), file=input2csv_file)
            
            # Escribir pérdidas eléctricas
            print("Elec,Elec Losses,%f,Elec_Losses,%s,TWh" % (elec_losses_value / 1000, color), file=input2csv_file)

        # Leer el archivo para modificar la línea de Elec Demand
        file_content = []
        with open(output_dir / file_name, "r") as input2csv_file:
            for line in input2csv_file:
                parts = line.strip().split(',')
                if parts[0] == "Elec" and parts[1] == "Elec Demand":
                    # Usar la demanda real calculada
                    line = "Elec,Elec Demand,%f,Elec Demand,%s,TWh\n" % (real_demand / 1000, color)
                file_content.append(line)

        # Escribir el contenido modificado de vuelta al archivo
        with open(output_dir / file_name, "w") as input2csv_file:
            input2csv_file.writelines(file_content)

RegroupElements = {
    "Nuclear": ["NUCLEAR"],
    "Gas Turbines": ["CCGT", "OCGT"],
    "CCGT_Ammonia": ["CCGT_AMMONIA"],
    "Coal_US": ["COAL_US"],
    "Coal_IGCC": ["COAL_IGCC"],
    "Solar PV": ["PV_ROOFTOP", "PV_UTILITY"],
    "CSP Power Block": ["PT_POWER_BLOCK", "ST_POWER_BLOCK"],
    "Wind Onshore": ["WIND_ONSHORE"],
    "Wind Offshore": ["WIND_OFFSHORE"],
    "Hydro": ["HYDRO_DAM", "HYDRO_RIVER"],
    # "Hydro Dam":        ["HYDRO_DAM"],
    # "Hydro River":      ["HYDRO_RIVER"],
    "Tidal Power": ["TIDAL_STREAM", "TIDAL_RANGE"],
    "Wave": ["WAVE"],
    "Geothermal": ["GEOTHERMAL"],
    "Diesel Genset": ["GENSET_DIESEL"],
    "Biomass power plant": ["ST_BIOMASS", "ST_SNG", "FB_ST_BIOMASS", "CFB_ST_BIOMASS", "BFB_ST_BIOMASS"],
    "Fuel cell": ["FUEL_CELL"],
    "Biofuels": ["FERMENTATION_TO_BIOETHANOL", "ESTERIFICATION_TO_BIODIESEL"],
    "Ind Cogen": ["IND_COGEN_GAS", "IND_COGEN_WOOD", "IND_COGEN_WASTE"],
    "Ind Boiler": ["IND_BOILER_GAS", "IND_BOILER_WOOD", "IND_BOILER_BIOWASTE",
                   "IND_BOILER_OIL", "IND_BOILER_COAL", "IND_BOILER_WASTE", "IND_BOILER_DIESEL", "IND_BOILER_LPG"],
    "Ind Direct Elec": ["IND_DIRECT_ELEC"],
    "Mechanical Energy": ["COMM_MACHINERY_DIESEL", "COMM_MACHINERY_EL", "IND_MACHINERY_EL", "TRACTOR_DIESEL", "TRACTOR_EL", "AGR_MACHINERY_DIESEL", "AGR_MACHINERY_EL", "MIN_MACHINERY_DIESEL", "MIN_MACHINERY_EL", "FISH_MACHINERY_DIESEL", "FISH_MACHINERY_EL"],
    "HPs": ["DHN_HP_ELEC", "DEC_HP_ELEC"],
    "DHN Tech": ["DHN_COGEN_GAS", "DHN_COGEN_WOOD", "DHN_COGEN_WASTE",
                 #"DHN_COGEN_WET_BIOMASS", "DHN_COGEN_BIO_HYDROLYSIS",
                 "DHN_BOILER_GAS", "DHN_BOILER_WOOD", "DHN_BOILER_OIL",
                 "DHN_DEEP_GEO", "DHN_SOLAR"],
    "DEC Heat": ["DEC_THHP_GAS", "DEC_COGEN_GAS", "DEC_COGEN_OIL", "DEC_ADVCOGEN_GAS",
                 "DEC_ADVCOGEN_H2", "DEC_BOILER_GAS", "DEC_BOILER_WOOD", "DEC_BOILER_OIL",
                 "DEC_SOLAR", "DEC_DIRECT_ELEC"],
    "Stoves": ["STOVE_WOOD", "STOVE_LPG", "STOVE_NG", "STOVE_OIL", "STOVE_ELEC"],
    "Lighting": ["CONVENTIONAL_BULB", "LED_BULB", "CONVENTIONAL_LIGHT", "LED_LIGHT"],
    "Food preservation": ["REFRIGERATOR_EL"],
    "Cooling": ["DEC_ELEC_COLD", "IND_ELEC_COLD"
                     #"DEC_THHp_GAS_COLD", 
                     ],
    # "Big Split":        ["BIG_SPLIT"],
    # "Chiller":          ["CHILLER_WC"],

    #"Mobility": ["TRAMWAY_TROLLEY", "BUS_COACH_DIESEL", "BUS_COACH_HYDIESEL",
                 #"BUS_COACH_CNG_STOICH", "BUS_COACH_FC_HYBRIDH2", "TRAIN_PUB",
                 #"PLANE_JETFUEL",
                 #"CAR_GASOLINE", "CAR_DIESEL", "CAR_NG",
                 #"CAR_METHANOL", "CAR_HEV", "CAR_PHEV", "CAR_BEV", "CAR_FUEL_CELL"],
    "Public Mob": ["TRAMWAY_TROLLEY", "BUS_COACH_DIESEL", "BUS_COACH_HYDIESEL", "BUS_COACH_CNG_STOICH", "BUS_COACH_FC_HYBRIDH2", "TRAIN_PUB", "CAR_FG_PUBLIC", "BUS_FG_PUBLIC", "PICKUP_TRUCK_FG_PUBLIC", "SUV_FG_PUBLIC", "MOTORCYCLE_FG_PUBLIC", "CAR_GASOLINE_PUBLIC", "BUS_GASOLINE_PUBLIC", "PICKUP_TRUCK_GASOLINE_PUBLIC", "SUV_GASOLINE_PUBLIC", "MOTORCYCLE_GASOLINE_PUBLIC", "CAR_DIESEL_PUBLIC", "BUS_DIESEL_PUBLIC", "PICKUP_TRUCK_DIESEL_PUBLIC", "SUV_DIESEL_PUBLIC", "MOTORCYCLE_DIESEL_PUBLIC", "BUS_ELEC_PUBLIC", "PICKUP_TRUCK_ELEC_PUBLIC", "SUV_ELEC_PUBLIC", "MOTORCYCLE_ELEC_PUBLIC"],
    "Private Mob": ["CAR_GASOLINE", "CAR_DIESEL", "CAR_NG", "CAR_METHANOL", "CAR_HEV", "CAR_PHEV", "CAR_BEV", "CAR_FUEL_CELL", "CAR_FG_PRIVATE", "BUS_FG_PRIVATE", "PICKUP_TRUCK_FG_PRIVATE", "SUV_FG_PRIVATE", "MOTORCYCLE_FG_PRIVATE", "CAR_GASOLINE_PRIVATE", "BUS_GASOLINE_PRIVATE", "PICKUP_TRUCK_GASOLINE_PRIVATE", "SUV_GASOLINE_PRIVATE", "MOTORCYCLE_GASOLINE_PRIVATE", "CAR_DIESEL_PRIVATE", "BUS_DIESEL_PRIVATE", "PICKUP_TRUCK_DIESEL_PRIVATE", "SUV_DIESEL_PRIVATE", "MOTORCYCLE_DIESEL_PRIVATE", "BUS_ELEC_PRIVATE", "PICKUP_TRUCK_ELEC_PRIVATE", "SUV_ELEC_PRIVATE", "MOTORCYCLE_ELEC_PRIVATE"],

    "Freight": ["TRAIN_FREIGHT", "TRAIN_FREIGHT_ELEC", "BOAT_FREIGHT_DIESEL", "BOAT_FREIGHT_NG",
                "BOAT_FREIGHT_METHANOL", "BOAT_FREIGHT_ELEC", "TRUCK_DIESEL", "TRUCK_FUEL_CELL",
                "TRUCK_NG", "TRUCK_METHANOL", "TRUCK_ELEC", "TRUCK_FG", "VAN_FG", "CARGO_MOTORCYCLE_FG", "TRUCK_GASOLINE", "VAN_GASOLINE", "CARGO_MOTORCYCLE_GASOLINE", "TRUCK_DIESEL_P", "VAN_DIESEL", "CARGO_MOTORCYCLE_DIESEL", "VAN_ELEC", "CARGO_MOTORCYCLE_ELEC", "PLANE"],
    # "Rail Freight":    ["TRAIN_FREIGHT"],
    # "Boat Freight":     ["BOAT_FREIGHT_DIESEL", "BOAT_FREIGHT_NG", "BOAT_FREIGHT_METHANOL"],
    # "Road Freight":     ["TRUCK_DIESEL", "TRUCK_METHANOL", "TRUCK_FUEL_CELL",
    #                     "TRUCK_ELEC", "TRUCK_NG"],
    # "Int. Freight": ["CONTAINER_CARGO_DIESEL", "CONTAINER_CARGO_LNG",
    #                  "CONTAINER_CARGO_METHANOL", "CONTAINER_CARGO_AMMONIA",
    #                  "CONTAINER_CARGO_FUELCELL_AMMONIA", "CONTAINER_CARGO_RETRO_METHANOL",
    #                  "CONTAINER_CARGO_RETRO_AMMONIA", "CONTAINER_CARGO_FUELCELL_LH2"],
    "CSP collector": ["PT_COLLECTOR", "ST_COLLECTOR"],
    "H2.": ["H2_ELECTROLYSIS", "H2_NG", "H2_BIOMASS", "AMMONIA_TO_H2"],
    "Gasifi SNG": ["GASIFICATION_SNG"],
    "To Methane": ["SYN_METHANATION", "BIOMETHANATION_WET_BIOMASS", "BIOMETHANATION_BIOWASTE"],
    "Pyrolise": ["PYROLYSIS_BIOWASTE_TO_FUELS"
                 #"PYROLYSIS_LIGNO_TO_LFO", "PYROLYSIS_LIGNO_TO_FUELS", "PYROLYSIS_BIOWASTE_TO_LFO"
                 ],
    "To Methanol": ["SYN_METHANOLATION", "METHANE_TO_METHANOL", "BIOMASS_TO_METHANOL"],
    "Haber Bosch": ["HABER_BOSCH"],
    #"Fischer-Tropsch": ["FISCHER_TROPSCH_DIESEL", "FISCHER_TROPSCH_GASOLINE", "FISCHER_TROPSCH_JETFUEL"],
    "Non-Energy Demand": ["OIL_TO_HVC", "GAS_TO_HVC", "BIOMASS_TO_HVC", "METHANOL_TO_HVC"],
    #"Elec in/out": ["ELECTRICITY"],
    "Prod. & Imp. Gasoline": ["GASOLINE"],
    #"Oil RE imports": ["GASOLINE_RE"],
    #"Jet Fuel imports": ["JET_FUEL"],
    #"Jet Fuel RE imports": ["JET_FUEL_RE"],
    "Prod. & Imp. Diesel": ["DIESEL"],
    "Bioethanol imports": ["BIOETHANOL"],
    "Biodiesel imports": ["BIODIESEL"],
    "LFO": ["LFO"],
    "Fossil Gas": ["GAS"],
    "Gas RE imports": ["GAS_RE"],
    "Ligno. biomass": ["WOOD", "ENERGY_CROPS_2"],
    "Biowaste": ["BIOWASTE", "BIOMASS_RESIDUES"],
    "Wet biomass": ["WET_BIOMASS"],
    "Coal imports": ["COAL"],
    "Waste": ["WASTE"],
    "H2 imports": ["H2"],
    "H2 RE imports": ["H2_RE"],
    "Ammonia imports": ["AMMONIA"],
    "Methanol imports": ["METHANOL"],
    "Ammonia RE imports": ["AMMONIA_RE"],
    "Methanol RE imports": ["METHANOL_RE"],
    "DEC Sto.": ["TS_DEC_DIRECT_ELEC", "TS_DEC_HP_ELEC", "TS_DEC_THHP_GAS",
                "TS_DEC_COGEN_GAS", "TS_DEC_COGEN_OIL", "TS_DEC_ADVCOGEN_GAS",
                "TS_DEC_ADVCOGEN_H2", "TS_DEC_BOILER_GAS", "TS_DEC_BOILER_WOOD",
                "TS_DEC_BOILER_OIL"],
    "DHN Sto.": ["TS_DHN_DAILY", "TS_DHN_SEASONAL"],
    "Cold Sto.": ["TS_COLD"],
    "End Use": ["END_USES"]
}

# Names of the layers.
RegroupLayers = {
    "Elec": ["ELECTRICITY"],
    "Gasoline": ["GASOLINE"],
    "Jet Fuel": ["JET_FUEL"],
    "Diesel": ["DIESEL"],
    "LPG": ["LPG"],
    "LFO": ["LFO"],
    "Fossil Gas": ["GAS"],
    "Ligno. biomass": ["WOOD", 
                      #"ENERGY_CROPS_2"
                      ],
    "Biowaste": ["BIOWASTE", 
                #"BIOMASS_RESIDUES"
                ],
    "Wet biomass": ["WET_BIOMASS"],
    "Coal": ["COAL"],
    "Waste": ["WASTE"],
    "H2.": ["H2"],
    "Ammonia": ["AMMONIA"],
    "Methanol": ["METHANOL"],
    "Non-Energy Demand": ["HVC"],
    "Heat HT": ["HEAT_HIGH_T"],
    "Heat LT DHN": ["HEAT_LOW_T_DHN"],
    "Heat LT": ["HEAT_LOW_T_DECEN"],
    "Cooling": ["SPACE_COOLING", "PROCESS_COOLING"],
    # "Space Cool": ["SPACE_COOLING"],
    # "Process Cool": ["PROCESS_COOLING"],
    # "Mobility": ["MOB_PUBLIC", "MOB_PRIVATE"],
    "Public Mob": ["MOB_PUBLIC"],
    "Private Mob": ["MOB_PRIVATE"],
    "Freight": ["MOB_FREIGHT_RAIL", "MOB_FREIGHT_ROAD", "MOB_FREIGHT_BOAT"],
    "Lighting": ["LIGHTING_R_C", "LIGHTING_P"],
    "Cooking": ["COOKING"],
    "Food preservation": ["FOOD_PRESERVATION"],
    "Mechanical Energy": ["MECHANICAL_ENERGY_COMM", "MECHANICAL_ENERGY_IND", "MECHANICAL_ENERGY_MOV_AGR", "MECHANICAL_ENERGY_FIX_AGR", "MECHANICAL_ENERGY_MIN", "MECHANICAL_ENERGY_FISH_OTHERS"],
    # "Rail Freight": ["MOB_FREIGHT_RAIL"],
    # "Road Freight": ["MOB_FREIGHT_ROAD"],
    # "Boat Freight": ["MOB_FREIGHT_BOAT"],
    #"Int. Freight": ["CONTAINER_FREIGHT"],
    "CSP": ["PT_HEAT", "ST_HEAT"],
    "Solar": ["RES_SOLAR"]
}

# Link the tech with the layer name (as in LayerColor) of its output. If several
# output, they must be listed in order of appearance of the layer in the "Layers_in_out.scv"
# file (from left to right).


TechLayer = {
    "Solar PV": "RES Solar",
    "Wind Onshore": "RES Wind",
    "Wind Offshore": "RES Wind",
    "Hydro": "RES Hydro",
    # "Hydro Dam":        "RES Hydro",
    # "Hydro River":      "RES Hydro",
    "Tidal Power": "RES Hydro",
    "Wave": "RES Hydro",
    "Geothermal": "Geothermal",
    "Gas Turbines": "Fossil Gas",
    "Diesel Genset": "Diesel",
    "Biomass power plant": "Biowaste",
    "Nuclear": "Uranium"
}

# It is only relevant to show storage of layer linked to end use only.

EndUseStorage = {
    "DEC Sto.": ["TS_DEC_DIRECT_ELEC", "TS_DEC_HP_ELEC", "TS_DEC_THHP_GAS",
                "TS_DEC_COGEN_GAS", "TS_DEC_COGEN_OIL", "TS_DEC_ADVCOGEN_GAS",
                "TS_DEC_ADVCOGEN_H2", "TS_DEC_BOILER_GAS", "TS_DEC_BOILER_WOOD",
                "TS_DEC_BOILER_OIL"],
    "DHN Sto.": ["TS_DHN_DAILY", "TS_DHN_SEASONAL"],
    "Cold Sto.": ["TS_COLD"]
}

StorageField = {
    "Year energy flux": ["Year_energy_flux"]
}

StorageLayer = {
    "DEC Sto.": ["Heat LT DEC"],
    "DHN Sto.": ["Heat LT DHN"],
    # "Cold Sto.":         ["Space Cool"]
}

# Name of the end use of the end use layer. If the layer is end use only
# (There is no "-" is the layer column of "Layers_in_out" matrix), its
# name is the same as the name of the layer itself.

EndUseLayer = [
    "Heat LT DHN",
    "Heat LT",
    "Cooling",
    # "Space Cool",
    # "Process Cool",
    #"Mobility",
    "Freight",
    "Public Mob",
    "Private Mob",
    "Mech. Energy",
    "Lighting"
    # "Road Freight",
    # "Boat Freight",
    # "Rail Freight",
    #"Int. Freight"
]

EndUseName = {
    "Elec": "Elec Demand",
    "Ammonia": "Non-Energy Demand",
    "Methanol": "Non-Energy Demand",
    #"HVC": "Non-Energy Demand",
    #"Heat HT": "Ind Heat Demand",
}

LayerColor = {
    "Elec": "#00BFFF",
    "Gasoline": "#808080",
    "Jet Fuel": "#deb887",
    "Diesel": "#D3D3D3",
    "LFO": "#8B008B",
    "LPG": "#da70d6",
    "Fossil Gas": "#FFD700",
    "Ligno. biomass": "#CD853F",
    "Biowaste": "#3F5901",
    "Wet biomass": "#336600",
    "Coal": "#A0522D",
    "Uranium": "#FFC0CB",
    "Waste": "#808000",
    "H2.": "#FF00FF",
    "Ammonia": "#000ECD",
    "Methanol": "#CC0066",
    "Non-Energy Demand": "#00FFFF",
    # "HVC": "#00FFFF",
    "Heat HT": "#DC143C",
    "Heat LT DHN": "#FA8072",
    "Heat LT": "#6a5acd",
    "Cooling": "#00CED1",
    "Cooking": "#00CED1",
    "Lighting": "#00CED1",
    "Mech. Energy": "#00CED1",
    # "Space Cool": "#00CED1",
    # "Process Cool": "#00CED1",
    "RES Wind": "#27AE34",
    "RES Solar": "#FF7F50",
    "RES Hydro": "#00CED1",
    "RES Geo": "#FF0000",
    "Solar": "#FF7F50",
    "CSP": "#FF7F50",
    "Geothermal": "#FF0000"
}
