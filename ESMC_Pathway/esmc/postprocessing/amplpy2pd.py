"""Helpers to convert amplpy objects to pandas and JSON.

ESMC Pathway - Pablo Jimenez Zabalaga
Based on EnergyScope Multi-Cell (original author: Paolo Thiran).
"""
import json

def print_json(my_sets, file):  # printing the dictionary containing all the sets into directory/sets.json
        with open(file, 'w') as fp:
            json.dump(my_sets, fp, indent=4, sort_keys=True)
        return

def read_json(file):
        # reading the saved dictionary containing all the sets from directory/sets.json
        with open(file, 'r') as fp:
            data = json.load(fp)
        return data

