import toml

def open_toml(path_to_tomlFile):
    with open(path_to_tomlFile) as f:
        config = toml.load(f)
    
    return config
      
# cutoff_pp2_upper = config['Prediction_tools']['pp2hvar']['damaging']
# cutoff_pp2_middle = config['Prediction_tools']['pp2hvar']['possibly']