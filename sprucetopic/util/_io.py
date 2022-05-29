
def read_config(config_file):
    import yaml
    with open(config_file) as f:
        params = yaml.safe_load(f)
    return params
