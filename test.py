import configparser
config=configparser.ConfigParser()
config.read('coordinates.txt','ISO-8859-1')
Default_config=config['Default']
print(Default_config)
