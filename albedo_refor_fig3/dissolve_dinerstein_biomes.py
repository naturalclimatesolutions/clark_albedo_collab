import geopandas as gpd
biomes = gpd.read_file('./Ecoregions2017.shp').loc[:, ['BIOME_NUM',
                                                       'BIOME_NAME',
                                                       'geometry',
                                                      ]].dissolve(by='BIOME_NUM')
biomes = biomes.reset_index()
biomes.columns = ['BIOME_NUM', 'BIOME_NAME', 'geometry']
biomes.to_file('./Ecoregions2017_BIOME_DISSOLVE.shp', index=False)
