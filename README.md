# meteo_interpolations_wsl
Scripts in order to interpolate meteorological variables and validate them with LWF data

# Interpolation files
interpolate_precipitation_hourly: Interpolate hourly precipitation data of MeteoSwiss stations

interpolate_precipitation_10min: Interpolate 10min precipitation data of MeteoSwiss stations

interpolate_temperature_10min: Interpolate 10min temperature data of MeteoSwiss stations

interpolate_relativehumidity_10min: Interpolate 10min relative humidity data of MeteoSwiss stations

interpolate_radiation_10min: Interpolate 10min radiation data of MeteoSwiss stations

# Validation files
precip_statistics_hourly_anywslsite: Statistics of hourly precipitation interpolations, CombiPrecip and LWF precipitation data

precip_statistics_hourly_vertical: Statistics of hourly precipitation interpolations including a vertical interpolation, CombiPrecip and LWF precipitation data

precip_validation_hourly_anywslsite: Bar plot comparison of hourly precipitation interpolations, CombiPrecip and LWF precipitation data

precip_temporal_clustering: Event identification for hourly precipitation interpolations, CombiPrecip and LWF precipitation data

precip_statistics_10minres_anywslsite_part1: Statistics (part 1) of 10min precipitation interpolations, CombiPrecip and LWF precipitation data

precip_statistics_10minres_anywslsite_part2: Statistics (part 2) of 10min precipitation interpolations, CombiPrecip and LWF precipitation data

precip_validation_10minres_anywslsite: Bar plot comparison of 10min precipitation interpolations, CombiPrecip and LWF precipitation data

temp_validation_10minres_anywslsite: Plot comparison of 10min temperature interpolations and LWF temperature data

relhum_validation_10minres_anywslsite: Plot comparison of 10min relative humidity interpolations and LWF relative humidity data

globalrad_diurnalcycles_monthly_anywslsite: Plot comparison of mean monthly diurnal cycles of 10min global radiation interpolations and LWF global radiation data

globalrad_validation_10minres_anywslsite: Plot comparison of 10min global radiation interpolations and LWF global radiation data

# Functions
import_data: Import routines for MeteoSwiss station data, LWF station data, CombiPrecip

functions: Other functions

# Writing routines
write_precip_to_csv: Write hourly precipitation interpolations and CombiPrecip to csv files

write_to_csv: Write 10min interpolations to csv files
