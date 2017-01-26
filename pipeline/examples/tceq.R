source('pipeline.R')

pipeline.load_TCEQ_libs()

dataframe <- load.file_as_dataframe("examples/raw_data/TCEQ/2007_raw.csv", maximum_rows=100000)

stats <- TCEQ.generate_air_pollution_statistics(dataframe, region=NA, pollutants=NA, AQS_codes=NA, stats=NA)

stats$hourly_averages
stats$hourly_maximums
stats$eight_hour_rolling_averages
stats$eight_hour_maximums