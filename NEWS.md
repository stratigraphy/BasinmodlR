# Changes in version 0.1.0

first release

added functions:

run_basin_model
runs a one-dimensional burial, compaction and thermal history model from a
table of stratigraphic units and a table of erosion events

plot_basin_model
composite figure of the burial path, the temperature raster and the hydrocarbon
generation windows (also available as plot method)

plot_event_chart
petroleum system event chart based on the result of run_basin_model

added methods:

print.basin_model
plot.basin_model

added data sets:

`brabant`, `condroz` and `ardennes`
example basins for the Brabant Massif, the Condroz Inlier (Huy) and the
Cambro-Ordovician inliers of the Ardennes
