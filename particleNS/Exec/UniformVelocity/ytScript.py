# import yt

# ds = yt.load("detonation/plt00450")
# my_sphere = ds.sphere([0.5, 0.5, 0.5], (100, "kpc"))
# plot = yt.ProfilePlot(
#     my_sphere, "x", [("rho")], weight_field=None, accumulation=True
# )
# plot.save()

import yt

# # Load the dataset
# ds = yt.load("detonation/plt00450")

# # Create a line plot of the variables 'u' and 'v' with 1000 sampling points evenly
# plot = yt.LinePlot(
#     ds, [("rho")], (0.0, 0.0, 0.0), (1.0, 1.0, 0.0), 1000
# )

# # Add a legend
# plot.annotate_legend(("rho"))

# # Save the line plot
# plot.save()

es = yt.loaders.load_simulation("detonation", "AMReX")
es.get_time_series(redshifts=[5, 4, 3, 2, 1, 0])