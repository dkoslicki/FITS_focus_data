from fits_stats import *
dir_name = "/Users/dmk333/Dropbox/Astrophotography/FocusTest/"
df = get_stats_dir(dir_name)
df.to_csv(os.path.join(dir_name, "stats.csv"))

