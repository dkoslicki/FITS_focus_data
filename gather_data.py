from fits_stats import *

# Just process a single file to test for speed
from time import time
fits_file_name = "/Users/dmk333/Dropbox/Astrophotography/FocusTest/L_Ha_2509_300s.fit"
t1 = time()
res = get_star_stats_fits(fits_file_name, plot_stats=False)
t2 = time()
print(f"Time to process a single file:{t2-t1}")
input("Press Enter to continue...")

dir_name = "/Users/dmk333/Dropbox/Astrophotography/FocusTest/"
df = get_stats_dir(dir_name)
df.to_csv(os.path.join(dir_name, "stats.csv"))

