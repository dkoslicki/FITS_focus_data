import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
# Analysis
dir_name = "/Users/dmk333/Dropbox/Astrophotography/FocusTest/"
df = pd.read_csv(os.path.join(dir_name, "stats.csv"))
# filter out weird focus positions
df_filtered = df[df['focus_pos'] > 31250]
# plot focuses
fig, axes = plt.subplots(1, 2)
sns.histplot(data=df_filtered, x='my_FWHM', hue='filter', stat='count', edgecolor=None, ax=axes[0])
sns.histplot(data=df_filtered, x='their_FWHM', hue='filter', stat='count', edgecolor=None, ax=axes[1])
plt.show()
# plot focus positions
sns.histplot(data=df_filtered, x='focus_pos', stat='count', edgecolor=None)
plt.show()
sns.histplot(data=df_filtered, x='temp', stat='count', edgecolor=None)
plt.show()
sns.scatterplot(x=df_filtered['focus_pos'], y=df_filtered["my_FWHM"], hue=df_filtered['filter'])
plt.show()

from mpl_toolkits.mplot3d import Axes3D
color_map = {'B': 'blue',
             "G": "green",
             "Ha": "darkred",
             "OIII": 'teal',
             "OII": "teal",
             "R": "red",
             "SII": "orange"}
c = list(map(lambda x: color_map[x], df_filtered['filter']))
fig = plt.figure(figsize=(6,6))
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)
for filter in color_map.keys():
    df_temp = df_filtered[df_filtered['filter']==filter]
    sc = ax.scatter(df_temp['focus_pos'], df_temp['my_FWHM'], df_temp['temp'], c=color_map[filter], label=filter)
ax.set_xlabel('Focus position')
ax.set_ylabel('FWHM')
ax.set_zlabel('Temp')
plt.legend()
plt.show()
