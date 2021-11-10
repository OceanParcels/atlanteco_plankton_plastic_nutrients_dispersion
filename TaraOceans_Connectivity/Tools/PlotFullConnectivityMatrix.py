import numpy as np
import matplotlib.pyplot as plt

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/Full_connectivity_output/'

# np.savez_compressed(home_folder + 'Stations_min-T_connectivity.npz', codes=final_stations_code, matrix=min_T_matrix)
data = np.load(home_folder + 'Stations_MTTcount_connectivity_nan.npz', allow_pickle=True)
codes = data['codes']
con_matrix = data['matrix']
print('maximum time: ', np.nanmax(con_matrix))

# Source: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
fig = plt.figure(figsize=(16, 14), dpi=100)
plt.margins(0, 0)
ax = plt.gca()
ax.set_title(
    "Minimum Travel Time (MTT) Connectivity between all stations (No constraints)",
    pad=30)  # min/max Temperature: 15$^\circ$C to 35$^\circ$C Temperature range: 2$^\circ$C
ax.set_xlabel("Destination")
ax.set_ylabel("Source", labelpad=20)
ax.set_xticks(np.arange(len(codes)))
ax.set_yticks(np.arange(len(codes)))
ax.set_xticklabels(codes)
ax.set_yticklabels(codes)
plt.tick_params(axis='both', which='major', labelsize=7)

# Y axis labels on top
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
plt.xticks(rotation=90)
# plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")

# Turn spines off and create white grid.
ax.spines[:].set_visible(False)
ax.set_xticks(np.arange(len(codes) + 1) - .5, minor=True)
ax.set_yticks(np.arange(len(codes) + 1) - .5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
ax.tick_params(which="minor", bottom=False, left=False)

plt.imshow(con_matrix / 12, cmap='turbo')
cbar = plt.colorbar(orientation='vertical')
cbar.set_label('Years')
plt.savefig(home_folder + "Plots/Stations_MTTcount_connectivity_nan1.pdf", bbox_inches='tight',
            pad_inches=0.2)
