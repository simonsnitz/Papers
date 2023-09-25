
# PI CHART #

'''
import matplotlib.pyplot as plt
import numpy as np

y = np.array([74,26])

plt.pie(y)
plt.show() 
'''





# VENN DIAGRAM #

'''
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
 
# First way to call the 2 group Venn diagram:
venn2(subsets = (16, 14, 12), set_labels = ('Group A', 'Group B'))
plt.show()
'''





# VERTICAL GROUPED BAR CHART #

'''
import matplotlib.pyplot as plt
import numpy as np

species = ("Information not in database", "No chemical association", "Operon doesn't fit model", "Other")
penguin_means = {
    'Ligify': (42, 20, 5, 3),
    'TFBMiner': (34, 20, 17, 3),
}

x = np.arange(len(species))  # the label locations
width = 0.25  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in penguin_means.items():
    offset = width * multiplier
    rects = ax.barh(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_xlabel('Count')
ax.set_yticks(x + width, species)

ax.legend(loc='upper left', ncols=3)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlim(0, 50)

plt.show()
'''
