I collected a few routines that are useful to work with XPS files from SpecsLab in Python.

This folder is a Python Module, so copy it to a place where your interpreter can find it,
e.g. C:\Users\XYZ\.ipython\XPS



you can then call the functions like this:


from XPS.get_region import get_region
from XPS.shirley import remove_shirley_background

#load data from SpecsLab xml file and select region
region = get_region(filename, groupNo, regionNo)

E = region.x_be    		#binding energy scale
s = region.y_avg_counts_mcd     #signal