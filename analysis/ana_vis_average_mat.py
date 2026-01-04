import ana_utils
import numpy as np
import os

# Mean SC and FC 
ana_utils.plot_group_matrix('sc', *ana_utils.get_mean_matrix('sc-wei'))
ana_utils.plot_group_matrix('fc', *ana_utils.get_mean_matrix('fc'))
