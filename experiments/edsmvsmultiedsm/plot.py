#
#    Copyright (C) 2017 Ahmad Retha and Solon P. Pissis.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import numpy
import matplotlib.pyplot as plt

x_axis_labels = []
multiedsm_points = []
edsm_points = []
with open('results.csv') as f:
    f.readline()
    for line in f:
        row = line.split(',')
        x_axis_labels.append(int(row[0]))
        multiedsm_points.append(float(row[1]))
        edsm_points.append(float(row[2]))

# Plot the data for left axis
plt.xticks(x_axis_labels, x_axis_labels)
plt.plot(x_axis_labels, multiedsm_points, '-bs', label='MultiEDSM')
plt.plot(x_axis_labels, edsm_points, '-ro', label='EDSM-BV')
# apply labels and legend
plt.xlabel('No. of patterns')
plt.ylabel('Time (s)')
plt.legend(loc='upper left')
#plt.show()
plt.savefig('edsmvsmultiedsm.pdf')
