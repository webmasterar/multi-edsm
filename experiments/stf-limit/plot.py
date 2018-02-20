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

import matplotlib.pyplot as plt

x_axis_labels = []
preprocessing_points = []
processing_points = []
memory_points = []
duration_points = []
with open('results.csv') as f:
    f.readline()
    for line in f:
        row = line.split(',')
        x_axis_labels.append(int(row[0]))
        preprocessing_points.append(float(row[1]))
        processing_points.append(float(row[2]))
        memory_points.append(float(row[3]))
        duration_points.append(int(row[4]))

# Plot the data for left axis
fig, ax1 = plt.subplots()
ax1.set_xticks(x_axis_labels)
ax1.set_ylim([4500,7000])
ax1.plot(x_axis_labels, processing_points, '-ro', label='Processing Time')
#ax1.plot(x_axis_labels, duration_points, '-bs', label='Duration')
#plot the data for right axis
ax2 = ax1.twinx()
ax2.set_ylim([100,275])
ax2.plot(x_axis_labels, memory_points, '-g^', label='Memory')
# apply labels and legend
ax1.set_xlabel('Parameter $\\tau$')
ax1.set_ylabel('Time (s)')
ax2.set_ylabel('Memory (MB)')
ax2.yaxis.set_label_position('right')
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
plt.legend(h1+h2, l1+l2, loc='upper right')
#plt.show()
plt.savefig('stflimit.pdf')
