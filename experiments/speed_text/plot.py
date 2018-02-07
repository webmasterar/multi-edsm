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
processing_points = []
logprocessing_points = []
with open('results.csv') as f:
    f.readline()
    for line in f:
        row = line.split(',')
        x_axis_labels.append(int(row[0]))
        processing_points.append(float(row[1]))
        logprocessing_points.append(float(row[2]))

plt.xticks(range(len(x_axis_labels)), x_axis_labels)
plt.plot(range(len(x_axis_labels)), logprocessing_points, '-ro')
#plt.title('Speed of MultiEDSM with increasing text size')
plt.xlabel('Total text size (N)')
plt.ylabel('Log processing time (s)')
#plt.show()
plt.savefig('speed_text.pdf')

