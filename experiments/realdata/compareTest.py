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

import os
import sys
import argparse
import subprocess
import multiprocessing

def edsm(pattFile, outFile, edsFile):
    cmd = './../../../edsm/edsm ' + edsFile + ' ' + pattFile
    output = subprocess.Popen([cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    with open(outFile, 'a') as f:
        f.write(output)
        f.write('\n')

def multiedsm(pattFile, outFile, edsFile):
    cmd = './../../multiedsm -m 4g -s ' + edsFile + ' -p ' + pattFile
    output = subprocess.Popen([cmd], shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
    with open(outFile, 'a') as f:
        f.write(output)
        f.write('\n')

def main():
    parser = argparse.ArgumentParser(description='Compare EDSM and Multi-EDSM results')
    parser.add_argument('patterns_file', type=argparse.FileType('r'), help='The patterns file to read')
    parser.add_argument('eds_file', help='The EDS file to read')
    args = parser.parse_args()
    if not os.path.isfile(args.eds_file):
        print 'EDS File missing!'
        sys.exit(1)

    i = 0
    for pattern in args.patterns_file:
        pattFile = 'pattern_' + str(i) + '.pat'
        outFile = 'compare_results/result_' + str(i) + '.txt'
        with open(pattFile, 'w') as f:
            f.write(pattern.strip())
        if os.path.isfile(outFile):
            os.remove(outFile)
        p1 = multiprocessing.Process(target=edsm, args=(pattFile, outFile, args.eds_file,))
        p2 = multiprocessing.Process(target=multiedsm, args=(pattFile, outFile, args.eds_file,))
        p1.start()
        p2.start()
        p1.join()
        p2.join()
        os.remove(pattFile)
        i = i + 1

if __name__ == '__main__':
    main()
