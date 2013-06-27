'''
Created on Mar 21, 2013

@author: tiago
'''

from plotter import Plotter
import argparse
import re
import os
import cPickle as cp
import numpy as np

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='Plot the cell systems.')
	parser.add_argument('folder', metavar='Folder', type=str, nargs=1,
					help='The folder with simulation data')
	parser.add_argument('--ps', dest='ps', metavar='End frame', nargs='?',
					type=int, const=-1, default=0,
					help="plot simulation results until specified frame, or all frames if not specified (default: no plot)")
	parser.add_argument('-c', dest='c', action="store_true",
					help="Choose this option if the source of the data is the C version of the sim.")
	
	args = parser.parse_args()
	if args.folder[0][-1] == "/" :
		folder = args.folder[0]
	else :
		folder = args.folder[0] + "/"
	
	if not os.path.isdir(folder):
		print "You made a mistake with the folder path. Please double check it. Quitting..."
		exit()
	
	try:
		with open(folder+"params.txt", 'rb') as f:
			params = f.read()
			print params
	except EnvironmentError:
		pass
	
	if args.ps == 0:
		exit()
		
	d = float(re.search("rest length: (.+)\n", params).group(1))
	pt = float(re.search("pressure threshold: (.+)\n", params).group(1))
	mx = float(re.search("max x: (.+)\n", params).group(1))
	my = float(re.search("max y: (.+)\n", params).group(1))
	if args.c:
		steps = int(re.search("ca steps run: (.+)\n", params).group(1))
	else:
		steps = int(re.search("number of ca steps: (.+)\n", params).group(1))
	
	if args.ps > 0 and args.ps <= steps:
		steps = args.ps #plot only until requested frame
	else:
		steps += 1 #plot all states, including final state
			
	pl = Plotter(d, pt, folder)
	pl.setLims(mx, my)
	
	for i in xrange(steps):
		if args.c:
			try:
				filep = open(folder+"state"+str(i).zfill(4)+'.cs')
			except EnvironmentError:
				pl.next()
				continue
			with filep:
				x = np.array(filep.readline().split("; ")[:-1], dtype=float)
				y = np.array(filep.readline().split("; ")[:-1], dtype=float)
				s = np.array(filep.readline().split("; ")[:-1], dtype=float)
				v = np.array(filep.readline().split("; ")[:-1], dtype=float)
				neigh = filep.readline().split("; ")[:-1]
				nnlist = []
				for n in neigh:
					nnlist.append(set(n[1:-1].split(",")[:-1]))
			pl.setupFromArrays(x, y, s, v)
		else:
			try:
				filep = open(folder+"step"+str(i).zfill(4)+'.cs')
			except EnvironmentError:
				pl.next()
				continue
			with filep:
				cs = cp.load(filep)
			pl.setupCellSystem(cs)
		pl.plotCellSystem()
		pl.drawCellStates()
		pl.drawCellPressure()
		pl.drawPressureValues()
		pl.next()
		print "Done with step ", i