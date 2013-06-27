'''
Created on Apr 24, 2013

@author: tiago
'''
import cPickle as cp
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
	f = open("data/cluster.pickle", 'r')
	fw = open("data/cluster.csv", 'w')
	fw.write("rid, eqtime, avsize, finalsize, stateent, runtime, netent\n")
	data = cp.load(f)
	nRules = len(data)
	atEntropy = np.zeros((nRules))
	aSEntropy = np.zeros((nRules))
	fsEntropy = np.zeros((nRules))
	seEntropy = np.zeros((nRules))
	rsEntropy = np.zeros((nRules))
	neEntropy = np.zeros((nRules))
	for rid in sorted(data):
		fw.write('%d, %f, %f, %d, %f, %f, %f\n' % 
				(data[rid][0], data[rid][1], data[rid][2], data[rid][3], data[rid][4], data[rid][5], data[rid][6]))
		atEntropy[rid] = data[rid][1]
		aSEntropy[rid] = data[rid][2]
		fsEntropy[rid] = data[rid][3]
		seEntropy[rid] = data[rid][4]
		rsEntropy[rid] = data[rid][5]
		neEntropy[rid] = data[rid][6]
	fw.close()
	n, bins, patches = plt.hist(atEntropy, 50, facecolor='green', alpha=0.75)
	plt.title("Average Equilibration Time Histogram")
	plt.savefig("avEqTime.pdf")
	plt.clf()
	n, bins, patches = plt.hist(aSEntropy, 50, facecolor='green', alpha=0.75)
	plt.title("Average System Size Histogram")
	plt.savefig("avSysSize.pdf")
	plt.clf()
	n, bins, patches = plt.hist(fsEntropy, 50, facecolor='green', alpha=0.75)
	plt.title("Final System Size Histogram")
	plt.savefig("finSysSize.pdf")
	plt.clf()
	n, bins, patches = plt.hist(seEntropy, 50, facecolor='green', alpha=0.75)
	plt.title("State Entropy Histogram")
	plt.savefig("stateEnt.pdf")
	plt.clf()
	n, bins, patches = plt.hist(rsEntropy, 50, facecolor='green', alpha=0.75)
	plt.title("Runtime Histogram")
	plt.savefig("runTime.pdf")
	plt.clf()
	n, bins, patches = plt.hist(neEntropy, 50, facecolor='green', alpha=0.75)
	plt.title("Network Entropy Histogram")
	plt.savefig("netEnt.pdf")
	plt.clf()
	
	plt.scatter(seEntropy, neEntropy, marker='.', lw = 0, s=2, facecolor='0.5')
	plt.xlabel("State Entropy")
	plt.ylabel("Network Entropy")
	plt.savefig("entScatter.png")
	plt.clf()
	
	plt.scatter(neEntropy, fsEntropy, marker='.', lw = 0, s=2, facecolor='0.5')
	plt.xlabel("Network Entropy")
	plt.ylabel("Final Size")
	plt.savefig("sizeEntScatter.png")
	plt.clf()
	
	plt.scatter(seEntropy, fsEntropy, marker='.', lw = 0, s=2, facecolor='0.5')
	plt.xlabel("State Entropy")
	plt.ylabel("Final Size")
	plt.savefig("sizeStateScatter.png")
	plt.clf()
	