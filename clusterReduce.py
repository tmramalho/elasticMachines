import os
import re
import numpy as np
import cPickle as cp

home = os.path.expanduser("~")
dir = home + "/Downloads/data"
data = {}
i = 0
for f in os.listdir(dir):
	with open(dir+"/"+f+"/"+"params.txt") as fp:
		params = fp.read()
		ar = int(re.search("Automaton rule: (.+)\n", params).group(1))
		at = float(re.search("average equilibration time: (.+)\n", params).group(1))
		aS = float(re.search("average system size: (.+)\n", params).group(1))
		fs = float(re.search("final system size: (.+)\n", params).group(1))
		se = float(re.search("state entropy: (.+)\n", params).group(1))
		rs = float(re.search("runtime in seconds: (.+)\n", params).group(1))
		ne = float(re.search("network entropy: (.+)\n", params).group(1))
		data[ar] = np.array([ar,at,aS,fs,se,rs,ne])
	i += 1
	print i
with open(home+"/Dropbox/workspace/elasticMachines/data/cluster.pickle", "w") as fp:
	cp.dump(data, fp, cp.HIGHEST_PROTOCOL)