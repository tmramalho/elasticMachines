'''
Created on Feb 6, 2013

@author: tiago
'''

from tissue import Tissue
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Simulate a morphogenetic cellular automaton.')
	parser.add_argument('--k', dest='k',
					type=float, default=3.0,
					help="elastic constant for the cells")
	parser.add_argument('--l', dest='l',
					type=float, default=1.0,
					help="dampening constant")
	parser.add_argument('--d', dest='d',
					type=float, default=0.01,
					help="equilibrium intra nuclei distance")
	parser.add_argument('--nx', dest='nx',
					type=int, default=100,
					help="initial number of cells in x dim")
	parser.add_argument('--ny', dest='ny',
					type=int, default=3,
					help="initial number of cells in x dim")
	parser.add_argument('--ns', dest='ns',
					type=int, default=10,
					help="number of cellular automata iterations")
	parser.add_argument('--rid', dest='rid',
					type=int, default=19574,
					help="integer describing the cellular automata transition table.Must be < 2^16, else mod will be taken")
	parser.add_argument('--scan', dest='scan', nargs=2, metavar=('beginID', 'endID'),
					type=int, default=(0,0),
					help="Scan a range of rules, if selected --rid will be ignored.")
	parser.add_argument('--xe', dest='xe',
					type=float, default=1.0,
					help="Cells will be initially uniformely distributed with positions in range [0, x_end]")
	parser.add_argument('--ye', dest='ye',
					type=float, default=0.02,
					help="Cells will be initially uniformely distributed with positions in range [0, x_end]")
	
	
	args = parser.parse_args()
	if args.scan[0] == args.scan[1]:
		t = Tissue(args.k, args.l, args.d, args.nx, args.ny, args.ns, args.rid)
		t.setupDevelopment(args.xe, args.ye)
		t.run()
	else:
		for i in range(args.scan[0], args.scan[1]+1):
			t = Tissue(args.k, args.l, args.d, args.nx, args.ny, args.ns, i)
			t.setupDevelopment(args.xe, args.ye)
			t.run()