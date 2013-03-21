'''
Created on Feb 6, 2013

@author: tiago
'''

from tissue import Tissue

if __name__ == '__main__':
	t = Tissue(3.0, 1, 0.01, 100, 3)
	t.setupDevelopment(1,0.02)
	t.run()
	
	
	'''
	self.plotter.plotCellStates(self.cs, self.xp, self.yp, self.edges, self.tt.circumcenters, self.tt.triangle_nodes)
		self.plotter.drawCells(self.cs, self.tt.triangle_nodes, self.xp, self.yp)
		self.plotter.next()
		self.plotter = Plotter(d, self.pt)
		self.plotter.setLims(maxX, self.my)
		'''