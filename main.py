'''
Created on Feb 6, 2013

@author: tiago
'''

from tissue import Tissue

if __name__ == '__main__':
	t = Tissue(3.0, 1, 0.1, 10, 3)
	t.setupDevelopment(1,0.2)
	t.run()