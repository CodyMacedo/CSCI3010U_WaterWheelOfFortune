'''
each particle will be contained within a list of all particles
the list of all particles is processed with a number of threads
threads will each handle 100 out of n total particles, so 
for i in range(0, n/100):
	for j in range(i*100, i*100+100):
		thread ODE particle[i*100+j]
screen will be split up into grid, particles will place self
into grid based on position
collision handler will process collisions in each grid space, also threaded 