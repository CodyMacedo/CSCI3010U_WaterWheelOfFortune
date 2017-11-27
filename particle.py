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
'''


import pygame, sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode
import random as rand
import math


# set up the colors
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)
g = -9.8


def normalize(v):
    return v / np.linalg.norm(v)


class Particle(pygame.sprite.Sprite):

    def __init__(self, imgfile, radius, mass=1.0):
        pygame.sprite.Sprite.__init__(self)


        self.image = pygame.image.load(imgfile)
        self.image = pygame.transform.scale(self.image, (radius, radius)) 
        self.state = [0, 0, 0, 0]
        self.mass = mass
        self.t = 0
        self.radius = radius


        self.solver = ode(self.f)
        self.solver.set_integrator('dop853')
        self.solver.set_initial_value(self.state, self.t)


    def f(self, t, y):
        y[3] = g * self.mass
        return [y[2], y[3], 0, 0]


    def set_pos(self, pos):
        self.state[0:2] = pos
        self.solver.set_initial_value(self.state, self.t)
        return self


    def set_vel(self, vel):
        #print vel, 'fsdf'

        self.state[2:] = vel
        self.solver.set_initial_value(self.state, self.t)
        return self


    def update(self, dt):
        self.t += dt
        self.state = self.solver.integrate(self.t)


    def move_by(self, delta):
        self.state[0:2] = np.add(self.pos, delta)
        return self


    def draw(self, surface):
        rect = self.image.get_rect()
        rect.center = (self.state[0], 600-self.state[1]) # Flipping y
        surface.blit(self.image, rect)


    def pprint(self):
        print 'Particle', self.state

class Wheel(pygame.sprite.Sprite):
    def __init__(self, center, radius, mass=1000):
        pygame.sprite.Sprite.__init__(self)


        self.state = [0,0]
        self.lines = []
        self.mass = mass
        self.t = 0
        self.center = center
        self.radius = radius
        self.angle = 0


        self.solver = ode(self.f)
        self.solver.set_integrator('dop853')
        self.solver.set_initial_value(self.state, self.t)


    def f(self, t, y):
        return [y[0], y[1]]


    def set_vel(self, vel):
        self.state[0] = vel
        self.solver.set_initial_value(self.state, self.t)
        return self


    def update(self, dt):
        self.t += dt
        self.state = self.solver.integrate(self.t)


    def draw(self, surface):
        self.angle = (self.angle - 1) % 360

        self.angle +=1.5

        for i in range(0,359, 45):
            x = self.center[0] + math.cos(math.radians(self.angle + i)) * self.radius
            y = self.center[1] + math.sin(math.radians(self.angle + i)) * self.radius
            self.lines.append(pygame.draw.line(surface, BLACK, self.center, (x,y)))

        self.circle = pygame.draw.circle(surface, BLACK, self.center, (int)(self.radius*.7), 1)
        

    def pprint(self):
        print 'Wheel', self.state

class World:

    def __init__(self, height, width):
        self.particles = []
        self.wheels =[]
        self.height = height
        self.width = width
        self.e = 1. # Coefficient of restitution


    def add(self, imgfile, radius, mass=1.0):
        particle = Particle(imgfile, radius, mass)
        self.particles.append(particle)
        return particle

    def addWheel(self, centre, radius):
        wheel = Wheel(centre, radius)
        self.wheels.append(wheel)
        return wheel


    def pprint(self):
        print '#particles', len(self.particles)
        for d in self.particles:
            d.pprint()


    def draw(self, screen):
        for d in self.particles:
            d.draw(screen)
        for w in self.wheels:
            w.draw(screen)


    def update(self, dt):
        self.check_for_collision()

        for d in self.particles:
            d.update(dt)


    def compute_collision_response(self, i, j):
        pass


    def check_for_collision(self):
        for i in range(0, len(self.particles)):
            if (self.particles[i].state[0] - self.particles[i].radius <= 0 or
                self.particles[i].state[0] + self.particles[i].radius >= 800):
                self.particles[i].state[2] *= -1
            elif (self.particles[i].state[1] - self.particles[i].radius <= 0 or
                self.particles[i].state[1] + self.particles[i].radius >= 600):
                self.particles[i].state[3] *= -1
            else:
                for j in range(i+1, len(self.particles)):
                    if i == j:
                        continue
                    #print 'Checking particles', i, 'and', j
                    pos_i = np.array(self.particles[i].state[0:2])
                    pos_j = np.array(self.particles[j].state[0:2])
                    dist_ij = np.sqrt(np.sum((pos_i - pos_j)**2))


                    #print pos_i, pos_j, dist_ij


                    radius_i = self.particles[i].radius
                    radius_j = self.particles[j].radius
                    if dist_ij > radius_i + radius_j:
                        continue


                    # May be a collision
                    vel_i = np.array(self.particles[i].state[2:])
                    vel_j = np.array(self.particles[j].state[2:])
                    relative_vel_ij = vel_i - vel_j
                    n_ij = normalize(pos_i - pos_j)


                    #print relative_vel_ij, n_ij


                    if np.dot(relative_vel_ij, n_ij) >= 0:
                        continue


                    # Ouch!
                    #print 'Collision between particles', i, 'and', j, '!!'
                    mass_i = self.particles[i].mass
                    mass_j = self.particles[j].mass

                    # Don't confuse this J with j
                    J = -(1+self.e) * np.dot(relative_vel_ij, n_ij) / ((1./mass_i) + (1./mass_j))


                    vel_i_aftercollision = vel_i + n_ij * J / mass_i
                    vel_j_aftercollision = vel_j - n_ij * J / mass_j


                    #print 'Response'
                    #print vel_i_aftercollision.shape, vel_j_aftercollision.shape

                    self.particles[i].set_vel(vel_i_aftercollision)
                    self.particles[j].set_vel(vel_j_aftercollision)
                    break # Only handle a single collision per instance




def main():

   # initializing pygame
    pygame.init()


    clock = pygame.time.Clock()


    # top left corner is (0,0)
    win_width = 800     # 500 cm = 5 m
    win_height = 600
    screen = pygame.display.set_mode((win_width, win_height))
    pygame.display.set_caption('Water Wheel of Fortune')


    world = World(win_height, win_width)
    
    for i in range(100):
        newRadius = 3 # smallest radius is 3, anything smaller is invisible
        newMass = 1
        # make sure particle is within the walls
        newPos = np.array([rand.uniform(0 + newRadius, 800 - newRadius), rand.uniform(0 + newRadius, 600 - newRadius)])
        newVel = np.array([0, 0])
        
        world.add('sandparticle.png', newRadius, newMass).set_pos(newPos).set_vel(newVel)

    world.addWheel([400, 300], 200)

   
    dt = 0.1
    pause = False

    while True:
        # 30 fps
        if not pause:
            clock.tick(30)
        
        event = pygame.event.poll()
        if event.type == pygame.QUIT:
            sys.exit(0)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_q:
            pygame.quit()
            sys.exit(0)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_p:
            pause = not pause
        else:
            pass

        if not pause:
            # Clear the background, and draw the sprites
            screen.fill(WHITE)
            world.draw(screen)
            world.update(dt)

            pygame.display.update()


if __name__ == '__main__':
    main()