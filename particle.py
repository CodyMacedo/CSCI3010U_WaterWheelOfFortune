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

win_width = 800     # 500 cm = 5 m
win_height = 600

# set up the colors
BLACK = (0, 0, 0)
GREY = (150,150,150)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)


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
        self.gravity = -9.8


        self.solver = ode(self.f)
        self.solver.set_integrator('dop853')
        self.solver.set_initial_value(self.state, self.t)
        self.solver.set_f_params(self.gravity)


    def f(self, t, y, arg1):
        return [y[2], y[3], 0, arg1]


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
        rect.center = (self.state[0], win_height-self.state[1]) # Flipping y
        surface.blit(self.image, rect)


    def pprint(self):
        print 'Particle', self.state




class Wheel(pygame.sprite.Sprite):
    def __init__(self, center, radius, mass=1000):
        pygame.sprite.Sprite.__init__(self)


        self.state = np.zeros(6)
        self.state[0:2] = np.zeros(2) # position
        self.state[2:4] = np.zeros(2) # angular velocity
        self.state[4:6] = np.zeros(2) # angular momentum
        self.lines = []
        self.mass = mass
        self.t = 0
        self.center = center
        self.radius = radius
        self.angle = 0

        self.torque = 0
        

        self.solver = ode(self.f)
        self.solver.set_integrator('dop853')
        self.solver.set_initial_value(self.state, self.t)


    def f(self, t, y):
        return [y[0], y[1]]


    def set_vel(self, vel):
        self.state[2:4] = vel
        self.solver.set_initial_value(self.state, self.t)
        return self


    def update(self, dt):
        self.t += dt
        self.state = self.solver.integrate(self.t)


    def draw(self, surface):
        self.angle = (self.angle - 1) % 360

        self.angle +=1.5

        for i in range(0,316, 45):
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
        self.e = 0.1 # Coefficient of restitution


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
        self.check_outside_screen()

        for d in self.particles:
            d.update(dt)


    def compute_collision_response(self, i, j):
        pass


    def check_outside_screen(self):
        self.particles = [x for x in self.particles if self.outside_screen(x)]


    def outside_screen(self, particle):
        if (particle.state[0] < -particle.radius):
            return False
        elif (particle.state[0] > win_width + particle.radius):
            return False
        elif (particle.state[1] < -particle.radius):
            return False
        else:
            return True


    def check_for_collision(self):
        for i in range(0, len(self.particles)):
            if (self.particles[i].state[0] - self.particles[i].radius <= 0 or
                self.particles[i].state[0] + self.particles[i].radius >= 800):
                self.particles[i].state[2] *= -1*self.e
            elif (self.particles[i].state[1] - self.particles[i].radius <= 0):
                self.particles[i].state[3] = 0

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
                # break # Only handle a single collision per instance


            # check for particle - wheel collission
            for j in range(0, len(self.wheels)):
                #print 'Checking particles', i, 'and wheel', j
                pos_i = np.array(self.particles[i].state[0:2])
                pos_j = np.array(self.wheels[j].center)
                dist_ij = np.sqrt(np.sum((pos_i - pos_j)**2))

                radius_i = self.particles[i].radius
                radius_j = self.wheels[j].radius*.7
                if dist_ij > radius_i + radius_j:
                    continue


                # May be a collision
                vel_i = np.array(self.particles[i].state[2:])
                vel_j = 0
                relative_vel_ij = vel_i - vel_j
                n_ij = normalize(pos_i - pos_j)


                #print relative_vel_ij, n_ij


                if np.dot(relative_vel_ij, n_ij) >= 0:
                    continue


                # Ouch!
                #print 'Collision between particles', i, 'and', j, '!!'
                mass_i = self.particles[i].mass
                mass_j = self.wheels[j].mass

                # Don't confuse this J with j
                J = -(1+self.e) * np.dot(relative_vel_ij, n_ij) / ((1./mass_i) + (1./mass_j))


                vel_i_aftercollision = vel_i + n_ij * J / mass_i


                #print 'Response'
                #print vel_i_aftercollision.shape, vel_j_aftercollision.shape

                self.particles[i].set_vel(vel_i_aftercollision)
                # break # Only handle a single collision per instance
                
                
                # ANGULAR COLISION #
                
                # detect collision with lines on wheel
                for x in range(len(wheels[j].lines)):
                    line = self.wheels.lines[x]
                    A = self.wheels[j].center
                    C = self.particles[i].state[0:2]
                    
                    if A == line.topleft:
                        B = line.bottomright
                    elif A == line.bottomright:
                        B = line.topleft
                    elif A == line.topright:
                        B = line.bottomleft
                    else:
                        B = line.topright
                    
                    dist = np.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2)
                    
                    Dx = (B[0]-A[0])/dist
                    Dy = (B[1]-A[1])/dist
                    
                    t = Dx*(C[0]-A[0])+Dy*(C[1]-A[1])
                    
                    Ex = t*Dx+A[0]
                    Ey = t*Dy+A[1]
                    
                    dist2 = np.sqrt((Ex-C[0])**2+(Ey-C[1])**2)
                    
                    if (dist2 < self.particles[i].radius):
                        #Do conservation of momentum for angular momentum
                    
                        # Nothing I tried worked, I'll try to input the equations tomorrow morning
                

               



def main():

   # initializing pygame
    pygame.init()


    clock = pygame.time.Clock()


    # top left corner is (0,0)
    screen = pygame.display.set_mode((win_width, win_height))
    pygame.display.set_caption('Water Wheel of Fortune')


    world = World(win_height, win_width)

    world.addWheel([400, 300], 200)

   
    pause = False
    rain = False

    dt = 0.3
    pRadius = 10 # smallest radius is 3, anything smaller is invisible
    pMass = 1

    pygame.time.set_timer(pygame.USEREVENT + 1, 100)

    if rain:
        range = [0 + pRadius, win_width - pRadius]
    else:
        spoutPos = 380
        spoutWidth = 40
        range = [spoutPos + pRadius, spoutPos + spoutWidth + pRadius]
    
    print "\n\nPress P key to pause or resume"
    print "Press R key to toggle rain or spout"
    print "Press A or D keys to move spout left or right\n\n"



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
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_r:
            rain = not rain
            if rain: 
                range = [0 + pRadius, win_width - pRadius]
            else:
                range = [spoutPos + pRadius, spoutPos + spoutWidth + pRadius]
        elif event.type == pygame.USEREVENT + 1:
            # new particle
            if (len(world.particles) < 100):
                # make sure particle is within the walls
                newPos = np.array([rand.uniform(range[0],range[1]), win_height])
                newVel = np.array([0, 0])
                
                world.add('sandparticle.png', pRadius, pMass).set_pos(newPos).set_vel(newVel)
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_a and not rain:
            spoutPos -= 10
            range = [spoutPos + pRadius, spoutPos + spoutWidth + pRadius]
        elif event.type == pygame.KEYDOWN and event.key == pygame.K_d and not rain:
            spoutPos += 10
            range = [spoutPos + pRadius, spoutPos + spoutWidth + pRadius]
        else:
            pass

        if not pause:
            # Clear the background, and draw the sprites
            screen.fill(WHITE)
            world.draw(screen)
            world.update(dt)

            if not rain:
                pygame.draw.rect(screen, GREY, (spoutPos, 0, spoutWidth*1.5, 30))

            pygame.display.update()


if __name__ == '__main__':
    main()