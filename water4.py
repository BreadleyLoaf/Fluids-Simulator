from water3 import Fluid
import pygame

def main():
    resolution = 20
    width = 60
    height = 20
    gravity = 2
    particleRadius = 0.2
    dt = 0.5
    numParticles = 3000
    incompressibilityIters = 30
    f = Fluid(width, height, gravity, dt, numParticles, incompressibilityIters)

    
    pygame.init()
    screen = pygame.display.set_mode((width * resolution, height * resolution))
    pygame.display.set_caption("WE LOVE LULU")

    running = True
    while running:
        screen.fill((0,0,0))

        f.sim()
        for p in range(numParticles):
            pygame.draw.circle(screen, (150,200,255), (int(f.particleX[p] * resolution), int(f.particleY[p] * resolution)), particleRadius * resolution)
        pygame.display.flip()

        #print("Frame")        
        #if input("frame:") == "q":
        #   running = False
    pygame.quit()

main()