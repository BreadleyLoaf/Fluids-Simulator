from water3 import Fluid
import pygame

def main():
    resolution = 5
    width = 100
    height = 100
    gravity = 5
    particleRadius = 0.5
    dt = 1/2
    numParticles = 5000
    incompressibilityIters = 30
    f = Fluid(width, height, gravity, particleRadius, dt, numParticles, incompressibilityIters)

    
    pygame.init()
    screen = pygame.display.set_mode((width * resolution, height * resolution))
    pygame.display.set_caption("WE LOVE LULU")

    running = True
    while running:
        print("Frame")
        screen.fill((0,0,0))

        f.sim()
        for p in range(numParticles):
            #print(int(f.particleX[p] * resolution), int(f.particleY[p] * resolution))
            pygame.draw.circle(screen, (150,200,255), (int(f.particleX[p] * resolution), int(f.particleY[p] * resolution)), particleRadius * resolution)

        pygame.display.flip()

        if input("frame:") == "q":
           running = False
    pygame.quit()

main()