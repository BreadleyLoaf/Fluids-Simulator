import numpy as np
import pygame
import time

def main():
    rows = 100
    cols = 200
    u = np.zeros((rows+1, cols), dtype=np.float64)
    v = np.zeros((rows, cols+1), dtype=np.float64)

    gravity = np.array([0.98, 0], dtype=np.float64)
    timeStep = 1

    nonSolid = np.ones((rows, cols), dtype=int)
    nonSolid[:,0] = 0
    nonSolid[:,cols-1] = 0
    nonSolid[0,:] = 0
    nonSolid[rows-1,:] = 0
    nonSolidMask = nonSolid.astype(bool)
    nonSolidU = np.ones((rows+1, cols), dtype=int)
    nonSolidV = np.ones((rows, cols+1), dtype=int)
    nonSolidNeighbours = np.zeros((rows, cols), dtype=int)
    def updateSolids(nonSolidNeighbours, nonSolidMask, nonSolidU, nonSolidV):
        nonlocal nonSolid

        nonSolidMask[:,:] = nonSolid.astype(bool)

        nonSolidNeighbours[1:,:] = nonSolid[:-1:]
        nonSolidNeighbours[:-1,:] += nonSolid[1:,:]
        nonSolidNeighbours[:,1:] += nonSolid[:,:-1]
        nonSolidNeighbours[:,:-1] += nonSolid[:,1:]

        nonSolidU[:,:] = 1
        nonSolidV[:,:] = 1
        nonSolidU[1:,:] = nonSolid
        nonSolidU[:-1,:] *= nonSolid
        nonSolidV[:,1:] = nonSolid
        nonSolidV[:,:-1] *= nonSolid

    density = np.zeros((rows, cols), dtype=np.float64)
    def stopSolidBorders(u, v):
        nonlocal nonSolidU, nonSolidV

        u[~nonSolidU.astype(bool)] = 0
        v[~nonSolidV.astype(bool)] = 0
        density[~nonSolidMask] = 0

    def handleGravity(u, v):
        u += gravity[0] * timeStep
        v += gravity[1] * timeStep

    #red-black gauss seidel method for divergence handling
    I, J = np.indices(nonSolid.shape)
    redM = ((I + J) % 2 == 1).astype(bool) #    0,0 is false
    blackM = ((I + J) % 2 == 0).astype(bool) #  0,0 is true
    def handleDivergence(u, v):
        nonlocal nonSolid, redM, blackM, nonSolidNeighbours, nonSolidMask, nonSolidU, nonSolidV
        divergence = np.zeros((rows, cols), dtype=np.float64)
        redMask = redM & nonSolidMask
        blackMask = blackM & nonSolidMask

        for i in range(100): #increase this number to reduce divergence

            print("u original 1", u[98,1], u[99,1], v[98,2], v[98,1], nonSolidNeighbours[98,1])
            print(                          (u[98,1]-u[99,1]+v[98,2]-v[98,1])/nonSolidNeighbours[98,1])
            print(u[98,1]+(nonSolidU[98,1]* (u[98,1]-u[99,1]+v[98,2]-v[98,1])/nonSolidNeighbours[98,1]))
            print(u[98,1]-(nonSolidU[98,1]* (u[98,1]-u[99,1]+v[98,2]-v[98,1])/nonSolidNeighbours[98,1]))
            print(       -(nonSolidU[98,1]* (u[98,1]-u[99,1]+v[98,2]-v[98,1])/nonSolidNeighbours[98,1]))
            print(nonSolidU[98,1])
            print("")
            divergence[redMask] = ( u[1:,:][redMask] - u[:-1,:][redMask] + v[:,1:][redMask] - v[:,:-1][redMask] ) / nonSolidNeighbours[redMask]

            print(u)
            u[:-1,:][redMask] += divergence[redMask] * nonSolidU[:-1,:][redMask]
            print(divergence)
            print(u)

            u[1:,:][redMask] -= divergence[redMask] * nonSolidU[1:,:][redMask]

            v[:,:-1][redMask] += divergence[redMask] * nonSolidV[:,:-1][redMask]
            v[:,1:][redMask] -= divergence[redMask] * nonSolidV[:,1:][redMask]
            print("u original 2", u[98,1], u[99,1], v[98,2], v[98,1], nonSolidNeighbours[98,1])
            print(                          (u[98,1]-u[99,1]+v[98,2]-v[98,1])/nonSolidNeighbours[98,1])
            print(u[98,1]+(nonSolidU[98,1]* (u[98,1]-u[99,1]+v[98,2]-v[98,1])/nonSolidNeighbours[98,1]))
            print(u[98,1]-(nonSolidU[98,1]* (u[98,1]-u[99,1]+v[98,2]-v[98,1])/nonSolidNeighbours[98,1]))
            print(nonSolidU[98,1])
            print("")

            divergence[blackMask] = ( u[1:,:][blackMask] - u[:-1,:][blackMask] + v[:,1:][blackMask] - v[:,:-1][blackMask] ) / nonSolidNeighbours[blackMask]
            u[:-1,:][blackMask] += divergence[blackMask] * nonSolidU[:-1,:][blackMask]
            u[1:,:][blackMask] -= divergence[blackMask] * nonSolidU[1:,:][blackMask]
            v[:,:-1][blackMask] += divergence[blackMask] * nonSolidV[:,:-1][blackMask]
            v[:,1:][blackMask] -= divergence[blackMask] * nonSolidV[:,1:][blackMask]

    def testDivergence():
        #Test case for handleDivergence
        u = np.random.randn(rows + 1, cols) * 0.1   # x-face velocities
        v = np.random.randn(rows, cols + 1) * 0.1   # y-face velocities
        updateSolids(nonSolidNeighbours, nonSolidMask, nonSolidU, nonSolidV)
        stopSolidBorders(u, v)
        div_before = (u[1:,:] - u[:-1,:]) + (v[:,1:] - v[:,:-1])
        handleDivergence(u, v)
        div_after = (u[1:,:] - u[:-1,:]) + (v[:,1:] - v[:,:-1])
        print("Mean absolute divergence before:", np.sum(np.abs(div_before)) / np.sum(nonSolid))
        print("Mean absolute divergence after: ", np.sum(np.abs(div_after)) / np.sum(nonSolid))
        print("Reduction factor:", np.sum(np.abs(div_before)) / np.sum(np.abs(div_after)))
    
    # semi-legrangian advection method
    def handleAdvection(density):
        nonlocal u, v

        #face centered velocities
        velocityU = (u[1:,:] + u[:-1,:]) * 0.5
        velocityV = (v[:,1:] + v[:,:-1]) * 0.5

        X, Y = np.meshgrid(np.arange(cols), np.arange(rows))
        oldX = X - velocityU * timeStep
        oldY = Y - velocityV * timeStep

        oldX = np.clip(oldX, 0, cols - 1.001)
        oldY = np.clip(oldY, 0, rows - 1.001)

        oldDensity = density.copy()

        offsetX, oldXLeft = np.modf(oldX)
        oldXLeft = oldXLeft.astype(int)
        oldXRight = oldXLeft + 1

        offsetY, oldYBottom = np.modf(oldY)
        oldYBottom = oldYBottom.astype(int)
        oldYTop = oldYBottom + 1

        density[:,:] =  (1-offsetX) * (1-offsetY) * oldDensity[oldYBottom, oldXLeft] + \
                        offsetX * (1-offsetY) *     oldDensity[oldYBottom, oldXRight] + \
                        (1-offsetX) * offsetY *     oldDensity[oldYTop, oldXLeft] + \
                        offsetX * offsetY *         oldDensity[oldYTop, oldXRight]
    
    pygame.init()
    screen = pygame.display.set_mode((rows, cols))
    
    densityMultiplier = 100
    density = np.random.randn(rows, cols) * densityMultiplier

    updateSolids(nonSolidNeighbours, nonSolidMask, nonSolidU, nonSolidV)

    running = True
    count = 0
    while running:
        
        handleGravity(u, v)
        stopSolidBorders(u, v)
        oldU = u.copy()
        handleDivergence(u, v)

        if False:
            print("\nframe\n")
            print(oldU)
            print("")
            print(u-oldU)
        
        handleAdvection(density)
        densitySurface = np.uint8(255 / (1 + np.exp(-density/densityMultiplier)))

        surface = pygame.surfarray.make_surface(np.stack([densitySurface]*3, axis=-1))
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        screen.blit(surface, (0, 0))
        pygame.display.flip()
        
        time.sleep(0.2)
        count += 1

        #for j in range(95,101):
        #    for i in range(0,3):
        #        print(u[j][i], end='\t')
        #    print("")
        #print("\n\n")
        if count > 50:
            running = False

    pygame.quit()

main()