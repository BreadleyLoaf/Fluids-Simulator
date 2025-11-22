import math
import random

class Fluid():
    def __init__(self, width, height, gravity, particleRadius, dt, numParticles, incompressibilityIters, overCompression=1.0):
        self.numCells = width * height
        self.gravity = gravity

        self.width = width
        self.height = height
        self.particleRadius = particleRadius
        self.dt = dt
        self.numParticles = numParticles
        self.incompressibilityIters = incompressibilityIters
        self.overCompression = overCompression

        self.u = [0] * ((width+1) * height)
        self.v = [0] * (width * (height+1))
        self.density = [0] * self.numCells
        
        self.notSolid = [True] * self.numCells
        self.notSolidU = [0] * ((width+1) * height)
        self.notSolidV = [0] * (width * (height+1))

        self.isWater = [False] * self.numCells
        self.isWaterU = [0] * ((width+1) * height)
        self.isWaterV = [0] * (width * (height+1))

        self.particleX = [random.randint(3, width-4) for _ in range(numParticles)]
        self.particleY = [random.randint(3, height-4) for _ in range(numParticles)]
        self.particleU = [0] * numParticles
        self.particleV = [0] * numParticles

        self.initSolids()
    
    def initSolids(self):
        self.notSolid = [True] * self.numCells
        self.notSolidU = [1] * ((self.width+1) * self.height)
        self.notSolidV = [1] * (self.width * (self.height+1))

        for i in range(self.width):
            self.notSolid[self.xy(i, 0)] = False
            self.notSolid[self.xy(i, self.height-1)] = False
        for i in range(self.height):
            self.notSolid[self.xy(0, i)] = False
            self.notSolid[self.xy(self.width-1, i)] = False
        
        for i in range(self.width):
            for j in range(self.height):
                if not self.notSolid[self.xy(i, j)]:
                    self.notSolidU[self.xy(i  , j  , 1)] = 0
                    self.notSolidU[self.xy(i+1, j  , 1)] = 0
                    self.notSolidV[self.xy(i  , j     )] = 0
                    self.notSolidV[self.xy(i  , j+1   )] = 0

    def printAll(self, n):
        if False:
            print("\n step", n)
            print("x", self.particleX)
            print("y", self.particleY)
            print("u", self.particleU)
            print("v", self.particleV)
            print("ugrid", self.u)
            print("vgrid", self.v)

    def sim(self):
        self.printAll(0)

        self.move()

        self.enforceOutOfBounds()

        self.updateWaterUV()

        self.particlesToGrid()
        self.enforceSolidUV()
        self.printAll(1)

        self.enforceIncompressability()
        self.printAll(2)

        self.gridToParticles()
        self.printAll(3)


    def move(self):
        for p in range(self.numParticles):
            self.particleV[p] += (self.gravity * self.dt)
            self.particleX[p] += (self.particleU[p] * self.dt)
            self.particleY[p] += (self.particleV[p] * self.dt)
   
    def enforceOutOfBounds(self):
        for p in range(self.numParticles):
            if self.particleX[p] < 0:
                self.particleX[p] = 0
                self.particleU[p] = 0
            if self.particleX[p] >= self.width:
                self.particleX[p] = self.width - 0.001
                self.particleU[p] = 0

            if self.particleY[p] < 0:
                self.particleY[p] = 0
                self.particleV[p] = 0
            if self.particleY[p] >= self.height:
                self.particleY[p] = self.height - 0.001
                self.particleV[p] = 0

    def xy(self, x, y, isU=0):
        return x + y * (self.width + isU)

    def clamp(self, z, small, big):
        return max(small, min(big, z))

    def clampFloor(self, z, small, big):
        return max(small, min(big, math.floor(z)))

    def enforceSolidUV(self):
        for i in range(self.width):
            for j in range(self.height):
                if not self.notSolid[self.xy(i, j)]:
                    self.u[self.xy(i  , j  , 1)] = 0
                    self.u[self.xy(i+1, j  , 1)] = 0
                    self.v[self.xy(i  , j     )] = 0
                    self.v[self.xy(i  , j+1   )] = 0

    def updateWaterUV(self):
        self.isWater = [False] * self.numCells
        self.isWaterU = [0] * ((self.width+1) * self.height)
        self.isWaterV = [0] * (self.width * (self.height+1))

        for p in range(self.numParticles):
            ix = self.clampFloor(self.particleX[p], 0, self.width-1)
            iy = self.clampFloor(self.particleY[p], 0, self.height-1)
            self.isWater[self.xy(ix, iy)] = True
        
        for i in range(self.width):
            for j in range(self.height):
                if self.isWater[self.xy(i, j)]:
                    self.isWaterU[self.xy(i  , j  , 1)] = 1
                    self.isWaterU[self.xy(i+1, j  , 1)] = 1
                    self.isWaterV[self.xy(i  , j     )] = 1
                    self.isWaterV[self.xy(i  , j+1   )] = 1

    def interp(self, x, y, val, weights, grid, waterGrid, particleToGrid, isU=0):
        x = self.clamp(x, 0, self.width-2)
        y = self.clamp(y, 0, self.height-2)
        ix = self.clampFloor(x, 0, self.width-2)
        iy = self.clampFloor(y, 0, self.height-2)

        tx = x - ix
        ty = y - iy
        sx = 1 - tx
        sy = 1 - ty

        w0 = sx * sy
        w1 = tx * sy
        w2 = tx * ty
        w3 = sx * ty

        if particleToGrid:
            grid[self.xy(ix  , iy  )] += val * w0
            grid[self.xy(ix+1, iy  )] += val * w1
            grid[self.xy(ix+1, iy+1)] += val * w2
            grid[self.xy(ix  , iy+1)] += val * w3

            weights[self.xy(ix  , iy  , isU)] += w0
            weights[self.xy(ix+1, iy  , isU)] += w1
            weights[self.xy(ix+1, iy+1, isU)] += w2
            weights[self.xy(ix  , iy+1, isU)] += w3
        else:
            valid0 = waterGrid[self.xy(ix  , iy  , isU)]
            valid1 = waterGrid[self.xy(ix+1, iy  , isU)]
            valid2 = waterGrid[self.xy(ix+1, iy+1, isU)]
            valid3 = waterGrid[self.xy(ix  , iy+1, isU)]

            d = valid0 * w0 + valid1 * w1 + valid2 * w2 + valid3 * w3

            if d > 0:
                return (valid0 * w0 * grid[self.xy(ix  , iy  , isU)] + 
                        valid1 * w1 * grid[self.xy(ix+1, iy  , isU)] + 
                        valid2 * w2 * grid[self.xy(ix+1, iy+1, isU)] + 
                        valid3 * w3 * grid[self.xy(ix  , iy+1, isU)]) / d
            return 0

    def updateDensity(self):
        self.density = [0] * self.numCells
        useless = [0] * self.numCells
        for p in range(self.numParticles):
            self.interp(self.particleX[p], self.particleY[p], 0, self.density, useless, 0, True)
        
        waterCount = 0
        densityCount = 0
        for i in range(self.width):
            for j in range(self.height):
                if self.isWater[self.xy(i, j)]:
                    waterCount += 1
                    densityCount += self.density[self.xy(i, j)]
        
        return densityCount / waterCount

    def particlesToGrid(self):
        self.u = [0] * ((self.width+1) * self.height)
        self.v = [0] * (self.width * (self.height+1))
        
        weights = [0] * ((self.width+1) * self.height)
        for p in range(self.numParticles):
            self.interp(self.particleX[p] + 0.5, self.particleY[p], self.particleU[p], weights, self.u, 0, True, 1)
        for i in range(self.width+1):
            for j in range(self.height):
                if weights[self.xy(i, j, 1)] > 0:
                    self.u[self.xy(i, j, 1)] /= weights[self.xy(i, j, 1)]
        
        weights = [0] * (self.width * (self.height+1))
        for p in range(self.numParticles):
            self.interp(self.particleX[p], self.particleY[p] + 0.5, self.particleV[p], weights, self.v, 0, True)
        for i in range(self.width):
            for j in range(self.height+1):
                if weights[self.xy(i, j)] > 0:
                    self.v[self.xy(i, j)] /= weights[self.xy(i, j)]

    def gridToParticles(self):
        for p in range(self.numParticles):
            self.particleU[p] = self.interp(self.particleX[p] + 0.5, self.particleY[p], 0, 0, self.u, self.isWaterU, False, 1)
        for p in range(self.numParticles):
            self.particleV[p] = self.interp(self.particleX[p], self.particleY[p] + 0.5, 0, 0, self.v, self.isWaterV, False)

    def calculateTotalDivergence(self):
        totalDivergence = 0
        numSolids = 0
        for x in range(self.width):
            for y in range(self.height):
                if (self.notSolid[self.xy(x, y)]):
                    numSolids += 1
                    totalDivergence += abs(self.u[self.xy(x+1, y  , 1)] - 
                                           self.u[self.xy(x  , y  , 1)] + 
                                           self.v[self.xy(x  , y+1   )] - 
                                           self.v[self.xy(x  , y     )])
        return totalDivergence / numSolids

    def enforceIncompressability(self):
        averageDensity = self.updateDensity()
        for _ in range(self.incompressibilityIters):
            print("Average Divergence: ", self.calculateTotalDivergence())
            for x in range(self.width):
                for y in range(self.height):
                    if (self.notSolid[self.xy(x, y)]):
                        sur = self.notSolidU[self.xy(x+1, y  , 1)]
                        sul = self.notSolidU[self.xy(x  , y  , 1)]
                        svd = self.notSolidV[self.xy(x  , y+1   )]
                        svu = self.notSolidV[self.xy(x  , y     )]
                    
                        divergence = (self.u[self.xy(x+1, y  , 1)] -
                                      self.u[self.xy(x  , y  , 1)] +
                                      self.v[self.xy(x  , y+1   )] - 
                                      self.v[self.xy(x  , y     )])
                                     
                                       
                        s = sul + sur + svd + svu
                        if (s == 0):
                            continue
                        
                        #if averageDensity > 0:
                        #    compression = self.density[self.xy(x, y)] - averageDensity
                        #    if compression > 0:
                        #        divergence -= compression

                        divergence = -divergence / s

                        self.u[self.xy(x+1, y  , 1)] += (divergence * sur)
                        self.u[self.xy(x  , y  , 1)] -= (divergence * sul)
                        self.v[self.xy(x  , y+1   )] += (divergence * svd)
                        self.v[self.xy(x  , y     )] -= (divergence * svu)
