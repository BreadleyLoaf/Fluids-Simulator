import math
import random

class Fluid():
    def __init__(self, width, height, gravity, dt, numParticles, incompressibilityIters, overCompression=1.0):
        self.numCells = width * height
        self.numCellsU = (width+1) * height
        self.numCellsV = width * (height+1)
        self.gravity = gravity * dt

        self.width = width
        self.height = height
        self.dt = dt
        self.numParticles = numParticles
        self.incompressibilityIters = incompressibilityIters
        self.overCompression = overCompression

        self.u = [0] * self.numCellsU
        self.v = [0] * self.numCellsV
        self.density = [0] * self.numCells
        
        self.notSolid = [True] * self.numCells
        self.notSolidU = [0] * self.numCellsU
        self.notSolidV = [0] * self.numCellsV

        self.isWater = [False] * self.numCells
        self.isWaterU = [0] * self.numCellsU
        self.isWaterV = [0] * self.numCellsV

        self.particleX = [random.uniform(3, int(width/3*2)-4) for _ in range(numParticles)]
        self.particleY = [random.uniform(3, height-4) for _ in range(numParticles)]
        self.particleU = [0] * numParticles
        self.particleV = [0] * numParticles

        self.initSolids()
    
    def initSolids(self):
        self.notSolid = [True] * self.numCells
        self.notSolidU = [1] * self.numCellsU
        self.notSolidV = [1] * self.numCellsV

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

    def sim(self):
        self.move(self.dt, self.width, self.height)
        self.updateWaterUV(self.width, self.height)
        self.particlesToGrid(self.width, self.height)
        self.enforceSolidUV()
        averageDensity = self.updateDensity(self.width, self.height)
        self.enforceIncompressability(averageDensity)
        self.gridToParticles(self.width, self.height)

    def move(self, dt, width, height):
        for p in range(self.numParticles):
            self.particleV[p] += (self.gravity)
            self.particleX[p] += (self.particleU[p] * dt)
            self.particleY[p] += (self.particleV[p] * dt)

        width -= 1
        height -= 1
        for p in range(self.numParticles):
            if self.particleX[p] < 1:
                self.particleX[p] = 1
                self.particleU[p] = 0
            elif self.particleX[p] > width:
                self.particleX[p] = width
                self.particleU[p] = 0

            if self.particleY[p] < 0:
                self.particleY[p] = 0
                self.particleV[p] = 0
            elif self.particleY[p] > height:
                self.particleY[p] = height
                self.particleV[p] = 0

    def xy(self, x, y, isU=0):
        return x + y * (self.width + isU)

    def enforceSolidUV(self):
        for i in range(self.numCellsU):
            if not self.notSolidU[i]:
                self.u[i] = 0
        for i in range(self.numCellsV):
            if not self.notSolidV[i]:
                self.v[i] = 0

    def updateWaterUV(self, width, height):
        isWater = [False] * self.numCells
        isWaterU = [0] * self.numCellsU
        isWaterV = [0] * self.numCellsV

        for p in range(self.numParticles):
            ix = int(self.particleX[p])
            iy = int(self.particleY[p])
            
            locV = ix + iy * width
            if not isWater[locV]:
                locU = ix + iy * (width + 1)

                isWater[locV] = True
                isWaterU[locU      ] = 1
                isWaterU[locU+1    ] = 1
                isWaterV[locV      ] = 1
                isWaterV[locV+width] = 1
        
        self.isWater = isWater
        self.isWaterU = isWaterU
        self.isWaterV = isWaterV

    def interp(self, width, height, x, y, val, weights, grid, waterGrid, particleToGrid, isU=0, fillGrid=True):
        x = min(x, width-2+isU)
        y = min(y, height-2)
        ix = int(x)
        iy = int(y)

        tx = x - ix
        ty = y - iy
        sx = 1 - tx
        sy = 1 - ty

        w0 = sx * sy
        w1 = tx * sy
        w2 = tx * ty
        w3 = sx * ty

        i0 = ix + iy * (width + isU)
        i1 = i0 + 1
        i2 = i1 + width + isU
        i3 = i2 - 1

        if particleToGrid:
            if fillGrid:
                grid[i0] += val * w0
                grid[i1] += val * w1
                grid[i2] += val * w2
                grid[i3] += val * w3
            weights[i0] += w0
            weights[i1] += w1
            weights[i2] += w2
            weights[i3] += w3
        else:
            valid0 = waterGrid[i0]
            valid1 = waterGrid[i1]
            valid2 = waterGrid[i2]
            valid3 = waterGrid[i3]

            d = valid0 * w0 + valid1 * w1 + valid2 * w2 + valid3 * w3

            if d > 0:
                return (valid0 * w0 * grid[i0] + 
                        valid1 * w1 * grid[i1] + 
                        valid2 * w2 * grid[i2] + 
                        valid3 * w3 * grid[i3]) / d
            return 0

    def updateDensity(self, width, height):
        density = [0] * self.numCells
        for p in range(self.numParticles):
            self.interp(width, height, self.particleX[p], self.particleY[p], 0, density, 0, 0, 1, 0, False)
        
        waterCount = 0
        densityCount = 0
        for i in range(self.numCells):
            if self.isWater[i]:
                waterCount += 1
                densityCount += density[i]
        
        self.density = density
        return densityCount / waterCount

    def particlesToGrid(self, width, height):
        u = [0] * self.numCellsU
        v = [0] * self.numCellsV
        
        weights = [0] * self.numCellsU
        for p in range(self.numParticles):
            self.interp(width, height, self.particleX[p] + 0.5, self.particleY[p], self.particleU[p], weights, u, 0, True, 1)
        for i in range(self.numCellsU):
            if weights[i] > 0:
                u[i] /= weights[i]
        
        weights = [0] * self.numCellsV
        for p in range(self.numParticles):
            self.interp(width, height, self.particleX[p], self.particleY[p] + 0.5, self.particleV[p], weights, v, 0, True)
        for i in range(self.numCellsV):
            if weights[i] > 0:
                v[i] /= weights[i]

        self.u = u
        self.v = v

    def gridToParticles(self, width, height):
        for p in range(self.numParticles):
            self.particleU[p] = self.interp(width, height, self.particleX[p] + 0.5, self.particleY[p], 0, 0, self.u, self.isWaterU, False, 1)
            self.particleV[p] = self.interp(width, height, self.particleX[p], self.particleY[p] + 0.5, 0, 0, self.v, self.isWaterV, False)

    def calculateTotalDivergence(self):
        totalDivergence = 0
        numWater = 0
        for x in range(self.width):
            for y in range(self.height):
                if (self.notSolid[self.xy(x, y)]):
                    numWater += 1
                    totalDivergence += abs(self.u[self.xy(x+1, y  , 1)] - 
                                           self.u[self.xy(x  , y  , 1)] + 
                                           self.v[self.xy(x  , y+1   )] - 
                                           self.v[self.xy(x  , y     )])
        return totalDivergence / numWater

    def enforceIncompressability(self, averageDensity):
        for _ in range(self.incompressibilityIters):
            #if _ % 5 == 0:
            #    print("Average Divergence in", _, "runs:", self.calculateTotalDivergence())
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
                        
                        if averageDensity > 0:
                            compression = self.density[self.xy(x, y)] - averageDensity
                            if compression > 0:
                                divergence -= 1 * compression

                        divergence = -divergence / s

                        self.u[self.xy(x+1, y  , 1)] += (divergence * sur)
                        self.u[self.xy(x  , y  , 1)] -= (divergence * sul)
                        self.v[self.xy(x  , y+1   )] += (divergence * svd)
                        self.v[self.xy(x  , y     )] -= (divergence * svu)
