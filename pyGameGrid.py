# ...existing code...
import pygame
import numpy as np
import math
import colorsys
import eulerianGrid as eg
import sys
import time

# --- Pygame Setup ---
pygame.init()
WIDTH, HEIGHT = 600, 600
screen = pygame.display.set_mode((WIDTH, HEIGHT))
clock = pygame.time.Clock()

# --- Grid Settings ---
rows, cols = 100, 100
cell_w = WIDTH / cols
cell_h = HEIGHT / rows

# initialize class-level grid in eulerianGrid (matches eg.gridPoint.initializeGrid signature)
eg.gridPoint.initializeGrid(rows, cols)

# --- Velocity Field (vx, vy) (kept for fast drawing) ---
# IMPORTANT vx is y and vy is x (shoutout I setup the grid sideways)
vx = np.zeros((cols, rows), dtype=np.float32)
vy = np.zeros((cols, rows), dtype=np.float32)
vd = np.zeros((cols, rows), dtype=np.float32)

def somethingToColor(x):
    x = max(-60, min(60, x))
    r = int((1 - 2**(-x/1000))*255)
    return r,r,r

def velocityToColor(vx_val, vy_val):
    """
    Convert velocity vector to RGB color based on direction (HSV hue) and magnitude.
    """
    angle = math.atan2(vy_val, vx_val)  # range: [-pi, pi]
    hue = (angle + math.pi) / (2 * math.pi)  # map to [0, 1]
    mag = math.sqrt(vx_val**2 + vy_val**2)
    brightness = min(1.0, mag * 2.0)  # adjust scaling if needed
    r, g, b = colorsys.hsv_to_rgb(hue, 1.0, brightness)
    return int(r * 255), int(g * 255), int(b * 255)

def drawVelocityField():
    for i in range(rows):
        for j in range(cols):
            #color = velocityToColor(vx[i, j], vy[i, j])
            color = somethingToColor(vd[i, j])
            if i == 2 and j == 80:
                color = 100,5,50
            rect = pygame.Rect(int(j * cell_w), int(i * cell_h), int(cell_w) + 1, int(cell_h) + 1)
            pygame.draw.rect(screen, color, rect)

def updateVelocityField():
    # try to advance the physics in the eulerianGrid if available
    eg.gridPoint.nextFrame()

    # copy velocities from the class-level grid into local arrays for drawing
    
    for i in range(rows):
        for j in range(cols):
            v = eg.gridPoint.grid[i][j]
            vx[cols-j-1, i] = float(v.velocity[1])
            vy[cols-j-1, i] = float(v.velocity[0])
            vd[cols-j-1, i] = float(v.density)


# --- Main Loop ---
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
    updateVelocityField()
    drawVelocityField()
    pygame.display.flip()
    time.sleep(0.01)
    print("\n-----------------------------------\n")
    #sys.quit



pygame.quit()
# ...existing code...