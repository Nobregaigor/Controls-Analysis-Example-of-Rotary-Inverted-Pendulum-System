# Import standard modules.
import sys

# Import non-standard modules.
import pygame as pyg
from pygame.locals import *
import numpy as np
import numpy.linalg as la
import scipy.linalg as spla
import matplotlib.pyplot as plt
from pprint import pprint as pp
from classes.controls import controls
from classes.ctr_sys import ctr_sys


class Background(pyg.sprite.Sprite):
    def __init__(self, image_file, location, w=1280,h=720):
        pyg.sprite.Sprite.__init__(self)  #call Sprite initializer
        self.image = pyg.image.load(image_file)
        self.image = pyg.transform.scale(self.image, (w, h))
        self.rect = self.image.get_rect()
        self.location = location
        self.rect.left, self.rect.top = location

    def draw(self,screen,moved=False):
        if moved == True:
            self.rect.left, self.rect.top = self.location
        screen.blit(self.image, self.rect)

class Rocket():
    def __init__(self,x0,y0,theta0):
        self.x = x0
        self.y = y0
        self.theta = theta0

        self.img_h = 135
        self.img_w = 110
        self.img_x = x0 + 0.5 * self.img_w
        self.img_y = y0 + 0.5 * self.img_h

        self.img = Background('rocket_simulation/rocket.png',[self.img_x,self.img_y],self.img_w,self.img_h)

    def update_location(self,x,y):
        self.x = x
        self.y = y
        self.img_x = x + 0.5 * self.img_w
        self.img_y = y + 0.5 * self.img_h
        self.img.location = [self.img_x,self.img_y]
        self.img.rect.left = self.img_x
        self.img.rect.top = self.img_y
        print('location: ' + str(self.img.location))

    def draw(self,screen):
        self.img.draw(screen,True)


def update(dt):
  """
  Update game. Called once per frame.
  dt is the amount of time passed since last frame.
  """

  # Go through events that are passed to the script by the window.
  for event in pyg.event.get():

    # We need to handle these events. Initially the only one you'll want to care
    # about is the QUIT event, because if you don't handle it, your game will crash
    # whenever someone tries to exit.
    if event.type == QUIT:
      pyg.quit() # Opposite of pygame.init
      sys.exit() # Not including this line crashes the script on Windows. Possibly
      # on other operating systems too, but I don't know for sure.
    # Handle other events as you wish.

def draw(screen,BackGround,Platform,Rocket):
  """
  Draw things to the window. Called once per frame.
  """
  screen.fill((0, 0, 0))
  BackGround.draw(screen)
  Platform.draw(screen)
  Rocket.draw(screen)

  pyg.display.flip()

def runPyGame():
  # Initialise PyGame.
  pyg.init()

  # Set up the clock. This will tick every frame and thus maintain a relatively constant framerate. Hopefully.
  fps = 60.0
  fpsClock = pyg.time.Clock()


  # Set up the window.
  width, height = 1280, 720
  screen = pyg.display.set_mode((width, height))

  # screen is the surface representing the window.

  BackGround = Background('rocket_simulation/desert-nights-moon.jpg', [0,0])
  Platform = Background('rocket_simulation/platform.png', [500,720-300],300,100)
  Rkt = Rocket(0,0,0)

  # Main game loop.
  dt = 1/fps # dt is the time since last frame.

  x = 0
  y = 0
  while True: # Loop forever!
    update(dt) # You can update/draw here, I've just moved the code for neatness.
    Rkt.update_location(x,y)
    draw(screen,BackGround,Platform,Rkt)

    x += 1
    y += 1

    dt = fpsClock.tick(fps)

if __name__ == '__main__':
    runPyGame()
