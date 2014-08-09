#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""module mercury qui définit des objets pour correspondant à chaque fichier de paramètre du code mercury (Chambers, 1999). Ces 
objets ne définissent pas de manière plus simple de définir une simulation. C'est simplement un module qui donne la possibilité 
de définir et d'écrire les fichiers de paramètres dans un script python."""
__author__ = "Autiwa <autiwa@gmail.com>"
__date__ = "2012-07-26"
__version__ = "2.2.2"

import os
import pdb # usefull to debug with pdb.set_trace()

AN = 365.25  # nombre de jours dans un an, c'est plus simple ensuite pour calculer T

class MetaClass(object):
  """class that define general methods usefull for inheritance for other class. To sum up, this class doesn't do anything and is defined to be used for inheritance only"""
  
  def __init__(self):
    """ for the moment, do nothing. Nothing to initialize"""
  
  def doc(self, method="class"):
    """method to display the documentation of a class. 
    
    method="class" : if nothing is given, by default, it prints the documentation of the class. Else, the string must be the name of a method
    
    return : nothing. It prints the documentation of the class, by default, or of the method whose name is given in parameter
    """

    if (method == "class"):
      print("Documentation of the Class of this object :",type(self))
      print(self.__doc__)
    elif (method in dir(self)):  # test if a method with the name contained in the var 'method' exist within the current class
      print("Documentation of '"+method+"'")
      print(eval("self."+method+".__doc__"))
    else:
      raise Warning("There is no method named '"+method+"' in this class")

class Body(MetaClass):
  """
  class that define a planet and his characteristics. This is a meta class, used by BodyAst, BodyCom and BodyCart
  
  Optional parameters common to all types of bodies:
  m = X    where X is a real number, to indicate the body's mass in
      Solar masses. If you don't specify a value the mass is
      assumed to be 0.

  r = X    where X is a real number, to indicate the maximum
      distance from the body (in Hill radii) that constitutes
      a close encounter. If you don't include this the default
      is r=1

  d = X    where X is a real number, to indicate the density of the
      body in g/cm^3. If you don't include this the default is d=1

  a1 = X   where X is a real number, to indicate the A1 non-gravitational
      force parameter for this body. Realistically this should be
      zero for Big bodies (the default is 0).

  a2 = X   where X is a real number, to indicate the A2 non-gravitational
      force parameter for this body (the default is 0).

  a3 = X   where X is a real number, to indicate the A1 non-gravitational
      force parameter for this body (the default is 0).
            
  name=None : The name of the planet. If none is given, a default name will be set, equal to PLANETEi where i is the number of instance of the class Planet
  
  method :
  you can display properties of the planet via print name_instance since __str__ have been overload
  you can test if a planet is the same as another via (a_planet == another_planet) (m, a, e and I must be the same in order have True)
  """
  
  NB_SMALL = 0  # The number of small bodies. Usefull for the default name of the planets
  NB_BIG = 0 # The number of big bodies
  BIG_PREFIX = "BIG_" # must take 4 characters
  SMALL_PREFIX = "AST_" # must take 4 characters

  def __init__(self, type, name=None, m=None, r=None, d=None, a1=None, a2=None, a3=None, ep=None):
    """set the optional parameters of the body"""
    
    MetaClass.__init__(self)
    
    self.m, self.r, self.d, self.a1, self.a2, self.a3 = m, r, d, a1, a2, a3
    
    if (type == "big"):
      self.type = type
      Body.NB_BIG += 1
    elif (type == "small"):
      self.type = type
      Body.NB_SMALL += 1
    else:
      raise TypeError("the current type does not exist, you must specify either 'big' or 'small'")
    
    if ((self.type == "small") and (ep != None)):
      self.ep = ep
      self.isEPOCH = True
    else:
      self.ep = None
      self.isEPOCH = False
    
    # By default, we give a name based on the index of the planet (assuming all the instances will be on the planetary system. 
    # That means that if we define two planetary systems, the second one will have planet with indexes that follow indexes of the 
    # first planetary system.
    if (name == None):
      if (self.type == "small"):
        self.name = Body.SMALL_PREFIX+number_fill(Body.NB_SMALL, 4)
      elif (self.type == "big"):
        self.name = Body.BIG_PREFIX+number_fill(Body.NB_BIG, 4)
    else:
      self.name = name
    
  
  def __str__(self):
    """overload the str method. As a consequence, you can print the object via print name_instance
    
    return : a string that represent the properties of the object
    """
    
    texte = ""
          
    texte += self.name
    
    if (self.m != None):
      texte += "\t m="+str(self.m)
    if (self.r != None):
      texte += "\t r="+str(self.r)
    if (self.d != None):
      texte += "\t d="+str(self.d)
    if (self.a1 != None):
      texte += "\t a1="+str(self.a1)
    if (self.a2 != None):
      texte += "\t a2="+str(self.a2)
    if (self.a3 != None):
      texte += "\t a3="+str(self.a3)
    if (self.isEPOCH):
      texte += "\t Ep="+str(self.ep)

    return texte
  
  def format(self):
    """method that return a formatted output usefull to write properties of the Body in a file (either big.in or small.in)
    
    return : a string that represent the properties of the object and that fit requirements for mercury
    """
    
    texte = ""
          
    texte += self.name
    
    if (self.m != None):
      texte += " m="+str(self.m)
    if (self.r != None):
      texte += " r="+str(self.r)
    if (self.d != None):
      texte += " d="+str(self.d)
    if (self.a1 != None):
      texte += " a1="+str(self.a1)
    if (self.a2 != None):
      texte += " a2="+str(self.a2)
    if (self.a3 != None):
      texte += " a3="+str(self.a3)
    if (self.isEPOCH):
      texte += " Ep="+str(self.ep)

    return texte
  
  def get_info(self, line):
    """method that get info from a line formatted by a mercury file (either "small.in" or "big.in"). 
    A lot of arguments will be optional and set to None if nothing is given in the line.
    """
    line = line.rstrip("\n")
    elements = line.split()
    self.name = elements[0]
    
    correspondance = {"m":self.m , "r":self.r , "d":self.d , "a1":self.a1 , "a2":self.a2 , "a3":self.a3 , "ep":self.ep}
    
    for element in elements[1:]:
      (key, value) = element.split("=")
      if (key == 'm'):
        self.m = float(value)
      elif (key == 'r'):
        self.r = float(value)
      elif (key == 'd'):
        self.d = float(value)
      elif (key == 'a1'):
        self.a1 = float(value)
      elif (key == 'a2'):
        self.a2 = float(value)
      elif (key == 'a3'):
        self.a3 = float(value)
      elif (key == 'ep'):
        self.ep = float(value)
    
    if (self.ep != None):
      self.isEPOCH = True
    else:
      self.isEPOCH = False
    
  
  @classmethod
  def resetCounter(cls):
    """method to reset the counters of bodies. 
    
    The class parameters NB_SMALL and NB_BIG will be set to 0. 
    
    WARNING: Problems could occurs if you use the method outside the 'Body' class (in inherited class for instance)"""
    
    cls.NB_SMALL = 0  # The number of small bodies. Usefull for the default name of the planets
    cls.NB_BIG = 0 # The number of big bodies
  
class BodyAst(Body):
  """Class that define a body with asteroïdal coordinates
  
  Parameters:
  a = semi-major axis (in AU)
  e = eccentricity
  I = inclination (degrees)
  g = argument of pericentre (degrees)
  n = longitude of the ascending node (degrees)
  M = mean anomaly (degrees)
  sx=0, sy=0, sz=0 : the 3 components of spin angular momentum for the body,
               in units of solar masses AU^2 per day
  """
  
  def __init__(self, type, a, e, I, g, n, M, sx=0, sy=0, sz=0, name=None, m=None, r=None, d=None, a1=None, a2=None, a3=None, ep=None):
    """Initialisation of the object"""
    
    Body.__init__(self, type=type, name=name, m=m, r=r, d=d, a1=a1, a2=a2, a3=a3, ep=ep)
    
    self.style = "Asteroidal"
    
    self.a, self.e, self.I, self.g, self.n, self.M = a, e, I, g, n, M
    
    self.sx, self.sy, self.sz = sx, sy, sz
  
  def __str__(self):
    """to redefine the str() function, and add information to the one in the parent class 'Body'
    """
    
    texte = Body.__str__(self)
    texte += "\n"
    
    texte += "a = " + str(self.a) +" ua\t" +"e = " + str(self.e) + "\t" +"I = " + str(self.I) +"°\n"
    texte += "g = " + str(self.g) +"°\t" +"n = " + str(self.n) +"°\t" +"M = " + str(self.M) +"°\n"
    texte += "sx = " + str(self.sx) +" ms.au^2/day\t" +"sy = " + str(self.sy) + " ms.au^2/day\t" +"sz = " + str(self.sz) +" ms.au^2/day"
    return texte
  
  def format(self):
    """method that return a formatted output usefull to write properties of the Body in a file (either big.in or small.in)
    
    return : a string that represent the properties of the object and that fit requirements for mercury
    """
    
    texte = Body.format(self)
    texte += "\n"
    
    texte += str(self.a) +" " + str(self.e) + " " + str(self.I) +" " + str(self.g) +" " + str(self.n) +" " + str(self.M)+" " + str(self.sx)+" " + str(self.sy)+" " + str(self.sz)+"\n"
    
    return texte
  
  def get_info(self, lines):
    """method that get info from a line formatted by a mercury file (either "small.in" or "big.in"). 
    A lot of arguments will be optional and set to None if nothing is given in the line.
    """
    Body.get_info(self, lines[0])
    
    line = lines[1].rstrip("\n")
    elements = line.split()
    
    (self.a, self.e, self.I, self.g, self.n, self.M, self.sx, self.sy, self.sz) = map(float, elements)
    
  
    
class BodyCom(Body):
  """Class that define a body with cometary coordinates
  
  Parameters:
  q = periastron (in AU)
  e = eccentricity
  I = inclination (degrees)
  g = argument of pericentre (degrees)
  n = longitude of the ascending node (degrees)
  T = epoch of pericentre (days)
  sx, sy, sz : the 3 components of spin angular momentum for the body,
               in units of solar masses AU^2 per day
  """
  
  def __init__(self, type, q, e, I, g, n, T, sx=0, sy=0, sz=0, name=None, m=None, r=None, d=None, a1=None, a2=None, a3=None, ep=None):
    """Initialisation of the object"""
    
    Body.__init__(self, type=type, name=name, m=m, r=r, d=d, a1=a1, a2=a2, a3=a3, ep=ep)
    
    self.style = "Cometary"
    
    self.q, self.e, self.I, self.g, self.n, self.T = q, e, I, g, n, T
    
    self.sx, self.sy, self.sz = sx, sy, sz
  
  def __str__(self):
    """to redefine the str() function, and add information to the one in the parent class 'Body'
    """
    
    texte = Body.__str__(self)
    texte += "\n"
    texte += "q = " + str(self.q) +" ua\t" +"e = " + str(self.e) + "\t" +"I = " + str(self.I) +"°\n"
    texte += "g = " + str(self.g) +"°\t" +"n = " + str(self.n) +"°\t" +"M = " + str(self.T) +" days\n"
    texte += "sx = " + str(self.sx) +" ms.au^2/day\t" +"sy = " + str(self.sy) + " ms.au^2/day\t" +"sz = " + str(self.sz) +" ms.au^2/day"
    
    return texte
  
  def format(self):
    """method that return a formatted output usefull to write properties of the Body in a file (either big.in or small.in)
    
    return : a string that represent the properties of the object and that fit requirements for mercury
    """
    
    texte = Body.format(self)
    texte += "\n"

    texte += str(self.q) +" " + str(self.e) + " " + str(self.I) +" " + str(self.g) +" " + str(self.n) +" " + str(self.T)+" " + str(self.sx)+" " + str(self.sy)+" " + str(self.sz)+"\n"
    
    return texte
  
  def get_info(self, lines):
    """method that get info from a line formatted by a mercury file (either "small.in" or "big.in"). 
    A lot of arguments will be optional and set to None if nothing is given in the line.
    """
    Body.get_info(self, lines[0])
    
    line = lines[1].rstrip("\n")
    elements = line.split()
    
    (self.q, self.e, self.I, self.g, self.n, self.T, self.sx, self.sy, self.sz) = map(float, elements)


class BodyCart(Body):
  """Class that define a body with cartesian coordinates
  
  Parameters:
  x, y, z : the 3 components of position in AU
  vx, vy, vz : the 3 components of velocity in AU/day
  sx=0, sy=0, sz=0 : the 3 components of spin angular momentum for the body,
               in units of solar masses AU^2 per day
  """
  
  def __init__(self, type, x ,y ,z ,vx ,vy ,vz, sx=0, sy=0, sz=0, name=None, m=None, r=None, d=None, a1=None, a2=None, a3=None, ep=None):
    """Initialisation of the object"""
    
    Body.__init__(self, type=type, name=name, m=m, r=r, d=d, a1=a1, a2=a2, a3=a3, ep=ep)
    
    self.style = "Cartesian"
    
    self.x, self.y, self.z, self.vx, self.vy, self.vz = x, y, z, vx, vy, vz
    
    self.sx, self.sy, self.sz = sx, sy, sz

  def __str__(self):
    """to redefine the str() function, and add information to the one in the parent class 'Body'
    """
    
    texte = Body.__str__(self)
    texte += "\n"

    texte += "x = " + str(self.x) +" ua\t" +"y = " + str(self.y) + " ua\t" +"z = " + str(self.z) +" ua\n"
    texte += "vx = " + str(self.vx) +" ua/day\t" +"vy = " + str(self.vy) + " ua/day\t" +"vz = " + str(self.vz) +" ua/day\n"
    texte += "sx = " + str(self.sx) +" ms.au^2/day\t" +"sy = " + str(self.sy) + " ms.au^2/day\t" +"sz = " + str(self.sz) +" ms.au^2/day"
    
    return texte
  
  def format(self):
    """method that return a formatted output usefull to write properties of the Body in a file (either big.in or small.in)
    
    return : a string that represent the properties of the object and that fit requirements for mercury
    """
    
    texte = Body.format(self)
    texte += "\n"
    
    texte += str(self.x) +" " + str(self.y) + " " + str(self.z) +"\n"
    texte += str(self.vx) +" " + str(self.vy) + " " + str(self.vz) + "\n"
    texte += str(self.sx) +" " + str(self.sy)+" " + str(self.sz)+"\n"
    
    return texte
  
  def get_info(self, lines):
    """method that get info from a line formatted by a mercury file (either "small.in" or "big.in"). 
    A lot of arguments will be optional and set to None if nothing is given in the line.
    """
    Body.get_info(self, lines[0])
    
    line = lines[1].rstrip("\n")
    elements = line.split()
    
    (self.x, self.y, self.z, self.vx, self.vy, self.vz, self.sx, self.sy, self.sz) = map(float, elements)

class PlanetarySystem(MetaClass):
  """
  class that define a planetary system, i.e orbitals elements of each member of the system. For each parameters we can give a list of values that represent the planets or a tuple that give the minimum and maximum values for random generation of the parameters. In addition, you may set the nb_planets parameter. If none is given it will be assumed to be the length of the longest list. If no list is given, you must give in addition the number of planets by setting nb_planets. The nb_planets is defined in init only. By set_X, you can only set list that have the same size. 
  
  Parameters
  planets=[] : a list of object of type Bodyast, BodyCart or BodyCom
  m_star : the mass of the central body in solar mass
  epoch : epoch of start of integration in days
  
  Methods
  __str__ : simply print on the screen the informations available on the planetary system when you do "print name_instance". It overrides the existing __str__ method that is called by the function "print"
  à compléter.
  append : allow to append one or several 'Planet' object(s) to the PlanetarySystem
  set_epoch : allow to change the epoch
  set_m_star : allow to change the mass of the star (in solar masses)
  """

  def __init__(self, bodies=[], m_star=1.0, epoch=0):
    """
    We set the list of orbitals elements as parameters of the instance
    """
    
    MetaClass.__init__(self)
    
    self.small = []
    self.big = []
    
    for body in bodies:
      if body.type == "small":
        self.small.append(body)
      elif body.type == "big":
        self.big.append(body)
      else:
        Warning("Invalid type for:"+str(body))
    
    # Big and small bodies must have respectively only one type 
    # of parameters (in Cartesian, Cometary and Asteroidal), so 
    # we take the style of the first element, it will be tested further out.
    try:
      self.BigStyle = self.big[0].style
    except:
      self.BigStyle = "Asteroidal" # If there is no small bodies, we define a default value for an empty big.in file.
    
    try:
      self.SmallStyle = self.small[0].style
    except:
      self.SmallStyle = "Asteroidal" # If there is no small bodies, we define a default value for an empty small.in file.
    
    if not(self.__nonzero__()):
      Warning("The system is not valid")
    
    self.m_star = m_star
    self.epoch = epoch
  
  def __nonzero__(self):
    """method that test if the planetary system is valid. This method is called when you use bool(instance)
    
    Parameters : None
    
    Return : True if the system is valid. Return False if not. 
    """
    if (len(self.small) != 0):
      isEPOCH = False
      for body in self.small:
        if (body.style != self.SmallStyle):
          Warning("style of this object:"+str(body)+" is different from the expected one ("+self.SmallStyle+")")
          return False
        if body.isEPOCH:
          isEPOCH = True
      
      if isEPOCH:
        ep = self.small[0]
        for body in self.small:
          if (body.ep != ep):
            return False
    
    for body in self.big:
      if (body.style != self.BigStyle):
        return False
      
    return True
  
  def __str__(self):
    """Display the current parameters of the planetary system
    """

    if not(bool(self)):
      print("No valid planetary system is set.")
      return
    
    texte = ""
    if ((len(self.big)) > 0):
      texte += "mass of the central body :"+str(self.m_star)+"ms\n"
      texte += "epoch :"+str(self.epoch)+"\n"
      
      texte += "Big Bodies:\n"
      
      for body in self.big:
        texte += str(body)+"\n"
      
      if (len(self.small) != 0):
        texte += "Small Bodies:\n"
      for body in self.small:
        texte += str(body)+"\n"
    
    return texte
  
  def set_m_star(self, m_star):
    """ set the mass
    
    m_star : the masse of the central body in solar masses
    """
    
    if (m_star < 0):
      raise ValueError("m_star = "+str(m_star)+"\nm_star must be a positive floating number.")
    
    self.m_star = m_star
    
  def set_epoch(self, epoch):
    """ set the epoch
    
    epoch : the epoch of osculation in days (i.e. the time at which the initial coordinates/elements are valid).
    """
    
    self.epoch = epoch
  
  def append(self, *bodies):
    """function that add planets to the planetary system
    
    Parameters
    *bodies : unknown number of arguments that must be of type 'BodyAst', 'BodyCom' or 'BodyCart' 
    """
    # we insert theses planets in beginning of the lists because of the possible random generation that extend the list at the end. To avoid errors, in particuliar on names, we insert them at the beginning (the random generation doesn't touch the beginning of the list)
    for body in bodies:
      if (body.type == "small"):
        self.small.append(body)
      elif (body.type == "big"):
        self.big.append(body)
      else:
        Warning("The type of "+str(body)+" is incorrect")

class Big(object):
  """class that define an object equivalent to big.in, a parameter file of a mercury simulation. 
  If nothing is given, an empty object and planetary system is created.
  
  Parameters:
  system : an object of type 'PlanetarySystem'
  """
  
  BIG_START = ")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n" + \
        ") Lines beginning with `)' are ignored.\n" + \
        ")---------------------------------------------------------------------\n"
  BIG_INT = ")---------------------------------------------------------------------\n"
  
  def __init__(self, system=None):
    
    if (type(system) == PlanetarySystem):
      self.system = system
    elif (system == None):
      self.system = PlanetarySystem()
    else:
      raise TypeError("'system must be a 'PlanetarySystem' object")
      
    
    
  def write(self):
    """write all the data in a file named 'big.in' in the current working directory"""
    
    #########################
    #on crée le fichier et We write the header
    #########################
    bigin = open('big.in','w')

    #on recopie l'entête du fichier type
    bigin.write(Big.BIG_START)
    bigin.write(" style (Cartesian, Asteroidal, Cometary) = "+str(self.system.BigStyle)+"\n")
    bigin.write(" epoch (in days) = "+str(self.system.epoch)+"\n")
    bigin.write(Big.BIG_INT)

    ##########################
    # écriture des paramètres orbitaux des planètes
    ##########################
    for body in self.system.big:
      bigin.write(body.format())

    bigin.close()
    
  def __str__(self):
    """to overwrite the str() method"""
    
    string = Big.BIG_START
    string += " style (Cartesian, Asteroidal, Cometary) = "+str(self.system.BigStyle)+"\n"
    string += " epoch (in days) = "+str(self.system.epoch)+"\n"
    string += Big.BIG_INT
    for body in self.system.big:
      string += body.format()
      
    return string
    
  def read(self):
    """method to read properties from a 'big.in' file in the current working directory
    
    /!\ We cannot set the mass of the star in the planetary system with this method. 
    Since the mass of the star is not used by 'big.in', it is not really important, 
    but it could be if the mass was used somewhere else from the object embedded in the Big object.
    """
    
    bigin = open('big.in','r')
    
    lines = []
    for line in bigin:
      if (line[0] != ")"):
        lines.append(line)
    bigin.close()
    
    (key, value) = lines[0].split("=")
    if (key.count("style") > 0):
      BigStyle = value.split()[0]
    else:
      raise TypeError("'style' parameter expected at this line")
      
    (key, value) = lines[1].split("=")
    if (key.count("epoch") > 0):
      epoch = float(value.split()[0])
    else:
      raise TypeError("'epoch' parameter expected at this line")
    
    # We delete the two first lines
    del(lines[0:2])
    
    if (BigStyle == 'Asteroidal'):
      body_type = BodyAst
    elif (BigStyle == 'Cometary'):
      body_type = BodyCom
    elif (BigStyle == 'Cartesian'):
      body_type = BodyCart
    else:
      raise ValueError("the style of coordinate do not match with anything : "+BigStyle)
    
    # We merge the lines by two because all the informations of one planet are contained in two consecutive lines
    bodies = []
    i = 0
    while (i < len(lines)):
      line1 = lines[i]
      
      i += 1
      # We append lines while the total number of elements is not 9
      nb_el = 0
      els = []
      while (nb_el < 9):
        els.extend(lines[i].split())
        nb_el = len(els)
        i += 1
      line2 = " ".join(els)
      
      #~ pdb.set_trace()
      body = body_type("big", 0, 0, 0, 0, 0, 0)
      body.get_info([line1, line2])
      bodies.append(body)
    
    self.system = PlanetarySystem(bodies=bodies, epoch=epoch)

def readBig():
  """function that return an object "Big" by reading a 'big.in' file in the current working directory
  """
  
  bigin = Big()
  bigin.read()
  return bigin

class Small(object):
  """class that define an object equivalent to small.in, a parameter file of a mercury simulation
  Parameters:
  system : an object of type 'PlanetarySystem'
  """
  
  SMALL_START = ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n" + \
        ") Lines beginning with `)' are ignored.\n" + \
        ")---------------------------------------------------------------------\n"
  SMALL_INT = ")---------------------------------------------------------------------\n"
  
  def __init__(self, system=None):
    
    self.system = system
    
  def write(self):
    """write all the data in a file named 'big.in' in the current working directory"""
    
    #########################
    #on crée le fichier et We write the header
    #########################
    smallin = open('small.in','w')

    #on recopie l'entête du fichier type
    smallin.write(Small.SMALL_START)
    smallin.write(" style (Cartesian, Asteroidal, Cometary) = "+self.system.SmallStyle+"\n")
    smallin.write(Small.SMALL_INT)
    

    ##########################
    # écriture des paramètres orbitaux des planètes/petits corps
    ##########################
    for body in self.system.small:
      smallin.write(body.format())

    smallin.close()

class Element(object):
  """class that define an object equivalent to element.in, a parameter file of a mercury simulation
  
  Parameters
  output_interval = "365.25" : (in days) interval between two outputs (i.e two instants when we write orbital elements)
  format_sortie = " a21e e21e i8.4 g8.4 n8.4 l8.4 m21e "  : variables we want to write and their associated formats (fortran notation)
          The code letters are:
          a = semi-major axis (in AU)
          b = apocentre distance (in AU, b is short for Big q)
          d = density (g per cm^3)
          e = eccentricity
          f = true anomaly (degrees)
          g = argument of perihelion (degrees)
          i = inclination (degrees)
          l = mean anomaly (degrees)
          m = mass (solar masses)
          n = longitude of ascending node
          o = obliquity (degrees)
          p = longitude of perihelion (degrees)
          q = pericentre distance (AU)
          r = radial distance (AU)
          s = spin period (days)
          x, y or z = Cartesian coordinates x, y or z
          u, v or w = Cartesian velocities vx, vy or vz
          Note that if you choose to express an element using a large number
          of significant figures, the last few digits might not be meaningful
          if the output precision of the original integation was low or medium.
  coord = "Cen" : (central body, barycentric, Jacobi) coordonate type we want for speed and position. 
  time_format = "years" : (years, days) in what format is written time
  relative_time = "yes" : (yes/no) is the time during the simulation must be written with respect to the initial time (yes) or not (no) (usefull when the time have a signification, like in our solar system)
  """
  
  ELEMENT_START = ")O+_06 element  (WARNING: Do not delete this line!!)\n" + \
          ") Lines beginning with `)' are ignored.\n" + \
          ")---------------------------------------------------------------------\n" + \
          " number of input files = 1\n" + \
          ")---------------------------------------------------------------------\n" + \
          ") List the input files, one per line\n" + \
          " xv.out\n" + \
          ")---------------------------------------------------------------------\n"
  ELEMENT_INT = ")---------------------------------------------------------------------\n" + \
          ") Output format? (e.g. a8.4 => semi-major axis with 8 digits & 4 dec. places)\n"
  ELEMENT_END = ")---------------------------------------------------------------------\n" + \
          ") Which bodies do you want? (List one per line or leave blank for all bodies)\n" + \
          ")\n"
  
  def __init__(self, output_interval="365.25", format_sortie=" a21e e21e i8.4 g8.4 n8.4 l8.4 m21e ", coord="Cen", time_format="years", relative_time="yes"):
    """initialisation of the class"""
    
    self.output_interval = output_interval
    self.format_sortie = format_sortie
    self.coord = coord
    self.time_format = time_format
    self.relative_time = relative_time
  
  def set_format_sortie(self, format_sortie):
    """explicitely set the format_sortie variable
    
    Parameter
    format_sortie : (in days) interval between two outputs (i.e two instants when we write orbital elements). The default format is :
    " a21e e21e i8.4 g8.4 n8.4 l8.4 m21e "
    """
    
    self.format_sortie = str(format_sortie)
  
  def set_output_interval(self, output_interval):
    """explicitely set the output_interval variable
    
    Parameter
    output_interval : (in days) interval between two outputs (i.e two instants when we write orbital elements)
    """
    
    self.output_interval = str(output_interval)
  
  def set_time_format(self, time_format):
    """explicitely set the time_format variable. 
    
    Parameter
    time_format : sous quel format est écrit le temps (years, days)
    """
    
    self.time_format = time_format
  
  def set_relative_time(self, relative_time):
    """explicitely set the relative_time variable. 
    
    Parameter
    relative_time : (yes/no) est-ce que le temps au cours de la simulation doit être exprimé avec pour référence le début de la simulation (yes) ou pas (dans le cas où la date a une signification dans la simulation)
    """
    self.relative_time = relative_time
  
  def read(self, filename="element.in"):
    """method that allow to get, from a file 'element.in' all the parameters at once
    
    Parameter 
    filename="element.in" : by default, the name is the regular name. But to continue the integration, we must define 
    the file in param.dmp, so we must be able to specify a filename"""
    
    paramin = open(filename, 'r')
    tab = paramin.readlines()
    paramin.close()
    
    parameters = []
    for line in tab:
      if (line[0] != ")"):
        parameters.append(line)
    
    for i in [0, 2, 3, 4, 5]:
        tmp = parameters[i].split("=")[-1]
        tmp = tmp[:-1] # We take off the ending '\n'
        tmp = tmp.split()[0] # we take of extra spaces before and after
        parameters[i] = tmp 
    
    self.coord = parameters[2]
    self.output_interval = float(parameters[3])
    self.time_format = parameters[4]
    self.relative_time = parameters[5]
    self.format_sortie = parameters[6]
    # if perculiar planets only are used, you must add some lines here to get the 
    # last lines of the parameters array that are not currently stored
  
  def write(self):
    """write all the data in a file named 'element.in' in the current working directory"""
    
    ## We generate the file "element.in" with values passed in parameter.s.
    element = open('element.in','w')
    ## We write the header
    element.write(Element.ELEMENT_START)
    ## On écrit dans quel système de coordonnées on veut écrire les données
    element.write(" type of elements (central body, barycentric, Jacobi) = "+str(self.coord)+"\n")
    ## On écrit l'intervalle minimum en jours entre les outputs.
    element.write(" minimum interval between outputs (days) = "+str(self.output_interval)+"\n")
    ## Dans quel format on veut écrire les temps (jours ou années)
    element.write(" express time in days or years = "+str(self.time_format)+"\n")
    ## Est-ce que le temps commence à 0 ("yes") en début de simulation ou pas ("no").
    element.write(" express time relative to integration start time = "+str(self.relative_time)+"\n")
    ## Format de sortie et nom des variables que l'on souhaite avoir.
    # (par défaut, on veut afficher a, e, I, g, n, l, m)
    element.write(Element.ELEMENT_INT)
    element.write(self.format_sortie + "\n")
    ## On écrit la fin du fichier et on le ferme
    element.write(Element.ELEMENT_END)
    element.close()

class Param(object):
  """class that define an object equivalent to param.in, a parameter file of a mercury simulation

  Parameters
  algorithme : ("MVS", "BS", "BS2", "RADAU", "HYBRID")
  start_time : (in days) time when the simulation begin
  stop_time : (in days) time when the simulation stop
  h : timestep in days
  nb_points=None : number of outputs we want. If not set, then data_dump AND output_interval MUST be set
  accuracy = 1.e-12 : precision we want for energy between timesteps
  stop_integration = "no" : (yes/no) if we want to stop integration after a close encounter
  collisions = "yes" : (yes/no) if we want to allow collisions
  fragmentation = "no" : (yes/no) if we allow fragmentation during collision (curently not implemented, so no effects)
  time_format = "years" : sous quel format est écrit le temps (years, days)
  relative_time = "yes" : (yes/no) est-ce que le temps au cours de la simulation doit être exprimé avec pour référence le début de la simulation (yes) ou pas (dans le cas où la date a une signification dans la simulation)
  output_precision = "high" : (low, medium, high) for the precision we want for the outputs (4, 9, 15)
  relativity = "no" : (yes/no) if we want to include relativity (currently not implemented)
  user_force = "no" : (yes/no) if we want to take into account the user-defined force put in the associated sub-routine
  ejection_distance = 1000 : (in au) distance from the star where objets will be treated as ejected (and then erased for the next timestep)
  radius_star = 0.005 : (in au) radius of the central body
  central_mass = 1.0 : (in solar mass) the mass of the central body
  J2 = 0 : J2 moment of the central body in units of its radius^2
  J4 = 0 : J4 moment of the central body in units of its radius^4
  J6 = 0 : J6 moment of the central body in units of its radius^6
  changeover = 3. : (in hill radii) distance from where the HYBRID algorithm will change from MVS to BS
  periodic effect = 100 : The number of timesteps between other periodic effects. At present this controls how often mercury6_1.for checks for ejections and recomputes objects' Hill radii.
  data_dump=None : MUST be set if nb_points is not set
  output_interval=None : MUST be set if nb_points is not set
  
  Remarks : 
  If you want to do the integration back in time (integration from the future to the past), you just have to invert values of start_time and stop_time. Thus, instead of 
    start_stime = 0.
    stop_time = 365.25
  do
    start_stime = 365.25
    stop_time = 0
  leaving the other parameters unchanged.
  """
  
  PARAM_START = ")O+_06 Integration parameters  (WARNING: Do not delete this line!!)\n" + \
          ") Lines beginning with `)' are ignored.\n" + \
          ")---------------------------------------------------------------------\n" + \
          ") Important integration parameters:\n" + \
          ")---------------------------------------------------------------------\n"

  PARAM_INT = ")---------------------------------------------------------------------\n" + \
        ") Integration options:\n" + \
        ")---------------------------------------------------------------------\n"

  PARAM_PAR = ")---------------------------------------------------------------------\n" + \
        ") These parameters do not need to be adjusted often:\n" + \
        ")---------------------------------------------------------------------\n"
  
  def __init__(self, algorithme, start_time, stop_time, h, accuracy=1.e-12, 
  stop_integration="no", collisions="yes", fragmentation="no", time_format="years", 
  relative_time="yes", output_precision="high", relativity="no", user_force="no", 
  ejection_distance=1000, radius_star=0.005, central_mass=1.0, J2=0, J4=0, J6=0, 
  changeover=3., periodic_effect=100, data_dump=500, output_interval=365.25):
    """initialise the class and store the datas
    
    time_format and relative_time are common information with element and param. Do something to combine thoses data
    """
    
    self.algorithme = algorithme
    self.start_time, self.stop_time = start_time, stop_time
    self.output_interval = output_interval
    self.h = h
    self.accuracy = accuracy
    self.stop_integration = stop_integration
    self.collisions = collisions
    self.fragmentation = fragmentation
    self.time_format = time_format
    self.relative_time = relative_time
    self.output_precision = output_precision
    self.relativity = relativity
    self.user_force = user_force
    self.ejection_distance = ejection_distance
    self.radius_star = radius_star
    self.central_mass = central_mass
    self.J2, self.J4, self.J6 = J2, J4, J6
    self.changeover = changeover
    self.data_dump = data_dump
    self.periodic_effect = periodic_effect
    
  
  @property
  def integration_time(self):
    """return the integration time, regards to the start and stop time"""
    
    return self.stop_time - self.start_time
  
  def get_output_interval(self):
    """return the value contained in 'self.output_interval'
    
    Return : output_interval
    """
    
    return self.output_interval
  
  def set_time_format(self, time_format):
    """explicitely set the time_format variable. 
    
    Parameter
    time_format : in what format is written the time (years, days)
    """
    
    self.time_format = time_format
  
  def set_stop_time(self, stop_time):
    """explicitely set the stop_time variable. 
    
    Parameter
    stop_time : stop time, in days
    """
    
    self.stop_time = stop_time
  
  def set_algorithme(self,algorithm):
    """set the algorithm to the value given in parameter"""
    
    self.algorithme = algorithm
    
  def get_start_time(self):
    """return the value contained in 'self.start_time'
    
    Return : start_time (in days)
    """
    
    return self.start_time
  
  def get_stop_time(self):
    """return the value contained in 'self.stop_time'
    
    Return : stop_time (in days)
    """
    
    return self.stop_time
  
  def get_time_format(self):
    """return the value contained in 'self.time_format'
    
    Return : time_format
    """
    
    return self.time_format
  
  def set_relative_time(self, relative_time):
    """explicitely set the relative_time variable. 
    
    Parameter
    relative_time : (yes/no) est-ce que le temps au cours de la simulation doit être exprimé avec pour référence le début de la simulation (yes) ou pas (dans le cas où la date a une signification dans la simulation)
    """
    self.relative_time = relative_time
  
  def get_relative_time(self):
    """return the value contained in 'self.relative_time'
    
    Return : relative_time
    """
    
    return self.relative_time
  
  def read(self, filename="param.in"):
    """method that allow to get, from a file 'param.in' all the parameters at once
    
    Parameter 
    filename="param.in" : by default, the name is the regular name. But to continue the integration, we must define 
    the file in param.dmp, so we must be able to specify a filename"""
    
    paramin = open(filename, 'r')
    tab = paramin.readlines()
    paramin.close()
    
    parameters = []
    for line in tab:
      if (line[0] != ")"):
        tmp = line.split("=")[-1]
        tmp = tmp[:-1] # We take off the ending '\n'
        tmp = tmp.split()[0] # we take of extra spaces before and after
        parameters.append(tmp) # we take the last element [-1] of the split, but without the last caracter (which is '\n')
    
    self.algorithme = parameters[0]
    self.start_time = float(parameters[1])
    self.stop_time = float(parameters[2])
    self.output_interval = float(parameters[3])
    self.h = float(parameters[4])
    self.accuracy = float(parameters[5])
    self.stop_integration = parameters[6]
    self.collisions = parameters[7]
    self.fragmentation = parameters[8]
    self.time_format = parameters[9]
    self.relative_time = parameters[10]
    self.output_precision = parameters[11]
    # Not used
    self.relativity = parameters[13]
    self.user_force = parameters[14]
    self.ejection_distance = float(parameters[15])
    self.radius_star = float(parameters[16])
    self.central_mass = float(parameters[17])
    self.J2 = float(parameters[18])
    self.J4 = float(parameters[19])
    self.J6 = float(parameters[20])
    # Not used
    # Not used
    self.changeover = float(parameters[23])
    self.data_dump = int(parameters[24])
    self.periodic_effect = int(parameters[25])
    
  
  def write(self, filename="param.in"):
    """write all the data in a file named 'param.in' in the current working directory
    
    Parameter
    filename="param.in" : by default, the name is the regular name. But to continue the integration, we must define the file in param.dmp, so we must be able to specify a filename"""
    
    param = open(filename,'w')
    param.write(Param.PARAM_START)
    param.write(" algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = "+self.algorithme+"\n")
    param.write(" start time (days) = "+str(self.start_time)+"\n")
    param.write(" stop time (days) = %e\n" % self.stop_time)
    param.write(" output interval (days) = %e\n" % self.output_interval)
    param.write(" timestep (days) = "+str(self.h)+"\n")
    param.write(" accuracy parameter = "+str(self.accuracy)+"\n")
    param.write(Param.PARAM_INT)

    param.write(" stop integration after a close encounter = "+self.stop_integration+"\n")
    param.write(" allow collisions to occur = "+self.collisions+"\n")
    param.write(" include collisional fragmentation = "+self.fragmentation+"\n")
    param.write(" express time in days or years = "+self.time_format+"\n")
    param.write(" express time relative to integration start time = "+self.relative_time+"\n")
    param.write(" output precision = "+self.output_precision+"\n")
    param.write(" < not used at present >\n")
    param.write(" include relativity in integration = "+self.relativity+"\n")
    param.write(" include user-defined force = "+self.user_force+"\n")
    param.write(Param.PARAM_PAR)
    param.write(" ejection distance (AU) = "+str(self.ejection_distance)+"\n")
    param.write(" radius of central body (AU) = "+str(self.radius_star)+"\n")
    param.write(" central mass (solar) = "+str(self.central_mass)+"\n")
    param.write(" central J2 = "+str(self.J2)+"\n")
    param.write(" central J4 = "+str(self.J4)+"\n")
    param.write(" central J6 = "+str(self.J6)+"\n")
    param.write(" < not used at present >\n")
    param.write(" < not used at present >\n")
    param.write(" Hybrid integrator changeover (Hill radii) = "+str(self.changeover)+"\n")
    param.write(" number of timesteps between data dumps = "+str(self.data_dump)+"\n")
    param.write(" number of timesteps between periodic effects = "+str(self.periodic_effect)+"\n")
    param.close()

class Close(object):
  """class that define an object equivalent to element.in, a parameter file of a mercury simulation
  
  Parameters:
  time_format = "years" : (years, days) in what format is written time
  relative_time = "yes" : (yes/no) is the time during the simulation must be written with respect to the initial time (yes) or not (no) (usefull when the time have a signification, like in our solar system)
  """
  
  CLOSE_START = ")O+_06 close  (WARNING: Do not delete this line!!)\n" + \
          ") Lines beginning with `)' are ignored.\n" + \
          ")---------------------------------------------------------------------\n" + \
          " number of input files = 1\n" + \
          ")---------------------------------------------------------------------\n" + \
          ") List the input files, one per line\n" + \
          " ce.out\n" + \
          ")---------------------------------------------------------------------\n"
  CLOSE_END = ")---------------------------------------------------------------------\n" + \
          ") Which bodies do you want? (List one per line or leave blank for all bodies)\n" + \
          ")\n"
  
  def __init__(self, time_format="years", relative_time="yes"):
    """initialisation of the class"""

    self.time_format = time_format
    self.relative_time = relative_time
  
  
  def write(self):
    """write all the data in a file named 'element.in' in the current working directory"""
    
    ## We generate the file "element.in" with values passed in parameter.s.
    close = open('close.in','w')
    ## We write the header
    close.write(Close.CLOSE_START)
    ## Dans quel format on veut écrire les temps (jours ou années)
    close.write(" express time in days or years = "+str(self.time_format)+"\n")
    ## Est-ce que le temps commence à 0 ("yes") en début de simulation ou pas ("no").
    close.write(" express time relative to integration start time = "+str(self.relative_time)+"\n")
    ## On écrit la fin du fichier et on le ferme
    close.write(Close.CLOSE_END)
    close.close()

class Disk(object):
  """class that define an object equivalent to disk.in, a parameter file of a mercury simulation (implemented by my user_module)
  
  Parameters:
  b_over_h=None : the smoothing length for the planet's potential [No dim]
  adiabatic_index=None : the adiabatic index for the gas equation of state [No dim]
  mean_molecular_weight=None : the mean molecular weight in mass of a proton [m_proton]
  surface_density=None : a tuple of 2 values. The first must be the surface density at R=1AU (in g/cm^2). The second is the negative slope of the surface density power law
  temperature=None : a tuple of 2 values. (radius_min, radius_max), in AU, boundary values for the radius sample of the temperature profile. The number of point is given directly in mfo_user
  viscosity=None : constant viscosity of the disk [cm^2/s]
  """
  
  TORQUE_TYPES = ['real', 'mass_dependant', 'linear_indep', 'tanh_indep', 'manual']
  
  DAMPING_TYPES = ['cossou', 'pierens', 'fendyke', 'none']
  
  OPACITY_TYPES = ['bell', 'zhu', 'chambers', 'hure']
  
  VISCOSITY_TYPES = ['constant', 'alpha', 'alpha_dz']
  
  BOUNDARIES = ["open", "closed"]
  
  DISK_START = "!# ------------------------------------------------\n" + \
          "!# Parameter file for various properties of the disk. \n" + \
          "!# ------------------------------------------------\n" + \
          "!# blanck line or with spaces will be skipped. \n" + \
          "!# In fact, the only lines that matter are non commented lines with a\n" + \
          "!# '=' character to distinguish the identificator and the value(s)\n" + \
          "!# (each value must be separated with at least one space. \n" + \
          "!# Line must not be longer than 80 character, but comments can be far\n" + \
          "!# bigger than that, even on line with a parameter to read.\n\n"
  
  DISK_COMMENT = {'b/h':"!# the smoothing parameter for the planet's potential  \n" +\
                        "!#  two values : first for Lindblad and second for corotation torque \n" +\
                        "!#  one value  : apply for both", 
            'adiabatic_index':"!# the adiabatic index for the gas equation of state", 
            'mean_molecular_weight':"!# the mean molecular weight in mass of a proton", 
            'surface_density':"!# Here we define the power law for surface density \n" +\
                                "!#  sigma(R) = sigma_0 * R^(-sigma_index) \n" + \
                                "!#  where sigma_0 is the surface density at (R=1AU) [g/cm^2] and sigma_index\n" +\
                                "!#  is the negative slope of the surface density power law (alpha in the paper)\n" +\
                                "!# If 'manual' is specified, then the surface density profile will instead be\n" +\
                                "!#  read from 'surface_density_profile.dat', two columns, the first being \n" +\
                                "!#  orbital distance and the second the surface density in g/cm^2", 
            'is_irradiation':"!# (0, False) if there is no stellar irradiation to compute temperature profile\n" + \
                             "!, (1, True) if there is stellar irradiation", 
            'r_star':"!# in solar radii, the radius of the star for irradiation",
            't_star':"!# in K, the temperature of the star for irradiation",
            'disk_albedo':"!# [0..1] the disk albedo for irradiation",
            'disk_edges':"!# Here we define the radius_min and radius_max for the radius sample of the disk \n" +\
                         "!# (used for temperature profile for instance)", 
            'inner_smoothing_width':"!# The width (in unit of the inner boundary radius) of the region right after\n" + \
                                    "!# the inner edge where the surface density is damped so\n" + \
                                    "!# that the surface density at the inner edge will be 0",
            'viscosity_type':"!# %s define the viscosity\n" % VISCOSITY_TYPES +\
                             "!# constant : constant viscosity (being defined with the 'viscosity' parameter)\n"+\
                             "!# alpha : alpha prescription, (alpha being defined with the 'alpha' parameter)\n"+\
                             "!# alpha_dz : (alphas are 3 values stored in 'alpha_dz', and separation radius \n"+\
                             "!#            are 2 values stored in radius_dz)", 
            'viscosity':"!# Constant viscosity of the disk [cm^2/s]", 
            'alpha':"!# Alpha value used for alpha prescription of the viscosity (adimensioned)",
            'alpha_dz':"!# the alpha value for the alpha-prescription in the 3 regions",
            'radius_dz':"!# the two radius that separate the 3 different alpha-regions",
            'opacity_type':"!# %s define the torque type.\n" % OPACITY_TYPES +\
                          "!# bell : from bell & lin 1994\n" +\
                          "!# chambers : from chambers 2009\n" +\
                          "!# hure : opacity table from (hure, 2000)\n" + \
                          "!# zhu : From zhu & hartmann 2009",
            'is_turbulence':"!# (0, False) if there is no turbulence, (1, True) if there is turbulence", 
            'turbulent_forcing':"!# The value of the adimensioned parameter that control the strength of the resonance. \n" +\
                                "!# If not specified, an auto value, based on the value of the viscosity is used.",
            'dissipation_type':"!# integer to tell if there is dissipation of the disk or not. \n" +\
                                             "!# 0 for no dissipation\n" +\
                                             "!# (1) for viscous dissipation (Not implemented anymore)\n" +\
                                             "!# 2 for exponential decay of the initial profile\n" +\
                                             "!# 3 for mixed exponential decay (viscous then photoevap)", 
            'tau_viscous':"!# (years) the exponential decay for the viscous phase (dissipation_type = 3)", 
            'tau_photoevap':"!# (years) exponential decay for photoevaporation phase (dissipation_type = 3)", 
            'dissipation_time_switch':'! (years) the time when we switch from viscous to photoevap (dissipation_type = 3)',
            'disk_exponential_decay':'! Value of the exponential decay timescale for the dissipation of the disk\n' +\
                                     '! (only if dissipation_type = 2)',
            'sample':"!# number of point to the 1D radial grid of the disk"}
  INTERACTIONS_COMMENT = {'damping_type':"!# %s define the corotation damping type.\n" % DAMPING_TYPES +\
                          "!# cossou : from (cossou & raymond, 2013)\n" +\
                          "!# pierens : from (pierens & cossou, 2013)\n" +\
                          "!# fendyke : from (fendyke & nelson, 2013)\n" +\
                          "!# none : no corotation damping",
  'torque_type':"!# %s define the torque type.\n" % TORQUE_TYPES +\
                          "!# real : The torque from (pardekooper et al., 2011)\n" +\
                          "!# linear_indep : Mass independant convergence zone with a linear torque profile\n" +\
                          "!# tanh_indep : Mass independant convergence zone with a \n" +\
                          "!#              tanh torque profile that saturate at a given value\n" +\
                          "!# mass_dependant : Mass dependant convergence zone where for a \n" +\
                          "!#                  given mass the torque profile is linear with distance\n" +\
                          "!# manual : the code will read the file 'torque_profile.dat' that must contain 2\n" +\
                          "!# columns, the first being the semi major axis in AU, and the second the torque",
            'torque_profile_steepness':"!# Gamma = a * x + b. Here is the steeness 'a' of the linear torque profile, both mass-(in)dependant", 
            'saturation_torque':"!# the assymptot for the arctan mass indep convergence zone", 
            'indep_cz':"!# The position of the convergence zone in the 'mass_independant' torque case", 
            'mass_dep_m_min':"!# lower mass for the 'mass_dependant' convergence zone", 
            'mass_dep_m_max':"!# top mass for the 'mass_dependant' convergence zone", 
            'mass_dep_cz_m_min':"!# position of the convergence zone for the lower mass ('mass_dependant' case)", 
            'mass_dep_cz_m_max':"!# position of the convergence zone for the top mass ('mass_dependant' case)"}
  
  DEFAULT = {'b/h':0.4, 
          'adiabatic_index':1.4, 
          'mean_molecular_weight':2.35, 
          'surface_density':(500,0.5), 
          'is_irradiation':1, 
          'r_star':2.5,
          't_star':4000.,
          'disk_albedo':0.5,
          'disk_edges':(0.1, 100.), 
          'inner_smoothing_width':0.05,
          'viscosity_type':'constant', 
          'viscosity':1e15,
          'alpha':1e-3,
          'alpha_dz':(1e-2, 1e-4, 1e-2),
          'radius_dz':(1., 10.),
          'opacity_type':"hure",
          'opacity_type':"cossou",
          'is_turbulence':0, 
          'turbulent_forcing':1.3e-4,
          'dissipation_type':0, 
          'tau_viscous':1e7, 
          'tau_photoevap':3e4, 
          'dissipation_time_switch':2e6,
          'disk_exponential_decay':1e6,
          'sample':400,
          'torque_type':"real",
          'torque_profile_steepness':1., 
          'saturation_torque':1., 
          'indep_cz':3.0, 
          'mass_dep_m_min':1, 
          'mass_dep_m_max':60, 
          'mass_dep_cz_m_min':4, 
          'mass_dep_cz_m_max':30}
  
  def __init__(self, b_over_h=None, adiabatic_index=None, mean_molecular_weight=None, surface_density=None, disk_edges=None, 
               viscosity_type=None, viscosity=None, alpha=None, alpha_dz=None, radius_dz=None, damping_type=None,
               is_turbulence=None, turbulent_forcing=None, inner_smoothing_width=None, tau_viscous=None, tau_photoevap=None, 
               dissipation_time_switch=None, is_irradiation=None, r_star=None, t_star=None, disk_albedo=None, opacity_type=None,
               sample=None, dissipation_type=None, disk_exponential_decay=None, 
               torque_type=None, torque_profile_steepness=None, saturation_torque=None, indep_cz=None, 
               mass_dep_m_min=None, mass_dep_m_max=None, mass_dep_cz_m_min=None, mass_dep_cz_m_max=None):
    """initialisation of the class"""

    self.disk_parameter = {}
    self.interaction_parameter = {}
    
    # either a float or a tuple of two floats
    self.disk_parameter['b/h'] = b_over_h
    
    self.disk_parameter['adiabatic_index'] = adiabatic_index
    
    self.disk_parameter['mean_molecular_weight'] = mean_molecular_weight
    
    self.disk_parameter['surface_density'] = surface_density
    
    if (is_irradiation != None):
      self.disk_parameter['is_irradiation'] = int(is_irradiation)
    else:
      self.disk_parameter['is_irradiation'] = None
    
    self.disk_parameter['r_star'] = r_star
    self.disk_parameter['t_star'] = t_star
    self.disk_parameter['disk_albedo'] = disk_albedo
    
    self.disk_parameter['disk_edges'] = disk_edges
    
    self.disk_parameter['inner_smoothing_width'] = inner_smoothing_width
    
    self.disk_parameter['tau_viscous'] = tau_viscous
    self.disk_parameter['tau_photoevap'] = tau_photoevap
    self.disk_parameter['dissipation_time_switch'] = dissipation_time_switch
    
    if (viscosity_type != None):
      if (viscosity_type in Disk.VISCOSITY_TYPES):
        self.disk_parameter['viscosity_type'] = viscosity_type
      else:
        raise NameError("'viscosity_type' does not correspond to an existing type : '%s'.\nAvailable viscosity_type types : %s" % (viscosity_type, Disk.VISCOSITY_TYPES))
    
    self.disk_parameter['viscosity'] = viscosity
    self.disk_parameter['alpha'] = alpha
    self.disk_parameter['alpha_dz'] = alpha_dz
    self.disk_parameter['radius_dz'] = radius_dz
    
    if (opacity_type != None):
      if (opacity_type in Disk.OPACITY_TYPES):
        self.disk_parameter['opacity_type'] = opacity_type
      else:
        raise NameError("'opacity_type' does not correspond to an existing type : '%s'.\nAvailable opacity_type types : %s" % (opacity_type, Disk.OPACITY_TYPES))
    
    if (is_turbulence != None):
      self.disk_parameter['is_turbulence'] = int(is_turbulence)
    else:
      self.disk_parameter['is_turbulence'] = None
    self.disk_parameter['turbulent_forcing'] = turbulent_forcing
    
    self.disk_parameter['sample'] = sample
    
    self.disk_parameter['dissipation_type'] = dissipation_type
    
    self.disk_parameter['disk_exponential_decay'] = disk_exponential_decay
    
    #####################################################################
    
    if (damping_type != None):
      if (damping_type in Disk.DAMPING_TYPES):
        self.interaction_parameter['damping_type'] = damping_type
      else:
        raise NameError("'damping_type' does not correspond to an existing type : '%s'.\nAvailable damping types : %s" % (damping_type, Disk.DAMPING_TYPES))
      
    if (torque_type != None):
      if (torque_type in Disk.TORQUE_TYPES):
        self.interaction_parameter['torque_type'] = torque_type
      else:
        raise NameError("'torque_type' does not correspond to an existing type : '%s'.\nAvailable torque types : %s" % (torque_type, Disk.TORQUE_TYPES))
      
    self.interaction_parameter['torque_profile_steepness'] = torque_profile_steepness
    
    self.interaction_parameter['saturation_torque'] = saturation_torque
    
    self.interaction_parameter['indep_cz' ] = indep_cz 
    
    self.interaction_parameter['mass_dep_m_min' ] = mass_dep_m_min 
    self.interaction_parameter['mass_dep_m_max' ] = mass_dep_m_max 
      
    self.interaction_parameter['mass_dep_cz_m_min' ] = mass_dep_cz_m_min 
    self.interaction_parameter['mass_dep_cz_m_max' ] = mass_dep_cz_m_max 
    
    
  def demo(self):
    """will create and write a demo file"""
    
    self.disk_parameter['b/h'] = 0.4
  
    self.disk_parameter['adiabatic_index'] = 1.4
  
    self.disk_parameter['mean_molecular_weight'] = 2.35
  
    self.disk_parameter['surface_density'] = (500, 0.5)
  
    self.disk_parameter['is_irradiation'] = 1
  
    self.disk_parameter['disk_edges'] = (0.1, 100)
  
    self.disk_parameter['inner_smoothing_width'] = 0.05
  
    self.disk_parameter['tau_viscous'] = None
    self.disk_parameter['tau_photoevap'] = None
    self.disk_parameter['dissipation_time_switch'] = None
  
    self.disk_parameter['viscosity_type'] = 'constant'
    self.disk_parameter['viscosity'] = 1e15
    self.disk_parameter['alpha'] = None

    self.disk_parameter['opacity_type'] = 'bell'
    
    self.disk_parameter['is_turbulence'] = 0
    self.disk_parameter['turbulent_forcing'] = None
  
    self.disk_parameter['sample'] = 400
  
    self.disk_parameter['dissipation_type'] = 0
  
    self.disk_parameter['disk_exponential_decay'] = None
  
  #####################################################################
    self.interaction_parameter['torque_type'] = 'real'
    
    self.interaction_parameter['torque_profile_steepness'] = None
  
    self.interaction_parameter['saturation_torque'] = None
  
    self.interaction_parameter['indep_cz' ] = None
  
    self.interaction_parameter['mass_dep_m_min' ] = None
    self.interaction_parameter['mass_dep_m_max' ] = None
    
    self.interaction_parameter['mass_dep_cz_m_min' ] = None
    self.interaction_parameter['mass_dep_cz_m_max' ] = None 
    
    self.write()
  
  def write(self):
    """write all the data in a file named 'element.in' in the current working directory"""
    
    ## We generate the file "disk.in" with values passed in parameter.
    disk = open('disk.in','w')
    ## We write the header
    disk.write(Disk.DISK_START)
    
    disk.write("!*****************************\n")
    disk.write("!*      Disk Parameters      *\n")
    disk.write("!*****************************\n")
    for (key, value) in sorted(self.disk_parameter.items()):
      disk.write("\n")
      disk.write(Disk.DISK_COMMENT[key]+"\n") # the comment character is directly in the COMMENT element because some elements might be multilines.
      if (type(value) in (list, tuple)):
        disk.write(key+" = "+" ".join(map(str, value))+"\n")
      elif (value == None):
        disk.write("!%s = %s\n" % (key, Disk.DEFAULT[key]))
      else:
        disk.write(key+" = "+str(value)+"\n")
    
    disk.write("!*****************************\n")
    disk.write("!* Interactions disk/planets *\n")
    disk.write("!*****************************\n")
    for (key, value) in sorted(self.interaction_parameter.items()):
      disk.write("\n")
      disk.write(Disk.INTERACTIONS_COMMENT[key]+"\n") # the comment character is directly in the COMMENT element because some elements might be multilines.
      if (type(value) in (list, tuple)):
        disk.write(key+" = "+" ".join(map(str, value))+"\n")
      elif (value == None):
        disk.write("!%s = %s\n" % (key, Disk.DEFAULT[key]))
      else:
        disk.write(key+" = "+str(value)+"\n")
    
    disk.close()


class Message(object):
  """class that creates and write the message.in parameter file, needed by element6,mercury and maybe close"""
  
  FILE = """  1  6  days 
  2  6  years
  3 13  solar masses
  4  3  AU
  5  3 no 
  6  3 yes
  7  3 low
  8  6 medium
  9  4 high
 10  0 
 11 33            Integration parameters
 12 33            ----------------------
 13 14    Algorithm:
 14 38 Second-order mixed-variable symplectic
 15 24 Bulirsch-Stoer (general)
 16 37 Bulirsch-Stoer (conservative systems)
 17 16 15th-order RADAU
 18  0 
 19  0
 20  0 
 21  0 
 22  5 Test
 23 48 Hybrid symplectic integrator (mixed coordinates)
 24 44 Hybrid symplectic (close binary coordinates)
 25 43 Hybrid symplectic (wide binary coordinates)
 26 32    Integration start epoch:
 27 32    Integration stop  epoch:
 28 32    Output interval:
 29 32    Element origin:
 30 31    Initial timestep:
 31 36    Accuracy parameter:
 32 36    Central mass:
 33 36    J_2:
 34 36    J_4:
 35 36    J_6:
 36 36    Ejection distance:
 37 36    Radius of central body:
 38 29    Number of Big bodies:
 39 29    Number of Small bodies:
 40 37    Output precision: 
 41 40    Includes collisions:
 42 40    Includes fragmentation: 
 43  0 
 44  0 
 45 40    Includes relativity: 
 46 40    Includes user-defined force routine: 
 47 10 barycentre 
 48 12 central body
 49  0 
 50  0 
 51 30            Integration details
 52 30            -------------------
 53 29    Initial energy:
 54 29    Initial angular momentum:
 55 65    Integrating massive bodies and particles up to the same epoch.
 56 34    Beginning the main integration.
 57 24    Integration complete.
 58 48    Fractional energy change due to integrator: 
 59 48    Fractional angular momentum change:
 60 57    Fractional energy change due to collisions/ejections: 
 61 57    Fractional angular momentum change:
 62 47    Continuing integration from dump files at 
 63  6 Time: 
 64  6 Date: 
 65  9    dE/E: 
 66  9    dL/L: 
 67 35  collided with the central body at 
 68 12  ejected at 
 69 12  was hit by 
 70 34  removed due to an encounter with 
 71  4  at 
 72 26  solar masses AU^2 day^-2
 73 26  solar masses AU^2 day^-1
 74 36  lost mass due to rotational breakup
 75 24  removed due to small a
 76  0 
 77  0 
 78  0 
 79  0 
 80  0 
 81  8  ERROR:
 82 49        Modify mercury.inc and recompile Mercury.
 83 62        Check the file containing initial data for Big bodies.
 84 64        Check the file containing initial data for Small bodies.
 85 57        Check the file containing integration parameters.
 86 22        Check files.in
 87 27 This file already exists:  
 88 34 This file is needed to continue:  
 89 30 This filename is duplicated: 
 90 40 The total number of bodies exceeds NMAX.
 91 68 Data style on first line must be Cartesian, Asteroidal or Cometary
 92 68 You cannot integrate non-gravitational forces using this algorithm.
 93 64 You cannot integrate a user-defined force using this algorithm.
 94 64 You cannot integrate massive Small bodies using this algorithm.
 95 66 Massive Small bodies must have the same epoch as the Big bodies.
 96 49 Check character implies input file is corrupted.
 97 62 Mass, density, encounter limit must be >= 0 for this object:
 98 46 This integration algorithm is not available: 
 99 50 A problem occurred reading the parameter on line
100 50 A problem occurred reading data for this object: 
101 56 A problem occured reading the epoch for the Big bodies.
102 67 You cannot use non-zero J2,J4,J6 using the close-binary algorithm.
103 34 Two objects both have this name: 
104 36         is corrupted at line number: 
105 42 Central-body radius exceeds maximum radius. 
106 68 Maximum/Central radius is large. Output precision will be degraded. 
107 58 Coordinate origin must be Central, Barycentric or Jacobi.
108  0 
109  0 
110  0 
111  0 
112  0 
113  0 
114  0 
115  0 
116  0 
117  0 
118  0 
119  0 
120  0 
121 10  WARNING:
122 53 Truncating the name of this object to 8 characters: 
123 30 Main integration is backwards.
124 26 No Big bodies are present.
125 28 No Small bodies are present.
126 50 Stopping integration due to an encounter between 
127 45 Throwing this object into the central body: 
128 42 Setting output threshhold DA to infinity.
129 42 Setting output threshhold DE to infinity.
130 42 Setting output threshhold DI to infinity.
131 43 Increasing the radius of the central body.
132 56 Total number of current close encounters exceeds CMAX.
133  0 
134  0 
135  0 
136  0 
137  0 
138  0 
139  0 
140  0 
141  0 
142  0 
143  0 
144  0 
145  0 
146  0 
147  0 
148  0 
149  0 
150  0 
151 67 )O+_05 Integration parameters  (WARNING: Do not delete this line!!)
152 66 )O+_05 Big-body initial data  (WARNING: Do not delete this line!!)
153 68 )O+_05 Small-body initial data  (WARNING: Do not delete this line!!)
154 39 ) Lines beginning with `)' are ignored.
155 70 )---------------------------------------------------------------------
156 43  style (Cartesian, Asteroidal, Cometary) = 
157 20  epoch (in days) = 
158 35 ) Important integration parameters:
159 48  algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = 
160 21  start time (days) = 
161 20  stop time (days) = 
162 26  output interval (days) = 
163 19  timestep (days) = 
164 22  accuracy parameter = 
165 22 ) Integration options:
166 44  stop integration after a close encounter = 
167 29  allow collisions to occur = 
168 37  include collisional fragmentation = 
169 33  express time in days or years = 
170 51  express time relative to integration start time = 
171 20  output precision = 
172 24  < Not used at present > 
173 37  include relativity in integration = 
174 30  include user-defined force = 
175 52 ) These parameters do not need to be adjusted often:
176 26  ejection distance (AU) = 
177 31  radius of central body (AU) = 
178 31  central mass (solar masses) = 
179 14  central J2 = 
180 14  central J4 = 
181 14  central J6 = 
182 24  < Not used at present > 
183 24  < Not used at present > 
184 45  Hybrid integrator changeover (Hill radii) = 
185 42  number of timesteps between data dumps = 
186 48  number of timesteps between periodic effects = 
187 41  origin (Central, Barycentric, Jacobi) = 
188  0 
189  0 
190  0 
191  0 
192  0 
193  0 
194  0 
195  0 
196  0 
197  0 
198  0 
199  0 
200  0 """
  def __init__(self):
    """Initialisation of the class"""
    
    #Nothing to initialize for the moment
  
  def write(self):
    """Create the file 'message.in'"""
    
    fichier = open("message.in", "w")
    fichier.write(Message.FILE)
    fichier.close()

class Files(object):
  """function that creates and write the files.in parameter file, needed by mercury6"""
  
  FILE = """ big.in
 small.in
 param.in
 xv.out
 ce.out
 info.out
 big.dmp
 small.dmp
 param.dmp
 restart.dmp
 """
  
  def __init__(self):
    """Initialisation of the class"""
    
    # Nothing to initialize for the moment
  
  def write(self):
    """Create the file 'files.in'"""
    fichier = open("files.in", "w")
    fichier.write(Files.FILE)
    fichier.close()

# If we launch the module as a script, we run a series of tests on the defined objects
if __name__=='__main__':
  test = BodyAst("small", 1, 2, 3, 4, 5, 6)
  print(test)
  
  mercury = BodyCart("big",name="MERCURY", x=-3.83966017419175965E-01, y=-1.76865300855700736E-01, z=2.07959213998758705E-02,
  vx=5.96286238644834141E-03, vy=-2.43281292146216750E-02, vz=-2.53463209848734695E-03, m=1.66013679527193009E-07,
  r=20.0e0, d=5.43)
  venus = BodyCart("big",name="VENUS", x=6.33469157915745540E-01, y=3.49855234102151691E-01, z=-3.17853172088953667E-02,
  vx=-9.84258038001823571E-03, vy=1.76183746921837227E-02, vz=8.08822351013463794E-04, m=2.44783833966454430E-06,
  r=20.0e0, d=5.24)
  earthmoo = BodyCart("big", name="EARTHMOO", x=2.42093942183383037E-01, y=-9.87467766698604366E-01, z=-4.54276292555233496E-06,
  vx=1.64294055023289365E-02, vy=4.03200725816140870E-03, vz=1.13609607260006795E-08, m=3.04043264264672381E-06,
  r=20.0e0, d=5.52)
  mars = BodyCart("big", name="MARS", x=2.51831018120174499E-01, y=1.52598983115984788E+00, z=2.57781137811807781E-02,
  vx=-1.32744166042475433E-02, vy=3.46582959610421387E-03, vz=3.98930013246952611E-04, m=3.22715144505386530E-07,
  r=20.0e0, d=3.94)
  jupiter = BodyCart("big", name="JUPITER", x=4.84143144246472090E+00, y=-1.16032004402742839E+00, z=-1.03622044471123109E-01,
  vx=1.66007664274403694E-03, vy=7.69901118419740425E-03, vz=-6.90460016972063023E-05, m=9.54791938424326609E-04,
  r=3.0e0, d=1.33)
  saturn = BodyCart("big", name="SATURN", x=8.34336671824457987E+00, y=4.12479856412430479E+00, z=-4.03523417114321381E-01,
  vx=-2.76742510726862411E-03, vy=4.99852801234917238E-03, vz=2.30417297573763929E-05, m=2.85885980666130812E-04,
  r=3.0e0, d=0.70)
  uranus = BodyCart("big", name="URANUS", x=1.28943695621391310E+01, y=-1.51111514016986312E+01, z=-2.23307578892655734E-01,
  vx=2.96460137564761618E-03, vy=2.37847173959480950E-03, vz=-2.96589568540237556E-05, m=4.36624404335156298E-05,
  r=3.0e0, d=1.30)
  neptune = BodyCart("big", name="NEPTUNE", x=1.53796971148509165E+01, y=-2.59193146099879641E+01, z=1.79258772950371181E-01,
  vx=2.68067772490389322E-03, vy=1.62824170038242295E-03, vz=-9.51592254519715870E-05, m=5.15138902046611451E-05,
  r=3.0e0, d=1.76)
  pluto = BodyCart("big", name="PLUTO", x=-1.15095623952731607E+01, y=-2.70779438829451422E+01, z=6.22871533567077229E+00,
  vx=2.97220056963797431E-03, vy=-1.69820233395912967E-03, vz=-6.76798264809371094E-04, m=7.39644970414201173E-09,
  r=3.0e0, d=1.1)
  
  apollo = BodyAst("small", name="APOLLO", a=1.4710345, e=.5600245, I=6.35621, 
  g=285.63908, n=35.92313, M=15.77656, ep=2450400.5)
  jason = BodyAst("small", name="JASON", a=2.2157309, e=.7644575, I=4.84834, 
  g=336.49610, n=169.94137, M=293.37226, ep=2450400.5)
  khufu = BodyAst("small", name="KHUFU", a=0.9894948, e=.4685310, I=9.91298, 
  g=54.85927, n=152.64772, M=66.69818, ep=2450600.5)
  minos = BodyAst("small", name="MINOS", a=1.1513383, e=.4127106, I=3.93863, 
  g=239.50170, n=344.85893, M=8.93445, ep=2450400.5)
  orpheus = BodyAst("small", name="ORPHEUS", a=1.2091305, e=.3226805, I=2.68180, 
  g=301.55128, n=189.79654, M=28.31467, ep=2450400.5)
  toutatis = BodyAst("small", name="TOUTATIS", a=2.5119660, e=.6335854, I=0.46976, 
  g=274.82273, n=128.20968, M=50.00728, ep=2450600.5)
  
  solarSystem = PlanetarySystem(bodies=[mercury, venus, earthmoo, mars, jupiter, saturn, uranus, 
  neptune, pluto, apollo, jason, khufu, minos, orpheus, toutatis], m_star=1.0, epoch=2451000.5)
  
  bigin = Big(solarSystem)
  bigin.write()
  
  smallin = Small(solarSystem)
  smallin.write()
  
  elementin = Element(format_sortie=" a8.5 e8.6 i8.4 g8.4 n8.4 l8.4 m13e ", coord="Cen", 
  output_interval=365.2e1, time_format="years", relative_time="yes")
  elementin.write()
  
  closein = Close(time_format="years", relative_time="yes")
  closein.write()
  
  paramin = Param(algorithme="HYBRID", start_time=2451179.5, stop_time=2462502.5, output_interval=365.25e0, 
  h=8, accuracy=1.e-12, stop_integration="no", collisions="no", fragmentation="no", 
  time_format="years", relative_time="no", output_precision="medium", relativity="no", 
  user_force="no", ejection_distance=100, radius_star=0.005, central_mass=1.0, 
  J2=0, J4=0, J6=0, changeover=3., data_dump=500, periodic_effect=100)
  paramin.write()
  
  filesin = Files()
  filesin.write()
  
  messagein = Message()
  messagein.write()
  
  bigin2 = readBig()
  bigin2.write()
  
  diskin = Disk(b_over_h=0.4, adiabatic_index=1.4, mean_molecular_weight=2.35, surface_density=(500, 0.5), 
            viscosity=1e15, sample=200, dissipation_type=0)
  diskin.write()
  
  pdb.set_trace()
  
  for file in ["big.in", "small.in", "param.in", "element.in", "close.in", "message.in", "files.in"]:
    os.remove(file)
  
  print(mercury)
#TODO
# Tester si les smalls ont la même époque. Si ce n'est pas le cas, alors aucun ne doit avoir de masse. 

# Diskin() : In earlier versions, there were a "temperature=(500,1)", a "dissipation=False" (replaced by dissipation_type=0") and a "nb_sample=200" (replaced by sample=200) options

