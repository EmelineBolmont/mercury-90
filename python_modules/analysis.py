#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""module that contains generic functions that help analyse simulations of planetary dynamic."""

__author__ = "Autiwa <autiwa@gmail.com>"
__date__ = "2011-09-21"
__version__ = "1.1"

import numpy as np
from fractions import Fraction

def get_x_s(mass):
  """with the mass of the planet in solar mass, 
  and for planets around 3AU in my disk (for the fixed value of 'h')"""
  mstar = 1.0 # the stellar mass in solar mass
  b_over_h = 0.4
  h = 0.045
  q = mass / mstar
  adiabatic_index = 1.4
  
  x_s = (1.1 * (0.4 / b_over_h)**0.25 / adiabatic_index**0.25) * np.sqrt(q / h)
  
  return x_s

def get_possible_resonances(periodRatio, uncertainty=0.05, denominator_limit=12, numerator_limit=20, sampling=10):
  """Give a list of 'Fraction' objects that correspond to possible Mean Motion Resonances for a given period ratio. 
  Optional argument are the uncertainty that will determine the range of perio ratios to test around the given value. 
  The denominator and numerator limit above which all the resonances will be skipped
  
  Parameter :
  periodRatio : [float] the periodRatio between the two considered parameters
  
  Optional parameters :
  uncertainty=0.05 : the period range to test will be [p*(1-0.05) ; p * (1+0.05)]
  denominator_limit=12 : 13:12 resonance will be ok, but 14:13 resonance will be skipped
  numerator_limit=20 : 20:12 resonance will be ok, but 21:12 resonance will be skipped
  sampling=10 : The number of period ratio values for the range of values of period ratio to be tested around the nominal value 
                (given the uncertainty)
  
  Return :
  list of 'Fraction' objects, each one representing a Mean Motion Resonance to be tested  
  """
  
  periodMin = periodRatio * (1 - uncertainty)
  periodMax = periodRatio * (1 + uncertainty)
  
  # We do not want period ratios less than 1 (this only happens for coorbitals I think)
  if (periodMin < 1.):
    periodMin = 1.
  
  periodWidth = periodMax - periodMin
  deltaPeriod = periodWidth/sampling
  
  
  periods = [periodMin + deltaPeriod * i for i in range(sampling)]

  resonances = []
  for period_i in periods:
    fraction = Fraction(period_i).limit_denominator(denominator_limit)
    resonances.append(fraction)

  # We exclude all values that appears several time to only keep one occurence of each value
  resonances = list(set(resonances))

  # We sort the resonances to get the more interesting first (3:2 before 32:27 for instance)
  tmp = [(res.numerator, res) for res in resonances]
  tmp.sort()
  resonances = [element[1] for element in tmp if element[1].numerator < numerator_limit]
  
  #~ print(uncertainty)
  #~ print(periodMin)
  #~ print(periodMax)
  #~ print(resonances)
  #~ exit()
  
  return resonances

def isResonance(res, g_inner, n_inner, M_inner, g_outer, n_outer, M_outer, nb_points=50, angle_center_value=0, std_threshold=20.):
  """Given a resonance as a Fraction object, and g, n M for inner and
  outer planet, the function return if there is the resonance between
  the two planets
  
  Parameters : 
  res : a Fraction object (for instance Fraction(3,2))
  g_inner, n_inner, M_inner : g, n, M for the inner planet (in degrees)
  g_outer, n_outer, M_outer : g, n, M for the outer planet (in degrees)
  
  Optional parameters :
  nb_points : [50] The number of points, at the end the each time span, that we will use to test resonance
  angler_center_value : [0.] (in degrees) we will overlap all the angle, through congruence, to have angle centered over this value. 
                        For instance, with 0, angle will be between -180 and +180 degrees.
  std_threshold = [20.] in degrees, the value below which we will consider that a resonant angle prove the existence of a MMR.
  
  Return : 
  True or False
  """
  outer_period_nb = res.denominator
  inner_period_nb = res.numerator
  
  angle_min = angle_center_value - 180.
  angle_max = angle_center_value + 180.
  
  # Resonances are usually displayed as (p+q):p where q is the order of
  # the resonance. We retreive thoses parameters
  p = outer_period_nb
  q = inner_period_nb - outer_period_nb

  # We calculate the resonant angles
  long_of_peri_inner = g_inner + n_inner
  mean_longitude_inner = M_inner + long_of_peri_inner

  long_of_peri_outer = g_outer + n_outer
  mean_longitude_outer = M_outer + long_of_peri_outer

  phi = np.empty((q+1, nb_points)) # we create an array, empty for the moment, that will contain all the resonant angles associated with the supposed resonance.

  temp_value = inner_period_nb * mean_longitude_outer - outer_period_nb * mean_longitude_inner

  for i in range(q+1):
    phi[i] = temp_value - i * long_of_peri_inner - (q - i) * long_of_peri_outer

  # We take modulo 2*pi of the resonant angle
  phi = phi%(360.)
  too_low = phi < angle_min
  too_high = phi > angle_max
  phi[too_low] = phi[too_low] + 360.
  phi[too_high] = phi[too_high] - 360.

  delta_longitude = long_of_peri_outer - long_of_peri_inner
  delta_longitude = delta_longitude%(360.)
  too_low = delta_longitude < angle_min
  too_high = delta_longitude > angle_max
  delta_longitude[too_low] = delta_longitude[too_low] + 360.
  delta_longitude[too_high] = delta_longitude[too_high] - 360.
  
  # If one of the std's is small (Typically, around 25, but 
  # I had once a 80 that was not a resonance, so I think a threshold around 40 is a good one)
  standard_deviation = min(phi.std(1))
    
  if (standard_deviation < std_threshold):
    return True
  else:
    return False
