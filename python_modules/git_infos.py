#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Module librement utilisable
"""module to retrieve information about a git repository and create a source code containing thoses informations.


HOW TO USE : 
import git_infos
git_infos.write_infos_in_f90_file()
"""
from __future__ import print_function

__author__ = "Autiwa <autiwa@gmail.com>"
__date__ = "04 october 2012"
__version__ = "$Revision: 1.0 $"
__credits__ = """Thanks to Bastien for his help when I was learning Python and his usefull tutorial"""

import subprocess
import pdb

def run_command(commande):
  """lance une commande qui sera typiquement soit une liste, soit une 
  commande seule. La fonction renvoit un tuple avec la sortie, 
  l'erreur et le code de retour"""
  if (type(commande)==list):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  elif (type(commande)==str):
    process = subprocess.Popen(commande, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
  else:
    raise TypeError("The command is neither a string nor a list.")
  (process_stdout, process_stderr) = process.communicate()
  returncode = process.poll()
  # there is .poll() or .wait() but I don't remember the difference. For some kind of things, one of the two was not working
  return (process_stdout, process_stderr, returncode)

def get_current_branch():
  """function that return as a string the current branch of the git repository"""
  (stdout, stderr, returnCode) = run_command("git branch")
  
  if (returnCode != 0):
    return None
  
  lines = stdout.split("\n")
  for line in lines:
    if (line[0] == '*'):
      return line[2:]

def get_current_revision():
  """function that return as a string the current revision of the git repository"""
  (stdout, stderr, returnCode) = run_command("git log|head -1")
  
  if (returnCode != 0):
    return None
  
  commit = stdout.split()[1]
  return commit

def is_non_committed_modifs():
  """function that return as a boolean if theere is non committed modifications in the repository"""
  (stdout, stderr, returnCode) = run_command("git diff|wc -l")
  
  if (returnCode != 0):
    return None
  
  nbLines = int(stdout)
  
  return (nbLines != 0)
  
def list_tag(commit):
  """list the tags that exists linking towar the considered commit
  
  Return :
  The list of tags corresponding to 'commit'. If none, an empty list is returned.
  """
  (stdout, stderr, returnCode) = run_command("git tag -l --contains %s" % commit)
  
  tags = stdout.split("\n")[0:-1] # We do not include the extra "" in the end.
  
  return tags 

def write_infos_in_f90_file(main_branch='master'):
  """This function will create a fortran file that will store, as variable, some infos about a git repository"""
  
  F90_BEGIN = "module git_infos\n" + \
              "! Automatically generated file through Makefile.py, do not modify manually !\n" + \
              "implicit none\n\n"

  F90_END = "\nend module git_infos"

  
  branch = get_current_branch()
  commit = get_current_revision()
  isModifs = is_non_committed_modifs()
  tags = list_tag(commit)
  
  if (branch != main_branch):
    print("Warning: The current branch is %s" % branch)
  
  f90source = open("git_infos.f90", 'w')
  f90source.write(F90_BEGIN)
  f90source.write("character(len=40), parameter :: commit = '%s'\n" % commit)
  f90source.write("character(len=%d), parameter :: branch = '%s'\n" % (len(branch), branch))
  if (isModifs):
    f90source.write("character(len=80), parameter :: modifs = '/!\ There is non committed modifications'\n")
  else:
    f90source.write("character(len=80), parameter :: modifs = 'This is a pure version (without any local modifs)'\n")
    
  
  if (len(tags)==0):
    tag_text = "There is no tag"
  else:
    tag_text = " ; ".join(tags)
  
  f90source.write("character(len=%d), parameter :: tags = '%s'\n" % (len(tag_text), tag_text))
  f90source.write(F90_END)
  f90source.close()
  
  
