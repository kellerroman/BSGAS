#!/usr/bin/env python3
# -*- coding: utf-8 -*-
######################################################################################################
# DEPENDENCY MAKEFILE CREATER FOR FORTRAN PROGRAMS
# AUTHOR: ROMAN KELLER
# START DATE : 12.03.2015
# LAST CHANGE: 12.03.2015, RK
# CREATES A DEPENDENCY TREE FROM A SOURCE FOLDER AND WRITES IT TO A MAKEFILE
# OPTIONAL: EXTRACT DEPENDENCY TREE
# OPTIONAL: CREATE COMPILER FLAG DEPENDENCIES
######################################################################################################
#
# USAGE:
#
# ./make_dependencies.py [SOURCE_FOLDER]
#
# Es kann auch komplett ohne KOMMANDOZEILENPARAMETER gestartet werden. 
# Dann erfolgen die Arbeitschritte interaktiv.
#
######################################################################################################
#
# CHANGELOG:
# 12.02.2015,RK:     START DES CHANGELOGS:
# 13.02.2015,RK: DEPENDENCIES FÜR CFLAGS HINZUGEFÜGT
# 17.04.2015,RK: CONVERTED TO PYTHON3
######################################################################################################
#
# TODO:
#
######################################################################################################


import sys
import os

def get_filepaths(directory):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = filename
            filepath = os.path.join(root, filename)
            filepath = os.path.relpath(filepath,directory)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.

######################################################################################################
################## CONFIG

endungen = [".F90",".F",".f",".f90"]
# DATEIENDUNGE DER SOURCEODE_DATEI WELCHE DURCHSUCHT WERDEN SOLLEN
source_name = "_SOURCES"
# MAKEFILE VARIABLEN NAME ZU DME DIE EINZELNEN CODEDATEIEN HINZUGEFÜGT WERDEN
object_prefix=""
# PREFIX, WELCHER VOR DIE EINZELNEN OBJECTE GESTELLT WIRD
cflags_dir=""
# ORDNER IN WELCHEM DIE CFLAGS DATEINE GESPEICHERT WERDEN
cflags_exclude=["comments","use_mpi"]
# CFLAGS DIE IGNORIERT WERDEN SOLLEN
source_makefile="inc_Makefile_src"
# DATEINAME FÜR DAS MAKEFILE MIT DEN SOURCE FILES
depend_makefile="inc_Makefile_dep"
# DATEINAME FÜR DAS MAKEFILE MIT DEN DEPENCIES (OBJECTS UND CFLAGS)

print_dep = False
print_cflags = False
include_cflags_dependencies = False


################## END  CONFIG
######################################################################################################

print("DEPENDENCY MAKEFILE CREATER FOR FORTRAN PROGRAMS")
print("V 0001 - 13.03.2015")
folder = ""
files = []
routine_links = []
lines = []
structure = {}
link = []
cflag_link = []
cflags = []

for i in sys.argv:
    il = i.lower()
    if il ==  "-debug":
        do_debug = True
        print("-------------------------------------------------")
        print("------------RUNNING IN DEBUG MODE----------------")
        print("-------------------------------------------------")
    elif os.path.isdir(i):
        folder = i
    elif il == "-print_cflags":
        print_cflags = True
    elif il == "-print_dep":
        print_dep = True
    elif il == "-exclude_cflags":
        include_cflags_dependencies = False

if folder == "":
    print("Folder to work on:")
    while 1==1:
        folder = input()
        if os.path.isdir(folder):
            break 
        else:
            print("Folder not found! Try again:")
    
print("Working on Folder:",folder)
for f in get_filepaths(folder):
  for end in endungen:
    if f.endswith(end):
      files.append(f)
      routine_links.append([])
      cflag_link.append([])
for i,f in enumerate(files):
   #print f
   fobj = open(os.path.join(folder,f), "r") 
   for line in fobj:
      lines.append(line.rstrip())
      l = line.split("!")[0].strip().lower()
      if l.find("module ") == 0 or l.find("subroutine ") == 0 or l.find("function ") == 0 or l.find("program ") == 0:
#         if l.find("function ") == 0:
#            print l
         l = l.split(" ")[1].split("(")[0]
         #print l,"is in",files[i]

         structure[l] = i
      if "call " in l:
         l = l.split("call ")[1].strip().split(" ")[0].split("(")[0]
         #print l,"Called in",files[i]
         if l not in routine_links[i]:
             routine_links[i].append(l)
      if "use " in l:
         l = l.split("use ")[1].strip().split(" ")[0].split("(")[0].split(",")[0]
         #print l,"Called in",files[i]
         if l not in routine_links[i]:
             routine_links[i].append(l)
      if l.find("#") == 0 and include_cflags_dependencies:
         if "ifdef" in l or "ifndef" in l:
            l = l.split(" ")[1]
            #print l
            if l not in cflags_exclude:
               if l not in cflag_link[i]:
                 cflag_link[i].append(l)
                 if l not in cflags:
                    cflags.append(l)
   fobj.close()

file_dep = []


if print_cflags:
   for l in cflags:
      print(l, end=' ')
      for i,f in enumerate(cflag_link):
         if l in f:
            print(f)


for i,f in enumerate(files):
   if print_dep:
      print(f)
   file_dep.append([])
   for l in routine_links[i]:
      if structure.get(l,-1) != -1 and structure.get(l,-1) != i:
         if print_dep:
            print("   ->",l,"form",files[structure.get(l)])
         if structure.get(l) not in file_dep[i]:
            file_dep[i].append(structure.get(l))


fsrc = open(source_makefile,"w+")
print(source_name+"+=", end=' ', file=fsrc)
for i,f in enumerate(files):
  print(os.path.basename(f), end=' ', file=fsrc)
print(file=fsrc)
fsrc.close()

fdep = open(depend_makefile,"w+")
for i,f in enumerate(files):
  if len(file_dep[i]) > 0 or len(cflag_link[i]) > 0:
   print(object_prefix+os.path.splitext(os.path.basename(f))[0]+".o"+":", end=' ', file=fdep)

   for l in file_dep[i]:
      print(object_prefix+os.path.basename(files[l]).split(".")[0]+".o", end=' ', file=fdep)

   for l in cflag_link[i]:
      print(cflags_dir+l, end=' ', file=fdep)
   print(file=fdep)
fdep.close()

quit()

scan_lines.scan_lines(myglobal.input_lines)




