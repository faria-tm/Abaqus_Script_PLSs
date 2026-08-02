# ---------------------------------------------------------------------------
# Model-creator-module-loadcellless
# ---------------------------------------------------------------------------
# Builds a parametric unit-cell + load-cell Abaqus/CAE model for each
# "fingerprint" (a string of 0/1 digits describing which of the 20 wall
# positions in the unit cell are populated) listed in FINGERPRINT_FILE.
#
# For every fingerprint the script:
#   1. Constructs the unit-cell geometry (normal, diagonal, and triangular
#      wall parts), including drainage holes for resin removal.
#   2. Assembles the walls indicated by the fingerprint's binary digits,
#      merges them into a single part, and mirrors it to complete the cell.
#   3. Adds a simple steel "load-cell" fixture, materials, contact
#      interactions, an implicit dynamics step, boundary conditions, and
#      a mesh.
#   4. Computes the part mass and appends it to a results file.
#   5. Saves the finished model as its own .cae file.
#
# Requirements:
#   - Abaqus/CAE (run via `abaqus cae noGUI=Model-creator-module-loadcellless_0.25.py`
#     or from within the Abaqus/CAE Python kernel).
#   - A FINGERPRINT_FILE (default: fingerprint.txt) in the working directory,
#     with one fingerprint string per line, e.g.:
#       11000000000000000000
#       10100000000000000000
#
# Output:
#   - OUTPUT_DIR/mass_results.txt  - "fingerprint, mass_kg" for each cell
#   - OUTPUT_DIR/<fingerprint>.cae - one CAE file per successfully built cell
#
# Units follow Abaqus' unit-agnostic convention consistent with mm, tonne,
# N, MPa, s (hence density given as e.g. 2.85e-09 tonne/mm^3).
# ---------------------------------------------------------------------------

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
from array import *
import math
import numpy
import os, glob
import shutil
from abaqus import *
from abaqusConstants import *
import multiprocessing
import ctypes
import mesh
from math import pi

# ---------------------------------------------------------------------------
# Configuration - edit these for your local setup
# ---------------------------------------------------------------------------
FINGERPRINT_FILE = 'fingerprint.txt'
OUTPUT_DIR = './output'

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

session.journalOptions.setValues(replayGeometry=COORDINATE)
mdb = Mdb()

with open(FINGERPRINT_FILE, 'r') as f:
    input_files = [line.strip() for line in f.readlines()]
print (input_files)

ff = open(os.path.join(OUTPUT_DIR, 'mass_results.txt'), 'w')
ff.write("fingerprint, mass_kg\n")

for fingerprint in enumerate(input_files):
    try:
        fingerprint = fingerprint[1]
        #Normal_Wall_Properties---------------------------------------------------
        #field_inputs = (('Enter the width (mm):','10'),('Enter wall thickness (mm):','0.25'),('Enter Drainage holes' radius (mm):','0.7'),('Enter your structure fingerprint :','11000000000000000000'))
        #width,thickness,hole_radius,fingerprint = getInputs(fields=field_inputs,label='Specimen Properties',dialogTitle = 'Create Specimen')
        width = 10
        thickness = 0.25
        hole_radius = 0.7
        print (fingerprint)
        walls = [int(digit) for digit in fingerprint]

        n = 0
        for i in range (len(walls)):
           if walls[i]==1:
              n =n+1
        #Diagonal_Wall_Parameters-------------------------------------------------------------
        x = thickness/(2*cos(pi/4))
        width_D = sqrt(2*(width-x))
        side_D = width-2*thickness

        #Triangular_Wall_Parameters-------------------------------------------------------------
        Width_T = sqrt(2)*(width-0.0)

        #Sketch_Normal_Wall--------------------------------------------------------------
        mymodel = mdb.models['Model-1']
        mysketch1 = mymodel.ConstrainedSketch('Firstsketch',500)
        g1, v1, d1, c1 = mysketch1.geometry, mysketch1.vertices, mysketch1.dimensions, mysketch1.constraints
        mysketch1.rectangle(point1=(0.0, 0.0), point2=(width, width))
        part1 = mymodel.Part(name='Part-1', dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        part1.BaseSolidExtrude(sketch=mysketch1, depth=thickness)

        f, e = part1.faces, part1.edges
        t = part1.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(3.333333, 3.333333,
            0.25)), sketchUpEdge=e.findAt(coordinates=(10.0, 7.5, 0.25)),
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(5.0, 5.0,
            0.25))
        s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
            sheetSize=28.28, gridSpacing=0.7, transform=t)
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        s1.setPrimaryObject(option=SUPERIMPOSE)
        part1.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
        s1.ArcByCenterEnds(center=(0.0, -5.0), point1=(-hole_radius, -5.0), point2=(hole_radius, -5.0),
            direction=CLOCKWISE)
        s1.Line(point1=(-hole_radius, -5.0), point2=(hole_radius, -5.0))
        f1, e1 = part1.faces, part1.edges
        part1.CutExtrude(sketchPlane=f1.findAt(coordinates=(3.333333, 3.333333, 0.25)),
            sketchUpEdge=e1.findAt(coordinates=(10.0, 7.5, 0.25)),
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s1,
            flipExtrudeDirection=OFF)
        s1.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']


        #Sketch_Diagonal_Wall_NEW--------------------------------------------------------------
        mysketch2 = mymodel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        g2, v2, d2, c2 = mysketch2.geometry, mysketch2.vertices, mysketch2.dimensions, mysketch2.constraints
        mysketch2.setPrimaryObject(option=STANDALONE)
        mysketch2.Line(point1=(0.0, 0.0), point2=(x, 0.0))
        mysketch2.Line(point1=(x, 0.0), point2=(width, width-x))
        mysketch2.Line(point1=(width, width-x), point2=(width, width))
        mysketch2.Line(point1=(width, width), point2=(width-x, width))
        mysketch2.Line(point1=(width-x, width), point2=(0.0, x))
        mysketch2.Line(point1=(0.0, x), point2=(0.0, 0.0))
        part2 = mymodel.Part(name='Part-2', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part2.BaseSolidExtrude(sketch=mysketch2, depth=width)
        del mdb.models['Model-1'].sketches['__profile__']

        f, e = part2.faces, part2.edges
        t = part2.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(3.274408, 3.451184,
            6.666667)), sketchUpEdge=e.findAt(coordinates=(9.823223, 10.0, 2.5)),
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(4.911612,
            5.088388, 5.0))
        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
            sheetSize=34.43, gridSpacing=0.86, transform=t)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=SUPERIMPOSE)
        part2.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, -hole_radius))
        f1, e1 = part2.faces, part2.edges
        part2.CutExtrude(sketchPlane=f1.findAt(coordinates=(3.274408, 3.451184, 6.666667)),
            sketchUpEdge=e1.findAt(coordinates=(9.823223, 10.0, 2.5)),
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s,
            flipExtrudeDirection=OFF)
        s.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']


        #Sketch_Triangular_Wall--------------------------------------------------------------
        mysketch3 = mymodel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        g3, v3, d3, c3 = mysketch3.geometry, mysketch3.vertices, mysketch3.dimensions, mysketch3.constraints
        mysketch3.setPrimaryObject(option=STANDALONE)
        mysketch3.Line(point1=(0.0, 0.0), point2=((width-0.0)*cos(pi/4), (width-0.0)*1.22474487))
        mysketch3.Line(point1=((width-0.0)*cos(pi/4), (width-0.0)*1.22474487), point2=(Width_T, 0.0))
        mysketch3.Line(point1=(Width_T, 0.0), point2=(0.0, 0.0))
        part3 = mymodel.Part(name='Part-3', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part3.BaseSolidExtrude(sketch=mysketch3, depth=thickness)

        #Sketch_Triangular_Wall_Inverse--------------------------------------------------------------
        mysketch4 = mymodel.ConstrainedSketch(name='__profile__', sheetSize=200.0)
        g4, v4, d4, c4 = mysketch4.geometry, mysketch4.vertices, mysketch4.dimensions, mysketch4.constraints
        mysketch4.setPrimaryObject(option=STANDALONE)
        mysketch4.Line(point1=((width-0.0)*cos(pi/4), 0.0), point2=(Width_T, (width-0.0)*1.22474487))
        mysketch4.Line(point1=(Width_T, (width-0.0)*1.22474487), point2=(0.0, (width-0.0)*1.22474487))
        mysketch4.Line(point1=(0.0, (width-0.0)*1.22474487), point2=((width-0.0)*cos(pi/4), 0.0))
        part4 = mymodel.Part(name='Part-4', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part4.BaseSolidExtrude(sketch=mysketch4, depth=thickness)

        #Unit Cell Assembly---------------------------------------------------------------
        #Plate_13----------------------------------------------------------------
        if walls[12]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-13', part=part3, dependent=ON)
            myassembly.rotate(instanceList=('Part-plate-13', ), axisPoint=(0.0, 0.0, 0.0),
                axisDirection=(0.0, -1.0, 0.0), angle=45.0)
            myassembly.rotate(instanceList=('Part-plate-13', ), axisPoint=(0.0, 0.0, 0.0),
                axisDirection=(cos(pi/4), 0.0, cos(pi/4)), angle=-35.2644)

        #Plate_14----------------------------------------------------------------
        if walls[13]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-14', part=part3, dependent=ON)
            myassembly.translate(instanceList=('Part-plate-14', ), vector=(0.0, 0.0, -thickness*0.75))
            myassembly.rotate(instanceList=('Part-plate-14', ), axisPoint=(0.0, 0.0, 0.0),
                axisDirection=(0.0, -1.0, 0.0), angle=45.0)
            myassembly.rotate(instanceList=('Part-plate-14', ), axisPoint=(0.0, 0.0, 0.0),
                axisDirection=(cos(pi/4), 0.0, cos(pi/4)), angle=35.2644)

        #Plate_15----------------------------------------------------------------
        if walls[14]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-15', part=part3, dependent=ON)
            myassembly.translate(instanceList=('Part-plate-15', ), vector=(-(Width_T-width), 0.0, 0.0))
            myassembly.rotate(instanceList=('Part-plate-15', ), axisPoint=(width, 0.0, 0.0),
                axisDirection=(0.0, -1.0, 0.0), angle=-45.0)
            myassembly.rotate(instanceList=('Part-plate-15', ), axisPoint=(width, 0.0, 0.0),
                axisDirection=(-cos(pi/4), 0.0, cos(pi/4)), angle=35.2644)

        #Plate_16----------------------------------------------------------------
        if walls[15]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-16', part=part3, dependent=ON)
            myassembly.translate(instanceList=('Part-plate-16', ), vector=(-(Width_T-width), 0.0, 0.0))
            myassembly.translate(instanceList=('Part-plate-16', ), vector=(0.0, 0.0, -thickness*0.75))
            myassembly.rotate(instanceList=('Part-plate-16', ), axisPoint=(width, 0.0, 0.0),
                axisDirection=(0.0, -1.0, 0.0), angle=-45.0)
            myassembly.rotate(instanceList=('Part-plate-16', ), axisPoint=(width, 0.0, 0.0),
                axisDirection=(-cos(pi/4), 0.0, cos(pi/4)), angle=-35.2644)

         #Plate_17----------------------------------------------------------------
        if walls[16]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-17', part=part4, dependent=ON)
            myassembly.translate(instanceList=('Part-plate-17', ), vector=(0.0, -(((width-0.0)*1.22474487)-width), -thickness*0.75))
            myassembly.rotate(instanceList=('Part-plate-17', ), axisPoint=(0.0, width, 0.0),
                axisDirection=(0.0, -1.0, 0.0), angle=45.0)
            myassembly.rotate(instanceList=('Part-plate-17', ), axisPoint=(0.0, width, 0.0),
                axisDirection=(cos(pi/4), 0.0, cos(pi/4)), angle=-35.2644)


        #Plate_18----------------------------------------------------------------
        if walls[17]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-18', part=part4, dependent=ON)
            myassembly.translate(instanceList=('Part-plate-18', ), vector=(0.0, -(((width-0.0)*1.22474487)-width), 0.0))
            myassembly.rotate(instanceList=('Part-plate-18', ), axisPoint=(0.0, width, 0.0),
                axisDirection=(0.0, -1.0, 0.0), angle=45.0)
            myassembly.rotate(instanceList=('Part-plate-18', ), axisPoint=(0.0, width, 0.0),
                axisDirection=(cos(pi/4), 0.0, cos(pi/4)), angle=35.2644)

        #Plate_19----------------------------------------------------------------
        if walls[18]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-19', part=part4, dependent=ON)
            myassembly.translate(instanceList=('Part-plate-19', ), vector=(-(Width_T-width), -(((width-0.0)*1.22474487)-width), -thickness*0.75))
            myassembly.rotate(instanceList=('Part-plate-19', ), axisPoint=(width, width, 0.0),
                axisDirection=(0.0, -1.0, 0.0), angle=-45.0)
            myassembly.rotate(instanceList=('Part-plate-19', ), axisPoint=(width, width, 0.0),
                axisDirection=(-cos(pi/4), 0.0, cos(pi/4)), angle=35.2644)

        #Plate_20----------------------------------------------------------------
        if walls[19]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-20', part=part4, dependent=ON)
            myassembly.translate(instanceList=('Part-plate-20', ), vector=(-(Width_T-width), -(((width-0.0)*1.22474487)-width), 0.0))
            myassembly.rotate(instanceList=('Part-plate-20', ), axisPoint=(width, width, 0.0),
                axisDirection=(0.0, -1.0, 0.0), angle=-45.0)
            myassembly.rotate(instanceList=('Part-plate-20', ), axisPoint=(width, width, 0.0),
                axisDirection=(-cos(pi/4), 0.0, cos(pi/4)), angle=-35.2644)

        #Plate_1----------------------------------------------------------------
        if walls[0]==1:
           myassembly = mdb.models['Model-1'].rootAssembly
           myassembly.Instance(name='Part-plate-1', part=part1, dependent=ON)

        #Plate_2----------------------------------------------------------------
        if walls[1]==1:
             myassembly = mdb.models['Model-1'].rootAssembly
             myassembly.Instance(name='Part-plate-2', part=part1, dependent=ON)
             myassembly.rotate(instanceList=('Part-plate-2', ), axisPoint=(width, 0.0, 0.0),
                 axisDirection=(0.0, 1.0, 0.0), angle=90.0)
             myassembly.translate(instanceList=('Part-plate-2', ), vector=(-thickness, 0.0, 0.0))

        #Plate_3----------------------------------------------------------------
        if walls[2]==1:
             myassembly = mdb.models['Model-1'].rootAssembly
             myassembly.Instance(name='Part-plate-3', part=part1, dependent=ON)
             myassembly.translate(instanceList=('Part-plate-3', ), vector=(0.0, 0.0, width-thickness))

        #Plate_4----------------------------------------------------------------
        if walls[3]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-4', part=part1, dependent=ON)
            myassembly.rotate(instanceList=('Part-plate-4', ), axisPoint=(0.0, 0.0, 0.0),
                axisDirection=(0.0, 1.0, 0.0), angle=-90.0)
            myassembly.translate(instanceList=('Part-plate-4', ), vector=(thickness, 0.0, 0.0))

        #Plate_5----------------------------------------------------------------
        if walls[4]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-5', part=part1, dependent=ON)
            myassembly.rotate(instanceList=('Part-plate-5', ), axisPoint=(0.0, 0.0, 0.0),
                axisDirection=(1.0, 0.0, 0.0), angle=90.0)
            myassembly.translate(instanceList=('Part-plate-5', ), vector=(0.0, thickness, 0.0))

        #Plate_6----------------------------------------------------------------
        if walls[5]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-6', part=part1, dependent=ON)
            myassembly.rotate(instanceList=('Part-plate-6', ), axisPoint=(0.0, width, 0.0),
                axisDirection=(1.0, 0.0, 0.0), angle=-90.0)
            myassembly.translate(instanceList=('Part-plate-6', ), vector=(0.0, -thickness, 0.0))

        #Plate_7----------------------------------------------------------------
        if walls[6]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-7', part=part2, dependent=ON)
        #Plate_8----------------------------------------------------------------
        if walls[7]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-8', part=part2, dependent=ON)
            myassembly.rotate(instanceList=('Part-plate-8', ), axisPoint=(thickness, thickness, thickness),
                axisDirection=(0.0, 0.0, 1.0), angle=45.0)
            myassembly.rotate(instanceList=('Part-plate-8', ), axisPoint=(thickness, thickness, 0.0),
                    axisDirection=(0.0, 1.0, 0.0), angle=90.0)
            myassembly.rotate(instanceList=('Part-plate-8', ), axisPoint=(thickness, thickness, 0.0),
                    axisDirection=(1.0, 0.0, 0.0), angle=45.0)
            myassembly.translate(instanceList=('Part-plate-8', ), vector=(-thickness, 0.0, thickness))
            #myassembly.translate(instanceList=('Part-plate-8', ), vector=(0.0, 0.0, thickness-x))
        #Plate_9----------------------------------------------------------------
        if walls[8]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-9', part=part2, dependent=ON)
            myassembly.rotate(instanceList=('Part-plate-9', ), axisPoint=(thickness, thickness, thickness),
                axisDirection=(0.0, 0.0, 1.0), angle=45.0)
            myassembly.rotate(instanceList=('Part-plate-9', ), axisPoint=(thickness, thickness, width-thickness),
                    axisDirection=(0.0, 1.0, 0.0), angle=-90.0)
            myassembly.rotate(instanceList=('Part-plate-9', ), axisPoint=(thickness, thickness, width-thickness),
                    axisDirection=(1.0, 0.0, 0.0), angle=-45.0)
            #myassembly.translate(instanceList=('Part-plate-9', ), vector=(0.0, 0.0, thickness-x))
        #Plate_10----------------------------------------------------------------
        if walls[9]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-10', part=part2, dependent=ON)
            myassembly.translate(instanceList=('Part-plate-10', ), vector=(width, 0.0, 0.0))
            myassembly.rotate(instanceList=('Part-plate-10', ), axisPoint=(width, 0.0, 0.0),
                axisDirection=(0.0, 0.0, 1.0), angle=90.0)

        #Plate_11----------------------------------------------------------------
        if walls[10]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-11', part=part2, dependent=ON)
            myassembly.rotate(instanceList=('Part-plate-11', ), axisPoint=(thickness, thickness, thickness),
                axisDirection=(0.0, 0.0, 1.0), angle=-45.0)
            myassembly.rotate(instanceList=('Part-plate-11', ), axisPoint=(thickness, thickness, thickness),
                axisDirection=(1.0, 0.0, 0.0), angle=-90.0)
            myassembly.rotate(instanceList=('Part-plate-11', ), axisPoint=(thickness, thickness, thickness),
                axisDirection=(0.0, 1.0, 0.0), angle=-45.0)
            #myassembly.translate(instanceList=('Part-plate-11', ), vector=(0.0, 0.0, thickness-x))

        #Plate_12----------------------------------------------------------------
        if walls[11]==1:
            myassembly = mdb.models['Model-1'].rootAssembly
            myassembly.Instance(name='Part-plate-12', part=part2, dependent=ON)
            #myassembly.translate(instanceList=('Part-plate-12', ), vector=(0.0, 0.0, thickness))
            myassembly.rotate(instanceList=('Part-plate-12', ), axisPoint=(thickness, thickness, thickness),
                axisDirection=(0.0, 0.0, 1.0), angle=-45.0)
            myassembly.rotate(instanceList=('Part-plate-12', ), axisPoint=(thickness, thickness, width-thickness),
                axisDirection=(1.0, 0.0, 0.0), angle=90.0)
            myassembly.rotate(instanceList=('Part-plate-12', ), axisPoint=(thickness, thickness, width-thickness),
                axisDirection=(0.0, 1.0, 0.0), angle=45.0)

        #Hole on part3------------------------------------------------------------
        e = part3.edges
        part3.Round(radius=0.125, edgeList=(e[0], e[2], e[4], e[6], e[8], e[7]))
        e1 = part3.edges
        part3.Round(radius=0.125, edgeList=(e1[3], e1[4], e1[6]))

        f, e = part3.faces, part3.edges
        t = part3.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(7.071068, 4.082483,
            0.25)), sketchUpEdge=e.findAt(coordinates=(8.784708, 9.029337, 0.25)),
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(7.071068,
            4.082483, 0.25))
        s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
            sheetSize=36.76, gridSpacing=0.91, transform=t)
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        s1.setPrimaryObject(option=SUPERIMPOSE)
        part3.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
        s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, -hole_radius))
        f1, e1 = part3.faces, part3.edges
        part3.CutExtrude(sketchPlane=f1.findAt(coordinates=(7.071068, 4.082483, 0.25)),
            sketchUpEdge=e1.findAt(coordinates=(8.784708, 9.029337, 0.25)),
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s1,
            flipExtrudeDirection=OFF)
        s1.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']

        #Hole on part4------------------------------------------------------------
        e = part4.edges
        part4.Round(radius=0.125, edgeList=(e[0], e[2], e[4], e[6], e[8], e[7]))
        e1 = part4.edges
        part4.Round(radius=0.125, edgeList=(e1[3], e[4], e[6]))

        f, e = part4.faces, part4.edges
        t = part4.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(7.071068, 8.164966,
            0.25)), sketchUpEdge=e.findAt(coordinates=(12.211989, 9.154337, 0.25)),
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(7.071068,
            8.164966, 0.25))
        s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
            sheetSize=36.92, gridSpacing=0.92, transform=t)
        g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
        s1.setPrimaryObject(option=SUPERIMPOSE)
        part4.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
        s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, -hole_radius))
        f1, e1 = part4.faces, part4.edges
        part4.CutExtrude(sketchPlane=f1.findAt(coordinates=(7.071068, 8.164966, 0.25)),
            sketchUpEdge=e1.findAt(coordinates=(12.211989, 9.154337, 0.25)),
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s1,
            flipExtrudeDirection=OFF)
        s1.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']



        # Merging the unit cell-----------------------------------------------
        SingleInstances_List = myassembly.instances.keys()
        myassembly.InstanceFromBooleanMerge(name='UnitCell', instances=([myassembly.instances[SingleInstances_List[i]]
            for i in range(len(SingleInstances_List))] ), keepIntersections=OFF,
            originalInstances=DELETE, mergeNodes=BOUNDARY_ONLY,
            nodeMergingTolerance=1e-06, domain=BOTH)

        # Cutting the unit cell---------------------------------------------------
        part4 = mdb.models['Model-1'].parts['UnitCell']
        v2, e1 = part4.vertices, part4.edges
        part4.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=width)
        part4.DatumAxisByPrincipalAxis(principalAxis=YAXIS)
        part4.DatumPointByCoordinate(coords=(20.0, 20.0, 0.0))
        v2, d2 = part4.vertices, part4.datums
        part4.DatumAxisByParToEdge(edge=d2[3], point=d2[4])
        del mdb.models['Model-1'].sketches['Firstsketch']
        f, e, d2 = part4.faces, part4.edges, part4.datums
        Transform1 = part4.MakeSketchTransform(sketchPlane=d2[2], sketchUpEdge=d2[5], sketchPlaneSide=SIDE1,
            sketchOrientation=RIGHT, origin=(width, 0.0, 0.0))
        sketch4 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
            sheetSize=36.98, gridSpacing=0.92, transform=Transform1)
        g4, v4, d4, c4 = sketch4.geometry, sketch4.vertices, sketch4.dimensions, sketch4.constraints
        sketch4.Line(point1=(-thickness/2, width-0.025), point2=(-thickness/2, thickness/2))
        sketch4.Line(point1=(-thickness/2, thickness/2), point2=(-width+0.025, thickness/2))
        sketch4.Line(point1=(-width+0.025, thickness/2), point2=(-width+0.025, width-0.025))
        sketch4.Line(point1=(-width+0.025, width-0.025), point2=(-thickness/2, width-0.025))
        sketch4.Line(point1=(width, -width), point2=(width, 2*width))
        sketch4.Line(point1=(width, 2*width), point2=(-2*width, 2*width))
        sketch4.Line(point1=(-2*width, 2*width), point2=(-2*width, -width))
        sketch4.Line(point1=(-2*width, -width), point2=(width, -width))
        f1, e1, d1 = part4.faces, part4.edges, part4.datums
        part4.CutExtrude(sketchPlane=d1[2], sketchUpEdge=d1[5], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=sketch4,
            flipExtrudeDirection=OFF)
        sketch4.unsetPrimaryObject()
        e = part4.edges
        part4.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=width)
        part4.DatumPointByCoordinate(coords=(width, thickness/2, width))
        v2, d2 = part4.vertices, part4.datums
        part4.DatumAxisByParToEdge(edge=d2[5], point=d2[8])
        f, e1, d2 = part4.faces, part4.edges, part4.datums
        Transform2 = part4.MakeSketchTransform(sketchPlane=d2[7], sketchUpEdge=d2[9], sketchPlaneSide=SIDE1,
            sketchOrientation=RIGHT, origin=(width, 0.0, width))
        sketch5 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
            sheetSize=37.47, gridSpacing=0.93, transform=Transform2)
        g, v, d, c = sketch5.geometry, sketch5.vertices, sketch5.dimensions, sketch5.constraints
        sketch5.rectangle(point1=(-thickness/2,width ), point2=(-width+0.025, 0.0))
        sketch5.rectangle(point1=(2*width,2*width ), point2=(-2*width, -2*width))
        f1, e, d1 = part4.faces, part4.edges, part4.datums
        part4.CutExtrude(sketchPlane=d1[7], sketchUpEdge=d1[9], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
            sketch=sketch5, flipExtrudeDirection=OFF)
        del mdb.models['Model-1'].sketches['__profile__']

        #Mirroring-------------------------------------------------------------------
        part4 = mdb.models['Model-1'].parts['UnitCell']
        part4.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=width-thickness/2)
        f, e1, d2 = part4.faces, part4.edges, part4.datums
        part4.Mirror(mirrorPlane=d2[11], keepOriginal=ON)

        part4 = mdb.models['Model-1'].parts['UnitCell']
        part4.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=thickness/2)
        f1, e, d1 = part4.faces, part4.edges, part4.datums
        part4.Mirror(mirrorPlane=d1[13], keepOriginal=ON)

        part4 = mdb.models['Model-1'].parts['UnitCell']
        part4.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=thickness/2)
        f, e1, d2 = part4.faces, part4.edges, part4.datums
        part4.Mirror(mirrorPlane=d2[15], keepOriginal=ON)

        #LoadCell--------------------------------------------------------------------
        sketch6 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
            sheetSize=200.0)
        g6, v6, d6, c6 = sketch6.geometry, sketch6.vertices, sketch6.dimensions, sketch6.constraints
        sketch6.rectangle(point1=(-1.5*width, -1.5*width), point2=(1.5*width, 1.5*width))
        part5 = mdb.models['Model-1'].Part(name='Part-loadcell', dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        part5.BaseSolidExtrude(sketch=sketch6, depth=width/5)
        del mdb.models['Model-1'].sketches['__profile__']

        #LoadCell_Assembly-----------------------------------------------------------
        myassembly.Instance(name='Part-loadcell-1', part=part5, dependent=ON)
        myassembly.rotate(instanceList=('Part-loadcell-1', ), axisPoint=(0.0, 0, 0),
            axisDirection=(1, 0.0, 0.0), angle=90.0)
        myassembly.translate(instanceList=('Part-loadcell-1', ), vector=((2*width-0.5)/2, width+width/5-0.025, 0.0))

        myassembly.Instance(name='Part-loadcell-2', part=part5, dependent=ON)
        myassembly.rotate(instanceList=('Part-loadcell-2', ), axisPoint=(0.0, 0, 0),
            axisDirection=(1, 0.0, 0.0), angle=-90.0)
        myassembly.translate(instanceList=('Part-loadcell-2', ), vector=((2*width-0.5)/2, -((width-thickness)+width/5)+0.025, 0.0))

        # RP-----------------------------------------------------------------
        #myassembly.ReferencePoint(point=(-5.25, 11.975, 15.0))
        #myassembly.ReferencePoint(point=(-5.25, -11.725, 15.0))
        part5 = mdb.models['Model-1'].parts['Part-loadcell']
        v1, e1, d1, n1 = part5.vertices, part5.edges, part5.datums, part5.nodes
        part5.ReferencePoint(point=v1.findAt(coordinates=(-15.0, 15.0, 2.0)))

        # UnitCell_Material_Properties------------------------------------------------
        mdb.models['Model-1'].Material(name='Material-Unitcell')
        mdb.models['Model-1'].materials['Material-Unitcell'].Elastic(table=((1760.0, 0.35), ))
        mdb.models['Model-1'].materials['Material-Unitcell'].Density(table=((2.85e-09, ), ))
        mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', material='Material-Unitcell', thickness=None)
        part4 = mdb.models['Model-1'].parts['UnitCell']
        c = part4.cells
        cells = (c[0],)
        part4.SectionAssignment(region=cells, sectionName='Section-1', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

        #Loadcell_Material_Properties------------------------------------------------
        mdb.models['Model-1'].Material(name='Steel')
        mdb.models['Model-1'].materials['Steel'].Elastic(table=((210000.0, 0.3), ))
        mdb.models['Model-1'].materials['Steel'].Density(table=((7.85e-09, ), ))
        mdb.models['Model-1'].HomogeneousSolidSection(name='Section-2',
            material='Steel', thickness=None)
        session.viewports['Viewport: 1'].setValues(displayedObject=part5)
        c = part5.cells
        cells = (c[0],)
        part5.SectionAssignment(region=cells, sectionName='Section-2', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

        # #Step------------------------------------------------------------------------
        # mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial',
            # maxNumInc=100000, initialInc=0.001, minInc=1e-15, nlgeom=ON)
        mdb.models['Model-1'].ImplicitDynamicsStep(name='Step-1',
             previous='Initial', maxNumInc=100000, initialInc=0.1, minInc=1e-15)

        mdb.models['Model-1'].steps['Step-1'].setValues(timePeriod=0.1,
                initialInc=0.01)

        r1 = myassembly.instances['Part-loadcell-1'].referencePoints
        refPoints1=(r1[2], )
        myassembly.Set(referencePoints=refPoints1, name='Set-RP')

        #Surfaces---------------------------------------------------------------------

        s1 = myassembly.instances['UnitCell-1'].faces
        side1Faces1 = s1.findAt(((6.890767, 19.890798, 13.432006), ),((19.583333, 19.975, 13.483334),), ((26.18601,
            19.91439, 6.974446), ), ((13.31399, 19.91439, -6.474446), ), ((
            32.609233, 19.890798, -12.932006), ), ((19.603531, 19.875602, 19.780691), ), ((0.285759, 19.975, 0.331106), ), ((19.771102, 19.750235, 19.958278), ), ((
                39.088842, 19.975, 0.331106), ), ((19.728898, 19.750235, -19.458278),
                ), ((0.411159, 19.975, 0.168894), ), ((19.896469, 19.875602,
                -19.280691), ), ((39.21424, 19.975, 0.168894), ), ((26.075518, 19.883788, -6.406871), ), ((36.126409,
                19.727203, -15.832677), ), ((22.902999, 19.879727, -15.880535), ), ((
                35.83735, 19.943158, -3.053078), ), ((16.579026, 19.908122, -3.219471),
                ), ((3.373591, 19.727203, -15.832677), ), ((16.597348, 19.846799,
                -15.849407), ), ((6.905921, 19.724041, -6.112154), ), ((22.920974,
                19.908122, 3.719471), ), ((36.126409, 19.727203, 16.332677), ), ((
                22.902652, 19.846799, 16.349407), ), ((32.594079, 19.724041, 6.612155),
                ), ((13.424482, 19.883788, 6.906871), ), ((3.373591, 19.727203,
                16.332677), ), ((16.597001, 19.879727, 16.380535), ), ((3.66265,
                19.943158, 3.553078), ), ((26.16831, 19.521977, -6.555232), ), ((13.342082,
                19.721444, -6.58192), ), ((26.157918, 19.721444, 7.081919), ), ((
                13.33169, 19.521978, 7.055232), ), ((26.098117, 19.80799, -12.710671), ), ((13.401883,
                19.80799, -12.710671), ), ((26.098117, 19.80799, 13.210671), ), ((
                13.401883, 19.80799, 13.210671), ), ((39.107462, 19.784518, -19.409296), ), ((0.351738,
                19.746539, -19.447275), ), ((39.148262, 19.746539, 19.947275), ), ((
                0.392538, 19.784518, 19.909296), ), ((6.6, 19.975, -19.357149), ), ((6.6, 19.975,
                19.739298), ), ((0.064989, 9.975, 0.041667), ))

        myassembly.Surface(side1Faces=side1Faces1, name='Surf-1')
        s1 = myassembly.instances['Part-loadcell-1'].faces
        side1Faces1 = s1.findAt(((4.75, 9.975, -5.0), ))
        myassembly.Surface(side1Faces=side1Faces1, name='Surf-T')

        s1 = myassembly.instances['UnitCell-1'].faces
        side1Faces1 = s1.findAt(((26.108256, -19.091722, -12.651264), ), ((32.698032,
                -19.344736, -13.006969), ), ((26.216335, -19.398633, -6.35463), ), ((
                10.090001, -19.41946, -9.429875), ), ((3.770499, -19.027159,
                -16.120744), ), ((13.354364, -19.018369, -6.529406), ), ((29.409999,
                -19.41946, 9.929875), ), ((35.729501, -19.02716, 16.620744), ), ((
                26.145636, -19.018369, 7.029406), ), ((19.583333, -19.475, 13.483334),
                ), ((13.391744, -19.091722, 13.151264), ), ((6.801968, -19.344736,
                13.506969), ), ((13.283665, -19.398633, 6.85463), ), ((19.666667, -19.475, -12.666667), ),
                ((6.6, -19.475, 0.416667), ), ((26.226129, -19.448603, -6.440331), ), ((13.327378,
                -19.380892, -6.507338), ), ((26.172622, -19.380892, 7.007338), ), ((
                13.27387, -19.448603, 6.940331), ), ((26.098117, -19.30799, -12.710671), ), ((13.401883,
                -19.30799, -12.710671), ), ((26.098117, -19.30799, 13.210671), ), ((
                13.401883, -19.30799, 13.210671), ), ((19.771102, -19.250235, -19.458278), ), ((19.603531,
                -19.375602, -19.280691), ), ((19.896469, -19.375602, 19.780691), ), ((
                19.728898, -19.250235, 19.958278), ), ((39.148262, -19.246539, -19.447275), ), ((0.392538,
                -19.284518, -19.409296), ), ((39.107462, -19.284518, 19.909296), ), ((
                0.351738, -19.246539, 19.947275), ), ((26.324999, -19.475, -19.357149), ), ((6.6, -19.475,
                19.857149), ), ((19.75, -19.475, -19.282149), ), ((19.685011, -9.725, 0.041667), ))
        myassembly.Surface(side1Faces=side1Faces1, name='Surf-2')

        s1 = myassembly.instances['Part-loadcell-2'].faces
        side1Faces1 = s1.findAt(((4.75, -9.725, 5.0), ))
        myassembly.Surface(side1Faces=side1Faces1, name='Surf-B')

        #Interactions---------------------------------------------------------------
        mdb.models['Model-1'].ContactProperty('IntProp-1')
        mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
            formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
            pressureDependency=OFF, temperatureDependency=OFF, dependencies=0,
            table=((0.2, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION,
            fraction=0.005, elasticSlipStiffness=None)
        mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=ON,
            constraintEnforcementMethod=DEFAULT)

        mdb.models['Model-1'].ContactStd(name='Int-1', createStepName='Initial')
        mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
              stepName='Initial', useAllstar=ON)
        mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
              stepName='Initial', assignments=((GLOBAL, SELF, 'IntProp-1'), ))

        c1 = myassembly.instances['Part-loadcell-1'].cells
        cells1 = (c1[0], )
        r1 = myassembly.instances['Part-loadcell-1'].referencePoints
        refPoints1=(r1[2], )
        mdb.models['Model-1'].RigidBody(name='Constraint-1', refPointRegion=refPoints1,
           bodyRegion=cells1)

        c1 = myassembly.instances['Part-loadcell-2'].cells
        cells1 = (c1[0], )
        r1 = myassembly.instances['Part-loadcell-2'].referencePoints
        refPoints1=(r1[2], )
        mdb.models['Model-1'].RigidBody(name='Constraint-2', refPointRegion=refPoints1,
            bodyRegion=cells1)

        #Surface Contact-------------------------------------------------------------------

        region1=myassembly.surfaces['Surf-T']
        region2=myassembly.surfaces['Surf-1']

        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-2',
                createStepName='Initial', main=region1, secondary=region2,
                sliding=FINITE, thickness=ON, interactionProperty='IntProp-1',
                adjustMethod=NONE, initialClearance=OMIT, datumAxis=None,
                clearanceRegion=None)

        region1=myassembly.surfaces['Surf-B']
        region2=myassembly.surfaces['Surf-2']
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-3',
            createStepName='Initial', main=region1, secondary=region2,
            sliding=FINITE, thickness=ON, interactionProperty='IntProp-1',
            adjustMethod=NONE, initialClearance=OMIT, datumAxis=None,
            clearanceRegion=None)

        #Field_output--------------------------------------------------------------
        mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(
           frequency=1)

        #HistoryOutput-------------------------------------------------------------
        regionDef=mdb.models['Model-1'].rootAssembly.sets['Set-RP']
        mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(variables=(
                'U2', 'RF2'), region=regionDef, sectionPoints=DEFAULT, rebar=EXCLUDE)
        mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(
              frequency=1)

        # #BoundaryCondition---------------------------------------------------------

        mdb.models['Model-1'].TabularAmplitude(name='Amp-1', timeSpan=STEP,
            smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
        r1 = myassembly.instances['Part-loadcell-1'].referencePoints
        refPoints1=(r1[2], )
        mdb.models['Model-1'].DisplacementBC(name='BC-2', createStepName='Step-1',
            region=refPoints1, u1=0.0, u2=-1, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0,
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='',
            localCsys=None)

        r1 = myassembly.instances['Part-loadcell-2'].referencePoints
        refPoints1=(r1[2], )
        mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Step-1',
              region=refPoints1, localCsys=None)
        #Mesh----------------------------------------------------------------------
        c = part4.cells[0]
        part4.setMeshControls(regions=(c,), elemShape=TET, technique=FREE)
        elemType1 = mesh.ElemType(elemCode=C3D8R)
        elemType2 = mesh.ElemType(elemCode=C3D6)
        elemType3 = mesh.ElemType(elemCode=C3D4,
                secondOrderAccuracy=OFF, distortionControl=DEFAULT)
        part4.setElementType(regions=(c,), elemTypes=(elemType1, elemType2,
              elemType3))
        part4.seedPart(size=1.2, deviationFactor=0.1, minSizeFactor=0.1)
        part4.generateMesh()

        part5.seedPart(size=1000.0, deviationFactor=0.1, minSizeFactor=0.1)
        part5.generateMesh()

        # Calculate the mass of the part------------------------------------------
        mass_props = part4.getMassProperties()
        mass = mass_props['mass'] * 0.001
        print("Mass of the unit cell: {} Kg".format(mass))

        # Write the mass to the results file---------------------------------------
        ff.write("{}, {}\n".format(fingerprint, mass))
        ff.flush()

        # Save this structure's CAE file--------------------------------------------
        mdb.saveAs(pathName=os.path.join(OUTPUT_DIR, fingerprint))

        # Close and reset the model database before the next iteration-------------
        mdb.close()
        mdb = Mdb()

    except Exception as e:
        print("Error occurred while creating {}: {}".format(fingerprint, e))
        try:
            mdb.close()
        except Exception:
            pass
        mdb = Mdb()
        continue

ff.close()
