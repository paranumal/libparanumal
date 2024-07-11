# This is simpel scipt to convert cgns mesh files to gmsh msh files
# Usage: 
# python cgnsToMsh --cgnsFile </full/path/to/cgnsFile> --mshFile </full/path/to/mshFile>
# To check help type: python cgnsToMsh.py --h
import gmsh 
import argparse
import os
import numpy as np
import sys

gmsh.initialize()
parser = argparse.ArgumentParser(prog='CGNS_Converter',
																description='Convert CGNS mesh files to GMSH v2.2',
                    						epilog='--------------------------------')

parser.add_argument('--cgnsFile', help='cgns file to read')
parser.add_argument('--mshFile', default='out2.msh', help='mesh file to write')
parser.add_argument('--outDir', default='./outMshDir', help='mesh file directory to write')
args = parser.parse_args()

# Create a folder for new mesh file
os.system('mkdir '+ args.outDir)

# gmsh handles the correct file extensions automatically 
gmsh.open(args.cgnsFile)
gmsh.option.setNumber("Mesh.MshFileVersion",2.2)   
gmsh.write(args.outDir+"/"+ args.mshFile)

gmsh.finalize()
