# This is simpel scipt to convert cgns mesh files to gmsh msh files
# Usage: 
# python cgnsToMsh --cgnsFile </full/path/to/cgnsFile> --mshFile </full/path/to/mshFile>
# To check help type: python cgnsToMsh.py --h


from paraview.simple import *
import numpy as np
import os

"""Function that counts number of files
with specific extension in directory """
def directory(path,extension):
  list_dir = []
  list_dir = os.listdir(path)
  count = 0
  for file in list_dir:
    if file.endswith(extension): # eg: '.txt'
      count += 1
  return count

"""Choose input/output directory and filenames"""
vtu_input_directory = "."
cgns_output_directory = "."
input_filename_root = "."
output_filename_root = "."

""" Create output directory """
os.system('mkdir '+ cgns_output_directory)

file_dir  = os.listdir('.')
file_data = [i for i in file_dir if i.endswith('.vtu')] 
files_vtu = np.array(file_data)

# for file in files_vtu:
# 	reader = simple.OpenDataFile(file)
# 	print(file[:file.rfind(".")]+'.cgns')
# 	writer = simple.CreateWriter(file[:file.rfind(".")]+'.cgns', reader)

"""Start Loop over all files """
number_of_vtu = directory(vtu_input_directory,'.vtu')
for index in range(1,number_of_vtu):
    in_filename = files_vtu[index]; 
    out_filename = in_filename[:in_filename.rfind(".")]+'.cgns' 
    loadfile = in_filename
    writefile = out_filename
    r = XMLPartitionedUnstructuredGridReader( FileName=loadfile)
    writer = CreateWriter(writefile,r)
    writer.FieldAssociation = "Points"
    writer.UpdatePipeline()



# import os
# import numpy as np
# import sys
# from paraview import simple
# import argparse
# parser = argparse.ArgumentParser(prog='CGNS_Converter',
# 																description='Convert CGNS mesh files to GMSH v2.2',
#                     						epilog='--------------------------------')

# parser.add_argument('--outDir', default='.', help='mesh file to write')
# args = parser.parse_args()


# file_dir  = os.listdir(args.outDir)
# file_data = [i for i in file_dir if i.endswith('.vtu')] 
# files_vtu = np.array(file_data)

# for file in files_vtu:
# 	reader = simple.OpenDataFile(file)
# 	print(file[:file.rfind(".")]+'.cgns')
# 	writer = simple.CreateWriter(file[:file.rfind(".")]+'.cgns', reader)


# for file in os.listdir(args.outDir):
#     if file.endswith(".vtu"):
#         print(os.path.join(args.outDir, file))




# reader = simple.OpenDataFile("cns_0000_0000.vtu")
# writer = simple.CreateWriter("cns_000_0000.cgns", reader)
# writer.WriteAllTimeSteps = 1
# writer.FieldAssociation = "Points"
# writer.UpdatePipeline()
