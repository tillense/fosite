import vtk
import os

def compute_spect(basename, sfiles, nfiles):
   print "VTS files are joined together via PVTS file into one global file."
   for i in range(sfiles,nfiles+1):
	if (i < 10):
		ext = "_000" + str(i)
        elif (i < 100):
		ext = "_00" + str(i)
        elif (i < 1000):
		ext = "_0" + str(i)
        else:
		ext = "_" + str(i);

	filename = basename + ext + ".pvts"
	if not(os.path.isfile(filename)):
		print "ERROR in file: ", filename
		print "file not found!"
		break
	reader = vtk.vtkXMLPStructuredGridReader()
	reader.SetFileName(filename)
        reader.Update()
	mesh = reader.GetOutput()

	#nx,ny,nz =  mesh.GetDimensions()
	#print nx, ny, nz
	print "processing file: ", filename


	writer = vtk.vtkXMLStructuredGridWriter()
	writer.SetInputData(mesh)
	ofilename = "merged" + basename + ext + ".vts"
	writer.SetFileName(ofilename)
	writer.SetDataModeToAppended()
	writer.Write()
   print "finished!"


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description='Merge many vts to one global file.')
	parser.add_argument('basename', metavar='basename', type=str, help="basename of input file")
	parser.add_argument('--first','-f',type=int, default=0, help="first index")
	parser.add_argument('--last','-l',type=int,default=10000, help="last index")

	args = parser.parse_args()
	compute_spect(args.basename, args.first, args.last)

