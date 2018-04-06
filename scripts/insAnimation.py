#### insAnimation.py
####   usage: pvbatch pvbatch insAnimation.py 1 ../examples/cnsTri2D/foo ./ foo


#### import the simple module from the paraview
from paraview.simple import *
import glob 
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

Ndatasets = int(sys.argv[1])
datasetIn = sys.argv[2]
directoryOut = sys.argv[3]
imageFilesOut = sys.argv[4]

dataset = range(0,Ndatasets)
datasetDisplay = range(0,Ndatasets)
for n in range(0,Ndatasets):
	files = glob.glob(datasetIn+'_%04d_*.vtu' % n)

	# create a new 'XML Unstructured Grid Reader'
	dataset[n] = XMLUnstructuredGridReader(FileName=files)
	dataset[n].PointArrayStatus = ['Pressure', 'Divergence', 'Vorticity', 'Velocity']

	# get animation scene
	animationScene1 = GetAnimationScene()

	# update animation scene based on data timesteps
	animationScene1.UpdateAnimationUsingDataTimeSteps()

	# set active source
	SetActiveSource(dataset[n])


# create a new 'Group Datasets'
groupDatasets1 = GroupDatasets(Input=dataset)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1599, 803]

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
# renderView1.InteractionMode = '2D'
# renderView1.CameraPosition = [-0.39905500411987305, 0.0, 10000.0]
# renderView1.CameraFocalPoint = [-0.39905500411987305, 0.0, 0.0]

# # show data in view
groupDatasets1Display = Show(groupDatasets1, renderView1)
# trace defaults for the display properties.
groupDatasets1Display.Representation = 'Surface'
groupDatasets1Display.ColorArrayName = [None, '']
groupDatasets1Display.OSPRayScaleArray = 'Divergence'
groupDatasets1Display.OSPRayScaleFunction = 'PiecewiseFunction'
groupDatasets1Display.SelectOrientationVectors = 'Divergence'
groupDatasets1Display.ScaleFactor = 1.8
groupDatasets1Display.SelectScaleArray = 'Divergence'
groupDatasets1Display.GlyphType = 'Arrow'
groupDatasets1Display.GlyphTableIndexArray = 'Divergence'
groupDatasets1Display.DataAxesGrid = 'GridAxesRepresentation'
groupDatasets1Display.PolarAxes = 'PolarAxesRepresentation'
groupDatasets1Display.ScalarOpacityUnitDistance = 0.1199621937768176
groupDatasets1Display.GaussianRadius = 0.9
groupDatasets1Display.SetScaleArray = ['POINTS', 'Divergence']
groupDatasets1Display.ScaleTransferFunction = 'PiecewiseFunction'
groupDatasets1Display.OpacityArray = ['POINTS', 'Divergence']
groupDatasets1Display.OpacityTransferFunction = 'PiecewiseFunction'

# # hide data in view
# #Hide(foo_0001_000, renderView1)

# # hide data in view
# #Hide(foo_0000_000, renderView1)

# # update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(groupDatasets1Display, ('POINTS', 'Vorticity'))

# Hide the scalar bar for this color map if no visible data is colored by it.
#HideScalarBarIfNotNeeded(vtkBlockColorsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
groupDatasets1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
#groupDatasets1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Vorticity'
vorticityLUT = GetColorTransferFunction('Vorticity')

# Rescale transfer function
vorticityLUT.RescaleTransferFunction(-10.0, 10.0)

# get opacity transfer function/opacity map for 'Vorticity'
vorticityPWF = GetOpacityTransferFunction('Vorticity')

# Rescale transfer function
vorticityPWF.RescaleTransferFunction(-10.0, 10.0)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [2.742803327974873, 0.05383143533708823, 8.0]
renderView1.CameraFocalPoint = [2.742803327974873, 0.05383143533708823, 0.0]
renderView1.CameraParallelScale = 2.409167061817433

# save animation
SaveAnimation(imageFilesOut+'.avi', renderView1, ImageResolution=[1596, 800], FrameRate=15, FrameWindow=[0, len(files)])

#### uncomment the following to render all views
RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

# save screenshot
#SaveScreenshot(imageFilesOut+'.png', renderView1, ImageResolution=[1599, 803])