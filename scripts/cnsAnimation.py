#### cnsAnimation.py
####   usage: pvbatch pvbatch cnsAnimation.py 1 ../examples/cnsTri2D/foo foo


#### import the simple module from the paraview
from paraview.simple import *
import glob 
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

Ndatasets = int(sys.argv[1])
datasetIn = sys.argv[2]
imageFilesOut = sys.argv[3]

dataset = range(0,Ndatasets)
datasetDisplay = range(0,Ndatasets)
for n in range(0,Ndatasets):
  Nfiles = len(glob.glob(datasetIn+'_%04d_*.vtu' % n))
  files = range(0,Nfiles)
  for m in range(0,Nfiles):
    files[m] = glob.glob(datasetIn+'_%04d_%04d.vtu' % (n,m))[0]

  # create a new 'XML Unstructured Grid Reader'
  dataset[n] = XMLUnstructuredGridReader(FileName=files)
  dataset[n].PointArrayStatus = ['Density', 'Velocity', 'Vorticity']

  # set active source
  #SetActiveSource(dataset[n])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

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
groupDatasets1Display.OSPRayScaleArray = 'Density'
groupDatasets1Display.OSPRayScaleFunction = 'PiecewiseFunction'
groupDatasets1Display.SelectOrientationVectors = 'Density'
groupDatasets1Display.ScaleFactor = 1.0527610063552857
groupDatasets1Display.SelectScaleArray = 'Density'
groupDatasets1Display.GlyphType = 'Arrow'
groupDatasets1Display.GlyphTableIndexArray = 'Density'
groupDatasets1Display.DataAxesGrid = 'GridAxesRepresentation'
groupDatasets1Display.PolarAxes = 'PolarAxesRepresentation'
groupDatasets1Display.ScalarOpacityUnitDistance = 0.09258228982369943
groupDatasets1Display.GaussianRadius = 0.5263805031776428
groupDatasets1Display.SetScaleArray = ['POINTS', 'Density']
groupDatasets1Display.ScaleTransferFunction = 'PiecewiseFunction'
groupDatasets1Display.OpacityArray = ['POINTS', 'Density']
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
vorticityLUT.RescaleTransferFunction(-4.0, 4.0)

# get opacity transfer function/opacity map for 'Vorticity'
vorticityPWF = GetOpacityTransferFunction('Vorticity')

# Rescale transfer function
vorticityPWF.RescaleTransferFunction(-4.0, 4.0)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0
groupDatasets1Display.SetScalarBarVisibility(renderView1, False)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [3.8648142362157496, -0.12321016819782263, 25.0]
renderView1.CameraFocalPoint = [3.8648142362157496, -0.12321016819782263, 0.0]
renderView1.CameraParallelScale = 6.210760565455133

# save animation
#SaveAnimation(imageFilesOut+'.avi', renderView1, ImageResolution=[1596, 800], FrameRate=15, FrameWindow=[0, len(files)])
SaveAnimation(imageFilesOut+'.png', renderView1, ImageResolution=[1855, 1163],
    FrameWindow=[0, 10])

#### uncomment the following to render all views
#RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).

# save screenshot
#SaveScreenshot(imageFilesOut+'.png', renderView1, ImageResolution=[1599, 803])