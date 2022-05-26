# script-version: 2.0
# Catalyst state generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1602, 755]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [15.999999664485241, 16.000000135640732, 26.041901441337483]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [36.66411697921538, 48.55513931624985, 71.55408524634686]
renderView1.CameraFocalPoint = [15.99999966448525, 16.000000135640718, 26.041901441337455]
renderView1.CameraViewUp = [-0.32121385382146983, 0.8331843190266527, -0.4501394790964824]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 15.438727875657248
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #2'
layout2 = CreateLayout(name='Layout #2')
layout2.AssignView(0, renderView1)
layout2.SetSize(1602, 755)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
fluidpvd = PVDReader(registrationName='fluid', FileName='/home/jifutan/insituBloodFlow/BloodFlow/examples/singleCell/vtk_output/fluid.pvd')
fluidpvd.CellArrays = ['vtkTestType']
fluidpvd.PointArrays = ['velocity', 'vorticity', 'velocityNorm']

# create a new 'PVD Reader'
cellspvd = PVDReader(registrationName='cells', FileName='/home/jifutan/insituBloodFlow/BloodFlow/examples/singleCell/vtk_output/cells.pvd')

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from cellspvd
cellspvdDisplay = Show(cellspvd, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.InterpretValuesAsCategories = 1
vtkBlockColorsLUT.AnnotationsInitialized = 1
vtkBlockColorsLUT.Annotations = ['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11']
vtkBlockColorsLUT.ActiveAnnotatedValues = ['0', '2', '4', '6', '1', '3', '5', '7', '8', '9', '10', '11']
vtkBlockColorsLUT.IndexedColors = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5]

# trace defaults for the display properties.
cellspvdDisplay.Representation = 'Surface'
cellspvdDisplay.ColorArrayName = ['FIELD', 'vtkBlockColors']
cellspvdDisplay.LookupTable = vtkBlockColorsLUT
cellspvdDisplay.SelectTCoordArray = 'None'
cellspvdDisplay.SelectNormalArray = 'None'
cellspvdDisplay.SelectTangentArray = 'None'
cellspvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
cellspvdDisplay.SelectOrientationVectors = 'None'
cellspvdDisplay.ScaleFactor = 1.8000001233250211
cellspvdDisplay.SelectScaleArray = 'None'
cellspvdDisplay.GlyphType = 'Arrow'
cellspvdDisplay.GlyphTableIndexArray = 'None'
cellspvdDisplay.GaussianRadius = 0.09000000616625105
cellspvdDisplay.SetScaleArray = [None, '']
cellspvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
cellspvdDisplay.OpacityArray = [None, '']
cellspvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
cellspvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
cellspvdDisplay.PolarAxes = 'PolarAxesRepresentation'

# show data from fluidpvd
fluidpvdDisplay = Show(fluidpvd, renderView1, 'UniformGridRepresentation')

# get opacity transfer function/opacity map for 'vtkBlockColors'
vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

# trace defaults for the display properties.
fluidpvdDisplay.Representation = 'Outline'
fluidpvdDisplay.ColorArrayName = ['FIELD', 'vtkBlockColors']
fluidpvdDisplay.LookupTable = vtkBlockColorsLUT
fluidpvdDisplay.SelectTCoordArray = 'None'
fluidpvdDisplay.SelectNormalArray = 'None'
fluidpvdDisplay.SelectTangentArray = 'None'
fluidpvdDisplay.OSPRayScaleArray = 'velocity'
fluidpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
fluidpvdDisplay.SelectOrientationVectors = 'None'
fluidpvdDisplay.ScaleFactor = 8.6
fluidpvdDisplay.SelectScaleArray = 'None'
fluidpvdDisplay.GlyphType = 'Arrow'
fluidpvdDisplay.GlyphTableIndexArray = 'None'
fluidpvdDisplay.GaussianRadius = 0.43
fluidpvdDisplay.SetScaleArray = ['POINTS', 'velocity']
fluidpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
fluidpvdDisplay.OpacityArray = ['POINTS', 'velocity']
fluidpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
fluidpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
fluidpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
fluidpvdDisplay.ScalarOpacityUnitDistance = 1.882751655065373
fluidpvdDisplay.ScalarOpacityFunction = vtkBlockColorsPWF
fluidpvdDisplay.OpacityArrayName = ['POINTS', 'velocity']
fluidpvdDisplay.SliceFunction = 'Plane'
fluidpvdDisplay.Slice = 43

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
fluidpvdDisplay.ScaleTransferFunction.Points = [-3.0756343308182305e-06, 0.0, 0.5, 0.0, 3.075634330818228e-06, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
fluidpvdDisplay.OpacityTransferFunction.Points = [-3.0756343308182305e-06, 0.0, 0.5, 0.0, 3.075634330818228e-06, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
fluidpvdDisplay.SliceFunction.Origin = [15.0, 15.0, 40.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for vtkBlockColorsLUT in view renderView1
vtkBlockColorsLUTColorBar = GetScalarBar(vtkBlockColorsLUT, renderView1)
vtkBlockColorsLUTColorBar.Title = 'vtkBlockColors'
vtkBlockColorsLUTColorBar.ComponentTitle = ''

# set color bar visibility
vtkBlockColorsLUTColorBar.Visibility = 1

# show color legend
cellspvdDisplay.SetScalarBarVisibility(renderView1, True)

# show color legend
fluidpvdDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
# pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
# pNG1.Writer.FileName = 'RenderView1_%.6ts%cm.png'
# pNG1.Writer.ImageResolution = [1602, 755]
# pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
# SetActiveSource(pNG1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.EnableCatalystLive = 1
options.CatalystLiveURL = 'localhost:11111'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
