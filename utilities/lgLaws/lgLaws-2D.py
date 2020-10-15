#!pvpython
# Authors: Octavio Garcia Aguirre, David Harispe

# To "install" this program, locate it inside the "bin" folder in Paraview's
# installation. You WILL need this bin folder added to cmd's (or whatever your
# command line is) PATH enviroment variable
# You will also probably need to change this file in specific ways for
# different cases
# - Change the plot overline origin and endpoint depending on your domain (line 105 - 110).
# - Increase the plot overline resolution or increase the "eps" to detect
#   very small bands (line 199).


# trace generated using paraview version 5.8.0-RC1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import csv
import math
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
liesegangfoam = OpenFOAMReader(FileName='Liesegang.foam')
liesegangfoam.MeshRegions = ['internalMesh']
liesegangfoam.CellArrays = ['A', 'A_0', 'A_0_0', 'B', 'B_0', 'B_0_0', 'C', 'C_0', 'C_0_0', 'D', 'D_0', 'D_0_0', 'ddt0(A)', 'ddt0(B)', 'ddt0(C)', 'ddt0(D)']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1602, 757]

# get layout
layout1 = GetLayout()

# show data in view
liesegangfoamDisplay = Show(liesegangfoam, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
liesegangfoamDisplay.Representation = 'Surface'
liesegangfoamDisplay.ColorArrayName = [None, '']
liesegangfoamDisplay.OSPRayScaleArray = 'A'
liesegangfoamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
liesegangfoamDisplay.SelectOrientationVectors = 'A'
liesegangfoamDisplay.ScaleFactor = 0.0019999999552965165
liesegangfoamDisplay.SelectScaleArray = 'A'
liesegangfoamDisplay.GlyphType = 'Arrow'
liesegangfoamDisplay.GlyphTableIndexArray = 'A'
liesegangfoamDisplay.GaussianRadius = 9.999999776482583e-05
liesegangfoamDisplay.SetScaleArray = ['POINTS', 'A']
liesegangfoamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
liesegangfoamDisplay.OpacityArray = ['POINTS', 'A']
liesegangfoamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
liesegangfoamDisplay.DataAxesGrid = 'GridAxesRepresentation'
liesegangfoamDisplay.PolarAxes = 'PolarAxesRepresentation'
liesegangfoamDisplay.ScalarOpacityUnitDistance = 0.0008100527439950064
liesegangfoamDisplay.ExtractedBlockIndex = 1

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
liesegangfoamDisplay.ScaleTransferFunction.Points = [1.8688099961089299e-19, 0.0, 0.5, 0.0, 500.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
liesegangfoamDisplay.OpacityTransferFunction.Points = [1.8688099961089299e-19, 0.0, 0.5, 0.0, 500.0, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(liesegangfoamDisplay, ('POINTS', 'D'))

# rescale color and/or opacity maps used to include current data range
liesegangfoamDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
liesegangfoamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'D'
dLUT = GetColorTransferFunction('D')

# get opacity transfer function/opacity map for 'D'
dPWF = GetOpacityTransferFunction('D')

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(Input=liesegangfoam,
    Source='High Resolution Line Source')

# init the 'High Resolution Line Source' selected for 'Source'
origin = [0, 0, 0.001]
plotOverLine1.Source.Point1 = origin
plotOverLine1.Source.Point2 = [0.02, 0.02, 0.001]
plotOverLine1.Source.Resolution = 10000

# show data in view
plotOverLine1Display = Show(plotOverLine1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plotOverLine1Display.Representation = 'Surface'
plotOverLine1Display.ColorArrayName = ['POINTS', 'D']
plotOverLine1Display.LookupTable = dLUT
plotOverLine1Display.OSPRayScaleArray = 'A'
plotOverLine1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plotOverLine1Display.SelectOrientationVectors = 'A'
plotOverLine1Display.ScaleFactor = 0.0019999999552965165
plotOverLine1Display.SelectScaleArray = 'A'
plotOverLine1Display.GlyphType = 'Arrow'
plotOverLine1Display.GlyphTableIndexArray = 'A'
plotOverLine1Display.GaussianRadius = 9.999999776482583e-05
plotOverLine1Display.SetScaleArray = ['POINTS', 'A']
plotOverLine1Display.ScaleTransferFunction = 'PiecewiseFunction'
plotOverLine1Display.OpacityArray = ['POINTS', 'A']
plotOverLine1Display.OpacityTransferFunction = 'PiecewiseFunction'
plotOverLine1Display.DataAxesGrid = 'GridAxesRepresentation'
plotOverLine1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plotOverLine1Display.ScaleTransferFunction.Points = [1.8688099961089299e-19, 0.0, 0.5, 0.0, 499.1441650390625, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plotOverLine1Display.OpacityTransferFunction.Points = [1.8688099961089299e-19, 0.0, 0.5, 0.0, 499.1441650390625, 1.0, 0.5, 0.0]

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
# uncomment following to set a specific view size
# lineChartView1.ViewSize = [400, 400]

# show data in view
plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1, 'XYChartRepresentation')

# trace defaults for the display properties.
plotOverLine1Display_1.CompositeDataSetIndex = [0]
plotOverLine1Display_1.UseIndexForXAxis = 0
plotOverLine1Display_1.XArrayName = 'arc_length'
plotOverLine1Display_1.SeriesVisibility = ['D']
plotOverLine1Display_1.SeriesLabel = ['A', 'A', 'A_0', 'A_0', 'A_0_0', 'A_0_0', 'arc_length', 'arc_length', 'B', 'B', 'B_0', 'B_0', 'B_0_0', 'B_0_0', 'C', 'C', 'C_0', 'C_0', 'C_0_0', 'C_0_0', 'D', 'D', 'D_0', 'D_0', 'D_0_0', 'D_0_0', 'ddt0(A)', 'ddt0(A)', 'ddt0(B)', 'ddt0(B)', 'ddt0(C)', 'ddt0(C)', 'ddt0(D)', 'ddt0(D)', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotOverLine1Display_1.SeriesColor = ['A', '0', '0', '0', 'A_0', '0.89', '0.1', '0.11', 'A_0_0', '0.22', '0.49', '0.72', 'arc_length', '0.3', '0.69', '0.29', 'B', '0.6', '0.31', '0.64', 'B_0', '1', '0.5', '0', 'B_0_0', '0.65', '0.34', '0.16', 'C', '0', '0', '0', 'C_0', '0.89', '0.1', '0.11', 'C_0_0', '0.22', '0.49', '0.72', 'D', '0.3', '0.69', '0.29', 'D_0', '0.6', '0.31', '0.64', 'D_0_0', '1', '0.5', '0', 'ddt0(A)', '0.65', '0.34', '0.16', 'ddt0(B)', '0', '0', '0', 'ddt0(C)', '0.89', '0.1', '0.11', 'ddt0(D)', '0.22', '0.49', '0.72', 'vtkValidPointMask', '0.3', '0.69', '0.29', 'Points_X', '0.6', '0.31', '0.64', 'Points_Y', '1', '0.5', '0', 'Points_Z', '0.65', '0.34', '0.16', 'Points_Magnitude', '0', '0', '0']
plotOverLine1Display_1.SeriesPlotCorner = ['A', '0', 'A_0', '0', 'A_0_0', '0', 'arc_length', '0', 'B', '0', 'B_0', '0', 'B_0_0', '0', 'C', '0', 'C_0', '0', 'C_0_0', '0', 'D', '0', 'D_0', '0', 'D_0_0', '0', 'ddt0(A)', '0', 'ddt0(B)', '0', 'ddt0(C)', '0', 'ddt0(D)', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine1Display_1.SeriesLabelPrefix = ''
plotOverLine1Display_1.SeriesLineStyle = ['A', '1', 'A_0', '1', 'A_0_0', '1', 'arc_length', '1', 'B', '1', 'B_0', '1', 'B_0_0', '1', 'C', '1', 'C_0', '1', 'C_0_0', '1', 'D', '1', 'D_0', '1', 'D_0_0', '1', 'ddt0(A)', '1', 'ddt0(B)', '1', 'ddt0(C)', '1', 'ddt0(D)', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotOverLine1Display_1.SeriesLineThickness = ['A', '2', 'A_0', '2', 'A_0_0', '2', 'arc_length', '2', 'B', '2', 'B_0', '2', 'B_0_0', '2', 'C', '2', 'C_0', '2', 'C_0_0', '2', 'D', '2', 'D_0', '2', 'D_0_0', '2', 'ddt0(A)', '2', 'ddt0(B)', '2', 'ddt0(C)', '2', 'ddt0(D)', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotOverLine1Display_1.SeriesMarkerStyle = ['A', '0', 'A_0', '0', 'A_0_0', '0', 'arc_length', '0', 'B', '0', 'B_0', '0', 'B_0_0', '0', 'C', '0', 'C_0', '0', 'C_0_0', '0', 'D', '0', 'D_0', '0', 'D_0_0', '0', 'ddt0(A)', '0', 'ddt0(B)', '0', 'ddt0(C)', '0', 'ddt0(D)', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotOverLine1Display_1.SeriesMarkerSize = ['A', '4', 'A_0', '4', 'A_0_0', '4', 'arc_length', '4', 'B', '4', 'B_0', '4', 'B_0_0', '4', 'C', '4', 'C_0', '4', 'C_0_0', '4', 'D', '4', 'D_0', '4', 'D_0_0', '4', 'ddt0(A)', '4', 'ddt0(B)', '4', 'ddt0(C)', '4', 'ddt0(D)', '4', 'vtkValidPointMask', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Points_Magnitude', '4']

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView1, layout=layout1, hint=0)

animationScene1.GoToLast()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.009999999776482582, 0.009999999776482582, 0.05517515492749509]
renderView1.CameraFocalPoint = [0.009999999776482582, 0.009999999776482582, 0.0005000000237487257]
renderView1.CameraParallelScale = 0.01415097138302004

# Our code
writer = CreateWriter('plotoverline.csv',plotOverLine1)
writer.FieldAssociation = "Point Data"
writer.UpdatePipeline()

D = []
pos = []
with open('plotoverline.csv') as csvfile:
    csvread = csv.DictReader(csvfile)
    for row in csvread:
        val = row['D']
        if val != 'nan':
            D.append(float(val))
            x = float(row['Points:0'])
            y = float(row['Points:1'])
            z = float(row['Points:2'])
            pos.append((x,y,z))

maxval = -float('inf')
for val in D:
    if val > maxval:
        maxval = val

print('Calculated rho',maxval)
saturateds = []
vals = []
eps = 0.1
#print('Calculated rho',maxval)
for idx,val in enumerate(D):
    if val >= (maxval-eps):
        saturateds.append(idx)
        vals.append(val)



#print('Saturateds',saturateds)
#print('vals',vals)
#unrreal base case just to avoid unnecesary IFs
streaks = [[-5,-5]]
for idx in saturateds:
    last = streaks[-1]
    if last[1] == (idx-1):
        last[1] = idx
    else:
        streaks.append([idx,idx])

        
#Remove base case
streaks = streaks[1:]
#Remove "0" width streaks
#(causes division by 0 when checking width law...)
streaks = [x for x in streaks if x[0]!=x[1]]

def distance(a,b=(0,0,0)):
    x = a[0] - b[0]
    y = a[1] - b[1]
    z = a[2] - b[2]
    return math.sqrt(x*x+y*y+z*z)

print('{:>12},{:>12},{:>12},{:>12},{:>12}'.format('Begin','End','Width','Width law','Space law'))
obegin = float('nan')
oend   = float('nan')
owidth = float('nan')
maxwidth = -float('inf')
minwidth = float('inf')
maxspace = -float('inf')
minspace = float('inf')
avgwidth = 0
avgspace = 0
count = 0
for idx,s in enumerate(streaks):
    begin = distance(pos[s[0]],origin)
    end   = distance(pos[s[1]],origin)
    width = distance(pos[s[0]],pos[s[1]])
    widthlaw = width/owidth
    spacelaw = end/oend
    print(f'{begin:-12.6f},{end:-12.6f},{width:-12.6f},{widthlaw:-12.6f},{spacelaw:-12.6f}')
    obegin = begin
    oend   = end
    owidth = width
    if idx == 0:#ignore NaN
        continue
    avgwidth += widthlaw
    avgspace += spacelaw
    count+=1
    if(widthlaw > maxwidth):
        maxwidth = widthlaw
    if(widthlaw < minwidth):
        minwidth = widthlaw
    if(spacelaw > maxspace):
        maxspace = spacelaw
    if(spacelaw < minspace):
        minspace = spacelaw

avgwidth /= count
avgspace /= count
print('{:>12},{:>12},{:>12},{:-12.6f},{:-12.6f}'.format('Min','Min','Min',minwidth,minspace))
print('{:>12},{:>12},{:>12},{:-12.6f},{:-12.6f}'.format('Max','Max','Max',maxwidth,maxspace))
print('{:>12},{:>12},{:>12},{:-12.6f},{:-12.6f}'.format('Average','Average','Average',avgwidth,avgspace))
print('Fin primer banda',distance(pos[streaks[0][1]],origin))
print('Fin ultima banda',distance(pos[streaks[-1][1]],origin))
print('Numero de bandas:',len(streaks))

#### uncomment the following to render all views
#RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
Interact()
