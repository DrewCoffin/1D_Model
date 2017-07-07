/coyote                                                   January 2011

The procedures in this directory are taken from the IDL program library (the
"Coyote library") developed by David Fanning ( http://www.dfanning.com/ ).  
They are duplicated here because they are used by other procedures in the
Astronomy library.    Users who have already have the Coyote library installed
can delete this directory.     Please inform  wayne.landsman@nasa.gov if there
are any discrepancies between these procedures and those in the current Coyote
library.  


Contents of /pro/coyote

ASPECT() - Calculate position coordinates to return a plot with given aspect ratio
CENTERTLB - Position a widget program on the display at an arbitrary location.
COLOR24() - Convert a RGB color triple into the equivalent 24-bit long integer.
COLORSAREIDENTICAL() - Returns a 1 if the two input colors refer to the same color
DECOMPOSEDCOLOR = Determine if current graphics device is using color decomposition
ERROR_MESSAGE() - Device independent error messaging function
FSC_COLOR() - Obtain drawing colors by name 
FSC_COLORFILL - Wrapper to PolyFill to fill a polygon 
FSC_CONTOUR - Wrapper to Contour to with extra features 
FSC_ERASE - Wrapper to ERASE to erase a graphics window with a particular color.
FSC_PLOT - Wrapper to PLOT that provide device independence
FSC_PLOTS - Wrapper to PLOTS that provide device independence
FSC_TEXT - Wrapper to XYOUTS that provide device independence
FSC_WINDOW - Implement a "smart" resizeable graphics window
GETDECOMPOSEDSTATE() - Get the color decomposition state of the current device. 
PICKCOLORNAME - Provide a blocking widget for selecting a color name
PROGRAMROOTDIR - Portable way of finding root directory of program distribution
PS_BACKGROUND - Device-independent way to set the background color in  PostScript
PSWINDOW - Calculate size of a PostScript window with specified aspect ratio
SETDECOMPOSEDSTATE - Set the color decomposition state of the current  device
SETDEFAULTVALUE  - Set default values for positional and keyword arguments
SYMCAT() - Provide plotting symbols not in the standard PSYM definitions.
TVREAD() - Screen dumps with proper color decomposition
XCOLORS - Enhanced version of XLOADCT to interactively load color tables
