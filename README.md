# ANUGA-erosion
(Experimental) Operators for adding Erosion to ANUGA hydro, and scenarios for erosion in natural dam-breaching floods on Mars.

<b>Minor bug fixes, April 2021.</b>

<b>Last major update, Oct. 4, 2019</b>

<b>1) Notes on installing ANUGA hydro:</b><br>
See https://github.com/GeoscienceAustralia/anuga_core/blob/master/INSTALL.rst<br>

<p>I've recently had consistent success on  multiple Windows and Linux systems installing the ANUGA dependencies via the conda/miniconda/anaconda, following the instructions.  Create an ANUGA conda environment (make sure to choose python=2.7), using the conda instructions in the INSTALL file above to get the packages that ANUGA needs. Then, when I want to run code that calls anuga, I use "conda activate anuga".  {Could be source activate anuga, depending on vagaries of your system setup}.  If all is well, the full ANUGA package will  import in python. 

<b>2) Using the New Operators:</b><br>
<p>Warning 1: Obviously, these erosion operators (and other modules) are experimental and may have bugs or be broken in unexpected ways.  
<p>Warning 2: Sediment mass conservation is -close- to okay, but not 100%. The issue appears when sediment crosses between triangles of different sizes. This is a bear to solve and has bedeviled me for >2 years; improvements are encouraged.
<p>Warning 3: Some of the empirical expressions are potentially based on observations of experiments extremely far from the conditions of interest.  Oh, and these empirical expressions were on Earth, not Mars. They are based on nondimensional scaling, so should at least be reasonable, but there is much that could be going wrong. <br>

<p>Using these operators ... I include a copy of the relevant module file in whatever scenario I'm running.  Presumably you could install them somewhere else on the path that python looks at and they would work as well.

Note that now the model_param.py file is now used to set parameters across modules.

<b>3) Example scenarios:</b><br>
<p>a) Idealized terrain
<p>b) Jezero
<p>c) Coholich flume experiment, LC/LR

<b>4) Viewing the results:</b> 
<p>For fast visualization of results, I use AnugaViewer: https://sourceforge.net/projects/anuga/files/anuga_viewer_windows/<br>
Elevations or other quantities in the domain can also be converted into ArcMap-readable ascii files. SWW output can be viewed in QGIS with the crayfish extension.
