# ANUGA-erosion
(Experimental) Operators for adding Erosion to ANUGA hydro, and scenarios for erosion in natural dam-breaching floods on Mars.

<b>updated, Oct. 2, 2019</b>

1) Notes on installing ANUGA hydro:<br>
See https://github.com/GeoscienceAustralia/anuga_core/blob/master/INSTALL.rst<br>

<p>I've recently had consistent success on  multiple Windows and Linux systems installing the ANUGA dependencies via the conda/anaconda (https://www.anaconda.com/download/) instructions.  I create an ANUGA conda environment (make sure to choose python=2.7), following the conda instructions in the INSTALL file above for all the packages that ANUGA needs. Then, when I want to run code callign anuga, I use "conda activate anuga".  {It could be source activate anuga, depending on vagaries of your system setup}.  If all is well, the full ANUGA package will  import in python. 

2) Using the New Operators:<br>
<p> Warning 1: These operators (and other modules) are experimental and may have bugs or be broken in unexpected ways.  <br>
Warning 2: Sediment mass conservation is close to okay, although there is still a slight bug that can occur at interfaces where mesh size changes.  The physics, of course, could still be really wrong. In particular, some of the empirical expressions are based on observations of experiments extremely far from the conditions I'm interested in.  Oh, and these empirical expressions were on Earth, not Mars. They are based on nondimensional scaling, so should at least be reasonable, but there is much that could be going wrong. <br>
<Br>
Installing these operators ... I include a copy of the relevant module file in whatever scenario I'm running.  Presumably you could install them somewhere else on the path that python looks at and they would work as well.
<br>
Note that now the model_param.py file is now used to set parameters across modules.

3) Example scenarios:<br>
a) Idealized terrain
b) Jezero

4) Viewing the results:  
<p>For fast visualization of results, I use AnugaViewer: https://sourceforge.net/projects/anuga/files/anuga_viewer_windows/<br>
Elevations or other quantities in the domain can also be converted into ArcMap-readable ascii files. 
