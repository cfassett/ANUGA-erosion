# ANUGA-erosion
(Experimental) Operators for adding Erosion to ANUGA hydro, and scenarios for erosion in natural dam-breaching floods on Mars.

1) Notes on installing ANUGA hydro:
See https://github.com/GeoscienceAustralia/anuga_core/blob/master/INSTALL.rst<br>

I've recently had consistent sucess in Windows and Linux installing the ANUGA dependencies via conda or Anaconda (https://www.anaconda.com/download/).  I create an ANUGA Anaconda environment, say "anuga" using Anaconda navigator, then at an Anaconda prompt, then "activate anuga" where "anuga here" is whatever the anuga environment is called.  The first time you need to follow the conda instructions in the INSTALL file above.  If all is well, ANUGA will then import in python.  Any future time you want to use ANUGA, it is basically a matter of opening an Anaconda Prompt and activating the ANUGA environment "activate anuga".

2) Using the New Operators:
Warning 1: These operators (and other modules) are experimental and may have bugs or be broken in unexpected ways.  <br>
Warning 2: Sediment mass conservation has been verified to be okay, but otherwise the physics could still be really wrong. In particular, some of the empirical expressions are based on observations of experiments extremely far from the conditions in the big floods being studied here.  Oh, and they were on Earth, not Mars.<br>
<Br>&nbsp
Installing these operators ... I usually just include a copy of the relevant module in whatever scenario I'm running.  Presumably you could install them somewhere else on the path that python looks at and they would work as well.

3) Example scenarios

4) Viewing the results:  
For fast visualization of results, I use AnugaViewer: https://sourceforge.net/projects/anuga/files/anuga_viewer_windows/
The domain can also be converted into ArcMap 
