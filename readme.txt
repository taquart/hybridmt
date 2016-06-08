hybridMT
Seismic moment tensor inversion and moment tensor refinement package
version 1.1.1
last update: 2016.06.08

AUTHORS

Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
Patricia Martinez-Garzon <patricia@gfz-potsdam.de>

DOWNLOAD

The newest version of the hybridMT package together with 
installation instructions, quick start guide, documentation and examples 
are available at: http://www.induced.pl/hybridmt and 
http://www.gfz-potdam.de/hybridmt.

INSTALLATION

Unpack the downloaded toolbox to either the directory that exists in your 
MATLAB search path or to the newly created directory. Run MATLAB and 
proceed to the directory containing the toolbox files. Execute the 
hybridmt_install.m script from command window. The script will add/refresh
path to the directory where the toolbox is placed and create the searchable
help help system. 

ACQUIRING HELP

1) The package is provided with HTML help system that integrates 
   with MATLAB environment once hybridmt_install.m script is called. To 
   access the package help system press F1 (this will open the MATLAB  
   documentation HOME page. Proceed downwards through the page and click  
   the "Supplemental Software" at the bottom of the page. This should  
   result in opening of hybridMT package help system. In older versions of 
   MATLAB (<2014a) the help system of hybriMT package may be integrated 
   within MATLAB documentation. In this case search for "hybridMT" package
   in the list of topics of MATLAB documentation or write "hybridMT" in the 
   search box.

2) Short description of hybridMT package is available online through
   http://www.induced.pl/hybridmt and http://www.induced.pl/focimt web 
   pages.

3) The HTML help files are stored in /html subdirectory of your local 
   installation of the toolbox. You can open index.html file to open the 
   help in your web browser. 
   
4) The package contains documentation in PDF format.

WEBSITE

   http://www.induced.pl/hybridmt
   http://www.induced.pl/focimt
   http://www.gfz-potsdam.de/hybridmt

LICENSE

The package is delivered under GPL v3 license. See the license.txt file for
details.

ACKNOWLEDGEMENTS

focMT uses portions of FORTRAN code by Pawel Wiejacz (Institute of 
Geophysics, Polish Academy of Sciences, Warsaw, Poland) related to the 
seismic moment tensor inversion. fociMT binaries uses CAIRO library 
(http://cairographics.org/). The source code contains routines from PSMECA 
program, which is a part of Generic Mapping Tools (GMT) software package 
available to download from http://gmt.soest.hawaii.edu/. The 1D velocity 
model ray-tracing routines were translated and upgraded from FORTRAN from 
hypoDD v1.3 package by Felix Waldhauser, see 
http://www.ldeo.columbia.edu/~felixw/hypoDD.html for details.

