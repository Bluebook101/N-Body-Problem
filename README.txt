This code was written and compiled on a Windows machine and so the terminals and 
commands for Gnuplot were written to plot on a Windows computer with X11 separately 
installed. 

I have attempted to write a code which I think will compile on Mac, but this has not
been tested and so I do not know if it will compile properly. It was written by
looking at how Gnuplot_From_C was written as that ran on a Mac and reverse-engineering 
it. 


The difference involves adding an extra definition in the main code which is then 
called before the plot. The definition of GNUPLOT_SCRIPT is also changed. The 
definition differences are shown below:

----------------------------------------------------------------------------------------
Mac:

#define GNUPLOT.EXE "gnuplot"
#define GNUPLOT_SCRIPT "SEM.script" 

Later on in main()

snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );

----------------------------------------------------------------------------------------
Windows:

#define GNUPLOT_SCRIPT "./SEM.script" 

Later on in main()

snprintf(command, sizeof(command), "%s", GNUPLOT_SCRIPT );

----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------





There is a difference in the .script file too, with the text below added for 
Windows at the beginning of the file:
----------------------------------------------------------------------------------------

#!/bin/gnuplot -persist
set terminal x11

----------------------------------------------------------------------------------------

Hopefully the Mac version will compile and run as it does on Windows, and as is 
demonstrated by the output seen in the report for this assignment.