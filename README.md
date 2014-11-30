Biophysics
==========
Code with biophysical models written for the subject [Biophysics](http://materias.df.uba.ar/biofisa2014c2/) at the [Physics Department](http://df.uba.ar/) in the University of Buenos Aires. 

In models.py you will find ode that you will able to integrate using the scipy.ode integrator or the scipy.odeint integrator.
* *HH* stands for the Hodgkin-Huxley model
* *cable* is the pde of the cable equations
* *soma* and dendrite are HH like equations
* *nmda* are the dynamic equations of a nmda channel

The *HH*, *LR* and *llinas* models have been rewritten in python from the orignal MatLab version by Jacobo Sitt that you can find in [here](http://materias.df.uba.ar/biofisa2014c2/guias/)

nature.ipynb is an IPython notebook for reproducing the fgure 1d in [this](http://www.nature.com/nature/journal/v464/n7293/full/nature08947.html) paper.  
