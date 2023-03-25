.. _UsefulTools:
Useful Tools
----------------

Setting your workspace
^^^^^^^^^^^^^^^^^^^^^^^^


A.1- Download and install the PIC code Smilei
^^^^^^^^^^^^^^^^^^^^^^^^

Before you start your practical work,
you must install the PIC code Smile on your computer. 
To do this, go to the Smilei webpage https://smileipic.github.io/Smilei/Use/installation.html, select your operating system 
and carefully follow all the instructions. 
If you have already installed the PIC code Smilei, skip this part.
At the end of this step, you must have the files called ``smilei`` and ``smilei_test`` in the folder you intalled Smilei.

Besides install the PIC code `Smilei`, install postprocessing library `happi` and keep the location of your executable filles `Smilei` and `happi`
in mind. (`e.g.`  ``path/to/executable``)


A.2- Prepare your simulation
^^^^^^^^^^^^^^^^^^^^^^^^

- Enter your working space (folder where you will place your simulations)

.. code-block:: bash

    cd path/to/workdir

where ``path/to/workdir`` is the path to your workdir.
- In your workdir, create a new simulation folder, for example a folder called ``sim1``, where you will run your simulation 1:

.. code-block:: bash

    mkdir sim1
    cd sim1

Each time you do it, choose a convenient name of the folder to
remember which simulation it contains. In order to avoid overwriting data, it is recommended to create 
a new simulation folder for each simulation.

- Inside the simulation folder, create a symbolic link to the Smilei executables. The comand line to create a symbolic link to Smilei is
.. code-block:: bash

    ln -s path/to/executable/smilei
    ln -s path/to/executable/smilei_test

As antecipated ``path/to/executable`` is the path to your executable filles. So you need to insert the actual path where your files ``smilei``
and ``smilei_test`` are.


- Inside the simulation folder, you also need a file to submit your simulation job to the job scheduler, `e.g.`` ``submission_script.sh``.
.. code-block:: bash
   
   #!/bin/bash
   NP=$1  # argument (nb of MPI processes)
   FTR=$2 # argument (file to read = .py file)
   DIAG=$3 # argument (put -d if you want the script to run python diags) 

   PS=/path/to/Smilei_Folder #you must change this path according.
   # input/results/smilei(.exe) directories
   INPDIR=$PS/Input #you must create the folders Input and Results before
   RESDIR=$PS/Results
   SMIDIR=$PS/Smilei

   DIAGDIR=$PS/Smilei/smilei-laser_from_file/scripts/TPUPMC # script for the diagnostics

   # preparing the code
   START=`pwd`
   rm -rf $RESDIR/$FTR
   mkdir -p $RESDIR/$FTR
   cp $INPDIR/$FTR.py $RESDIR/$FTR
   cd $RESDIR/$FTR
   pwd

   # running the code
   DYLD_LIBRARY_PATH=$HDF5_ROOT_DIR/lib mpirun -np $NP $SMIDIR/smilei $FTR.py
   #mpirun -np $NP $SMIDIR/Smilei $FTR.py

   cd $START
   # running the diagnostics
   if [ "$DIAG" == "-d" ]
      then
         $DIAGDIR/SmileiQt.py $RESDIR/$FTR
   elif [ "$DIAG" == "-davg" ]
      then
         $DIAGDIR/SmileiQt_avg.py $RESDIR/$FTR		
   fi

Check the script ``submission_script.sh`` above and modify accoding.





Finally, you need your input file with instructrions for `Smilei`. In this pratical work, we will call it ``InputNamelist.py``. All instructions to run a simulation on SPW is in this fille.

A.3- Run your simulation
^^^^^^^^^^^^^^^^^^^^^^^^
- Check if you have all the required files (executables, submission script, input namelist) through the command:

.. code-block:: bash
   
    ls
in your prompt. In this present case, you must have the symbolic links ``smilei`` and ``smilei_test`` , the submission file ``submission_script.sh`` and the Input namelist ``InputNamelist.py``.


- Now you check your namelist does not contain syntax errors, use the ``smilei_test`` executable on the namelist (you will need to load the same libraries used for the code compilation): ``./smilei_test InputNamelist.py``. If you see the line ``END TEST MODE``, the namelist does not contain syntax errors and can be run.

- Launch your simulation job:

.. code-block:: bash
   
    ./run_new.sh 2 InputNamelist




A.4- Postprocess your simulation results
^^^^^^^^^^^^^^^^^^^^^^^^

- Open ``IPython`` (before, you will need to load the Python modules and define variables like how you did to compile the code, and be sure you have compiled ``happi``):

.. code-block:: bash
   
    ipython

- Import the libraries you need:

.. code-block:: bash
   
    import happi
    import numpy as np
    import matplotlib.pyplot as plt 

The output files have the extension ``.h5`` and can be opened  with the postprocessing library ``happi``. You will need also the 
file ``smilei.py``, generated at the start of your simulation.

- Open your simulation:

.. code-block:: bash
   
    S = happi.Open("path/to/my/results")

again, ``"path/to/my/results"`` is an example, you need to put the path of your simulation. 
If you use simply ``S = happi.Open()``, the library ``happi`` open the results inside the current working directory.

-  Now you can use the commands in the section postprocessing.

A.5- Command line cheatsheet
^^^^^^^^^^^^^^^^^^^^^^^^

- ``pwd``: shows the path of the current working directory.

- ``cd path``: go to ``path``

- ``ls``: shows the content of the current directory.

- ``ls path``: shows the content in ``path``.

- ``rm file``: removes ``file``. To remove a folder, you will need an additional flag: ``rm -r folder`` (be careful).

- ``cp source_file destination_path``: copies ``source_file`` to the ``destination_path``.

- ``scp source_file destination_path`` : same as ``cp``, but you can also transfer folders and files to a different machine, e.g. from the cluster to your computer and vice versa. You have to provide your username, the server address and your password, e.g. ``scp source_file username@server:/destination_path/``. This command can be used to transfer output files from the cluster to your computer for later postprocessing if so you prefer (of course larger data files will need more time to transfer).

- ``mv source destination``: move ``source`` (can be a file or directory) to a ``destination``. If the ``destination`` does not specify a path, the command renames ``source`` with the name
``destination``.

- ``ipython``: opens ``Ipython``, where also the previous commands can be used. To run a Python script inside this interface, use ``%run script_name.py``.

.. _Postprocessing:
Postprocessing
^^^^^^^^^^^^^^^^^^^^^^^^

A fundamental part of working with simulation codes is the 
postprocessing of the results. Smilei includes an entire ``Python`` library 
for postprocessing. 
However, to plot your first results and make quantitative evaluations 
you do not need to be an expert of ``Python``.

For your convenience and quick reference, here we include only the commands 
you will need for this practical. Do not hesitate to copy and paste 
the following commands in ``IPython`` and adapt them to the problem you are solving.

Remember that the results are in normalized units. 
The library ``happi`` also allows to convert to SI units, but this will not be taught in this practical 
(details in the `documentation <https://smileipic.github.io/Smilei/Use/post-processing.html>`_).


B.1- Compilation of happi library
^^^^^^^^^^^^^^^^^^^^^^^^

It is sufficient to use the command ``make happi`` in the code folder 
(after you have loaded the Python modules, see the file ``ClusterEnvironment.pdf``). 
Then, to analyze the results of your simulation, open the ``IPython`` interface 
(just use the command ``ipython`` in the command line terminal).

B.2- Open a simulation
^^^^^^^^^^^^^^^^^^^^^^^^
To import the library ``happi`` in ``IPython`` and open a simulation in the folder, use::

   import happi; S = happi.Open("path/to/simulation")

In this specific example the folder’s path is called for example ``"path/to/simulation"`` 
(use the path of your simulation instead!). 

The last command will create an object called ``S``, our simulation, 
which contains all the necessary data, taken from the input namelist and from the 
output files. 

You can easily access parameters from the input namelist, for example::

   S.namelist.Lx
   S.namelist.dx

In general, if you tap ``S.`` or add the name of the blocks and then use the tab key, 
you will see the available blocks and variables.

B.3- Plot diagnostics
^^^^^^^^^^^^^^^^^^^
To open a specific diagnostic, like the ``Probe1`` defined in the namelist, 
and plot the longitudinal electric field ``Ex`` contained in that diagnostic, use::

   S.Probe.Probe1("Ex").plot()

Other physical fields defined on the grid that you can plot are for example ``Ey``
(the electric field component in the `y` direction), 
``Rho`` (the charge density). Remember that you can also specify operations 
on the fields, like ``Ey*Ey+Ex*Ex``, when you declare your variable.

By default, the last command will only plot the requested field obtained 
in the last simulation output available for that diagnostic. 
You may instead be interested in a specific iteration of the simulation (in code units), 
like iteration `1200`. To plot only that timestep, just specify it inside the diagnostic block::

   S.Probe.Probe1("Ex", timesteps=1200).plot()

Remember that this timestep corresponds to physical time ``1200*dt``, where ``dt`` 
is the simulation timestep, which can be found with ``dt=S.namelist.Main.timestep``.

To know which iterations are available in your diagnostic, you can use::

   S.Probe.Probe1("Ex").getAvailableTimesteps()

B.4- Visualize multiple timesteps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Normally you have a sequence of outputs, so you may want to see an animation 
of the outputs or to be able to slide between the saved timesteps. 
It is possible to do it with these commands respectively::

    S.Probe.Probe1("Ex").animate()
    S.Probe.Probe1("Ex").slide()

In the last case, just slide with the horizontal bar to see the evolution of the plotted quantity at
different iterations.

B.5- Modify elements of the plot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Like in Python, you may be interested into specifying the figure number, 
or change the colormap, or specifying a maximum or minimum value plotted. 
You can include the same corresponding keywords inside the plot/animate/slide command. 
As an example where all these elements are specified::

   S.Probe.Probe1("Ex").plot(figure=2, vmin = -0.1, vmax = 0.1 , cmap = "seismic")

B.6- Plot multiple lines
^^^^^^^^^^^^^^^^^^^^^^^^^
You may be interested in visualizing multiple curves in the same plot window. 
Then the command ``happi.multiPlot`` is what you need.

For example, if you want to plot two quantities from the same simulation, 
scaling them through multiplying factors::

   import happi
   S = happi.Open("path/to/simulation")
   E = S.Probe.Probe1("0.1*Ex", timesteps=1000, label = "E")
   rho = S.Probe.Probe1("-10.*Rho", timesteps=1000, label="charge density")
   happi.multiPlot(E, rho, figure = 1)

The previous example draws two curves, but you can use multiPlot to plot more curves.

Note that you can plot also different timesteps from the same simulation with the same procedure. 
Similarly, you can plot two quantities from two or more simulations::

   import happi
   S1 = happi.Open("path/to/simulation1")
   Ex1 = S1.Probe.Probe0("Ex",timesteps=1000)
   S2 = happi.Open("path/to/simulation2")
   Ex2 = S2.Probe.Probe0("Ex",timesteps=1000)
   happi.multiPlot(Ex1,Ex2)

B.7- Export the data
^^^^^^^^^^^^^^^^^^^^
Those shown above are all the ``happi`` commands you may need for this practical. 
If you prefer instead to analyze your results with ``numpy`` arrays in Python, 
you can easily export your diagnostic to a ``numpy`` array, for example::

   import happi
   import numpy as np
   S = happi.Open("path/to/simulation")
   myArrayVariable = S.Probe.Probe1("Ex").getData()
   myArrayVariable = S.Probe.Probe1("Ex", timesteps=1200).getData()
   myArrayVariable = np.asarray(myArrayVariable)

In case you want to export the data to a text file ``.txt`` and read it with 
another language, you can write this array on a text file using::

   np.savetxt("file_name.txt", myArrayVariable)
