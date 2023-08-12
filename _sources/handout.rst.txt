Handout
---------------------


Before starting the practical work, we will present the general theory behind surface plasma waves. 
Then, we will discuss the normalized system of units commonly used in PIC codes. 
After this brief introduction, we will begin building the PIC code namelist and so, the practical work.

Specifically, we will start the practical work by creating a namelist that defines the simulation
box and simulates a laser pulse propagating in the vacuum. Once we confirm the simulation is running well,
we will add a plasma target to the system consisting of ions and free electrons, which emulate a metallic target. 
We will initially assume the plasma target has a flat target surface. 
Next, we will assuming periodic structure across the target surface and
finally, optimize the target surface to improve the laser-plasma interaction. 
After each of these steps, we will study and analyze relevant diagnostic tools 
to understand the relevant physical phenomena of surface plasma waves.

To complete the practical work, read the instructions in the following sections, 
complete the related simulation exercises and ask the instructor if you have any doubts. 

Prelude
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


0.1- Introduction to Surface Plasma Waves
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The interaction of an intense laser with an overdense plasma with a
sharp density gradient provides ways to accelerate charged particles up
to the relativistic limit. A specific optimal set up that might improve
the effectiveness of the laser-plasma interaction is provided by the use
of periodic grooves (called gratings) at the metal surface. Whenever
appropriately chosen, it can induce the laser to excite the called
Surface Plasma Wave (SPW). Specifically, the periodic groves are
engraved in a solid metal, which, as soon as it starts interacting with
the intense and short laser pulse (:math:`>10^{18}` W/cm\ :math:`^2` and
:math:`<100 f`\ s) at high contrast becomes a sharp-edged over-dense
plasma. If the frequency and wavelength of the interacting laser wave
match the ones described by the dispersion relation of the surface
plasma waves, up to :math:`70\%` of the energy from the laser can be
transferred to the plasma. Moreover, a significant percentage of the
electrons can be accelerated along the surface at relativistic speed.

Exploiting these properties at ultra-high intensity could drive surface
plasma waves with large amplitude, leading to generation of
unprecedented current of energetic electrons as well as bright coherent
radiation. Resonant excitation of SPW implies the knowledge of
the dispersion relation.


0.2- Surface Plasma Waves dispersion relation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As anticipated, Surface Plasma Wave (SPW), also known as surface plasmon
polariton, is a propagating electromagnetic wave (let us denote by
:math:`\omega` its angular frequency) confined at the interface between
a dielectric/transparent medium (with dielectric function
:math:`\epsilon_d(\omega)>0`) and a conducting/opaque medium (with
dielectric function :math:`\epsilon_c(\omega)<0`). Of particular
importance for this study are SPW at the interface between a vacuum
(:math:`\epsilon_d(\omega)=1` located in the region :math:`x<0`) and a
collisionless plasma, with density :math:`n_0` located in the region
:math:`x>0`. In the nonrelativistic and cold plasma limit, we can
describe it using the Drude model, that is taking
:math:`\epsilon_c(\omega) = 1-\omega_p^2/\omega^2`, with
:math:`\omega_p = \sqrt{e^2 n_0/(\epsilon_0 m_e)}` the electron plasma
frequency, :math:`e` the elementary charge, :math:`m_e` the electron
mass and :math:`\epsilon_0` the vacuum permittivity. To the laser be
able to excite the SPW, the plasma needs to be opaque to electromagnetic
waves with angular frequency :math:`\omega`, what implies in
:math:`n_0>n_c` where :math:`n_c= \epsilon_0 m_e \omega^2/e^2` is the
so-called critical density (correspondingly, :math:`\omega_p > \omega`).

The SPW has the following interesting properties: (i) it propagates
along the interface (we denote here by :math:`y` its propagation
direction) with a phase-velocity smaller than :math:`c`, the speed of
light in vacuum; (ii) it is a transverse-magnetic mode (TM-mode) with
the magnetic field along the direction :math:`z` perpendicular to its
direction of propagation, (iii) it supports both transverse
(:math:`E_x`) and longitudinal (:math:`E_y`) electric fields, and
(iv) all fields decay exponentially with the distance
:math:`\vert x\vert` perpendicular to the plasma-dielectric interface.

The dispersion relation for such a wave reads:

.. math:: \frac{c^2k^2}{\omega^2} = \frac{\epsilon_d(\omega)\,\epsilon_c(\omega)}{\epsilon_d(\omega)+\epsilon_c(\omega)} = \frac{\omega_p^2/\omega^2-1}{\omega_p^2/\omega^2-2},

and the wave only exists for
:math:`\vert \epsilon_c(\omega)\vert > \epsilon_d(\omega)`, which in the
case of a simple Drude model requires :math:`\omega_p^2/\omega^2>2`
(correspondingly :math:`n_0 > 2 n_c`). Even though SPW can exists in the
presence of a flat surface, the easiest way to excite them in
laser-plasma experiments is by analogy with solid state case by
impinging a laser pulse onto, for example, a plasma grating, *i.e.* a
plasma which surface is modulated so that the vacuum-plasma interface is
located at :math:`x_g(y) = (h/2)\,\sin(2\pi y/d)`, with :math:`h` the
grating depth and :math:`d` its period. This *grating-coupling* approach
allows to excite a surface plasma wave at the same frequency
(:math:`\omega\equiv\omega_0`) than the driving laser pulse when the
laser angle of incidence :math:`\theta_{\rm inc}` satisfies the
phase-matching condition
:math:`\sin(\theta_{\rm inc}) = \sqrt{(\omega_p^2/\omega_0^2-1)/(\omega_p^2/\omega_0^2-2)} - l \lambda_0/d`
with :math:`\lambda_0=2\pi c/\omega_0` the laser wavelength and
:math:`l \in \mathbb{Z}`. Throughout this practical work, we will focus
our attention onto the resonance condition for :math:`l=1`.

When dealing with high-intensity lasers interaction, and in particular
when the laser electric field :math:`E_0` becomes of the order of
:math:`m_e c \omega_0/e` (*i.e.* for a normalized vector potential
:math:`a_0 \equiv e E_0/(m_e c\,\omega_0) \gtrsim 1`), the electron
dynamics becomes relativistic. It tells us that, 
if we want to generate super energetic electrons (MeV), we have to consider 
:math:`a_0 \gt 1`.

The phase-matching condition introduced above leads to define the
optimal angle of incidence for the laser to excite the surface plasma
wave:

.. math:: \theta_{\rm opt} = \arcsin\left( \sqrt{\frac{n_0/n_c-1}{n_0/n_c-2}} - \frac{\lambda_0}{d}\right)\,.

Thus, we have to irradiate the laser pulse on the plasma target surface 
so that it has an angle :math:`\theta_{\rm opt}` in relation to the normal surface.
In what follows, we will put this simple model to the test by studying
the excitation of SPW in the relativistic regime using PIC simulations.

**Exercise 1:** Use the definition for the critical
plasma density :math:`n_c=\varepsilon_0 m_e \omega_r^2/e^2`, the electron plasma
frequency :math:`\omega_p = \sqrt{e^2 n_0/(\epsilon_0 m_e)}` and assume that :math:`\omega=\omega_0` to show the relation :math:`\omega_p^2/\omega_0^2=n_0/n_c` used in 
the equation above.

**Exercise 2:** Assuming that the diffraction grating has a period equal
to twice the wavelength of the laser, :math:`d=2\lambda_0` what is the resonance angle (in degrees)?


0.3- Physical units normalization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the PIC code Smilei (and also in most PIC codes) all input and output quantities are in normalized units (also called code units) since such a 
system of units facilitates computational calculations. It is crucial to learn how to work 
with this kind of unit since it is widely used in computer simulations.
You can perform conversions from normalized units to SI quantities, following :ref:`Table. 1 <table1>`:


.. _table1:
.. list-table:: Normalizations
   :widths: 25 25
   :header-rows: 1

   * - Physical quantity 
     - Reference for normalization
   * - Time
     - :math:`\omega_r^{-1}`   
   * - Velocity
     - :math:`c`
   * - Charge
     - :math:`e`
   * - Mass
     - :math:`m_e`
   * - Momentum
     - :math:`m_e c`
   * - Energy, Temperature
     - :math:`m_e c^2`
   * - Length
     - :math:`c/\omega_r`
   * - Number density
     - :math:`n_r=\varepsilon_0 m_e \omega_r^2/e^2`
   * - Electric field
     - :math:`m_e \omega_r c /e`
   * - Magnetic field
     - :math:`m_e \omega_r /e`


As shown in :ref:`Table. 1 <table1>`, charge and mass are normalized to the elementary 
charge :math:`e` and electron mass :math:`m_e`, respectively. Additionally, all velocities 
are normalized to the speed of light in vacuum :math:`c`, which naturally arises
from Maxwell’s equations. In this set of units, the normalized
speed of light in vacuum is :math:`1`. 

Now, the unit of time - here defined as :math:`\omega_r^{-1}`, 
with :math:`\omega_r` being the reference angular frequency- it is not defined a priori and is instead
chosen by the user. Once the unit of time (or equivalently the unit of length) 
is chosen, all other units are uniquely defined as detailed in :ref:`Table. 1 <table1>`.

It is worth noting that the density associated with each plasma species is not 
in units of :math:`(c/\omega_r)^{-3}` but instead in units of :math:`n_r=\varepsilon_0 m_e \omega_r^2/e^2`. In the following 
exercises, where you will be required to make conversions from code units, 
assume that the laser wavelength is :math:`\lambda_0=0.8\mu m` (that of a Ti:Sa laser system). 
This choice implies that a unit of length is :math:`c/\omega_r=\lambda_0/2\pi` and a unit of time is
:math:`\omega_r^{-1}=\lambda_0/2\pi c`, what is in accordance with :ref:`Table. 1 <table1>`.


It is  also possible to perform these conversions automatically using the postprocessing library ``happi``, 
but this is not covered in this practical work. Instead, `see tutorials, basic units <https://smileipic.github.io/tutorials/basics_units.html>`_ for more information.


**- Exercise 3:** 
Assuming :math:`\lambda_0=0.8\mu m` (a typical Ti:Sa laser system), 
*i)* what is the value of the critical density :math:`n_c=`? *ii)* What is the value of the normalizing 
electric field :math:`E_0=m_e \omega_r c /e`? 
This choice of :math:`\lambda_0` shall be utilized throughout *all* subsequent exercises. 


0.4- Simulation virtual workspace 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before run the simulations, it is necessary to prepare the interface and the computing environment. 
To do so, follow the guidelines provided in the `UsefulTools`. This includes install Smilei PIC code, 
compiling the code, creating simulation folders. 
It is advisable to create a folder for each simulation to avoid data loss. 
Within each folder, copy the ``InputNamelist.py`` file you with write and create the symbolic links to Smilei executable files (``smilei`` and ``smilei_test``). 
You can check for errors in your input file by using the command ``./smilei_test InputNamelist.py``. 
If everything goes well, the word "END" should appear on your screen (it might be interesting to test after changes in the namelist). 
Then, launch your simulation using the ``submission_script.sh``. It would take some minutes in your computer. After run your first simulation,
compile the ``happi`` library for post-processing and data analysis, if you have not done.

0.5- Simulation general setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To put the simple model presented in previous sections to the test, we
will performed a series of 2D3V PIC simulations with the open-source
code SMILEI, considering the simplest interaction configuration, where a
laser pulse impinges at an incidence angle :math:`\theta_{\rm inc}` onto
a plasma grating. In this way, periodic conditions are used in the
transverse :math:`y`-direction (for both fields and particles) while
absorbing/injecting (so-called Silver-Muller) boundary conditions are
used in the :math:`x` direction for the fields. Particles are either
reflected at the :math:`x_{\rm min}` border of the simulation domain
(even though we have ensured that no particle reaches this border) or
thermalized at the :math:`x_{\rm max}` border.



The driven laser pulse will be represented by a P-polarized Gaussian pulse
with angular frequency :math:`\omega_0`, correspondingly wavelength
:math:`\lambda_0 = 2\pi c/\omega_0=0.8\mu`\ m, waist equal to
:math:`6\lambda_0` (:math:`=4.8\mu`\ m) and pulse duration equal to
:math:`\tau_L=10\lambda_0/c` (:math:`\simeq27f`\ s) FHWM. The laser pulse impinges over the plasma interface, through an angle
:math:`\theta_ {inc}` in relation to the normal surface. As discussed,
depending on the laser and plasma characteristics, as well as the
grooves etched on the target interface, both the SPW may be excited, and
electrons may be accelerated.

This laser impinges onto a plasma target of constant density
:math:`n_0` with three different interfaces plasma-vacuum interfaces (the three cases that we are going to analyze in this TP) 
(i) The first simply consists of a continuous interface located at :math:`x=0`. 
This will be our reference case, to compare with the results observed when grattings are present. (ii)
The second interface to be analyzed consists of a sinusoidal modulation present on the entire surface.
The interface is located at :math:`x_g(y) = (h/2)\,\sin(2\pi y/d)` where :math:`d=2\lambda_0` (:math:`=1.6\mu`\ m)
is the grating period, and :math:`h=0.3 \lambda_0` (:math:`=0.08\mu`\ m) is the
grating depth. This later value is small enough to keep valid the plane
interface condition and large enough to diffract the impinging laser wave, enabling so the surface
plasma wave excitation. (iii) As a third case to be analyzed, we will remove the modulation where there is no laser-plasma interaction.
This will reduce any radiation losses and consequently increase the efficiency of the particle-SPW coupling.
For this study, the plasma grating consists of electrons with a very small initial temperature (cold) of
:math:`T_e -> 0` as well as a neutralizing background of
immobile ions.

In all simulations, we will assume the box extends over :math:`15\lambda_0`
(:math:`=12\mu`\ m) in the :math:`x`-direction (roughly
:math:`9.6\lambda_0` (:math:`=12.8\mu`\ m) of vacuum and
:math:`3\lambda_0` (:math:`=2.4\mu`\ m) of plasma), and
:math:`72\lambda_0` (:math:`=57.6\mu`\ m) in the :math:`y`-direction.
The spatial resolution was set to
:math:`\Delta x = \Delta y = \lambda_0/64` (:math:`=0.00625\mu`\ m).
Analogously to previous case, the simulation timestep is chosen to be
:math:`\Delta t = 0.95~\Delta x/\sqrt{2}` that corresponds to the Courant-Friedrich-Lewy (CFL) condition. Every cells contains
initially :math:`16` randomly distributed particles of each species.

0.6- Simulation input namelist
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To reproduce the physical experiment just described, we need
write a script (here called input namelist file) to be interpreted by the PIC code Smilei.
The complet ``InputNamelist.py`` script with all the instructions can be found at the end of this Handout, 
however, we are going to build it, step-by-step.

To write the ``InputNamelist.py`` file, you can use any text editor, *e.g.,* ``nano``, ``gedit`` or ``vim``, still, our recomendation is use ``nano``. 
Ok, let us start it, but before...


.. The file starts with the definition of physical constants and units,
.. mesh points, integration timestep, etc. Some of these parameters are inserted in the first block of the simulation, called ``Main`` block. 
.. Others may be useful for conversions between units or to define other variables in the file.

.. In the ``Main`` block, you also find the geometry of the simulation, which is ``Cartesian2D``.
.. The grid is defined on cartesian space with coordinates `x`, `y`, but the particles velocity is defined in the 3D space 
.. `x`, `y`, `z` (See :ref:`Fig. 1 <Mesh_and_Reference>`). 

.. Note, all the lines starting by the hashtag simbol are comments to the code

.. At the end of the InputNamelist.py file,
.. there are blocks starting with the word ``Diag``. As you can imagine, these blocks are for the diagnostic of the code. 
.. The first Diag is a ``DiagProbe`` defined on a line 
.. (so a 1D diagnostic), on the propagation axis of the laser (the `x` axis). 
.. This diagnostic returns the value of some physical fields 
.. along that axis. We call this probe ``Probe0`` (the ``0`` because 
.. it is the first ``Probe`` in the namelist). The second diagnostic block is 
.. a ``DiagProbe`` defined on the plane `xy` (so a 2D diagnostic). 
.. This is the second probe of the namelist, so it is called ``Probe1`` (Python starts counting from zero.)
 
.. Inside these ``DiagProbe`` blocks, the physical quantities 
.. we are interested in are specified. For this practical, you do not have 
.. to write new postprocessing scripts, since Smilei includes a postprocessing 
.. library, called ``happi``, to analyze the code results. For the purposes of 
.. the practical, the happi commands explained in Postprocessing section of the :ref:`Useful Tools <Postprocessing>` 
.. and in scripts are enough to start.

.. _Exercise4:
**- Exercise 4:** The script we want to build considering all proposals above, *i)* What are the trasnversal size ``Lx`` and the longitudinal size ``Ly`` (Longitudinal with reference to the SPW wave we want to excite) 
of the simulation box? *ii)* How many cells there are in the simulation box? 





.. See Figures :ref:`1 <Mesh_and_Reference>` for reference, and find these lengths in the InputNamelist.py.





.. .. _laserpulseinvacuum:
.. First part: Laser pulse
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^

.. 1.1- Laser pulse propagating in vacuum
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^

.. Make sure that in  `InputNamelist.py <https://github.com/>`_, only the commands defining the simulation box
.. and the laser pulse are uncommented. If everything is correct and you have 
.. followed all the steps up to this point, launch your simulation. In this simulation,
.. we will only see the laser pulse propagating in vacuum.


.. 1.2- Laser pulse propagation diagnostics
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^


.. Second part: Laser-plasma interaction
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^


.. 2.1- Laser pulse interacting with a flat target
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^
.. Nessa primeira parte do nosso trabalho pratico, nós vamos considerar um alvo de plasma 
.. cuja superficie é plana. Nesse caso, haverá uma interação entre o pulso laser e o plasma que 
.. será suficiente para acelerar alguns eletrons. Para realizar essa simulação, descomente os blocos indicados 
.. como seção 2.1.


.. 2.2- Flat target disgnostics
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^


.. 2.3- Laser pulse interacting with a structured target
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^
.. Nós podemos aprimorar o acoplamento laser-plasma para a excitação das ondas de superficie de plasma 
.. assumindo gratings. Assim, nessa parte do nosso trabalho pratico, vamos modificar a função densidade do plasma 
.. para incluir gratings. 


.. 2.4- Grating target disgnostics
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^

.. 2.5- Comparation
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^


.. 2.6- Laser pulse interacting with a optimized target
.. ^^^^^^^^^^^^^^^^^^^^^^
.. ^^^^^^^^^^^^^^^^^^^^^^



Laser Pulse
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


1.1- Laser pulse propagating in vacuum
^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^


As antecipated, as a starting point in this practical work we will simulating a laser propagating through vacuum. 
To accomplish this, we will begin by defining a simulation box and including the necessary requirements for a Smilei code.
We will be using the script provided below as a starting point for our simulation. 
It is important to carefully review the code and its implementation of box conditions, as it will serve as the foundation for next simulation 
and your future works.

This is the InputNamelist.py::

  import numpy as np

  ### NORMALIZATIONS
  t0 = 2.0*np.pi      # optical cycle [in code units]
  l0 = 2.0*np.pi      # laser wavelength [in code units]

  ### SIMULATION BOX
  Lx   = 15.*l0       #simulation box size Lx
  Ly   = 72.*l0       #simulation box size Ly

  #### RESOLUTION
  res  = 64           #number of divisions each wavelenght

  dx   = dy = l0/res  #dx,dy cell size
  dt   = 0.95 * dx/np.sqrt(2.) #time increment


  Tsim = 30.*t0       #number the cycles the simulation will run

  ###############################################################################
  ####################### MAIN CODE #############################################

  Main(
      geometry = "2Dcartesian",
      interpolation_order = 2 ,

      cell_length = [dx,dy],
      grid_length = [Lx,Ly],
      number_of_patches = [1,16],  #patches in x,y-directions. The product x*y must be larger than the number of processor you are using

      timestep = dt,
      simulation_time = Tsim,

      EM_boundary_conditions = [['silver-muller'],['periodic']],

      random_seed = smilei_mpi_rank,
      print_every = int(Tsim/dt/10.)
  )


  ################## LASER PULSE ######################################

  ### LASER PROPERTIES
  aL=0.1                          #normalized laser amplitude a0=0.86 sqrt(I/10^18)
  theta_deg = 30.                 #angle of incidence in degrees ArcSin[Sqrt[(n/gammal - 1)/(n/gammal - 2)] - 1 l/d]*180/Pi?
  theta  = theta_deg*np.pi/180.   #angle of incidence in rads

  waist  =  6.*l0     # laser waist
  Tfwhm  = 10.*t0     # laser fwhm

  xfoc = Lx/2      # the laser is focalized over the surface
  yfoc = Ly/2.     # the laser is focalized at the target's center

  ### LASER PULSE
  LaserGaussian2D(
      a0              = aL,
      omega           = 1.,
      focus           = [xfoc,yfoc],
      waist           = waist,
      incidence_angle = theta,
      time_envelope   = tgaussian(start=0., duration=4.*Tfwhm, fwhm=Tfwhm, center=2.*Tfwhm, order=2)
  )

  #################### DIAGNOSTICS ####################################
  sub_grid=2 #to reduce the number of points in the output file
  Period = int(t0/dt)

  DiagFields(
      every = Period,
      fields = ['Bz','Ey','Ex' ],
      subgrid = np.s_[::sub_grid, ::sub_grid]
  )


Now let us plot and analyse the result from our simulation.
For this we can use the commands below. 
Specifically, you can either create a python script and include all the commands simultaneously, 
either you can open the python interface and type each of the lines (recommended for simple figures).
Inside the result folder (`i.e.`` sim1), use the script bellow to observe the `Bz` field ::

  import happi
  S=happi.Open('.',verbose=False)
  S.Field(0,"Bz").slide(cmap='bwr', aspect=1)



Using comand `.plot()` you can plot the field in a given instant of time, for example, to check the instant
of time :math:`t=10\tau_0`, we can use the following script ::

  import happi
  import matplotlib.pyplot as plt

  S=happi.Open('.',verbose=False)
  timesteps=S.Field(0,"Bz").getTimesteps()
  tcycle=10 #this corresponds to the time the laser oscilated 10 times.
  S.Field(0,"Bz", timesteps=timesteps[tcycle]).plot(cmap='bwr', aspect=1, saveAs="Bz_t=10.png")
  plt.show()

.. _Exercise5:
**- Exercise 5:** What are the trasnversal size ``Lx`` and the longitudinal size ``Ly`` obsereved in the Figure? Is it consistent with our script? Explain.

.. _Exercise6:
**- Exercise 6:** Modify the script above to plot the electric field `Ey` at the instant `Tsim`, used in the Smile script.
What is the maximum amplitude of the `Ey` electric field? In which units?


Laser-plasma interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


2.1- Laser pulse interacting with a flat target
^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^

Now we have seen that our simulation is working correctly, it means the
laser propagate correctly inside the simulation box, we are then ready to add the plasma target.
As our first approach, we are going to assume that the plasma interface is flat.
In the Smilei code, we can define the plasma function using a python function, as follows:

.. code-block:: bash

  def density_profile(x,y):
      if(x>Xsurface):
          return n0
      else:
          return 0.

Equally, we have to define the species that compose the
plasma. In this practical work, we will assume the plasma is an ionized gas composed of free electrons and ions,
with temperature going to zero. We define the species that make up the plasma as follows:

.. code-block:: bash

  Species(
      name = "ion",
      position_initialization = "random",
      momentum_initialization = "cold",
      particles_per_cell = nppc,
      mass = 1836.,
      charge = 1.,
      number_density = density_profile,
      boundary_conditions = [["remove","thermalize"],["periodic"]],
      thermal_boundary_temperature = [T0],
      time_frozen = 0.0
    )

  Species(
      name = "eon",
      position_initialization = "ion",
      momentum_initialization = "mj",
      particles_per_cell = nppc,
      mass = 1.,
      charge = -1.,
      number_density = density_profile,
      temperature = [T0],
      boundary_conditions = [["remove","thermalize"],["periodic"]],
      thermal_boundary_temperature = [T0],
      time_frozen = 0.0
    )

Before running this script, we are also going to change to over the plasma surface, the laser focal position, `xfoc` and add a diagnostic
to measure the electrons phase space.

.. code-block:: bash

  DiagParticleBinning(
      deposited_quantity = "weight",
      every = [0, Period],
      time_average = 1,
      species = ["eon"],
      axes = [
          ["y", 0, Ly, 400],
          ["py",-0.3,0.3, 120]
      ]
  )

Following the steps described, we arrive at the following script to describe the interaction of the laser with a flat surface plasma:

.. code-block:: bash

  import numpy as np

  ### NORMALIZATIONS
  t0 = 2.0*np.pi      # optical cycle [in code units]
  l0 = 2.0*np.pi      # laser wavelength [in code units]

  ### SIMULATION BOX
  Lx   = 15.*l0       #simulation box size Lx
  Ly   = 72.*l0       #simulation box size Ly

  #### RESOLUTION
  res  = 64           #number of divisions each wavelenght

  dx   = dy = l0/res  #dx,dy cell size
  dt   = 0.95 * dx/np.sqrt(2.) #time increment


  Tsim = 72.*t0       #number the cycles the simulation will run

  ###############################################################################
  ####################### MAIN CODE #############################################

  Main(
      geometry = "2Dcartesian",
      interpolation_order = 2 ,

      cell_length = [dx,dy],
      grid_length = [Lx,Ly],
      number_of_patches = [1,16],  #patches in x,y-directions. The product x*y must be larger than the number of processor you are using

      timestep = dt,
      simulation_time = Tsim,

      EM_boundary_conditions = [['silver-muller'],['periodic']],

      random_seed = smilei_mpi_rank,
      print_every = int(Tsim/dt/10.)
  )

  ################## PLASMA ##########

  one_eV = 1/511.e3   # eV definition
  n0   = 20           # plasma density
  T0 = 0.1*one_eV     # plasma temperature- cold
  nppc=16             # number of particles per cell

  Xsurface=Lx-3*l0   #position of the plasma surface, 3*l0 is the plasma size

  def density_profile(x,y):
      if(x>Xsurface):
          return n0
      else:
          return 0.

  Species(
      name = "ion",
      position_initialization = "random",
      momentum_initialization = "cold",
      particles_per_cell = nppc,
      mass = 1836.,
      charge = 1.,
      number_density = density_profile,
      boundary_conditions = [["remove","thermalize"],["periodic"]],
      thermal_boundary_temperature = [T0],
      time_frozen = 0.0
    )

  Species(
      name = "eon",
      position_initialization = "ion",
      momentum_initialization = "mj",
      particles_per_cell = nppc,
      mass = 1.,
      charge = -1.,
      number_density = density_profile,
      temperature = [T0],
      boundary_conditions = [["remove","thermalize"],["periodic"]],
      thermal_boundary_temperature = [T0],
      time_frozen = 0.0
    )


  ################## LASER PULSE ######################################

  ### LASER PROPERTIES
  aL=0.1                          #normalized laser amplitude a0=0.86 sqrt(I/10^18)
  theta_deg = 30.                 #angle of incidence in degrees ArcSin[Sqrt[(n/gammal - 1)/(n/gammal - 2)] - 1 l/d]*180/Pi?
  theta  = theta_deg*np.pi/180.   #angle of incidence in rads

  waist  =  6.*l0     # laser waist
  Tfwhm  = 10.*t0     # laser fwhm

  xfoc = Xsurface    # the laser is focalized over the surface
  yfoc = Ly/2.       # the laser is focalized at the target's center

  ### LASER PULSE
  LaserGaussian2D(
      a0              = aL,
      omega           = 1.,
      focus           = [xfoc,yfoc],
      waist           = waist,
      incidence_angle = theta,
      time_envelope   = tgaussian(start=0., duration=4.*Tfwhm, fwhm=Tfwhm, center=2.*Tfwhm, order=2)
  )

  #################### DIAGNOSTICS ####################################
  sub_grid=2 #to reduce the number of points in the output file
  Period = int(t0/dt)

  DiagFields(
      every = Period,
      fields = ['Bz','Ey','Ex', 'Rho_eon' ],
      subgrid = np.s_[::sub_grid, ::sub_grid]
  )

  DiagParticleBinning(
      deposited_quantity = "weight",
      every = [0, Period],
      time_average = 1,
      species = ["eon"],
      axes = [
             ["y", 0, Ly, 400],
             ["py",-0.3,0.3, 120]]
  )

To check the plasma surface, we can analyze the quantity `Rho_eon`, which refers to the density of electrons-`eon` ::

  import happi
  import matplotlib.pyplot as plt
  S=happi.Open('.',verbose=False)
  timesteps=S.Field(0,"-Rho_eon").getTimesteps()
  tcycle=1
  S.Field(0,"-Rho_eon", timesteps=timesteps[tcycle]).plot(cmap='binary', aspect=1, vmin=0,vmax=20)
  plt.show()

Note that we use the minus sign in front of the quantity `Rho_eon`, since we are interested in the modulus of the density only.
Likewise, we can multiply this quantity by any scalar and even by other quantities, for example, we can measure the excess charge,
ploting the quantity `Rho_ion-Rho_eon`.
Furthermore, we can also use the commands below to see the laser propagation in this new setup ::

  import happi
  S=happi.Open('.',verbose=False)
  S.Field(0,"Bz").slide(cmap='bwr', aspect=1)

Using a more elaborate script and working with the matplotlib package, we can simultaneously plot the electric field and the plasma
For that we can use the following script 

.. code-block:: bash

  import numpy as np
  from matplotlib import *
  import matplotlib.pyplot as plt
  from matplotlib.pyplot import *
  #import sys
  #sys.path.insert(0, ".")
  import happi


  ####### load simulation
  S = happi.Open('.', verbose=False)

  for tver in [20]: #checking the instant of time, tver=20
      plt.clf() 
      fs=11
      sizex=4/2.5
      sizey=15.2/2.5
      plt.figure(1, figsize=(sizex,sizey))
      plt.rc('text',usetex=True)
      plt.rc('font', **{'family': 'serif'  })
      plt.rc('xtick',labelsize=fs)
      plt.rc('ytick',labelsize=fs)
      plt.rc('axes',linewidth=0.4)
      plt.subplots_adjust(bottom=0.13,left=0.12,right=0.8,top=0.98)
      aspect0=0.5
      ax=gca()

      #-------------------------------------------------------------------------------
      lambda0   = 1  #use 0.8 if you want the plot in micrometers.

      diag  = S.Field(0,"Ey")
      timesteps=S.Field(0,"Ey").getTimesteps()
      #print(timesteps)

      diag  = S.Field(0,"Ey", timesteps=timesteps[tver])
      xx  = diag.getAxis('x')
      yy  = diag.getAxis('y') 
      Ey = np.array(diag.getData()[0] )

      diagr  = S.Field(0,'-Rho_eon', timesteps=timesteps[tver])
      xxr  = diagr.getAxis('x')
      yyr  = diagr.getAxis('y')
      Rho  = np.array(diagr.getData()[0] )


      xx=lambda0*xx/(2*np.pi)
      yy=lambda0*yy/(2*np.pi)

      ymin = np.min(yy)
      ymax = np.max(yy)
      xmin = np.min(xx)
      xmax = np.max(xx)

      plt.xlabel('$x/\lambda_0$')
      plt.ylabel('$y/\lambda_0$')

      #to define the range and the ticks
      extent = xmin, xmax,ymin,ymax

      minorLocatorx   =  matplotlib.ticker.IndexLocator(base=3, offset=0)
      ax.xaxis.set_minor_locator(minorLocatorx)
      
      minorLocatory   =  matplotlib.ticker.IndexLocator(base=9, offset=0)
      ax.yaxis.set_minor_locator(minorLocatory)

      ax.xaxis.set_ticks_position('both')

      for line in ax.xaxis.get_ticklines():
          line.set_color('k')
          line.set_markersize(2.5)
          line.set_markeredgewidth(0.3)

      for line in ax.xaxis.get_ticklines():
          line.set_color('k')
          line.set_markersize(5)
          line.set_markeredgewidth(0.6)

      for line in ax.yaxis.get_ticklines():
          line.set_color('k')
          line.set_markersize(2.5)
          line.set_markeredgewidth(0.3)

      for line in ax.yaxis.get_ticklines():
          line.set_color('k')
          line.set_markersize(5)
          line.set_markeredgewidth(0.6)

      ay=gca()
      ay.yaxis.labelpad=-10

      ay=gca()
      ay.xaxis.labelpad=-10


      ax.set_ylim( ymin, ymax )
      ax.set_xlim( xmin, xmax )
      ax.set_yticks([ymin,18,54,ymax])
      ax.set_xticks([xmin,12,xmax])


      #-----------------------------------------------------------------------------
      #mask to let transparent the vacuum region
      rcut=0.1
      Rho[abs(Rho) <  rcut] = np.nan
      #-----------------------------------------------------------------------------
      #parameter
      vmax0=S.namelist.n0
      v0=S.namelist.aL
      #-----------------------------------------------------------------------------

      #plot field
      im1 = plt.imshow(np.real(Ey).T, cmap='seismic', aspect=aspect0, vmin=-v0, vmax=v0, alpha=1, origin='lower', interpolation='bilinear', extent=extent)
      #cb1 = plt.colorbar(im1, shrink=0.5, label='$E_y/E_0$', pad=0.1)

      #plot density
      im2 = plt.imshow(np.real(Rho).T, cmap='binary', aspect=aspect0, vmin=0.00000, vmax=vmax0, alpha=1, origin='lower', interpolation='bilinear', extent=extent)
      #cb2 = plt.colorbar(im2, shrink=0.5, label='$n_e/n_c$', pad=0.09)
  

      plt.savefig('field_%3i.png'%(tver), bbox_inches='tight', dpi=100)



.. _Exercise7:
**- Exercise 7:** Modify the script above to plot the Poynting vector `Sy` at the instant `Tsim`, used in the Smilei script.
What is the maximum amplitude of the  Poynting vector? In which units is it presented? What are the values in the S.I units?

.. _Exercise8:
**- Exercise 8:**
Modify the loop on `tver` to build a plot at each instant in time.
Use the outcome figures to create a video. There are several ways to create the video from pictures. One way is to use the command ::
  
  convert -delay 1 -loop 0 *.png video.gif

In an alternative approach, you could use the `happi` library to create the video ::

  S.Field(0,'Bz').animate(movie="Bz.gif", dpi=120,fps=20,vmin=-0.1,vmax=0.1,ymax=300)


2.2- Laser pulse interacting with a sinusoidal target
^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^

Now we are going to modify our script so that the surface of the plasma will be described by a sinusoidal function. One way to do this consists of modifying the density_profile(x,y) function such as ::

  d=2*l0      # period
  depth=0.1*l0  # depth
  kg= l0/d    # normalized wavevector

  def density_profile(x,y):
      if(x>Xsurface+0.5*depth*np.sin(kg*y) ):
          return n0
      else:
          return 0.

You can check the plasma surface using the following command from `happi` library ::

  import happi
  import matplotlib.pyplot as plt

  S=happi.Open('.',verbose=False)
  timesteps=S.Field(0,"Ey").getTimesteps()
  tcycle=1 
  S.Field(0,"-Rho_eon", timesteps=timesteps[tcycle]).plot(cmap='Blues',xmin=60,xmax=90,ymin=20,ymax=60)
  plt.show()

Likewise, you can change the `x` and the `y` limits using the script used in the previous section, such as ::

    ax.set_ylim( 30, 38 )
    ax.set_xlim( 10, 14 )
    ax.set_yticks([30,32,36,38])
    ax.set_xticks([10,11,13,14])


.. _Exercise 9:
**- Exercise 9:** Create a figure to verify the period of the diffraction grating over the plasma. It is correct according to our demand?


We can then analyze the phase space of the particles. For this, we will use the library `happi` as follows::

  import sys
  sys.path.insert(0, "/Users/samuelmarini/Samuel/Smilei/smilei")
  import happi
  import numpy as np
  import matplotlib.pyplot as plt

  S = happi.Open('/Users/samuelmarini/Samuel/Smilei/Results/spw1', verbose=False)
  tver=70
  timesteps=S.ParticleBinning(0).getTimesteps()
  diag  =  S.ParticleBinning(0, data_log=True,timesteps=timesteps[tver])
  diag.plot()

.. _Exercise 10:
**- Exercise 10:** 
Analyze the longitudinal phase space at different instants of time. What can you interpret with it?
To go further, suppose that ``aL=5``, change the script and describe what happens with the electrons now.
Note that you have to change the diagnostics ranges, in order to make a correct analysis.


2.3- Laser pulse interacting with a optimized target
^^^^^^^^^^^^^^^^^^^^^^
^^^^^^^^^^^^^^^^^^^^^^

We finally reached the last stage of this TP and now we are going to simulate the interaction of the laser pulse with 
the plasma target whose interface is optimized. Find bellow the full InputNamelist.py ::

  import numpy as np

  ### NORMALIZATIONS
  t0 = 2.0*np.pi      # optical cycle [in code units]
  l0 = 2.0*np.pi      # laser wavelength [in code units]

  ### SIMULATION BOX
  Lx   = 15.*l0       #simulation box size Lx
  Ly   = 72.*l0       #simulation box size Ly

  #### RESOLUTION
  res  = 64           #number of divisions each wavelenght

  dx   = dy = l0/res  #dx,dy cell size
  dt   = 0.95 * dx/np.sqrt(2.) #time increment


  Tsim = 72.*t0       #number the cycles the simulation will run

  ###############################################################################
  ####################### MAIN CODE #############################################

  Main(
      geometry = "2Dcartesian",
      interpolation_order = 2 ,

      cell_length = [dx,dy],
      grid_length = [Lx,Ly],
      number_of_patches = [1,16],  #patches in x,y-directions. The product x*y must be larger than the number of processor you are using

      timestep = dt,
      simulation_time = Tsim,

      EM_boundary_conditions = [['silver-muller'],['periodic']],

      random_seed = smilei_mpi_rank,
      print_every = int(Tsim/dt/10.)
  )

  ################## PLASMA ##########

  one_eV = 1/511.e3   # eV definition
  n0   = 20           # plasma density
  T0 = 0.1*one_eV     # plasma temperature- cold
  nppc=16              # number of particles per cell

  Xsurface=Lx-3*l0   #position of the plasma surface, 3*l0 is the plasma size

  waist0=6*l0
  d=2*l0      # period
  depth=1*l0  # depth
  kg= l0/d    # normalized wavevector

  lysize = round(Ly/(2*d))
  nri    = round(waist0/(l0))+1 #three times more ripples than the area illuminated by the laser waist
  phase  = 0 #to change the ripples phase-- check in the begining of the simulation if it is correct.

  def density_profile(x,y):
      for nl in range(0,int(Ly/d)+1):
          if( nl>=lysize-nri and nl<=lysize+nri and y>=lysize*d-nri*l0 and y<=lysize*d+nri*l0 and x>Xsurface+0.5*depth*np.sin(kg*y)):
              return n0
          elif(nl>lysize+nri and y>lysize*d+nri*l0 and x>Xsurface):
              return n0
          elif(nl<lysize-nri and y<lysize*d-nri*l0 and x>Xsurface):
              return n0
      else:
          return 0.


  Species(
      name = "ion",
      position_initialization = "random",
      momentum_initialization = "cold",
      particles_per_cell = nppc,
      mass = 1836.,
      charge = 1.,
      number_density = density_profile,
      boundary_conditions = [["remove","thermalize"],["periodic"]],
      thermal_boundary_temperature = [T0],
      time_frozen = 0.0
    )

  Species(
      name = "eon",
      position_initialization = "ion",
      momentum_initialization = "mj",
      particles_per_cell = nppc,
      mass = 1.,
      charge = -1.,
      number_density = density_profile,
      temperature = [T0],
      boundary_conditions = [["remove","thermalize"],["periodic"]],
      thermal_boundary_temperature = [T0],
      time_frozen = 0.0
    )



  ################## LASER PULSE ######################################

  ### LASER PROPERTIES
  aL=0.1                          #normalized laser amplitude a0=0.86 sqrt(I/10^18)
  theta_deg = 30.                 #angle of incidence in degrees ArcSin[Sqrt[(n/gammal - 1)/(n/gammal - 2)] - 1 l/d]*180/Pi?
  theta  = theta_deg*np.pi/180.   #angle of incidence in rads

  waist  = waist0     # laser waist
  Tfwhm  = 10.*t0     # laser fwhm

  xfoc = Xsurface    # the laser is focalized over the surface
  yfoc = Ly/2.       # the laser is focalized at the target's center

  ### LASER PULSE
  LaserGaussian2D(
      a0              = aL,
      omega           = 1.,
      focus           = [xfoc,yfoc],
      waist           = waist,
      incidence_angle = theta,
      time_envelope   = tgaussian(start=0., duration=4.*Tfwhm, fwhm=Tfwhm, center=2.*Tfwhm, order=2)
  )


  #################### DIAGNOSTICS ####################################
  sub_grid=2 #to reduce the number of points in the output file
  Period = int(t0/dt)

  DiagFields(
      every = Period,
      fields = ['Bz','Ey','Ex', 'Rho_eon' ],
      subgrid = np.s_[::sub_grid, ::sub_grid]
  )

  DiagParticleBinning(
      deposited_quantity = "weight",
      every = [0, Period],
      time_average = 1,
      species = ["eon"],
      axes = [
          ["y", 0, Ly, 400],
          ["py",-0.3,0.3, 120]
      ]
  )

.. _Exercise 11:
**- Exercise 11:** Compare the longitudinal electric field at instant of time `TSim`,
assuming the three different configurations studied in this TP. 
Hint: To analyze the field on the surface, you can use the command ::

  import happi
  import numpy as np
  import matplotlib.pyplot as plt

  S = happi.Open('.', verbose=False)
  tver=60
  timesteps=S.Field(0,"Bz").getTimesteps()
  Xsurface=S.namelist.Xsurface-3*S.namelist.dx
  v0=S.namelist.aL
  diag  =  S.Field(0,"Bz", timesteps=timesteps[tver], average = {"x":Xsurface}).plot(vmin=-v0,vmax=v0)

And to put the three curves on the same panel, you can use the multiplot command from the happi library (see useful tools).

.. _Exercise 12:
**- Exercise 12:** Run three new simulations ( one for each target configuration ) assuming `aL=5`. Compare the transverse and longitudinal momentum of electrons at the instant of time `TSim`,
assuming the three different configurations studied in this TP. Use the phase space of each case for comparison and explain the result you obtain. Before run the simulations, do not forget modify the observation range in the Diagnostics.


Conclusion
==========

In this pratical work, we have explored the laser-plasma coupling in the relativistic regime of interaction. 
In particular, we have investigated conditions to improve laser-plasma energy transfer as well
as to generate relativistic electrons through computational
simulations using the Smilei PIC code.



References
^^^^^^^^^^

.. [smilei] `J. Derouillat, A. Beck, F. Pérez, T. Vinci, M. Chiaramello, A. Grassi, M. Flé, G. Bouchard, I. Plotnikov, N. Aunai, J. Dargent, C. Riconda, and M. Grech, Smilei: A collaborative, open-source, multi-purpose particle-in-cell code for plasma simulation, Comput. Phys. Commun. 222, 351 (2018).`

.. [Raynaud] `M. Raynaud, J. Kupersztych, C. Riconda, J. C. Adam, A.Heron,  Strongly enhanced laser absorption and electron acceleration via resonant excitation of surface plasma waves, Physics of Plasmas 14, 092702 (2007).``

.. [Marini] `S. Marini, P. Kleij, F. Amiranoff, M. Grech, M. Raynaud, and C. Riconda, Key parameters for surface plasma wave excitation in the ultra-high intensity regime. Physics of Plasmas 28, 073104 (2021).`


