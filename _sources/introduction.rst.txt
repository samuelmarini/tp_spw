=======================================================================================
Introduction
=======================================================================================




.. :Author: Samuel
:Date: February 2023

Observations of surface plasma waves (SPW) using the particle-in-cell (PIC) code SMILEI
===============

The particle-in-cell (PIC) method is a powerful computational tool for
studying plasma physics, allowing scientists to validate theoretical
predictions and understand experimental results. Its importance has
increased significantly in recent years, due to advancements in
computational systems and algorithm optimization. In this  practical work,
with the aim of understanding PIC computational algorithms, we will investigate the
excitation of a surface plasma wave (SPW) using the PIC code Smilei.

Surface plasma waves are electromagnetic waves confined at the interface
of two media with opposite permittivity, `e.g.` the vacuum-plasma
interface. These waves are excited when an electromagnetic wave like a
laser pulse is irradiated onto an overdense plasma target containing a corrugated
surface, known as gratings. When the laser pulse is directed onto the
plasma target at the resonant angle, the resulting diffracted wave can
excite the surface plasma wave. 
The surface plasma waves have been studied particularly in the context of
particle (electron) acceleration, as they exhibit intense longitudinal electric
field and subluminal phase velocity, which enables the electrons
acceleration to energies of tens of MeV, when assuming a laser pulse with
intensity of approximately :math:`\sim 10^{19} W/m^2`.

In this way, to observe the generation of energetic electrons
through surface plasma waves, we will collaborate to
build a Python script to be interpreted by the Smilei PIC code. In the
script, we need to define the characteristics of our simulation box and
describe how the laser-plasma setup should be arranged in order to
excite the surface plasma wave. It all based on the physical theory.

To build the script to be interpreted by the Smilei PIC code, 
we will go through several steps that will help
you even in the future construction and modeling of different physical
phenomena using PIC simulations. These steps include understanding and describing the laser
pulse, adding and understanding the plasma target, and generating
diagnostics that allow for the identification of the physical phenomena
we are interested in.

On the Useful Tools tab, you will find everything you need to know to prepare
your computer to run the computational simulations with the Smilei PIC
code. There you will also find some useful commands and scripts for
Diagnostics generation on SPW. On the Handout tab, you will find the
exercises developed for this TP so you can understand plasma waves as
well as you can understand PIC algorithms to use in your future
projects.

Get ready!
