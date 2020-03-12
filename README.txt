-*-Mode: text;-*-

$HeadURL$
$Id$

=============================================================================
=============================================================================

This is a pre-release of the PERFECT Suite.  It should be distributed
only to participants in the DARPA PERFECT program.

The Suite has been jointly developed by Pacific Northwest National
Laboratory (PNNL) and Georgia Tech Research Institute (GTRI) with
support from the DARPA PERFECT program.

For more information on the PERFECT Suite see:
  http://hpc.pnl.gov/projects/PERFECT

People:
  Kevin Barker (PNNL)
  Thomas Benson (GTRI)
  Dan Campbell (GTRI)
  David Ediger (GTRI)
  Roberto Gioiosa (PNNL)
  Adolfy Hoisie (PNNL)
  Darren Kerbyson (PNNL)
  Joseph Manzano (PNNL)
  Andres Marquez (PNNL)
  Leon Song (PNNL)
  Nathan Tallent (PNNL)
  Antonino Tumeo (PNNL)

=============================================================================
=============================================================================

In this pre-release, the following directories are populated:
  doc/           : PERFECT Suite manual
  suite/pa1      : "PERFECT Application 1" kernels
  suite/stap/    : "Space-Time Adaptive Processing" kernels
  suite/sar/     : "Synthetic Aperture Radar" kernels
  suite/required/: Required kernels

The make system is based on GNU make.  Type "make help" to see a list of
supported targets.

To build and test the suite's kernels, issue the command
    make check

The make system defaults to GCC, but that can be changed by overriding
CC and CFLAGS, e.g.,
     # environment override
    CC=clang make
or
     # make override
    make CC=clang


=============================================================================
=============================================================================

PNNL's disclaimer appears below:

This computer software was prepared by Battelle Memorial Institute,
hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with
the Department of Energy (DOE).  All rights in the computer software
are reserved by DOE on behalf of the United States Government and the
Contractor as provided in the Contract.  You are authorized to use
this computer software for Governmental purposes but it is not to be
released or distributed to the public.
 
This material was prepared as an account of work sponsored by an
agency of the United States Government.  Neither the United States
Government nor the United States Department of Energy, nor the
Contractor, nor any of their employees, nor any jurisdiction or
organization that has cooperated in the development of these
materials, makes any warranty, express or implied, or assumes any
legal liability or responsibility for the accuracy, completeness, or
usefulness or any information, apparatus, product, software, or
process disclosed, or represents that its use would not infringe
privately owned rights.
 
Reference herein to any specific commercial product, process, or
service by trade name, trademark, manufacturer, or otherwise does not
necessarily constitute or imply its endorsement, recommendation, or
favoring by the United States Government or any agency thereof, or
Battelle Memorial Institute. The views and opinions of authors
expressed herein do not necessarily state or reflect those of the
United States Government or any agency thereof.

		PACIFIC NORTHWEST NATIONAL LABORATORY
			     operated by
			       BATTELLE
			       for the
		  UNITED STATES DEPARTMENT OF ENERGY
		   under Contract DE-AC05-76RL01830

