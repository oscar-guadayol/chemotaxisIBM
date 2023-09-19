# chemotaxisIBM
is a collection of two matlab scripts, chemotaxisIBM_squareDomain.m and chemotaxisIBM_phycosphere.m, that simulate the behaviour of a population of Escherichia coli cellsswimming in 2D space in the presence of a gradient of non-metabolizable alpha-methyl-aspartate. The response of individual cells to the nutrients is modelled with a  that codifies the E. coli's chemosensory circuit. 

The chemotactic behaviour is modelled using a coarse-grain model of the chemosensory circuit of E. coli for aspartate, following Tu et al 2008 and Kalinin et al 2009.
The motility pattern is defined by the rotational diffusivity coefficient of the cells and by their tumble time and angle, running speed and run length. The model allows for multiple modes motility patterns, where the cells alternate between different reorientation behaviours.

 Copyright (C) 2023,  Oscar Guadayol <oscar_at_guadayol.cat>

 LICENSE:
  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License, version 3.0, as
  published by the Free Software Foundation.
 
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.
 
  You should have received a copy of the GNU General Public License along
  with this program.  If not, see <http://www.gnu.org/licenses/>.
