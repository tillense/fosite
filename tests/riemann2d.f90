!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# module: riemann2d.f90                                                     #
!#                                                                           #
!# Copyright (C) 2006-2012                                                   #
!# Tobias Illenseer <tillense@astrophysik.uni-kiel.de>                       #
!#                                                                           #
!# This program is free software; you can redistribute it and/or modify      #
!# it under the terms of the GNU General Public License as published by      #
!# the Free Software Foundation; either version 2 of the License, or (at     #
!# your option) any later version.                                           #
!#                                                                           #
!# This program is distributed in the hope that it will be useful, but       #
!# WITHOUT ANY WARRANTY; without even the implied warranty of                #
!# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, GOOD TITLE or        #
!# NON INFRINGEMENT.  See the GNU General Public License for more            #
!# details.                                                                  #
!#                                                                           #
!# You should have received a copy of the GNU General Public License         #
!# along with this program; if not, write to the Free Software               #
!# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
!#                                                                           #
!#############################################################################
!>
!! \test 19 different 2D Riemann problems for the inviscid Euler equations
!! \author Tobias Illenseer
!!
!! <div class="row"> <div class="col-md-6">
!! ### Test description
!!
!! Numerical solutions of the inviscid Euler equations in 2D for the
!! 19 tests proposed by Schulz-Rinne et al. \cite schulz1993 . These 2D
!! Riemann problems were also used by several other authors (\cite lax1998 ,
!! \cite kurganov2002 , \cite liska2003 ) and became a standard collection of tests
!! for compressible flow solvers. \n 
!! The test examines the shock/contact discontinuity and rarefraction wave interaction
!! in two dimensions subdividing a quadratic computational domain into four equally
!! sized quadrants.  The computations were
!! carried out on the unit square \f$ \left[-0.5,0.5 \right]\times \left[ -0.5,0.5 \right] \f$ 
!! centered on the origin of a planar cartesian grid.
!! In all tests we assume an ideal gas equation of state with a constant ratio of specific heats
!! (adiabatic index \f$\gamma \f$).
!! </div> <div class="col-md-6">
!! Simulation parameters     ||
!! --------------------------|------------------
!! adiabatic index           | \f$ \gamma = 1.4 \f$
!! reconstruction            | conservative
!! limiter                   | monocent
!! diffusivity               | \f$ \theta = 1.2 \f$
!! numerical flux functions  | KT (Kurganov-Tadmor)
!! boundary conditions       | absorbing
!!
!! </div> </div>
!! <div class="row"> <div class="col-md-4">
!! ### Initial conditions
!!
!! The initial setup provides a distinct and constant state vector in each quadrant.
!! To distinguish the different configurations we use the numbering scheme given in
!! \cite kurganov2002 .
!! <div class="row"> <div class="col-md-6" align="center">
!!  | Quadrant | numbers |
!!  |:-----: |:-----:|
!!  |   2    |   1   |
!!  |   3    |   4   |
!! </div> </div>  
!!  The panel below shows the results for the density obtained for the 19 different
!!  test cases. The resolution was always set to \f$ 400 \times 400 \f$. The diagrams
!!  were also genereted in SVG format. Thus if you want to zoom in and study the details
!!  klick on the link below.
!! </div> <div class="col-md-4">
!! <div class="row"> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_01_thumb.png ""
!!   \ref conf01 "Configuration No. 1"
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_05_thumb.png ""
!!   \ref conf05 "Configuration No. 5"
!! </div> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_02_thumb.png ""
!!   \ref conf02 "Configuration No. 2"
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_06_thumb.png ""
!!   \ref conf06 "Configuration No. 6"
!! </div> </div>
!! </div> <div class="col-md-4">
!! <div class="row"> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_03_thumb.png ""
!!   \ref conf03 "Configuration No. 3"
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_07_thumb.png ""
!!   \ref conf07 "Configuration No. 7"
!! </div> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_04_thumb.png ""
!!   \ref conf04 "Configuration No. 4"
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_08_thumb.png ""
!!   \ref conf08 "Configuration No. 8"
!! </div> </div>
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-4" align="center">
!! <div class="row"> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_01_thumb.png ""
!!   \ref conf09 "Configuration No. 9"
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_15_thumb.png ""
!!   \ref conf15 "Configuration No. 15"
!! </div> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_10_thumb.png ""
!!   \ref conf10 "Configuration No. 10"
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_16_thumb.png ""
!!   \ref conf16 "Configuration No. 16"
!! </div> </div>
!! </div> <div class="col-md-4">
!! <div class="row"> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_11_thumb.png ""
!!   \ref conf11 "Configuration No. 11"
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_17_thumb.png ""
!!   \ref conf17 "Configuration No. 17"
!! </div> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_12_thumb.png ""
!!   \ref conf12 "Configuration No. 12"
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_18_thumb.png ""
!!   \ref conf18 "Configuration No. 18"
!! </div> </div>
!! </div> <div class="col-md-4">
!! <div class="row"> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_13_thumb.png ""
!!   \ref conf13 "Configuration No. 13"
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_19_thumb.png ""
!!   \ref conf19 "Configuration No. 19"
!! </div> <div class="col-md-6" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_14_thumb.png ""
!!   \ref conf14 "Configuration No. 14"
!! </div> </div>
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 1 \anchor conf01
!!
!! The initial condition yields four rarefactions along the
!! interfaces of the four quadrants.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.0         |0.0         |1.0  
!!     2    | 0.5197  |-0.7259     |0.0         |0.4  
!!     3    | 0.1072  |-0.7259     |-1.4045     |0.0439
!!     4    | 0.2579  |0.0         |-1.4045     |0.15 
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_01_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_01_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 2 \anchor conf02
!!
!! The initial condition yields four rarefactions along the
!! interfaces of the four quadrants.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.0         |0.0         |1.0  
!!     2    | 0.5197  |-0.7259     |0.0         |0.4  
!!     3    | 1.0     |-0.7259     |-0.7259     |1.0
!!     4    | 0.5197  |0.0         |-0.7259     |0.4  
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_02_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_02_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 3 \anchor conf03
!!
!! The initial condition yields four shocks along the
!! interfaces of the four quadrants.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.5     |0.0         |0.0         |1.5
!!     2    | 0.5323  |1.2060      |0.0         |0.3
!!     3    | 0.1380  |1.2060      |1.2060      |0.029
!!     4    | 0.5323  |0.0         |1.2060      |0.3
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_03_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_03_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 4 \anchor conf04
!!
!! The initial condition yields four shocks along the
!! interfaces of the four quadrants.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.1     |0.0         |0.0         |1.1  
!!     2    | 0.5065  |0.8939      |0.0         |0.35
!!     3    | 1.1     |0.8939      |0.8939      |1.1  
!!     4    | 0.5065  |0.0         |0.8939      |0.35
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_04_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_04_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 5 \anchor conf05
!!
!! The initial condition yields four contacts along the
!! interfaces of the four quadrants.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |-0.75       |-0.5        |1.0  
!!     2    | 2.0     |-0.75       |0.5         |1.0  
!!     3    | 1.0     |0.75        |0.5         |1.0  
!!     4    | 3.0     |0.75        |-0.5        |1.0  
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_05_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_05_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 6 \anchor conf06
!!
!! The initial condition yields four contacts along the
!! interfaces of the four quadrants.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.75        |-0.5        |1.0  
!!     2    | 2.0     |0.75        |0.5         |1.0  
!!     3    | 1.0     |-0.75       |0.5         |1.0  
!!     4    | 3.0     |-0.75       |-0.5        |1.0  
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_06_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_06_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 7 \anchor conf07
!!
!! The initial condition contains two rarefactions at the upper and the
!! right interfaces and two contacts at the lower and left interfaces.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.1         |0.1         |1.0  
!!     2    | 0.5197  |-0.6259     |0.1         |0.4  
!!     3    | 0.8     |0.1         |0.1         |0.4   
!!     4    | 0.5197  |0.1         |-0.6259     |0.4   
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_07_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_07_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 8 \anchor conf08
!!
!! The initial condition ontains two rarefactions at the upper and the
!! right interfaces and two contacts at the lower and left interfaces.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 0.5197  |0.1         |0.1         |0.4  
!!     2    | 1.0     |-0.6259     |0.1         |1.0  
!!     3    | 0.8     |0.1         |0.1         |1.0   
!!     4    | 1.0     |0.1         |-0.6259     |1.0  
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_08_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_08_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 9 \anchor conf09
!!
!! The initial condition ycontains two rarefactions at the left and right
!! interfaces and two contacts at the upper and lower interfaces. 
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.0         |0.3         |1.0  
!!     2    | 1.0     |0.0         |-0.3        |1.0  
!!     3    | 1.0390  |0.0         |-0.8133     |0.4
!!     4    | 0.5197  |0.0         |-0.4259     |0.4
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_09_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_09_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 10 \anchor conf10
!!
!! The initial condition contains two rarefactions at the left and right
!! interfaces and two contacts at the upper and lower interfaces. 
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.0         |0.4297      |1.0  
!!     2    | 0.5     |0.0         |0.6076      |1.0     
!!     3    | 0.2281  |0.0         |-0.6076     |0.3333  
!!     4    | 0.4562  |0.0         |-0.4297     |0.3333  
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_10_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_10_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 11 \anchor conf11
!!
!! The initial condition contains two shocks at the upper and right
!! interfaces and two contacts at the lower and left interfaces. 
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.1         |0.0         |1.0  
!!     2    | 0.5313  |0.8276      |0.0         |0.4  
!!     3    | 0.8     |0.1         |0.0         |0.4    
!!     4    | 0.5313  |0.1         |0.7276      |0.4    
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_11_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_11_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 12 \anchor conf12
!!
!! The initial condition contains two shocks at the upper and right
!! interfaces and two contacts at the lower and left interfaces.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.0         |0.0         |0.4  
!!     2    | 0.5313  |0.7276      |0.0         |1.0  
!!     3    | 0.8     |0.0         |0.0         |1.0   
!!     4    | 1.0     |0.0         |0.7276      |1.0   
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_12_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_12_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 13 \anchor conf13
!!
!! The initial condition contains two shocks at the left and right
!! interfaces and two contacts at the upper and lower interfaces.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.0         |-0.3        |1.0  
!!     2    | 2.0     |0.0         |0.3         |1.0  
!!     3    | 1.0625  |0.0         |0.8145      |0.4   
!!     4    | 0.5313  |0.0         |0.4276      |0.4  
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_13_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_13_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 14 \anchor conf14
!!
!! The initial condition contains two shocks at the left and right
!! interfaces and two contacts at the upper and lower interfaces.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 2.0     |0.0         |-0.5606     |8.0  
!!     2    | 1.0     |0.0         |-1.2172     |8.0  
!!     3    | 0.4736  |0.0         |1.2172      |2.6667
!!     4    | 0.9474  |0.0         |1.1606      |2.6667
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_14_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_14_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 15 \anchor conf15
!!
!! The initial condition contains two contacts at the lower and left interfaces,
!! a rarefaction at the upper and shock at the right interface.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.1         |-0.3        |1.0  
!!     2    | 0.5197  |-0.6259     |-0.3        |0.4  
!!     3    | 0.8     |0.1         |-0.3        |0.4   
!!     4    | 0.5313  |0.1         |0.4276      |0.4   
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_15_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_15_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 16 \anchor conf16
!!
!! The initial condition contains two contacts at the lower and left interfaces,
!! a rarefaction at the upper and shock at the right interface.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 0.5313  |0.1         |0.1         |0.4  
!!     2    | 1.02222 |-0.6179     |0.1         |1.0  
!!     3    | 0.8     |0.1         |0.1         |1.0   
!!     4    | 1.0     |0.1         |0.8276      |1.0  
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_16_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_16_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 17 \anchor conf17
!!
!! The initial conditioncontains two contacts at the upper and lower interfaces,
!! a shock at the left and rarefaction at the right interface.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.0         |-0.4        |1.0  
!!     2    | 2.0     |0.0         |-0.3        |1.0  
!!     3    | 1.0625  |0.0         |0.2145      |0.4   
!!     4    | 0.5197  |0.0         |-1.1259     |0.4  
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_17_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_17_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 18 \anchor conf18
!!
!! The initial condition contains two contacts at the upper and lower interfaces,
!! a shock at the left and rarefaction at the right interface.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.0         |1.0         |1.0  
!!     2    | 2.0     |0.0         |-0.3        |1.0    
!!     3    | 1.0625  |0.0         |0.2145      |0.4    
!!     4    | 0.5197  |0.0         |0.2741      |0.4    
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_18_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_18_contour.svg)
!! </div> </div>
!!
!! <div class="row"> <div class="col-md-6">
!! #### Configuration No. 19 \anchor conf19
!!
!! The initial condition contains two contacts at the upper and lower interfaces,
!! a shock at the left and rarefaction at the right interface.
!! Quadrant | Density | X-Velocity | Y-Velocity | Pressure 
!! :-------:|:--------|:-----------|:-----------|:--------
!!     1    | 1.0     |0.0         |0.3         |1.0  
!!     2    | 2.0     |0.0         |-0.3        |1.0    
!!     3    | 1.0625  |0.0         |0.2145      |0.4    
!!     4    | 0.5197  |0.0         |-0.4259     |0.4    
!!
!! </div> <div class="col-md-6" align = "center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_19_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_19_contour.svg)
!! </div> </div>
!!
!! ### References
!!
!! - \cite schulz1993 C. W. Schulz-Rinne et al.: Numerical Solution of the Riemann Problem
!!     for Gas Dynamics, SIAM J. Sci. Comp. 14 (1993), 1394-1414
!! - \cite lax1998  P. Lax, X.-D. Liu: Solution of Two-dimensional Riemann Problems of
!!     Gas Dynamics by Positive Schemes, SIAM J. Sci. Comp. 19 (1998), 319-340
!! - \cite kurganov2002 A. Kurganov, E. Tadmor: Solution of Two-Dimensional Riemann Problems for
!!     Gas Dynamics without Riemann Problem Solvers, NMPDE 18 (2002), 561-588
!!
PROGRAM riemann2d
  USE fosite_mod
#include "tap.h"
  IMPLICIT NONE
  !--------------------------------------------------------------------------!
  ! simulation parameter
  INTEGER, PARAMETER :: ICNUM = 19        ! initial condition (see ref. [3])
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry of the mesh
!!$  INTEGER, PARAMETER :: MGEO = POLAR
!!$  INTEGER, PARAMETER :: MGEO = LOGPOLAR
!!$  INTEGER, PARAMETER :: MGEO = TANPOLAR
!!$  INTEGER, PARAMETER :: MGEO = SINHPOLAR
  INTEGER, PARAMETER :: XRES = 100         ! resolution
  INTEGER, PARAMETER :: YRES = 100
  INTEGER, PARAMETER :: ZRES = 1
  REAL, PARAMETER    :: RMIN = 1.0E-4      ! inner radius for polar grids
  ! output file parameter
  INTEGER, PARAMETER :: ONUM = 1           ! number of output data sets
  CHARACTER(LEN=256), PARAMETER &          ! output data dir
                     :: ODIR = './'
  CHARACTER(LEN=256), PARAMETER &          ! output data file name
                     :: OFNAME = 'riemann2d'
  !--------------------------------------------------------------------------!
  CLASS(fosite),Dimension(:), ALLOCATABLE      :: Sim
  INTEGER                         :: ic
  !--------------------------------------------------------------------------!

  TAP_PLAN(ICNUM)

  ALLOCATE(Sim(ICNUM))
  

  
  DO ic=1,ICNUM

    CALL Sim(ic)%InitFosite()

     CALL MakeConfig(Sim(ic), Sim(ic)%config, ic)

!  CALL PrintDict(config)

     CALL Sim(ic)%Setup()

     ! set initial condition
     CALL InitData(Sim(ic)%Mesh, Sim(ic)%Physics, Sim(ic)%Timedisc, ic)
     CALL Sim(ic)%Run()

     CALL Sim(ic)%Finalize()
     TAP_CHECK(.TRUE.,"Simulation finished")

  END DO

  DEALLOCATE(Sim)

  TAP_DONE

CONTAINS

  SUBROUTINE MakeConfig(Sim, config, ic)
    USE functions, ONLY : Asinh
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Fosite)     :: Sim
    TYPE(Dict_TYP),POINTER :: config
    INTEGER           :: ic
    !------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER           :: bc(6)
    TYPE(Dict_TYP), POINTER :: mesh, physics, boundary, datafile, &
                               timedisc, fluxes
    REAL              :: x1,x2,y1,y2,z1,z2,sc
    REAL              :: test_stoptime
    CHARACTER(LEN=3)  :: fext
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    INTENT(IN)        :: ic
    !------------------------------------------------------------------------!
    ! mesh settings and boundary conditions
    SELECT CASE(MGEO)
    CASE(CARTESIAN)
       sc = 1.0
       x1 = -0.5
       x2 =  0.5
       y1 = -0.5
       y2 =  0.5
       z1 = -0.0
       z2 =  0.0
       bc(WEST)  = NO_GRADIENTS
       bc(EAST)  = NO_GRADIENTS
       bc(SOUTH) = NO_GRADIENTS
       bc(NORTH) = NO_GRADIENTS
       bc(BOTTOM)= NO_GRADIENTS
       bc(TOP)   = NO_GRADIENTS
!    CASE(POLAR)
!       sc = 1.0
!       x1 = RMIN
!       x2 = 0.5*SQRT(2.0)
!       y1 = 0.0
!       y2 = 2*PI
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
!       bc(SOUTH) = PERIODIC
!       bc(NORTH) = PERIODIC
!    CASE(LOGPOLAR)
!       sc = 0.3
!       x1 = LOG(RMIN/sc)
!       x2 = LOG(0.5*SQRT(2.0)/sc)
!       y1 = 0.0
!       y2 = 2*PI
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
!       bc(SOUTH) = PERIODIC
!       bc(NORTH) = PERIODIC
!    CASE(TANPOLAR)
!       sc = 0.3
!       x1 = ATAN(RMIN/sc)
!       x2 = ATAN(0.5*SQRT(2.0)/sc)
!       y1 = 0.0
!       y2 = 2*PI
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
!       bc(SOUTH) = PERIODIC
!       bc(NORTH) = PERIODIC
!    CASE(SINHPOLAR)
!       sc = 0.3
!       x1 = Asinh(RMIN/sc)
!       x2 = Asinh(0.5*SQRT(2.0)/sc)
!       y1 = 0.0
!       y2 = 2*PI
!       bc(WEST)  = NO_GRADIENTS
!       bc(EAST)  = NO_GRADIENTS
!       bc(SOUTH) = PERIODIC
!       bc(NORTH) = PERIODIC
    CASE DEFAULT
       CALL Sim%Physics%Error("riemann2d::MakeConfig", &
            " geometry should be one of cartesian,polar,logpolar,tanpolar,sinhpolar")
    END SELECT

    ! mesh settings
    mesh => Dict("meshtype" / MIDPOINT, &
           "geometry" / MGEO, &
           "inum"     / XRES, &
           "jnum"     / YRES, &
           "knum"     / ZRES, &
           "xmin"     / x1, &
           "xmax"     / x2, &
           "ymin"     / y1, &
           "ymax"     / y2, &
           "zmin"     / z1, &
           "zmax"     / z2, &
           "gparam"   / sc)

    ! boundary conditions
    boundary => Dict("western" / bc(WEST), &
               "eastern" / bc(EAST), &
               "southern" / bc(SOUTH), &
               "northern" / bc(NORTH), &
               "bottomer" / bc(BOTTOM),&
               "topper"   / bc(TOP))

    ! physics settings
    physics => Dict("problem" / EULER, &
              "gamma"   / GAMMA)         ! ratio of specific heats        !

    ! flux calculation and reconstruction method
    fluxes => Dict("order"     / LINEAR, &
             "fluxtype"  / KT, &
             "variables" / CONSERVATIVE, &        ! vars. to use for reconstruction!
             "limiter"   / MONOCENT, &    ! one of: minmod, monocent,...   !
             "theta"     / 1.2)          ! optional parameter for limiter !

    ! runtime of the test problem
    SELECT CASE(ic)
    CASE(1)  ! KT test 1
       test_stoptime = 0.2
    CASE(2)  ! Riemann problem no. 2
       test_stoptime = 0.2
    CASE(3)  ! Riemann problem no. 3
       test_stoptime = 0.3
    CASE(4)  ! Riemann problem no. 4
       test_stoptime = 0.25
    CASE(5)  ! Riemann problem no. 5
       test_stoptime = 0.23
    CASE(6)  ! Riemann problem no. 6
       test_stoptime = 0.3
    CASE(7)  ! Riemann problem no. 7
       test_stoptime = 0.25
    CASE(8)  ! Riemann problem no. 8
       test_stoptime = 0.25
    CASE(9)  ! Riemann problem no. 9
       test_stoptime = 0.3
    CASE(10)  ! Riemann problem no. 10
       test_stoptime = 0.15
    CASE(11)  ! Riemann problem no. 11
       test_stoptime = 0.3
    CASE(12) ! Riemann problem no. 12
       test_stoptime = 0.25
    CASE(13)  ! Riemann problem no. 13
       test_stoptime = 0.3
    CASE(14)  ! Riemann problem no. 14
       test_stoptime = 0.1
    CASE(15) ! Riemann problem no. 15
       test_stoptime = 0.2
    CASE(16)  ! Riemann problem no. 16
       test_stoptime = 0.2
    CASE(17) ! Riemann problem no. 17
       test_stoptime = 0.3
    CASE(18) ! Riemann problem no. 18
       test_stoptime = 0.2
    CASE(19) ! Riemann problem no. 19
       test_stoptime = 0.3
    CASE DEFAULT
       CALL Sim%Physics%Error("InitProgram", &
            "Sorry, this 2D Riemann problem is currently not supported!")
    END SELECT

    ! time discretization settings
    timedisc => Dict( &
           "method"   / MODIFIED_EULER, &
           "order"    / 3, &
           "cfl"      / 0.4, &
           "stoptime" / test_stoptime, &
           "dtlimit"  / 1.0E-10, &
           "maxiter"  / 10000000)

    ! initialize data input/output
    ! append test number to file names
    WRITE (fext, '(A,I2.2)') "_", ic
!    datafile => Dict("fileformat" / VTK, &
    datafile => Dict( &
             !   "fileformat" / GNUPLOT, "filecycles" / 0, &
                "fileformat" / VTK, &
                "filename"   / (TRIM(ODIR) // TRIM(OFNAME) // TRIM(fext)), &
                "count"      / ONUM)

    config => Dict("mesh" / mesh, &
             "physics"  / physics, &
             "boundary" / boundary, &
             "fluxes"   / fluxes, &
             "timedisc" / timedisc, &
!             "logfile"  / logfile, &
             "datafile" / datafile)
  END SUBROUTINE MakeConfig

  !> Set initial conditions
  SUBROUTINE InitData(Mesh,Physics,Timedisc,ic)
  IMPLICIT NONE
    !------------------------------------------------------------------------!
    CLASS(Physics_base) :: Physics
    CLASS(Mesh_base)    :: Mesh
    CLASS(Timedisc_base):: Timedisc
    INTEGER           :: ic
    !------------------------------------------------------------------------!
    ! Local variable declaration
    REAL, DIMENSION(Mesh%IGMIN:Mesh%IGMAX,Mesh%JGMIN:Mesh%JGMAX, &
                                          Mesh%KGMIN:Mesh%KGMAX,2) :: vxy
    REAL              :: xmin,ymin,xmax,ymax,x0,y0
#ifdef PARALLEL
    REAL              :: xmin_all,xmax_all,ymin_all,ymax_all
    INTEGER           :: ierr
#endif
    CHARACTER(LEN=64) :: teststr
    !------------------------------------------------------------------------!
    INTENT(IN)        :: Mesh,Physics,ic
    INTENT(INOUT)     :: Timedisc
    !------------------------------------------------------------------------!
    ! minima and maxima of _cartesian_ coordinates
    xmin = MINVAL(Mesh%bccart(:,:,:,1))
    ymin = MINVAL(Mesh%bccart(:,:,:,2))
    xmax = MAXVAL(Mesh%bccart(:,:,:,1))
    ymax = MAXVAL(Mesh%bccart(:,:,:,2))
#ifdef PARALLEL
    CALL MPI_Allreduce(xmin,xmin_all,1,DEFAULT_MPI_REAL,MPI_MIN,Mesh%comm_cart,ierr)
    xmin = xmin_all
    CALL MPI_Allreduce(ymin,ymin_all,1,DEFAULT_MPI_REAL,MPI_MIN,Mesh%comm_cart,ierr)
    ymin = ymin_all
    CALL MPI_Allreduce(xmax,xmax_all,1,DEFAULT_MPI_REAL,MPI_MAX,Mesh%comm_cart,ierr)
    xmax = xmax_all
    CALL MPI_Allreduce(ymax,ymax_all,1,DEFAULT_MPI_REAL,MPI_MAX,Mesh%comm_cart,ierr)
    ymax = ymax_all
#endif
    x0 = xmin + 0.5*ABS(xmax-xmin)
    y0 = ymin + 0.5*ABS(ymax-ymin)

    Timedisc%pvar%data4d(:,:,:,:) = 0.
    vxy(:,:,:,:) = 0.

    SELECT CASE(ic)
    CASE(1)
       teststr = "2D Riemann problem no. 1"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = .5197
          vxy(:,:,:,1) = -.7259
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = .4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = .1072
          vxy(:,:,:,1) = -.7259
          vxy(:,:,:,2) = -1.4045
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = .0439
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = .2579
          vxy(:,:,:,2) = -1.4045
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = .15
       END WHERE
    
    CASE(2)
       teststr = "2D Riemann problem no. 2"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = .5197
          vxy(:,:,:,1) = -.7259
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = .4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = -.7259
          vxy(:,:,:,2) = -.7259
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = .5197
          vxy(:,:,:,2) = -.7259
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = .4
       END WHERE

    CASE(3)
       teststr = "2D Riemann problem no. 3"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.5
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.5
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = .5323
          vxy(:,:,:,1) = 1.206
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = .3
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.138
          vxy(:,:,:,1) = 1.206
          vxy(:,:,:,2) = 1.206
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.029
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = .5323
          vxy(:,:,:,2) = 1.206
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = .3
       END WHERE

    CASE(4)
       teststr = "2D Riemann problem no. 4"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.1
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = .5065
          vxy(:,:,:,1) = .8939
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = .35
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.1
          vxy(:,:,:,1) = .8939
          vxy(:,:,:,2) = .8939
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.1
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = .5065
          vxy(:,:,:,2) = .8939
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = .35
       END WHERE

    CASE(5)
       teststr = "2D Riemann problem no. 5"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = -.75
          vxy(:,:,:,2) = -.5
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 2.
          vxy(:,:,:,1) = -.75
          vxy(:,:,:,2) = .5
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = .75
          vxy(:,:,:,2) = .5
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 3.
          vxy(:,:,:,1) = .75
          vxy(:,:,:,2) = -.5
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       END WHERE
    CASE(6)
       teststr = "2D Riemann problem no. 6"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.75
          vxy(:,:,:,2) = -0.5
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 2.
          vxy(:,:,:,1) = 0.75
          vxy(:,:,:,2) = 0.5
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = -0.75
          vxy(:,:,:,2) = 0.5
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 3.
          vxy(:,:,:,1) = -0.75
          vxy(:,:,:,2) = -0.5
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       END WHERE
    CASE(7)
       teststr = "2D Riemann problem no. 7"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,:,1) = -0.6259
          vxy(:,:,:,2) = 0.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.8
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = -0.6259
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       END WHERE
    CASE(8)
       teststr = "2D Riemann problem no. 8"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = -0.6259
          vxy(:,:,:,2) = 0.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.8
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = -0.6259
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       END WHERE
    CASE(9)
       teststr = "2D Riemann problem no. 9"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 2.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.039
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.8133
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.4259
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       END WHERE
    CASE(10)
       teststr = "2D Riemann problem no. 10"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.4297
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.6076
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.2281
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.6076
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.3333
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.4562
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.4297
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.3333
       END WHERE
    CASE(11)
       teststr = "2D Riemann problem no. 11"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5313
          vxy(:,:,:,1) = 0.8276
          vxy(:,:,:,2) = 0.
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.8
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5313
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.7276
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       END WHERE
    CASE(12)
       teststr = "2D Riemann problem no. 12"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5313
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.7276
          vxy(:,:,:,2) = 0.
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.8
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.7276
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       END WHERE
    CASE(13)
       teststr = "2D Riemann problem no. 13"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 2.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.0625
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.8145
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5313
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.4276
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       END WHERE
    CASE(14)
       teststr = "2D Riemann problem no. 14"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 2.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.5606
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 8.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -1.2172
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 8.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.4736
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 1.2172
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 2.6667
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.9474
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 1.1606
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 2.6667
       END WHERE
    CASE(15)
       teststr = "2D Riemann problem no. 15"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = -0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,:,1) = -0.6259
          vxy(:,:,:,2) = -0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.8
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = -0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5313
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.4276
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       END WHERE
    CASE(16)
       teststr = "2D Riemann problem no. 16"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5313
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.02222
          vxy(:,:,:,1) = -0.6179
          vxy(:,:,:,2) = 0.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.8
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.1
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.1
          vxy(:,:,:,2) = 0.8276
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       END WHERE
    CASE(17)
       teststr = "2D Riemann problem no. 17"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.4
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 2.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.0625
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.2145
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -1.1259
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       END WHERE

    CASE(18)
       teststr = "2D Riemann problem no. 18"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 1.
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 2.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.0625
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.2145
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.2741
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       END WHERE

    CASE(19)
       teststr = "2D Riemann problem no. 19"
       WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 1
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
          ! no. 2
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 2.
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.3
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 1.
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 3
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 1.0625
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = 0.2145
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
          ! no. 4
          Timedisc%pvar%data4d(:,:,:,Physics%DENSITY) = 0.5197
          vxy(:,:,:,1) = 0.
          vxy(:,:,:,2) = -0.4259
          Timedisc%pvar%data4d(:,:,:,Physics%PRESSURE) = 0.4
       END WHERE
  
    CASE DEFAULT
       CALL Mesh%Error("InitData", &
            "Sorry, this 2D Riemann problem is currently not supported!")
    END SELECT
   
!    CALL Mesh%geometry%Convert2Curvilinear(Mesh%bcenter,vxy,&
!         Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY:Physics%YVELOCITY))
   Timedisc%pvar%data4d(:,:,:,Physics%XVELOCITY) = vxy(:,:,:,1)
   Timedisc%pvar%data4d(:,:,:,Physics%YVELOCITY) = vxy(:,:,:,2)
    CALL Physics%Convert2Conservative(Mesh,Timedisc%pvar%data4d,Timedisc%cvar%data4d)
    CALL Mesh%Info(" DATA-----> initial condition: " // TRIM(teststr))

  END SUBROUTINE InitData

END PROGRAM riemann2d
!> <div class="row"> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_01_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_01_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_02_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_02_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_03_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_03_contour.svg) \n
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_04_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_04_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_05_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_05_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_06_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_06_contour.svg) \n
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_07_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_07_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_08_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_08_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_09_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_09_contour.svg) \n
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_10_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_10_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_11_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_11_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_12_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_12_contour.svg) \n
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_13_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_13_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_14_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_14_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_15_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_15_contour.svg) \n
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_16_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_16_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_17_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_17_contour.svg)
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_18_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_18_contour.svg) \n
!! </div> <div class="col-md-4" align="center">
!!   \image html http://www.astrophysik.uni-kiel.de/fosite/riemann2d_19_contour.png ""
!!   [view  SVG] (http://www.astrophysik.uni-kiel.de/fosite/riemann2d_19_contour.svg)
!! </div> </div>
