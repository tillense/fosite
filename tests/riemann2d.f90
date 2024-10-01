!#############################################################################
!#                                                                           #
!# fosite - 3D hydrodynamical simulation program                             #
!# module: riemann2d.f90                                                     #
!#                                                                           #
!# Copyright (C) 2006-2024                                                   #
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
  INTEGER, PARAMETER :: SELECT_TEST = 6   ! run this test, should be in 1..19
  LOGICAL, PARAMETER :: RUN_ALL = .TRUE.  ! set this if you want to run all 19 tests
  REAL, PARAMETER    :: GAMMA = 1.4        ! ratio of specific heats
  ! mesh settings
  INTEGER, PARAMETER :: MGEO = CARTESIAN   ! geometry of the mesh
!   INTEGER, PARAMETER :: MGEO = CYLINDRICAL
!   INTEGER, PARAMETER :: MGEO = LOGCYLINDRICAL
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
  ! do not change these parameters
  INTEGER, PARAMETER :: NUM_TESTS = 19     ! total number of 2D setups
  REAL, PARAMETER, DIMENSION(NUM_TESTS) :: & ! stop time for each test
                    TEST_STOPTIME = (/ 0.2,0.2,0.3,0.25,0.23,0.3,0.25, &
                                       0.25,0.3,0.15,0.3,0.25,0.3,0.1, &
                                       0.2,0.2,0.3,0.2,0.3/)
  !--------------------------------------------------------------------------!
  CLASS(fosite), ALLOCATABLE      :: Sim
  LOGICAL                         :: ok(NUM_TESTS), mpifinalize = .FALSE.
  INTEGER                         :: ic,err
  !--------------------------------------------------------------------------!

  IF (RUN_ALL) THEN
TAP_PLAN(NUM_TESTS)
  ELSE
TAP_PLAN(1)
  END IF

  DO ic=1,NUM_TESTS
    IF (.NOT.RUN_ALL.AND.ic.NE.SELECT_TEST) CYCLE ! skip if not selected

    IF (ALLOCATED(Sim)) DEALLOCATE(Sim)
    ALLOCATE(Sim,STAT=err)
    IF (err.NE.0) &
      CALL Sim%Error("riemann2d","cannot allocate memory for Sim")

    CALL Sim%InitFosite()
    CALL MakeConfig(Sim, Sim%config, ic)

!  CALL PrintDict(config)

    CALL Sim%Setup()

    ! set initial condition
    CALL InitData(Sim%Mesh, Sim%Physics, Sim%Timedisc, ic)
    CALL Sim%Run()
    ok(ic) = .NOT.Sim%aborted

    IF ((RUN_ALL.AND.ic.EQ.NUM_TESTS).OR.(.NOT.RUN_ALL)) mpifinalize=.TRUE.
    CALL Sim%Finalize(mpifinalize)

  END DO

  DO ic=1,NUM_TESTS
    IF (.NOT.RUN_ALL.AND.ic.NE.SELECT_TEST) CYCLE ! skip if not selected
TAP_CHECK(ok(ic),"Simulation finished")
  END DO

  IF (ALLOCATED(Sim)) DEALLOCATE(Sim)

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
    CHARACTER(LEN=16) :: fext
    !------------------------------------------------------------------------!
    INTENT(INOUT)     :: Sim
    INTENT(IN)        :: ic
    CHARACTER(LEN=8)  :: geo_str
    !------------------------------------------------------------------------!
    IF (ic.LT.1.OR.ic.GT.NUM_TESTS) &
       CALL Sim%Error("riemann2d::MakeConfig","invalid test number selected")

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
       geo_str = "-cart"
   CASE(CYLINDRICAL)
      sc = 1.0
      x1 = 0.0
      x2 = 0.5*SQRT(2.0)
      y1 = 0.0
      y2 = 2*PI
      z1 = -0.0
      z2 = 0.0
      bc(WEST)   = ABSORBING
      bc(EAST)   = NO_GRADIENTS
      bc(SOUTH)  = PERIODIC
      bc(NORTH)  = PERIODIC
      bc(BOTTOM) = REFLECTING
      bc(TOP)    = REFLECTING
      geo_str = "-cyl"
   CASE(LOGCYLINDRICAL)
      sc = 0.3
      x1 = LOG(RMIN/sc)
      x2 = LOG(0.5*SQRT(2.0)/sc)
      y1 = -PI
      y2 = PI
      z1 = 0.0
      z2 = 0.0
      bc(WEST)   = AXIS
      bc(EAST)   = NO_GRADIENTS
      bc(SOUTH)  = PERIODIC
      bc(NORTH)  = PERIODIC
      bc(BOTTOM) = NO_GRADIENTS
      bc(TOP)    = NO_GRADIENTS
      geo_str = "-logcyl"
    CASE DEFAULT
      CALL Sim%Error("riemann2d::MakeConfig","geometry not supported in this test")
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
           "output/commutator" / 1,&
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

    ! time discretization settings
    timedisc => Dict( &
           "method"   / MODIFIED_EULER, &
           "order"    / 3, &
           "cfl"      / 0.4, &
           "stoptime" / TEST_STOPTIME(ic), &
           "dtlimit"  / 1.0E-10, &
           "maxiter"  / 10000000, &
           "output/geometrical_sources" / 1)

    ! initialize data input/output
    ! append test number to file names
    WRITE (fext, '(A,I2.2)') TRIM(geo_str) // "_", ic
    datafile => Dict( &
!                "fileformat" / GNUPLOT, "filecycles" / 0, &
                "fileformat" / VTK, &
!                 "fileformat" / XDMF, &
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
    TYPE(marray_base) :: vcart,vcurv
    REAL              :: xmin,ymin,xmax,ymax,x0,y0
#ifdef PARALLEL
    REAL              :: xmin_all,xmax_all,ymin_all,ymax_all
    INTEGER           :: ierr
#endif
    REAL, DIMENSION(NUM_TESTS,4):: den,pre,vx,vy
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

    vcart = marray_base(3)
    vcurv = marray_base(3)

    ! set defaults
    den(:,:) = 1.
    vx(:,:)  = 0.
    vy(:,:)  = 0.
    pre(:,:) = 1.

    IF (ic.LT.1.OR.ic.GT.NUM_TESTS) &
      CALL Mesh%Error("InitData","Sorry, this 2D Riemann problem is currently not supported!")

    WRITE (teststr,'(A,I2)') "2D Riemann problem no. ", ic

    ! initial conditions

    ! Test configuration no. 1
    ! 1st quadrant
    den(1,1) = 1.
    pre(1,1) = 1.
    ! 2nd quadrant
    den(1,2) = .5197
    vx(1,2)  = -.7259
    pre(1,2) = .4
    ! 3rd quadrant
    den(1,3) = .1072
    vx(1,3)  = -.7259
    vy(1,3)  = -1.4045
    pre(1,3) = .0439
    ! 4th quadrant
    den(1,4) = .2579
    vy(1,4)  = -1.4045
    pre(1,4) = .15
    
    ! Test configuration no. 2
    ! 1st quadrant
    den(2,1) = 1.
    pre(2,1) = 1.
    ! 2nd quadrant
    den(2,2) = .5197
    vx(2,2) = -.7259
    pre(2,2) = .4
    ! 3rd quadrant
    den(2,3) = 1.
    vx(2,3) = -.7259
    vy(2,3) = -.7259
    pre(2,3) = 1.
    ! 4th quadrant
    den(2,4) = .5197
    vy(2,4) = -.7259
    pre(2,4) = .4

    ! Test configuration no. 3
    ! 1st quadrant
    den(3,1) = 1.5
    pre(3,1) = 1.5
    ! 2nd quadrant
    den(3,2) = .5323
    vx(3,2) = 1.206
    pre(3,2) = .3
    ! 3rd quadrant
    den(3,3) = 0.138
    vx(3,3) = 1.206
    vy(3,3) = 1.206
    pre(3,3) = 0.029
    ! 4th quadrant
    den(3,4) = .5323
    vy(3,4) = 1.206
    pre(3,4) = .3

    ! Test configuration no. 4
    ! 1st quadrant
    den(4,1) = 1.1
    pre(4,1) = 1.1
    ! 2nd quadrant
    den(4,2) = .5065
    vx(4,2) = .8939
    pre(4,2) = .35
    ! 3rd quadrant
    den(4,3) = 1.1
    vx(4,3) = .8939
    vy(4,3) = .8939
    pre(4,3) = 1.1
    ! 4th quadrant
    den(4,4) = .5065
    vy(4,4) = .8939
    pre(4,4) = .35

    ! 1st quadrant
    ! Test configuration no. 5
    den(5,1) = 1.
    vx(5,1) = -.75
    vy(5,1) = -.5
    pre(5,1) = 1.
    ! 2nd quadrant
    den(5,2) = 2.
    vx(5,2) = -.75
    vy(5,2) = .5
    pre(5,2) = 1.
    ! 3rd quadrant
    den(5,3) = 1.
    vx(5,3) = .75
    vy(5,3) = .5
    pre(5,3) = 1.
    ! 4th quadrant
    den(5,4) = 3.
    vx(5,4) = .75
    vy(5,4) = -.5
    pre(5,4) = 1.

    ! Test configuration no. 6
    ! 1st quadrant
    den(6,1) = 1.
    vx(6,1) = 0.75
    vy(6,1) = -0.5
    pre(6,1) = 1.
    ! 2nd quadrant
    den(6,2) = 2.
    vx(6,2) = 0.75
    vy(6,2) = 0.5
    pre(6,2) = 1.
    ! 3rd quadrant
    den(6,3) = 1.
    vx(6,3) = -0.75
    vy(6,3) = 0.5
    pre(6,3) = 1.
    ! 4th quadrant
    den(6,4) = 3.
    vx(6,4) = -0.75
    vy(6,4) = -0.5
    pre(6,4) = 1.

    ! Test configuration no. 7
    ! 1st quadrant
    den(7,1) = 1.
    vx(7,1) = 0.1
    vy(7,1) = 0.1
    pre(7,1) = 1.
    ! 2nd quadrant
    den(7,2) = 0.5197
    vx(7,2) = -0.6259
    vy(7,2) = 0.1
    pre(7,2) = 0.4
    ! 3rd quadrant
    den(7,3) = 0.8
    vx(7,3) = 0.1
    vy(7,3) = 0.1
    pre(7,3) = 0.4
    ! 4th quadrant
    den(7,4) = 0.5197
    vx(7,4) = 0.1
    vy(7,4) = -0.6259
    pre(7,4) = 0.4

    ! Test configuration no. 8
    ! 1st quadrant
    den(8,1) = 0.5197
    vx(8,1) = 0.1
    vy(8,1) = 0.1
    pre(8,1) = 0.4
    ! 2nd quadrant
    den(8,2) = 1.
    vx(8,2) = -0.6259
    vy(8,2) = 0.1
    pre(8,2) = 1.
    ! 3rd quadrant
    den(8,3) = 0.8
    vx(8,3) = 0.1
    vy(8,3) = 0.1
    pre(8,3) = 1.
    ! 4th quadrant
    den(8,4) = 1.
    vx(8,4) = 0.1
    vy(8,4) = -0.6259
    pre(8,4) = 1.

    ! Test configuration no. 9
    ! 1st quadrant
    den(9,1) = 1.
    vx(9,1) = 0.
    vy(9,1) = 0.3
    pre(9,1) = 1.
    ! 2nd quadrant
    den(9,2) = 2.
    vx(9,2) = 0.
    vy(9,2) = -0.3
    pre(9,2) = 1.
    ! 3rd quadrant
    den(9,3) = 1.039
    vx(9,3) = 0.
    vy(9,3) = -0.8133
    pre(9,3) = 0.4
    ! 4th quadrant
    den(9,4) = 0.5197
    vx(9,4) = 0.
    vy(9,4) = -0.4259
    pre(9,4) = 0.4

    ! Test configuration no. 10
    ! 1st quadrant
    den(10,1) = 1.
    vx(10,1) = 0.
    vy(10,1) = 0.4297
    pre(10,1) = 1.
    ! 2nd quadrant
    den(10,2) = 0.5
    vx(10,2) = 0.
    vy(10,2) = 0.6076
    pre(10,2) = 1.
    ! 3rd quadrant
    den(10,3) = 0.2281
    vx(10,3) = 0.
    vy(10,3) = -0.6076
    pre(10,3) = 0.3333
    ! 4th quadrant
    den(10,4) = 0.4562
    vx(10,4) = 0.
    vy(10,4) = -0.4297
    pre(10,4) = 0.3333

    ! Test configuration no. 11
    ! 1st quadrant
    den(11,1) = 1.
    vx(11,1) = 0.1
    vy(11,1) = 0.
    pre(11,1) = 1.
    ! 2nd quadrant
    den(11,2) = 0.5313
    vx(11,2) = 0.8276
    vy(11,2) = 0.
    pre(11,2) = 0.4
    ! 3rd quadrant
    den(11,3) = 0.8
    vx(11,3) = 0.1
    vy(11,3) = 0.
    pre(11,3) = 0.4
    ! 4th quadrant
    den(11,4) = 0.5313
    vx(11,4) = 0.1
    vy(11,4) = 0.7276
    pre(11,4) = 0.4

    ! Test configuration no. 12
    ! 1st quadrant
    den(12,1) = 0.5313
    vx(12,1) = 0.
    vy(12,1) = 0.
    pre(12,1) = 0.4
    ! 2nd quadrant
    den(12,2) = 1.
    vx(12,2) = 0.7276
    vy(12,2) = 0.
    pre(12,2) = 1.
    ! 3rd quadrant
    den(12,3) = 0.8
    vx(12,3) = 0.
    vy(12,3) = 0.
    pre(12,3) = 1.
    ! 4th quadrant
    den(12,4) = 1.
    vx(12,4) = 0.
    vy(12,4) = 0.7276
    pre(12,4) = 1.

    ! Test configuration no. 13
    ! 1st quadrant
    den(13,1) = 1.
    vx(13,1) = 0.
    vy(13,1) = -0.3
    pre(13,1) = 1.
    ! 2nd quadrant
    den(13,2) = 2.
    vx(13,2) = 0.
    vy(13,2) = 0.3
    pre(13,2) = 1.
    ! 3rd quadrant
    den(13,3) = 1.0625
    vx(13,3) = 0.
    vy(13,3) = 0.8145
    pre(13,3) = 0.4
    ! 4th quadrant
    den(13,4) = 0.5313
    vx(13,4) = 0.
    vy(13,4) = 0.4276
    pre(13,4) = 0.4

    ! Test configuration no. 14
    ! 1st quadrant
    den(14,1) = 2.
    vx(14,1) = 0.
    vy(14,1) = -0.5606
    pre(14,1) = 8.
    ! 2nd quadrant
    den(14,2) = 1.
    vx(14,2) = 0.
    vy(14,2) = -1.2172
    pre(14,2) = 8.
    ! 3rd quadrant
    den(14,3) = 0.4736
    vx(14,3) = 0.
    vy(14,3) = 1.2172
    pre(14,3) = 2.6667
    ! 4th quadrant
    den(14,4) = 0.9474
    vx(14,4) = 0.
    vy(14,4) = 1.1606
    pre(14,4) = 2.6667

    ! Test configuration no. 15
    ! 1st quadrant
    den(15,1) = 1.
    vx(15,1) = 0.1
    vy(15,1) = -0.3
    pre(15,1) = 1.
    ! 2nd quadrant
    den(15,2) = 0.5197
    vx(15,2) = -0.6259
    vy(15,2) = -0.3
    pre(15,2) = 0.4
    ! 3rd quadrant
    den(15,3) = 0.8
    vx(15,3) = 0.1
    vy(15,3) = -0.3
    pre(15,3) = 0.4
    ! 4th quadrant
    den(15,4) = 0.5313
    vx(15,4) = 0.1
    vy(15,4) = 0.4276
    pre(15,4) = 0.4

    ! Test configuration no. 16
    ! 1st quadrant
    den(16,1) = 0.5313
    vx(16,1) = 0.1
    vy(16,1) = 0.1
    pre(16,1) = 0.4
    ! 2nd quadrant
    den(16,2) = 1.02222
    vx(16,2) = -0.6179
    vy(16,2) = 0.1
    pre(16,2) = 1.
    ! 3rd quadrant
    den(16,3) = 0.8
    vx(16,3) = 0.1
    vy(16,3) = 0.1
    pre(16,3) = 1.
    ! 4th quadrant
    den(16,4) = 1.
    vx(16,4) = 0.1
    vy(16,4) = 0.8276
    pre(16,4) = 1.

    ! Test configuration no. 17
    ! 1st quadrant
    den(17,1) = 1.
    vx(17,1) = 0.
    vy(17,1) = -0.4
    pre(17,1) = 1.
    ! 2nd quadrant
    den(17,2) = 2.
    vx(17,2) = 0.
    vy(17,2) = -0.3
    pre(17,2) = 1.
    ! 3rd quadrant
    den(17,3) = 1.0625
    vx(17,3) = 0.
    vy(17,3) = 0.2145
    pre(17,3) = 0.4
    ! 4th quadrant
    den(17,4) = 0.5197
    vx(17,4) = 0.
    vy(17,4) = -1.1259
    pre(17,4) = 0.4

    ! Test configuration no. 18
    ! 1st quadrant
    den(18,1) = 1.
    vx(18,1) = 0.
    vy(18,1) = 1.
    pre(18,1) = 1.
    ! 2nd quadrant
    den(18,2) = 2.
    vx(18,2) = 0.
    vy(18,2) = -0.3
    pre(18,2) = 1.
    ! 3rd quadrant
    den(18,3) = 1.0625
    vx(18,3) = 0.
    vy(18,3) = 0.2145
    pre(18,3) = 0.4
    ! 4th quadrant
    den(18,4) = 0.5197
    vx(18,4) = 0.
    vy(18,4) = 0.2741
    pre(18,4) = 0.4

    ! Test configuration no. 19
    ! 1st quadrant
    den(19,1) = 1.
    vx(19,1) = 0.
    vy(19,1) = 0.3
    pre(19,1) = 1.
    ! 2nd quadrant
    den(19,2) = 2.
    vx(19,2) = 0.
    vy(19,2) = -0.3
    pre(19,2) = 1.
    ! 3rd quadrant
    den(19,3) = 1.0625
    vx(19,3) = 0.
    vy(19,3) = 0.2145
    pre(19,3) = 0.4
    ! 4th quadrant
    den(19,4) = 0.5197
    vx(19,4) = 0.
    vy(19,4) = -0.4259
    pre(19,4) = 0.4

    SELECT TYPE(pvar => Timedisc%pvar)
    TYPE IS(statevector_euler) ! non-isothermal HD
      WHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
        pvar%density%data3d(:,:,:)    = den(ic,1)
        pvar%pressure%data3d(:,:,:)   = pre(ic,1)
        vcart%data4d(:,:,:,1) = vx(ic,1)
        vcart%data4d(:,:,:,2) = vy(ic,1)
      ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).GE.y0) )
        pvar%density%data3d(:,:,:)    = den(ic,2)
        pvar%pressure%data3d(:,:,:)   = pre(ic,2)
        vcart%data4d(:,:,:,1) = vx(ic,2)
        vcart%data4d(:,:,:,2) = vy(ic,2)
      ELSEWHERE ( (Mesh%bccart(:,:,:,1).LT.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
        pvar%density%data3d(:,:,:)    = den(ic,3)
        pvar%pressure%data3d(:,:,:)   = pre(ic,3)
        vcart%data4d(:,:,:,1) = vx(ic,3)
        vcart%data4d(:,:,:,2) = vy(ic,3)
      ELSEWHERE ( (Mesh%bccart(:,:,:,1).GE.x0).AND.(Mesh%bccart(:,:,:,2).LT.y0) )
        pvar%density%data3d(:,:,:)    = den(ic,4)
        pvar%pressure%data3d(:,:,:)   = pre(ic,4)
        vcart%data4d(:,:,:,1) = vx(ic,4)
        vcart%data4d(:,:,:,2) = vy(ic,4)
      END WHERE
      vcart%data2d(:,3) = 0.0 ! z-velocity
      CALL Mesh%geometry%Convert2Curvilinear(Mesh%bcenter,vcart%data4d,vcurv%data4d)
      pvar%velocity%data2d(:,1:2) = vcurv%data2d(:,1:2)
    END SELECT
   
    CALL Physics%Convert2Conservative(Timedisc%pvar,Timedisc%cvar)
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
