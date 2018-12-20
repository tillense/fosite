!#############################################################################
!#                                                                           #
!# fosite - 2D hydrodynamical simulation program                             #
!# header: tap.h                                                             #
!#                                                                           #
!# Copyright (C) 2013                                                        #
!# Manuel Jung <mjunf@astrophysik.uni-kiel.de>                               #
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

!----------------------------------------------------------------------------!
! TAP: Test Anything Protocoll Header
! For a short definition visit
! http://podwiki.hexten.net/TAP/TAP13.html?page=TAP13
! The interface is loosly oriented at the c library libtap
! https://github.com/zorgnax/libtap
! and the boost testing framework
! http://www.boost.org/doc/libs/1_53_0/libs/test/doc/html/utf/testing-tools/reference.html
!----------------------------------------------------------------------------!

USE tap

#define TAP_PLAN(tests) CALL tap_plan(tests)
#define TAP_DONE CALL tap_done()
#define TAP_CHECK(expr,name) CALL tap_check_at_loc(__FILE__, __LINE__, expr, name)
#define TAP_CHECK_CLOSE(a,b,eps,name) CALL tap_check_close_at_loc(__FILE__, __LINE__, a, b, eps, name)
#define TAP_CHECK_SMALL(a,eps,name) CALL tap_check_small_at_loc(__FILE__, __LINE__, a, eps, name)
#define TAP_CHECK_EQV(a,b,name) CALL tap_check_op_at_loc(__FILE__, __LINE__, ".EQV.", a.EQV.b, a, b, name)
#define TAP_CHECK_EQ(a,b,name) CALL tap_check_op_at_loc(__FILE__, __LINE__, ".EQ.", a.EQ.b, a, b, name)
#define TAP_CHECK_NE(a,b,name) CALL tap_check_op_at_loc(__FILE__, __LINE__, ".NE.", a.NE.b, a, b, name)
#define TAP_CHECK_GT(a,b,name) CALL tap_check_op_at_loc(__FILE__, __LINE__, ".GT.", a.GT.b, a, b, name)
#define TAP_CHECK_GE(a,b,name) CALL tap_check_op_at_loc(__FILE__, __LINE__, ".GE.", a.GE.b, a, b, name)
#define TAP_CHECK_LT(a,b,name) CALL tap_check_op_at_loc(__FILE__, __LINE__, ".LT.", a.LT.b, a, b, name)
#define TAP_CHECK_LE(a,b,name) CALL tap_check_op_at_loc(__FILE__, __LINE__, ".LE.", a.LE.b, a, b, name)

