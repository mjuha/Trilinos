#@HEADER
# ************************************************************************
# 
#                     New_Package Example Package
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
#@HEADER
#np# This file is an example of a typical Makefile.am to build a test exe.  

EXEEXT = .exe

#np# List all executable names here
noinst_PROGRAMS = Object_test

#np# One SOURCES primary is needed for EACH exe.  (Notice for all of the following 
#np# primaries that the pattern is "(Exe_name)_(PRIMARY_NAME) = ").  List all source 
#np# files needed to compile the exe.
Object_test_SOURCES = $(srcdir)/cxx_main.cpp

#np# One DEPENDENCIES primary is needed for each exe.  List any libraries that the 
#np# exe depends on.  (This primary is not absolutely essential.)
Object_test_DEPENDENCIES = $(top_builddir)/src/libepetra.a

#np# One CXXFLAGS primary can be defined for each exe.  Include any -I's necessary for 
#np# building the exe.
Object_test_CXXFLAGS = -I$(top_srcdir)/src

#np# One LDADD primary can be defined for each exe.  Include any LDFLAGS necessary for 
#np# building the exe.
Object_test_LDADD = $(top_builddir)/src/libepetra.a \
	$(LIBS) $(BLAS_LIBS) $(FLIBS)


