// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Ioss_CodeTypes.h"           // for IntVector
#include "Ioss_ElementTopology.h"     // for ElementTopology
#include <Ioss_ElementVariableType.h> // for ElementVariableType
#include <Ioss_TriShell4.h>
#include <cassert> // for assert

//------------------------------------------------------------------------
// Define a variable type for storage of this elements connectivity
namespace Ioss {
  class St_TriShell4 : public ElementVariableType
  {
  public:
    static void factory() { static St_TriShell4 registerThis; }

  protected:
    St_TriShell4() : ElementVariableType("trishell4", 4) {}
  };
} // namespace Ioss
// ========================================================================
namespace {
  struct Constants
  {
    static const int nnode     = 4;
    static const int nedge     = 3;
    static const int nedgenode = 2;
    static const int nface     = 2;
    static const int nfacenode = 4;
    static const int nfaceedge = 3;
    static int       edge_node_order[nedge][nedgenode];
    static int       face_node_order[nface][nfacenode];
    static int       face_edge_order[nface][nfaceedge];
    static int       nodes_per_face[nface + 1];
    static int       edges_per_face[nface + 1];
  };
} // namespace

// Edge numbers are zero-based [0..number_edges)
int Constants::edge_node_order[nedge][nedgenode] = // [edge][edge_node]
    {{0, 1}, {1, 2}, {2, 0}};

// Face numbers are zero-based [0..number_faces)
int Constants::face_node_order[nface][nfacenode] = // [face][face_node]
    {{0, 1, 2, 3}, {0, 2, 1, 3}};

int Constants::face_edge_order[nface][nfaceedge] = // [face][face_edge]
    {{0, 1, 2}, {2, 1, 0}};

// face 0 returns number of nodes for all faces if homogenous
//        returns -1 if faces have differing topology
int Constants::nodes_per_face[nface + 1] = {4, 4, 4};

// face 0 returns number of edges for all faces if homogenous
//        returns -1 if faces have differing topology
int Constants::edges_per_face[nface + 1] = {3, 3, 3};

void Ioss::TriShell4::factory()
{
  static Ioss::TriShell4 registerThis;
  Ioss::St_TriShell4::factory();
}

Ioss::TriShell4::TriShell4() : Ioss::ElementTopology("trishell4", "ShellTriangle_4")
{
  Ioss::ElementTopology::alias("trishell4", "Shell_Tri_4_3D");
  Ioss::ElementTopology::alias("trishell4", "SHELL_TRIANGLE_4");
}

Ioss::TriShell4::~TriShell4() = default;

int Ioss::TriShell4::parametric_dimension() const { return 2; }
int Ioss::TriShell4::spatial_dimension() const { return 3; }
int Ioss::TriShell4::order() const { return 2; }

int Ioss::TriShell4::number_corner_nodes() const { return 3; }
int Ioss::TriShell4::number_nodes() const { return Constants::nnode; }
int Ioss::TriShell4::number_edges() const { return Constants::nedge; }
int Ioss::TriShell4::number_faces() const { return Constants::nface; }

int Ioss::TriShell4::number_nodes_edge(int /* edge */) const { return Constants::nedgenode; }

int Ioss::TriShell4::number_nodes_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::nodes_per_face[face];
}

int Ioss::TriShell4::number_edges_face(int face) const
{
  // face is 1-based.  0 passed in for all faces.
  assert(face >= 0 && face <= number_faces());
  return Constants::edges_per_face[face];
}

Ioss::IntVector Ioss::TriShell4::edge_connectivity(int edge_number) const
{
  assert(edge_number > 0 && edge_number <= Constants::nedge);
  Ioss::IntVector connectivity(Constants::nedgenode);

  for (int i = 0; i < Constants::nedgenode; i++) {
    connectivity[i] = Constants::edge_node_order[edge_number - 1][i];
  }

  return connectivity;
}

Ioss::IntVector Ioss::TriShell4::face_connectivity(int face_number) const
{
  assert(face_number > 0 && face_number <= number_faces());
  Ioss::IntVector connectivity(Constants::nodes_per_face[face_number]);

  for (int i = 0; i < Constants::nodes_per_face[face_number]; i++) {
    connectivity[i] = Constants::face_node_order[face_number - 1][i];
  }

  return connectivity;
}

Ioss::IntVector Ioss::TriShell4::element_connectivity() const
{
  Ioss::IntVector connectivity(number_nodes());
  for (int i = 0; i < number_nodes(); i++) {
    connectivity[i] = i;
  }
  return connectivity;
}

Ioss::ElementTopology *Ioss::TriShell4::face_type(int face_number) const
{
  assert(face_number >= 0 && face_number <= number_faces());
  //  return Ioss::ElementTopology::factory("triface4");
  return Ioss::ElementTopology::factory("tri4");
}

Ioss::ElementTopology *Ioss::TriShell4::edge_type(int edge_number) const
{
  assert(edge_number >= 0 && edge_number <= number_edges());
  return Ioss::ElementTopology::factory("edge2");
}

Ioss::IntVector Ioss::TriShell4::face_edge_connectivity(int face_number) const
{
  assert(face_number > 0 && face_number <= Constants::nface);

  int             nface_edge = number_edges_face(face_number);
  Ioss::IntVector fcon(nface_edge);

  for (int i = 0; i < nface_edge; i++) {
    fcon[i] = Constants::face_edge_order[face_number - 1][i];
  }

  return fcon;
}
