r"""
Differential form modules

The set of differential forms along an open subset `U` of some manifold `S`
with values in a open subset `V` of a manifold `M` (possibly `S=M` and `U=V`)
is a module over the algebra `C^\infty(U)` of differentiable scalar fields
on `U`. It is a free module iff `V` is parallelizable.
Accordingly, two classes are devoted to differential form modules:

- :class:`DiffFormModule` for differential forms with values in a generic (in
  practice, not parallelizable) open set `V`
- :class:`DiffFormFreeModule` for differential forms with values in a
  parallelizable open set `V`

AUTHORS:

- Eric Gourgoulhon (2015): initial version

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.modules import Modules
from diffform import DiffForm, DiffFormParal

class DiffFormModule(UniqueRepresentation, Parent):
    r"""
    Module of differential forms of a given degree `p` along an open subset `U`
    of some manifold `S` with values in a open subset `V` of a manifold `M`.

    This is a module over `C^\infty(U)`, the ring (algebra) of differentiable
    scalar fields on `U`.

    The standard case of differential forms *on* a manifold corresponds to
    `U=V` (and hence `S=M`). Another common case is `\Phi` being an
    immersion.

    If `V` is parallelizable, the class :class:`DiffFormFreeModule` should
    be used instead.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the mapping `\Phi:\; U \rightarrow V`.
    - ``degree`` -- positive integer; the degree `p` of the differential forms

    """

    Element = DiffForm

    def __init__(self, vector_field_module, degree):
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        name = "/\^" + str(degree) + "(" + domain._name
        latex_name = r"\Lambda^{" + str(degree) + r"}\left(" + \
                     domain._latex_name
        if dest_map is domain._identity_map:
            name += ")"
            latex_name += r"\right)"
        else:
            name += "," + dest_map._name + ")"
            latex_name += "," + dest_map._latex_name + r"\right)"
        self._vmodule = vector_field_module
        self._degree = degree
        self._name = name
        self._latex_name = latex_name
        # the member self._ring is created for efficiency (to avoid calls to
        # self.base_ring()):
        self._ring = domain.scalar_field_algebra()
        Parent.__init__(self, base=self._ring, category=Modules(self._ring))
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain
        # NB: self._zero_element is not constructed here, since no element
        # can be constructed here, to avoid some infinite recursion.
