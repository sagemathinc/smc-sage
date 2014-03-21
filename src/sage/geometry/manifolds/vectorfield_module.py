r"""
Vector field module. 


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version
"""

#******************************************************************************
#       Copyright (C) 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.finite_free_module import FiniteFreeModule
from sage.geometry.manifolds.scalarfield import ScalarField

class VectorFieldFreeModule(FiniteFreeModule):
    r"""
    Module of vector fields along an open subset `U` of some immersed 
    submanifold `S` of a manifold `M` with values in a parallelizable
    open subset `V` of `M`. 
    
    Since `V` is parallelizable, the module is a free module over `C^\infty(U)`,
    the ring of differentiable scalar fields on `U`. 
    
    The standard case of vector fields *on* a manifold corresponds to 
    `U=V` (and hence `S=M`).

    INPUT:
    
    - ``domain`` -- open subset `U` on which the vector fields are defined
    - ``ambient_domain`` -- (default: None) parallelizable open subset `V` 
      of the ambient manifold `M`; if None, it is set to ``domain``.
    
    """
    
    # Element = ???

    def __init__(self, domain, ambient_domain=None):
        name = "X(" + domain.name
        latex_name = r"\mathcal{X}\left(" + domain.latex_name
        if ambient_domain is None or ambient_domain is domain:
            ambient_domain = domain
            name += ")" 
            latex_name += r"\right)" 
        else:
            name += "," + ambient_domain.name + ")" 
            latex_name += "," + ambient_domain.latex_name + r"\right)" 
        manif = ambient_domain.manifold
        FiniteFreeModule.__init__(self, domain.scalar_field_ring(), 
                                  manif.dim, name=name, latex_name=latex_name, 
                                  start_index=manif.sindex,
                                  output_formatter=ScalarField.function_chart)
        self.domain = domain
        self.ambient_domain = ambient_domain

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "free module "
        if self.name is not None:
            description += self.name + " "
        description += "of vector fields "
        if self.domain == self.ambient_domain:
            description += "on the " + str(self.domain)
        else:
            description += "along the " + str(self.domain) + " within the " + \
                           str(self.ambient_domain)
        return description
