r"""
Tensor field module. 


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

from sage.tensor.modules.tensor_free_module import TensorFreeModule
from tensorfield import TensorFieldParal

class TensorFieldFreeModule(TensorFreeModule):
    r"""
    Module of tensor fields of a given type `(k,l)` along an open subset `U` 
    of some immersed  submanifold `S` of a manifold `M` with values in 
    a parallelizable open subset `V` of `M`. 
    
    Since `V` is parallelizable, the module is a free module over `C^\infty(U)`,
    the ring of differentiable scalar fields on `U`. 
    
    The standard case of tensor fields *on* a manifold corresponds to 
    `U=V` (and hence `S=M`).

    INPUT:
    
    - ``vector_field_module`` -- free module `X(U,V)` of vector fields along
      `U` with values on `V`
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank and 
      `l` the covariant rank
    
    """

    
    Element = TensorFieldParal

    def __init__(self, vector_field_module, tensor_type):
        domain = vector_field_module.domain
        ambient_domain = vector_field_module.ambient_domain
        kcon = tensor_type[0]
        lcov = tensor_type[1]
        name = "TF^(" + str(kcon) + "," + str(lcov) + ")(" + domain.name
        latex_name = "TF^(" + str(kcon) + "," + str(lcov) + r")\left(" + \
                     domain.latex_name
        if ambient_domain is domain:
            ambient_domain = domain
            name += ")" 
            latex_name += r"\right)" 
        else:
            name += "," + ambient_domain.name + ")" 
            latex_name += "," + ambient_domain.latex_name + r"\right)" 
        TensorFreeModule.__init__(self, vector_field_module, tensor_type, 
                                  name=name, latex_name=latex_name)
        self.domain = domain
        self.ambient_domain = ambient_domain

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "free module "
        if self.name is not None:
            description += self.name + " "
        description += "of type-(%s,%s)" % \
                           (str(self.tensor_type[0]), str(self.tensor_type[1]))
        description += " tensors fields "
        if self.domain == self.ambient_domain:
            description += "on the " + str(self.domain)
        else:
            description += "along the " + str(self.domain) + " within the " + \
                           str(self.ambient_domain)
        return description

