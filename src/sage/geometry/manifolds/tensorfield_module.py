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

from sage.tensor.modules.finite_free_module import FiniteFreeModule
from sage.tensor.modules.tensor_free_module import TensorFreeModule
from scalarfield import ScalarField
from vectorfield import VectorFieldParal
from tensorfield import TensorFieldParal

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
    
    Element = VectorFieldParal

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
        self.domain = domain
        self.ambient_domain = ambient_domain
        manif = ambient_domain.manifold
        FiniteFreeModule.__init__(self, domain.scalar_field_ring(), 
                                  manif.dim, name=name, latex_name=latex_name, 
                                  start_index=manif.sindex,
                                  output_formatter=ScalarField.function_chart)

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

    def tensor_module(self, k, l):
        r"""
        Return the free module of all tensors of type (k,l) defined on 
        ``self``. 
        
        INPUT: 
        
        - ``k`` -- (non-negative integer) the contravariant rank, the tensor type 
          being (k,l)
        - ``l`` -- (non-negative integer) the covariant rank, the tensor type 
          being (k,l)
        
        OUTPUT:

        - instance of 
          :class:`TensorFieldFreeModule` 
          representing the free module of type-`(k,l)` tensors on the 
          free module ``self``. 
        
        EXAMPLES:
        
        
        """
        if (k,l) not in self._tensor_modules:
            self._tensor_modules[(k,l)] = TensorFieldFreeModule(self, (k,l))
        return self._tensor_modules[(k,l)]

    def linear_form(self, name=None, latex_name=None):
        r"""
        Construct a linear form on the free module. 
        
        A *linear form* on a free module `M` over a ring `R` is a map
        `M\rightarrow R` that is linear. It can be viewed as a tensor field
        of type (0,1) on `M`. 

        INPUT:
    
        - ``name`` -- (string; default: None) name given to the linear 
          form
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
          linear form; if none is provided, the LaTeX symbol is set to 
          ``name``
          
        OUTPUT:
        
        - instance of 
          :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleLinForm` 

        EXAMPLES:
        

        """
        from tensorfield import TensorFieldParal
        #!# provisory
        return TensorFieldParal(self, (0,1), name=name, latex_name=latex_name)


#******************************************************************************

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

