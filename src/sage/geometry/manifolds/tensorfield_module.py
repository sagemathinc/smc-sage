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

from sage.modules.module import Module
from sage.structure.unique_representation import UniqueRepresentation
from sage.tensor.modules.tensor_free_module import TensorFreeModule
from tensorfield import TensorField, TensorFieldParal

class TensorFieldModule(UniqueRepresentation, Module):
    r"""
    Module of tensor fields of a given type `(k,l)` along an open subset `U` 
    of some manifold `S` with values in a open subset `V` of 
    a manifold `M`.
    
    This is a module over `C^\infty(U)`, the ring (algebra) of differentiable 
    scalar fields on `U`. 
    
    The standard case of tensor fields *on* a manifold corresponds to 
    `U=V` (and hence `S=M`). Another common case is `\Phi` being an 
    immersion.

    If `V` is parallelizable, the class :class:`TensorFieldFreeModule` should
    be used instead.
    
    INPUT:
    
    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector 
      fields along `U` associated with the mapping `\Phi:\; U \rightarrow V`. 
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank and 
      `l` the covariant rank
    
    """
    
    Element = TensorField

    def __init__(self, vector_field_module, tensor_type):
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        kcon = tensor_type[0]
        lcov = tensor_type[1]
        name = "TF^(" + str(kcon) + "," + str(lcov) + ")(" + domain._name
        latex_name = "TF^(" + str(kcon) + "," + str(lcov) + r")\left(" + \
                     domain._latex_name
        if dest_map is None:
            name += ")" 
            latex_name += r"\right)" 
        else:
            name += "," + dest_map._name + ")" 
            latex_name += "," + dest_map._latex_name + r"\right)" 
        self._vmodule = vector_field_module
        self._tensor_type = tensor_type
        self._name = name
        self._latex_name = latex_name
        # the member self._ring is created for efficiency (to avoid calls to 
        # self.base_ring()):
        self._ring = domain.scalar_field_algebra() 
        Module.__init__(self, self._ring)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain
        # NB: self._zero_element is not constructed here, since no element 
        # can be constructed here, to avoid some infinite recursion. 

    #### Methods required for any Parent 

    def _element_constructor_(self, comp=[], frame=None, name=None, 
                              latex_name=None, sym=None, antisym=None):
        r"""
        Construct a tensor field
        """
        if comp == 0:
            if not hasattr(self, '_zero_element'):
                self._zero_element = self._element_constructor_(name='zero', 
                                                                latex_name='0')
            return self._zero_element
        resu = self.element_class(self._vmodule, self._tensor_type, name=name, 
                                  latex_name=latex_name, sym=sym, 
                                  antisym=antisym)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) tensor field
        """
        resu = self.element_class(self._vmodule, self._tensor_type)
        return resu
            
    #### End of methods required for any Parent 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "module "
        if self._name is not None:
            description += self._name + " "
        description += "of type-(%s,%s)" % \
                           (str(self._tensor_type[0]), str(self._tensor_type[1]))
        description += " tensors fields "
        if self._dest_map is None:
            description += "on the " + str(self._domain)
        else:
            description += "along the " + str(self._domain) + \
                           " mapped into the " + str(self._ambient_domain)
        return description


#******************************************************************************

class TensorFieldFreeModule(TensorFreeModule):
    r"""
    Module of tensor fields of a given type `(k,l)` along an open subset `U` 
    of some manifold `S` with values in a parallelizable open subset `V` of 
    a manifold `M`.
    
    Since `V` is parallelizable, the module is a free module over `C^\infty(U)`,
    the ring (algebra) of differentiable scalar fields on `U`. 
    
    The standard case of tensor fields *on* a manifold corresponds to 
    `U=V` (and hence `S=M`). Another common case is `\Phi` being an 
    immersion.

    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector 
      fields along `U` associated with the mapping `\Phi:\; U \rightarrow V`. 
    - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant rank and 
      `l` the covariant rank
    
    """

    
    Element = TensorFieldParal

    def __init__(self, vector_field_module, tensor_type):
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        kcon = tensor_type[0]
        lcov = tensor_type[1]
        name = "TF^(" + str(kcon) + "," + str(lcov) + ")(" + domain._name
        latex_name = "TF^(" + str(kcon) + "," + str(lcov) + r")\left(" + \
                     domain._latex_name
        if dest_map is None:
            name += ")" 
            latex_name += r"\right)" 
        else:
            name += "," + dest_map._name + ")" 
            latex_name += "," + dest_map._latex_name + r"\right)" 
        TensorFreeModule.__init__(self, vector_field_module, tensor_type, 
                                  name=name, latex_name=latex_name)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "free module "
        if self._name is not None:
            description += self._name + " "
        description += "of type-(%s,%s)" % \
                           (str(self._tensor_type[0]), str(self._tensor_type[1]))
        description += " tensors fields "
        if self._dest_map is None:
            description += "on the " + str(self._domain)
        else:
            description += "along the " + str(self._domain) + \
                           " mapped into the " + str(self._ambient_domain)
        return description

