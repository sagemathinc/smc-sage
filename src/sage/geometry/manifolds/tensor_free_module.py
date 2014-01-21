r"""
Tensor free modules

The class :class:`TensorFreeModule` implements the free module `T^{(k,l)}(M)` 
of tensors of type (k,l) acting as multilinear forms on a free module `M` 
over a commutative ring `R`. The free module `M` itself is canonically 
identified with `T^{(1,0)}(M)`. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

EXAMPLES:

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
from free_module_tensor import FreeModuleTensor, FreeModuleVector

# From sage/modules/module.pyx:
# ----------------------------
### The new Module class that should be the base of all Modules
### The derived Module class must implement the element
### constructor:
#
# class MyModule(sage.modules.module.Module):
#     Element = MyElement
#     def _element_constructor_(self, x):
#         return self.element_class(x)
#


class TensorFreeModule(Module):
    r"""
    Free module `T^{(k,l)}(M)` of tensors of type (k,l) acting as multilinear 
    forms on a free module `M` over a commutative ring `R`.
    
    .. NOTE::
    
        The class :class:`TensorFreeModule` inherits from :class:`Module` and
        not from :class:`FreeModule_generic` since the latter is a derived class
        of :class:`Module_old`, which does not conform to the new coercion model
        
    INPUT:
    
    - ``fmodule`` -- free module `M` (must be an instance of 
      :class:`GenFreeModule`)
    - ``tensor_type`` -- pair (k,l) with k being the contravariant rank and l 
      the covariant rank
    - ``name`` -- (default: None) name given to the free module
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the free module; 
      if none is provided, it is set to ``name``
    
    """
    
    Element = FreeModuleTensor
    
    def __init__(self, fmodule, tensor_type, name=None, latex_name=None):
        self.tensor_type = tuple(tensor_type)
        Module.__init__(self, fmodule.ring)
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        self.rank = pow(fmodule.rank, self.tensor_type[0] + self.tensor_type[1])
        self.fmodule = fmodule
        # Unique representation:
        if self.tensor_type in self.fmodule.tensor_modules:
            raise TypeError("The module of tensors of type" + 
                            str(self.tensor_type) + 
                            " has already been created.")
        else:
            self.fmodule.tensor_modules[self.tensor_type] = self
    
    #### Methods required for any Parent 
    def _element_constructor_(self, name=None, latex_name=None):
        r"""
        Construct a tensor
        """
        return self.element_class(self.fmodule, self.tensor_type, name, latex_name)

    def _an_element_(self):
        r"""
        Construct some (unamed) tensor
        """
        return self.element_class(self.fmodule, self.tensor_type)
            
    #### End of methods required for any Parent 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "free module "
        if self.name is not None:
            description += self.name + " "
        description += "of type (%s,%s)" % \
                           (str(self.tensor_type[0]), str(self.tensor_type[1]))
        description += " tensors on the " + str(self.fmodule)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        if self.latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self.latex_name


#******************************************************************************

class GenFreeModule(TensorFreeModule):
    r"""
    Free module over a commutative ring `R`.
    
    INPUT:
    
    - ``ring`` -- commutative ring `R` over which the free module is 
      constructed.
    - ``rank`` -- rank of the free module
    - ``name`` -- (default: None) name given to the free module
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the free module; 
      if none is provided, it is set to ``name``
    - ``start_index`` -- (default: 0) lower bound of the range of indices in 
      basis defined on the free module

    """
    
    Element = FreeModuleVector
    
    def __init__(self, ring, rank, name=None, latex_name=None, start_index=0):
        self.ring = ring
        self.rank = rank 
        self.tensor_modules = {} # Dictionary of the tensor modules built on 
                                 # self (dict. keys = (k,l) --the tensor type)
        TensorFreeModule.__init__(self, self, (1,0), name, latex_name)
        self.sindex = start_index
        self.known_bases = []  # List of known bases on the free module
        self.def_basis = None # default basis
        self.basis_changes = {} # Dictionary of the changes of bases

    #### Methods required for any Parent 
    def _element_constructor_(self, name=None, latex_name=None):
        r"""
        Construct a element of the module
        """
        return self.element_class(self, name, latex_name)

    def _an_element_(self):
        r"""
        Construct some (unamed) element of the module
        """
        return self.element_class(self)
            
    #### End of methods required for any Parent 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "free module "
        if self.name is not None:
            description += self.name + " "
        description += "over the " + str(self.ring)
        return description
        
    def tensor_module(self, k, l):
        r"""
        Return the free module of all tensors of type (k,l) defined on 
        ``self``. 
        
        INPUT: 
        
        - ``k`` -- the contravariant rank, the tensor type being (k,l)
        - ``l`` -- the covariant rank, the tensor type being (k,l)
        
        OUTPUT:
        
        - the free module `T^{(k,l)}(M)` of tensors of type (k,l) acting as 
          multilinear forms on the free module ``self``. 
          
        """
        if (k,l) not in self.tensor_modules:
            self.tensor_modules[(k,l)] = TensorFreeModule(self, (k,l))
        return self.tensor_modules[(k,l)]

    def irange(self, start=None):
        r"""
        Single index generator, labelling the elements of a basis.
                
        INPUT:
        
        - ``start`` -- (default: None) initial value of the index; if none is 
          provided, ``self.sindex`` is assumed

        OUTPUT:
        
        - an iterable index, starting from ``start`` and ending at
          ``self.sindex + self.dim -1``

        EXAMPLES:
               
        """
        si = self.sindex
        imax = self.rank + si
        if start is None:
            i = si
        else:
            i = start
        while i < imax:
            yield i
            i += 1















