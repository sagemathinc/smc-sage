r"""
Tensors on free modules

The class :class:`FreeModuleTensor` implements tensors over a free module `M`,
i.e. elements of the free module `T^{(k,l)}(M)` of tensors of type (k,l) 
acting as multilinear forms on `M`. 

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

from sage.structure.element import ModuleElement  
#!# or from sage.structure.element import Element
# to avoid arithmetics defined in ModuleElement ??

class FreeModuleTensor(ModuleElement):
    r"""
    Tensor over a free module `M`.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` (must be an instance of 
      :class:`GenFreeModule`)
    - ``tensor_type`` -- pair (k,l) with k being the contravariant rank and l 
      the covariant rank
    - ``name`` -- (default: None) name given to the tensor
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor; 
      if none is provided, the LaTeX symbol is set to ``name``
    """
    def __init__(self, fmodule, tensor_type, name=None, latex_name=None):
        ModuleElement.__init__(self, fmodule.tensor_module(*tensor_type))
        self.fmodule = fmodule
        self.tensor_type = tuple(tensor_type)
        self.tensor_rank = self.tensor_type[0] + self.tensor_type[1]
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        self.components = {}    # components not set yet

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "tensor "
        if self.name is not None:
            description += self.name + " " 
        description += "of type (%s,%s)" % \
                           (str(self.tensor_type[0]), str(self.tensor_type[1]))
        description += " on the " + str(self.fmodule)
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

# From sage/modules/module.pyx:
#-----------------------------
### The Element should also implement _rmul_ (or _lmul_)
#
# class MyElement(sage.structure.element.ModuleElement):
#     def _rmul_(self, c):
#         ...


class FreeModuleVector(FreeModuleTensor):
    r"""
    Element (vector) of a free module `M`.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` (must be an instance of 
      :class:`GenFreeModule`)
    - ``name`` -- (default: None) name given to the vector
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the vector; 
      if none is provided, the LaTeX symbol is set to ``name``
    """
    def __init__(self, fmodule, name=None, latex_name=None):
        FreeModuleTensor.__init__(self, fmodule, (1,0), name, latex_name)

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "element "
        if self.name is not None:
            description += self.name + " " 
        description += "on the " + str(self.fmodule)
        return description

        

        
