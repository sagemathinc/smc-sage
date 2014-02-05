"""
Specific rank-2 tensors on free modules

Four derived classes of :class:`FreeModuleTensor` are devoted to rank-2 
tensors::

* :class:`FreeModuleEndomorphism` for endomorphisms (type-(1,1) tensors)
* :class:`FreeModuleAutomorphism` for invertible endomorphisms
* :class:`FreeModuleIdentityMap` for the identity map on a free module
* :class:`FreeModuleSymBilinForm` for symmetric bilinear forms (symmetric 
  type-(0,2) tensor 
  
Antisymmetric type-(0,2) tensors are dealt with by the class
:class:`FreeModuleAltForm`.

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

from free_module_tensor import FreeModuleTensor

class FreeModuleEndomorphism(FreeModuleTensor):
    r"""
    Endomorphism (considered as a type-(1,1) tensor) on a free module.

    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` 
      (must be an instance of :class:`FiniteFreeModule`)
    - ``name`` -- (default: None) name given to the endomorphism
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
      endomorphism; if none is provided, the LaTeX symbol is set to ``name``
      
    """
    def __init__(self, fmodule, name=None, latex_name=None):
        FreeModuleTensor.__init__(self, fmodule, (1,1), name=name, 
                                  latex_name=latex_name)
                                  
    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "endomorphism "
        if self.name is not None:
            description += self.name + " " 
        description += " on the " + str(self.fmodule)
        return description

    def _new_instance(self):
        r"""
        Create a :class:`FreeModuleEndomorphism` instance. 
        
        """
        return FreeModuleEndomorphism(self.fmodule)

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`FreeModuleTensor.__call__` to allow for a single 
        vector argument. 
        """
        if len(arg) > 1:
            # the endomorphism acting as a type (1,1) tensor on a pair 
            # (linear form, vector), returning a scalar:
            return FreeModuleTensor.__call__(self, *arg) 
        # the endomorphism acting as such, on a vector, returning a vector:
        vector = arg[0]
        if not isinstance(vector, FreeModuleVector):
            raise TypeError("The argument must be an element of a free module.")
        basis = self.common_basis(vector)
        t = self.components[basis]
        v = vector.components[basis]
        fmodule = self.fmodule
        result = FreeModuleVector(fmodule)
        for i in fmodule.irange():
            res = 0
            for j in fmodule.irange():
                res += t[[i,j]]*v[[j]]
            result.set_comp(basis)[i] = res
        # Name of the output:
        result.name = None
        if self.name is not None and vector.name is not None:
            result.name = self.name + "(" + vector.name + ")"
        # LaTeX symbol for the output:
        result.latex_name = None
        if self.latex_name is not None and vector.latex_name is not None:
            result.latex_name = self.latex_name + r"\left(" + \
                              vector.latex_name + r"\right)"
        return result

