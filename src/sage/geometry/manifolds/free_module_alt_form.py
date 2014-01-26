r"""
Alternating forms on free modules

The class :class:`FreeModuleAltForm` implement alternating forms on a free 
module over a commutative ring. 

It is a subclass of :class:`FreeModuleTensor`, alternating forms being a 
special type of tensors. 

A subclass of :class:`FreeModuleAltForm` is :class:`FreeModuleLinForm` for 
alternating forms of degree 1, i.e. linear forms. 

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

from free_module_tensor import FreeModuleTensor, FreeModuleVector
from comp import Components, CompFullyAntiSym

class FreeModuleAltForm(FreeModuleTensor):
    r"""
    Alternating form over a free module `M`.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` 
      (must be an instance of :class:`FiniteFreeModule`)
    - ``degree`` -- the degree of the alternating form (i.e. its tensor rank)
    - ``name`` -- (default: None) name given to the alternating form
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the alternating 
      form; if none is provided, the LaTeX symbol is set to ``name``

    """
    def __init__(self, fmodule, degree, name=None, latex_name=None):
        FreeModuleTensor.__init__(self, fmodule, (0,degree), name=name, 
                                  latex_name=latex_name, antisym=range(degree))
        FreeModuleAltForm._init_derived(self) # initialization of derived quantities

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "alternating form "
        if self.name is not None:
            description += self.name + " " 
        description += "of degree " + str(self.tensor_rank) + " on the " + \
                       str(self.fmodule)
        return description

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        FreeModuleTensor._init_derived(self)  

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        FreeModuleTensor._del_derived(self)  

    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 

        This method, which is already implemented in 
        :meth:`FreeModuleTensor._new_comp`, is redefined here for efficiency
        """
        fmodule = self.fmodule  # the base free module
        if self.tensor_rank == 1: 
            return Components(fmodule.ring, basis, 1,
                              start_index=fmodule.sindex,
                              output_formatter=fmodule.output_formatter)
        else:
            return CompFullyAntiSym(fmodule.ring, basis, self.tensor_rank, 
                                    start_index=fmodule.sindex,
                                    output_formatter=fmodule.output_formatter)

    def degree(self):
        r"""
        Return the degree of the alternating form. 
        """
        return self.tensor_rank

#******************************************************************************

class FreeModuleLinForm(FreeModuleAltForm):
    r"""
    Linear form over a free module `M`.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` 
      (must be an instance of :class:`FiniteFreeModule`)
    - ``name`` -- (default: None) name given to the linear form
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the linear 
      form; if none is provided, the LaTeX symbol is set to ``name``

    """
    def __init__(self, fmodule, name=None, latex_name=None):
        FreeModuleAltForm.__init__(self, fmodule, 1, name=name, 
                                   latex_name=latex_name)

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "linear form "
        if self.name is not None:
            description += self.name + " " 
        description += "on the " + str(self.fmodule)
        return description

    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 
              
        This method, which is already implemented in 
        :meth:`FreeModuleAltForm._new_comp`, is redefined here for efficiency
        """
        fmodule = self.fmodule  # the base free module
        return Components(fmodule.ring, basis, 1, start_index=fmodule.sindex,
                          output_formatter=fmodule.output_formatter)

    def __call__(self, vector):
        r"""
        The linear form acting on an element of the module.
        
        INPUT:
        
        - ``vector`` -- an element of the module (instance of 
          :class:`FreeModuleVector`)
        
        OUTPUT:
        
        - ring element `\langle \omega, v \rangle`
          
        """
        if not isinstance(vector, FreeModuleVector):
            raise TypeError("The argument must be a free module element.")
        basis = self.common_basis(vector)
        if basis is None:
            raise ValueError("No common basisfor the components.")
        omega = self.components[basis]
        vv = vector.components[basis]
        resu = 0
        for i in self.manifold.irange():
            resu += omega[[i]]*vv[[i]]
        # Name of the output:
        resu.name = None
        if self.name is not None and vector.name is not None:
            resu.name = self.name + "(" + vector.name + ")"
        # LaTeX symbol for the output:
        resu.latex_name = None
        if self.latex_name is not None and vector.latex_name is not None:
            resu.latex_name = self.latex_name + r"\left(" + \
                              vector.latex_name + r"\right)"
        return resu
















