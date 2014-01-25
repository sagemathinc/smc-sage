r"""
Tensor free modules

The class :class:`FiniteFreeModule` implements free modules `M` of finite rank
over a commutative ring `R`. It inherits (via :class:`TensorFreeModule`) from 
the generic class :class:`sage.modules.module.Module`. 

The class :class:`TensorFreeModule` implements the free modules `T^{(k,l)}(M)` 
of tensors of type `(k,l)` acting as multilinear forms on a free module `M` 
of finite rank over a commutative ring `R`. The module `T^{(k,l)}(M)` is
canonically identified with the following tensor product:

.. MATH::

    T^{(k,l)}(M) = \underbrace{M\otimes\cdots\otimes M}_{k\ \; \mbox{times}}
    \otimes \underbrace{M^*\otimes\cdots\otimes M^*}_{l\ \; \mbox{times}}
    
where `M^*` stands for the dual of the free module `M`. 
In particular, the free module `M` itself is canonically identified 
with `T^{(1,0)}(M)`. Accordingly, the class :class:`FiniteFreeModule` 
inherits from  :class:`TensorFreeModule`. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

EXAMPLES:

    Free module of rank 3 over `\ZZ`::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M') ; M
        rank-3 free module M over the Integer Ring
        sage: M.category()
        Category of modules over Integer Ring
        sage: M.rank()
        3

    Set of all tensors of type (1,2) defined on the module M::
    
        sage: T = M.tensor_module(1,2) ; T
        free module of type (1,2) tensors on the rank-3 free module M over the Integer Ring
        sage: T.category()
        Category of modules over Integer Ring
        sage: T.rank()
        27

    The module `M` itself is considered as the set of tensors of type (1,0)::
    
        sage: M is M.tensor_module(1,0)
        True

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
    Set `T^{(k,l)}(M)` of tensors of type `(k,l)` acting as multilinear 
    forms on a free module `M` over a commutative ring `R`. 
    `T^{(k,l)}(M)` is a free module over `R`, which can be canonically 
    identified with the following tensor product of modules:

    .. MATH::
    
        T^{(k,l)}(M) = \underbrace{M\otimes\cdots\otimes M}_{k\ \; \mbox{times}}
        \otimes \underbrace{M^*\otimes\cdots\otimes M^*}_{l\ \; \mbox{times}}
        
    where `M^*` stands for the dual of the free module `M`. 
    
    .. NOTE::
    
        The class :class:`TensorFreeModule` inherits from the generic class
        :class:`sage.modules.module.Module` and not from 
        :class:`sage.modules.free_module.FreeModule_generic` since the latter 
        is a derived class of :class:`sage.modules.module.Module_old`, which 
        does not conform to the new coercion model.
        Besides, the class :class:`TensorFreeModule` does not inherit 
        from the class :class:`CombinatorialFreeModule` since the latter is
        devoted to modules *with a basis*. 

    :class:`TensorFreeModule` is a Sage *Parent* class whose elements belong 
    to the class :class:`FreeModuleTensor`. 

    INPUT:
    
    - ``fmodule`` -- free module `M` of finite rank (must be an instance of 
      :class:`FiniteFreeModule`)
    - ``tensor_type`` -- pair (k,l) with k being the contravariant rank and l 
      the covariant rank
    - ``name`` -- (string; default: None) name given to the tensor module
    - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
      tensor module; if none is provided, it is set to ``name``
    
    EXAMPLES:
    
    Set of tensors of type (1,2) on a free module of rank 3 over `\ZZ`::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: T = TensorFreeModule(M, (1,2)) ; T
        free module of type (1,2) tensors on the rank-3 free module M over the Integer Ring
        
    T is a module (actually a free module) over `\ZZ`::
    
        sage: T.category()
        Category of modules over Integer Ring
        sage: T in Modules(ZZ)
        True
        sage: T.rank()
        27

    T is a *parent* object, whose elements are instances of 
    :class:`FreeModuleTensor`::

        sage: t = T.an_element() ; t
        tensor of type (1,2) on the rank-3 free module M over the Integer Ring
        sage: from sage.geometry.manifolds.free_module_tensor import FreeModuleTensor
        sage: isinstance(t, FreeModuleTensor)
        True
        sage: t in T
        True
        sage: T.is_parent_of(t)
        True
  
    The tensor modules over a given module `M` are unique::
    
        sage: T is M.tensor_module(1,2)
        True

    Instead of a direct call to the constructor of class 
    :class:`TensorFreeModule`, one should construct them by the method
    :meth:`FiniteFreeModule.tensor_module` of the base module::
    
        sage: T = M.tensor_module(1,2) ; T
        free module of type (1,2) tensors on the rank-3 free module M over the Integer Ring
        sage: F = M.tensor_module(0,1) ; F
        free module of type (0,1) tensors on the rank-3 free module M over the Integer Ring
        
    The base module itself is considered as the set of type (1,0) tensors::
    
        sage: M is M.tensor_module(1,0)
        True
            
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
        self._rank = pow(fmodule._rank, self.tensor_type[0] + self.tensor_type[1])
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
        return self.element_class(self.fmodule, self.tensor_type, name, 
                                  latex_name)

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

    def rank(self):
        r"""
        Return the rank of the free module ``self``.
        """
        return self._rank
        

#******************************************************************************

class FiniteFreeModule(TensorFreeModule):
    r"""
    Free module of finite rank over a commutative ring `R`.

    This class inherits (via :class:`TensorFreeModule`) from the generic 
    class :class:`sage.modules.module.Module`. 
    
    .. NOTE::
    
        The class :class:`FiniteFreeModule` does not inherit from 
        :class:`sage.modules.free_module.FreeModule_generic` since the latter 
        is a derived class of :class:`sage.modules.module.Module_old`, 
        which does not conform to the new coercion model.
        Besides, the class :class:`FiniteFreeModule` does not inherit 
        from the class :class:`CombinatorialFreeModule` since the latter is
        devoted to modules *with a basis*. 

    The class :class:`FiniteFreeModule` is a Sage *Parent* class whose elements 
    belong to the class :class:`FreeModuleVector`. 

    INPUT:
    
    - ``ring`` -- commutative ring `R` over which the free module is 
      constructed.
    - ``rank`` -- (positive integer) rank of the free module
    - ``name`` -- (string; default: None) name given to the free module
    - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the free 
      module; if none is provided, it is set to ``name``
    - ``start_index`` -- (integer; default: 0) lower bound of the range of 
      indices in bases defined on the free module

    EXAMPLES:
    
    Free module of rank 3 over `\ZZ`::
    
        sage: M = FiniteFreeModule(ZZ, 3) ; M
        rank-3 free module over the Integer Ring
        sage: M = FiniteFreeModule(ZZ, 3, name='M') ; M  # declaration with a name
        rank-3 free module M over the Integer Ring
        sage: M.category()
        Category of modules over Integer Ring
        sage: M.base_ring()
        Integer Ring
        sage: M.rank()
        3
        
    The LaTeX output is adjusted via the parameter ``latex_name``::
    
        sage: latex(M)  # the default is the symbol provided in the string ``name``
        M
        sage: M = FiniteFreeModule(ZZ, 3, name='M', latex_name=r'\mathcal{M}')
        sage: latex(M)
        \mathcal{M}

    M is a *parent* object, whose elements are instances of 
    :class:`FreeModuleVector`::
    
        sage: v = M.an_element() ; v
        element of the rank-3 free module M over the Integer Ring
        sage: from sage.geometry.manifolds.free_module_tensor import FreeModuleVector
        sage: isinstance(v, FreeModuleVector)
        True        
        sage: v in M
        True
        sage: M.is_parent_of(v)
        True
        
    M is in the category of modules over `\ZZ`::

        sage: M.category()
        Category of modules over Integer Ring
        sage: M in Modules(ZZ)
        True

    Indices on the free module, such as indices labelling the element of a
    basis, are provided by the generator method :meth:`irange`. By default, 
    they range from 0 to the module's rank minus one::
    
        sage: for i in M.irange(): print i,
        0 1 2

    This can be changed via the parameter ``start_index`` in the module 
    construction::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M', start_index=1)
        sage: for i in M.irange(): print i,
        1 2 3

    """
    
    Element = FreeModuleVector
    
    def __init__(self, ring, rank, name=None, latex_name=None, start_index=0):
        self.ring = ring
        self._rank = rank 
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
        description = "rank-" + str(self._rank) + " free module "
        if self.name is not None:
            description += self.name + " "
        description += "over the " + str(self.ring)
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
        
        - ``start`` -- (integer; default: None) initial value of the index; if none is 
          provided, ``self.sindex`` is assumed

        OUTPUT:
        
        - an iterable index, starting from ``start`` and ending at
          ``self.sindex + self.rank() -1``

        EXAMPLES:
               
        """
        si = self.sindex
        imax = self._rank + si
        if start is None:
            i = si
        else:
            i = start
        while i < imax:
            yield i
            i += 1

