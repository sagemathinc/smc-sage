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

from comp import Components, CompWithSym, CompFullySym, CompFullyAntiSym


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
    - ``sym`` -- (default: None) a symmetry or a list of symmetries among the 
      tensor arguments: each symmetry is described by a tuple containing 
      the positions of the involved arguments, with the convention position=0
      for the first argument. For instance:
        * sym=(0,1) for a symmetry between the 1st and 2nd arguments 
        * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
          arguments and a symmetry between the 2nd, 4th and 5th arguments.
    - ``antisym`` -- (default: None) antisymmetry or list of antisymmetries 
      among the arguments, with the same convention as for ``sym``. 

    """
    def __init__(self, fmodule, tensor_type, name=None, latex_name=None,
                 sym=None, antisym=None):
        ModuleElement.__init__(self, fmodule.tensor_module(*tensor_type))
        self.fmodule = fmodule
        self.tensor_type = tuple(tensor_type)
        self.tensor_rank = self.tensor_type[0] + self.tensor_type[1]
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        self.components = {}    # components on various bases (not set yet)
        # Treatment of symmetry declarations:
        self.sym = []
        if sym is not None and sym != []:
            if isinstance(sym[0], (int, Integer)):  
                # a single symmetry is provided as a tuple -> 1-item list:
                sym = [tuple(sym)]
            for isym in sym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self.tensor_rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                 " not in [0," + str(self.tensor_rank-1) + "]")
                    self.sym.append(tuple(isym))       
        self.antisym = []
        if antisym is not None and antisym != []:
            if isinstance(antisym[0], (int, Integer)):  
                # a single antisymmetry is provided as a tuple -> 1-item list:
                antisym = [tuple(antisym)]
            for isym in antisym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self.tensor_rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                " not in [0," + str(self.tensor_rank-1) + "]")
                    self.antisym.append(tuple(isym))
        # Final consistency check:
        index_list = []
        for isym in self.sym:
            index_list += isym
        for isym in self.antisym:
            index_list += isym
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("Incompatible lists of symmetries: the same " + 
                             "position appears more than once.")
        # Initialization of derived quantities:
        FreeModuleTensor._init_derived(self) 


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

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        pass # no derived quantities

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        pass # no derived quantities

    def symmetries(self):
        r"""
        Print the list of symmetries and antisymmetries.
        
        EXAMPLES:
        
        Various symmetries / antisymmetries for a rank-4 tensor::
        
            sage: from sage.geometry.manifolds.tensor_free_module import GenFreeModule
            sage: from sage.geometry.manifolds.free_module_tensor import FreeModuleTensor
            sage: M = GenFreeModule(ZZ, 3, name='M')
            sage: t = FreeModuleTensor(M, (4,0), 'T') # no symmetry declared
            sage: t.symmetries()
            no symmetry;  no antisymmetry
            sage: t = FreeModuleTensor(M, (4,0), 'T', sym=(0,1))
            sage: t.symmetries()
            symmetry: (0, 1);  no antisymmetry
            sage: t = FreeModuleTensor(M, (4,0), 'T', sym=[(0,1), (2,3)])
            sage: t.symmetries()
            symmetries: [(0, 1), (2, 3)];  no antisymmetry
            sage: t = FreeModuleTensor(M, (4,0), 'T', sym=(0,1), antisym=(2,3))
            sage: t.symmetries()
            symmetry: (0, 1);  antisymmetry: (2, 3)
            
        """
        if len(self.sym) == 0:
            s = "no symmetry; "
        elif len(self.sym) == 1:
            s = "symmetry: " + str(self.sym[0]) + "; "
        else:
            s = "symmetries: " + str(self.sym) + "; " 
        if len(self.antisym) == 0:
            a = "no antisymmetry"
        elif len(self.antisym) == 1:
            a = "antisymmetry: " + str(self.antisym[0])
        else:
            a = "antisymmetries: " + str(self.antisym)   
        print s, a
         
    def set_name(self, name, latex_name=None):
        r"""
        Set (or change) the text name and LaTeX name of the tensor.

        INPUT:
        
        - ``name`` -- name given to the tensor
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor; 
          if none is provided, the LaTeX symbol is set to ``name``

        """
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
       
    
    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 
        
        This method, to be called by :meth:`comp`, must be redefined by derived 
        classes to adapt the output to the relevant subclass of 
        :class:`Components`.
        
        """
        if self.sym == [] and self.antisym == []:
            return Components(self.fmodule.ring, basis, self.tensor_rank,
                              start_index=self.fmodule.sindex)
        for isym in self.sym:
            if len(isym) == self.tensor_rank:
                return CompFullySym(self.fmodule.ring, basis, self.tensor_rank,
                                    start_index=self.fmodule.sindex)
        for isym in self.antisym:
            if len(isym) == self.tensor_rank:
                return CompFullyAntiSym(self.fmodule.ring, basis, 
                                        self.tensor_rank, 
                                        start_index=self.fmodule.sindex)
        return CompWithSym(self.fmodule.ring, basis, self.tensor_rank, 
                           start_index=self.fmodule.sindex, sym=self.sym,
                           antisym=self.antisym)        


    def comp(self, basis=None, from_basis=None):
        r"""
        Return the components in a given basis.
        
        If the components are not known already, they are computed by the tensor
        change-of-basis formula from components in another basis. 
        
        INPUT:
        
        - ``basis`` -- (default: None)  vector basis in which the components 
          are required; if none is provided, the components are assumed to 
          refer to the domain's default basis
        - ``from_basis`` -- (default: None) basis from which the
          required components are computed, via the tensor change-of-basis 
          formula, if they are not known already in the basis ``basis``; 
          if none, a basis is picked in ``self.components``.
 
        OUTPUT: 
        
        - components in the basis ``basis``, as an instance of the 
          class :class:`Components` 
        
        EXAMPLES:
        
        """
        fmodule = self.fmodule
        if basis is None: 
            basis = fmodule.def_basis
        if basis not in self.components:
            # The components must be computed from 
            # those in the basis from_basis
            if from_basis is None: 
                for known_basis in self.components:
                    if (known_basis, basis) in self.fmodule.basis_changes \
                        and (basis, known_basis) in self.fmodule.basis_changes:
                        from_basis = known_basis
                        break
                if from_basis is None:
                    raise ValueError("No basis could be found for computing " + 
                                     "the components in the " + str(basis))
            elif from_basis not in self.components:
                raise ValueError("The tensor components are not known in the " +
                                 "basis "+ str(from_basis))
            (n_con, n_cov) = self.tensor_type
            if n_cov > 0:
                if (from_basis, basis) not in fmodule.basis_changes:
                    raise ValueError("The change-of-basis matrix from the " + 
                                     str(from_basis) + " to the " + str(basis) 
                                     + " has not been set.")
                pp = \
                  fmodule.basis_changes[(from_basis, basis)].comp(from_basis)
                # pp not used if n_cov = 0 (pure contravariant tensor)
            if n_con > 0:
                if (basis, from_basis) not in fmodule.basis_changes:
                    raise ValueError("The change-of-basis matrix from the " + 
                                     str(basis) + " to the " + str(from_basis) +
                                     " has not been set.")
                ppinv = \
                  fmodule.basis_changes[(basis, from_basis)].comp(from_basis)
                # ppinv not used if n_con = 0 (pure covariant tensor)
            old_comp = self.components[from_basis]
            new_comp = self._new_comp(basis)
            rank = self.tensor_rank
            # loop on the new components:
            for ind_new in new_comp.non_redundant_index_generator(): 
                # Summation on the old components multiplied by the proper 
                # change-of-basis matrix elements (tensor formula): 
                res = 0 
                for ind_old in fmodule.index_generator(rank): 
                    t = old_comp[[ind_old]]
                    for i in range(n_con): # loop on contravariant indices
                        t *= ppinv[[ind_new[i], ind_old[i]]]
                    for i in range(n_con,rank):  # loop on covariant indices
                        t *= pp[[ind_old[i], ind_new[i]]]
                    res += t
                new_comp[ind_new] = res
            self.components[basis] = new_comp
            # end of case where the computation was necessary
        return self.components[basis]

    def set_comp(self, basis=None):
        r"""
        Return the components in a given basis for assignment.
        
        The components with respect to other bases are deleted, in order to 
        avoid any inconsistency. To keep them, use the method :meth:`add_comp` 
        instead.
        
        INPUT:
        
        - ``basis`` -- (default: None) basis in which the components are
          defined; if none is provided, the components are assumed to refer to 
          the module's default basis.
         
        OUTPUT: 
        
        - components in the given basis, as an instance of the 
          class :class:`Components`; if such components did not exist
          previously, they are created.  
        
        EXAMPLES:
        

        """
        if basis is None: 
            basis = self.fmodule.def_basis
        if basis not in self.components:
            if basis not in self.fmodule.known_bases:
                raise ValueError("The " + str(basis) + " has not been " +
                                 "defined on the " + str(self.fmodule))
            self.components[basis] = self._new_comp(basis)
        self._del_derived() # deletes the derived quantities
        self.del_other_comp(basis)
        return self.components[basis]

    def add_comp(self, basis=None):
        r"""
        Return the components in a given basis for assignment, keeping the
        components in other bases. 
        
        To delete the components in other bases, use the method 
        :meth:`set_comp` instead. 
        
        INPUT:
        
        - ``basis`` -- (default: None) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis.
          
        .. WARNING::
        
            If the tensor has already components in other bases, it 
            is the user's responsability to make sure that the components
            to be added are consistent with them. 
         
        OUTPUT: 
        
        - components in the given basis, as an instance of the 
          class :class:`Components`; if such components did not exist
          previously, they are created.  
        
        EXAMPLES:
        
        """
        if basis is None: basis = self.fmodule.def_basis
        if basis not in self.components:
            if basis not in self.fmodule.known_bases:
                raise ValueError("The " + str(basis) + " has not been " +
                                 "defined on the " + str(self.fmodule))
            self.components[basis] = self._new_comp(basis)
        self._del_derived() # deletes the derived quantities
        return self.components[basis]


    def del_other_comp(self, basis=None):
        r"""
        Delete all the components but those corresponding to ``basis``.
        
        """
        if basis is None: basis = self.fmodule.def_basis
        if basis not in self.components:
            raise ValueError("The components w.r.t. the " + 
                             str(basis) + " have not been defined.")
        to_be_deleted = []
        for other_basis in self.components:
            if other_basis != basis:
                to_be_deleted.append(other_basis)
        for other_basis in to_be_deleted:
            del self.components[other_basis]


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

        

        
