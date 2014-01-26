r"""
Tensors on free modules

The class :class:`FreeModuleTensor` implements tensors over a free module `M`,
i.e. elements of the free module `T^{(k,l)}(M)` of tensors of type `(k,l)`
acting as multilinear forms on `M`. 

A tensor of type `(k,l)` is a multilinear map:

.. MATH::

    \underbrace{M^*\times\cdots\times M^*}_{k\ \; \mbox{times}}
    \times \underbrace{M\times\cdots\times M}_{l\ \; \mbox{times}}
    \longrightarrow R
    
where `M^*` stands for the dual of the free module `M` and `R` for the 
commutative ring over which `M` is defined. The integer `k+l`
is called the tensor rank. 

Various derived classes of :class:`FreeModuleTensor` are devoted to specific 
tensors:

* :class:`FreeModuleVector` for elements of `M` (vectors), considered as rank-1 
  contravariant tensors
* :class:`FreeModuleLinForm`for elements of `M` (linear forms), considered as 
  rank-1 covariant tensors
* :class:`FreeModuleAltForm` for alternating forms (fully antisymmetric 
  covariant tensors)

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

EXAMPLES:

    A tensor of type (1,1) on a rank-3 free module over `\ZZ`::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: t = M.tensor((1,1), name='t') ; t
        type-(1,1) tensor t on the rank-3 free module M over the Integer Ring
        sage: t.parent()
        free module of type-(1,1) tensors on the rank-3 free module M over the Integer Ring
        sage: t in M.tensor_module(1,1)
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

from sage.rings.integer import Integer
from sage.structure.element import ModuleElement  
#!# or from sage.structure.element import Element
# to avoid arithmetics defined in ModuleElement ??

from comp import Components, CompWithSym, CompFullySym, CompFullyAntiSym

class FreeModuleTensor(ModuleElement):
    r"""
    Tensor over a free module `M`.
    
    INPUT:
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` (must be an 
      instance of :class:`FiniteFreeModule`)
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

    EXAMPLES:

    A tensor of type (1,1) on a rank-3 free module over `\ZZ`::
    
        sage: M = FiniteFreeModule(ZZ, 3, name='M')
        sage: t = M.tensor((1,1), name='t') ; t
        type-(1,1) tensor t on the rank-3 free module M over the Integer Ring

    Tensors are *Element* objects whose parents are tensor free modules::
    
        sage: t.parent()
        free module of type-(1,1) tensors on the rank-3 free module M over the Integer Ring
        sage: t.parent() is M.tensor_module(1,1)
        True
        
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
        description = "type-(%s,%s) tensor" % \
                           (str(self.tensor_type[0]), str(self.tensor_type[1]))
        if self.name is not None:
            description += " " + self.name
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

    def view(self, basis=None, format_spec=None):
        r"""
        Displays the tensor in terms of its expansion onto a given basis.
        
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        INPUT:
                
        - ``basis`` -- (default: None) basis of the free module with respect to 
          which the tensor is expanded; if none is provided, the module's 
          default basis is assumed
        - ``format_spec`` -- (default: None) format specification passed to 
          ``self.fmodule.output_formatter`` to format the output.

        EXAMPLES:
                    
        """
        from sage.misc.latex import latex
        from utilities import is_atomic, FormattedExpansion
        if basis is None:
            basis = self.fmodule.def_basis
        cobasis = basis.dual_basis
        comp = self.comp(basis)
        terms_txt = []
        terms_latex = []
        n_con = self.tensor_type[0]
        for ind in comp.index_generator():
            ind_arg = ind + (format_spec,)
            coef = comp[ind_arg]
            if coef != 0:
                bases_txt = []
                bases_latex = []
                for k in range(n_con):
                    bases_txt.append(basis[ind[k]].name)
                    bases_latex.append(latex(basis[ind[k]]))
                for k in range(n_con, self.tensor_rank):
                    bases_txt.append(cobasis[ind[k]].name)
                    bases_latex.append(latex(cobasis[ind[k]]))
                basis_term_txt = "*".join(bases_txt)    
                basis_term_latex = r"\otimes ".join(bases_latex)    
                if coef == 1:
                    terms_txt.append(basis_term_txt)
                    terms_latex.append(basis_term_latex)
                elif coef == -1:
                    terms_txt.append("-" + basis_term_txt)
                    terms_latex.append("-" + basis_term_latex)
                else:
                    coef_txt = repr(coef)
                    coef_latex = latex(coef)
                    if is_atomic(coef_txt):
                        terms_txt.append(coef_txt + " " + basis_term_txt)
                    else:
                        terms_txt.append("(" + coef_txt + ") " + 
                                         basis_term_txt)
                    if is_atomic(coef_latex):
                        terms_latex.append(coef_latex + basis_term_latex)
                    else:
                        terms_latex.append(r"\left(" + coef_latex + r"\right)" + 
                                           basis_term_latex)

        if terms_txt == []:
            expansion_txt = "0"
        else:
            expansion_txt = terms_txt[0]
            for term in terms_txt[1:]:
                if term[0] == "-":
                    expansion_txt += " - " + term[1:]
                else:
                    expansion_txt += " + " + term
        if terms_latex == []:
            expansion_latex = "0"
        else:
            expansion_latex = terms_latex[0]
            for term in terms_latex[1:]:
                if term[0] == "-":
                    expansion_latex += term
                else:
                    expansion_latex += "+" + term
        result = FormattedExpansion(self)            
        if self.name is None:
            result.txt = expansion_txt
        else:
            result.txt = self.name + " = " + expansion_txt
        if self.latex_name is None:
            result.latex = expansion_latex
        else:
            result.latex = latex(self) + " = " + expansion_latex
        return result
    
    def symmetries(self):
        r"""
        Print the list of symmetries and antisymmetries.
        
        EXAMPLES:
        
        Various symmetries / antisymmetries for a rank-4 tensor::
        
            sage: from sage.geometry.manifolds.tensor_free_module import FiniteFreeModule
            sage: from sage.geometry.manifolds.free_module_tensor import FreeModuleTensor
            sage: M = FiniteFreeModule(ZZ, 3, name='M')
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
       
    def _new_instance(self):
        r"""
        Create a :class:`FreeModuleTensor` instance of the same tensor type and 
        with the same symmetries.

        This method must be redefined by derived classes of 
        :class:`FreeModuleTensor`.
        
        """
        return FreeModuleTensor(self.fmodule, self.tensor_type, sym=self.sym, 
                                antisym=self.antisym)

    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 
        
        This method, to be called by :meth:`comp`, must be redefined by derived 
        classes to adapt the output to the relevant subclass of 
        :class:`Components`.
        
        """
        fmodule = self.fmodule  # the base free module
        if self.sym == [] and self.antisym == []:
            return Components(fmodule.ring, basis, self.tensor_rank,
                              start_index=fmodule.sindex,
                              output_formatter=fmodule.output_formatter)
        for isym in self.sym:
            if len(isym) == self.tensor_rank:
                return CompFullySym(fmodule.ring, basis, self.tensor_rank,
                                    start_index=fmodule.sindex,
                                    output_formatter=fmodule.output_formatter)
        for isym in self.antisym:
            if len(isym) == self.tensor_rank:
                return CompFullyAntiSym(fmodule.ring, basis, self.tensor_rank, 
                                        start_index=fmodule.sindex,
                                     output_formatter=fmodule.output_formatter)
        return CompWithSym(fmodule.ring, basis, self.tensor_rank, 
                           start_index=fmodule.sindex, 
                           output_formatter=fmodule.output_formatter,
                           sym=self.sym, antisym=self.antisym)        

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
        if self is self.parent()._zero_element: #!# this is maybe not very efficient
            raise ValueError("The zero element cannot be changed.")
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

    def __getitem__(self, indices):
        r"""
        Return a component w.r.t. the free module's default basis.

        INPUT:
        
        - ``indices`` -- list of indices defining the component
    
        """
        return self.comp()[indices]
        
    def __setitem__(self, indices, value):
        r"""
        Set a component w.r.t. the free module's default basis.

        INPUT:
        
        - ``indices`` -- list of indices defining the component
    
        """        
        self.set_comp()[indices] = value


    def copy(self):
        r"""
        Returns an exact copy of ``self``.
        
        The name and the derived quantities are not copied. 
        
        EXAMPLES:
        """
        resu = self._new_instance()
        for basis, comp in self.components.items():
             resu.components[basis] = comp.copy()
        return resu

        
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
    
    - ``fmodule`` -- free module `M` over a commutative ring `R` (must be an 
      instance of :class:`FiniteFreeModule`)
    - ``name`` -- (default: None) name given to the vector
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the vector; 
      if none is provided, the LaTeX symbol is set to ``name``

    """
    def __init__(self, fmodule, name=None, latex_name=None):
        FreeModuleTensor.__init__(self, fmodule, (1,0), name=name, 
                                  latex_name=latex_name)

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "element "
        if self.name is not None:
            description += self.name + " " 
        description += "of the " + str(self.fmodule)
        return description

    def _new_comp(self, basis): 
        r"""
        Create some components in the given basis. 
              
        This method, which is already implemented in 
        :meth:`FreeModuleTensor._new_comp`, is redefined here for efficiency
        """
        fmodule = self.fmodule  # the base free module
        return Components(fmodule.ring, basis, 1, start_index=fmodule.sindex,
                          output_formatter=fmodule.output_formatter)


        
