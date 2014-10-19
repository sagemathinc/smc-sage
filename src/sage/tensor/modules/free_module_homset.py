r"""
Sets of morphisms between free modules

The class :class:`FreeModuleHomset` implements sets (actually free modules) of
homomorphisms between two free modules of finite rank over the same 
commutative ring. 

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

from sage.categories.homset import Homset
from free_module_morphism import FiniteRankFreeModuleMorphism

class FreeModuleHomset(Homset):
    r"""
    Set of homomorphisms between free modules of finite rank
    
    Given two free modules `M` and `N` of respective ranks `m` and `n` over a 
    commutative ring `R`, the class :class:`FreeModuleHomset` implements the 
    set `\mathrm{Hom}(M,N)` of homomorphisms `M\rightarrow N`. This is a 
    *parent* class, whose *elements* are instances of 
    :class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`

    The set `\mathrm{Hom}(M,N)` is actually a free module of rank `mn` over 
    `R`, but this aspect is not taken into account here. 

    INPUT:
    
    - ``fmodule1`` -- free module `M` (domain of the homomorphisms); must be
      an instance of 
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``fmodule2`` -- free module `N` (codomain of the homomorphisms); must be
      an instance of 
      :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    - ``name`` -- (string; default: None) name given to the hom-set; if None, 
      Hom(M,N) will be used
    - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
      hom-set; if None, `\mathrm{Hom}(M,N)` will be used

    EXAMPLES:
    
    Set of homomorphisms between two free modules over `\ZZ`::
    
        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
        sage: H = Hom(M,N) ; H
        Set of Morphisms from rank-3 free module M over the Integer Ring to rank-2 free module N over the Integer Ring in Category of modules over Integer Ring
        sage: type(H)
        <class 'sage.tensor.modules.free_module_homset.FreeModuleHomset_with_category_with_equality_by_id'>
        sage: H.category()
        Category of hom sets in Category of modules over Integer Ring

    Hom-sets are cached::
    
        sage: H is Hom(M,N)
        True

    The LaTeX formatting is::
    
        sage: latex(H)
        \mathrm{Hom}\left(M,N\right)

    As usual, the construction of an element is performed by the ``__call__`` 
    method; the argument can be the matrix representing the morphism in the 
    default bases of the two modules::
    
        sage: e = M.basis('e')
        sage: f = N.basis('f')
        sage: phi = H([[-1,2,0], [5,1,2]]) ; phi
        Generic morphism:
          From: rank-3 free module M over the Integer Ring
          To:   rank-2 free module N over the Integer Ring
        sage: phi.parent() is H
        True

    An example of construction from a matrix w.r.t. bases that are not the
    default ones::
    
        sage: ep = M.basis('ep', latex_symbol=r"e'")
        sage: fp = N.basis('fp', latex_symbol=r"f'")
        sage: phi2 = H([[3,2,1], [1,2,3]], bases=(ep,fp)) ; phi2
        Generic morphism:
          From: rank-3 free module M over the Integer Ring
          To:   rank-2 free module N over the Integer Ring

    The zero element::
    
        sage: z = H.zero() ; z
        Generic morphism:
          From: rank-3 free module M over the Integer Ring
          To:   rank-2 free module N over the Integer Ring
        sage: z.matrix(e,f)
        [0 0 0]
        [0 0 0]

    """

    Element = FiniteRankFreeModuleMorphism

    def __init__(self, fmodule1, fmodule2, name=None, latex_name=None):
        r"""
        TESTS::
        
            sage: from sage.tensor.modules.free_module_homset import FreeModuleHomset
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: FreeModuleHomset(M, N)
            Set of Morphisms from rank-3 free module M over the Integer Ring to rank-2 free module N over the Integer Ring in Category of modules over Integer Ring
            sage: H = FreeModuleHomset(M, N, name='L(M,N)', latex_name=r'\mathcal{L}(M,N)')
            sage: latex(H)
            \mathcal{L}(M,N)
        
        """
        from finite_rank_free_module import FiniteRankFreeModule
        if not isinstance(fmodule1, FiniteRankFreeModule):
            raise TypeError("fmodule1 = " + str(fmodule1) + " is not an " + 
                            "instance of FiniteRankFreeModule.")
        if not isinstance(fmodule2, FiniteRankFreeModule):
            raise TypeError("fmodule1 = " + str(fmodule2) + " is not an " + 
                            "instance of FiniteRankFreeModule.")
        if fmodule1.base_ring() != fmodule2.base_ring():
            raise ValueError("The domain and codomain are not defined over " + 
                            "the same ring.")
        Homset.__init__(self, fmodule1, fmodule2)
        if name is None:
            self._name = "Hom(" + fmodule1._name + "," + fmodule2._name + ")"
        else:
            self._name = name
        if latex_name is None:
            self._latex_name = \
                    r"\mathrm{Hom}\left(" + fmodule1._latex_name + "," + \
                    fmodule2._latex_name + r"\right)"
        else:
            self._latex_name = latex_name
        
    def _latex_(self):
        r"""
        LaTeX representation of the object.
        
        EXAMPLES::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: H = Hom(M,N)
            sage: H._latex_()
            '\\mathrm{Hom}\\left(M,N\\right)'
            sage: latex(H)  # indirect doctest
            \mathrm{Hom}\left(M,N\right)
        
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name
            
    def __call__(self, *args, **kwds):
        r"""
        To bypass Homset.__call__ and enforce the call to _element_constructor_
        """
        from sage.structure.parent import Parent
        return Parent.__call__(self, *args, **kwds)
        
    #### Methods required for any Parent 

    def _element_constructor_(self, matrix_rep, bases=None, name=None, 
                              latex_name=None):
        r"""
        Construct an element of ``self``, i.e. a homomorphism M --> N, where
        M is the domain of ``self`` and N its codomain. 
        
        INPUT:
        
        - ``matrix_rep`` -- matrix representation of the homomorphism with 
          respect to the bases ``basis1`` and ``basis2``; this entry can actually
          be any material from which a matrix of size rank(N)*rank(M) can be 
          constructed
        - ``bases`` -- (default: None) pair (basis_M, basis_N) defining the 
          matrix representation, basis_M being a basis of module `M` and
          basis_N a basis of module `N` ; if None the pair formed by the 
          default bases of each module is assumed. 
        - ``name`` -- (string; default: None) name given to the homomorphism
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
          homomorphism; if None, ``name`` will be used. 
          
        """
        return self.element_class(self, matrix_rep, bases=bases, name=name, 
                                  latex_name=latex_name)
    
    def _an_element_(self):
        r"""
        Construct some (unamed) element.
        
        EXAMPLE::
        
        """
        ring = self.base_ring()
        m = self.domain().rank()
        n = self.codomain().rank()
        matrix_rep = [[ring.an_element() for i in range(m)] for j in range(n)]
        return self.element_class(self, matrix_rep)
            
    #### End of methods required for any Parent 
