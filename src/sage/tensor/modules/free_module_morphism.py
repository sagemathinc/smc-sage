r"""
Morphisms between free modules

The class :class:`FiniteRankFreeModuleMorphism` implements homomorphisms
between two free modules of finite rank over the same commutative ring. 

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

from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from finite_rank_free_module import FiniteRankFreeModule

class FiniteRankFreeModuleMorphism(Morphism):
    r"""
    Homomorphism between free modules of finite rank
    
    This is an *element* class, whose *parent* class is
    :class:`~sage.tensor.modules.free_module_homset.FreeModuleHomset`.

    An instance of this class is a homomorphism
    
    .. MATH::

        \phi:\ M \longrightarrow N, 
        
    where `M` and `N` are two free modules of finite rank over the same 
    commutative ring `R`. 
    
    INPUT:
    
    - ``parent`` -- hom-set Hom(M,N) to which the homomorphism belongs
    - ``matrix_rep`` -- matrix representation of the homomorphism with 
      respect to the bases ``bases``; this entry can actually
      be any material from which a matrix of size rank(N)*rank(M) of 
      elements of `R` can be constructed
    - ``bases`` -- (default: None) pair (basis_M, basis_N) defining the 
      matrix representation, basis_M being a basis of module `M` and
      basis_N a basis of module `N` ; if None the pair formed by the 
      default bases of each module is assumed. 
    - ``name`` -- (string; default: None) name given to the homomorphism
    - ``latex_name`` -- (string; default: None) LaTeX symbol to denote the 
      homomorphism; if None, ``name`` will be used. 
    
    EXAMPLES:
    
    A homomorphism between two free modules over `\ZZ` is contructed 
    as an element of the corresponding hom-set, by means of the function
    ``__call__``::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
        sage: e = M.basis('e') ; f = N.basis('f')
        sage: H = Hom(M,N) ; H
        Set of Morphisms from rank-3 free module M over the Integer Ring to rank-2 free module N over the Integer Ring in Category of modules over Integer Ring
        sage: phi = H([[2,-1,3], [1,0,-4]], name='phi', latex_name=r'\phi') ; phi
        Generic morphism:
          From: rank-3 free module M over the Integer Ring
          To:   rank-2 free module N over the Integer Ring

    Since no bases have been specified in the argument list, the provided 
    matrix is relative to the default bases of modules M and N, so that
    the above is equivalent to::
    
        sage: phi = H([[2,-1,3], [1,0,-4]], bases=(e,f), name='phi', latex_name=r'\phi') ; phi
        Generic morphism:
          From: rank-3 free module M over the Integer Ring
          To:   rank-2 free module N over the Integer Ring

    An alternative way to construct a homomorphism is to call the method
    :meth:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule.hom`
    on the domain::
    
        sage: phi = M.hom(N, [[2,-1,3], [1,0,-4]], bases=(e,f), name='phi', latex_name=r'\phi') ; phi
        Generic morphism:
          From: rank-3 free module M over the Integer Ring
          To:   rank-2 free module N over the Integer Ring
    
    The parent of a homomorphism is of course the corresponding hom-set::

        sage: phi.parent() is H
        True
        sage: phi.parent() is Hom(M,N)
        True

    Due to Sage's category scheme, the class of the homomorphism phi is
    actually a derived class of :class:`FiniteRankFreeModuleMorphism`::
    
        sage: type(phi)
        <class 'sage.tensor.modules.free_module_morphism.FreeModuleHomset_with_category_with_equality_by_id.element_class'>
        sage: isinstance(phi, sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism)
        True

    The domain and codomain of the homomorphism are returned respectively by 
    the methods ``domain()`` and ``codomain()``, which is are implemented as
    Sage's constant functions::
    
        sage: phi.domain()
        rank-3 free module M over the Integer Ring
        sage: phi.codomain()
        rank-2 free module N over the Integer Ring
        sage: type(phi.domain)
        <type 'sage.misc.constant_function.ConstantFunction'>
    
    The matrix of the homomorphism with respect to a pair of bases is 
    returned by the method :meth:`matrix`::
    
        sage: phi.matrix(e,f)
        [ 2 -1  3]
        [ 1  0 -4]

    """
    def __init__(self, parent, matrix_rep, bases=None, name=None, 
                 latex_name=None):
        r"""
        TESTS::
        """
        from sage.matrix.constructor import matrix
        from sage.misc.constant_function import ConstantFunction
        Morphism.__init__(self, parent)
        fmodule1 = parent.domain()
        fmodule2 = parent.codomain()
        if bases is None:
            def_basis1 = fmodule1.default_basis()
            if def_basis1 is None:
                raise ValueError("The " + str(fmodule1) + " has no default " + 
                                 "basis.")
            def_basis2 = fmodule2.default_basis()
            if def_basis2 is None:
                raise ValueError("The " + str(fmodule2) + " has no default " + 
                                 "basis.")
            bases = (def_basis1, def_basis2)
        else:
            bases = tuple(bases)  # insures bases is a tuple
            if len(bases) != 2:
                raise TypeError("The argument bases must contain 2 bases.")
            if bases[0] not in fmodule1.bases():
                raise TypeError(str(bases[0]) + " is not a basis on the " + \
                                str(fmodule1) + ".")
            if bases[1] not in fmodule2.bases():
                raise TypeError(str(bases[1]) + " is not a basis on the " + \
                                str(fmodule2) + ".")

        if isinstance(matrix_rep, ConstantFunction):
            # Construction of the zero morphism
            if matrix_rep().is_zero():
                matrix_rep = 0
        ring = fmodule1.base_ring()
        n1 = fmodule1.rank()
        n2 = fmodule2.rank()
        self._matrices = {bases: matrix(ring, n2, n1, matrix_rep)}
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        
        EXAMPLES::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi', latex_name=r'\Phi')
            sage: phi._latex_()
            '\\Phi'
            sage: latex(phi)  # indirect doctest
            \Phi
        
        ::
        
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='F')
            sage: phi._latex_()
            'F'
            
        ::
        
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi._latex_()
            '\\mbox{Generic morphism:\n  From: rank-3 free module M over the Integer Ring\n  To:   rank-2 free module N over the Integer Ring}'

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def __nonzero__(self):
        r"""
        Return True if ``self`` is nonzero and False otherwise. 
        
        This method is called by self.is_zero(). 
        
        EXAMPLES::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[2,-1,3], [1,0,-4]])
            sage: phi.__nonzero__()
            True
            sage: phi.is_zero() # indirect doctest
            False
            sage: phi = M.hom(N, 0)
            sage: phi.__nonzero__()
            False
            sage: phi.is_zero() # indirect doctest
            True
            sage: Hom(M,N).zero().__nonzero__()
            False

        """
        # Some matrix representation is picked at random:
        matrix_rep = self._matrices.items()[0][1]
        return not matrix_rep.is_zero()         

    def matrix(self, basis1=None, basis2=None):
        r"""
        Return the matrix of ``self`` w.r.t to a pair of bases. 

        INPUT:
        
        - ``basis1`` -- (default: None) basis of the domain of ``self``; if 
          None, the domain's default basis is assumed
        - ``basis2`` -- (default: None) basis of the codomain of ``self``; if 
          None, the codomain's default basis is assumed
        
        OUTPUT:
        
        - the matrix representing representing the homomorphism ``self`` w.r.t
          to bases ``basis1`` and ``basis2``
        
        EXAMPLES:
        
        Matrix of a homomorphism between two `\ZZ`-modules::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e')
            sage: f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi.matrix()     # default bases
            [-1  2  0]
            [ 5  1  2]
            sage: phi.matrix(e,f)  # bases explicited
            [-1  2  0]
            [ 5  1  2]
            sage: type(phi.matrix())
            <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>

        """
        fmodule1 = self.domain()
        fmodule2 = self.codomain()
        if basis1 is None:
            basis1 = fmodule1.default_basis()
        elif basis1 not in fmodule1.bases():
            raise TypeError(str(basis1) + " is not a basis on the " + \
                            str(fmodule1) + ".")
        if basis2 is None:
            basis2 = fmodule2.default_basis()
        elif basis2 not in fmodule2.bases():
            raise TypeError(str(basis2) + " is not a basis on the " + \
                            str(fmodule2) + ".")
        if (basis1, basis2) not in self._matrices:
            raise KeyError("The matrix of " + str(self) + " is not known " + \
                           "w.r.t. " + str(basis1) + " and " + str(basis2) + \
                           ".")
        return self._matrices[(basis1, basis2)]
