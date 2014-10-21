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

    def _call_(self, element):
        r"""
        Action of the homomorphism ``self`` on some free module element

        INPUT:
        
        - ``element`` -- element of the domain of ``self``
        
        OUTPUT:
        
        - the image of ``element`` by ``self``
        
        EXAMPLE:
        
        Images of a homomorphism between two `\ZZ`-modules::
        
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: N = FiniteRankFreeModule(ZZ, 2, name='N')
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]], name='phi', latex_name=r'\phi')
            sage: v = M([1,2,3], basis=e, name='v')
            sage: w = phi(v) ; w
            element phi(v) of the rank-2 free module N over the Integer Ring
            sage: w.view()
            phi(v) = 3 f_0 + 13 f_1
            
        Tests::
        
            sage: for i in range(2):
            ....:     print w[i] == sum( phi.matrix()[i,j]*v[j] for j in range(3) ),
            ....:     
            True True
            sage: phi.matrix(e,f)
            [-1  2  0]
            [ 5  1  2]
            sage: phi(e[0]).view()
            phi(e_0) = -f_0 + 5 f_1
            sage: phi(e[1]).view()
            phi(e_1) = 2 f_0 + f_1
            sage: phi(e[2]).view()
            phi(e_2) = 2 f_1
            
        Image of an element that is not defined on the default basis::
        
            sage: a = M.automorphism()
            sage: a[0,2], a[1,0], a[2,1] = 1, -1, -1
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: v = M([1,2,3], basis=ep, name='v')
            sage: w = phi(v) ; w
            element phi(v) of the rank-2 free module N over the Integer Ring
            sage: w.view()
            phi(v) = -5 f_0 + 10 f_1
            sage: for i in range(2):
            ....:     print w[i] == sum( phi.matrix(ep,f)[i,j]*v[ep,j] for j in range(3) ),
            ....:     
            True True
            
        """
        dom = self.parent().domain()
        sindex = dom._sindex
        codom = self.parent().codomain()
        basis_codom = codom.default_basis()
        # Search for a common basis to compute the result
        for basis in element._components:
            try:
                self.matrix(basis, basis_codom)
                basis_dom = basis
                break
            except ValueError:
                continue
        else:
            raise ValueError("No common basis found to evaluate the image " + 
                             "of " + str(element) + " by " + str(self) + ".")
        # Components of the result obtained by matrix multiplication
        mat = self.matrix(basis_dom, basis_codom)
        vcomp = element._components[basis_dom]
        tresu = []
        for i in range(codom.rank()):
            s = 0
            for j in range(dom.rank()):
                s += mat[i,j] * vcomp[[j+sindex]]
            tresu.append(s)
        # Name of the result
        if self._name is not None and element._name is not None:
            resu_name = self._name + '(' + element._name + ')'
        else:
            resu_name = None
        if self._latex_name is not None and element._latex_name is not None:
            resu_latex_name = self._latex_name + r'\left(' + \
                              element._latex_name + r'\right)'
        else:
            resu_latex_name = None
        # Creation of the result
        return codom(tresu, basis=basis_codom, name=resu_name, 
                     latex_name=resu_latex_name)


    def matrix(self, basis1=None, basis2=None):
        r"""
        Return the matrix of ``self`` w.r.t to a pair of bases. 

        If the matrix is not known already, it is computed from the matrix in
        another pair of bases by means of the change-of-bases formula. 
        
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
            sage: e = M.basis('e') ; f = N.basis('f')
            sage: phi = M.hom(N, [[-1,2,0], [5,1,2]])
            sage: phi.matrix()     # default bases
            [-1  2  0]
            [ 5  1  2]
            sage: phi.matrix(e,f)  # bases explicited
            [-1  2  0]
            [ 5  1  2]
            sage: type(phi.matrix())
            <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>

        Matrix in bases different from those in which the homomorphism has
        been defined::
        
            sage: a = M.automorphism()
            sage: a[0,2], a[1,0], a[2,1] = 1, -1, -1
            sage: a[:]
            [ 0  0  1]
            [-1  0  0]
            [ 0 -1  0]
            sage: ep = e.new_basis(a, 'ep', latex_symbol="e'")
            sage: b = N.automorphism()
            sage: b[0,1], b[1,0] = -1, 1
            sage: b[:]
            [ 0 -1]
            [ 1  0]
            sage: fp = f.new_basis(b, 'fp', latex_symbol="f'")
            sage: phi.matrix(ep, fp)
            [-1 -2  5]
            [ 2  0  1]

        Check of the change-of-bases formula::
        
            sage: phi.matrix(ep, fp) == matrix(b.inverse()[:]) * phi.matrix(e,f) * matrix(a[:])
            True

        Single change of basis::
        
            sage: phi.matrix(ep, f)
            [-2  0 -1]
            [-1 -2  5]
            sage: phi.matrix(ep,f) == phi.matrix(e,f) * matrix(a[:])
            True
            sage: phi.matrix(e, fp)
            [ 5  1  2]
            [ 1 -2  0]
            sage: phi.matrix(e, fp) == matrix(b.inverse()[:]) * phi.matrix(e,f)
            True

        """
        from sage.matrix.constructor import matrix
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
            b1_list = [bases[0] for bases in self._matrices]
            b2_list = [bases[1] for bases in self._matrices]
            if basis1 in b1_list:
                for b2 in b2_list:
                    if (basis2, b2) in fmodule2._basis_changes:
                        nb2 = b2
                        break
                else:
                    raise ValueError("No start basis could be found for " + 
                                     "applying the change-of-bases formula.")
                change2 = fmodule2._basis_changes[(basis2, nb2)]
                mat2 = matrix( [[change2[[i,j]] for j in fmodule2.irange()] 
                                                  for i in fmodule2.irange()] )
                self._matrices[(basis1, basis2)] = \
                                            mat2 * self._matrices[(basis1,nb2)]
            elif basis2 in b2_list:
                for b1 in b1_list:
                    if (b1, basis1) in fmodule1._basis_changes:
                        nb1 = b1
                        break
                else:
                    raise ValueError("No start basis could be found for " + 
                                     "applying the change-of-bases formula.")
                change1 = fmodule1._basis_changes[(nb1, basis1)]
                mat1 = matrix( [[change1[[i,j]] for j in fmodule1.irange()]
                                                  for i in fmodule1.irange()] )
                self._matrices[(basis1, basis2)] = \
                                            self._matrices[(nb1,basis2)] * mat1
            else: # most general change-of-basis formula
                for (b1, b2) in self._matrices:
                    if (b1, basis1) in fmodule1._basis_changes and \
                       (basis2, b2) in fmodule2._basis_changes:
                        nb1, nb2 = b1, b2
                        break
                else:
                    raise ValueError("No start basis could be found for " + 
                                     "applying the change-of-bases formula.")
                change1 = fmodule1._basis_changes[(nb1, basis1)]
                change2 = fmodule2._basis_changes[(basis2, nb2)]
                mat1 = matrix( [[change1[[i,j]] for j in fmodule1.irange()]
                                                  for i in fmodule1.irange()] )
                mat2 = matrix( [[change2[[i,j]] for j in fmodule2.irange()] 
                                                  for i in fmodule2.irange()] )
                self._matrices[(basis1, basis2)] = \
                                        mat2 * self._matrices[(nb1,nb2)] * mat1
        return self._matrices[(basis1, basis2)]

    def _common_bases(self, other):
        r"""
        Return a pair of bases in which ``self`` and ``other`` have a known 
        matrix representation. 
        
        INPUT:
        
        - ``other`` -- another homomorphism in the same hom-set
        
        OUTPUT:
        
        - a pair of bases in which ``self`` and ``other`` have a known 
          matrix representation. 
          
        """
        resu = None
        for bases in self._matrices:
            try:
                other.matrix(bases)
                resu = bases
                break
            except ValueError:
                continue
        if resu is None:
            for bases in other._matrices:
                try:
                    self.matrix(bases)
                    resu = bases
                    break
                except ValueError:
                    continue
        return resu
