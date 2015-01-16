"""
Free module automorphisms

AUTHORS:

- Eric Gourgoulhon (2015): initial version

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.element import MultiplicativeGroupElement
from sage.tensor.modules.free_module_tensor import FreeModuleTensor

class FreeModuleAutomorphism(FreeModuleTensor, MultiplicativeGroupElement):
    r"""
    Automorphism (considered as a type-`(1,1)` tensor) on a free module.

    INPUT:

    - ``fmodule`` -- free module `M` over a commutative ring `R`
      (must be an instance of :class:`FiniteRankFreeModule`)
    - ``name`` -- (default: ``None``) name given to the automorphism
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      automorphism; if none is provided, the LaTeX symbol is set to ``name``
    - ``is_identity`` -- (default: ``False``) determines whether the
      constructed object is the identity automorphism, i.e. the identity map
      of `M` considered as an automorphism (the identity element of the
      general linear group)

    EXAMPLES:

    Automorphism of a rank-2 free module (vector space) on `\QQ`::

        sage: M = FiniteRankFreeModule(QQ, 2, name='M')
        sage: a = M.automorphism(name='A') ; a
        Automorphism A of the Rank-2 free module M over the Rational Field
        sage: a.parent()
        General linear group of the Rank-2 free module M over the Rational
         Field
        sage: a.parent() is M.general_linear_group()
        True

    Automorphisms are tensors of type `(1,1)`::

        sage: a.tensor_type()
        (1, 1)
        sage: a.tensor_rank()
        2

    Setting the components in a basis::

        sage: e = M.basis('e') ; e
        Basis (e_0,e_1) on the Rank-2 free module M over the Rational Field
        sage: a[:] = [[1, 2], [-1, 3]]
        sage: a[:]
        [ 1  2]
        [-1  3]
        sage: a.display(basis=e)
        A = e_0*e^0 + 2 e_0*e^1 - e_1*e^0 + 3 e_1*e^1

    The inverse automorphism is obtained via the method :meth:`inverse`,
    or the operator ~, or the exponent -1::

        sage: b = a.inverse() ; b
        Automorphism A^(-1) of the Rank-2 free module M over the Rational Field
        sage: b is ~a
        True
        sage: b is a^(-1)
        True
        sage: b.display(basis=e)
        A^(-1) = 3/5 e_0*e^0 - 2/5 e_0*e^1 + 1/5 e_1*e^0 + 1/5 e_1*e^1
        sage: b[:]
        [ 3/5 -2/5]
        [ 1/5  1/5]
        sage: a[:] * b[:]  # check that b is indeed the inverse of a
        [1 0]
        [0 1]

    """
    def __init__(self, fmodule, name=None, latex_name=None, is_identity=False):
        r"""
        TESTS::
    
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: from sage.tensor.modules.free_module_automorphism import FreeModuleAutomorphism
            sage: a = FreeModuleAutomorphism(M, name='a')
            sage: a[e,:] = [[-1,0,0],[0,1,2],[0,1,3]]
            sage: TestSuite(a).run(skip="_test_category") # see below
    
        In the above test suite, _test_category fails because a is not an
        instance of a.parent().category().element_class. Actually automorphism
        must be constructed via FreeModuleLinearGroup.element_class and
        not by a direct call to FreeModuleAutomorphism::

            sage: a = M.general_linear_group().element_class(M, name='a') 
            sage: a[e,:] = [[-1,0,0],[0,1,2],[0,1,3]]
            sage: TestSuite(a).run()

        Test suite on the identity map::

            sage: id = M.general_linear_group().one() 
            sage: TestSuite(id).run()

        Test suite on the automorphism obtained as GL.an_element()::

            sage: b = M.general_linear_group().an_element()
            sage: TestSuite(b).run()

        """
        if is_identity:
            if name is None:
                name = 'Id'
            if latex_name is None and name == 'Id':
                latex_name = r'\mathrm{Id}'
        FreeModuleTensor.__init__(self, fmodule, (1,1), name=name,
                                  latex_name=latex_name,
                                  parent=fmodule.general_linear_group())
        # MultiplicativeGroupElement attributes:
        # - none
        # Local attributes:
        self._is_identity = is_identity
        self._inverse = None    # inverse automorphism not set yet
        self._matrices = {}

    #### SageObject methods ####

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: M.automorphism()
            Automorphism of the Rank-3 free module M over the Rational Field
            sage: M.automorphism(name='a')
            Automorphism a of the Rank-3 free module M over the Rational Field
            sage: M.identity_map()
            Identity map of the Rank-3 free module M over the Rational Field

        """
        if self._is_identity:
            description = "Identity map "
        else:
            description = "Automorphism "
            if self._name is not None:
                description += self._name + " "
        description += "of the {}".format(self._fmodule)
        return description

    #### End of SageObject methods ####

    #### FreeModuleTensor methods ####

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.automorphism(name='a')
            sage: a._new_instance()
            Automorphism of the Rank-3 free module M over the Integer Ring
            sage: Id = M.identity_map()
            sage: Id._new_instance()
            Automorphism of the Rank-3 free module M over the Integer Ring

        """
        return self.__class__(self._fmodule)

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: a = M.automorphism(name='a')
            sage: e = M.basis('e')
            sage: a[e,:] = [[1,0,-1], [0,3,0], [0,0,2]]
            sage: b = a.inverse()
            sage: a._inverse
            Automorphism a^(-1) of the Rank-3 free module M over the Rational Field
            sage: a._del_derived()
            sage: a._inverse  # has been reset to None

        """
        # First delete the derived quantities pertaining to FreeModuleTensor:
        FreeModuleTensor._del_derived(self)
        # Then reset the inverse automorphism to None:
        if self._inverse is not None:
            self._inverse._inverse = None  # (it was set to self)
            self._inverse = None
        # and delete the matrices:
        self._matrices.clear()

    def _new_comp(self, basis):
        r"""
        Create some (uninitialized) components of ``self`` in a given basis.

        INPUT:

        - ``basis`` -- basis of the free module on which ``self`` is defined

        OUTPUT:

        - an instance of :class:`~sage.tensor.modules.comp.Components` or,
          if ``self`` is the identity, of the subclass
          :class:`~sage.tensor.modules.comp.KroneckerDelta`

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism()
            sage: a._new_comp(e)
            2-indices components w.r.t. Basis (e_0,e_1,e_2) on the Rank-3 free
             module M over the Integer Ring
            sage: id = M.identity_map()
            sage: id._new_comp(e)
            Kronecker delta of size 3x3
            sage: type(id._new_comp(e))
            <class 'sage.tensor.modules.comp.KroneckerDelta'>

        """
        from comp import KroneckerDelta
        if self._is_identity:
            fmodule = self._fmodule
            return KroneckerDelta(fmodule._ring, basis,
                                  start_index=fmodule._sindex,
                                  output_formatter=fmodule._output_formatter)
        return FreeModuleTensor._new_comp(self, basis)


    def components(self, basis=None, from_basis=None):
        r"""
        Return the components in a given basis.

        If the components are not known already, they are computed by the
        tensor change-of-basis formula from components in another basis.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          required; if none is provided, the components are assumed to refer
          to the module's default basis
        - ``from_basis`` -- (default: ``None``) basis from which the
          required components are computed, via the tensor change-of-basis
          formula, if they are not known already in the basis ``basis``;
          if none, a basis from which both the components and a change-of-basis
          to ``basis`` are known is selected.

        OUTPUT:

        - components in the basis ``basis``, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`,
          or, for the identity automorphism, of the subclass
          :class:`~sage.tensor.modules.comp.KroneckerDelta`

        EXAMPLES:

        Components of an automorphism on a rank-3 free `\ZZ`-module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: a = M.automorphism([[-1,0,0],[0,1,2],[0,1,3]], name='a')
            sage: a.components(e)
            2-indices components w.r.t. Basis (e_1,e_2,e_3) on the Rank-3 free
             module M over the Integer Ring
            sage: a.components(e)[:]
            [-1  0  0]
            [ 0  1  2]
            [ 0  1  3]

        Since e is the module's default basis, it can be omitted::

            sage: a.components() is a.components(e)
            True

        A shortcut is ``a.comp()``::

            sage: a.comp() is a.components()
            True
            sage: a.comp(e) is a.components()
            True

        Components in another basis::

            sage: f1 = -e[2]
            sage: f2 = 4*e[1] + 3*e[3]
            sage: f3 = 7*e[1] + 5*e[3]
            sage: f = M.basis('f', from_family=(f1,f2,f3))
            sage: a.components(f)
            2-indices components w.r.t. Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring
            sage: a.components(f)[:]
            [  1  -6 -10]
            [ -7  83 140]
            [  4 -48 -81]

        Some check of the above matrix::

            sage: a(f[1]).display(f)
            a(f_1) = f_1 - 7 f_2 + 4 f_3
            sage: a(f[2]).display(f)
            a(f_2) = -6 f_1 + 83 f_2 - 48 f_3
            sage: a(f[3]).display(f)
            a(f_3) = -10 f_1 + 140 f_2 - 81 f_3

        Components of the identity map::

            sage: id = M.identity_map()
            sage: id.components(e)
            Kronecker delta of size 3x3
            sage: id.components(e)[:]
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: id.components(f)
            Kronecker delta of size 3x3
            sage: id.components(f)[:]
            [1 0 0]
            [0 1 0]
            [0 0 1]

        """
        if self._is_identity:
            if basis is None:
                basis = self._fmodule._def_basis
            if basis not in self._components:
                self._components[basis] = self._new_comp(basis)
            return self._components[basis]
        else:
            return FreeModuleTensor.components(self, basis=basis,
                                               from_basis=from_basis)

    comp = components

    def set_comp(self, basis=None):
        r"""
        Return the components in a given basis for assignment.

        The components with respect to other bases are deleted, in order to
        avoid any inconsistency. To keep them, use the method :meth:`add_comp`
        instead.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis

        OUTPUT:

        - components in the given basis, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`; if such components did not exist
          previously, they are created.

        EXAMPLE:

        For the identity map, it is not permitted to set components::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.identity_map()
            sage: a.set_comp(e)
            Traceback (most recent call last):
            ...
            TypeError: the components of the identity map cannot be changed

        """
        if self._is_identity:
            raise TypeError("the components of the identity map cannot be " +
                            "changed")
        else:
            return FreeModuleTensor.set_comp(self, basis=basis)

    def add_comp(self, basis=None):
        r"""

        Return the components in a given basis for assignment, keeping the
        components in other bases.

        To delete the components in other bases, use the method
        :meth:`set_comp` instead.

        INPUT:

        - ``basis`` -- (default: ``None``) basis in which the components are
          defined; if none is provided, the components are assumed to refer to
          the module's default basis

        .. WARNING::

            If the tensor has already components in other bases, it
            is the user's responsability to make sure that the components
            to be added are consistent with them.

        OUTPUT:

        - components in the given basis, as an instance of the
          class :class:`~sage.tensor.modules.comp.Components`;
          if such components did not exist previously, they are created

        EXAMPLE:

        For the identity map, it is not permitted to set components::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.identity_map()
            sage: a.add_comp(e)
            Traceback (most recent call last):
            ...
            TypeError: the components of the identity map cannot be changed

        """
        if self._is_identity:
            raise TypeError("the components of the identity map cannot be " +
                            "changed")
        else:
            return FreeModuleTensor.add_comp(self, basis=basis)

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`FreeModuleTensor.__call__` to allow for a single
        argument (module element).

        EXAMPLES:

        Call with a single argument --> return a module element::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: a = M.automorphism([[-1,0,0],[0,1,2],[0,1,3]], name='a')
            sage: v = M([2,1,4], name='v')
            sage: s = a.__call__(v) ; s
            Element a(v) of the Rank-3 free module M over the Integer Ring
            sage: s.display()
            a(v) = -2 e_1 + 9 e_2 + 13 e_3
            sage: s == a(v)
            True
            sage: s == a.contract(v)
            True

        Call with two arguments (:class:`FreeModuleTensor` behaviour)
        --> return a scalar::

            sage: b = M.linear_form(name='b')
            sage: b[:] = 7, 0, 2
            sage: a.__call__(b,v)
            12
            sage: a(b,v) == a.__call__(b,v)
            True
            sage: a(b,v) == s(b)
            True

        Identity map with a single argument --> return a module element::

            sage: id = M.identity_map()
            sage: s = id.__call__(v) ; s
            Element v of the Rank-3 free module M over the Integer Ring
            sage: s == v
            True
            sage: s == id(v)
            True
            sage: s == id.contract(v)
            True

        Identity map with two arguments (:class:`FreeModuleTensor` behaviour)
        --> return a scalar::

            sage: id.__call__(b,v)
            22
            sage: id(b,v) == id.__call__(b,v)
            True
            sage: id(b,v) == b(v)
            True

        """
        from free_module_tensor import FiniteRankFreeModuleElement
        if len(arg) > 1:
            # The automorphism acting as a type-(1,1) tensor on a pair
            # (linear form, module element), returning a scalar:
            if self._is_identity:
                if len(arg) != 2:
                    raise TypeError("wrong number of arguments")
                linform = arg[0]
                if linform._tensor_type != (0,1):
                    raise TypeError("the first argument must be a linear form")                
                vector = arg[1]
                if not isinstance(vector, FiniteRankFreeModuleElement):
                    raise TypeError("the second argument must be a module" +
                                    " element")
                return linform(vector)
            else: # self is not the identity automorphism:
                return FreeModuleTensor.__call__(self, *arg)
        # The automorphism acting as such, on a module element, returning a
        # module element:
        vector = arg[0]
        if not isinstance(vector, FiniteRankFreeModuleElement):
            raise TypeError("the argument must be an element of a free module")
        if self._is_identity:
            return vector
        basis = self.common_basis(vector)
        t = self._components[basis]
        v = vector._components[basis]
        fmodule = self._fmodule
        result = vector._new_instance()
        for i in fmodule.irange():
            res = 0
            for j in fmodule.irange():
                res += t[[i,j]]*v[[j]]
            result.set_comp(basis)[i] = res
        # Name of the output:
        result._name = None
        if self._name is not None and vector._name is not None:
            result._name = self._name + "(" + vector._name + ")"
        # LaTeX symbol for the output:
        result._latex_name = None
        if self._latex_name is not None and vector._latex_name is not None:
            result._latex_name = self._latex_name + r"\left(" + \
                              vector._latex_name + r"\right)"
        return result

    #### End of FreeModuleTensor methods ####

    #### MultiplicativeGroupElement methods ####

    def __invert__(self):
        r"""
        Return the inverse automorphism.

        OUTPUT:

        - instance of :class:`FreeModuleAutomorphism` representing the
          automorphism that is the inverse of ``self``.

        EXAMPLES:

        Inverse of an automorphism on a rank-3 free module::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: a = M.automorphism(name='A')
            sage: e = M.basis('e')
            sage: a[:] = [[1,0,-1], [0,3,0], [0,0,2]]
            sage: b = a.__invert__() ; b
            Automorphism A^(-1) of the Rank-3 free module M over the Rational Field
            sage: b[:]
            [  1   0 1/2]
            [  0 1/3   0]
            [  0   0 1/2]

        We may check that ``b`` is the inverse of ``a`` by performing the
        matrix product of the components in the basis ``e``::

            sage: a[:] * b[:]
            [1 0 0]
            [0 1 0]
            [0 0 1]

        Another check is of course::

            sage: b.__invert__() == a
            True

        The inverse of the identity map is of course itself::

            sage: id = M.identity_map()
            sage: id.__invert__() is id
            True

        """
        from sage.matrix.constructor import matrix
        from comp import Components
        if self._is_identity:
            return self
        if self._inverse is None:
            if self._name is None:
                inv_name = None
            else:
                inv_name = self._name  + '^(-1)'
            if self._latex_name is None:
                inv_latex_name = None
            else:
                inv_latex_name = self._latex_name + r'^{-1}'
            fmodule = self._fmodule
            si = fmodule._sindex
            nsi = fmodule._rank + si
            self._inverse = self.__class__(fmodule, inv_name, inv_latex_name)
            for basis in self._components:
                try:
                    mat = self.matrix(basis)
                except (KeyError, ValueError):
                    continue
                mat_inv = mat.inverse()
                cinv = Components(fmodule._ring, basis, 2, start_index=si,
                                  output_formatter=fmodule._output_formatter)
                for i in range(si, nsi):
                    for j in range(si, nsi):
                        cinv[i, j] = mat_inv[i-si,j-si]
                self._inverse._components[basis] = cinv
            self._inverse._inverse = self
        return self._inverse

    inverse = __invert__

    def _mul_(self, other):
        r"""
        Automorphism composition.

        This implements the group law of GL(M), M being the module of ``self``. 

        INPUT:

        - ``other`` -- an automorphism of the same module as ``self``
        OUPUT:

        - the automorphism resulting from the composition of ``other`` and
        ``self.``

        EXAMPLES::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[1,2],[1,3]]) 
            sage: b = M.automorphism([[3,4],[5,7]])
            sage: c = a._mul_(b) ; c
            Automorphism of the Rank-2 free module M over the Integer Ring
            sage: c.matrix()
            [13 18]
            [18 25]

        TESTS::
        
            sage: c.parent() is a.parent()
            True
            sage: c.matrix() == a.matrix() * b.matrix()
            True
            sage: c(e[0]) == a(b(e[0]))
            True
            sage: c(e[1]) == a(b(e[1]))
            True
            sage: a.inverse()._mul_(c) == b
            True
            sage: c._mul_(b.inverse()) == a
            True
            sage: id = M.identity_map()
            sage: id._mul_(a) == a
            True
            sage: a._mul_(id) == a
            True
            sage: a._mul_(a.inverse()) == id
            True
            sage: a.inverse()._mul_(a) == id
            True

        """
        # No need for consistency check since self and other are guaranted
        # to have the same parent. In particular, they are defined on the same
        # free module.
        #
        # Special cases:
        if self._is_identity:
            return other
        if other._is_identity:
            return self
        if other is self._inverse or self is other._inverse:
            return self._fmodule.identity_map()
        # General case:
        fmodule = self._fmodule
        resu = self.__class__(fmodule)
        basis = self.common_basis(other)
        if basis is None:
            raise ValueError("no common basis for the composition")
        # The composition is performed as a tensor contraction of the last
        # index of self (position=1) and the first index of other (position=0):
        resu._components[basis] = self._components[basis].contract(1,
                                                    other._components[basis],0)
        return resu

    #### End of MultiplicativeGroupElement methods ####

    def __mul__(self, other):
        r"""
        Redefinition of
        :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.__mul__`
        so that * dispatch either to automorphism composition or to the tensor
        product.

        EXAMPLES:

        Automorphism composition::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[1,2],[1,3]])
            sage: b = M.automorphism([[3,4],[5,7]])
            sage: s = a*b ; s
            Automorphism of the Rank-2 free module M over the Integer Ring
            sage: s.matrix()
            [13 18]
            [18 25]
            sage: s.matrix() == a.matrix() * b.matrix()
            True
            sage: s(e[0]) == a(b(e[0]))
            True
            sage: s(e[1]) == a(b(e[1]))
            True
            sage: s.display()
            13 e_0*e^0 + 18 e_0*e^1 + 18 e_1*e^0 + 25 e_1*e^1

        Tensor product::

            sage: c =  M.tensor((1,1)) ; c
            Type-(1,1) tensor on the Rank-2 free module M over the Integer Ring
            sage: c[:] = [[3,4],[5,7]]
            sage: c[:] == b[:]  # c and b have the same components
            True
            sage: s = a*c ; s
            Type-(2,2) tensor on the Rank-2 free module M over the Integer Ring
            sage: s.display()
            3 e_0*e_0*e^0*e^0 + 4 e_0*e_0*e^0*e^1 + 6 e_0*e_0*e^1*e^0
             + 8 e_0*e_0*e^1*e^1 + 5 e_0*e_1*e^0*e^0 + 7 e_0*e_1*e^0*e^1
             + 10 e_0*e_1*e^1*e^0 + 14 e_0*e_1*e^1*e^1 + 3 e_1*e_0*e^0*e^0
             + 4 e_1*e_0*e^0*e^1 + 9 e_1*e_0*e^1*e^0 + 12 e_1*e_0*e^1*e^1
             + 5 e_1*e_1*e^0*e^0 + 7 e_1*e_1*e^0*e^1 + 15 e_1*e_1*e^1*e^0
             + 21 e_1*e_1*e^1*e^1

        """
        if isinstance(other, FreeModuleAutomorphism):
            return self._mul_(other)  # general linear group law
        else:
            return FreeModuleTensor.__mul__(self, other)  # tensor product

    def __imul__(self, other):
        r"""
        Redefinition of
        :meth:`sage.structure.element.ModuleElement.__imul__`

        """
        return self.__mul__(other)


    def matrix(self, basis1=None, basis2=None):
        r"""
        Return the matrix of ``self`` w.r.t to a pair of bases.

        If the matrix is not known already, it is computed from the matrix in
        another pair of bases by means of the change-of-bases formula.

        INPUT:

        - ``basis1`` -- (default: ``None``) basis of the free module on which
          ``self`` is defined ; if none is provided, the module's default basis
            is assumed
        - ``basis2`` -- (default: ``None``) basis of the free module on which
          ``self`` is defined ; if none is provided, ``basis2`` is set to
          ``basis1``.

        OUTPUT:

        - the matrix representing representing the automorphism ``self`` w.r.t
          to bases ``basis1`` and ``basis2``; more precisely, the columns of
          this matrix are formed by the components w.r.t. ``basis2`` of 
          the images of the elements of ``basis1``.

        EXAMPLES:

        Matrices of an automorphism of a rank-3 free `\ZZ`-module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M', start_index=1)
            sage: e = M.basis('e')
            sage: a = M.automorphism([[-1,0,0],[0,1,2],[0,1,3]], name='a')
            sage: a.matrix(e)
            [-1  0  0]
            [ 0  1  2]
            [ 0  1  3]
            sage: a.matrix()
            [-1  0  0]
            [ 0  1  2]
            [ 0  1  3]
            sage: f = M.basis('f', from_family=(-e[2], 4*e[1]+3*e[3],  7*e[1]+5*e[3])) ; f
            Basis (f_1,f_2,f_3) on the Rank-3 free module M over the Integer Ring
            sage: a.matrix(f)
            [  1  -6 -10]
            [ -7  83 140]
            [  4 -48 -81]

        Check of the above matrix::

            sage: a(f[1]).display(f)
            a(f_1) = f_1 - 7 f_2 + 4 f_3
            sage: a(f[2]).display(f)
            a(f_2) = -6 f_1 + 83 f_2 - 48 f_3
            sage: a(f[3]).display(f)
            a(f_3) = -10 f_1 + 140 f_2 - 81 f_3

        Check of the change-of-basis formula::

            sage: p = M.change_of_basis(e,f).matrix(e)
            sage: a.matrix(f) == ~p * a.matrix(e) * p
            True

        Check that the matrix of the product of two automorphisms is the
        product of their matrices::

            sage: b = M.change_of_basis(e,f) ; b
            Automorphism of the Rank-3 free module M over the Integer Ring
            sage: b.matrix(e)
            [ 0  4  7]
            [-1  0  0]
            [ 0  3  5]
            sage: (a*b).matrix(e) == a.matrix(e) * b.matrix(e)
            True

        Check that the matrix of the inverse automorphism is the inverse of the
        automorphism's matrix::

            sage: (~a).matrix(e)
            [-1  0  0]
            [ 0  3 -2]
            [ 0 -1  1]
            sage: (~a).matrix(e) == ~(a.matrix(e))
            True

        Matrices of the identity map::

            sage: id = M.identity_map()
            sage: id.matrix(e)
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: id.matrix(f)
            [1 0 0]
            [0 1 0]
            [0 0 1]

        """
        from sage.matrix.constructor import matrix
        fmodule = self._fmodule
        if basis1 is None:
            basis1 = fmodule.default_basis()
        elif basis1 not in fmodule.bases():
            raise TypeError("{} is not a basis on the {}".format(basis1,
                                                                 fmodule))
        if basis2 is None:
            basis2 = basis1
        elif basis2 not in fmodule.bases():
            raise TypeError("{} is not a basis on the {}".format(basis2,
                                                                 fmodule))
        if (basis1, basis2) not in self._matrices:
            if basis2 == basis1:
                comp = self.components(basis1)
                mat = [[comp[[i,j]] for j in fmodule.irange()] 
                                                     for i in fmodule.irange()]
                self._matrices[(basis1, basis1)] = matrix(mat)
            else:
                # 1/ determine the matrix w.r.t. basis1:
                self.matrix(basis1)
                # 2/ perform the change (basis1, basis1) --> (basis1, basis2):
                raise NotImplementedError("basis1 != basis2 not implemented yet")
        return self._matrices[(basis1, basis2)]

    def det(self):
        r"""
        Return the determinant of ``self``.

        OUTPUT:

        - element of the base ring of the module on which ``self`` is defined,
          equal to the determinant of ``self``.
          
        EXAMPLES:

        Determinant of an automorphism on a `\ZZ`-module of rank 2::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[4,7],[3,5]], name='a')
            sage: a.matrix(e)
            [4 7]
            [3 5]
            sage: a.det()
            -1
            sage: det(a)
            -1
            sage: ~a.det()  # determinant of the inverse automorphism
            -1
            sage: id = M.identity_map()
            sage: id.det()
            1

        """
        self.matrix() # forces the update of the matrix in the module's default
                      # basis, to make sure that the dictionary self._matrices
                      # is not empty
        return self._matrices.values()[0].det() # pick a random value in the
                                                # dictionary self._matrices
                                                # and compute the determinant


    def trace(self):
        r"""
        Return the trace of ``self``.

        OUTPUT:

        - element of the base ring of the module on which ``self`` is defined,
          equal to the trace of ``self``.
          
        EXAMPLES:

        Trace of an automorphism on a `\ZZ`-module of rank 2::

            sage: M = FiniteRankFreeModule(ZZ, 2, name='M')
            sage: e = M.basis('e')
            sage: a = M.automorphism([[4,7],[3,5]], name='a')
            sage: a.matrix(e)
            [4 7]
            [3 5]
            sage: a.trace()
            9
            sage: id = M.identity_map()
            sage: id.trace()
            2

        """
        self.matrix() # forces the update of the matrix in the module's default
                      # basis, to make sure that the dictionary self._matrices
                      # is not empty
        return self._matrices.values()[0].trace() # pick a random value in the
                                                  # dictionary self._matrices
                                                  # and compute the trace
