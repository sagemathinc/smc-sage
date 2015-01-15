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

    EXAMPLES:

    Automorphism of a rank-2 free module (vector space) on `\QQ`::

        sage: M = FiniteRankFreeModule(QQ, 2, name='M')
        sage: a = M.automorphism('A') ; a
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
    def __init__(self, fmodule, name=None, latex_name=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_automorphism import FreeModuleAutomorphism
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = FreeModuleAutomorphism(M, name='a')
            sage: a[e,:] = [[1,0,1],[0,2,0],[0,0,-3]]
            sage: #!# TestSuite(a).run(skip="_test_category") # see below

        In the above test suite, _test_category fails because a is not an
        instance of a.parent().category().element_class.

        """
        FreeModuleTensor.__init__(self, fmodule, (1,1), name=name,
                                  latex_name=latex_name,
                                  parent=fmodule.general_linear_group())
        self._inverse = None    # inverse automorphism not set yet
        # MultiplicativeGroupElement attributes:
        # - none
        # FiniteRankFreeModuleMorphism attributes:
        self._is_identity = False
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

        """
        description = "Automorphism "
        if self._name is not None:
            description += self._name + " "
        description += "of the {}".format(self._fmodule)
        return description

    #### FreeModuleTensor methods ####

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self``.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.automorphism(name='a')
            sage: a._new_instance()
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

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`FreeModuleTensor.__call__` to allow for a single
        argument (module element).

        EXAMPLES:

        Call with a single argument --> return a module element::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: a = M.automorphism(name='a')
            sage: e = M.basis('e')
            sage: a[0,1], a[1,1], a[2,1] = 2, 4, -5
            sage: v = M([2,1,4], name='v')
            sage: s = a.__call__(v) ; s
            Element a(v) of the Rank-3 free module M over the Integer Ring
            sage: s.display()
            a(v) = 2 e_0 + 4 e_1 - 5 e_2
            sage: s == a(v)
            True
            sage: s == a.contract(v)
            True

        Call with two arguments (:class:`FreeModuleTensor` behaviour)
        --> return a scalar::

            sage: b = M.linear_form(name='b')
            sage: b[:] = 7, 0, 2
            sage: a.__call__(b,v)
            4
            sage: a(b,v) == a.__call__(b,v)
            True
            sage: a(b,v) == s(b)
            True

        """
        from free_module_tensor import FiniteRankFreeModuleElement
        if len(arg) > 1:
            # the endomorphism acting as a type-(1,1) tensor on a pair
            # (linear form, module element), returning a scalar:
            return FreeModuleTensor.__call__(self, *arg)
        # the endomorphism acting as such, on a module element, returning a
        # module element:
        vector = arg[0]
        if not isinstance(vector, FiniteRankFreeModuleElement):
            raise TypeError("the argument must be an element of a free module")
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
            sage: a = M.automorphism('A')
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

        """
        from sage.matrix.constructor import matrix
        from comp import Components
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

    #### End of MultiplicativeGroupElement methods ####

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
                return FiniteRankFreeModuleMorphism.matrix(self, basis1,
                                                           basis2)
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
            sage: a = M.automorphism('a')
            sage: a[:] = [[4,7],[3,5]]
            sage: a.matrix(e)
            [4 7]
            [3 5]
            sage: a.det()
            -1
            sage: ~a.det()  # determinant of the inverse automorphism
            -1

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
            sage: a = M.automorphism('a')
            sage: a[:] = [[4,7],[3,5]]
            sage: a.matrix(e)
            [4 7]
            [3 5]
            sage: a.trace()
            9
        
        """
        self.matrix() # forces the update of the matrix in the module's default
                      # basis, to make sure that the dictionary self._matrices
                      # is not empty
        return self._matrices.values()[0].trace() # pick a random value in the
                                                  # dictionary self._matrices
                                                  # and compute the trace

        
#******************************************************************************        

class FreeModuleIdentityTensor(FreeModuleAutomorphism):
    r"""
    Identity map (considered as a type-(1,1) tensor) on a free module.

    INPUT:

    - ``fmodule`` -- free module `M` over a commutative ring `R`
      (must be an instance of :class:`FiniteRankFreeModule`)
    - ``name`` -- (default: 'Id') name given to the identity tensor.
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the identity
      tensor; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    Identity tensor on a rank-3 free module::

        sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
        sage: e = M.basis('e')
        sage: a = M.identity_tensor() ; a
        Identity tensor on the Rank-3 free module M over the Integer Ring

    The LaTeX symbol is set by default to `\mathrm{Id}`, but can be changed::

        sage: latex(a)
        \mathrm{Id}
        sage: a = M.identity_tensor(latex_name=r'\mathrm{1}')
        sage: latex(a)
        \mathrm{1}

    The identity is a tensor of type `(1,1)` on the free module::

        sage: a.parent()
        General linear group of the Rank-3 free module M over the Integer Ring
        sage: a.tensor_type()
        (1, 1)
        sage: a.tensor_rank()
        2

    Its components are Kronecker deltas in any basis::

        sage: a[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: a.comp() # components in the module's default basis (e)
        Kronecker delta of size 3x3
        sage: a.display()
        Id = e_0*e^0 + e_1*e^1 + e_2*e^2
        sage: f = M.basis('f')
        sage: a.comp(basis=f)
        Kronecker delta of size 3x3
        sage: a.comp(f)[:]
        [1 0 0]
        [0 1 0]
        [0 0 1]

    The components can be read, but cannot be set::

        sage: a[1,1]
        1
        sage: a[1,1] = 2
        Traceback (most recent call last):
        ...
        TypeError: the components of the identity map cannot be changed

    The identity tensor acting on a module element::

        sage: v = M([2,-3,1], basis=e, name='v')
        sage: v.display()
        v = 2 e_0 - 3 e_1 + e_2
        sage: u = a(v) ; u
        Element v of the Rank-3 free module M over the Integer Ring
        sage: u is v
        True

    The identity tensor acting as a type-`(1,1)` tensor on a pair (linear form,
    module element)::

        sage: w = M.tensor((0,1), name='w') ; w
        Linear form w on the Rank-3 free module M over the Integer Ring
        sage: w[:] = [0, 3, 2]
        sage: s = a(w,v) ; s
        -7
        sage: s == w(v)
        True

    The identity tensor is its own inverse::

        sage: a.inverse() == a
        True
        sage: a.inverse() is a
        True

    """
    def __init__(self, fmodule, name='Id', latex_name=None):
        r"""
        TESTS::

            sage: from sage.tensor.modules.free_module_automorphism import FreeModuleIdentityTensor
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: Id = FreeModuleIdentityTensor(M)
            sage: TestSuite(Id).run(skip="_test_category") # see below

        In the above test suite, _test_category fails because Id is not an
        instance of Id.parent().category().element_class.

        """
        if latex_name is None and name == 'Id':
            latex_name = r'\mathrm{Id}'
        FreeModuleAutomorphism.__init__(self, fmodule, name=name,
                                              latex_name=latex_name)
        self._inverse = self    # the identity is its own inverse
        self.comp() # Initializing the components in the module's default basis

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: M.identity_tensor()
            Identity tensor on the Rank-3 free module M over the Integer Ring

        """
        description = "Identity tensor "
        if self._name != 'Id':
            description += self._name + " "
        description += "on the " + str(self._fmodule)
        return description

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: id = M.identity_tensor()
            sage: id._del_derived()

        """
        # FreeModuleAutomorphism._del_derived is bypassed:
        FreeModuleTensor._del_derived(self)

    def _new_comp(self, basis):
        r"""
        Create some components in the given basis.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: id = M.identity_tensor()
            sage: id._new_comp(e)
            Kronecker delta of size 3x3
            sage: type(id._new_comp(e))
            <class 'sage.tensor.modules.comp.KroneckerDelta'>

        """
        from comp import KroneckerDelta
        fmodule = self._fmodule  # the base free module
        return KroneckerDelta(fmodule._ring, basis, start_index=fmodule._sindex,
                              output_formatter=fmodule._output_formatter)

    def components(self, basis=None, from_basis=None):
        r"""
        Return the components in a given basis as a Kronecker delta.

        INPUT:

        - ``basis`` -- (default: ``None``) module basis in which the components
          are required; if none is provided, the components are assumed to
          refer to the module's default basis
        - ``from_basis`` -- (default: ``None``) unused (present just for
          ensuring compatibility with ``FreeModuleTensor.comp`` calling list)

        OUTPUT:

        - components in the basis ``basis``, as an instance of the
          class :class:`~sage.tensor.modules.comp.KroneckerDelta`

        EXAMPLES:

        Components of the identity map on a rank-3 free module::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.identity_tensor()
            sage: a.components(basis=e)
            Kronecker delta of size 3x3

        For the module's default basis, the argument ``basis`` can be omitted::

            sage: a.components() is a.components(basis=e)
            True
            sage: a.components()[:]
            [1 0 0]
            [0 1 0]
            [0 0 1]

        A shortcut is ``a.comp()``::

            sage: a.comp() is a.components()
            True

        """
        if basis is None:
            basis = self._fmodule._def_basis
        if basis not in self._components:
            self._components[basis] = self._new_comp(basis)
        return self._components[basis]

    comp = components

    def set_comp(self, basis=None):
        r"""
        Redefinition of the generic tensor method
        :meth:`FreeModuleTensor.set_comp`: should not be called.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.identity_tensor()
            sage: a.set_comp(e)
            Traceback (most recent call last):
            ...
            TypeError: the components of the identity map cannot be changed

        """
        raise TypeError("the components of the identity map cannot be changed")

    def add_comp(self, basis=None):
        r"""
        Redefinition of the generic tensor method
        :meth:`FreeModuleTensor.add_comp`: should not be called.

        EXAMPLE::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.identity_tensor()
            sage: a.add_comp(e)
            Traceback (most recent call last):
            ...
            TypeError: the components of the identity map cannot be changed

        """
        raise TypeError("the components of the identity map cannot be changed")

    def __call__(self, *arg):
        r"""
        Redefinition of :meth:`FreeModuleEndomorphismTensor.__call__`.

        EXAMPLES:

        Call with a single argument --> return a module element::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: id = M.identity_tensor()
            sage: v = M([-1,4,3])
            sage: s = id.__call__(v) ; s
            Element of the Rank-3 free module M over the Integer Ring
            sage: s == v
            True
            sage: s == id(v)
            True
            sage: s == id.contract(v)
            True

        Call with two arguments (:class:`FreeModuleTensor` behaviour) -->
        return a scalar::

            sage: b = M.linear_form(name='b')
            sage: b[:] = 7, 0, 2
            sage: id.__call__(b,v)
            -1
            sage: id(b,v) == id.__call__(b,v)
            True
            sage: id(b,v) == b(v)
            True

        """
        from free_module_tensor import FiniteRankFreeModuleElement
        from free_module_alt_form import FreeModuleAltForm
        if len(arg) == 1:
            # the identity map acting as such, on a module element:
            vector = arg[0]
            if not isinstance(vector, FiniteRankFreeModuleElement):
                raise TypeError("the argument must be a module element")
            return vector
            #!# should it be return vector.copy() instead ?
        elif len(arg) == 2:
            # the identity map acting as a type-(1,1) tensor on a pair
            # (1-form, vector), returning a scalar:
            linform = arg[0]
            if not isinstance(linform, FreeModuleAltForm):
                raise TypeError("the first argument must be a linear form")
            if linform._tensor_type != (0,1):
                raise TypeError("the first argument must be a linear form")                
            vector = arg[1]
            if not isinstance(vector, FiniteRankFreeModuleElement):
                raise TypeError("the second argument must be a module element")
            return linform(vector)
        else:
            raise TypeError("wrong number of arguments")

