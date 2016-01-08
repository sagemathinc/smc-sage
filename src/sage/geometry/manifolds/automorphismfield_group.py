r"""
Group of tangent-space automorphism fields

Given an open subset `U` of a manifold `S` and a differentiable mapping
`\Phi: U \rightarrow V`, where `V` is an open subset of a manifold `M`,
the *group of tangent-space automorphism fields* associated with `U` and
`\Phi` is the general linear group `\mathrm{GL}(\mathcal{X}(U,\Phi))` of
the module `\mathcal{X}(U,\Phi)` of vector fields along `U` with values in
`V\subset M` (see
:class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldModule`).
Note that `\mathcal{X}(U,\Phi)` is a module over
`C^\infty(U)`, the algebra of smooth scalar fields on `U`.
Elements of `\mathrm{GL}(\mathcal{X}(U,\Phi))` are fields along `U` of
automorphisms of the tangent spaces to `M`.

Two classes implement `\mathrm{GL}(\mathcal{X}(U,\Phi))` depending whether
`V` is parallelizable or not:
:class:`AutomorphismFieldParalGroup` and :class:`AutomorphismFieldGroup`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

- Chap. 15 of R. Godement: *Algebra*, Hermann (Paris) / Houghton Mifflin
  (Boston) (1968)

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.groups import Groups
from sage.tensor.modules.free_module_linear_group import FreeModuleLinearGroup
from sage.geometry.manifolds.vectorfield_module import VectorFieldModule, \
                                                       VectorFieldFreeModule
from sage.geometry.manifolds.automorphismfield import AutomorphismField, \
                                                      AutomorphismFieldParal

class AutomorphismFieldGroup(UniqueRepresentation, Parent):
    r"""
    General linear group of the module of vector fields along an open subset
    `U` of some manifold `S` with values in an open subset `V` of a
    manifold `M`.

    Given an open subset `U` of a manifold `S` and a differentiable mapping
    `\Phi: U \rightarrow V`, where `V` is an open subset of a manifold `M`,
    the *group of tangent-space automorphism fields* associated with `U` and
    `\Phi` is the general linear group `\mathrm{GL}(\mathcal{X}(U,\Phi))` of
    the module `\mathcal{X}(U,\Phi)` of vector fields along `U` with values in
    `V\subset M` (see
    :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldModule`).
    Note that `\mathcal{X}(U,\Phi)` is a module over
    `C^\infty(U)`, the algebra of smooth scalar fields on `U`.
    Elements of `\mathrm{GL}(\mathcal{X}(U,\Phi))` are fields along `U` of
    automorphisms of the tangent spaces to `M`.

    If `V` is parallelizable, the class
    :class:`AutomorphismFieldParalGroup` must be used instead.

    This is a Sage *parent* class, the corresponding *element* class being
    :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField`.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector fields
      along `U` with values on `V`, as an instance of
      :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldModule`

    EXAMPLES:

    Group of tangent-space automorphism fields of the 2-sphere::

        sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                             intersection_name='W', restrictions1= x^2+y^2!=0, \
                                             restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: G = M.automorphism_field_group() ; G
        General linear group of the module X(M) of vector fields on the 2-dimensional manifold 'M'
        sage: type(G)
        <class 'sage.geometry.manifolds.automorphismfield_group.AutomorphismFieldGroup_with_category'>

    ``G`` is the general linear group of the vector field module
    `\mathcal{X}(M)`::

        sage: XM = M.vector_field_module() ; XM
        module X(M) of vector fields on the 2-dimensional manifold 'M'
        sage: G is XM.general_linear_group()
        True

    ``G`` is a non-abelian group::

        sage: G.category()
        Category of groups
        sage: G in Groups()
        True
        sage: G in CommutativeAdditiveGroups()
        False

    ``G`` is a *parent* object, whose elements are tangent-space automorphisms::

        sage: G.Element
        <class 'sage.geometry.manifolds.automorphismfield.AutomorphismField'>
        sage: a = G.an_element() ; a
        field of tangent-space identity maps on the 2-dimensional manifold 'M'
        sage: a.parent() is G
        True

    The identity element of the group ``G``::

        sage: e = G.one() ; e
        field of tangent-space identity maps on the 2-dimensional manifold 'M'
        sage: eU = U.default_frame() ; eU
        coordinate frame (U, (d/dx,d/dy))
        sage: eV = V.default_frame() ; eV
        coordinate frame (V, (d/du,d/dv))
        sage: e.display(eU)
        Id = d/dx*dx + d/dy*dy
        sage: e.display(eV)
        Id = d/du*du + d/dv*dv

    """

    Element = AutomorphismField

    def __init__(self, vector_field_module):
        r"""
        See :class:`AutomorphismfieldGroup` for documentation and examples.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U') ; V = M.open_subset('V')
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart()
            sage: transf = c_xy.transition_map(c_uv, (x+y, x-y),
            ....:  intersection_name='W', restrictions1= x>0,
            ....:  restrictions2= u+v>0)
            sage: inv = transf.inverse()
            sage: from sage.geometry.manifolds.automorphismfield_group import \
            ....:                                        AutomorphismFieldGroup
            sage: G = AutomorphismFieldGroup(M.vector_field_module()) ; G
            General linear group of the module X(M) of vector fields on the
             2-dimensional manifold 'M'

        """
        if not isinstance(vector_field_module, VectorFieldModule):
            raise TypeError("{} is not a module of vector fields".format(
                            vector_field_module))
        Parent.__init__(self, category=Groups())
        self._vmodule = vector_field_module
        self._one = None # to be set by self.one()


    #### Parent methods ####

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct a field of tangent-space automorphisms.

        OUTPUT:

        - instance of
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`

        """
        if comp == 1:
            return self.one()
        # standard construction
        resu = self.element_class(self._vmodule, name=name,
                                  latex_name=latex_name)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some specific field of tangent-space automorphisms.

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField`

        EXAMPLES:

        """
        return self.one() #!# not terrible...

    #### End of parent methods ####


    #### Monoid methods ####

    def one(self):
        r"""
        Return the group identity element of ``self``.

        The group identity element is the field of tangent-space identity maps.

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField`
          representing the identity element.

        EXAMPLE:

        Identity element of the group of tangent-space automorphism fields of
        the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                                 intersection_name='W', restrictions1= x^2+y^2!=0, \
                                                 restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: G = M.automorphism_field_group()
            sage: G.one()
            field of tangent-space identity maps on the 2-dimensional manifold 'M'
            sage: G.one().restrict(U)[:]
            [1 0]
            [0 1]
            sage: G.one().restrict(V)[:]
            [1 0]
            [0 1]

        """
        if self._one is None:
            self._one = self.element_class(self._vmodule, is_identity=True)
        return self._one

    #### End of monoid methods ####

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLE:
        """
        return "General linear group of the {}".format(self._vmodule)

    def _latex_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLE:

        """
        from sage.misc.latex import latex
        return r"\mathrm{GL}\left("+ latex(self._vmodule)+ r"\right)"


    def base_module(self):
        r"""
        Return the vector-field module of which ``self`` is the general linear
        group.

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldModule`

        EXAMPLE:

        Base module of the group of tangent-space automorphism fields of
        the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                                 intersection_name='W', restrictions1= x^2+y^2!=0, \
                                                 restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: G = M.automorphism_field_group()
            sage: G.base_module()
            module X(M) of vector fields on the 2-dimensional manifold 'M'
            sage: G.base_module() is M.vector_field_module()
            True

        """
        return self._vmodule


#******************************************************************************

class AutomorphismFieldParalGroup(FreeModuleLinearGroup):
    r"""
    General linear group of the module of vector fields along an open subset
    `U` of some manifold `S` with values in a parallelizable open subset `V`
    of a manifold `M`.

    Given an open subset `U` of a manifold `S` and a differentiable mapping
    `\Phi: U \rightarrow V`, where `V` is a parrallelizable open subset of
    a manifold `M`,
    the *group of tangent-space automorphism fields* associated with `U` and
    `\Phi` is the general linear group `\mathrm{GL}(\mathcal{X}(U,\Phi))` of
    the module `\mathcal{X}(U,\Phi)` of vector fields along `U` with values in
    `V\subset M` (see
    :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldFreeModule`).
    Note that `\mathcal{X}(U,\Phi)` is a module over
    `C^\infty(U)`, the algebra of smooth scalar fields on `U`.
    Elements of `\mathrm{GL}(\mathcal{X}(U,\Phi))` are fields along `U` of
    automorphisms of the tangent spaces to `M`.

    If `V` is not parallelizable, the class
    :class:`AutomorphismFieldGroup` must be used instead.

    This is a Sage *parent* class, the corresponding *element* class being
    :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismFieldParal`.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` with values on `V`, as an instance of
      :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldFreeModule`

    EXAMPLES:

    Group of tangent-space automorphism fields of a 2-dimensional
    parallelizable manifold::

        sage: Manifold._clear_cache_() # for doctests only
        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: XM = M.vector_field_module() ; XM
        free module X(M) of vector fields on the 2-dimensional manifold 'M'
        sage: from sage.geometry.manifolds.automorphismfield_group import \
        ....:                                       AutomorphismFieldParalGroup
        sage: G = AutomorphismFieldParalGroup(XM)  ; G
        General linear group of the free module X(M) of vector fields on the
         2-dimensional manifold 'M'
        sage: latex(G)
        \mathrm{GL}\left( \mathcal{X}\left(M\right) \right)
        sage: type(G)
        <class 'sage.geometry.manifolds.automorphismfield_group.AutomorphismFieldParalGroup_with_category'>

    Instead of importing ``AutomorphismFieldParalGroup`` in the global name
    space, it is recommended to use the method
    :meth:`~sage.geometry.manifolds.domain.ManifoldOpenSubset.automorphism_field_group`::

        sage: G = M.automorphism_field_group() ; G
        General linear group of the free module X(M) of vector fields on the
         2-dimensional manifold 'M'

    There is a unique instance of this group::

        sage: G is M.automorphism_field_group()
        True

    ``G`` is nothing but the general linear group of the module
    `\mathcal{X}(M)`::

        sage: G is XM.general_linear_group()
        True

    ``G`` is a group::

        sage: G.category()
        Category of groups
        sage: G in Groups()
        True

    It is not an abelian group::

        sage: G in CommutativeAdditiveGroups()
        False

    ``G`` is a *parent* object, whose elements are tangent-space automorphisms::

        sage: G.Element
        <class 'sage.geometry.manifolds.automorphismfield.AutomorphismFieldParal'>
        sage: a = G.an_element() ; a
        field of tangent-space automorphisms on the 2-dimensional manifold 'M'
        sage: a.parent() is G
        True

    As automorphisms of `\mathcal{X}(M)`, the elements of ``G`` map a vector
    field to a vector field::

        sage: v = XM.an_element() ; v
        vector field on the 2-dimensional manifold 'M'
        sage: v.display()
        2 d/dx + 2 d/dy
        sage: a(v)
        vector field on the 2-dimensional manifold 'M'
        sage: a(v).display()
        2 d/dx - 2 d/dy

    Indeed the matrix of ``a`` with respect to the frame
    `(\partial_x,\partial_y)` is::

        sage: a[X.frame(),:]
        [ 1  0]
        [ 0 -1]

    The elements of ``G`` can also be considered as tensor fields of
    type (1,1)::

        sage: a.tensor_type()
        (1, 1)
        sage: a.tensor_rank()
        2
        sage: a.domain()
        2-dimensional manifold 'M'
        sage: a.display()
        d/dx*dx - d/dy*dy

    The identity element of the group ``G`` is::

        sage: id = G.one() ; id
        field of tangent-space identity maps on the 2-dimensional manifold 'M'
        sage: id*a == a
        True
        sage: a*id == a
        True
        sage: a*a^(-1) == id
        True
        sage: a^(-1)*a == id
        True

    Construction of an element by providing its components w.r.t. to the
    manifold's default frame (frame associated to the coordinates `(x,y)`)::

        sage: b = G([[1+x^2,0], [0,1+y^2]]) ; b
        field of tangent-space automorphisms on the 2-dimensional manifold 'M'
        sage: b.display()
        (x^2 + 1) d/dx*dx + (y^2 + 1) d/dy*dy
        sage: (~b).display()  # the inverse automorphism
        1/(x^2 + 1) d/dx*dx + 1/(y^2 + 1) d/dy*dy

    Check of some group law::

        sage: (a*b)^(-1) == b^(-1) * a^(-1)
        True

    More generally, the full test suite of ``G`` is passed::

        sage: TestSuite(G).run()

    Invertible tensor fields of type (1,1) can be converted to elements of
    ``G``::

        sage: t = M.tensor_field(1, 1, name='t')
        sage: t[:] = [[1+exp(y), x*y], [0, 1+x^2]]
        sage: t1 = G(t) ; t1
        field of tangent-space automorphisms 't' on the 2-dimensional manifold 'M'
        sage: t1 in G
        True
        sage: t1.display()
        t = (e^y + 1) d/dx*dx + x*y d/dx*dy + (x^2 + 1) d/dy*dy
        sage: t1^(-1)
        field of tangent-space automorphisms 't^(-1)' on the 2-dimensional
         manifold 'M'
        sage: (t1^(-1)).display()
        t^(-1) = 1/(e^y + 1) d/dx*dx - x*y/(x^2 + (x^2 + 1)*e^y + 1) d/dx*dy
         + 1/(x^2 + 1) d/dy*dy

    Since any automorphism field can be considered as a tensor field of type
    (1,1) on ``M``, there is a coercion map from ``G`` to the module
    `T^{(1,1)}(M)` of type-(1,1) tensor fields::

        sage: T11 = M.tensor_field_module((1,1)) ; T11
        free module T^(1,1)(M) of type-(1,1) tensors fields on the 2-dimensional manifold 'M'
        sage: T11.has_coerce_map_from(G)
        True

    An explicit call of this coercion map is::

        sage: tt = T11(t1) ; tt
        tensor field 't' of type (1,1) on the 2-dimensional manifold 'M'
        sage: tt == t
        True

    An implicit call of the coercion map is performed to subtract an element of
    ``G`` from an element of `T^{(1,1)}(M)`::

        sage: s = t - t1 ; s
        tensor field 't-t' of type (1,1) on the 2-dimensional manifold 'M'
        sage: s.parent() is T11
        True
        sage: s.display()
        t-t = 0

    as well as for the reverse operation::

        sage: s = t1 - t ; s
        tensor field 't-t' of type (1,1) on the 2-dimensional manifold 'M'
        sage: s.display()
        t-t = 0

    """

    Element = AutomorphismFieldParal

    def __init__(self, vector_field_module):
        r"""
        See :class:`AutomorphismfieldParalGroup` for documentation and
        examples.

        TEST::

            sage: M = Manifold(2, 'M') ; M
            2-dimensional manifold 'M'
            sage: X.<x,y> = M.chart()
            sage: from sage.geometry.manifolds.automorphismfield_group import \
            ....:                                   AutomorphismFieldParalGroup
            sage: G = AutomorphismFieldParalGroup(M.vector_field_module()) ; G
            General linear group of the free module X(M) of vector fields on
             the 2-dimensional manifold 'M'
            sage: TestSuite(G).run()

        """
        if not isinstance(vector_field_module, VectorFieldFreeModule):
            raise TypeError("{} is not a free module of vector fields".format(
                            vector_field_module))
        FreeModuleLinearGroup.__init__(self, vector_field_module)