r"""
Vector field modules

The set of vector fields along an open subset `U` of some manifold `S`
with values in a open subset `V` of a manifold `M` (possibly `S=M` and `U=V`)
is a module over the algebra `C^\infty(U)` of differentiable scalar fields
on `U`. It is a free module iff `V` is parallelizable.
Accordingly, two classes are devoted to vector field modules:

- :class:`VectorFieldModule` for vector fields with values in a
  generic (in practice, not parallelizable) open set `V`
- :class:`VectorFieldFreeModule` for vector fields with values in a
  parallelizable open set `V`

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version

REFERENCES:

- S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*, vol. 1,
  Interscience Publishers (New York) (1963)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)
- B O'Neill : *Semi-Riemannian Geometry*, Academic Press (San Diego) (1983)

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.modules import Modules
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from vectorfield import VectorField, VectorFieldParal

class VectorFieldModule(UniqueRepresentation, Parent):
    r"""
    Module of vector fields along an open subset `U` of some manifold `S`
    with values in a open subset `V` of a manifold `M`.

    If `V` is parallelizable, the class :class:`VectorFieldFreeModule` should
    be used instead.

    Given a differentiable mapping

    .. MATH::

        \Phi:\ U\subset S \longrightarrow V\subset M

    the module `\mathcal{X}(U,\Phi)` is the set of all vector fields of
    the type

    .. MATH::

        v:\ U  \longrightarrow TM

    such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M


    The set `\mathcal{X}(U,\Phi)` is a module over `C^\infty(U)`, the ring
    (algebra) of differentiable scalar fields on `U` (see
    :class:`~sage.geometry.manifolds.scalarfield_algebra.ScalarFieldAlgebra`).

    The standard case of vector fields *on* a manifold corresponds to `S=M`,
    `U=V` and `\Phi = \mathrm{Id}_U`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `V` (`U` is then an open interval
    of `\RR`).

    This is a Sage *parent* class, the corresponding *element* class being
    :class:`~sage.geometry.manifolds.vectorfield.VectorField`.

    INPUT:

    - ``domain`` -- open subset `U` on which the vector fields are defined
    - ``dest_map`` -- (default: ``None``) destination map `\Phi:\ U \rightarrow V`
      (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`);
      if none is provided, the identity is assumed (case of vector fields *on*
      `U`)

    EXAMPLE:

    Module of vector fields on the 2-sphere::

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
        sage: XM = M.vector_field_module() ; XM
        module X(M) of vector fields on the 2-dimensional manifold 'M'

    `\mathcal{X}(M)` is a module over the algebra `C^\infty(M)`::

        sage: XM.category()
        Category of modules over algebra of scalar fields on the 2-dimensional manifold 'M'
        sage: XM.base_ring() is M.scalar_field_algebra()
        True

    `\mathcal{X}(M)` is not a free module::

        sage: isinstance(XM, FiniteRankFreeModule)
        False

    because `M = S^2` is not parallelizable::

        sage: M.is_manifestly_parallelizable()
        False

    On the contrary, the module of vector fields on `U` is a free module,
    since `U` is parallelizable (being a coordinate domain)::

        sage: XU = U.vector_field_module()
        sage: isinstance(XU, FiniteRankFreeModule)
        True
        sage: U.is_manifestly_parallelizable()
        True

    The zero element of the module::

        sage: z = XM.zero() ; z
        vector field 'zero' on the 2-dimensional manifold 'M'
        sage: z.display(c_xy.frame())
        zero = 0
        sage: z.display(c_uv.frame())
        zero = 0

    The module `\mathcal{X}(M)` coerces to any module of vector fields defined
    on a subdomain of `M`, for instance `\mathcal{X}(U)`::

        sage: XU.has_coerce_map_from(XM)
        True
        sage: XU.coerce_map_from(XM)
        Conversion map:
          From: module X(M) of vector fields on the 2-dimensional manifold 'M'
          To:   free module X(U) of vector fields on the open subset 'U' of the 2-dimensional manifold 'M'

    The conversion map is actually the restriction of vector fields defined
    on `M` to `U`.

    Sage test suite for modules is passed::

        sage: TestSuite(XM).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_nonzero_equal() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass

    """

    Element = VectorField

    def __init__(self, domain, dest_map=None):
        self._domain = domain
        name = "X(" + domain._name
        latex_name = r"\mathcal{X}\left(" + domain._latex_name
        if dest_map is None:
            dest_map = domain._identity_map
        self._dest_map = dest_map
        if dest_map is domain._identity_map:
            name += ")"
            latex_name += r"\right)"
        else:
            name += "," + self._dest_map._name + ")"
            latex_name += "," + self._dest_map._latex_name + r"\right)"
        self._ambient_domain = self._dest_map._codomain
        self._name = name
        self._latex_name = latex_name
        # the member self._ring is created for efficiency (to avoid calls to
        # self.base_ring()):
        self._ring = domain.scalar_field_algebra()
        Parent.__init__(self, base=self._ring, category=Modules(self._ring))
        # Dictionary of the tensor modules built on self
        #   (dict. keys = (k,l) --the tensor type)
        self._tensor_modules = {(1,0): self} # self is considered as the set of
                                            # tensors of type (1,0)
        # Dictionary of exterior powers of the dual of self
        #   (keys = p --the power degree) :
        self._dual_exterior_powers = {}
        # Zero element:
        if not hasattr(self, '_zero_element'):
            self._zero_element = self._element_constructor_(name='zero',
                                                            latex_name='0')
            for frame in self._domain._frames:
                if self._dest_map.restrict(frame._domain) == frame._dest_map:
                    self._zero_element.add_comp(frame)
                    # (since new components are initialized to zero)
        # Identity automorphism:
        self._identity_map = None # to be set by self.identity_map()
        # General linear group:
        self._general_linear_group = None # to be set by
                                          # self.general_linear_group()

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct an element of the module
        """
        if comp == 0:
            return self._zero_element
        if isinstance(comp, VectorField):
            if self._domain.is_subset(comp._domain) and \
                       self._ambient_domain.is_subset(comp._ambient_domain):
                return comp.restrict(self._domain)
            else:
                raise TypeError("Cannot coerce the " + str(comp) +
                                "to a vector field in " + str(self))
        resu = self.element_class(self, name=name, latex_name=latex_name)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) element of the module
        """
        resu = self.element_class(self)
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent
        """
        if isinstance(other, (VectorFieldModule, VectorFieldFreeModule)):
            return self._domain.is_subset(other._domain) and \
                   self._ambient_domain.is_subset(other._ambient_domain)
        else:
            return False

    #### End of parent methods

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "module "
        if self._name is not None:
            description += self._name + " "
        description += "of vector fields "
        if self._dest_map is self._domain._identity_map:
            description += "on the " + str(self._domain)
        else:
            description += "along the " + str(self._domain) + \
                           " mapped into the " + str(self._ambient_domain)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def tensor_module(self, k, l):
        r"""
        Return the module of all tensor fields of type (k,l) defined on
        ``self``.

        INPUT:

        - ``k`` -- (non-negative integer) the contravariant rank, the tensor type
          being (k,l)
        - ``l`` -- (non-negative integer) the covariant rank, the tensor type
          being (k,l)

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.tensor_field_module.TensorFieldModule`
          representing the free module
          `T^{(k,l)}(M)` of type-`(k,l)` tensors on the free module ``self``.

        EXAMPLES:
        """
        from tensorfield_module import TensorFieldModule
        if (k,l) not in self._tensor_modules:
            self._tensor_modules[(k,l)] = TensorFieldModule(self, (k,l))
        return self._tensor_modules[(k,l)]

    def dual_exterior_power(self, p):
        r"""
        Return the `p`-th exterior power of the dual of ``self``.

        If ``self`` is the vector field module `\mathcal{X}(U,\Phi)`, the
        `p`-th exterior power of its dual is the set `\Lambda^p(U,\Phi)` of
        `p`-forms along `U` with values in `\Phi(U)`. It is a module over
        `C^\infty(U)`, the ring (algebra) of differentiable scalar fields on
        `U`.

        INPUT:

        - ``p`` -- non-negative integer

        OUTPUT:

        - for `p\geq 1`, instance of
          :class:`~sage.geometry.manifolds.diffform_module.DiffFormModule`
          representing the module `\Lambda^p(U,\Phi)`; for `p=0`, the
          base ring, i.e. `C^\infty(U)`, is returned instead

        EXAMPLES:

        """
        from sage.geometry.manifolds.diffform_module import DiffFormModule
        if p == 0:
            return self._ring
        if p not in self._dual_exterior_powers:
            self._dual_exterior_powers[p] = DiffFormModule(self, p)
        return self._dual_exterior_powers[p]

    def dual(self):
        r"""
        Return the dual module.

        EXAMPLE:

        """
        return self.dual_exterior_power(1)

    def general_linear_group(self):
        r"""
        Return the general linear group of ``self``.

        If ``self`` is the module `\mathcal{X}(U,\Phi)`, the *general
        linear group* is the group `\mathrm{GL}(\mathcal{X}(U,\Phi))` of
        automorphisms of `\mathcal{X}(U,\Phi)`. Note that an automorphism of
        `\mathcal{X}(U,\Phi)` can also be viewed as a *field* along `U` of
        automorphisms of the tangent spaces of `V=\Phi(U)`.

        OUTPUT:

        - instance of class
          :class:`~sage.geometry.manifolds.automorphismfield_group.AutomorphismFieldGroup`
          representing `\mathrm{GL}(\mathcal{X}(U,\Phi))`

        EXAMPLES:

        """
        from sage.geometry.manifolds.automorphismfield_group import \
                                                         AutomorphismFieldGroup
        if self._general_linear_group is None:
            self._general_linear_group = AutomorphismFieldGroup(self)
        return self._general_linear_group

    def tensor(self, tensor_type, name=None, latex_name=None, sym=None,
               antisym=None, specific_type=None):
        r"""
        Construct a tensor on the vector field module.

        The tensor is actually a tensor field on the domain of ``self``.

        INPUT:

        - ``tensor_type`` -- pair (k,l) with k being the contravariant rank
          and l the covariant rank
        - ``name`` -- (string; default: ``None``) name given to the tensor
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          tensor; if none is provided, the LaTeX symbol is set to ``name``
        - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
          the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention position=0 for the first argument. For instance:

          * sym=(0,1) for a symmetry between the 1st and 2nd arguments
          * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments.

        - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
          among the arguments, with the same convention as for ``sym``.
        - ``specific_type`` -- (default: ``None``) specific subclass of
          :class:`~sage.geometry.manifolds.tensorfield.TensorField` for the
          output

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.tensorfield.TensorField`
          representing the tensor defined on ``self`` with the provided
          characteristics.

        EXAMPLES:

        """
        from automorphismfield import AutomorphismField
        from metric import Metric, RiemannMetric, LorentzMetric
        if tensor_type==(1,0):
            return self.element_class(self, name=name, latex_name=latex_name)
        elif tensor_type==(0,1):
            return self.linear_form(name=name, latex_name=latex_name)
        elif tensor_type==(1,1) and specific_type is not None:
            if issubclass(specific_type, AutomorphismField):
                return self.automorphism(name=name, latex_name=latex_name)
        elif tensor_type[0]==0 and tensor_type[1]>1 and antisym is not None \
                                                              and antisym !=[]:
            if isinstance(antisym, list):
                antisym0 = antisym[0]
            else:
                antisym0 = antisym
            if len(antisym0)==tensor_type[1]:
                return self.alternating_form(tensor_type[1], name=name,
                                             latex_name=latex_name)
            else:
                return self.tensor_module(*tensor_type).element_class(self,
                                 tensor_type, name=name, latex_name=latex_name,
                                 sym=sym, antisym=antisym)
        elif tensor_type==(0,2):
            if specific_type == Metric:
                return Metric(self, name, latex_name=latex_name)
                # NB: the signature is not passed
            elif specific_type == RiemannMetric:
                return RiemannMetric(self, name, latex_name=latex_name)
            elif specific_type == LorentzMetric:
                return LorentzMetric(self, name, latex_name=latex_name)
                # NB: the signature convention is not passed
            else:
                return self.tensor_module(0,2).element_class(self, (0,2),
                                              name=name, latex_name=latex_name,
                                              sym=sym, antisym=antisym)
        # Generic case
        return self.tensor_module(*tensor_type).element_class(self,
         tensor_type, name=name, latex_name=latex_name, sym=sym,
         antisym=antisym)

    def alternating_form(self, degree, name=None, latex_name=None):
        r"""
        Construct an alternating form on the module ``self``.

        INPUT:

        - ``degree`` -- the degree of the alternating form (i.e. its tensor rank)
        - ``name`` -- (string; default: ``None``) name given to the alternating
          form
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          alternating form; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.diffform.DiffForm`

        See
        :class:`~sage.geometry.manifolds.diffform.DiffForm`
        for further documentation.

        """
        return self.dual_exterior_power(degree).element_class(self, degree,
                                              name=name, latex_name=latex_name)

    def linear_form(self, name=None, latex_name=None):
        r"""
        Construct a linear form on the module ``self``.

        A linear form on the vector field module ``self`` is actually a field
        of linear forms (i.e. a 1-form) along the open subset `U` on which
        ``self`` is defined.

        INPUT:

        - ``name`` -- (string; default: ``None``) name given to the linear
          form
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          linear form; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.diffform.DiffForm`

        See
        :class:`~sage.geometry.manifolds.diffform.DiffForm`
        for further documentation.


        """
        return self.dual_exterior_power(1).element_class(self, 1, name=name,
                                                         latex_name=latex_name)

    def automorphism(self, name=None, latex_name=None):
        r"""
        Construct an automorphism of the module ``self``.

        An automorphism of the module ``self`` is actually a field
        of tangent-space automorphisms along the open subset `U` on which
        ``self`` is defined.

        INPUT:

        - ``name`` -- (string; default: ``None``) name given to the automorphism
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          automorphism; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField`

        See
        :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField`
        for more documentation.

        """
        return self.general_linear_group().element_class(self, name=name,
                                                         latex_name=latex_name)

    def identity_map(self, name='Id', latex_name=None):
        r"""
        Construct the identity map on the module ``self``.

        The identity map on the module ``self`` is actually a field
        of tangent-space identity maps along the open subset `U` on which
        ``self`` is defined.

        INPUT:

        - ``name`` -- (string; default: 'Id') name given to the identity map
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          identity map;  if none is provided, the LaTeX symbol is set to
          '\mathrm{Id}' if ``name`` is 'Id' and to ``name`` otherwise

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField`

        """
        if self._identity_map is None:
            self._identity_map = self.general_linear_group().one()
            if name != 'Id':
                if latex_name is None:
                    latex_name = name
                self._identity_map.set_name(name=name, latex_name=latex_name)
        return self._identity_map

    def metric(self, name, signature=None, latex_name=None):
        r"""
        Construct a pseudo-Riemannian metric (nondegenerate symmetric bilinear
        form) on the free module ``self``.

        A metric of the vector free module ``self`` is actually a field
        of tangent-space metrics along the open subset `U` on which
        ``self`` is defined.

        INPUT:

        - ``name`` -- (string) name given to the metric
        - ``signature`` -- (integer; default: ``None``) signature `S` of the
          metric: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the number of
          positive terms (resp. number of negative terms) in any diagonal writing
          of the metric components; if ``signature`` is not provided, `S` is set to
          the manifold's dimension (Riemannian signature)
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          metric; if none, it is formed from ``name``

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.metric.Metric`
          representing the defined pseudo-Riemannian metric.

        See :class:`~sage.geometry.manifolds.metric.Metric` for further
        documentation.

        """
        from metric import Metric
        return Metric(self, name, signature=signature, latex_name=latex_name)

    def riemann_metric(self, name, latex_name=None):
        r"""
        Construct a Riemannian metric (positive definite symmetric bilinear
        form) on the free module ``self``.

        A Riemannian metric of the vector free module ``self`` is actually a
        field of tangent-space Riemannian metrics along the open subset `U`
        on which ``self`` is defined.

        INPUT:

        - ``name`` -- (string) name given to the metric
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          metric; if none, it is formed from ``name``

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.metric.RiemannMetric`
          representing the defined Riemannian metric.

        See :class:`~sage.geometry.manifolds.metric.RiemannMetric` for
        further documentation.

        """
        from metric import RiemannMetric
        return RiemannMetric(self, name, latex_name=latex_name)

    def lorentz_metric(self, name, signature='positive', latex_name=None):
        r"""
        Construct a Lorentzian metric (symmetric bilinear
        form of signature (-,+,...,+) or (+,-,...,-)) on the free
        module ``self``.

        A Lorentzian metric of the vector free module ``self`` is actually a
        field of tangent-space Lorentzian metrics along the open subset `U` on
        which ``self`` is defined.

        INPUT:

        - ``name`` -- (string) name given to the metric
        - ``signature`` -- (string, default: 'positive') sign of the metric
          signature:

          * if set to 'positive', the signature is n-2, where n is the manifold's
            dimension, i.e. `(-,+,\cdots,+)`
          * if set to 'negative', the signature is -n+2, i.e. `(+,-,\cdots,-)`

        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          metric; if none, it is formed from ``name``

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.metric.LorentzMetric`
          representing the defined Loretnzian metric.

        See :class:`~sage.geometry.manifolds.metric.LorentzMetric` for
        further documentation.

        """
        from metric import LorentzMetric
        return LorentzMetric(self, name, signature=signature,
                             latex_name=latex_name)


#******************************************************************************

class VectorFieldFreeModule(FiniteRankFreeModule):
    r"""
    Free module of vector fields along an open subset `U` of some manifold `S`
    with values in a parallelizable open subset `V` of a manifold `M`.

    Given a differentiable mapping

    .. MATH::

        \Phi:\ U\subset S \longrightarrow V\subset M

    the module `\mathcal{X}(U,\Phi)` is the set of all vector fields of
    the type

    .. MATH::

        v:\ U  \longrightarrow TM

    such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M


    Since `V` is parallelizable, the set `\mathcal{X}(U,\Phi)` is a free module
    over `C^\infty(U)`, the ring (algebra) of differentiable scalar fields on
    `U` (see
    :class:`~sage.geometry.manifolds.scalarfield_algebra.ScalarFieldAlgebra`).
    Its rank is the dimension of `M`.

    The standard case of vector fields *on* a manifold corresponds to `S=M`,
    `U=V` and `\Phi = \mathrm{Id}_U`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `V` (`U` is then an open interval
    of `\RR`).

    This is a Sage *parent* class, the corresponding *element* class being
    :class:`~sage.geometry.manifolds.vectorfield.VectorFieldParal`.

    INPUT:

    - ``domain`` -- open subset `U` on which the vector fields are defined
    - ``dest_map`` -- (default: ``None``) destination map `\Phi:\ U \rightarrow V`
      (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`);
      if none is provided, the identity is assumed (case of vector fields *on*
      `U`)

    EXAMPLES:

    Module of vector fields on `\RR^2`::

        sage: M = Manifold(2, 'R^2')
        sage: cart.<x,y> = M.chart()  # Cartesian coordinates on R^2
        sage: XM = M.vector_field_module() ; XM
        free module X(R^2) of vector fields on the 2-dimensional manifold 'R^2'
        sage: XM.category()
        Category of modules over algebra of scalar fields on the 2-dimensional manifold 'R^2'
        sage: XM.base_ring() is M.scalar_field_algebra()
        True

    Since `\RR^2` is obviously parallelizable, ``XM`` is a free module::

        sage: isinstance(XM, FiniteRankFreeModule)
        True

    Some elements::

        sage: XM.an_element().display()
        2 d/dx + 2 d/dy
        sage: XM.zero().display()
        zero = 0
        sage: v = XM([-y,x]) ; v
        vector field on the 2-dimensional manifold 'R^2'
        sage: v.display()
        -y d/dx + x d/dy

    An example of module of vector fields with a destination map `\Phi`
    different from the identity map, namely a mapping
    `\Phi: I \rightarrow \RR^2`, where `I` is an open interval of `\RR`::

        sage: R.<t> = RealLine()
        sage: I = R.open_interval(0, 2*pi)
        sage: Phi = I.diff_mapping(M, coord_functions=[cos(t), sin(t)], name='Phi',
        ....:                      latex_name=r'\Phi') ; Phi
        Curve 'Phi' in the 2-dimensional manifold 'R^2'
        sage: Phi.display()
        Phi: (0, 2*pi) --> R^2
           t |--> (x, y) = (cos(t), sin(t))
        sage: XIM = I.vector_field_module(dest_map=Phi) ; XIM
        free module X((0, 2*pi),Phi) of vector fields along the Real interval
         (0, 2*pi) mapped into the 2-dimensional manifold 'R^2'
        sage: XIM.category()
        Category of modules over algebra of scalar fields on the Real interval
         (0, 2*pi)

    The rank of the free module `\mathcal{X}((0, 2\pi),\Phi)` is the dimension
    of the manifold `\RR^2`, namely two::

        sage: XIM.rank()
        2

    A basis of it is induced by the coordinate vector frame of `\RR^2`::

        sage: XIM.bases()
        [vector frame ((0, 2*pi), (d/dx,d/dy)) with values on the 2-dimensional
         manifold 'R^2']

    Some elements of this module::

        sage: XIM.an_element().display()
        2 d/dx + 2 d/dy
        sage: v = XIM([t, t^2]) ; v
        vector field along the Real interval (0, 2*pi) with values on the
         2-dimensional manifold 'R^2'
        sage: v.display()
        t d/dx + t^2 d/dy

    The test suite is passed::

        sage: TestSuite(XIM).run()

    Let us now consider the module of vector fields on the circle `S^1`; we
    start by constructing the `S^1` manifold::

        sage: M = Manifold(1, 'S^1')
        sage: U = M.open_subset('U')  # the complement of one point
        sage: c_t.<t> =  U.chart('t:(0,2*pi)') # the standard angle coordinate
        sage: V = M.open_subset('V') # the complement of the point t=pi
        sage: M.declare_union(U,V)   # S^1 is the union of U and V
        sage: c_u.<u> = V.chart('u:(0,2*pi)') # the angle t-pi
        sage: t_to_u = c_t.transition_map(c_u, (t-pi,), intersection_name='W', restrictions1 = t!=pi, restrictions2 = u!=pi)
        sage: u_to_t = t_to_u.inverse()
        sage: W = U.intersection(V)

    `S^1` cannot be covered by a single chart, so it cannot be covered by
    a coordinate frame. It is however parallelizable and we introduce a global
    vector frame as follows. We notice that on their common subdomain, `W`,
    the coordinate vectors `\partial/\partial t` and `\partial/\partial u`
    coincide, as we can check explicitely::

        sage: c_t.frame()[0].display(c_u.frame().restrict(W))
        d/dt = d/du

    Therefore, we can extend `\partial/\partial t` to all `V` and hence to all
    `S^1`, to form a vector field on `S^1` whose components w.r.t. both
    `\partial/\partial t` and `\partial/\partial u` are 1::

        sage: e = M.vector_frame('e')
        sage: U.set_frame_change(e.restrict(U), c_t.frame(), U.tangent_identity_field())
        sage: V.set_frame_change(e.restrict(V), c_u.frame(), V.tangent_identity_field())
        sage: e[0].display(c_t.frame())
        e_0 = d/dt
        sage: e[0].display(c_u.frame())
        e_0 = d/du

    Equipped with the frame `e`, the manifold `S^1` is manifestly
    parallelizable::

        sage: M.is_manifestly_parallelizable()
        True

    Consequently, the module of vector fields on `S^1` is a free module::

        sage: XM = M.vector_field_module() ; XM
        free module X(S^1) of vector fields on the 1-dimensional manifold 'S^1'
        sage: isinstance(XM, FiniteRankFreeModule)
        True
        sage: XM.category()
        Category of modules over algebra of scalar fields on the 1-dimensional manifold 'S^1'
        sage: XM.base_ring() is M.scalar_field_algebra()
        True

    The zero element::

        sage: z = XM.zero() ; z
        vector field 'zero' on the 1-dimensional manifold 'S^1'
        sage: z.display()
        zero = 0
        sage: z.display(c_t.frame())
        zero = 0

    The module `\mathcal{X}(S^1)` coerces to any module of vector fields
    defined on a subdomain of `S^1`, for instance `\mathcal{X}(U)`::

        sage: XU = U.vector_field_module() ; XU
        free module X(U) of vector fields on the open subset 'U' of the 1-dimensional manifold 'S^1'
        sage: XU.has_coerce_map_from(XM)
        True
        sage: XU.coerce_map_from(XM)
        Conversion map:
          From: free module X(S^1) of vector fields on the 1-dimensional manifold 'S^1'
          To:   free module X(U) of vector fields on the open subset 'U' of the 1-dimensional manifold 'S^1'

    The conversion map is actually the restriction of vector fields defined
    on `S^1` to `U`.

    Sage test suite for modules is passed::

        sage: TestSuite(XM).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_nonzero_equal() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass

    """

    Element = VectorFieldParal

    def __init__(self, domain, dest_map=None):
        from scalarfield import ScalarField
        self._domain = domain
        if dest_map is None:
            dest_map = domain._identity_map
        self._dest_map = dest_map
        self._ambient_domain = self._dest_map._codomain
        name = "X(" + domain._name
        latex_name = r"\mathcal{X}\left(" + domain._latex_name
        if self._dest_map == domain._identity_map:
            name += ")"
            latex_name += r"\right)"
        else:
            name += "," + self._dest_map._name + ")"
            latex_name += "," + self._dest_map._latex_name + r"\right)"
        manif = self._ambient_domain._manifold
        FiniteRankFreeModule.__init__(self, domain.scalar_field_algebra(),
                                  manif._dim, name=name, latex_name=latex_name,
                                  start_index=manif._sindex,
                                  output_formatter=ScalarField.function_chart)
        #
        # Special treatment when self._dest_map != identity:
        # bases of self are created from vector frames of the ambient domain
        #
        self._induced_bases = {}
        if self._dest_map != self._domain._identity_map:
            for frame in self._ambient_domain._top_frames:
                basis = self.basis(from_frame=frame)
                self._induced_bases[frame] = basis

        # Initialization of the components of the zero element:
        for frame in self._domain._frames:
            if frame._dest_map == self._dest_map:
                self._zero_element.add_comp(frame) # since new components are
                                                   # initialized to zero

    #### Parent methods

    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct an element of the module
        """
        if comp == 0:
            return self._zero_element
        if isinstance(comp, VectorField):
            if self._domain.is_subset(comp._domain) and \
                       self._ambient_domain.is_subset(comp._ambient_domain):
                return comp.restrict(self._domain)
            else:
                raise TypeError("Cannot coerce the " + str(comp) +
                                "to a vector field in " + str(self))
        resu = self.element_class(self, name=name, latex_name=latex_name)
        if comp != []:
            resu.set_comp(basis)[:] = comp
        return resu

    # Rem: _an_element_ is declared in the superclass FiniteRankFreeModule

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent
        """
        if isinstance(other, (VectorFieldModule, VectorFieldFreeModule)):
            return self._domain.is_subset(other._domain) and \
                   self._ambient_domain.is_subset(other._ambient_domain)
        else:
            return False

    #### End of parent methods

    #### Methods to be redefined by derived classes of FiniteRankFreeModule ####

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "free module "
        if self._name is not None:
            description += self._name + " "
        description += "of vector fields "
        if self._dest_map is self._domain._identity_map:
            description += "on the " + str(self._domain)
        else:
            description += "along the " + str(self._domain) + \
                           " mapped into the " + str(self._ambient_domain)
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

        - instance of
          :class:`~sage.geometry.manifolds.tensorfield_module.TensorFieldFreeModule`
          representing the free module of type-`(k,l)` tensors on the
          free module ``self``.

        EXAMPLES:


        """
        from tensorfield_module import TensorFieldFreeModule
        if (k,l) not in self._tensor_modules:
            self._tensor_modules[(k,l)] = TensorFieldFreeModule(self, (k,l))
        return self._tensor_modules[(k,l)]

    def dual_exterior_power(self, p):
        r"""
        Return the `p`-th exterior power of the dual of ``self``.

        If ``self`` is the vector field module `\mathcal{X}(U,\Phi)`, the
        `p`-th exterior power of its dual is the set `\Lambda^p(U,\Phi)` of
        `p`-forms along `U` with values in `\Phi(U)`. It is a module over
        `C^\infty(U)`, the ring (algebra) of differentiable scalar fields on
        `U`.

        INPUT:

        - ``p`` -- non-negative integer

        OUTPUT:

        - for `p\geq 1`, instance of
          :class:`~sage.geometry.manifolds.diffform_module.DiffFormModule`
          representing the module `\Lambda^p(U,\Phi)`; for `p=0`, the
          base ring, i.e. `C^\infty(U)`, is returned instead

        EXAMPLES:

        """
        from sage.geometry.manifolds.diffform_module import DiffFormFreeModule
        if p == 0:
            return self._ring
        if p not in self._dual_exterior_powers:
            self._dual_exterior_powers[p] = DiffFormFreeModule(self, p)
        return self._dual_exterior_powers[p]

    def general_linear_group(self):
        r"""
        Return the general linear group of ``self``.

        If ``self`` is the free module `\mathcal{X}(U,\Phi)`, the *general
        linear group* is the group `\mathrm{GL}(\mathcal{X}(U,\Phi))` of
        automorphisms of `\mathcal{X}(U,\Phi)`. Note that an automorphism of
        `\mathcal{X}(U,\Phi)` can also be viewed as a *field* along `U` of
        automorphisms of the tangent spaces of `V=\Phi(U)`.

        OUTPUT:

        - instance of class
          :class:`~sage.geometry.manifolds.automorphismfield_group.AutomorphismFieldParalGroup`
          representing `\mathrm{GL}(\mathcal{X}(U,\Phi))`

        EXAMPLES:

        """
        from sage.geometry.manifolds.automorphismfield_group import \
                                                    AutomorphismFieldParalGroup
        if self._general_linear_group is None:
            self._general_linear_group = AutomorphismFieldParalGroup(self)
        return self._general_linear_group

    def basis(self, symbol=None, latex_symbol=None, from_frame=None):
        r"""
        Define a basis (vector frame) of the free module.

        If the basis specified by the given symbol already exists, it is
        simply returned.
        If no argument is provided the module's default basis is returned.

        INPUT:

        - ``symbol`` -- (string; default: ``None``) a letter (of a few letters) to
          denote a generic element of the basis; if ``None`` and ``from_frame=None``
          the module's default basis is returned.
        - ``latex_symbol`` -- (string; default: ``None``) symbol to denote a
          generic element of the basis; if ``None``, the value of ``symbol`` is
          used.
        - ``from_frame`` -- (default: ``None``) vector frame `\tilde e` on the
          codomain `V` of the destination map `\Phi` of ``self``; the returned
          basis `e` is then such that
          `\forall p \in U, e(p) = \tilde e(\Phi(p))`

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.vectorframe.VectorFrame`
          representing a basis on ``self``.

        EXAMPLES:

        """
        from vectorframe import VectorFrame
        if symbol is None:
            if from_frame is None:
                return self.default_basis()
            else:
                symbol = from_frame._symbol
                latex_symbol = from_frame._latex_symbol
        for other in self._known_bases:
            if symbol == other._symbol:
                return other
        return VectorFrame(self, symbol=symbol, latex_symbol=latex_symbol,
                           from_frame=from_frame)

    def tensor(self, tensor_type, name=None, latex_name=None, sym=None,
               antisym=None, specific_type=None):
        r"""
        Construct a tensor on the free module.

        The tensor is actually a tensor field on the domain of ``self``.

        INPUT:

        - ``tensor_type`` -- pair (k,l) with k being the contravariant rank
          and l the covariant rank
        - ``name`` -- (string; default: ``None``) name given to the tensor
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          tensor; if none is provided, the LaTeX symbol is set to ``name``
        - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
          the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention position=0 for the first argument. For instance:

          * sym=(0,1) for a symmetry between the 1st and 2nd arguments
          * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments.

        - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
          among the arguments, with the same convention as for ``sym``.
        - ``specific_type`` -- (default: ``None``) specific subclass of
          :class:`~sage.geometry.manifolds.tensorfield.TensorFieldParal` for the
          output

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.tensorfield.TensorFieldParal`
          representing the tensor defined on ``self`` with the provided
          characteristics.

        EXAMPLES:


        """
        from automorphismfield import AutomorphismField, AutomorphismFieldParal
        from metric import Metric, RiemannMetric, LorentzMetric, MetricParal, \
                           RiemannMetricParal, LorentzMetricParal
        if tensor_type==(1,0):
            return self.element_class(self, name=name, latex_name=latex_name)
        elif tensor_type==(0,1):
            return self.linear_form(name=name, latex_name=latex_name)
        elif tensor_type==(1,1) and specific_type is not None:
            if issubclass(specific_type,
                          (AutomorphismField, AutomorphismFieldParal)):
                return self.automorphism(name=name, latex_name=latex_name)
        elif tensor_type[0]==0 and tensor_type[1]>1 and antisym is not None \
                                                              and antisym !=[]:
            if isinstance(antisym, list):
                antisym0 = antisym[0]
            else:
                antisym0 = antisym
            if len(antisym0)==tensor_type[1]:
                return self.alternating_form(tensor_type[1], name=name,
                                             latex_name=latex_name)
            else:
                return self.tensor_module(*tensor_type).element_class(self,
                                 tensor_type, name=name, latex_name=latex_name,
                                 sym=sym, antisym=antisym)
        elif tensor_type==(0,2):
            if specific_type == MetricParal or specific_type == Metric:
                return MetricParal(self, name, latex_name=latex_name)
                # NB: the signature is not passed
            elif specific_type == RiemannMetricParal or \
                                                specific_type == RiemannMetric:
                return RiemannMetricParal(self, name, latex_name=latex_name)
            elif specific_type == LorentzMetricParal or \
                                                specific_type == LorentzMetric:
                return LorentzMetricParal(self, name, latex_name=latex_name)
                # NB: the signature convention is not passed
            else:
                return self.tensor_module(0,2).element_class(self, (0,2),
                                              name=name, latex_name=latex_name,
                                              sym=sym, antisym=antisym)
        # Generic case
        return self.tensor_module(*tensor_type).element_class(self,
         tensor_type, name=name, latex_name=latex_name, sym=sym,
         antisym=antisym)

    def tensor_from_comp(self, tensor_type, comp, name=None, latex_name=None):
        r"""
        Construct a tensor on the free module from a set of components.

        The tensor is actually a tensor field on the domain of ``self``.
        The tensor symmetries are deduced from those of the components.

        INPUT:

        - ``tensor_type`` -- pair (k,l) with k being the contravariant rank and l
          the covariant rank
        - ``comp`` -- instance of :class:`~sage.tensor.modules.comp.Components`
          representing the tensor components in a given basis
        - ``name`` -- (string; default: ``None``) name given to the tensor
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the tensor;
          if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.tensorfield.TensorFieldParal`
          representing the tensor defined on ``self`` with the provided
          characteristics.

        EXAMPLES:

        """
        from sage.tensor.modules.comp import CompWithSym, CompFullySym, \
                                                               CompFullyAntiSym
        #
        # 0/ Compatibility checks:
        if comp._ring is not self._ring:
             raise TypeError("the components are not defined on the same" +
                            " ring as the module")
        if comp._frame not in self._known_bases:
            raise TypeError("the components are not defined on a basis of" +
                            " the module")
        if comp._nid != tensor_type[0] + tensor_type[1]:
            raise TypeError("number of component indices not compatible with "+
                            " the tensor type")
        #
        # 1/ Construction of the tensor:
        if tensor_type == (1,0):
            resu = self.element_class(self, name=name, latex_name=latex_name)
        elif tensor_type == (0,1):
            resu = self.linear_form(name=name, latex_name=latex_name)
        elif tensor_type[0] == 0 and tensor_type[1] > 1 and \
                                        isinstance(comp, CompFullyAntiSym):
            resu = self.alternating_form(tensor_type[1], name=name,
                                         latex_name=latex_name)
        else:
            resu = self.tensor_module(*tensor_type).element_class(self,
                                 tensor_type, name=name, latex_name=latex_name)
            # Tensor symmetries deduced from those of comp:
            if isinstance(comp, CompWithSym):
                resu._sym = comp._sym
                resu._antisym = comp._antisym
        #
        # 2/ Tensor components set to comp:
        resu._components[comp._frame] = comp
        #
        return resu

    def sym_bilinear_form(self, name=None, latex_name=None):
        r"""
        Construct a symmetric bilinear form on the free module ``self``.

        A symmetric bilinear form on the vector free module ``self`` is
        actually a field of tangent-space symmetric bilinear forms along
        the open subset `U` on which ``self`` is defined.

        INPUT:

        - ``name`` -- (string; default: ``None``) name given to the automorphism
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          automorphism; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.tensorfield.TensorFieldParal` of
          tensor type (0,2) and symmetric

        """
        return self.tensor((0,2), name=name, latex_name=latex_name, sym=(0,1))

    #### End of methods to be redefined by derived classes of FiniteRankFreeModule ####

    def metric(self, name, signature=None, latex_name=None):
        r"""
        Construct a pseudo-Riemannian metric (nondegenerate symmetric bilinear
        form) on the free module ``self``.

        A metric of the vector free module ``self`` is actually a field
        of tangent-space metrics along the open subset `U` on which
        ``self`` is defined.

        INPUT:

        - ``name`` -- (string) name given to the metric
        - ``signature`` -- (integer; default: ``None``) signature `S` of the
          metric: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the number of
          positive terms (resp. number of negative terms) in any diagonal writing
          of the metric components; if ``signature`` is not provided, `S` is set to
          the manifold's dimension (Riemannian signature)
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          metric; if none, it is formed from ``name``

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.metric.MetricParal`
          representing the defined pseudo-Riemannian metric.

        See :class:`~sage.geometry.manifolds.metric.MetricParal` for further
        documentation.

        """
        from metric import MetricParal
        return MetricParal(self, name, signature=signature,
                           latex_name=latex_name)

    def riemann_metric(self, name, latex_name=None):
        r"""
        Construct a Riemannian metric (positive definite symmetric bilinear
        form) on the free module ``self``.

        A Riemannian metric of the vector free module ``self`` is actually a
        field of tangent-space Riemannian metrics along the open subset `U`
        on which ``self`` is defined.

        INPUT:

        - ``name`` -- (string) name given to the metric
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          metric; if none, it is formed from ``name``

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.metric.RiemannMetricParal`
          representing the defined Riemannian metric.

        See :class:`~sage.geometry.manifolds.metric.RiemannMetricParal` for
        further documentation.

        """
        from metric import RiemannMetricParal
        return RiemannMetricParal(self, name, latex_name=latex_name)

    def lorentz_metric(self, name, signature='positive', latex_name=None):
        r"""
        Construct a Lorentzian metric (symmetric bilinear
        form of signature (-,+,...,+) or (+,-,...,-)) on the free
        module ``self``.

        A Lorentzian metric of the vector free module ``self`` is actually a
        field of tangent-space Lorentzian metrics along the open subset `U` on
        which ``self`` is defined.

        INPUT:

        - ``name`` -- (string) name given to the metric
        - ``signature`` -- (string, default: 'positive') sign of the metric
          signature:

          * if set to 'positive', the signature is n-2, where n is the manifold's
            dimension, i.e. `(-,+,\cdots,+)`
          * if set to 'negative', the signature is -n+2, i.e. `(+,-,\cdots,-)`

        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote the
          metric; if none, it is formed from ``name``

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.metric.LorentzMetricParal`
          representing the defined Loretnzian metric.

        See :class:`~sage.geometry.manifolds.metric.LorentzMetricParal` for
        further documentation.

        """
        from metric import LorentzMetricParal
        return LorentzMetricParal(self, name, signature=signature,
                                  latex_name=latex_name)
