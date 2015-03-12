r"""
Differentiable mappings between manifolds

The class :class:`DiffMapping` implements differentiable mappings from an open
subset `U` of a differentiable manifold `M` to some differentiable
manifold `N`:

.. MATH::

    \Phi: U\subset M \longrightarrow N


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version

REFERENCES:

- Chap. 1 of S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*,
  vol. 1, Interscience Publishers (New York) (1963)
- Chaps. 2 and 3 of J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed.,
  Springer (New York) (2013)


"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.morphism import Morphism
from sage.geometry.manifolds.domain import ManifoldSubset
from sage.geometry.manifolds.chart import FunctionChart, MultiFunctionChart

class DiffMapping(Morphism):
    r"""
    Differentiable mapping between two manifolds.

    This class implements differentiable mappings of the type

    .. MATH::

        \Phi: U\subset M \longrightarrow V\subset N

    where `M` and `N` are differentiable manifolds, `U` is an open subset
    of `M` and `V` is an open subset of `N`.

    In what follows, `M` is called the *start manifold* and
    `N` the *arrival manifold*.

    Differentiable mappings are the *morphisms* of the *category* of
    differentiable manifolds. The set of all differentiable mappings from
    `U` to `V` is therefore the homset between `U` and `V` and is denoted
    by `\mathrm{Hom}(U,V)`.

    The class :class:`DiffMapping` is a Sage *element* class, whose *parent*
    class is :class:`~sage.geometry.manifolds.manifold_homset.ManifoldHomset`.

    INPUT:

    - ``parent`` -- homset `\mathrm{Hom}(U,V)` to which the differentiable
      mapping belongs
    - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
      a dictionary of the coordinate expressions (as lists (or tuples) of the
      coordinates of the image expressed in terms of the coordinates of
      the considered point) with the pairs of charts (chart1, chart2)
      as keys (chart1 being a chart on `U` and chart2 a chart on `V`).
      If the dimension of the arrival manifold is 1, a single coordinate
      expression can be passed instead of a tuple with a single element
    - ``name`` -- (default: ``None``) name given to the differentiable mapping
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      differentiable mapping; if none is provided, the LaTeX symbol is set to
      ``name``
    - ``is_diffeomorphism`` -- (default: ``False``) determines whether the
      constructed object is a diffeomorphism; if set to ``True``,
      then the manifolds `M` and `N` must have the same dimension.
    - ``is_identity`` -- (default: ``False``) determines whether the
      constructed object is the identity map; if set to ``True``,
      then `V` must be `U` and the entry ``coord_functions`` is not used.

    .. NOTE::

        If the information passed by means of the argument ``coord_functions``
        is not sufficient to fully specify the differentiable mapping,
        further coordinate expressions, in other charts, can be subsequently
        added by means of the method :meth:`add_expr`

    EXAMPLES:

    The standard embedding of the sphere `S^2` into `\RR^3`::

        sage: M = Manifold(2, 'S^2') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                             intersection_name='W', restrictions1= x^2+y^2!=0, \
                                             restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: N = Manifold(3, 'R^3', r'\RR^3')  # R^3
        sage: c_cart.<X,Y,Z> = N.chart()  # Cartesian coordinates on R^3
        sage: Phi = M.diff_mapping(N, \
        ....: {(c_xy, c_cart): [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)],  \
        ....:  (c_uv, c_cart): [2*u/(1+u^2+v^2), 2*v/(1+u^2+v^2), (1-u^2-v^2)/(1+u^2+v^2)]}, \
        ....: name='Phi', latex_name=r'\Phi')
        sage: Phi
        differentiable mapping 'Phi' from the 2-dimensional manifold 'S^2' to
         the 3-dimensional manifold 'R^3'
        sage: Phi.parent()
        Set of Morphisms from 2-dimensional manifold 'S^2' to 3-dimensional
         manifold 'R^3' in Category of sets
        sage: Phi.parent() is Hom(M, N)
        True
        sage: type(Phi)
        <class 'sage.geometry.manifolds.diffmapping.ManifoldHomset_with_category.element_class'>
        sage: Phi.display()
        Phi: S^2 --> R^3
        on U: (x, y) |--> (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
        on V: (u, v) |--> (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

    The mapping can be initialized by the method
    :meth:`~sage.geometry.manifolds.domain.ManifoldOpenSubset.diff_mapping`
    only in a single pair of charts: the argument ``coord_functions`` is then
    a mere list of coordinate expressions (and not a dictionary) and the
    arguments ``chart1`` and ``chart2`` have to be provided if the charts
    differ from the default ones on the domain and/or the codomain::

        sage: Phi1 = M.diff_mapping(N, [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)], \
        ....: chart1=c_xy, chart2=c_cart, name='Phi', latex_name=r'\Phi')

    Since c_xy and c_cart are the default charts on respectively M and N, they
    can be omitted, so that the above declaration is equivalent to::

        sage: Phi1 = M.diff_mapping(N, [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)], \
        ....: name='Phi', latex_name=r'\Phi')

    With such a declaration, the differentiable mapping is only partially defined
    on the manifold `S^2`, being known in only one chart::

        sage: Phi1.display()
        Phi: S^2 --> R^3
        on U: (x, y) |--> (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))

    The definition can be completed by means of the method :meth:`add_expr`::

        sage: Phi1.add_expr(c_uv, c_cart, [2*u/(1+u^2+v^2), 2*v/(1+u^2+v^2), (1-u^2-v^2)/(1+u^2+v^2)])
        sage: Phi1.display()
        Phi: S^2 --> R^3
        on U: (x, y) |--> (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
        on V: (u, v) |--> (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

    At this stage, Phi1 and Phi are fully equivalent::

        sage: Phi1 == Phi
        True

    The test suite is passed::

        sage: TestSuite(Phi).run()
        sage: TestSuite(Phi1).run()

    The mapping acts on points::

        sage: np = M.point((0,0), chart=c_uv)  # the North pole
        sage: Phi(np)
        point on 3-dimensional manifold 'R^3'
        sage: Phi(np).coord() # Cartesian coordinates
        (0, 0, 1)
        sage: sp = M.point((0,0), chart=c_xy)  # the South pole
        sage: Phi(sp).coord() # Cartesian coordinates
        (0, 0, -1)

    Differential mappings can be composed by means of the operator ``*``: let
    us introduce the mapping `\RR^3\rightarrow \RR^2` corresponding to
    the projection from the point `(X,Y,Z)=(0,0,1)` onto the equatorial plane
    `Z=0`::

        sage: P = Manifold(2, 'R^2', r'\RR^2') # R^2 (equatorial plane)
        sage: cP.<xP, yP> = P.chart()
        sage: Psi = N.diff_mapping(P, (X/(1-Z), Y/(1-Z)), name='Psi',
        ....:                      latex_name=r'\Psi')
        sage: Psi
        differentiable mapping 'Psi' from the 3-dimensional manifold 'R^3' to
         the 2-dimensional manifold 'R^2'
        sage: Psi.display()
        Psi: R^3 --> R^2
           (X, Y, Z) |--> (xP, yP) = (-X/(Z - 1), -Y/(Z - 1))

    Then we compose ``Psi`` with ``Phi``, thereby getting a mapping
    `S^2\rightarrow \RR^2`::

        sage: ster = Psi*Phi ; ster
        differentiable mapping from the 2-dimensional manifold 'S^2' to the
         2-dimensional manifold 'R^2'

    Let us test on the South pole (``sp``) that ``ster`` is indeed the
    composite of ``Psi`` and ``Phi``::

        sage: ster(sp) == Psi(Phi(sp))
        True

    Actually ``ster`` is the stereographic projection from the North pole, as
    its coordinate expression reveals::

        sage: ster.display()
        S^2 --> R^2
        on U: (x, y) |--> (xP, yP) = (x, y)
        on V: (u, v) |--> (xP, yP) = (u/(u^2 + v^2), v/(u^2 + v^2))

    If the arrival manifold is 1-dimensional, a differentiable mapping must be
    defined by a single symbolic expression for each pair of charts, and not
    by a list/tuple with a single element::

        sage: N = Manifold(1, 'N')
        sage: c_N = N.chart('X')
        sage: Phi = M.diff_mapping(N, {(c_xy, c_N): x^2+y^2, \
        ....: (c_uv, c_N): 1/(u^2+v^2)})  # not ...[1/(u^2+v^2)] or (1/(u^2+v^2),)

    If the arrival manifold is the field of real numbers `\RR` (represented
    by :class:`~sage.geometry.manifolds.manifold.RealLine`), the action on a
    point returns a real number, i.e. the canonical coordinate of the image
    point, and not the image point itself::

        sage: Phi = M.diff_mapping(RealLine(), x^2+y^2)
        sage: Phi(M.point((1,2)))
        5

    An example of differentiable mapping `\RR \rightarrow \RR^2`::

        sage: R.<t> = RealLine()   # field R with canonical coordinate t
        sage: R2 = Manifold(2, 'R^2') # R^2
        sage: c_xy.<x,y> = R2.chart() # Cartesian coordinates on R^2
        sage: Phi = R.diff_mapping(R2, [cos(t), sin(t)], name='Phi') ; Phi
        Curve 'Phi' in the 2-dimensional manifold 'R^2'
        sage: Phi.parent()
        Set of Morphisms from field R of real numbers to 2-dimensional manifold
         'R^2' in Category of sets
        sage: Phi.parent() is Hom(R, R2)
        True
        sage: Phi.display()
        Phi: R --> R^2
           t |--> (x, y) = (cos(t), sin(t))

    An example of diffeomorphism between the unit open disk and the Euclidean
    plane `\RR^2`::

        sage: D = R2.open_subset('D', coord_def={c_xy: x^2+y^2<1}) # the open unit disk
        sage: Phi = D.diffeomorphism(R2, [x/sqrt(1-x^2-y^2), y/sqrt(1-x^2-y^2)],
        ....:                        name='Phi', latex_name=r'\Phi')
        sage: Phi
        diffeomorphism 'Phi' from the open subset 'D' of the 2-dimensional
         manifold 'R^2' to the 2-dimensional manifold 'R^2'
        sage: Phi.parent()
        Set of Morphisms from open subset 'D' of the 2-dimensional manifold
         'R^2' to 2-dimensional manifold 'R^2' in Category of facade sets
        sage: Phi.parent() is Hom(D, R2)
        True
        sage: Phi.display()
        Phi: D --> R^2
           (x, y) |--> (x, y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))

    The image of a point::

        sage: p = D.point((1/2,0))
        sage: q = Phi(p) ; q
        point on 2-dimensional manifold 'R^2'
        sage: q.coord()
        (1/3*sqrt(3), 0)

    The inverse diffeomorphism is computed by means of the method :meth:`inverse`::

        sage: Phi.inverse()
        diffeomorphism 'Phi^(-1)' from the 2-dimensional manifold 'R^2' to the
         open subset 'D' of the 2-dimensional manifold 'R^2'

    Equivalently, one may use the notations ``^(-1)`` or ``~`` to get the
    inverse::

        sage: Phi^(-1) is Phi.inverse()
        True
        sage: ~Phi is Phi.inverse()
        True

    Check that ``~Phi`` is indeed the inverse of ``Phi``::

        sage: (~Phi)(q) == p
        True
        sage: Phi * ~Phi == R2.identity_map()
        True
        sage: ~Phi * Phi == D.identity_map()
        True

    The coordinate expression of the inverse diffeomorphism::

        sage: (~Phi).display()
        Phi^(-1): R^2 --> D
           (x, y) |--> (x, y) = (x/sqrt(x^2 + y^2 + 1), y/sqrt(x^2 + y^2 + 1))

    A special case of diffeomorphism: the identity map of the open unit disk::

        sage: id = D.identity_map() ; id
        identity map 'Id_D' of the open subset 'D' of the 2-dimensional manifold 'R^2'
        sage: latex(id)
        \mathrm{Id}_{D}
        sage: id.parent()
        Set of Morphisms from open subset 'D' of the 2-dimensional manifold 'R^2' to open subset 'D' of the 2-dimensional manifold 'R^2' in Category of facade sets
        sage: id.parent() is Hom(D, D)
        True

    The identity map acting on a point::

        sage: id(p)
        point on 2-dimensional manifold 'R^2'
        sage: id(p) == p
        True
        sage: id(p) is p
        True

    The coordinate expression of the identity map::

        sage: id.display()
        Id_D: D --> D
           (x, y) |--> (x, y)

    The identity map is its own inverse::

        sage: id^(-1) is id
        True
        sage: ~id is id
        True

    """
    def __init__(self, parent, coord_functions=None, chart1=None, chart2=None,
                 name=None, latex_name=None, is_diffeomorphism=False,
                 is_identity=False):
        r"""
        Construct a differentiable mapping.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: f = Hom(M,N)({(X,Y): (x+y, x*y, x-y)}, name='f') ; f
            differentiable mapping 'f' from the 2-dimensional manifold 'M' to the 3-dimensional manifold 'N'
            sage: f.display()
            f: M --> N
               (x, y) |--> (u, v, w) = (x + y, x*y, x - y)
            sage: TestSuite(f).run()

        The identity map::

            sage: f = Hom(M,M)({}, is_identity=True) ; f
            identity map 'Id_M' of the 2-dimensional manifold 'M'
            sage: f.display()
            Id_M: M --> M
               (x, y) |--> (x, y)
            sage: TestSuite(f).run()

        """
        Morphism.__init__(self, parent)
        domain = parent.domain()
        codomain = parent.codomain()
        self._domain = domain
        self._codomain = codomain
        self._coord_expression = {}
        self._is_diffeo = False
        self._is_identity = False
        if is_identity:
            # Construction of the identity map
            self._is_identity = True
            self._is_diffeo = True
            if domain != codomain:
                raise ValueError("the domain and codomain must coincide " + \
                                 "for the identity map")
            if name is None:
                name = 'Id_' + domain._name
            if latex_name is None:
                latex_name = r'\mathrm{Id}_{' + domain._latex_name + r'}'
            self._name = name
            self._latex_name = latex_name
            for chart in domain.atlas():
                coord_funct = chart[:]
                self._coord_expression[(chart, chart)] = \
                                        MultiFunctionChart(chart, *coord_funct)
        else:
            # Construction of a generic differentiable mapping
            if is_diffeomorphism:
                self._is_diffeo = True
                if domain._manifold.dim() != codomain._manifold.dim():
                    raise ValueError("for a diffeomorphism, the source " +
                                     "manifold and target manifold must " +
                                     "have the same dimension")
            if coord_functions is not None:
                n2 = self._codomain._manifold._dim
                for chart_pair, expression in coord_functions.iteritems():
                    if chart_pair[0] not in self._domain._atlas:
                        raise ValueError("{} is not a chart ".format(
                                                              chart_pair[0]) +
                                     "defined on the {}".format(self._domain))
                    if chart_pair[1] not in self._codomain._atlas:
                        raise ValueError("{} is not a chart ".format(
                                                              chart_pair[1]) +
                                   " defined on the {}".format(self._codomain))
                    if n2 == 1:
                        # a single expression entry is allowed (instead of a
                        # tuple)
                        if not isinstance(expression, (tuple, list)):
                            expression = (expression,)
                    if len(expression) != n2:
                        raise ValueError("{} coordinate ".format(n2) +
                                         "functions must be provided")
                    self._coord_expression[chart_pair] = \
                             MultiFunctionChart(chart_pair[0], *expression)
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
            else:
                self._latex_name = latex_name
        # Initialization of derived quantities:
        DiffMapping._init_derived(self)

    #
    # SageObject methods
    #

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(2, 'N')
            sage: Y.<u,v> = N.chart()
            sage: f = Hom(M,N)({(X,Y): (x+y,x*y)})
            sage: f._repr_()
            "differentiable mapping from the 2-dimensional manifold 'M' to the
             2-dimensional manifold 'N'"
            sage: f = Hom(M,N)({(X,Y): (x+y,x*y)}, name='f')
            sage: f._repr_()
            "differentiable mapping 'f' from the 2-dimensional manifold 'M' to
             the 2-dimensional manifold 'N'"
            sage: f = Hom(M,N)({(X,Y): (x+y,x-y)}, name='f', is_diffeomorphism=True)
            sage: f._repr_()
            "diffeomorphism 'f' from the 2-dimensional manifold 'M' to the
             2-dimensional manifold 'N'"
            sage: f = Hom(M,M)({(X,X): (x+y,x-y)}, name='f', is_diffeomorphism=True)
            sage: f._repr_()
            "diffeomorphism 'f' of the 2-dimensional manifold 'M'"
            sage: f = Hom(M,M)({}, name='f', is_identity=True)
            sage: f._repr_()
            "identity map 'f' of the 2-dimensional manifold 'M'"

        """
        if self._is_identity:
            return "identity map '" + self._name + \
                   "' of the {}".format(self._domain)
        if self._is_diffeo:
            description = "diffeomorphism"
        else:
            description = "differentiable mapping"
        if self._name is not None:
            description += " '%s'" % self._name
        if self._domain == self._codomain:
            if self._is_diffeo:
                description += " of the {}".format(self._domain)
            else:
                description += " from the {} to itself".format(self._domain)
        else:
            description += " from the {} to the {}".format(self._domain,
                                                           self._codomain)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = Hom(M,M)({(X,X): (x+y,x*y)}, name='f')
            sage: f._latex_()
            'f'
            sage: f = Hom(M,M)({(X,X): (x+y,x*y)}, name='f', latex_name=r'\Phi')
            sage: f._latex_()
            '\\Phi'

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- another instance of :class:`DiffMapping` to compare with

        OUTPUT:

        - True if ``self`` is equal to ``other``,  or False otherwise

        """
        if not isinstance(other, DiffMapping):
            return False
        if self.parent() != other.parent():
            return False
        if self._is_identity:
            return other.is_identity()
        if other._is_identity:
            return self.is_identity()
        for charts, coord_functions in self._coord_expression.iteritems():
            try:
                if coord_functions.expr() != other.expr(*charts):
                    return False
            except ValueError:
                return False
        return True

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- another instance of :class:`DiffMapping` to compare with

        OUTPUT:

        - True if ``self`` is different from ``other``,  or False otherwise

        """
        return not self.__eq__(other)

    def __cmp__(self, other):
        r"""
        Old-style (Python 2) comparison operator.

        This is provisory, until migration to Python 3 is achieved.

        """
        if self.__eq__(other):
            return 0
        else:
            return -1

    #
    # Map methods
    #

    def _call_(self, point):
        r"""
        Compute the image of a point by ``self``.

        INPUT:

        - ``point`` -- point in the domain of ``self``, as an instance of
          :class:`~sage.geometry.manifolds.point.ManifoldPoint`

        OUTPUT:

        - image of the point by ``self`` (instance of
          :class:`~sage.geometry.manifolds.point.ManifoldPoint`)

        EXAMPLES:

        Planar rotation acting on a point::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'R^2', r'\RR^2') # Euclidean plane
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates
            sage: # A pi/3 rotation around the origin defined in Cartesian coordinates:
            sage: rot = M.diff_mapping(M, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
            sage: p = M.point((1,2), name='p')
            sage: q = rot(p) ; q
            point 'R(p)' on 2-dimensional manifold 'R^2'
            sage: q.coord()
            (-sqrt(3) + 1/2, 1/2*sqrt(3) + 1)

        Image computed after some change of coordinates::

            sage: c_spher.<r,ph> = M.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coord. on the plane
            sage: ch = c_spher.coord_change(c_cart, r*cos(ph), r*sin(ph))
            sage: p1 = M.point((sqrt(5), arctan(2)), chart=c_spher) # p1 is defined only in terms of c_spher
            sage: q1 = rot(p1) # but the computation of the action of rot is still possible
            sage: q1 == q
            True

        Image computed by means of spherical coordinates::

            sage: rot.add_expr(c_spher, c_spher, (r, ph+pi/3)) # now rot is known in terms of c_spher
            sage: p2 = M.point((sqrt(5), arctan(2)), chart=c_spher)
            sage: q2 = rot(p2) # computation on c_spher
            sage: q2 == q
            True

        """
        # NB: checking that ``point`` belongs to the mapping's domain has been
        # already performed by Map.__call__(); this check is therefore not
        # repeated here.
        from manifold import RealLine
        if self._is_identity:
            return point
        dom = self._domain
        chart1, chart2 = None, None
        for chart in point._coordinates:
            for chart_pair in self._coord_expression:
                if chart_pair[0] is chart:
                    chart1 = chart
                    chart2 = chart_pair[1]
                    break
            if chart1 is not None:
                break
        else:
            # attempt to perform a change of coordinate on the point
            for chart_pair in self._coord_expression:
                try:
                    point.coord(chart_pair[0])
                    chart1, chart2 = chart_pair
                except ValueError:
                    pass
                if chart1 is not None:
                    break
            else:
                raise ValueError("no pair of charts has been found to " +
                  "compute the action of the {} on the {}".format(self, point))
        coord_map = self._coord_expression[(chart1, chart2)]
        y = coord_map(*(point._coordinates[chart1]))
        if isinstance(self._codomain._manifold, RealLine):
            # special case of a mapping to R
            return y[0]
        else:
            if point._name is None or self._name is None:
                res_name = None
            else:
                res_name = self._name + '(' + point._name + ')'
            if point._latex_name is None or self._latex_name is None:
                res_latex_name = None
            else:
                res_latex_name = self._latex_name + r'\left(' + \
                                 point._latex_name + r'\right)'
            # The image point is created as an element of the domain of chart2:
            dom2 = chart2.domain()
            return dom2.element_class(dom2, coords=y, chart=chart2,
                                      name=res_name, latex_name=res_latex_name,
                                      check_coords=False)
    #
    # Morphism methods
    #

    def is_identity(self):
        r"""
        Check whether ``self`` is an identity map.

        EXAMPLES:

        Tests on differentiable mappings of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: M.identity_map().is_identity()  # obviously...
            True
            sage: Hom(M, M).one().is_identity()  # a variant of the obvious
            True
            sage: a = M.diff_mapping(M, coord_functions={(X,X): (x, y)})
            sage: a.is_identity()
            True
            sage: a = M.diff_mapping(M, coord_functions={(X,X): (x, y+1)})
            sage: a.is_identity()
            False

        Of course, if the codomain of ``self`` does not coincide with its
        domain, the outcome is ``False``::

            sage: N = Manifold(2, 'N')
            sage: Y.<u,v> = N.chart()
            sage: a = M.diff_mapping(N, {(X,Y): (x, y)})
            sage: a.display()
            M --> N
               (x, y) |--> (u, v) = (x, y)
            sage: a.is_identity()
            False

        """
        if self._is_identity:
            return True
        if self._codomain != self._domain:
            return False
        for chart in self._domain._top_charts:
            try:
                if chart[:] != self.expr(chart, chart):
                    return False
            except ValueError:
                return False
        # If this point is reached, ``self`` must be the identity:
        self._is_identity = True
        return True

    def _composition_(self, other, homset):
        r"""
        Composition of ``self`` with another morphism.

        The composition is performed on the right, i.e. the returned
        morphism is ``self*other``.

        INPUT:

        - ``other`` -- a differentiable mapping, whose codomain is the domain
          of ``self``
        - ``homset`` -- the homset of the differentiable mapping ``self*other``;
          this argument is required to follow the prototype of
          :meth:`~sage.categories.map.Map._composition_` and is determined by
          :meth:`~sage.categories.map.Map._composition` (single underscore),
          that is supposed to call the current method

        OUTPUT:

        - the composite mapping ``self*other``, as an instance of
          :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`

        """
        # This method is invoked by Map._composition (single underscore),
        # which is itself invoked by Map.__mul__ . The latter performs the
        # check other._codomain == self._domain. There is therefore no need
        # to perform it here.
        if self._is_identity:
            return other
        if other._is_identity:
            return self
        resu_funct = {}
        for chart1 in other._domain._top_charts:
            for chart2 in self._domain._top_charts:
                for chart3 in self._codomain._top_charts:
                    try:
                        self23 = self.multi_function_chart(chart2, chart3)
                        resu_funct[(chart1, chart3)] = \
                                self23(*(other.expr(chart1, chart2)),
                                       simplify=True)
                    except ValueError:
                        pass
        return homset(resu_funct)

    #
    # Monoid methods
    #

    def _mul_(self, other):
        r"""
        Composition of ``self`` with another morphism (endomorphism case).

        This applies only when the parent of ``self`` is a monoid, i.e. when
        ``self`` is an endomorphism of the category of differentiable manifolds,
        i.e. a differentiable mapping U --> U, where U is some open subset of
        a differentiable manifold.

        INPUT:

        - ``other`` -- a differentiable mapping, whose codomain is the domain
          of ``self``

        OUTPUT:

        - the composite mapping ``self*other``, as an instance of
          :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`

        """
        from sage.categories.homset import Hom
        dom = self._domain
        return self._composition_(other, Hom(dom, dom))

    #
    # Other methods
    #

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        self._restrictions = {} # dict. of restrictions to subdomains of
                                # self._domain
        if self._is_identity:
            self._inverse = self
        else:
            self._inverse = None
        self._diff = {} # dict. of the coord. expressions of the differential
                        # keys: pair of charts

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        self._restrictions.clear()
        if not self._is_identity:
            self._inverse = None
        self._diff.clear()

    def _display_expression(self, chart1, chart2, result):
        r"""
        Helper function for :meth:`display`.
        """
        from sage.misc.latex import latex
        try:
            expression = self.expr(chart1, chart2)
            coords1 = chart1[:]
            if len(coords1) == 1:
                coords1 = coords1[0]
            coords2 = chart2[:]
            if len(coords2) == 1:
                coords2 = coords2[0]
            if chart1._domain == self._domain:
                result._txt += "   "
                result._latex += " & "
            else:
                result._txt += "on " + chart1._domain._name + ": "
                result._latex += r"\mbox{on}\ " + latex(chart1._domain) + \
                                r": & "
            result._txt += repr(coords1) + " |--> "
            result._latex += latex(coords1) + r"& \longmapsto & "
            if chart2 == chart1:
                if len(expression) == 1:
                    result._txt += repr(expression[0]) + "\n"
                    result._latex += latex(expression[0]) + r"\\"
                else:
                    result._txt += repr(expression) + "\n"
                    result._latex += latex(expression) + r"\\"
            else:
                if len(expression) == 1:
                    result._txt += repr(coords2[0]) + " = " + \
                                  repr(expression[0]) + "\n"
                    result._latex += latex(coords2[0]) + " = " + \
                                    latex(expression[0]) + r"\\"
                else:
                    result._txt += repr(coords2) + " = " + \
                                  repr(expression) + "\n"
                    result._latex += latex(coords2) + " = " + \
                                    latex(expression) + r"\\"
        except (TypeError, ValueError):
            pass

    def display(self, chart1=None, chart2=None):
        r"""
        Display the expression of the differentiable mapping in one or more
        pair of charts.

        If the expression is not known already, it is computed from some
        expression in other charts by means of change-of-coordinate formulas.

        INPUT:

        - ``chart1`` -- (default: None) chart on the mapping's domain; if None,
          the display is performed on all the charts on the start manifold
          in which the mapping is known or computable via some change of
          coordinates
        - ``chart2`` -- (default: None) chart on the mapping's codomain; if
          None, the display is performed on all the charts on the codomain
          in which the mapping is known or computable via some change of
          coordinates

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        EXAMPLES:

        Standard embedding of the sphere `S^2` in `\RR^3`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'S^2') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: N = Manifold(3, 'R^3', r'\RR^3')  # R^3
            sage: c_cart.<X,Y,Z> = N.chart()  # Cartesian coordinates on R^3
            sage: Phi = M.diff_mapping(N, \
            ....: {(c_xy, c_cart): [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)],  \
            ....:  (c_uv, c_cart): [2*u/(1+u^2+v^2), 2*v/(1+u^2+v^2), (1-u^2-v^2)/(1+u^2+v^2)]}, \
            ....: name='Phi', latex_name=r'\Phi')
            sage: Phi.display(c_xy, c_cart)
            Phi: S^2 --> R^3
            on U: (x, y) |--> (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
            sage: Phi.display(c_uv, c_cart)
            Phi: S^2 --> R^3
            on V: (u, v) |--> (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

        The LaTeX output::

            sage: latex(Phi.display(c_xy, c_cart))
            \begin{array}{llcl} \Phi:& S^2 & \longrightarrow & \RR^3 \\ \mbox{on}\ U : & \left(x, y\right) & \longmapsto & \left(X, Y, Z\right) = \left(\frac{2 \, x}{x^{2} + y^{2} + 1}, \frac{2 \, y}{x^{2} + y^{2} + 1}, \frac{x^{2} + y^{2} - 1}{x^{2} + y^{2} + 1}\right) \end{array}

        If the argument ``chart2`` is not specified, the display is performed
        on all the charts on the arrival manifold in which the mapping is known
        or computable via some change of coordinates (here only one chart:
        c_cart)::

            sage: Phi.display(c_xy)
            Phi: S^2 --> R^3
            on U: (x, y) |--> (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))

        Similarly, if the argument ``chart1`` is omitted, the display is
        performed on all the charts on the start manifold in which the
        mapping is known or computable via some change of coordinates::

            sage: Phi.display(chart2=c_cart)
            Phi: S^2 --> R^3
            on U: (x, y) |--> (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
            on V: (u, v) |--> (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

        If neither ``chart1`` nor ``chart2`` is specified, the display is
        performed on all the pair of charts in which the mapping is known or
        computable via some change of coordinates::

            sage: Phi.display()
            Phi: S^2 --> R^3
            on U: (x, y) |--> (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
            on V: (u, v) |--> (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

        If a chart covers entirely the mapping's domain, the mention "on ..."
        is omitted::

            sage: Phi.restrict(U).display()
            Phi: U --> R^3
               (x, y) |--> (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))

        A shortcut of ``display()`` is ``disp()``::

            sage: Phi.disp()
            Phi: S^2 --> R^3
            on U: (x, y) |--> (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
            on V: (u, v) |--> (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

        """
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion
        result = FormattedExpansion()
        if self._name is None:
            symbol = ""
        else:
            symbol = self._name + ": "
        result._txt = symbol + self._domain._name + " --> " + \
                     self._codomain._name + "\n"
        if self._latex_name is None:
            symbol = ""
        else:
            symbol = self._latex_name + ":"
        result._latex = r"\begin{array}{llcl} " + symbol + r"&" + \
                       latex(self._domain) + r"& \longrightarrow & " + \
                       latex(self._codomain) + r"\\"
        if chart1 is None:
            if chart2 is None:
                for ch1 in self._domain._atlas:
                    for ch2 in self._codomain._atlas:
                        self._display_expression(ch1, ch2, result)
            else:
                for ch1 in self._domain._atlas:
                    self._display_expression(ch1, chart2, result)
        else:
            if chart2 is None:
                for ch2 in self._codomain._atlas:
                    self._display_expression(chart1, ch2, result)
            else:
                self._display_expression(chart1, chart2, result)
        result._txt = result._txt[:-1]
        result._latex = result._latex[:-2] + r"\end{array}"
        return result

    disp = display

    def view(self, chart1=None, chart2=None):
        r"""
        Deprecated method.

        Use method :meth:`display` instead.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: Phi = M.diff_mapping(M, (x+y, x-y), name='Phi')
            sage: Phi.view()
            doctest:...: DeprecationWarning: Use function display() instead.
            See http://trac.sagemath.org/15916 for details.
            Phi: M --> M
               (x, y) |--> (x + y, x - y)
            sage: Phi.display()
            Phi: M --> M
               (x, y) |--> (x + y, x - y)

        """
        from sage.misc.superseded import deprecation
        deprecation(15916, 'Use function display() instead.')
        return self.display(chart1=chart1, chart2=chart2)

    def multi_function_chart(self, chart1=None, chart2=None):
        r"""
        Return the functions of the coordinates representing the differentiable
        mapping in a given pair of charts.

        If these functions are not already known, they are computed from known
        ones by means of change-of-chart formulas.

        INPUT:

        - ``chart1`` -- (default: None) chart on the mapping's domain; if None,
          the domain's default chart is assumed
        - ``chart2`` -- (default: None) chart on the mapping's codomain; if
          None,  the codomain's default chart is assumed

        OUTPUT:

        - instance of class
          :class:`~sage.geometry.manifolds.chart.MultiFunctionChart`
          representing the differentiable mapping in the above two charts

        EXAMPLES:

        Differential mapping from a 2-dimensional manifold to a 3-dimensional
        one::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: N = Manifold(3, 'N')
            sage: c_uv.<u,v> = M.chart()
            sage: c_xyz.<x,y,z> = N.chart()
            sage: Phi = M.diff_mapping(N, (u*v, u/v, u+v), name='Phi', latex_name=r'\Phi')
            sage: Phi.display()
            Phi: M --> N
               (u, v) |--> (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.multi_function_chart(c_uv, c_xyz)
            functions (u*v, u/v, u + v) on the chart (M, (u, v))
            sage: Phi.multi_function_chart() # equivalent to above since 'uv' and 'xyz' are default charts
            functions (u*v, u/v, u + v) on the chart (M, (u, v))
            sage: type(Phi.multi_function_chart())
            <class 'sage.geometry.manifolds.chart.MultiFunctionChart'>

        Representation in other charts::

            sage: c_UV.<U,V> = M.chart()  # new chart on M
            sage: ch_uv_UV = c_uv.coord_change(c_UV, u-v, u+v)
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ.<X,Y,Z> = N.chart() # new chart on N
            sage: ch_xyz_XYZ = c_xyz.coord_change(c_XYZ, 2*x-3*y+z, y+z-x, -x+2*y-z)
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.multi_function_chart(c_UV, c_xyz)
            functions (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V) on the chart (M, (U, V))
            sage: Phi.multi_function_chart(c_uv, c_XYZ)
            functions (((2*u + 1)*v^2 + u*v - 3*u)/v, -((u - 1)*v^2 - u*v - u)/v, -((u + 1)*v^2 + u*v - 2*u)/v) on the chart (M, (u, v))
            sage: Phi.multi_function_chart(c_UV, c_XYZ)
            functions (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V), 1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V), 1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V)) on the chart (M, (U, V))

        """
        dom1 = self._domain; dom2 = self._codomain
        def_chart1 = dom1._def_chart; def_chart2 = dom2._def_chart
        if chart1 is None:
            chart1 = def_chart1
        if chart2 is None:
            chart2 = def_chart2
        if (chart1, chart2) not in self._coord_expression:
            # Some computation must be performed
            if self._is_identity and chart1 == chart2:
                # special case of the identity in a single chart:
                coord_functions = chart1[:]
                self._coord_expression[(chart1, chart1)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
                return self._coord_expression[(chart1, chart2)]
            # Some change of coordinates must be performed
            change_start = [] ; change_arrival = []
            for (ochart1, ochart2) in self._coord_expression:
                if chart1 == ochart1:
                    change_arrival.append(ochart2)
                if chart2 == ochart2:
                    change_start.append(ochart1)
            # 1/ Trying to make a change of chart only on the arrival domain:
            # the arrival default chart is privileged:
            sel_chart2 = None # selected chart2
            if def_chart2 in change_arrival \
                    and (def_chart2, chart2) in dom2._coord_changes:
                sel_chart2 = def_chart2
            else:
                for ochart2 in change_arrival:
                    if (ochart2, chart2) in dom2._coord_changes:
                        sel_chart2 = ochart2
                        break
            if sel_chart2 is not None:
                oexpr = self._coord_expression[(chart1, sel_chart2)]
                chg2 = dom2._coord_changes[(sel_chart2, chart2)]
                self._coord_expression[(chart1, chart2)] = \
                    MultiFunctionChart(chart1, *(chg2(*(oexpr.expr()))) )
                return self._coord_expression[(chart1, chart2)]

            # 2/ Trying to make a change of chart only on the start domain:
            # the start default chart is privileged:
            sel_chart1 = None # selected chart1
            if def_chart1 in change_start \
                    and (chart1, def_chart1) in dom1._coord_changes:
                sel_chart1 = def_chart1
            else:
                for ochart1 in change_start:
                    if (chart1, ochart1) in dom1._coord_changes:
                        sel_chart1 = ochart1
                        break
            if sel_chart1 is not None:
                oexpr = self._coord_expression[(sel_chart1, chart2)]
                chg1 = dom1._coord_changes[(chart1, sel_chart1)]
                self._coord_expression[(chart1, chart2)] = \
                    MultiFunctionChart(chart1,
                                       *(oexpr( *(chg1._transf.expr()) )) )
                return self._coord_expression[(chart1, chart2)]

            # 3/ If this point is reached, it is necessary to perform some
            # coordinate change both on the start domain and the arrival one
            # the default charts are privileged:
            if (def_chart1, def_chart2) in self._coord_expression \
                    and (chart1, def_chart1) in dom1._coord_changes \
                    and (def_chart2, chart2) in dom2._coord_changes:
                sel_chart1 = def_chart1
                sel_chart2 = def_chart2
            else:
                for (ochart1, ochart2) in self._coord_expression:
                    if (chart1, ochart1) in dom1._coord_changes \
                        and (ochart2, chart2) in dom2._coord_changes:
                        sel_chart1 = ochart1
                        sel_chart2 = ochart2
                        break
            if (sel_chart1 is not None) and (sel_chart2 is not None):
                oexpr = self._coord_expression[(sel_chart1, sel_chart2)]
                chg1 = dom1._coord_changes[(chart1, sel_chart1)]
                chg2 = dom2._coord_changes[(sel_chart2, chart2)]
                self._coord_expression[(chart1, chart2)] = \
                     MultiFunctionChart(chart1,
                                *(chg2( *(oexpr(*(chg1._transf.expr()))) )) )
                return self._coord_expression[(chart1, chart2)]

            # 4/ If this point is reached, the demanded value cannot be
            # computed
            raise ValueError("The expression of the mapping in the pair of " +
                "charts (" + str(chart1) + ", " + str(chart2) + ") cannot " +
                "be computed by means of known changes of charts.")

        return self._coord_expression[(chart1, chart2)]

    def expr(self, chart1=None, chart2=None):
        r"""
        Return the expression of the differentiable mapping in terms of
        specified coordinates.

        If the expression is not already known, it is computed from some known
        expression by means of change-of-chart formulas.

        INPUT:

        - ``chart1`` -- (default: None) chart on the mapping's domain; if None,
          the domain's default chart is assumed
        - ``chart2`` -- (default: None) chart on the mapping's codomain; if
          None, the codomain's default chart is assumed

        OUTPUT:

        - symbolic expression representing the differentiable mapping in the
          above two charts

        EXAMPLES:

        Differential mapping from a 2-dimensional manifold to a 3-dimensional
        one::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: N = Manifold(3, 'N')
            sage: c_uv.<u,v> = M.chart()
            sage: c_xyz.<x,y,z> = N.chart()
            sage: Phi = M.diff_mapping(N, (u*v, u/v, u+v), name='Phi', latex_name=r'\Phi')
            sage: Phi.display()
            Phi: M --> N
               (u, v) |--> (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.expr(c_uv, c_xyz)
            (u*v, u/v, u + v)
            sage: Phi.expr()  # equivalent to above since 'uv' and 'xyz' are default charts
            (u*v, u/v, u + v)
            sage: type(Phi.expr()[0])
            <type 'sage.symbolic.expression.Expression'>

        Expressions in other charts::

            sage: c_UV.<U,V> = M.chart()  # new chart on M
            sage: ch_uv_UV = c_uv.coord_change(c_UV, u-v, u+v)
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ.<X,Y,Z> = N.chart() # new chart on N
            sage: ch_xyz_XYZ = c_xyz.coord_change(c_XYZ, 2*x-3*y+z, y+z-x, -x+2*y-z)
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.expr(c_UV, c_xyz)
            (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V)
            sage: Phi.expr(c_uv, c_XYZ)
            (((2*u + 1)*v^2 + u*v - 3*u)/v,
             -((u - 1)*v^2 - u*v - u)/v,
             -((u + 1)*v^2 + u*v - 2*u)/v)
            sage: Phi.expr(c_UV, c_XYZ)
             (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V), 1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V), 1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V))

        A rotation in some Euclidean plane::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M') # the plane (minus a segment to have global regular spherical coordinates)
            sage: c_spher.<r,ph> = M.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates on the plane
            sage: rot = M.diff_mapping(M, (r, ph+pi/3), name='R') # pi/3 rotation around r=0
            sage: rot.expr()
            (r, 1/3*pi + ph)

        Expression of the rotation in terms of Cartesian coordinates::

            sage: c_cart.<x,y> = M.chart() # Declaration of Cartesian coordinates
            sage: ch_spher_cart = c_spher.coord_change(c_cart, r*cos(ph), r*sin(ph)) # relation to spherical coordinates
            sage: ch_spher_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x))
            Check of the inverse coordinate transformation:
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
               x == x
               y == y
            sage: rot.expr(c_cart, c_cart)
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        """
        return self.multi_function_chart(chart1, chart2).expr()

    def set_expr(self, chart1, chart2, coord_functions):
        r"""
        Set a new coordinate representation of the mapping.

        The expressions with respect to other charts are deleted, in order to
        avoid any inconsistency. To keep them, use :meth:`add_expr` instead.

        INPUT:

        - ``chart1`` -- chart for the coordinates on the mapping's domain
        - ``chart2`` -- chart for the coordinates on the mapping's codomain
        - ``coord_functions`` -- the coordinate symbolic expression of the
          mapping in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single
          expression is expected (not a list with a single element)

        EXAMPLES:

        Polar representation of a planar rotation initally defined in
        Cartesian coordinates::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'R^2', r'\RR^2')  # the Euclidean plane R^2
            sage: c_xy.<x,y> = M.chart() # Cartesian coordinate on R^2
            sage: U = M.open_subset('U', coord_def={c_xy: (y!=0, x<0)}) # the complement of the segment y=0 and x>0
            sage: c_cart = c_xy.restrict(U) # Cartesian coordinates on U
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates on U
            sage: # Links between spherical coordinates and Cartesian ones:
            sage: ch_cart_spher = c_cart.coord_change(c_spher, sqrt(x*x+y*y), atan2(y,x))
            sage: ch_cart_spher.set_inverse(r*cos(ph), r*sin(ph))
            Check of the inverse coordinate transformation:
               x == x
               y == y
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
            sage: rot = U.diff_mapping(U, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
            sage: rot.display(c_cart, c_cart)
            R: U --> U
               (x, y) |--> (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        Let us use the method :meth:`set_expr` to set the
        spherical-coordinate expression by hand::

            sage: rot.set_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.display(c_spher, c_spher)
            R: U --> U
               (r, ph) |--> (r, 1/3*pi + ph)

        The expression in Cartesian coordinates has been erased::

            sage: rot._coord_expression
            {(chart (U, (r, ph)),
              chart (U, (r, ph))): functions (r, 1/3*pi + ph) on the chart (U, (r, ph))}

        It is recovered (thanks to the known change of coordinates) by a call
        to :meth:`display`::

            sage: rot.display(c_cart, c_cart)
            R: U --> U
               (x, y) |--> (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
            sage: rot._coord_expression  # random (dictionary output)
            {(chart (U, (x, y)),
              chart (U, (x, y))): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart (U, (x, y)),
             (chart (U, (r, ph)),
              chart (U, (r, ph))): functions (r, 1/3*pi + ph) on the chart (U, (r, ph))}

        """
        if self._is_identity:
            raise NotImplementedError("set_expr() must not be used for the " +
                                      "identity map")
        if chart1 not in self._domain._atlas:
            raise ValueError("The " + str(chart1) +
               " has not been defined on the " + str(self._domain))
        if chart2 not in self._codomain._atlas:
            raise ValueError("The " + str(chart2) +
              " has not been defined on the " + str(self._codomain))
        self._coord_expression.clear()
        self._del_derived()
        n2 = self._codomain._manifold._dim
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError(str(n2) +
                                 " coordinate function must be provided.")
            self._coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
        else:
            self._coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, coord_functions)

    def add_expr(self, chart1, chart2, coord_functions):
        r"""
        Set a new coordinate representation of the mapping.

        The previous expressions with respect to other charts are kept. To
        clear them, use :meth:`set_expr` instead.

        INPUT:

        - ``chart1`` -- chart for the coordinates on the mapping's domain
        - ``chart2`` -- chart for the coordinates on the mapping's codomain
        - ``coord_functions`` -- the coordinate symbolic expression of the
          mapping in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single
          expression is expected (not a list with a single element)

        .. WARNING::

            If the mapping has already expressions in other charts, it
            is the user's responsability to make sure that the expression
            to be added is consistent with them.

        EXAMPLES:

        Polar representation of a planar rotation initally defined in
        Cartesian coordinates::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'R^2', r'\RR^2')  # the Euclidean plane R^2
            sage: c_xy.<x,y> = M.chart() # Cartesian coordinate on R^2
            sage: U = M.open_subset('U', coord_def={c_xy: (y!=0, x<0)}) # the complement of the segment y=0 and x>0
            sage: c_cart = c_xy.restrict(U) # Cartesian coordinates on U
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates on U
            sage: # Links between spherical coordinates and Cartesian ones:
            sage: ch_cart_spher = c_cart.coord_change(c_spher, sqrt(x*x+y*y), atan2(y,x))
            sage: ch_cart_spher.set_inverse(r*cos(ph), r*sin(ph))
            Check of the inverse coordinate transformation:
               x == x
               y == y
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
            sage: rot = U.diff_mapping(U, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
            sage: rot.display(c_cart, c_cart)
            R: U --> U
               (x, y) |--> (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        If we make Sage calculate the expression in terms of spherical
        coordinates, via the method :meth:`display`, we notice some difficulties
        in arctan2 simplifications::

            sage: rot.display(c_spher, c_spher)
            R: U --> U
               (r, ph) |--> (r, arctan2(1/2*(sqrt(3)*cos(ph) + sin(ph))*r, -1/2*(sqrt(3)*sin(ph) - cos(ph))*r))

        Therefore, we use the method :meth:`add_expr` to set the
        spherical-coordinate expression by hand::

            sage: rot.add_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.display(c_spher, c_spher)
            R: U --> U
               (r, ph) |--> (r, 1/3*pi + ph)

        The call to :meth:`add_expr` has not deleted the expression in
        terms of Cartesian coordinates, as we can check by printing the
        dictionary :attr:`_coord_expression`, which stores the various internal
        representations of the differentiable mapping::

            sage: rot._coord_expression # random (dictionary output)
            {(chart (U, (x, y)),
              chart (U, (x, y))): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart (U, (x, y)),
             (chart (U, (r, ph)),
              chart (U, (r, ph))): functions (r, 1/3*pi + ph) on the chart (U, (r, ph))}

        If, on the contrary, we use :meth:`set_expr`, the expression in
        Cartesian coordinates is lost::

            sage: rot.set_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot._coord_expression
            {(chart (U, (r, ph)),
              chart (U, (r, ph))): functions (r, 1/3*pi + ph) on the chart (U, (r, ph))}

        It is recovered (thanks to the known change of coordinates) by a call
        to :meth:`display`::

            sage: rot.display(c_cart, c_cart)
            R: U --> U
               (x, y) |--> (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
            sage: rot._coord_expression  # random (dictionary output)
            {(chart (U, (x, y)),
              chart (U, (x, y))): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart (U, (x, y)),
             (chart (U, (r, ph)),
              chart (U, (r, ph))): functions (r, 1/3*pi + ph) on the chart (U, (r, ph))}

        The rotation can be applied to a point by means of either coordinate
        system::

            sage: p = M.point((1,2))  #  p defined by its Cartesian coord.
            sage: q = rot(p)  # q is computed by means of Cartesian coord.
            sage: p1 = M.point((sqrt(5), arctan(2)), chart=c_spher) # p1 is defined only in terms of c_spher
            sage: q1 = rot(p1) # computation by means of spherical coordinates
            sage: q1 == q
            True

        """
        if self._is_identity:
            raise NotImplementedError("add_expr() must not be used for the " +
                                      "identity map")
        if chart1 not in self._domain._atlas:
            raise ValueError("The " + str(chart1) +
               " has not been defined on the " + str(self._domain))
        if chart2 not in self._codomain._atlas:
            raise ValueError("The " + str(chart2) +
              " has not been defined on the " + str(self._codomain))
        n2 = self._codomain._manifold._dim
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError(str(n2) +
                                 " coordinate function must be provided.")
            self._coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
        else:
            self._coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, coord_functions)

    def restrict(self, subdomain, subcodomain=None):
        r"""
        Restriction of the differentiable mapping to some subdomain of its
        domain of definition.

        INPUT:

        - ``subdomain`` -- the subdomain of ``self._domain`` (instance of
          :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`)
        - ``subcodomain`` -- (default: None) subdomain of ``self._codomain``;
          if None, ``self._codomain`` is assumed.

        OUTPUT:

        - the restriction of ``self`` to ``dom``, as an instance of
          class :class:`DiffMapping`

        EXAMPLE:

        Restriction to an annulus of a diffeomorphism between the open unit
        disk and `\RR^2`::

            sage: M = Manifold(2, 'R^2')  # R^2
            sage: c_xy.<x,y> = M.chart()  # Cartesian coord. on R^2
            sage: D = M.open_subset('D', coord_def={c_xy: x^2+y^2<1}) # the open unit disk
            sage: Phi = D.diff_mapping(M, [x/sqrt(1-x^2-y^2), y/sqrt(1-x^2-y^2)], name='Phi', latex_name=r'\Phi')
            sage: Phi.display()
            Phi: D --> R^2
               (x, y) |--> (x, y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))
            sage: c_xy_D = c_xy.restrict(D)
            sage: U = D.open_subset('U', coord_def={c_xy_D: x^2+y^2>1/2}) # the annulus 1/2 < r < 1
            sage: Phi.restrict(U)
            differentiable mapping 'Phi' from the open subset 'U' of the
             2-dimensional manifold 'R^2' to the 2-dimensional manifold 'R^2'
            sage: Phi.restrict(U).parent()
            Set of Morphisms from open subset 'U' of the 2-dimensional manifold
             'R^2' to 2-dimensional manifold 'R^2' in Category of facade sets
            sage: Phi.domain()
            open subset 'D' of the 2-dimensional manifold 'R^2'
            sage: Phi.restrict(U).domain()
            open subset 'U' of the 2-dimensional manifold 'R^2'
            sage: Phi.restrict(U).display()
            Phi: U --> R^2
               (x, y) |--> (x, y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))

        The result is cached::

            sage: Phi.restrict(U) is Phi.restrict(U)
            True

        The restriction of the identity map::

            sage: id = D.identity_map() ; id
            identity map 'Id_D' of the open subset 'D' of the 2-dimensional manifold 'R^2'
            sage: id.restrict(U)
            identity map 'Id_U' of the open subset 'U' of the 2-dimensional manifold 'R^2'
            sage: id.restrict(U) is U.identity_map()
            True

        """
        from sage.categories.homset import Hom
        if subdomain == self._domain:
            return self
        if subcodomain is None:
            if self._is_identity:
                subcodomain = subdomain
            else:
                subcodomain = self._codomain
        if (subdomain, subcodomain) not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("The specified domain is not a subset " +
                                 "of the domain of definition of the diff. " +
                                 "mapping.")
            if not subcodomain.is_subset(self._codomain):
                raise ValueError("The specified codomain is not a subset " +
                                 "of the codomain of the diff. mapping.")
            # Special case of the identity map:
            if self._is_identity:
                self._restrictions[(subdomain, subcodomain)] = \
                                                       subdomain.identity_map()
                return self._restrictions[(subdomain, subcodomain)]
            # Generic case:
            homset = Hom(subdomain, subcodomain)
            resu = homset.element_class(homset, name=self._name,
                                        latex_name=self._latex_name)
            for charts in self._coord_expression:
                for ch1 in charts[0]._subcharts:
                    if ch1._domain.is_subset(subdomain):
                        for ch2 in charts[1]._subcharts:
                            if ch2._domain.is_subset(subcodomain):
                                for sch2 in ch2._supercharts:
                                    if (ch1, sch2) in resu._coord_expression:
                                        break
                                else:
                                    for sch2 in ch2._subcharts:
                                        if (ch1, sch2) in resu._coord_expression:
                                            del resu._coord_expression[(ch1, sch2)]
                                    coord_functions = \
                                          self._coord_expression[charts].expr()
                                    resu._coord_expression[(ch1, ch2)] = \
                                      MultiFunctionChart(ch1, *coord_functions)
            self._restrictions[(subdomain, subcodomain)] = resu
        return self._restrictions[(subdomain, subcodomain)]


    def pullback(self, tensor):
        r"""
        Pullback operator associated with the differentiable mapping.

        INPUT:

        - ``tensor`` -- instance of class
          :class:`~sage.geometry.manifolds.tensorfield.TensorField`
          representing a fully covariant tensor field `T` on the mapping's
          codomain, i.e. a tensor field of type (0,p), with p a positive or
          zero integer. The case p=0 corresponds to a scalar field.

        OUTPUT:

        - instance of class
          :class:`~sage.geometry.manifolds.tensorfield.TensorField`
          representing a fully covariant tensor field on the mapping's domain
          that is the pullback of `T` given by ``self``.

        EXAMPLES:

        Pullback on `S^2` of a scalar field defined on `R^3`::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') # the complement of a meridian (domain of spherical coordinates)
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coord. on U
            sage: N = Manifold(3, 'R^3', r'\RR^3', start_index=1)
            sage: c_cart.<x,y,z> = N.chart() # Cartesian coord. on R^3
            sage: Phi = U.diff_mapping(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
            sage: f = N.scalar_field(x*y*z, name='f') ; f
            scalar field 'f' on the 3-dimensional manifold 'R^3'
            sage: f.display()
            f: R^3 --> R
               (x, y, z) |--> x*y*z
            sage: pf = Phi.pullback(f) ; pf
            scalar field 'Phi_*(f)' on the open subset 'U' of the 2-dimensional manifold 'S^2'
            sage: pf.display()
            Phi_*(f): U --> R
               (th, ph) |--> cos(ph)*cos(th)*sin(ph)*sin(th)^2

        Pullback on `S^2` of the standard Euclidean metric on `R^3`::

            sage: g = N.sym_bilin_form_field('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: g.display()
            g = dx*dx + dy*dy + dz*dz
            sage: pg = Phi.pullback(g) ; pg
            field of symmetric bilinear forms 'Phi_*(g)' on the open subset 'U' of the 2-dimensional manifold 'S^2'
            sage: pg.display()
            Phi_*(g) = dth*dth + sin(th)^2 dph*dph

        Pullback on `S^2` of a 3-form on `R^3`::

            sage: a = N.diff_form(3, 'A')
            sage: a[1,2,3] = f
            sage: a.display()
            A = x*y*z dx/\dy/\dz
            sage: pa = Phi.pullback(a) ; pa
            3-form 'Phi_*(A)' on the open subset 'U' of the 2-dimensional manifold 'S^2'
            sage: pa.display() # should be zero (as any 3-form on a 2-dimensional manifold)
            Phi_*(A) = 0

        """
        from sage.geometry.manifolds.tensorfield import TensorFieldParal
        from sage.geometry.manifolds.chart import FunctionChart
        # Special case of the identity map:
        if self._is_identity:
            return tensor  # no test for efficiency
        # Generic case:
        dom1 = self._domain
        dom2 = self._codomain
        tdom = tensor._domain
        if not tdom.is_subset(dom2):
            raise TypeError("The tensor field is not defined on the mapping " +
                            "codomain.")
        (ncon, ncov) = tensor._tensor_type
        if ncon != 0:
            raise TypeError("The pullback cannot be taken on a tensor " +
                            "with some contravariant part.")
        resu_name = None ; resu_latex_name = None
        if self._name is not None and tensor._name is not None:
            resu_name = self._name + '_*(' + tensor._name + ')'
        if self._latex_name is not None and tensor._latex_name is not None:
            resu_latex_name = self._latex_name + '_*' + tensor._latex_name
        if ncov == 0:
            # Case of a scalar field
            resu_fc = []
            for chart2 in tensor._express:
                for chart1 in dom1._atlas:
                    if (chart1, chart2) in self._coord_expression:
                        phi = self._coord_expression[(chart1, chart2)]
                        coord1 = chart1._xx
                        ff = tensor._express[chart2]
                        resu_fc.append( FunctionChart(chart1,
                                                       ff(*(phi(*coord1)))) )
            dom_resu = resu_fc[0]._chart._domain
            for fc in resu_fc[1:]:
                dom_resu = dom_resu.union(fc._chart._domain)
            resu = dom_resu.scalar_field(name=resu_name,
                                         latex_name=resu_latex_name)
            for fc in resu_fc:
                resu._express[fc._chart] = fc
        else:
            # Case of tensor field of rank >= 1
            if tensor._vmodule._dest_map is not tdom._identity_map:
                raise TypeError("The pullback in defined only for tensors " +
                                "on " + str(dom2) + ".")
            resu_rst = []
            for chart_pair in self._coord_expression:
                chart1 = chart_pair[0] ; chart2 = chart_pair[1]
                ch2dom = chart2._domain
                if ch2dom.is_subset(tdom):
                    self_r = self.restrict(chart1._domain, subcodomain=ch2dom)
                    tensor_r = tensor.restrict(ch2dom)
                    resu_rst.append( self_r._pullback_paral(tensor_r) )
            dom_resu = resu_rst[0]._domain
            for rst in resu_rst[1:]:
                dom_resu = dom_resu.union(rst._domain)
            resu = dom_resu.tensor_field(0, ncov, name=resu_name,
                                         latex_name=resu_latex_name,
                                         sym=resu_rst[0]._sym,
                                         antisym=resu_rst[0]._antisym)
            for rst in resu_rst:
                if rst._domain is not resu._domain:
                    resu._restrictions[rst._domain] = rst
            if isinstance(resu, TensorFieldParal):
                for rst in resu_rst:
                    for frame, comp in rst._components.iteritems():
                        resu._components[frame] = comp
        return resu

    def _pullback_paral(self, tensor):
        r"""
        Pullback on parallelizable domains.
        """
        from vectorframe import CoordFrame
        from sage.tensor.modules.comp import Components, CompWithSym, \
                                                 CompFullySym, CompFullyAntiSym
        dom1 = self._domain
        dom2 = self._codomain
        ncov = tensor._tensor_type[1]
        resu_name = None ; resu_latex_name = None
        if self._name is not None and tensor._name is not None:
            resu_name = self._name + '_*(' + tensor._name + ')'
        if self._latex_name is not None and tensor._latex_name is not None:
            resu_latex_name = self._latex_name + '_*' + tensor._latex_name
        fmodule1 = dom1.vector_field_module()
        ring1 = fmodule1._ring
        si1 = fmodule1._sindex
        of1 = fmodule1._output_formatter
        si2 = dom2._manifold._sindex
        resu = fmodule1.tensor((0,ncov), name=resu_name,
                               latex_name=resu_latex_name, sym=tensor._sym,
                               antisym=tensor._antisym)
        for frame2 in tensor._components:
            if isinstance(frame2, CoordFrame):
                chart2 = frame2._chart
                for chart1 in dom1._atlas:
                    if (chart1, chart2) in self._coord_expression:
                        # Computation at the component level:
                        frame1 = chart1._frame
                        tcomp = tensor._components[frame2]
                        if isinstance(tcomp, CompFullySym):
                            ptcomp = CompFullySym(ring1, frame1, ncov,
                                                  start_index=si1,
                                                  output_formatter=of1)
                        elif isinstance(tcomp, CompFullyAntiSym):
                            ptcomp = CompFullyAntiSym(ring1, frame1, ncov,
                                                      start_index=si1,
                                                      output_formatter=of1)
                        elif isinstance(tcomp, CompWithSym):
                            ptcomp = CompWithSym(ring1, frame1, ncov,
                                                 start_index=si1,
                                                 output_formatter=of1,
                                                 sym=tcomp.sym,
                                                 antisym=tcomp.antisym)
                        else:
                            ptcomp = Components(ring1, frame1, ncov,
                                                start_index=si1,
                                                output_formatter=of1)
                        phi = self._coord_expression[(chart1, chart2)]
                        jacob = phi.jacobian()
                        # X2 coordinates expressed in terms of X1 ones via the
                        # mapping:
                        coord2_1 = phi(*(chart1._xx))
                        for ind_new in ptcomp.non_redundant_index_generator():
                            res = 0
                            for ind_old in dom2._manifold.index_generator(ncov):
                                ff = tcomp[[ind_old]].function_chart(chart2)
                                t = FunctionChart(chart1, ff(*coord2_1))
                                for i in range(ncov):
                                    t *= jacob[ind_old[i]-si2][ind_new[i]-si1]
                                res += t
                            ptcomp[ind_new] = res
                        resu._components[frame1] = ptcomp
            return resu

    def __invert__(self, chart1=None, chart2=None):
        r"""
        Return the inverse of ``self`` if ``self`` is a diffeomorphism.

        INPUT:

        - ``chart1`` -- (default: None) chart in which the computation of the
          inverse is performed if necessary; if none is provided, the default
          chart of the domain of ``self`` will be used
        - ``chart2`` -- (default: None) chart in which the computation of the
          inverse is performed if necessary; if none is provided, the default
          chart of the codomain of ``self`` will be used

        OUTPUT:

        - the inverse diffeomorphism

        EXAMPLES:

        The inverse of a rotation in the Euclidean plane::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'R^2', r'\RR^2')
            sage: c_cart.<x,y> = M.chart()
            sage: # A pi/3 rotation around the origin:
            sage: rot = M.diffeomorphism(M, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
            sage: rot.inverse()
            diffeomorphism 'R^(-1)' of the 2-dimensional manifold 'R^2'
            sage: rot.inverse().display()
            R^(-1): R^2 --> R^2
               (x, y) |--> (1/2*sqrt(3)*y + 1/2*x, -1/2*sqrt(3)*x + 1/2*y)

        Checking that applying successively the diffeomorphism and its
        inverse results in the identity::

            sage: (a, b) = var('a b')
            sage: p = M.point((a,b)) # a generic point on M
            sage: q = rot(p)
            sage: p1 = rot.inverse()(q)
            sage: p1 == p
            True

        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        from sage.geometry.manifolds.utilities import simplify_chain
        if self._inverse is not None:
            return self._inverse
        if not self._is_diffeo:
            raise ValueError("the {} is not a diffeomorphism".format(self))
        if chart1 is None: chart1 = self._domain._def_chart
        if chart2 is None: chart2 = self._codomain._def_chart
        coord_map = self._coord_expression[(chart1, chart2)]
        n1 = len(chart1._xx)
        n2 = len(chart2._xx)
        # New symbolic variables (different from chart2._xx to allow for a
        #  correct solution even when chart2 = chart1):
        x2 = [ SR.var('xxxx' + str(i)) for i in range(n2) ]
        equations = [ x2[i] == coord_map._functions[i]._express
                      for i in range(n2) ]
        solutions = solve(equations, chart1._xx, solution_dict=True)
        if len(solutions) == 0:
            raise ValueError("No solution found")
        if len(solutions) > 1:
            raise ValueError("Non-unique solution found")
        #!# This should be the Python 2.7 form:
        # substitutions = {x2[i]: chart2._xx[i] for i in range(n2)}
        #
        # Here we use a form compatible with Python 2.6:
        substitutions = dict([(x2[i], chart2._xx[i]) for i in range(n2)])

        inv_functions = [solutions[0][chart1._xx[i]].subs(substitutions)
                           for i in range(n1)]
        for i in range(n1):
            x = inv_functions[i]
            try:
                inv_functions[i] = simplify_chain(x)
            except AttributeError:
                pass
        if self._name is None:
            name = None
        else:
            name = self._name + '^(-1)'

        if self._latex_name is None:
            latex_name = None
        else:
            latex_name = self._latex_name + r'^{-1}'
        self._inverse = self._codomain.diffeomorphism(self._domain,
                               coord_functions=inv_functions, chart1=chart2,
                               chart2=chart1, name=name, latex_name=latex_name)
        return self._inverse

    inverse = __invert__

    def differential(self, point):
        r"""
        Return the differential of ``self`` at a given point.

        If ``self`` is the differentiable mapping

        .. MATH::

            \Phi: U\subset M \longrightarrow N

        where `M` and `N` are differentiable manifolds and `U` is an open
        subset of `M`, the *differential* of `\Phi` at a point `p\in U` is the
        tangent space linear map:

        .. MATH::

            \mathrm{d}\Phi_p: T_p M \longrightarrow T_{\Phi(p)} N

        defined by

        .. MATH::

            \begin{array}{rccc}
            \forall v\in T_p M,\quad \mathrm{d}\Phi_p(v) : & C^\infty(N) &
                                                \longrightarrow & \mathbb{R} \\
                                & f & \longmapsto & v(f\circ \Phi)
            \end{array}

        INPUT:

        - ``point`` -- point `p` in the domain of ``self``

        OUTPUT:

        - `\mathrm{d}\Phi_p`, the differential of ``self`` at `p`, as an
          instance of
          :class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`

        EXAMPLES:

        Differential of a mapping between a 2-dimensional manifold and a
        3-dimensional one::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: Phi = M.diff_mapping(N, {(X,Y): (x-2*y, x*y, x^2-y^3)}, name='Phi',
            ....:                      latex_name = r'\Phi')
            sage: p = M.point((2,-1), name='p')
            sage: dPhip = Phi.differential(p) ; dPhip
            Generic morphism:
              From: tangent space at point 'p' on 2-dimensional manifold 'M'
              To:   tangent space at point 'Phi(p)' on 3-dimensional manifold 'N'
            sage: latex(dPhip)
            \mathrm{d}\Phi_{p}
            sage: dPhip.parent()
            Set of Morphisms from tangent space at point 'p' on 2-dimensional
             manifold 'M' to tangent space at point 'Phi(p)' on 3-dimensional
             manifold 'N' in Category of vector spaces over Symbolic Ring

        The matrix of `\mathrm{d}\Phi_p` w.r.t. to the default bases of
        `T_p M` and `T_{\Phi(p)} N`::

            sage: dPhip.matrix()
            [ 1 -2]
            [-1  2]
            [ 4 -3]

        """
        image_point = self(point)
        tsp_image = image_point.tangent_space()
        tsp_source = point.tangent_space()
        # Search for a common chart to perform the computation
        chartp = None
        # 1/ Search without any extra computation
        for chart in point._coordinates:
            for chart_pair in self._diff:
                if chart == chart_pair[0]:
                    chartp = chart_pair
                    break
            if chartp is not None:
                break
        else:
            # 2/ Search with a coordinate transformation on the point
            for chart_pair in self._diff:
                try:
                    point.coord(chart_pair[0])
                    chartp = chart_pair
                except ValueError:
                    pass
        if chartp is None:
            # 3/ Search with a coordinate evaluation of self
            for chart1 in point._coordinates:
                for chart2 in self._codomain.atlas():
                    try:
                        self.differential_functions(chart1, chart2)
                        chartp = (chart1, chart2)
                        break
                    except ValueError:
                        pass
                if chartp is not None:
                    break
        if chartp is None:
            raise ValueError("no common chart have been found for the " +
                     "coordinate expressions of {} and {}".format(self, point))
        diff_funct = self.differential_functions(*chartp)
        chart1 = chartp[0]
        chart2 = chartp[1]
        coord_point = point.coord(chart1)
        n1 = self._domain._manifold.dim()
        n2 = self._codomain._manifold.dim()
        matrix = [[diff_funct[i][j](*coord_point) for j in range(n1)]
                                                            for i in range(n2)]
        bases = (chart1.frame().at(point), chart2.frame().at(image_point))
        if self._name is not None and point._name is not None:
            name = 'd' + self._name + '_' + point._name
        else:
            name = None
        if self._latex_name is not None and point._latex_name is not None:
            latex_name = r'\mathrm{d}' + self._latex_name + r'_{' + \
                         point._latex_name + '}'
        else:
            latex_name = None
        return tsp_source.hom(tsp_image, matrix, bases=bases, name=name,
                              latex_name=latex_name)

    def differential_functions(self, chart1=None, chart2=None):
        r"""
        Return the coordinate expression of the differential of ``self``
        w.r.t. a pair of charts.

        If ``self`` is the differentiable mapping

        .. MATH::

            \Phi: U\subset M \longrightarrow V\subset N

        where `U` and `V` are two open subsets of the differentiable manifolds
        `M` and `N`, the *differential* of `\Phi` at a point `p\in U` is the
        tangent space linear map:

        .. MATH::

            \mathrm{d}\Phi_p: T_p M \longrightarrow T_{\Phi(p)} N

        defined by

        .. MATH::

            \begin{array}{rccc}
            \forall v\in T_p M,\quad \mathrm{d}\Phi_p(v) : & C^\infty(N) &
                                                \longrightarrow & \mathbb{R} \\
                                & f & \longmapsto & v(f\circ \Phi)
            \end{array}

        If the coordinate expression of `\Phi` is

        .. MATH::

            y^i = Y^i(x^1,\ldots,x^n)  \quad 1\leq i \leq m

        where $(x^1,\ldots,x^n)$ are coordinates of a chart on `U` and
        $(y^1,\ldots,y^m)$ are coordinates of a chart on `V`, the expression
        of the differential of `\Phi` w.r.t to these coordinates is

        .. MATH::

            J_{ij} = \frac{\partial Y^i}{\partial x^j} \quad 1\leq i \leq m,
                            \quad 1\leq j \leq n

        `\left. J_{ij} \right|_p` is then the matrix of the linear map
        `\mathrm{d}\Phi_p` with respect to the bases of `T_p M` and
        `T_{\Phi(p)} N` associated to the above charts:

        .. MATH::

            \mathrm{d}\Phi_p\left(  \left. \frac{\partial}{\partial x^j} \right| _p
                    \right) = \left. J_{ij} \right|_p \;
             \left. \frac{\partial}{\partial y^i} \right| _{\Phi(p)}

        INPUT:

        - ``chart1`` -- (default: ``None``) chart on the domain of ``self``
          (coordinates denoted by `(x^j)` above); if none is provided, the
          domain's default chart is assumed
        - ``chart2`` -- (default: ``None``) chart on the codomain of ``self``
          (coordinates denoted by `(y^i)` above); if none is provided, the
          codomain's default chart is assumed

        OUTPUT:

        - the functions `J_{ij}` as a double array, `J_{ij}` being the element
          ``[i][j]``, represented by an instance of
          :class:`~sage.geometry.manifolds.chart.FunctionChart`.
          To get symbolic expressions, use the method
          :meth:`jacobian_matrix` instead.

        EXAMPLES:

        Differential functions of a mapping between a 2-dimensional manifold
        and a 3-dimensional one::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: Phi = M.diff_mapping(N, {(X,Y): (x-2*y, x*y, x^2-y^3)}, name='Phi',
            ....:                      latex_name = r'\Phi')
            sage: J = Phi.differential_functions(X, Y) ; J
            [[1, -2], [y, x], [2*x, -3*y^2]]

        The elements of ``J`` are functions of the coordinates of chart ``X``::

            sage: J[2][0]
            2*x
            sage: type(J[2][0])
            <class 'sage.geometry.manifolds.chart.FunctionChart'>
            sage: J[2][0].display()
            (x, y) |--> 2*x

        In contrast, the method :meth:`jacobian_matrix` leads directly to
        symbolic expressions::

            sage: JJ = Phi.jacobian_matrix(X,Y) ; JJ
            [     1     -2]
            [     y      x]
            [   2*x -3*y^2]
            sage: JJ[2,0]
            2*x
            sage: type(JJ[2,0])
            <type 'sage.symbolic.expression.Expression'>
            sage: bool( JJ[2,0] == J[2][0].expr() )
            True

        """
        dom1 = self._domain; dom2 = self._codomain
        if chart1 is None:
            chart1 = dom1._def_chart
        if chart2 is None:
            chart2 = dom2._def_chart
        if (chart1, chart2) not in self._diff:
            # Some computation must be performed
            manif1 = dom1._manifold
            n2 = dom2._manifold.dim()
            funct = self.multi_function_chart(chart1, chart2)
            self._diff[(chart1, chart2)] = [[funct[i].diff(j) for j in
                                           manif1.irange()] for i in range(n2)]
        return self._diff[(chart1, chart2)]

    def jacobian_matrix(self, chart1=None, chart2=None):
        r"""
        Return the Jacobian matrix resulting from the coordinate expression of
        ``self`` w.r.t. a pair of charts.

        If the coordinate expression of ``self`` is

        .. MATH::

            y^i = Y^i(x^1,\ldots,x^n)  \quad 1\leq i \leq m

        where $(x^1,\ldots,x^n)$ are coordinates of a chart `X` on the
        domain of ``self`` and $(y^1,\ldots,y^m)$ are coordinates of a chart
        `Y` on the codomain of ``self``, the *Jacobian matrix* of the
        differentiable mapping ``self`` w.r.t. to charts `X` and `Y` is

        .. MATH::

            J = \left( \frac{\partial Y^i}{\partial x^j}
              \right) _{{1\leq i \leq m\atop 1\leq j \leq n}},

        where `i` is the row index and `j` the column one.

        INPUT:

        - ``chart1`` -- (default: ``None``) chart `X` on the domain of
          ``self``; if none is provided, the domain's default chart is assumed
        - ``chart2`` -- (default: ``None``) chart `Y` on the codomain of
          ``self``; if none is provided, the codomain's default chart is
          assumed

        OUTPUT:

        - the matrix `J` defined above

        EXAMPLES:

        Jacobian matrix of a mapping between a 2-dimensional manifold
        and a 3-dimensional one::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: Phi = M.diff_mapping(N, {(X,Y): (x-2*y, x*y, x^2-y^3)}, name='Phi',
            ....:                      latex_name = r'\Phi')
            sage: Phi.display()
            Phi: M --> N
               (x, y) |--> (u, v, w) = (x - 2*y, x*y, -y^3 + x^2)
            sage: J = Phi.jacobian_matrix(X, Y) ; J
            [     1     -2]
            [     y      x]
            [   2*x -3*y^2]
            sage: J.parent()
            Full MatrixSpace of 3 by 2 dense matrices over Symbolic Ring

        """
        from sage.matrix.constructor import matrix
        diff_funct = self.differential_functions(chart1, chart2)
        n1 = self._domain._manifold.dim()
        n2 = self._codomain._manifold.dim()
        return matrix( [[diff_funct[i][j].expr() for j in range(n1)]
                                                          for i in range(n2)] )
