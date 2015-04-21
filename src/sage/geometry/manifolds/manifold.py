r"""
Differentiable manifolds

The class :class:`Manifold` implements differentiable manifolds over `\RR`.

Ideally this class should inherit from a class describing topological
manifolds or at least topological spaces. Since such classes do not
exist in Sage yet, the class :class:`Manifold` inherits from the generic Sage
class :class:`~sage.structure.parent.Parent` (via the
class :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`)
and is declared to belong to the category of sets (Sage category
:class:`~sage.categories.sets_cat.Sets`).
The corresponding Sage :class:`~sage.structure.element.Element`'s are
implemented via the class :class:`~sage.geometry.manifolds.point.ManifoldPoint`.

The derived class :class:`RealLine` implements the field of real numbers
`\RR` as a manifold of dimension one, while the class
:class:`~sage.geometry.manifolds.manifold.OpenInterval` implements
open intervals of it.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version

REFERENCES:

- M. Berger & B. Gostiaux: *Geometrie differentielle, varietes, courbes et
  surfaces*, Presses Universitaires de France (Paris) (1987)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)

EXAMPLES:

    The sphere `S^2` as a 2-dimensional manifold::

        sage: M = Manifold(2, 'S^2')
        sage: M
        2-dimensional manifold 'S^2'
        sage: dim(M)
        2

    Let us consider the complement of the North pole; it is an open subset
    of `S^2`, which we call U::

        sage: U = M.open_subset('U') ; U
        open subset 'U' of the 2-dimensional manifold 'S^2'

    A standard chart on U is provided by the stereographic projection from the
    North pole to the equatorial plane::

        sage: stereoN.<x,y> = U.chart() ; stereoN
        chart (U, (x, y))

    Thanks to the operator ``<x,y>`` on the left-hand side, the coordinates
    declared in a chart (here x and y), are accessible by their names; they are
    Sage's symbolic variables::

        sage: y
        y
        sage: type(y)
        <type 'sage.symbolic.expression.Expression'>

    The South pole is the point of coordinates `(x,y)=(0,0)` in the above
    chart::

        sage: S = U.point((0,0), name='S') ; S
        point 'S' on 2-dimensional manifold 'S^2'

    Let us call V the subset that is the complement of the South pole and let
    us introduce on it the chart induced by the stereographic projection from
    the South pole to the equatorial plane::

        sage: V = M.open_subset('V') ; V
        open subset 'V' of the 2-dimensional manifold 'S^2'
        sage: stereoS.<u,v> = V.chart() ; stereoS
        chart (V, (u, v))

    The North pole is the point of coordinates `(u,v)=(0,0)` in this chart::

        sage: N = V.point((0,0), name='N') ; N
        point 'N' on 2-dimensional manifold 'S^2'

    To fully construct the manifold, we declare that it is the union of U
    and V::

        sage: M.declare_union(U,V)

    At this stage, the manifold's atlas contains two charts::

        sage: M.atlas()
        [chart (U, (x, y)), chart (V, (u, v))]

    To finalize things, we must declare the transition map between these two
    charts: calling W the intersection of U and V, (W is the subset of U
    defined by `x^2+y^2\not=0`, as well as the subset of V defined by
    `u^2+v^2\not=0`), we set::

        sage: transf = stereoN.transition_map(stereoS, (x/(x^2+y^2), y/(x^2+y^2)), \
                        intersection_name='W', restrictions1= x^2+y^2!=0, restrictions2= u^2+v^2!=0)
        sage: transf
        coordinate change from chart (W, (x, y)) to chart (W, (u, v))
        sage: W = U.intersection(V)
        sage: W.atlas()
        [chart (W, (x, y)), chart (W, (u, v))]
        sage: stereoN_W = W.atlas()[0]
        sage: stereoS_W = W.atlas()[1]

    The inverse of the transition map is computed by the method inverse()::

        sage: transf.inverse()(u,v)
        (u/(u^2 + v^2), v/(u^2 + v^2))

    At this stage, we have four open subsets on `S^2`::

        sage: M.list_of_subsets()
        [2-dimensional manifold 'S^2',
         open subset 'U' of the 2-dimensional manifold 'S^2',
         open subset 'V' of the 2-dimensional manifold 'S^2',
         open subset 'W' of the 2-dimensional manifold 'S^2']

    W is the open subset that is the complement of the two poles::

        sage: N in W
        False
        sage: S in W
        False

    The North pole lies in `V` and the South pole in `U`::

        sage: N in V, N in U
        (True, False)
        sage: S in U, S in V
        (True, False)

    Four charts have been defined on the manifold::

        sage: M.atlas()
        [chart (U, (x, y)), chart (V, (u, v)), chart (W, (x, y)), chart (W, (u, v))]

    The first defined chart is considered as the default chart on the
    manifold (unless it is changed by the method
    :meth:`~sage.geometry.manifolds.domain.ManifoldSubset.set_default_chart`)::

        sage: M.default_chart()
        chart (U, (x, y))

    Being the *default chart* means that its mention can be omitted when
    specifying some point coordinates::

        sage: p = M.point((1,2), name='p')  # a point is created with coordinates (1,2) in the default chart
        sage: p = M.point((1,2), chart=stereoN, name='p') # the full declaration, equivalent to the above one
        sage: p._coordinates # random (dictionary output):
        {chart (W, (x, y)): (1, 2), chart (U, (x, y)): (1, 2)}
        sage: p.coord() # if the chart is not specified, the default chart coordinates are returned:
        (1, 2)
        sage: p.coord(stereoS_W) # the coordinates in the chart stereoS_W are computed by means of the transition map:
        (1/5, 2/5)

    Manifolds are Sage *parent* objects, whose *elements* are points::

        sage: isinstance(M, Parent)
        True
        sage: M.category()
        Category of sets
        sage: p.parent()
        2-dimensional manifold 'S^2'
        sage: M.is_parent_of(p)
        True
        sage: p in M
        True
        sage: p == M((1,2))
        True

    The tangent vector space at point p::

        sage: Tp = p.tangent_space() ; Tp
        tangent space at point 'p' on 2-dimensional manifold 'S^2'
        sage: Tp.category()
        Category of vector spaces over Symbolic Ring
        sage: dim(Tp)
        2

    A scalar field on the sphere::

        sage: f = M.scalar_field({stereoN: atan(x^2+y^2), stereoS: pi/2-atan(u^2+v^2)}, name='f')
        sage: f
        scalar field 'f' on the 2-dimensional manifold 'S^2'
        sage: f.display()
        f: S^2 --> R
        on U: (x, y) |--> arctan(x^2 + y^2)
        on V: (u, v) |--> 1/2*pi - arctan(u^2 + v^2)
        sage: f(p)
        arctan(5)
        sage: f.parent()
        algebra of scalar fields on the 2-dimensional manifold 'S^2'
        sage: f.parent().category()
        Category of commutative algebras over Symbolic Ring

    A manifold has a default vector frame, which, unless otherwise specified,
    is the coordinate frame associated with the first defined chart::

        sage: M.default_frame()
        coordinate frame (U, (d/dx,d/dy))
        sage: latex(M.default_frame())
        \left(U ,\left(\frac{\partial}{\partial x },\frac{\partial}{\partial y }\right)\right)
        sage: M.default_frame() is stereoN.frame()
        True

    A vector field on the manifold::

        sage: w = M.vector_field('w')
        sage: w[stereoN.frame(), :] = [x, y]
        sage: w.add_comp_by_continuation(stereoS.frame(), W, stereoS)
        sage: w.display() # display in the default frame (stereoN.frame())
        w = x d/dx + y d/dy
        sage: w.display(stereoS.frame())
        w = -u d/du - v d/dv
        sage: w.parent()
        module X(S^2) of vector fields on the 2-dimensional manifold 'S^2'
        sage: w.parent().category()
        Category of modules over algebra of scalar fields on the 2-dimensional manifold 'S^2'

    Vector fields act on scalar fields::

        sage: w(f)
        scalar field 'w(f)' on the 2-dimensional manifold 'S^2'
        sage: w(f).display()
        w(f): S^2 --> R
        on U: (x, y) |--> 2*(x^2 + y^2)/(x^4 + 2*x^2*y^2 + y^4 + 1)
        on V: (u, v) |--> 2*(u^2 + v^2)/(u^4 + 2*u^2*v^2 + v^4 + 1)
        sage: w(f) == f.differential()(w)
        True

    The value of the vector field at point p::

        sage: w.at(p)
        tangent vector w at point 'p' on 2-dimensional manifold 'S^2'
        sage: w.at(p).display()
        w = d/dx + 2 d/dy
        sage: w.at(p).parent()
        tangent space at point 'p' on 2-dimensional manifold 'S^2'

    A 1-form on the manifold::

        sage: df = f.differential() ; df
        1-form 'df' on the 2-dimensional manifold 'S^2'
        sage: df.display()
        df = 2*x/(x^4 + 2*x^2*y^2 + y^4 + 1) dx + 2*y/(x^4 + 2*x^2*y^2 + y^4 + 1) dy
        sage: df.display(stereoS.frame())
        df = -2*u/(u^4 + 2*u^2*v^2 + v^4 + 1) du - 2*v/(u^4 + 2*u^2*v^2 + v^4 + 1) dv
        sage: df.parent()
        Module /\^1(S^2) of 1-forms on the 2-dimensional manifold 'S^2'
        sage: df.parent().category()
        Category of modules over algebra of scalar fields on the 2-dimensional
         manifold 'S^2'

    The value of the 1-form at point p::

        sage: df.at(p)
        Linear form df on the tangent space at point 'p' on 2-dimensional
         manifold 'S^2'
        sage: df.at(p).display()
        df = 1/13 dx + 2/13 dy
        sage: df.at(p).parent()
        Dual of the tangent space at point 'p' on 2-dimensional manifold 'S^2'

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

from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.infinity import infinity, minus_infinity
from sage.symbolic.ring import SR
from domain import ManifoldOpenSubset

class Manifold(ManifoldOpenSubset):
    r"""
    Differentiable manifold over `\RR`.

    Ideally this class
    should inherit from a class describing topological manifolds, or at
    least, topological spaces. Since such classes do not exist in Sage yet,
    the class :class:`Manifold` inherits from the generic Sage
    class :class:`~sage.structure.parent.Parent` (via the
    class :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`)
    and is declared to belong to the category of sets (Sage category
    :class:`~sage.categories.sets_cat.Sets`). Accordingly, the class
    :class:`Manifold` is a Sage *parent* class, the corresponding *element*
    class being :class:`~sage.geometry.manifolds.point.ManifoldPoint`.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``latex_name`` -- (default: None) string; LaTeX symbol to denote the
      manifold; if none is provided, it is set to ``name``
    - ``start_index`` -- (default: 0) integer; lower bound of the range of
      indices used for "indexed objects" on the manifold, e.g. coordinates
      in a chart or elements of a vector frame and the corresponding tensor
      components.

    EXAMPLES:

    A 4-dimensional manifold::

        sage: M = Manifold(4, 'M', latex_name=r'\mathcal{M}')
        sage: M
        4-dimensional manifold 'M'
        sage: latex(M)
        \mathcal{M}

    The input parameter ``start_index`` defines the range of indices on the
    manifold::

        sage: M = Manifold(4, 'M')
        sage: list(M.irange())
        [0, 1, 2, 3]
        sage: M = Manifold(4, 'M', start_index=2)
        sage: list(M.irange())
        [2, 3, 4, 5]

    A manifold is a Sage *parent* object, in the category of sets::

        sage: isinstance(M, Parent)
        True
        sage: M.category()
        Category of sets
        sage: M in Sets()
        True

    The corresponding Sage *elements* are points::

        sage: X.<t, x, y, z> = M.chart()
        sage: p = M.an_element(); p
        point on 4-dimensional manifold 'M'
        sage: p.parent()
        4-dimensional manifold 'M'
        sage: M.is_parent_of(p)
        True
        sage: p in M
        True

    The manifold's points are instances of class
    :class:`~sage.geometry.manifolds.point.ManifoldPoint`::

        sage: isinstance(p, sage.geometry.manifolds.point.ManifoldPoint)
        True

    A manifold has a predefined zero scalar field, mapping all the points
    to 0::

        sage: M.zero_scalar_field()
        scalar field 'zero' on the 4-dimensional manifold 'M'
        sage: M.zero_scalar_field()(p)
        0

    The manifold passes all the tests of the test suite relative to the
    category of Sets::

        sage: TestSuite(M).run(verbose=True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
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


    """
    def __init__(self, n, name, latex_name=None, start_index=0):
        from sage.rings.integer import Integer
        if not isinstance(n, (int, Integer)):
            raise TypeError("The manifold dimension must be an integer.")
        if n<1:
            raise ValueError("The manifold dimension must be strictly " +
                             "positive.")
        self._dim = n
        ManifoldOpenSubset.__init__(self, self, name, latex_name)
        self._sindex = start_index
        #!# self._subsets = [self]

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        return str(self._dim) + "-dimensional manifold '%s'" % self._name

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        return self._latex_name

    def dimension(self):
        r"""
        Return the dimension of the manifold.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: M.dimension()
            2

        A shortcut is ``dim()``::

            sage: M.dim()
            2

        The Sage global function ``dim`` can also be used::

            sage: dim(M)
            2

        """
        return self._dim

    dim = dimension

    def irange(self, start=None):
        r"""
        Single index generator.

        INPUT:

        - ``start`` -- (default: None) initial value of the index; if none is
          provided, ``self._sindex`` is assumed

        OUTPUT:

        - an iterable index, starting from ``start`` and ending at
          ``self._sindex + self._dim -1``

        EXAMPLES:

        Index range on a 4-dimensional manifold::

            sage: M = Manifold(4, 'M')
            sage: for i in M.irange():
            ...       print i,
            ...
            0 1 2 3
            sage: for i in M.irange(2):
            ...       print i,
            ...
            2 3
            sage: list(M.irange())
            [0, 1, 2, 3]

        Index range on a 4-dimensional manifold with starting index=1::

            sage: M = Manifold(4, 'M', start_index=1)
            sage: for i in M.irange():
            ...       print i,
            ...
            1 2 3 4
            sage: for i in M.irange(2):
            ...      print i,
            ...
            2 3 4

        """
        si = self._sindex
        imax = self._dim + si
        if start is None:
            i = si
        else:
            i = start
        while i < imax:
            yield i
            i += 1


    def index_generator(self, nb_indices):
        r"""
        Generator of index series.

        INPUT:

        - ``nb_indices`` -- number of indices in a series

        OUTPUT:

        - an iterable index series for a generic component with the specified
          number of indices

        EXAMPLES:

        Indices on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: for ind in M.index_generator(2):
            ...       print ind
            ...
            (1, 1)
            (1, 2)
            (2, 1)
            (2, 2)

        Loops can be nested::

            sage: for ind1 in M.index_generator(2):
            ...       print ind1, " : ",
            ...       for ind2 in M.index_generator(2):
            ...           print ind2,
            ...       print ""
            ...
            (1, 1)  :  (1, 1) (1, 2) (2, 1) (2, 2)
            (1, 2)  :  (1, 1) (1, 2) (2, 1) (2, 2)
            (2, 1)  :  (1, 1) (1, 2) (2, 1) (2, 2)
            (2, 2)  :  (1, 1) (1, 2) (2, 1) (2, 2)

        """
        si = self._sindex
        imax = self._dim - 1 + si
        ind = [si for k in range(nb_indices)]
        ind_end = [si for k in range(nb_indices)]
        ind_end[0] = imax+1
        while ind != ind_end:
            yield tuple(ind)
            ret = 1
            for pos in range(nb_indices-1,-1,-1):
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = si
                        ret = 1

    def submanifold(self, dim, name, latex_name=None, start_index=0):
        r"""
        Construct a submanifold of ``self``

        See class :class:`~sage.geometry.manifolds.submanifolds.Submanifold`
        for a complete documentation.

        INPUT:

        - ``dim`` -- dimension of the submanifold
        - ``name`` -- name given to the submanifold
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the
          submanifold
        - ``start_index`` -- (default: 0) lower bound of the range of indices
          on the submanifold

        OUTPUT:

        - instance of class
          :class:`~sage.geometry.manifolds.submanifolds.Submanifold`

        """
        from submanifold import Submanifold
        return Submanifold(self, dim, name, latex_name=latex_name,
                           start_index=start_index)

#******************************************************************************

class RealLine(Manifold):
    r"""
    Field of real numbers, as a manifold of dimension 1 (real line) with a
    canonical coordinate chart.

    INPUT:

    - ``coordinate`` -- (default: ``None``) string defining the symbol of the
      canonical coordinate set on the real line; if none is provided and
      ``names`` is ``None``, the symbol 't' is used
    - ``name`` -- (default: 'R') name given to the real line
    - ``latex_name`` -- (default: r'\\RR') LaTeX symbol to denote the real line
    - ``start_index`` -- (default: 0) unique value of the index for vectors and
      forms on the real line.
    - ``names`` -- (default: ``None``) used only when ``coordinate`` is
      ``None``: it must be a single-element tuple containing the canonical
      coordinate symbol (this is guaranted if the shortcut operator ``<>`` is
      used, see examples below).

    EXAMPLES:

    Constructing the real line without any argument::

        sage: R = RealLine() ; R
        field R of real numbers
        sage: latex(R)
        \RR

    R is a 1-dimensional manifold::

        sage: isinstance(R, Manifold)
        True
        sage: dim(R)
        1

    It is endowed with a canonical chart::

        sage: R.canonical_chart()
        chart (R, (t,))
        sage: R.canonical_chart() is R.default_chart()
        True
        sage: R.atlas()
        [chart (R, (t,))]

    The instance is unique (as long as the constructor arguments are the
    same)::

        sage: R is RealLine()
        True
        sage: R is RealLine(latex_name='R')
        False

    The canonical coordinate is returned by the method
    :meth:`canonical_coordinate`::

        sage: R.canonical_coordinate()
        t
        sage: t = R.canonical_coordinate()
        sage: type(t)
        <type 'sage.symbolic.expression.Expression'>

    However, it can be obtained in the same step as the real line construction
    by means of the shortcut operator ``<>``::

        sage: R.<t> = RealLine()
        sage: t
        t
        sage: type(t)
        <type 'sage.symbolic.expression.Expression'>

    The trick is performed by Sage preparser::

        sage: preparse("R.<t> = RealLine()")
        "R = RealLine(names=('t',)); (t,) = R._first_ngens(1)"

    In particular the ``<>`` operator is to be used to set a canonical
    coordinate symbol different from 't'::

        sage: R.<x> = RealLine()
        sage: R.canonical_chart()
        chart (R, (x,))
        sage: R.atlas()
        [chart (R, (x,))]
        sage: R.canonical_coordinate()
        x

    The LaTeX symbol of the canonical coordinate can be adjusted via the same
    syntax as a chart declaration (see
    :class:`~sage.geometry.manifolds.chart.Chart`)::

        sage: R.<x> = RealLine(coordinate=r'x:\xi')
        sage: latex(x)
        {\xi}
        sage: latex(R.canonical_chart())
        \left(\RR,({\xi})\right)

    The LaTeX symbol of the real line itself can also be customized::

        sage: R.<x> = RealLine(latex_name=r'\mathbb{R}')
        sage: latex(R)
        \mathbb{R}

    Elements of the real line can be constructed directly from a number::

        sage: p = R(2) ; p
        point on field R of real numbers
        sage: p.coord()
        (2,)
        sage: p = R(1.742) ; p
        point on field R of real numbers
        sage: p.coord()
        (1.74200000000000,)

    Symbolic variables can also be used::

        sage: p = R(pi, name='pi') ; p
        point 'pi' on field R of real numbers
        sage: p.coord()
        (pi,)
        sage: a = var('a')
        sage: p = R(a) ; p
        point on field R of real numbers
        sage: p.coord()
        (a,)

    """
    def __init__(self, coordinate=None, name='R', latex_name=r'\RR',
                 start_index=0, names=None):
        r"""
        Construct the real line manifold.

        TESTS::

            sage: R = RealLine() ; R
            field R of real numbers
            sage: R.category()
            Category of sets
            sage: TestSuite(R).run()

        """
        from chart import Chart
        Manifold.__init__(self, 1, name, latex_name=latex_name,
                          start_index=start_index)
        if coordinate is None:
            if names is None:
                coordinate = 't'
            else:
                coordinate = names[0]
        self._canon_chart = Chart(self, coordinates=coordinate)
        self._lower = minus_infinity  # for compatibility with OpenInterval
        self._upper = infinity        # idem

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: R = RealLine()
            sage: R._repr_()
            'field R of real numbers'
            sage: R = RealLine(name='r')
            sage: R._repr_()
            'field r of real numbers'

        """
        return "field " + self._name + " of real numbers"

    def _first_ngens(self, n):
        r"""
        Return the coordinate of the canonical chart

        This is useful only for the use of Sage preparser.

        TESTS::

            sage: R = RealLine()
            sage: R._first_ngens(1)
            (t,)
            sage: R = RealLine(coordinate='x')
            sage: R._first_ngens(1)
            (x,)
            sage: R = RealLine(names=('x',))
            sage: R._first_ngens(1)
            (x,)

        """
        return self._canon_chart[:]

    def _element_constructor_(self, coords=None, chart=None, name=None,
                              latex_name=None, check_coords=True):
        r"""
        Construct an element of ``self``.

        This is a redefinition of
        :meth:`sage.geometry.manifolds.domain.ManifoldSubset._element_constructor_`
        to allow for construction from a number (considered as the canonical
        coordinate)

        EXAMPLES::

            sage: R = RealLine()
            sage: R._element_constructor_((pi,)) # standard use of ManifoldSubset._element_constructor_
            point on field R of real numbers
            sage: R._element_constructor_(pi) # specific use with a single coordinate as argument
            point on field R of real numbers
            sage: R._element_constructor_(pi).coord()
            (pi,)
            sage: R._element_constructor_(pi) == R._element_constructor_((pi,))
            True

        """
        if coords in SR:
            coords = (coords,)
        return super(RealLine, self)._element_constructor_(coords=coords,
                                 chart=chart, name=name, latex_name=latex_name,
                                 check_coords=check_coords)

    def _Hom_(self, other, category=None):
        r"""
        Construct the set of curves in ``other`` with parameter in ``self``.

        INPUT:

        - ``other`` -- an open subset of some manifold
        - ``category`` -- (default: ``None``) not used here (to ensure
          compatibility with generic hook ``_Hom_``)

        OUTPUT:

        - the set of curves R  --> V,  where R is ``self`` and V is ``other``

        See class
        :class:`~sage.geometry.manifolds.manifold_homset.ManifoldCurveSet`
        for more documentation.

        """
        from sage.geometry.manifolds.manifold_homset import ManifoldCurveSet
        return ManifoldCurveSet(self, other)


    def canonical_chart(self):
        r"""
        Return the canonical chart defined on ``self``.

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.chart.Chart`

        EXAMPLES::

            sage: R = RealLine()
            sage: R.canonical_chart()
            chart (R, (t,))
            sage: R.<x> = RealLine()
            sage: R.canonical_chart()
            chart (R, (x,))

        """
        return self._canon_chart

    def canonical_coordinate(self):
        r"""
        Return the canonical coordinate defined on ``self``.

        OUTPUT:

        - the symbolic variable representing the canonical coordinate
          defined on ``self``.

        EXAMPLES::

            sage: R = RealLine()
            sage: R.canonical_coordinate()
            t
            sage: type(R.canonical_coordinate())
            <type 'sage.symbolic.expression.Expression'>
            sage: R.canonical_coordinate().is_real()
            True
            sage: R.<x> = RealLine()
            sage: R.canonical_coordinate()
            x

        """
        return self._canon_chart._xx[0]

    def lower_bound(self):
        r"""
        Return the lower bound (infimum) of ``self`` (considered as an
        interval)

        OUTPUT:

        - always ``-Infinity``.

        EXAMPLE::

            sage: R = RealLine()
            sage: R.lower_bound()
            -Infinity

        An alias of :meth:`lower_bound` is :meth:`inf`::

            sage: R.inf()
            -Infinity

        """
        return self._lower

    inf = lower_bound

    def upper_bound(self):
        r"""
        Return the upper bound (supremum) of ``self`` (considered as an
        interval)

        OUTPUT:

        - always ``+Infinity``.

        EXAMPLE::

            sage: R = RealLine()
            sage: R.upper_bound()
            +Infinity

        An alias of :meth:`upper_bound` is :meth:`sup`::

            sage: R.sup()
            +Infinity

        """
        return self._upper

    sup = upper_bound

    def open_interval(self, lower, upper):
        r"""
        Define an open interval of the field of real numbers represented
        by ``self``.

        INPUT:

        - ``lower`` -- lower bound of the interval (possibly ``-Infinity``)
        - ``upper`` -- upper bound of the interval (possibly ``+Infinity``)

        OUTPUT:

        - instance of class
          :class:`~sage.geometry.manifolds.manifold.OpenInterval`
          representing the open interval (``lower``, ``upper``).

        EXAMPLES:

        The interval `(0,\pi)`::

            sage: R.<t> = RealLine()
            sage: I = R.open_interval(0, pi) ; I
            Real interval (0, pi)

        The interval `(-\infty,1)`::

            sage: I = R.open_interval(-oo, 1) ; I
            Real interval (-Infinity, 1)

        The interval `(0, +\infty)`::

            sage: I = R.open_interval(0, +oo) ; I
            Real interval (0, +Infinity)

        The interval `(-\infty,+\infty)` is `\RR`::

            sage: I = R.open_interval(-oo, +oo) ; I
            field R of real numbers
            sage: I is R
            True

        See :class:`~sage.geometry.manifolds.manifold.OpenInterval` for more
        examples and documentation.

        """
        if lower == minus_infinity and upper == infinity:
            return self
        return OpenInterval(self, lower, upper)

    def diff_mapping(self, codomain, coord_expression, chart=None, name=None,
                     latex_name=None):
        r"""
        Define a differentiable mapping between ``self`` and an open subset of
        some manifold (possibly ``self``).

        This is a redefinition of
        :meth:`sage.geometry.manifolds.domain.ManifoldOpenSubset.diff_mapping`
        in order to have a curve as output.

        See :class:`~sage.geometry.manifolds.curve.ManifoldCurve` for a
        complete documentation.

        INPUT:

        - ``codomain`` -- mapping's codomain (the target manifold or some
          subset of it)
        - ``coord_expression`` -- the coordinate symbolic
          expression of the mapping: list (or tuple) of the coordinates of the
          image expressed
          in terms of the coordinates of the considered point; if the
          dimension of the target manifold is 1, a single expression is
          expected (not a list with a single element)
        - ``chart`` -- (default: None) chart of ``codomain`` in which the
          coordinates are given on the codomain; if none is provided, the
          coordinates are assumed to refer to default chart of ``codomain``
        - ``name`` -- (default: None) name given to the differentiable mapping
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the
          differentiable mapping; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - the differentiable mapping, as an instance of
          :class:`~sage.geometry.manifolds.curve.ManifoldCurve`

        """
        curve_set = self._Hom_(codomain)
        if not isinstance(coord_expression, dict):
            # Turn coord_expression into a dictionary:
            if chart is None:
                chart = codomain.default_chart()
            elif chart not in codomain.atlas():
                raise ValueError("the {} has not been".format(chart) +
                                     " defined on the {}".format(codomain))
            coord_expression = {chart: coord_expression}
        return curve_set(coord_expression, name=name, latex_name=latex_name)

#******************************************************************************

class OpenInterval(ManifoldOpenSubset):
    r"""
    Open real interval.

    This class implements open intervals as open subsets of the real line
    manifold (see :class:`RealLine`).

    INPUT:

    - ``real_line`` -- instance of
      :class:`~sage.geometry.manifolds.manifold.RealLine` representing the
      1-dimensional manifold of real numbers in which the open interval is
      included.
    - ``lower`` -- lower bound of the interval (possibly ``-Infinity``)
    - ``upper`` -- upper bound of the interval (possibly ``+Infinity``)

    EXAMPLES:

    The interval `(0,\pi)`::

        sage: R.<t> = RealLine() ; R
        field R of real numbers
        sage: from sage.geometry.manifolds.manifold import OpenInterval
        sage: I = OpenInterval(R, 0, pi) ; I
        Real interval (0, pi)
        sage: type(I)
        <class 'sage.geometry.manifolds.manifold.OpenInterval_with_category'>

    Instead of importing the class ``OpenInterval`` in the global namespace, it
    is recommended to use the method
    :meth:`~sage.geometry.manifolds.manifold.RealLine.open_interval` of the
    real line instead::

        sage: I = R.open_interval(0, pi) ; I
        Real interval (0, pi)
        sage: latex(I)
        \left( 0 , \pi \right)

    The result is cached (unique representation property)::

        sage: I is OpenInterval(R, 0, pi)
        True
        sage: I is R.open_interval(0, pi)
        True

    ``I`` is an open subset of the field of real numbers, considered as
    a 1-dimensional manifold::

        sage: I.is_subset(R)
        True
        sage: isinstance(I, sage.geometry.manifolds.domain.ManifoldOpenSubset)
        True
        sage: I.manifold()
        field R of real numbers

    An element of ``I``::

        sage: x = I.an_element() ; x
        point on field R of real numbers
        sage: x.coord() # coordinates in the default chart = canonical chart
        (1/2*pi,)

    As for any manifold subset, a specific element of ``I`` can be created
    by providing a tuple containing its coordinate(s) in a given chart::

        sage: x = I((2,)) # (2,) = tuple of coordinates in the canonical chart
        sage: x
        point on field R of real numbers

    But for convenience, it can also be created directly from the coordinate::

        sage: x = I(2) ; x
        point on field R of real numbers
        sage: x.coord()
        (2,)
        sage: I(2) == I((2,))
        True

    ``I`` is endoved with a canonical chart, which is inherited from the
    canonical chart of ``R``::

        sage: I.canonical_chart()
        chart ((0, pi), (t,))

    By default, the coordinates passed for the element ``x`` are those relative
    to the canonical chart::

        sage: I(2) ==  I((2,), chart=I.canonical_chart())
        True

    Since open intervals are facade sets (see
    :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`), the
    parent of their elements is the whole real-line manifold on which
    they are defined, i.e. ``R``::

        sage: I.category()
        Category of facade sets
        sage: x.parent()
        field R of real numbers
        sage: x.parent() is R
        True

    We may check whether a real number (considered as an element of ``R``)
    belongs to ``I``::

        sage: R(2) in I
        True
        sage: R(4) in I
        False

    The lower and upper bounds of the interval ``I``::

        sage: I.lower_bound()
        0
        sage: I.upper_bound()
        pi

    One of the endpoint can be infinite::

        sage: J = R.open_interval(1, +oo) ; J
        Real interval (1, +Infinity)
        sage: R(10^8) in J
        True
        sage: J.an_element().coord()
        (2,)

    """
    def __init__(self, real_line, lower, upper):
        r"""
        Construct an open inverval.

        TESTS::

            sage: from sage.geometry.manifolds.manifold import OpenInterval
            sage: R = RealLine()
            sage: I = OpenInterval(R, -1, 1) ; I
            Real interval (-1, 1)
            sage: TestSuite(I).run()
            sage: J = OpenInterval(R, -oo, 2) ; J
            Real interval (-Infinity, 2)
            sage: TestSuite(J).run()

        """
        from sage.misc.latex import latex
        if not isinstance(real_line, RealLine):
            raise TypeError("{} is not an instance of RealLine".format(
                                                                    real_line))
        name = "({}, {})".format(lower, upper)
        latex_name = r"\left(" + latex(lower) + ", " + latex(upper) + r"\right)"
        ManifoldOpenSubset.__init__(self, real_line, name,
                                    latex_name=latex_name)
        t = real_line.canonical_coordinate()
        if lower != minus_infinity:
            if upper != infinity:
                restrictions = [t>lower, t<upper]
            else:
                restrictions = t>lower
        else:
            if upper != infinity:
                restrictions = t<upper
            else:
                restrictions = None
        self._lower = lower
        self._upper = upper
        self._canon_chart = real_line.canonical_chart().restrict(self,
                                                     restrictions=restrictions)
        self._canon_chart._bounds = (((lower, False), (upper, False)),)

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: R = RealLine()
            sage: I = R.open_interval(-1,pi)
            sage: I._repr_()
            'Real interval (-1, pi)'
            sage: I = R.open_interval(-1,+oo)
            sage: I._repr_()
            'Real interval (-1, +Infinity)'
            sage: I = R.open_interval(-oo,0)
            sage: I._repr_()
            'Real interval (-Infinity, 0)'

        """
        return "Real interval " + self._name

    def _element_constructor_(self, coords=None, chart=None, name=None,
                              latex_name=None, check_coords=True):
        r"""
        Construct an element of ``self``.

        This is a redefinition of
        :meth:`sage.geometry.manifolds.domain.ManifoldSubset._element_constructor_`
        to allow for construction from a number (considered as the canonical
        coordinate)

        EXAMPLES::

            sage: R = RealLine()
            sage: I = R.open_interval(-1, 4)
            sage: I._element_constructor_((2,)) # standard used of ManifoldSubset._element_constructor_
            point on field R of real numbers
            sage: I._element_constructor_(2)  # specific use with a single coordinate
            point on field R of real numbers
            sage: I._element_constructor_(2).coord()
            (2,)
            sage: I._element_constructor_(2) == I._element_constructor_((2,))
            True
            sage: I._element_constructor_(pi)
            point on field R of real numbers
            sage: I._element_constructor_(pi).coord()
            (pi,)
            sage: I._element_constructor_(8)
            Traceback (most recent call last):
            ...
            ValueError: The coordinates (8,) are not valid on the chart
             ((-1, 4), (t,))

        """
        if coords in SR:
            coords = (coords,)
        return super(OpenInterval, self)._element_constructor_(coords=coords,
                                 chart=chart, name=name, latex_name=latex_name,
                                 check_coords=check_coords)


    def _Hom_(self, other, category=None):
        r"""
        Construct the set of curves in ``other`` with parameter in ``self``.

        INPUT:

        - ``other`` -- an open subset of some manifold
        - ``category`` -- (default: ``None``) not used here (to ensure
          compatibility with generic hook ``_Hom_``)

        OUTPUT:

        - the set of curves I  --> V,  where I is ``self`` and V is ``other``

        See class
        :class:`~sage.geometry.manifolds.manifold_homset.ManifoldCurveSet`
        for more documentation.

        """
        from sage.geometry.manifolds.manifold_homset import ManifoldCurveSet
        return ManifoldCurveSet(self, other)



    def canonical_chart(self):
        r"""
        Return the canonical chart defined on ``self``.

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.chart.Chart`

        EXAMPLES:

        Canonical chart on the interval `(0,\pi)`::

            sage: R.<t> = RealLine()
            sage: I = R.open_interval(0, pi)
            sage: I.canonical_chart()
            chart ((0, pi), (t,))

        The canonical chart is the restriction of the canonical chart of
        `\RR` (represented by ``R``) to ``I``::

            sage: I.canonical_chart() == R.canonical_chart().restrict(I)
            True

        The information about the coordinate bounds has been passed to the
        chart (``False`` means that the bound is not part of the interval)::

            sage: I.canonical_chart().coord_bounds()
            (((0, False), (pi, False)),)

        The symbol used for the coordinate of the canonical chart is that
        defined during the construction of the real line::

            sage: R.<x> = RealLine()
            sage: I = R.open_interval(0, pi)
            sage: I.canonical_chart()
            chart ((0, pi), (x,))

        """
        return self._canon_chart

    def canonical_coordinate(self):
        r"""
        Return the canonical coordinate defined on ``self``.

        OUTPUT:

        - the symbolic variable representing the canonical coordinate
          defined on ``self``.

        EXAMPLES:

        Canonical coordinate on the interval `(0,\pi)`::

            sage: R.<t> = RealLine()
            sage: I = R.open_interval(0, pi)
            sage: I.canonical_coordinate()
            t
            sage: type(I.canonical_coordinate())
            <type 'sage.symbolic.expression.Expression'>
            sage: t.is_real()
            True

        The canonical coordinate is the first (unique) coordinate of the
        canonical chart::

            sage: I.canonical_coordinate() is I.canonical_chart()[0]
            True

        Its symbol is inherited from that declared during the construction
        of the real line::

            sage: R.<x> = RealLine()
            sage: I = R.open_interval(0, pi)
            sage: I.canonical_coordinate()
            x

        """
        return self._canon_chart._xx[0]

    def lower_bound(self):
        r"""
        Return the lower bound (infimum) of ``self``.

        EXAMPLES::

            sage: R.<t> = RealLine()
            sage: I = R.open_interval(1/4, 3)
            sage: I.lower_bound()
            1/4
            sage: J = R.open_interval(-oo, 2)
            sage: J.lower_bound()
            -Infinity

        An alias of :meth:`lower_bound` is :meth:`inf`::

            sage: I.inf()
            1/4
            sage: J.inf()
            -Infinity

        """
        return self._lower

    inf = lower_bound

    def upper_bound(self):
        r"""
        Return the upper bound (supremum) of ``self``.

        EXAMPLES::

            sage: R.<t> = RealLine()
            sage: I = R.open_interval(1/4, 3)
            sage: I.upper_bound()
            3
            sage: J = R.open_interval(1, +oo)
            sage: J.upper_bound()
            +Infinity

        An alias of :meth:`upper_bound` is :meth:`sup`::

            sage: I.sup()
            3
            sage: J.sup()
            +Infinity

        """
        return self._upper

    sup = upper_bound
