r"""
Scalar fields

The class :class:`ScalarField` implements scalar fields on differentiable
manifolds over `\RR`, i.e. differentiable mappings of the form

.. MATH::

    f: U\subset M \longrightarrow \RR

where `U` is an open subset of the differentiable manifold `M`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version

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

from sage.structure.element import CommutativeAlgebraElement
from sage.rings.integer import Integer
from domain import ManifoldSubset
from chart import FunctionChart, ZeroFunctionChart, MultiFunctionChart

class ScalarField(CommutativeAlgebraElement):
    r"""
    Scalar field on a differentiable manifold.

    This class implements scalar fields on differentiable
    manifolds over `\RR`, i.e. differentiable mappings of the form

    .. MATH::

        f: U\subset M \longrightarrow \RR

    where `U` is an open subset of the differentiable manifold `M`.

    The class :class:`ScalarField`  inherits from the class  :class:`~sage.structure.element.CommutativeAlgebraElement`, a scalar field
    on `U` being an element of the commutative algebra `C^\infty(U)`, see
    :class:`~sage.geometry.manifolds.scalarfield_algebra.ScalarFieldAlgebra`.

    INPUT:

    - ``domain`` -- the manifold open subset `U` on which the scalar field is
      defined (must be an instance of class
      :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`)
    - ``coord_expression`` -- (default: ``None``) coordinate expression(s) of
      the scalar field; this can be either

      - a dictionary of coordinate expressions in various charts on the domain,
        with the charts as keys;
      - a single expression, which is the same in all charts defined on the
        domain (constant scalar field).

      NB: If ``coord_expression`` is ``None`` or incomplete, coordinate
      expressions can be added after the creation of the object, by means of
      the methods :meth:`add_expr`, :meth:`add_expr_by_continuation` and
      :meth:`set_expr`
    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      scalar field
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      scalar field; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A scalar field on the 2-sphere::

        sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                                intersection_name='W',
        ....:                                restrictions1= x^2+y^2!=0,
        ....:                                restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2), c_uv: (u^2+v^2)/(1+u^2+v^2)},
        ....:                    name='f') ; f
        scalar field 'f' on the 2-dimensional manifold 'M'
        sage: f.display()
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    For scalar fields defined by a single coordinate expression, the latter
    can be passed instead of the dictionary over the charts::

        sage: g = U.scalar_field(x*y, chart=c_xy, name='g') ; g
        scalar field 'g' on the open subset 'U' of the 2-dimensional manifold 'M'

    The above is indeed equivalent to::

        sage: g = U.scalar_field({c_xy: x*y}, name='g') ; g
        scalar field 'g' on the open subset 'U' of the 2-dimensional manifold 'M'

    Since ``c_xy`` is the default chart of ``U``, the argument ``chart`` can
    be skipped::

        sage: g = U.scalar_field(x*y, name='g') ; g
        scalar field 'g' on the open subset 'U' of the 2-dimensional manifold 'M'

    The scalar field `g` is defined on `U` and has an expression in terms of
    the coordinates `(u,v)` on `W=U\cap V`::

        sage: g.display()
        g: U --> R
           (x, y) |--> x*y
        on W: (u, v) |--> u*v/(u^4 + 2*u^2*v^2 + v^4)

    Scalar fields on `M` can also be declared with a single chart::

        sage: f = M.scalar_field(1/(1+x^2+y^2), chart=c_xy, name='f') ; f
        scalar field 'f' on the 2-dimensional manifold 'M'

    Their definition must then be completed by providing the expressions on
    other charts, via the method :meth:`add_expr`, to get a global cover of
    the manifold::

        sage: f.add_expr((u^2+v^2)/(1+u^2+v^2), chart=c_uv)
        sage: f.display()
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    We can even first declare the scalar field without any coordinate
    expression and provide them subsequently::

        sage: f = M.scalar_field(name='f')
        sage: f.add_expr(1/(1+x^2+y^2), chart=c_xy)
        sage: f.add_expr((u^2+v^2)/(1+u^2+v^2), chart=c_uv)
        sage: f.display()
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    We may also use the method :meth:`add_expr_by_continuation` to complete
    the coordinate definition using the analytic continuation from domains in
    which charts overlap::

        sage: f = M.scalar_field(1/(1+x^2+y^2), chart=c_xy, name='f') ; f
        scalar field 'f' on the 2-dimensional manifold 'M'
        sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
        sage: f.display()
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    A scalar field can also be defined by some unspecified function of the
    coordinates::

        sage: h = U.scalar_field(function('H', x, y), name='h') ; h
        scalar field 'h' on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: h.display()
        h: U --> R
           (x, y) |--> H(x, y)
        on W: (u, v) |--> H(u/(u^2 + v^2), v/(u^2 + v^2))

    We may use the argument ``latex_name`` to specify the LaTeX symbol denoting
    the scalar field if the latter is different from ``name``::

        sage: latex(f)
        f
        sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2), c_uv: (u^2+v^2)/(1+u^2+v^2)},
        ....:                    name='f', latex_name=r'\mathcal{F}')
        sage: latex(f)
        \mathcal{F}

    The coordinate expression in a given chart is obtained via the method
    :meth:`expr`, which returns a symbolic expression::

        sage: f.expr(c_uv)
        (u^2 + v^2)/(u^2 + v^2 + 1)
        sage: type(f.expr(c_uv))
        <type 'sage.symbolic.expression.Expression'>

    The method :meth:`function_chart` returns instead a function of the
    chart coordinates, i.e. an instance of
    :class:`~sage.geometry.manifolds.chart.FunctionChart`::

        sage: f.function_chart(c_uv)
        (u^2 + v^2)/(u^2 + v^2 + 1)
        sage: type(f.function_chart(c_uv))
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        sage: f.function_chart(c_uv).display()
        (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    The value returned by the method :meth:`expr` is actually the coordinate
    expression of the chart function::

        sage: f.expr(c_uv) is f.function_chart(c_uv).expr()
        True

    A constant scalar field is declared by setting the argument ``chart`` to
    ``'all'``::

        sage: c = M.scalar_field(2, chart='all', name='c') ; c
        scalar field 'c' on the 2-dimensional manifold 'M'
        sage: c.display()
        c: M --> R
        on U: (x, y) |--> 2
        on V: (u, v) |--> 2

    The constant value can be some unspecified parameter::

        sage: var('a')
        a
        sage: c = M.scalar_field(a, chart='all', name='c') ; c
        scalar field 'c' on the 2-dimensional manifold 'M'
        sage: c.display()
        c: M --> R
        on U: (x, y) |--> a
        on V: (u, v) |--> a

    A special case of constant field is the zero scalar field::

        sage: zer = M.scalar_field(0, chart='all') ; zer
        scalar field 'zero' on the 2-dimensional manifold 'M'
        sage: zer.display()
        zero: M --> R
        on U: (x, y) |--> 0
        on V: (u, v) |--> 0

    It can be obtained directly by means of the function
    :meth:`~sage.geometry.manifolds.domain.ManifoldOpenSubset.zero_scalar_field`::

        sage: zer is M.zero_scalar_field()
        True

    A third way is to get the zero scalar field as the zero element of the
    algebra `C^\infty(M)` of scalar fields on `M` (see below)::

        sage: zer is M.scalar_field_algebra().zero()
        True

    By definition, a scalar field acts on the manifold's points, sending
    them to real numbers::

        sage: N = M.point((0,0), chart=c_uv) # the North pole
        sage: S = M.point((0,0), chart=c_xy) # the South pole
        sage: E = M.point((1,0), chart=c_xy) # a point at the equator
        sage: f(N)
        0
        sage: f(S)
        1
        sage: f(E)
        1/2
        sage: h(E)
        H(1, 0)
        sage: c(E)
        a
        sage: zer(E)
        0

    A scalar field can be compared to another scalar field::

        sage: f == g
        False

    ...to a symbolic expression::

        sage: f == x*y
        False
        sage: g == x*y
        True
        sage: c == a
        True

    ...to a number::

        sage: f == 2
        False
        sage: zer == 0
        True

    ...to anything else::

        sage: f == M
        False

    .. RUBRIC:: Arithmetics of scalar fields

    Scalar fields on `M` (resp. `U`) belong to the algebra `C^\infty(M)`
    (resp. `C^\infty(U)`)::

        sage: f.parent()
        algebra of scalar fields on the 2-dimensional manifold 'M'
        sage: f.parent() is M.scalar_field_algebra()
        True
        sage: g.parent()
        algebra of scalar fields on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: g.parent() is U.scalar_field_algebra()
        True

    Consequently, scalar fields can be added::

        sage: s = f + c ; s
        scalar field 'f+c' on the 2-dimensional manifold 'M'
        sage: s.display()
        f+c: M --> R
        on U: (x, y) |--> (a*x^2 + a*y^2 + a + 1)/(x^2 + y^2 + 1)
        on V: (u, v) |--> ((a + 1)*u^2 + (a + 1)*v^2 + a)/(u^2 + v^2 + 1)

    and subtracted::

        sage: s = f - c ; s
        scalar field 'f-c' on the 2-dimensional manifold 'M'
        sage: s.display()
        f-c: M --> R
        on U: (x, y) |--> -(a*x^2 + a*y^2 + a - 1)/(x^2 + y^2 + 1)
        on V: (u, v) |--> -((a - 1)*u^2 + (a - 1)*v^2 + a)/(u^2 + v^2 + 1)

    Some tests::

        sage: f + zer == f
        True
        sage: f - f == zer
        True
        sage: f + (-f) == zer
        True
        sage: (f+c)-f == c
        True
        sage: (f-c)+c == f
        True

    We may add a number (interpretted as a constant scalar field) to a scalar
    field::

        sage: s = f + 1 ; s
        scalar field on the 2-dimensional manifold 'M'
        sage: s.display()
        M --> R
        on U: (x, y) |--> (x^2 + y^2 + 2)/(x^2 + y^2 + 1)
        on V: (u, v) |--> (2*u^2 + 2*v^2 + 1)/(u^2 + v^2 + 1)
        sage: (f+1)-1 == f
        True
        sage: s = a + f ; s
        scalar field on the 2-dimensional manifold 'M'
        sage: s == c + f
        True

    The addition of two scalar fields with different domains is possible if
    the domain of one of them is a subset of the domain of the other; the
    domain of the result is then this subset::

        sage: f.domain()
        2-dimensional manifold 'M'
        sage: g.domain()
        open subset 'U' of the 2-dimensional manifold 'M'
        sage: s = f + g ; s
        scalar field on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.domain()
        open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.display()
        U --> R
        (x, y) |--> (x*y^3 + (x^3 + x)*y + 1)/(x^2 + y^2 + 1)
        on W: (u, v) |--> (u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6 + u*v^3 + (u^3 + u)*v)/(u^6 + v^6 + (3*u^2 + 1)*v^4 + u^4 + (3*u^4 + 2*u^2)*v^2)

    The operation actually performed is `f|_U + g`::

        sage: s == f.restrict(U) + g
        True

    In Sage framework, the addition of `f` and `g` is permitted because
    there is a *coercion* of the parent of `f`, namely `C^\infty(M)`, to
    the parent of `g`, namely `C^\infty(U)`::

        sage: CM = M.scalar_field_algebra()
        sage: CU = U.scalar_field_algebra()
        sage: CU.has_coerce_map_from(CM)
        True

    The coercion map is nothing but the restriction to domain `U`::

        sage: CU.coerce(f) == f.restrict(U)
        True

    Since the algebra `C^\infty(M)` is a vector space over `\RR`, scalar fields
    can be multiplied by a number, either an explicit one::

        sage: s = 2*f ; s
        scalar field on the 2-dimensional manifold 'M'
        sage: s.display()
        M --> R
        on U: (x, y) |--> 2/(x^2 + y^2 + 1)
        on V: (u, v) |--> 2*(u^2 + v^2)/(u^2 + v^2 + 1)

    or a symbolic one::

        sage: s = a*f ; s
        scalar field on the 2-dimensional manifold 'M'
        sage: s.display()
        M --> R
        on U: (x, y) |--> a/(x^2 + y^2 + 1)
        on V: (u, v) |--> (a*u^2 + a*v^2)/(u^2 + v^2 + 1)

    Some tests::

        sage: 0*f == 0
        True
        sage: 0*f == zer
        True
        sage: 1*f == f
        True
        sage: (-2)*f == - f - f
        True

    The ring multiplication of the algebras `C^\infty(M)` and `C^\infty(U)`
    is the pointwise multiplication of functions::

        sage: s = f*f ; s
        scalar field 'f*f' on the 2-dimensional manifold 'M'
        sage: s.display()
        f*f: M --> R
        on U: (x, y) |--> 1/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1)
        on V: (u, v) |--> (u^4 + 2*u^2*v^2 + v^4)/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1)
        sage: s = g*h ; s
        scalar field 'g*h' on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.display()
        g*h: U --> R
           (x, y) |--> x*y*H(x, y)
        on W: (u, v) |--> u*v*H(u/(u^2 + v^2), v/(u^2 + v^2))/(u^4 + 2*u^2*v^2 + v^4)

    Thanks to the coercion `C^\infty(M)\rightarrow C^\infty(U)` mentionned
    above, it is possible to multiply a scalar field defined on `M` by a
    scalar field defined on `U`, the result being a scalar field defined on
    `U`::

        sage: f.domain(), g.domain()
        (2-dimensional manifold 'M', open subset 'U' of the 2-dimensional manifold 'M')
        sage: s = f*g ; s
        scalar field on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.display()
        U --> R
        (x, y) |--> x*y/(x^2 + y^2 + 1)
        on W: (u, v) |--> u*v/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2)
        sage: s == f.restrict(U)*g
        True

    Scalar fields can be divided (pointwise division)::

        sage: s = f/c ; s
        scalar field 'f/c' on the 2-dimensional manifold 'M'
        sage: s.display()
        f/c: M --> R
        on U: (x, y) |--> 1/(a*x^2 + a*y^2 + a)
        on V: (u, v) |--> (u^2 + v^2)/(a*u^2 + a*v^2 + a)
        sage: s = g/h ; s
        scalar field 'g/h' on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.display()
        g/h: U --> R
           (x, y) |--> x*y/H(x, y)
        on W: (u, v) |--> u*v/((u^4 + 2*u^2*v^2 + v^4)*H(u/(u^2 + v^2), v/(u^2 + v^2)))
        sage: s = f/g ; s
        scalar field on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.display()
        U --> R
        (x, y) |--> 1/(x*y^3 + (x^3 + x)*y)
        on W: (u, v) |--> (u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6)/(u*v^3 + (u^3 + u)*v)
        sage: s == f.restrict(U)/g
        True

    For scalar fields defined on a single chart domain, we may perform some
    arithmetics with symbolic expressions involving the chart coordinates::

        sage: s = g + x^2 - y ; s
        scalar field on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.display()
        U --> R
        (x, y) |--> x^2 + (x - 1)*y
        on W: (u, v) |--> -(v^3 - u^2 + (u^2 - u)*v)/(u^4 + 2*u^2*v^2 + v^4)

    ::

        sage: s = g*x ; s
        scalar field on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.display()
        U --> R
        (x, y) |--> x^2*y
        on W: (u, v) |--> u*v*x/(u^4 + 2*u^2*v^2 + v^4)

    ::

        sage: s = g/x ; s
        scalar field on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.display()
        U --> R
        (x, y) |--> y
        on W: (u, v) |--> u*v/((u^4 + 2*u^2*v^2 + v^4)*x)
        sage: s = x/g ; s
        scalar field on the open subset 'U' of the 2-dimensional manifold 'M'
        sage: s.display()
        U --> R
        (x, y) |--> 1/y
        on W: (u, v) |--> (u^2 + v^2)/v

    """
    def __init__(self, domain, coord_expression=None, name=None,
                 latex_name=None):
        CommutativeAlgebraElement.__init__(self, domain.scalar_field_algebra())
        self._manifold = domain._manifold
        self._domain = domain
        self._tensor_type = (0,0)
        self._is_zero = False # a priori, may be changed below or via
                              # method __nonzero__()
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._express = {} # dict of coordinate expressions (FunctionChart
                           # instances) with charts as keys
        if coord_expression is not None:
            if isinstance(coord_expression, FunctionChart):
                self._express[coord_expression._chart] = coord_expression
            elif isinstance(coord_expression, dict):
                for chart, expression in coord_expression.iteritems():
                    if isinstance(expression, FunctionChart):
                        self._express[chart] = expression
                    else:
                        self._express[chart] = FunctionChart(chart, expression)
            elif coord_expression == 0:
                self._is_zero = True
                for chart in self._domain._atlas:
                    self._express[chart] = chart._zero_function
            else:
                # coord_expression is independent of the chart (constant scalar
                # field)
                for chart in self._domain._atlas:
                    self._express[chart] = FunctionChart(chart,
                                                         coord_expression)
        self._init_derived()   # initialization of derived quantities

    ####### Required methods for an algebra element (beside arithmetic) #######

    def __nonzero__(self):
        r"""
        Return True if ``self`` is nonzero and False otherwise.

        This method is called by self.is_zero().

        EXAMPLES:

        Tests on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y)
            sage: f.is_zero()
            False
            sage: f.set_expr(0)
            sage: f.is_zero()
            True
            sage: g = M.scalar_field(0)
            sage: g.is_zero()
            True
            sage: M.zero_scalar_field().is_zero()
            True

        """
        if self._is_zero:
            return False
        if not self._express:
            # undefined scalar field
            return True
        iszero = True
        for funct in self._express.itervalues():
            iszero = iszero and funct.is_zero()
        self._is_zero = iszero
        return not iszero

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a scalar field

        OUTPUT:

        - True if ``self`` is equal to ``other``,  or False otherwise

        """
        if not isinstance(other, ScalarField):
            # We try a conversion of other to a scalar field, except if
            # other is None (since this would generate an undefined scalar
            # field)
            if other is None:
                return False
            try:
                other = self.parent()(other)  # conversion to a scalar field
            except TypeError:
                return False
        if other._domain != self._domain:
            return False
        if other.is_zero():
            return self.is_zero()
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the comparison.")
        resu = True
        for chart in com_charts:
            resu = resu and (self._express[chart] == other._express[chart])
        return resu

    def __ne__(self, other):
        r"""
        Non-equality operator.
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

    ####### End of required methods for an algebra element (beside arithmetic) #######

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        self._restrictions = {} # dict. of restrictions of self on subsets
                                # of self._domain, with the subsets as keys
        self._differential = None # differential
        self._lie_derivatives = {} # collection of Lie derivatives of self

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        self._restrictions.clear()
        self._differential = None
        # First deletes any reference to self in the vectors' dictionary:
        for vid, val in self._lie_derivatives.iteritems():
            del val[0]._lie_der_along_self[id(self)]
        self._lie_derivatives.clear()

    def _repr_(self):
        r"""
        String representation of ``self``.
        """
        description = "scalar field"
        if self._name is not None:
            description += " '%s'" % self._name
        description += " on the " + str(self._domain)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def set_name(self, name=None, latex_name=None):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``name`` -- (string; default: None) name given to the scalar field
        - ``latex_name`` -- (string; default: None) LaTeX symbol to denote
          the scalar field; if None while ``name`` is provided, the LaTeX
          symbol is set to ``name``.

        """
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` and on the same domain.

        """
        return self.__class__(self._domain)

    def domain(self):
        r"""
        Return the open subset on which the scalar field is defined.

        OUTPUT:

        - instance of class :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`
          representing the manifold's open subset on which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x+2*y)
            sage: f.domain()
            2-dimensional manifold 'M'
            sage: U = M.open_subset('U', coord_def={c_xy: x<0})
            sage: g = f.restrict(U)
            sage: g.domain()
            open subset 'U' of the 2-dimensional manifold 'M'

        """
        return self._domain

    def tensor_type(self):
        r"""
        Return the tensor type of ``self``, i.e. (0,0), when ``self`` is
        considered as a tensor field on the manifold.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x+2*y)
            sage: f.tensor_type()
            (0, 0)

        """
        return self._tensor_type

    def copy(self):
        r"""
        Return an exact copy of ``self``.

        EXAMPLES:

        Copy on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: g = f.copy()
            sage: type(g)
            <class 'sage.geometry.manifolds.scalarfield.ScalarFieldAlgebra_with_category.element_class'>
            sage: g.expr()
            x*y^2
            sage: g == f
            True
            sage: g is f
            False

        """
        result = self.__class__(self._domain)  #!# what about the name ?
        for chart, funct in self._express.iteritems():
            result._express[chart] = funct.copy()
        return result

    def function_chart(self, chart=None, from_chart=None):
        r"""
        Return the function of the coordinates representing the scalar field
        in a given chart.

        INPUT:

        - ``chart`` -- (default: None) chart with respect to which the
          coordinate expression is to be returned; if None, the
          default chart of the domain of ``self`` will be used
        - ``from_chart`` -- (default: None) chart from which the
          required expression is computed if it is not known already in the
          chart ``chart``; if None, a chart is picked in ``self._express``

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.chart.FunctionChart`
          representing the coordinate function of the scalar field in the
          given chart.

        EXAMPLES:

        Coordinate function on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: f.function_chart()
            x*y^2
            sage: f.function_chart(c_xy)  # equivalent form (since c_xy is the default chart)
            x*y^2
            sage: type(f.function_chart())
            <class 'sage.geometry.manifolds.chart.FunctionChart'>

        Expression via a change of coordinates::

            sage: c_uv.<u,v> = M.chart()
            sage: c_uv.coord_change(c_xy, u+v, u-v)
            coordinate change from chart (M, (u, v)) to chart (M, (x, y))
            sage: f._express # at this stage, f is expressed only in terms of (x,y) coordinates
            {chart (M, (x, y)): x*y^2}
            sage: f.function_chart(c_uv) # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: f.function_chart(c_uv) == (u+v)*(u-v)^2  # check
            True
            sage: f._express  # random (dict. output); f has now 2 coordinate expressions:
            {chart (M, (x, y)): x*y^2, chart (M, (u, v)): u^3 - u^2*v - u*v^2 + v^3}

        Usage in a physical context (simple Lorentz transformation - boost in
        x direction, with relative velocity v between o1 and o2 frames)::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: o1.<t,x> = M.chart()
            sage: o2.<T,X> = M.chart()
            sage: f = M.scalar_field(x^2 - t^2)
            sage: f.function_chart(o1)
            -t^2 + x^2
            sage: v = var('v'); gam = 1/sqrt(1-v^2)
            sage: o2.coord_change(o1, gam*(T - v*X), gam*(X - v*T))
            coordinate change from chart (M, (T, X)) to chart (M, (t, x))
            sage: f.function_chart(o2)
            -T^2 + X^2

        """
        if chart is None:
            chart = self._domain._def_chart
        else:
            if chart not in self._domain._atlas:
                raise TypeError("The " + str(chart) + " has not " + \
                      " been defined on the domain " + str(self._domain))
        if chart not in self._express:
            # Check whether chart corresponds to a subchart of a chart
            # where the expression of self is known:
            for known_chart in self._express:
                if chart in known_chart._subcharts:
                    new_expr = self._express[known_chart].expr()
                    self._express[chart] = FunctionChart(chart, new_expr)
                    return self._express[chart]
            # If this point is reached, the expression must be computed
            # from that in the chart from_chart, by means of a
            # change-of-coordinates formula:
            if from_chart is None:
                # from_chart in searched among the charts of known expressions
                # and subcharts of them
                known_express = self._express.copy()
                found = False
                for kchart in known_express:
                    for skchart in kchart._subcharts:
                        if (chart, skchart) in self._domain._coord_changes:
                            from_chart = skchart
                            found = True
                            if skchart not in self._express:
                                self._express[skchart] = FunctionChart(skchart,
                                                  self._express[kchart].expr())
                            break
                    if found:
                        break
                if not found:
                    raise ValueError("No starting chart could be found to " +
                           "compute the expression in the " + str(chart))
            change = self._domain._coord_changes[(chart, from_chart)]
            # old coordinates expressed in terms of the new ones:
            coords = [ change._transf._functions[i]._express
                       for i in range(self._manifold._dim) ]
            new_expr = self._express[from_chart](*coords)
            self._express[chart] = FunctionChart(chart, new_expr)
            self._del_derived()
        return self._express[chart]


    def expr(self, chart=None, from_chart=None):
        r"""
        Return the coordinate expression of the scalar field in a given
        chart.

        INPUT:

        - ``chart`` -- (default: None) chart with respect to which the
          coordinate expression is required; if None, the default
          chart of the domain of ``self`` will be used
        - ``from_chart`` -- (default: None) chart from which the
          required expression is computed if it is not known already in the
          chart ``chart``; if None, a chart is picked in ``self._express``

        OUTPUT:

        - symbolic expression representing the coordinate
          expression of the scalar field in the given chart.

        EXAMPLES:

        Expression of a scalar field on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: f.expr()
            x*y^2
            sage: f.expr(c_xy)  # equivalent form (since c_xy is the default chart)
            x*y^2
            sage: type(f.expr())
            <type 'sage.symbolic.expression.Expression'>

        Expression via a change of coordinates::

            sage: c_uv.<u,v> = M.chart()
            sage: c_uv.coord_change(c_xy, u+v, u-v)
            coordinate change from chart (M, (u, v)) to chart (M, (x, y))
            sage: f._express # at this stage, f is expressed only in terms of (x,y) coordinates
            {chart (M, (x, y)): x*y^2}
            sage: f.expr(c_uv) # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: bool( f.expr(c_uv) == (u+v)*(u-v)^2 ) # check
            True
            sage: f._express  # random (dict. output); f has now 2 coordinate expressions:
            {chart (M, (x, y)): x*y^2, chart (M, (u, v)): u^3 - u^2*v - u*v^2 + v^3}

        """
        return self.function_chart(chart, from_chart)._express

    def set_expr(self, coord_expression, chart=None):
        r"""
        Set the coordinate expression of the scalar field.

        The expressions with respect to other charts are deleted, in order to
        avoid any inconsistency. To keep them, use :meth:`add_expr` instead.

        INPUT:

        - ``coord_expression`` -- coordinate expression of the scalar field
        - ``chart`` -- (default: None) chart in which ``coord_expression`` is
          defined; if None, the default chart of the domain of ``self`` is
          assumed

        EXAMPLES:

        Setting scalar field expressions on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2 + 2*x*y +1)
            sage: f._express
            {chart (M, (x, y)): x^2 + 2*x*y + 1}
            sage: f.set_expr(3*y)
            sage: f._express  # the (x,y) expression has been changed:
            {chart (M, (x, y)): 3*y}
            sage: c_uv.<u,v> = M.chart()
            sage: f.set_expr(cos(u)-sin(v), c_uv)
            sage: f._express # the (x,y) expression has been lost:
            {chart (M, (u, v)): cos(u) - sin(v)}
            sage: f.set_expr(3*y)
            sage: f._express # the (u,v) expression has been lost:
            {chart (M, (x, y)): 3*y}

        """
        if chart is None:
            chart = self._domain._def_chart
        self._is_zero = False # a priori
        self._express.clear()
        self._express[chart] = FunctionChart(chart, coord_expression)
        self._del_derived()

    def add_expr(self, coord_expression, chart=None):
        r"""
        Add some coordinate expression to the scalar field.

        The previous expressions with respect to other charts are kept. To
        clear them, use :meth:`set_expr` instead.

        INPUT:

        - ``coord_expression`` -- coordinate expression of the scalar field
        - ``chart`` -- (default: None) chart in which ``coord_expression``
          is defined; if None, the default chart of the domain of ``self`` is
          assumed

        .. WARNING::

            If the scalar field has already expressions in other charts, it
            is the user's responsability to make sure that the expression
            to be added is consistent with them.

        EXAMPLES:

        Adding scalar field expressions on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2 + 2*x*y +1)
            sage: f._express
            {chart (M, (x, y)): x^2 + 2*x*y + 1}
            sage: f.add_expr(3*y)
            sage: f._express  # the (x,y) expression has been changed:
            {chart (M, (x, y)): 3*y}
            sage: c_uv.<u,v> = M.chart()
            sage: f.add_expr(cos(u)-sin(v), c_uv)
            sage: f._express # random (dict. output); f has now 2 expressions:
            {chart (M, (x, y)): 3*y, chart (M, (u, v)): cos(u) - sin(v)}

        """
        if chart is None:
            chart = self._domain._def_chart
        self._express[chart] = FunctionChart(chart, coord_expression)
        self._is_zero = False # a priori
        self._del_derived()

    def add_expr_by_continuation(self, chart, subdomain):
        r"""
        Set coordinate expression in a chart by continuation of the
        coordinate expression in a subchart.

        The continuation is performed by demanding that the coordinate
        expression is identical to that in the restriction of the chart to
        a given subdomain.

        INPUT:

        - ``chart`` -- coordinate chart `(U,(x^i))` in which the expression of
          the scalar field is to set
        - ``subdomain`` -- open subset `V\subset U` in which the expression
          in terms of the restriction of the coordinate chart `(U,(x^i))` to
          `V` is already known or can be evaluated by a change of coordinates.

        EXAMPLE:

        Scalar field on the sphere `S^2`::

            sage: M = Manifold(2, 'S^2')
            sage: U = M.open_subset('U') ; V = M.open_subset('V') # the complement of resp. N pole and S pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coordinates
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                             intersection_name='W', restrictions1= x^2+y^2!=0, \
                                             restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W =  U.intersection(V)  # S^2 minus the two poles
            sage: f = M.scalar_field(atan(x^2+y^2), chart=c_xy, name='f')

        The scalar field has been defined only on the domain covered by the
        chart c_xy, i.e. `U`::

            sage: f.display()
            f: S^2 --> R
            on U: (x, y) |--> arctan(x^2 + y^2)

        We note that on `W = U\cap V`, the expression of `f` in terms of
        coordinates `(u,v)` can be deduced from that in the coordinates
        `(x,y)` thanks to the transition map between the two charts::

            sage: f.display(c_uv.restrict(W))
            f: S^2 --> R
            on W: (u, v) |--> arctan(1/(u^2 + v^2))

        We use this fact to extend the definition of `f` to the open subset `V`,
        covered by the chart c_uv::

            sage: f.add_expr_by_continuation(c_uv, W)

        Then, `f` is known on the whole sphere::

            sage: f.display()
            f: S^2 --> R
            on U: (x, y) |--> arctan(x^2 + y^2)
            on V: (u, v) |--> arctan(1/(u^2 + v^2))

        """
        if not chart._domain.is_subset(self._domain):
            raise ValueError("The chart is not defined on a subset of " +
                             "the scalar field domain.")
        schart = chart.restrict(subdomain)
        self._express[chart] = FunctionChart(chart, self.expr(schart))
        self._is_zero = False # a priori
        self._del_derived()

    def _display_expression(self, chart, result):
        r"""
        Helper function for :meth:`display`.
        """
        from sage.misc.latex import latex
        try:
            expression = self.function_chart(chart)
            coords = chart[:]
            if len(coords) == 1:
                coords = coords[0]
            if chart._domain == self._domain:
                if self._name is not None:
                    result._txt += "   "
                result._latex += " & "
            else:
                result._txt += "on " + chart._domain._name + ": "
                result._latex += r"\mbox{on}\ " + latex(chart._domain) + r": & "
            result._txt += repr(coords) + " |--> " + repr(expression) + "\n"
            result._latex += latex(coords) + r"& \longmapsto & " + \
                            latex(expression) + r"\\"
        except (TypeError, ValueError):
            pass

    def display(self, chart=None):
        r"""
        Display the expression of the scalar field in a given chart.

        Without any argument, this function displays the expressions of the
        scalar field in all the charts defined on the scalar field's domain
        that are not restrictions of another chart to some subdomain
        (the "top charts").

        INPUT:

        - ``chart`` -- (default: None) chart with respect to which the
          coordinate expression is to be displayed; if None, the display is
          performed in all the top charts in which the coordinate expression is
          known.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        EXAMPLES:

        Various displays::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(sqrt(x+1), name='f')
            sage: f.display()
            f: M --> R
               (x, y) |--> sqrt(x + 1)
            sage: latex(f.display())
            \begin{array}{llcl} f:& M & \longrightarrow & \RR \\ & \left(x, y\right) & \longmapsto & \sqrt{x + 1} \end{array}
            sage: g = M.scalar_field(function('G', x, y), name='g')
            sage: g.display()
            g: M --> R
               (x, y) |--> G(x, y)
            sage: latex(g.display())
            \begin{array}{llcl} g:& M & \longrightarrow & \RR \\ & \left(x, y\right) & \longmapsto & G\left(x, y\right) \end{array}

        A shortcut of ``display()`` is ``disp()``::

            sage: f.disp()
            f: M --> R
               (x, y) |--> sqrt(x + 1)

        """
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion
        result = FormattedExpansion()
        if self._name is None:
            symbol = ""
        else:
            symbol = self._name + ": "
        result._txt = symbol + self._domain._name + " --> R\n"
        if self._latex_name is None:
            symbol = ""
        else:
            symbol = self._latex_name + ":"
        result._latex = r"\begin{array}{llcl} " + symbol + r"&" + \
                       latex(self._domain) + r"& \longrightarrow & \RR \\"
        if chart is None:
            for ch in self._domain._top_charts:
                self._display_expression(ch, result)
        else:
            self._display_expression(chart, result)
        result._txt = result._txt[:-1]
        result._latex = result._latex[:-2] + r"\end{array}"
        return result

    disp = display

    def view(self, chart=None):
        r"""
        Deprecated method.

        Use method :meth:`display` instead.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field((x+1)^2, name='f')
            sage: f.view()
            doctest:...: DeprecationWarning: Use function display() instead.
            See http://trac.sagemath.org/15916 for details.
            f: M --> R
               (x, y) |--> (x + 1)^2
            sage: f.display()
            f: M --> R
               (x, y) |--> (x + 1)^2

        """
        from sage.misc.superseded import deprecation
        deprecation(15916, 'Use function display() instead.')
        return self.display(chart=chart)


    def restrict(self, subdomain):
        r"""
        Restriction of the scalar field to a subdomain of its domain of
        definition.

        INPUT:

        - ``subdomain`` -- the subdomain (instance of
          :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`)

        OUTPUT:

        - instance of :class:`ScalarField` representing the restriction of
          ``self`` to ``subdomain``.

        EXAMPLE:

        Restriction of a scalar field defined on `\RR^2` to the unit open
        disc::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # Cartesian coordinates
            sage: U = M.open_subset('U', coord_def={X: x^2+y^2 < 1}) # U unit open disc
            sage: f = M.scalar_field(cos(x*y), name='f')
            sage: f_U = f.restrict(U) ; f_U
            scalar field 'f' on the open subset 'U' of the 2-dimensional manifold 'M'
            sage: f_U.display()
            f: U --> R
               (x, y) |--> cos(x*y)
            sage: f.parent()
            algebra of scalar fields on the 2-dimensional manifold 'M'
            sage: f_U.parent()
            algebra of scalar fields on the open subset 'U' of the 2-dimensional manifold 'M'

        The restriction to the whole domain is the identity::

            sage: f.restrict(M) is f
            True
            sage: f_U.restrict(U) is f_U
            True

        Restriction of the zero scalar field::

            sage: M.zero_scalar_field().restrict(U)
            scalar field 'zero' on the open subset 'U' of the 2-dimensional manifold 'M'
            sage: M.zero_scalar_field().restrict(U) is U.zero_scalar_field()
            True

        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("The specified domain is not a subset " +
                                 "of the domain of definition of the scalar " +
                                 "field.")
            # Special case of the zero scalar field:
            if self._is_zero:
                return subdomain._zero_scalar_field
            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.iteritems():
                if subdomain.is_subset(dom):
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    break
            else:
            # If this fails, the restriction is obtained via coercion
                resu = subdomain.scalar_field_algebra()(self)
                resu._name = self._name
                resu._latex_name = self._latex_name
                self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def common_charts(self, other):
        r"""
        Find common charts for the expressions of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a scalar field

        OUPUT:

        - list of common charts; if no common chart is found, None is
          returned (instead of an empty list).

        EXAMPLES:

        Search for common charts on a 2-dimensional manifold with 2
        overlapping domains::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: c_xy.<x,y> = U.chart()
            sage: V = M.open_subset('V')
            sage: c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: f = U.scalar_field(x^2)
            sage: g = M.scalar_field(x+y)
            sage: f.common_charts(g)
            [chart (U, (x, y))]
            sage: g.add_expr(u, c_uv)
            sage: f._express
            {chart (U, (x, y)): x^2}
            sage: g._express  # random (dictionary output)
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u}
            sage: f.common_charts(g)
            [chart (U, (x, y))]

        Common charts found as subcharts: the subcharts are introduced via
        a transition map between charts c_xy and c_uv on the intersecting
        subdomain `W = U\cap V`::

            sage: trans = c_xy.transition_map(c_uv, (x+y, x-y), 'W', x<0, u+v<0)
            sage: M.atlas()
            [chart (U, (x, y)), chart (V, (u, v)), chart (W, (x, y)), chart (W, (u, v))]
            sage: c_xy_W = M.atlas()[2]
            sage: c_uv_W = M.atlas()[3]
            sage: trans.inverse()
            coordinate change from chart (W, (u, v)) to chart (W, (x, y))
            sage: f.common_charts(g)
            [chart (U, (x, y))]
            sage: f.expr(c_xy_W)
            x^2
            sage: f._express  # random (dictionary output)
            {chart (U, (x, y)): x^2, chart (W, (x, y)): x^2}
            sage: g._express  # random (dictionary output)
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u}
            sage: g.common_charts(f)  # c_xy_W is not returned because it is subchart of 'xy'
            [chart (U, (x, y))]
            sage: f.expr(c_uv_W)
            1/4*u^2 + 1/2*u*v + 1/4*v^2
            sage: f._express  # random (dictionary output)
            {chart (U, (x, y)): x^2, chart (W, (x, y)): x^2, chart (W, (u, v)): 1/4*u^2 + 1/2*u*v + 1/4*v^2}
            sage: g._express  # random (dictionary output)
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u}
            sage: f.common_charts(g)
            [chart (U, (x, y)), chart (W, (u, v))]
            sage: # the expressions have been updated on the subcharts
            sage: g._express #  random (dictionary output)
            {chart (U, (x, y)): x + y, chart (V, (u, v)): u, chart (W, (u, v)): u}

        Common charts found by computing some coordinate changes::

            sage: W = U.intersection(V)
            sage: f = W.scalar_field(x^2, c_xy_W)
            sage: g = W.scalar_field(u+1, c_uv_W)
            sage: f._express
            {chart (W, (x, y)): x^2}
            sage: g._express
            {chart (W, (u, v)): u + 1}
            sage: f.common_charts(g)
            [chart (W, (u, v)), chart (W, (x, y))]
            sage: f._express # random (dictionary output)
            {chart (W, (u, v)): 1/4*u^2 + 1/2*u*v + 1/4*v^2, chart (W, (x, y)): x^2}
            sage: g._express # random (dictionary output)
            {chart (W, (u, v)): u + 1, chart (W, (x, y)): x + y + 1}

        """
        if not isinstance(other, ScalarField):
            raise TypeError("The second argument must be a scalar field.")
        dom1 = self._domain
        dom2 = other._domain
        coord_changes = self._manifold._coord_changes
        resu = []
        #
        # 1/ Search for common charts among the existing expressions, i.e.
        #    without performing any expression transformation.
        #    -------------------------------------------------------------
        for chart1 in self._express:
            if chart1 in other._express:
                resu.append(chart1)
        # Search for a subchart:
        known_expr1 = self._express.copy()
        known_expr2 = other._express.copy()
        for chart1 in known_expr1:
            if chart1 not in resu:
                for chart2 in known_expr2:
                    if chart2 not in resu:
                        if chart2 in chart1._subcharts:
                            self.expr(chart2)
                            resu.append(chart2)
                        if chart1 in chart2._subcharts:
                            other.expr(chart1)
                            resu.append(chart1)
        #
        # 2/ Search for common charts via one expression transformation
        #    ----------------------------------------------------------
        for chart1 in known_expr1:
            if chart1 not in resu:
                for chart2 in known_expr2:
                    if chart2 not in resu:
                        if (chart1, chart2) in coord_changes:
                            self.function_chart(chart2, from_chart=chart1)
                            resu.append(chart2)
                        if (chart2, chart1) in coord_changes:
                            other.function_chart(chart1, from_chart=chart2)
                            resu.append(chart1)
        if resu == []:
            return None
        else:
            return resu

    def __call__(self, p, chart=None):
        r"""
        Compute the value of the scalar field at a given point.

        INPUT:

        - ``p`` -- point in the scalar field's domain (type:
          :class:`~sage.geometry.manifolds.point.ManifoldPoint`)
        - ``chart`` -- (default: None) chart in which the coordinates of p
          are to be considered; if none is provided, a chart in which both p's
          coordinates and the expression of ``self`` are known is searched,
          starting from the default chart of self._domain

        OUTPUT:

        - value at p

        EXAMPLES:

        """
        #!# it should be "if p not in self_domain:" instead, but this test is
        # skipped for efficiency
        if p not in self._manifold:
            raise ValueError("The point " + str(p) +
                             " does not belong to the " + str(self._manifold))
        if self._is_zero:
            return 0
        if chart is None:
            # A common chart is searched:
            def_chart = self._domain._def_chart
            if def_chart in p._coordinates and def_chart in self._express:
                chart = def_chart
            else:
                for chart_p in p._coordinates:
                    if chart_p in self._express:
                        chart = chart_p
                        break
        if chart is None:
            # A change of coordinates is attempted for p:
            for chart_s in self._express:
                try:
                    p.coord(chart_s)
                    chart = chart_s
                    break
                except ValueError:
                    pass
            else:
                # A change of coordinates is attempted on the scalar field
                # expressions:
                for chart_p in p._coordinates:
                    try:
                        self.function_chart(chart_p)
                        chart = chart_p
                        break
                    except (TypeError, ValueError):
                        pass
        if chart is None:
            raise ValueError("No common chart has been found to evaluate " \
                "the action of " + str(self) + " on the " + str(p) + ".")
        return self._express[chart](*(p._coordinates[chart]))

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of ``self``

        """
        result = self._new_instance()
        for chart in self._express:
            result._express[chart] = + self._express[chart]
        if self._name is not None:
            result._name = '+' + self._name
        if self._latex_name is not None:
            result._latex_name = '+' + self._latex_name
        return result

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the negative of ``self``

        """
        result = self._new_instance()
        for chart in self._express:
            result._express[chart] = - self._express[chart]
        if self._name is not None:
            result._name = '-' + self._name
        if self._latex_name is not None:
            result._latex_name = '-' + self._latex_name
        return result


    #########  CommutativeAlgebraElement arithmetic operators ########

    def _add_(self, other):
        r"""
        Scalar field addition.

        INPUT:

        - ``other`` -- a scalar field (in the same algebra as self)

        OUPUT:

        - the scalar field resulting from the addition of ``self`` and
          ``other``

        """
        dom = self._domain
        # Special cases:
        if self._is_zero:
            return other.copy()
        if other._is_zero:
            return self.copy()
        # Generic case:
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the addition.")
        result = self.__class__(dom)
        for chart in com_charts:
            # FunctionChart addition:
            result._express[chart] = self._express[chart] + other._express[chart]
        if result.is_zero():
            return dom._zero_scalar_field
        if self._name is not None and other._name is not None:
            result._name = self._name + '+' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            result._latex_name = self._latex_name + '+' + other._latex_name
        return result

    def _sub_(self, other):
        r"""
        Scalar field subtraction.

        INPUT:

        - ``other`` -- a scalar field (in the same algebra as self)

        OUPUT:

        - the scalar field resulting from the subtraction of ``other`` from
          ``self``

        """
        dom = self._domain
        # Special cases:
        if self._is_zero:
            return -other
        if other._is_zero:
            return self.copy()
        # Generic case:
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the subtraction.")
        result = self.__class__(dom)
        for chart in com_charts:
            # FunctionChart subtraction:
            result._express[chart] = self._express[chart] - other._express[chart]
        if result.is_zero():
            return dom._zero_scalar_field
        if self._name is not None and other._name is not None:
            result._name = self._name + '-' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            result._latex_name = self._latex_name + '-' + other._latex_name
        return result


    def _mul_(self, other):
        r"""
        Scalar field multiplication.

        INPUT:

        - ``other`` -- a scalar field (in the same algebra as self)

        OUPUT:

        - the scalar field resulting from the multiplication of ``self`` by
          ``other``

        """
        from sage.tensor.modules.format_utilities import format_mul_txt, \
                                                         format_mul_latex
        dom = self._domain
        # Special cases:
        if self._is_zero or other._is_zero:
            return dom._zero_scalar_field
        # Generic case:
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the multiplication.")
        result = self.__class__(dom)
        for chart in com_charts:
            # FunctionChart multiplication:
            result._express[chart] = self._express[chart] * other._express[chart]
        result._name = format_mul_txt(self._name, '*', other._name)
        result._latex_name = format_mul_latex(self._latex_name, ' ',
                                             other._latex_name)
        return result

    def _div_(self, other):
        r"""
        Scalar field division.

        INPUT:

        - ``other`` -- a scalar field (in the same algebra as self)

        OUPUT:

        - the scalar field resulting from the division of ``self`` by
          ``other``

        """
        from sage.tensor.modules.format_utilities import format_mul_txt, \
                                                         format_mul_latex
        dom = self._domain
        # Special cases:
        if other._is_zero:
            raise ZeroDivisionError("Division of a scalar field by zero.")
        if self._is_zero:
            return dom._zero_scalar_field
        # Generic case:
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("No common chart for the division.")
        result = self.__class__(dom)
        for chart in com_charts:
            # FunctionChart division:
            result._express[chart] = self._express[chart] / other._express[chart]
        result._name = format_mul_txt(self._name, '/', other._name)
        result._latex_name = format_mul_latex(self._latex_name, '/',
                                             other._latex_name)
        return result

    def _lmul_(self, number):
        r"""
        Multiplication on the left of a scalar field by a real number.

        INPUT:

        - ``number`` -- an element of the ring on which the algebra is defined;
          mathematically, this should be a real number; here it is a member of
          the symbolic ring SR.

        OUPUT:

        - the scalar field ``number*self``

        """
        if number == 0:
            return self._domain._zero_scalar_field
        result = self.__class__(self._domain)
        for chart, expr in self._express.iteritems():
            result._express[chart] = number * expr
        return result

    def _rmul_(self, number):
        r"""
        Multiplication on the right of a scalar field by a real number.

        INPUT:

        - ``number`` -- an element of the ring on which the algebra is defined;
          mathematically, this should be a real number; here it is a member of
          the symbolic ring SR.

        OUPUT:

        - the scalar field ``number*self``

        """
        return self._lmul_(number) # since the algebra is commutative


    #########  End of CommutativeAlgebraElement arithmetic operators ########

    def differential(self):
        r"""
        Return the differential of the scalar field.

        OUTPUT:

        - the 1-form that is the differential of ``self``.

        EXAMPLES:

        Differential 1-form on a 3-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: f = M.scalar_field(cos(x)*z^3 + exp(y)*z^2, name='f')
            sage: df = f.differential() ; df
            1-form 'df' on the 3-dimensional manifold 'M'
            sage: df.display()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f

        Since the exterior derivative of a scalar field (considered a 0-form)
        is nothing but its differential, ``exterior_der()`` is an alias of
        ``differential()``::

            sage: df = f.exterior_der() ; df
            1-form 'df' on the 3-dimensional manifold 'M'
            sage: df.display()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f

        Differential computed on a chart that is not the default one::

            sage: c_uvw.<u,v,w> = M.chart()
            sage: g = M.scalar_field(u*v^2*w^3, c_uvw, name='g')
            sage: dg = g.differential() ; dg
            1-form 'dg' on the 3-dimensional manifold 'M'
            sage: dg._components
            {coordinate frame (M, (d/du,d/dv,d/dw)): 1-index components w.r.t.
             coordinate frame (M, (d/du,d/dv,d/dw))}
            sage: dg.comp(c_uvw.frame())[:, c_uvw]
            [v^2*w^3, 2*u*v*w^3, 3*u*v^2*w^2]
            sage: dg.display(c_uvw.frame(), c_uvw)
            dg = v^2*w^3 du + 2*u*v*w^3 dv + 3*u*v^2*w^2 dw

        The exterior derivative is nilpotent::

            sage: ddf = df.exterior_der() ; ddf
            2-form 'ddf' on the 3-dimensional manifold 'M'
            sage: ddf == 0
            True
            sage: ddf[:] # for the incredule
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: ddg = dg.exterior_der() ; ddg
            2-form 'ddg' on the 3-dimensional manifold 'M'
            sage: ddg == 0
            True

        """
        from sage.tensor.modules.format_utilities import format_unop_txt, \
                                                         format_unop_latex
        if self._differential is None:
            # A new computation is necessary:
            rname = format_unop_txt('d', self._name)
            rlname = format_unop_latex(r'\mathrm{d}', self._latex_name)
            self._differential = self._domain.one_form(name=rname,
                                                             latex_name=rlname)
            if self._is_zero:
                for chart in self._domain._atlas:
                    self._differential.add_comp(chart._frame) # since a newly
                                            # created set of components is zero
            else:
                for chart, func in self._express.iteritems():
                    diff_func = self._differential.add_comp(chart._frame)
                    for i in self._manifold.irange():
                        diff_func[i, chart] = func.diff(i)
        return self._differential

    exterior_der = differential

    def lie_der(self, vector):
        r"""
        Computes the Lie derivative with respect to a vector field.

        The Lie derivative is stored in the dictionary
        :attr:`_lie_derivatives`, so that there is no need to
        recompute it at the next call if neither ``self`` nor ``vector``
        have been modified meanwhile.

        In the present case (scalar field), the Lie derivative is equal to
        the scalar field resulting from the action of the vector field on
        ``self``.

        INPUT:

        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken

        OUTPUT:

        - the scalar field that is the Lie derivative of ``self`` with
          respect to ``vector``

        EXAMPLES:

        Lie derivative on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2*cos(y))
            sage: v = M.vector_field(name='v')
            sage: v[:] = (-y, x)
            sage: f.lie_der(v)
            scalar field on the 2-dimensional manifold 'M'
            sage: f.lie_der(v).expr()
            -x^3*sin(y) - 2*x*y*cos(y)

        Alternative expressions of the Lie derivative of a scalar field::

            sage: f.lie_der(v) == v(f)  # the vector acting on f
            True
            sage: f.lie_der(v) == f.differential()(v)  # the differential of f acting on the vector
            True

        A vanishing Lie derivative::

            sage: f.set_expr(x^2 + y^2)
            sage: f.lie_der(v).display()
            M --> R
            (x, y) |--> 0

        """
        if id(vector) not in self._lie_derivatives:
            # A new computation must be performed
            res = vector(self)
            self._lie_derivatives[id(vector)] = (vector, res)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]


    def hodge_star(self, metric):
        r"""
        Compute the Hodge dual of the scalar field with respect to some
        pseudo-Riemannian metric.

        If `f` is ``self``, the Hodge dual is the `n`-form
        `*f` defined by (`n` being the manifold's dimension)

        .. MATH::

            *f = f \epsilon

        where `\epsilon` is the volume form associated with some
        pseudo-Riemannian metric `g` on the manifold.

        INPUT:

        - ``metric``: the pseudo-Riemannian metric `g` defining the Hodge dual,
          via the volume form `\epsilon`; must be an instance of
          :class:`~sage.geometry.manifolds.metric.Metric`

        OUTPUT:

        - the `n`-form `*f`

        EXAMPLES:

        Hodge star of a scalar field in the Euclidean space `R^3`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: f = M.scalar_field(function('F',x,y,z), name='f')
            sage: sf = f.hodge_star(g) ; sf
            3-form '*f' on the 3-dimensional manifold 'M'
            sage: sf.display()
            *f = F(x, y, z) dx/\dy/\dz
            sage: ssf = sf.hodge_star(g) ; ssf
            scalar field '**f' on the 3-dimensional manifold 'M'
            sage: ssf.display()
            **f: M --> R
               (x, y, z) |--> F(x, y, z)
            sage: ssf == f # must hold for a Riemannian metric
            True

        """
        from sage.tensor.modules.format_utilities import format_unop_txt, \
                                                         format_unop_latex
        eps = metric.volume_form()
        dom_resu = self._domain.intersection(eps._domain) # result domain
        resu = self.restrict(dom_resu) * eps.restrict(dom_resu)
        if self._name is None:
            resu_name = None
        else:
            resu_name = format_unop_txt('*', self._name)
        if self._latex_name is None:
            resu_latex_name = None
        else:
            resu_latex_name = format_unop_latex(r'\star ', self._latex_name)
        resu.set_name(name=resu_name, latex_name=resu_latex_name)
        return resu

    def _function_name(self, func, func_latex):
        r"""
        Helper function to set the symbol of a function applied to the
        scalar field.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f', latex_name=r"\Phi")
            sage: f._function_name("cos", r"\cos")
            ('cos(f)', '\\cos\\left(\\Phi\\right)')
            sage: f = M.scalar_field({X: x+y})  # no name given to f
            sage: f._function_name("cos", r"\cos")
            (None, None)

        """
        if self._name is None:
            name = None
        else:
            name = func + "(" + self._name + ")"
        if self._latex_name is None:
            latex_name = None
        else:
            latex_name = func_latex + r"\left(" + self._latex_name + r"\right)"
        return name, latex_name

    def exp(self):
        r"""
        Exponential of the scalar field.

        OUTPUT:

        - the scalar field `\exp f`, where `f` is the current scalar field.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f', latex_name=r"\Phi")
            sage: g = exp(f) ; g
            scalar field 'exp(f)' on the 2-dimensional manifold 'M'
            sage: g.display()
            exp(f): M --> R
               (x, y) |--> e^(x + y)
            sage: latex(g)
            \exp\left(\Phi\right)

        Automatic simplifications occur::

            sage: f = M.scalar_field({X: 2*ln(1+x^2)}, name='f')
            sage: exp(f).display()
            exp(f): M --> R
               (x, y) |--> x^4 + 2*x^2 + 1

        The inverse function is :meth:`log`::

            sage: log(exp(f)) == f
            True

        """
        name, latex_name = self._function_name("exp", r"\exp")
        resu = self.__class__(self._domain, name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.exp()
        return resu

    def log(self):
        r"""
        Natural logarithm of the scalar field

        OUTPUT:

        - the scalar field `\ln f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f', latex_name=r"\Phi")
            sage: g = log(f) ; g
            scalar field 'ln(f)' on the 2-dimensional manifold 'M'
            sage: g.display()
            ln(f): M --> R
               (x, y) |--> log(x + y)
            sage: latex(g)
            \ln\left(\Phi\right)

        The inverse function is :meth:`exp`::

            sage: exp(log(f)) == f
            True

        """
        name, latex_name = self._function_name("ln", r"\ln")
        resu = self.__class__(self._domain, name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.log()
        return resu

    def __pow__(self, exponent):
        r"""
        The scalar field to a given power.

        INPUT:

        - ``exponent`` -- the exponent

        OUTPUT:

        - the scalar field `f^a`, where `f` is the current scalar field and
          `a` the exponent

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f', latex_name=r'\Phi')
            sage: g = f.__pow__(pi) ; g
            scalar field 'f^pi' on the 2-dimensional manifold 'M'
            sage: latex(g)
            {\Phi}^{ \pi }
            sage: g.display()
            f^pi: M --> R
               (x, y) |--> (x + y)^pi

        The global function ``pow`` can be used::

            sage: pow(f, pi) == f.__pow__(pi)
            True

        as well as the exponent notation::

            sage: f^pi == f.__pow__(pi)
            True

        Some checks::

            sage: pow(f, 2) == f*f
            True
            sage: pow(pow(f, 1/2), 2) == f
            True

        """
        from sage.misc.latex import latex
        if self._name is None:
            name = None
        else:
            name = self._name + "^{}".format(exponent)
        if self._latex_name is None:
            latex_name = None
        else:
            latex_name = r"{" + self._latex_name + r"}^{" + \
                         latex(exponent) + r"}"
        resu = self.__class__(self._domain, name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.__pow__(exponent)
        return resu
