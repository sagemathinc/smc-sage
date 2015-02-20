r"""
Curves in manifolds

Given a differentiable manifold `M`, a *differentiable curve* curve in
`M` is a differential mapping

.. MATH::

    \gamma: I \longrightarrow N

where `I` is an interval of `\RR`. 

Differentiable curves are implemented by the class :class:`ManifoldCurve`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

- Chap. 1 of S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*,
  vol. 1, Interscience Publishers (New York) (1963)
- Chap. 3 of J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed.,
  Springer (New York) (2013)

"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.symbolic.expression import Expression
from sage.misc.latex import latex
from sage.geometry.manifolds.point import ManifoldPoint
from sage.geometry.manifolds.diffmapping import DiffMapping
from sage.geometry.manifolds.utilities import simplify_chain

class ManifoldCurve(DiffMapping):
    r"""
    Curve in a differentiable manifold.


    EXAMPLES:

    LaTeX representation of curves::
    
        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: R.<t> = RealLine()
        sage: c = M.curve([cos(t), sin(2*t)], (t, 0, 2*pi))
        sage: latex(c)
        "\\mbox{Curve in the 2-dimensional manifold 'M'}"
        sage: c = M.curve([cos(t), sin(2*t)], (t, 0, 2*pi), name='c')
        sage: latex(c)
        'c'
        sage: c = M.curve([cos(t), sin(2*t)], (t, 0, 2*pi), name='c',
        ....:             latex_name=r'\gamma')
        sage: latex(c)
        '\\gamma'

    Curves `\RR\longrightarrow\RR` can be composed: the operator `\circle` is
    denoted by ``*``::

        sage: f = R.curve(t^2, (t,-oo,+oo))
        sage: g = R.curve(cos(t), (t,-oo,+oo))
        sage: s = g*f ; s
        Curve in the field R of real numbers
        sage: s.display()
        R --> R
           t |--> cos(t^2)
        sage: s = g*f ; s
        sage: s = f*g ; s
        Curve in the field R of real numbers
        sage: s.display()
        R --> R
           t |--> cos(t)^2

    """
    def __init__(self, parent, coord_expression=None, name=None,
                 latex_name=None, is_diffeomorphism=False, is_identity=False):
        r"""
        Construct a curve.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: I = R.open_interval(0, 2*pi)
            sage: from sage.geometry.manifolds.curve import ManifoldCurve
            sage: c = ManifoldCurve(Hom(I,M), coord_expression={X: (cos(t), sin(2*t))},
            ....:                   name='c') ; c
            Curve 'c' in the 2-dimensional manifold 'M'

        To pass the test suite, the curve must be constructed from the parent
        and not from a direct call to ManifoldCurve::

            sage: c = Hom(I,M)({X: (cos(t), sin(2*t))},  name='c') ; c
            Curve 'c' in the 2-dimensional manifold 'M'
            sage: TestSuite(c).run()

        The identity of interval I::

            sage: c = Hom(I,I)({}, is_identity=True) ; c
            Curve 'Id_(0, 2*pi)' in the Real interval (0, 2*pi)
            sage: TestSuite(c).run()

        """
        domain = parent.domain()
        codomain = parent.codomain()
        if coord_expression is None:
            coord_functions = None
        else:
            if not isinstance(coord_expression, dict):
                raise TypeError("{} is not a dictionary".format(
                                                             coord_expression))
            param_chart = domain.canonical_chart()
            codom_atlas = codomain.atlas()
            n = codomain.manifold().dim()
            coord_functions = {}
            for chart, expr in coord_expression.iteritems():
                if isinstance(chart, tuple):
                    # a pair of charts is passed:
                    coord_functions[chart] = expr
                else:
                    coord_functions[(param_chart, chart)] = expr
        DiffMapping.__init__(self, parent, coord_functions=coord_functions,
                             name=name, latex_name=latex_name,
                             is_diffeomorphism=is_diffeomorphism,
                             is_identity=is_identity)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: c = M.curve([cos(t), sin(2*t)], (t, 0, 2*pi))
            sage: c._repr_()
            "Curve in the 2-dimensional manifold 'M'"
            sage: c = M.curve([cos(t), sin(2*t)], (t, 0, 2*pi), name='c')
            sage: c._repr_()
            "Curve 'c' in the 2-dimensional manifold 'M'"

        """
        if self._codomain._manifold._dim == 1:
            return DiffMapping._repr_(self)
        description = "Curve "
        if self._name is not None:
            description += "'" + self._name + "' "
        description += "in the {}".format(self._codomain)
        return description

    def coord_functions(self, chart=None):
        r"""
        Return the coordinate functions expressing ``self`` in a given chart.

        INPUT:

        - ``chart`` -- (default: None) chart on the curve's codomain; if
          ``None``, the codomain's default chart is assumed

        OUTPUT:

        - symbolic expression representing the curve in the above chart

        EXAMPLES:

        Cartesian and polar expression of a curve in the Euclidean plane::

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
            sage: R.<t> = RealLine()
            sage: c = U.curve({c_spher: (1,t)}, (t, 0, 2*pi), name='c')
            sage: c.coord_functions(c_spher)
            (1, t)
            sage: c.coord_functions(c_cart)
            (cos(t), sin(t))

        Since ``c_cart`` is the default chart on ``U``, it can be omitted::

            sage: c.coord_functions()
            (cos(t), sin(t))

        Cartesian expression of a cardiod::
        
            sage: c = U.curve({c_spher: (2*(1+cos(t)), t)}, (t, 0, 2*pi), name='c')
            sage: c.coord_functions(c_cart)
            (2*cos(t)^2 + 2*cos(t), 2*(cos(t) + 1)*sin(t))

        """
        return self.expr(chart1=self._domain.canonical_chart(), chart2=chart)

    def __call__(self, t, simplify=True):
        r"""
        Image for a given value of the curve parameter.

        This is a redefinition of :meth:`sage.categories.map.Map.__call__`
        to allow for the direct call with some value of the parameter
        (numerical value or symbolic expression) instead
        of the element (ManifoldPoint) of the domain corresponding to that
        value. 

        EXAMPLES:

        Points on circle in the Euclidean plane::
        
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: c = M.curve([cos(t), sin(t)], (t, 0, 2*pi), name='c')
            sage: c(0)
            point 'c(0)' on 2-dimensional manifold 'M'
            sage: c(0) in M
            True
            sage: c(0).coord(X)
            (1, 0)
            sage: c(pi/4).coord(X)
            (1/2*sqrt(2), 1/2*sqrt(2))
            sage: c(t)
            point 'c(t)' on 2-dimensional manifold 'M'
            sage: c(t).coord(X)
            (cos(t), sin(t))
        
        """
        # Case of a point in the domain:
        if isinstance(t, ManifoldPoint):
            return DiffMapping.__call__(t)
        # Case of a value of the canonical coordinate in the domain:
        codom = self._codomain
        canon_chart = self._domain._canon_chart
        canon_coord = canon_chart._xx[0]
        if (canon_chart, codom._def_chart) in self._coord_expression:
            chart_pair = (canon_chart, codom._def_chart)
        else:
            chart_pair = self._coord_expression.keys()[0]  # a chart is picked
                                                           # at random
        coord_functions = self._coord_expression[chart_pair]._functions
        n = codom._manifold._dim
        dict_subs = {canon_coord: t}
        coords = [coord_functions[i].expr().substitute(dict_subs)
                  for i in range(n)]
        if simplify:
            coords = [simplify_chain(coords[i]) for i in range(n)]
        if self._name is not None:
            name = self._name + "({})".format(t)
        else:
            name = None
        if self._latex_name is not None:
            latex_name = self._latex_name + r"\left(" + latex(t) + r"\right)"
        else:
            latex_name = None
        return codom.element_class(codom, coords=coords, chart=chart_pair[1],
                                   name=name, latex_name=latex_name,
                                   check_coords=False)
