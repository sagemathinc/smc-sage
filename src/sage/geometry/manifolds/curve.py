r"""
Curves on manifolds

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
from sage.geometry.manifolds.domain import ManifoldOpenSubset
from sage.geometry.manifolds.diffmapping import DiffMapping
from sage.geometry.manifolds.utilities import simplify_chain

class ManifoldCurve(DiffMapping):
    r"""
    Curve in a differentiable manifold.
    
    """
    def __init__(self, parent, coord_expression=None, name=None,
                 latex_name=None, is_diffeomorphism=False, is_identity=False):
        r"""
        Construct a curve.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: from sage.geometry.manifolds.curve import ManifoldCurve
            sage: c = ManifoldCurve(M, {X: (cos(t), sin(2*t))}, tmin=0,
            ....:                   tmax=2*pi, name='c') ; c
            Curve 'c' in the 2-dimensional manifold 'M'
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
                if chart not in codom_atlas:
                    raise ValueError("{} is not a chart defined on {}".format(
                                                              chart, codomain))
                if len(expr) != n:
                    raise ValueError("{} coordinate functions ".format(n) + 
                                     "must be provided")
                coord_functions[(param_chart, chart)] = expr
        DiffMapping.__init__(self, parent, coord_functions=coord_functions,
                             name=name, latex_name=latex_name,
                             is_diffeomorphism=is_diffeomorphism,
                             is_identity=is_identity)
        # self._param = self._domain.canonical_coordinate()
        # self._tmin = self._domain.lower_bound()
        # self._tmax = self._domain.upper_bound()

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
        description = "Curve "
        if self._name is not None:
            description += "'" + self._name + "' "
        description += "in the {}".format(self._codomain)
        return description

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: R.<t> = RealLine()
            sage: c = M.curve([cos(t), sin(2*t)])
            sage: c._latex_()
            "\\mbox{Curve in the 2-dimensional manifold 'M'}"
            sage: c = M.curve([cos(t), sin(2*t)], name='c')
            sage: c._latex_()
            'c'
            sage: c = M.curve([cos(t), sin(2*t)], name='c', latex_name=r'\gamma')
            sage: c._latex_()
            '\\gamma'

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def expr(self, chart=None):
        r"""
        Return the expression of ``self`` in a given chart.

        INPUT:

        - ``chart`` -- (default: None) chart on the curve's codomain; if
          ``None``, the codomain's default chart is assumed

        OUTPUT:

        - symbolic expression representing the curve in the above chart

        EXAMPLES:

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
            sage: c = U.curve({c_spher: (1,t)}, tmin=0, tmax=2*pi, name='c')
            sage: c.expr(c_spher)
            (1, t)
            sage: c.expr(c_cart)
            (cos(t), sin(t))

        Since ``c_cart`` is the default chart on ``U``, it can be omitted::

            sage: c.expr()
            (cos(t), sin(t))

        Cartesian expression of a cardiod::
        
            sage: c = U.curve({c_spher: (2*(1+cos(t)), t)}, tmin=0, tmax=2*pi, name='c')
            sage: c.expr(c_cart)
            (2*cos(t)^2 + 2*cos(t), 2*(cos(t) + 1)*sin(t))

        """
        codom = self._codomain
        def_chart = codom._def_chart
        if chart is None:
            chart = def_chart
        if chart not in self._coord_expression:
            # Some computation must be performed
            sel_chart = None # selected chart
            if def_chart in self._coord_expression \
                                and (def_chart, chart) in codom._coord_changes:
                sel_chart = def_chart
            else:
                for ochart in self._coord_expression:
                    if (ochart, chart) in codom._coord_changes:
                        sel_chart = ochart
                        break
            if sel_chart is None:
                raise ValueError("No start chart has been found to compute " +
                      "the expression of {} in the {}".format(self, chart))
            oexpr = self._coord_expression[sel_chart]
            chg = codom._coord_changes[(sel_chart, chart)]
            self._coord_expression[chart] = chg(*oexpr)
        return self._coord_expression[chart]

    def __call__(self, t, simplify=True):
        r"""
        Image of a given value of the parameter.

        EXAMPLES:

        Points on circle in the Euclidean plane::
        
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: c = M.curve([cos(t), sin(t)], name='c')
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
        codom = self._codomain
        def_chart = codom._def_chart
        if def_chart in self._coord_expression:
            chart = def_chart
        else:
            chart = self._coord_expression.keys()[0]  # a chart is picked at
                                                      # random
        expr = self._coord_expression[chart]
        n = len(expr)
        dict_subs = {self._param: t}
        coords = [expr[i].substitute(dict_subs) for i in range(n)]
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
        return codom.point(coords=coords, chart=chart, name=name,
                           latex_name=latex_name)





        
