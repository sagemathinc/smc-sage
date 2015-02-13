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

from sage.structure.sage_object import SageObject
from sage.rings.infinity import Infinity
from sage.geometry.manifolds.domain import ManifoldOpenSubset
from sage.geometry.manifolds.chart import MultiFunctionChart

class ManifoldCurve(SageObject):
    r"""
    Curve in a differentiable manifold.
    
    """
    def __init__(self, codomain, coord_expression, tmin=-Infinity,
                 tmax=+Infinity, name=None, latex_name=None):
        r"""
        Construct a curve.
        """
        if not isinstance(codomain, ManifoldOpenSubset):
            raise TypeError(
                 "{} is not an open subset of some manifold".format(codomain))
        if not isinstance(coord_expression, dict):
            raise TypeError("{} is not a dictionary".format(coord_expression))
        self._codomain = codomain
        self._tmin = tmin
        self._tmax = tmax
        self._name = name
        if latex_name is None:
            self._latex_name = name
        else:
            self._latex_name = latex_name
        self._coord_expression = {}
        codom_atlas = codomain.atlas()
        n = codomain._manifold._dim
        for chart, expr in coord_expression.iteritems():
            if chart not in codom_atlas:
                raise ValueError("{} is not a chart defined on {}".format(
                                                              chart, codomain))
            if len(expr) != n:
                raise ValueError("{} coordinate functions must be ".format(n) + 
                                                                    "provided")
            self._coord_expression[chart] = MultiFunctionChart(chart, *expr)


    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        description = "Curve "
        if self._name is not None:
            description += "'" + self._name + "' "
        description += "in the {}".format(self._codomain)
        return description

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.
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
        if not isinstance(other, ManifoldCurve):
            return False
        if self._codomain != other._codomain:
            return False
        if self._tmin != other._tmin:
            return False
        if self._tmax != other._tmax:
            return False
        for chart, functions in self._coord_expression.iteritems():
            try:
                if functions.expr() != other.expr(chart):
                    return False
            except ValueError:
                return False
        return True


