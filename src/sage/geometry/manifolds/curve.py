r"""
Curves in manifolds

Given a differentiable manifold `M`, a *differentiable curve* curve in
`M` is a differentiable mapping

.. MATH::

    \gamma: I \longrightarrow M

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

    Given a differentiable manifold `M`, a *differentiable curve* curve in
    `M` is a differentiable mapping

    .. MATH::

        \gamma: I \longrightarrow M

    where `I` is an interval of `\RR`.

    This is a Sage *element* class, whose *parent* class is
    :class:`~sage.geometry.manifolds.manifold_homset.ManifoldCurveSet`.

    INPUT:

    - ``parent`` -- the set of curves `\mathrm{Hom}(I, M)` to which the curve
      belongs; this must be an instance of
      :class:`~sage.geometry.manifolds.manifold_homset.ManifoldCurveSet`
    - ``coord_expression`` -- (default: ``None``) dictionary (possibly empty)
      of the functions of the curve parameter `t` expressing the curve in
      various charts of `M`, the keys of the dictionary being the charts and
      the values being lists or tuples of `n` symbolic expressions of `t`,
      where `n` is the dimension of `M`
    - ``name`` -- (default: ``None``) string; symbol given to the curve
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      the curve; if none is provided, ``name`` will be used
    - ``is_diffeomorphism`` -- (default: ``False``) determines whether the
      constructed object is a diffeomorphism; if set to ``True``,
      then `M` must have dimension one.
    - ``is_identity`` -- (default: ``False``) determines whether the
      constructed object is the identity map; if set to ``True``,
      then `M` must be the interval `I`.

    EXAMPLES:

    The lemniscate of Gerono in the 2-dimensional Euclidean plane::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: R.<t> = RealLine()
        sage: c = M.curve({X: [sin(t), sin(2*t)/2]}, (t, 0, 2*pi), name='c') ; c
        Curve 'c' in the 2-dimensional manifold 'M'
        sage: type(c)
        <class 'sage.geometry.manifolds.curve.ManifoldCurveSet_with_category.element_class'>

    Curves are considered as (manifold) morphisms from real intervals to
    differentiable manifolds::

        sage: c.parent()
        Set of Morphisms from Real interval (0, 2*pi) to 2-dimensional manifold 'M' in Category of facade sets
        sage: I = R.open_interval(0, 2*pi)
        sage: c.parent() is Hom(I, M)
        True
        sage: c.domain()
        Real interval (0, 2*pi)
        sage: c.domain() is I
        True
        sage: c.codomain()
        2-dimensional manifold 'M'

    Accordingly, all methods of
    :class:`~sage.geometry.manifolds.diffmapping.DiffMapping` are available
    for them. In particular, the method
    :meth:`~sage.geometry.manifolds.diffmapping.DiffMapping.display`
    shows the coordinate representations in various charts of manifold ``M``::

        sage: c.display()
        c: (0, 2*pi) --> M
           t |--> (x, y) = (sin(t), 1/2*sin(2*t))

    Another mapping method is of course ``__call__``, which returns the image
    of a point in the curve's domain::

        sage: t0 = pi/2
        sage: I(t0)
        point on field R of real numbers
        sage: c(I(t0))
        point on 2-dimensional manifold 'M'
        sage: c(I(t0)).coord(X)
        (1, 0)

    For curves, the value of the parameter, instead of the corresponding point
    in the real line manifold, can be passed directly to the method
    ``__call__``::

        sage: c(t0)
        point 'c(1/2*pi)' on 2-dimensional manifold 'M'
        sage: c(t0).coord(X)
        (1, 0)
        sage: c(t0) == c(I(t0))
        True

    Instead of a dictionary of coordinate expressions, the curve can be
    defined by a single coordinate expression in a given chart::

        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), chart=X, name='c') ; c
        Curve 'c' in the 2-dimensional manifold 'M'
        sage: c.display()
        c: (0, 2*pi) --> M
           t |--> (x, y) = (sin(t), 1/2*sin(2*t))

    Since ``X`` is the default chart on ``M``, it can be omitted::

        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c') ; c
        Curve 'c' in the 2-dimensional manifold 'M'
        sage: c.display()
        c: (0, 2*pi) --> M
           t |--> (x, y) = (sin(t), 1/2*sin(2*t))

    Note that a curve in `M` can also be created as a differentiable mapping
    `I\rightarrow M`::

        sage: c1 = I.diff_mapping(M, coord_functions={X: [sin(t), sin(2*t)/2]},
        ....:                     name='c') ; c1
        Curve 'c' in the 2-dimensional manifold 'M'
        sage: c1.parent() is c.parent()
        True
        sage: c1 == c
        True

    LaTeX symbols representing a curve::

        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi))
        sage: latex(c)
        \mbox{Curve in the 2-dimensional manifold 'M'}
        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
        sage: latex(c)
        c
        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c',
        ....:             latex_name=r'\gamma')
        sage: latex(c)
        \gamma

    The curve's tangent vector field (velocity vector)::

        sage: v = c.tangent_vector_field() ; v
        vector field 'c'' along the Real interval (0, 2*pi) with values on the
         2-dimensional manifold 'M'
        sage: v.display()
        c' = cos(t) d/dx + (2*cos(t)^2 - 1) d/dy

    Its value at `t=\pi`::

        sage: v.at(R(pi))
        tangent vector c' at point on 2-dimensional manifold 'M'
        sage: v.at(R(pi)) in c(R(pi)).tangent_space()
        True
        sage: v.at(R(pi)).display()
        c' = -d/dx + d/dy

    Curves `\RR\rightarrow\RR` can be composed: the operator `\circ` is
    denoted by ``*``::

        sage: f = R.curve(t^2, (t,-oo,+oo))
        sage: g = R.curve(cos(t), (t,-oo,+oo))
        sage: s = g*f ; s
        differentiable mapping from the field R of real numbers to itself
        sage: s.display()
        R --> R
           t |--> cos(t^2)
        sage: s = f*g ; s
        differentiable mapping from the field R of real numbers to itself
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
            identity map 'Id_(0, 2*pi)' of the Real interval (0, 2*pi)
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
            return DiffMapping.__call__(self, t)
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

    def tangent_vector_field(self, name=None, latex_name=None):
        r"""
        Return the tangent vector field to ``self`` (velocity vector).

        INPUT:

        - ``name`` -- (default: ``None``) string; symbol given to the tangent
          vector field; if none is provided, the primed curve symbol (if any)
          will be used
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the tangent vector field; if ``None`` then (i) if ``name`` is
          ``None`` as well, the primed curve LaTeX symbol (if any) will be
          used or (ii) if ``name`` is not ``None``, ``name`` will be used

        OUTPUT:

        - the tangent vector field, as an instance of
          :class:`~sage.geometry.manifolds.vectorfield.VectorField`

        EXAMPLES:

        Tangent vector field to a circle curve in `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: c = M.curve([cos(t), sin(t)], (t, 0, 2*pi), name='c')
            sage: v = c.tangent_vector_field() ; v
            vector field 'c'' along the Real interval (0, 2*pi) with values on
             the 2-dimensional manifold 'R^2'
            sage: v.display()
            c' = -sin(t) d/dx + cos(t) d/dy
            sage: latex(v)
            {c'}
            sage: v.parent()
            free module X((0, 2*pi),c) of vector fields along the Real interval
             (0, 2*pi) mapped into the 2-dimensional manifold 'R^2'

        Value of the tangent vector field for some specific value of the
        curve parameter (`t=\pi`)::

            sage: R(pi) in c.domain()  # pi in (0, 2*pi)
            True
            sage: vp = v.at(R(pi)) ; vp
            tangent vector c' at point on 2-dimensional manifold 'R^2'
            sage: vp.parent() is c(R(pi)).tangent_space()
            True
            sage: vp.display()
            c' = -d/dy


        """
        vmodule = self._domain.vector_field_module(dest_map=self)
        if latex_name is None:
            if name is None:
                if self._latex_name is not None:
                    latex_name = r"{" + self._latex_name + r"'}"
            else:
                latex_name = name
        if name is None and self._name is not None:
            name = self._name + "'"
        resu = vmodule.element_class(vmodule, name=name, latex_name=latex_name)
        canon_chart = self._domain.canonical_chart()
        codom = self._codomain
        dim = codom._manifold._dim
        codom_top_charts = codom._top_charts
        for chart in codom_top_charts:
            try:
                jacob = self.differential_functions(canon_chart, chart)
            except ValueError:
                continue
            frame = vmodule.basis(from_frame=chart.frame())
            resu.add_comp(frame)[:, canon_chart] = [jacob[i][0]
                                                          for i in range(dim)]
        return resu


    def plot(self, chart=None, prange=None, include_end_point=(True, True),
             end_point_offset=(0.001, 0.001), max_value=8, plot_coords=None,
             mapping=None, parameters=None, color='red',  style='-',
             thickness=1, plot_points=75, label_axes=True,
             aspect_ratio='automatic'):
        r"""
        Plot of the curve in terms of a given chart of its codomain.

        INPUT:

        - ``chart`` -- (default: ``None``) the chart in terms of which the plot
          is performed; if ``None``, the default chart of the codomain of
          curve (or of the curve composed by ``mapping``) is used
        - ``prange`` -- (default: ``None``) range of the curve parameter for
          the plot; if ``None``, the entire parameter range declared during the
          curve construction is considered (with -Infinity
          replaced by ``-max_value`` and +Infinity by ``max_value``)
        - ``include_end_point`` -- (default: ``(True, True)``) determines
          whether the end points of ``prange`` are included in the plot
        - ``end_point_offset`` -- (default: ``(0.001, 0.001)``) offsets from
          the end points when they are not included in the plot: if
          ``include_end_point[0] == False``, the minimal value of the curve
          parameter used for the plot is ``prange[0] + end_point_offset[0]``,
          while if ``include_end_point[1] == False``, the maximal value is
          ``prange[1] - end_point_offset[1]``.
        - ``max_value`` -- (default: 8) numerical value substituted to
          +Infinity if the latter is the upper bound of the parameter range;
          similarly ``-max_value`` is the numerical valued substituted for
          -Infinity
        - ``ambient_coords`` -- (default: ``None``) tuple containing the 2 or 3
          coordinates of ``chart`` in terms of which the plot is
          performed; if ``None``, all the coordinates of ``chart`` are
          considered
        - ``mapping`` -- (default: ``None``) differentiable mapping (instance
          of :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`)
          providing the link between ``self`` and the chart ``chart`` (cf.
          above); if ``None``, chart is supposed to be defined on the
          codomain of the curve ``self``.
        - ``parameters`` -- (default: ``None``) dictionary giving the numerical
          values of the parameters that may appear in the coordinate expression
          of ``self``
        - ``color`` -- (default: 'red') color of the drawn curve
        - ``style`` -- (default: '-') color of the drawn curve; NB: ``style``
          is effective only for 2D plots
        - ``thickness`` -- (default: 1) thickness of the drawn curve
        - ``plot_points`` -- (default: 75) number of points to plot the curve
        - ``label_axes`` -- (default: ``True``) boolean determining whether the
          labels of the coordinate axes of ``chart`` shall be added to the
          graph; can be set to ``False`` if the graph is 3D and must be
          superposed with another graph.
        - ``aspect_ratio`` -- (default: 'automatic') aspect ratio of the plot;
          the default value ('automatic') applies only for 2D plots; for
          3D plots, the default value is ``1`` instead.

        OUTPUT:

        - a graphic object, either an instance of
          :class:`~sage.plot.graphics.Graphics` for a 2D plot (i.e. based on
          2 coordinates of ``chart``) or an instance of
          :class:`~sage.plot.plot3d.base.Graphics3d` for a 3D plot (i.e.
          based on 3 coordinates of ``chart``)

        EXAMPLES:

        Plot of the lemniscate of Gerono::

            sage: R2 = Manifold(2, 'R^2')
            sage: X.<x,y> = R2.chart()
            sage: R.<t> = RealLine()
            sage: c = R2.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
            sage: c.plot()  # 2D plot
            Graphics object consisting of 1 graphics primitive

        Plot for a subinterval of the curve's domain::

            sage: c.plot(prange=(0,pi))
            Graphics object consisting of 1 graphics primitive

        Plot with various options::

            sage: c.plot(color='green', style=':', thickness=3, aspect_ratio=1)
            Graphics object consisting of 1 graphics primitive

        Plot via a mapping to another manifold: loxodrome of a sphere viewed
        in `\RR^3`::

            sage: S2 = Manifold(2, 'S^2')
            sage: U = S2.open_subset('U')
            sage: XS.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: R3 = Manifold(3, 'R^3')
            sage: X3.<x,y,z> = R3.chart()
            sage: F = S2.diff_mapping(R3, {(XS, X3): [sin(th)*cos(ph),
            ....:                     sin(th)*sin(ph), cos(th)]}, name='F')
            sage: F.display()
            F: S^2 --> R^3
            on U: (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: c = S2.curve([2*atan(exp(-t/10)), t], (t, -oo, +oo), name='c')
            sage: graph_c = c.plot(mapping=F, max_value=40,
            ....:                  plot_points=200, thickness=2, label_axes=False)  # 3D plot
            sage: graph_S2 = XS.plot(X3, mapping=F, nb_values=11, color='black') # plot of the sphere
            sage: show(graph_c + graph_S2) # the loxodrome + the sphere

        Example of use of the argument ``parameters``: we define a curve with
        some symbolic parameters ``a`` and ``b``::

            sage: a, b = var('a b')
            sage: c = R2.curve([a*cos(t) + b, a*sin(t)], (t, 0, 2*pi), name='c')

        To make a plot, we set spectific values for ``a`` and ``b`` by means
        of the Python dictionary ``parameters``::

            sage: c.plot(parameters={a: 2, b: -3}, aspect_ratio=1)
            Graphics object consisting of 1 graphics primitive

        """
        from sage.rings.infinity import Infinity
        from sage.misc.functional import numerical_approx
        from sage.plot.graphics import Graphics
        from sage.plot.line import line
        from sage.geometry.manifolds.chart import Chart
        from sage.geometry.manifolds.utilities import set_axes_labels
        #
        # The "effective" curve to be plotted
        #
        if mapping is None:
            eff_curve = self
        else:
            eff_curve = mapping * self
        #
        # The chart w.r.t. which the curve is plotted
        #
        if chart is None:
            chart = eff_curve._codomain.default_chart()
        elif not isinstance(chart, Chart):
            raise TypeError("{} is not a chart".format(chart))
        #
        # Coordinates of the above chart w.r.t. which the curve is plotted
        #
        if plot_coords is None:
            plot_coords = chart[:]  # all chart coordinates are used
        n_pc = len(plot_coords)
        if n_pc != 2 and n_pc !=3:
            raise ValueError("Bad number of plot coordinates: {}".format(n_cp))
        ind_pc = [chart[:].index(pc) for pc in plot_coords] # indices of plot
                                                            # coordinates
        #
        # Parameter range for the plot
        #
        if prange is None:
            prange = (self._domain.lower_bound(), self._domain.upper_bound())
        elif not isinstance(prange, (tuple, list)):
            raise TypeError("{} is neither a tuple nor a list".format(prange))
        elif len(prange) != 2:
            raise ValueError("the argument prange must be a tuple/list " +
                             "of 2 elements")
        tmin = prange[0]
        tmax = prange[1]
        if tmin == -Infinity:
            tmin = -max_value
        elif not include_end_point[0]:
            tmin = tmin + end_point_offset[0]
        if tmax == Infinity:
            tmax = max_value
        elif not include_end_point[1]:
            tmax = tmax - end_point_offset[1]
        tmin = numerical_approx(tmin)
        tmax = numerical_approx(tmax)
        #
        # The coordinate expression of the effective curve
        #
        canon_chart = self._domain.canonical_chart()
        transf = None
        for chart_pair in eff_curve._coord_expression:
            if chart_pair == (canon_chart, chart):
                transf = eff_curve._coord_expression[chart_pair]
                break
        else:
            # Search for a subchart
            for chart_pair in eff_curve._coord_expression:
                for schart in chart._subcharts:
                    if chart_pair == (canon_chart, schart):
                        transf = eff_curve._coord_expression[chart_pair]
        if transf is None:
            raise ValueError("No expression has been found for " +
                              "{} in terms of {}".format(self, format))
        #
        # List of points for the plot curve
        #
        plot_curve = []
        dt = (tmax - tmin) / (plot_points - 1)
        t = tmin
        if parameters is None:
            for i in range(plot_points):
                x = transf(t, simplify=False)
                plot_curve.append( [numerical_approx(x[j]) for j in ind_pc] )
                t += dt
        else:
             for i in range(plot_points):
                x = transf(t, simplify=False)
                plot_curve.append(
                               [numerical_approx( x[j].substitute(parameters) )
                                for j in ind_pc] )
                t += dt
        #
        # The plot
        #
        resu = Graphics()
        resu += line(plot_curve, color=color, linestyle=style,
                     thickness=thickness)
        if n_pc==2:  # 2D graphic
            resu.set_aspect_ratio(aspect_ratio)
            if label_axes:
                resu.axes_labels([r'$'+latex(pc)+r'$' for pc in plot_coords])
        else: # 3D graphic
            if aspect_ratio == 'automatic':
                aspect_ratio = 1
            resu.aspect_ratio(aspect_ratio)
            if label_axes:
                labels = [str(pc) for pc in plot_coords]
                resu = set_axes_labels(resu, *labels)
        return resu
