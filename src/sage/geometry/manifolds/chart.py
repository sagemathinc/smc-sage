r"""
Coordinate charts

Five classes deal with coordinates on a differentiable manifold over `\RR`:

* :class:`Chart`: charts on a manifold
* :class:`FunctionChart`: real-valued functions of the coordinates of a
  given chart
* :class:`ZeroFunctionChart`: the null function of the coordinates of a
  given chart
* :class:`MultiFunctionChart`: functions from a chart codomain to `\RR^m`
* :class:`CoordChange`: transition maps between charts

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version

REFERENCES:

- Chap. 1 of J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer
  (New York) (2013)

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

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import RingElement
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.misc.latex import latex
from sage.geometry.manifolds.domain import ManifoldOpenSubset
from sage.geometry.manifolds.utilities import simplify_chain, ExpressionNice

class Chart(UniqueRepresentation, SageObject):
    r"""
    Chart on a manifold.

    Given a manifold `M` of dimension `n`, a *chart* is a pair `(U,\varphi)`,
    where `U` is an open subset of `M` and
    `\varphi: U \rightarrow V \subset \RR^n` is a homeomorphism from `U` to
    an open subset `V` of `\RR^n`.

    The components `(x^1,\ldots,x^n)` of `\varphi`, defined by
    `\varphi(p) = (x^1(p),\ldots,x^n(p))`, are called the *coordinates* of the
    chart `(U,\varphi)`.

    INPUT:

    - ``domain`` -- open subset `U` on which the chart is defined (must be
      an instance of :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`)
    - ``coordinates`` -- (default: '' (empty string)) single string defining
      the coordinate symbols and ranges: the coordinates are separated by ' '
      (space) and each coordinate has at most three fields, separated by ':':

        1. The coordinate symbol (a letter or a few letters)
        2. (optional) The interval `I` defining the coordinate range: if not
           provided, the coordinate is assumed to span all `\RR`; otherwise
           `I` must be provided in the form ``(a,b)`` (or equivalently
           ``]a,b[``). The bounds ``a`` and ``b`` can be ``+/-Infinity``,
           ``Inf``, ``infinity``, ``inf`` or ``oo``.
           For *singular* coordinates, non-open intervals such as ``[a,b]`` and
           ``(a,b]`` (or equivalently ``]a,b]``) are allowed.
           Note that the interval declaration must not contain any space
           character.
        3. (optional) The LaTeX spelling of the coordinate; if not provided the
           coordinate symbol given in the first field will be used.

      The order of the fields 2 and 3 does not matter and each of them can be
      omitted.
      If it contains any LaTeX expression, the string ``coordinates`` must be
      declared with the prefix 'r' (for "raw") to allow for a proper treatment
      of the backslash character (see examples below).
      If no interval range and no LaTeX spelling is to be set for any
      coordinate, the argument ``coordinates`` can be omitted when the
      shortcut operator ``<,>`` is used via Sage preparser (see examples below)
    - ``names`` -- (default: ``None``) unused argument, except if
      ``coordinates`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranted if the shortcut operator
      ``<,>`` is used).

    EXAMPLES:

    Cartesian coordinates on `\RR^3`::

        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart = M.chart('x y z') ; c_cart
        chart (R^3, (x, y, z))
        sage: type(c_cart)
        <class 'sage.geometry.manifolds.chart.Chart'>

    To have the coordinates accessible as global variables, one has to set::

        sage: (x,y,z) = c_cart[:]

    However, a shortcut is to use the declarator ``<x,y,z>`` in the left-hand
    side of the chart declaration (there is then no need to pass the string
    'x y z' to chart())::

        sage: Manifold._clear_cache_() # for doctests only
        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = M.chart() ; c_cart
        chart (R^3, (x, y, z))

    The coordinates are then immediately accessible::

        sage: y
        y
        sage: y is c_cart[2]
        True

    The trick is performed by Sage preparser::

        sage: preparse("c_cart.<x,y,z> = M.chart()")
        "c_cart = M.chart(names=('x', 'y', 'z',)); (x, y, z,) = c_cart._first_ngens(3)"

    Note that x, y, z declared in ``<x,y,z>`` are mere Python variable names
    and do not have to coincide with the coordinate symbols; for instance,
    one may write::

        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x1,y1,z1> = M.chart('x y z') ; c_cart
        chart (R^3, (x, y, z))

    Then y is not known as a global variable and the corresponding coordinate
    is accessible only through the global variable y1::

        sage: y1
        y
        sage: y1 is c_cart[2]
        True

    However, having the name of the Python variable coincide with the
    coordinate symbol is quite convenient; so it is recommended to declare::

        sage: Manifold._clear_cache_() # for doctests only
        sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)
        sage: c_cart.<x,y,z> = M.chart()

    Spherical coordinates on the subset `U` of `\RR^3` that is the
    complement of the half-plane `\{y=0, x\geq 0\}`::

        sage: U = M.open_subset('U')
        sage: c_spher.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi') ; c_spher
        chart (U, (r, th, ph))

    Note the prefix 'r' for the string defining the coordinates in the arguments of ``Chart``.

    Coordinates are some Sage symbolic variables::

        sage: type(th)
        <type 'sage.symbolic.expression.Expression'>
        sage: latex(th)
        {\theta}
        sage: assumptions(th)
        [th is real, th > 0, th < pi]

    Coordinate are also accessible by their indices::

        sage: x1 = c_spher[1]; x2 = c_spher[2]; x3 = c_spher[3]
        sage: print x1, x2, x3
        r th ph
        sage: (x1, x2, x3) == (r, th, ph)
        True

    The full set of coordinates is obtained by means of the operator [:]::

        sage: c_cart[:]
        (x, y, z)
        sage: c_spher[:]
        (r, th, ph)

    Let us check that the declared coordinate ranges have been taken into
    account::

        sage: c_cart.coord_range()
        x: (-oo, +oo); y: (-oo, +oo); z: (-oo, +oo)
        sage: c_spher.coord_range()
        r: (0, +oo); th: (0, pi); ph: (0, 2*pi)
        sage: bool(th>0 and th<pi)
        True
        sage: assumptions()  # list all current symbolic assumptions
        [x is real, y is real, z is real, r is real, r > 0, th is real,
         th > 0, th < pi, ph is real, ph > 0, ph < 2*pi]

    The coordinate ranges are used for simplifications::

        sage: simplify(abs(r)) # r has been declared to lie in the interval (0,+oo)
        r
        sage: simplify(abs(x)) # no positive range has been declared for x
        abs(x)
        sage: from sage.geometry.manifolds.utilities import simplify_abs_trig
        sage: simplify_abs_trig(abs(sin(th)))  # sin(th) is always positive
        sin(th)
        sage: simplify_abs_trig(abs(sin(ph)))  # sin(ph) is not
        abs(sin(ph))

    Each constructed chart is automatically added to the manifold's atlas::

        sage: M.atlas()
        [chart (R^3, (x, y, z)), chart (U, (r, th, ph))]

    and to the atlas of its domain::

        sage: U.atlas()
        [chart (U, (r, th, ph))]

    Manifold subsets may have a default chart, which, unless changed via the
    method
    :meth:`~sage.geometry.manifolds.domain.ManifoldSubset.set_default_chart`,
    is the first defined chart on the subset (or on a open subset of it)::

        sage: M.default_chart()
        chart (R^3, (x, y, z))
        sage: U.default_chart()
        chart (U, (r, th, ph))

    The chart map `\varphi` acting on a point is obtained by means of the
    call operator, i.e. the operator ``()``::

        sage: p = M.point((1,0,-2)) ; p
        point on 3-dimensional manifold 'R^3'
        sage: c_cart(p)
        (1, 0, -2)
        sage: c_cart(p) == p.coord(c_cart)
        True
        sage: q = M.point((2,pi/2,pi/3), c_spher) # point defined by its spherical coordinates
        sage: c_spher(q)
        (2, 1/2*pi, 1/3*pi)
        sage: c_spher(q) == q.coord(c_spher)
        True
        sage: a = U.point((1,pi/2,pi)) # the default coordinates on U are the spherical ones
        sage: c_spher(a)
        (1, 1/2*pi, pi)
        sage: c_spher(a) == a.coord(c_spher)
        True

    Cartesian coordinates on U as an example of chart construction with
    coordinate restrictions: since U is the complement of the half-plane
    `\{y=0, x\geq 0\}`, we must have `y\not=0` or `x<0` on U. Accordingly,
    we set::

        sage: c_cartU.<x,y,z> = U.chart()
        sage: c_cartU.add_restrictions((y!=0, x<0)) # the tuple (y!=0, x<0) means y!=0 or x<0
        sage: # c_cartU.add_restrictions([y!=0, x<0]) would have meant y!=0 AND x<0
        sage: U.atlas()
        [chart (U, (r, th, ph)), chart (U, (x, y, z))]
        sage: M.atlas()
        [chart (R^3, (x, y, z)), chart (U, (r, th, ph)), chart (U, (x, y, z))]
        sage: c_cartU.valid_coordinates(-1,0,2)
        True
        sage: c_cartU.valid_coordinates(1,0,2)
        False
        sage: c_cart.valid_coordinates(1,0,2)
        True

    Each constructed chart has its zero function, mapping the coordinates to 0;
    this zero function is an instance of :class:`ZeroFunctionChart`::

        sage: c_spher.zero_function()
        0
        sage: c_spher.zero_function().display()
        (r, th, ph) |--> 0
        sage: type(c_spher.zero_function())
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        sage: c_cart.zero_function()
        0
        sage: c_cart.zero_function().display()
        (x, y, z) |--> 0

    Of course, zero function defined on different charts are not equal::

        sage: c_cart.zero_function() == c_spher.zero_function()
        False

    Chart grids can be drawn in 2D or 3D graphics thanks to the method
    :meth:`plot`.

    """
    def __init__(self, domain, coordinates='', names=None):
        from sage.symbolic.ring import SR
        from sage.symbolic.assumptions import assume
        from vectorframe import CoordFrame
        if not isinstance(domain, ManifoldOpenSubset):
            raise TypeError("The first argument must be an open subset.")
        if coordinates == '':
            for x in names:
                coordinates += x + ' '
        self._manifold = domain._manifold
        self._domain = domain
        # Treatment of the coordinates:
        if ' ' in coordinates:
            coord_list = coordinates.split()
        else:
            coord_list = [coordinates]
        n = self._manifold._dim
        if len(coord_list) != n:
            raise ValueError("The list of coordinates must contain " + \
                              str(n) + " elements.")
        xx_list = [] # will contain the coordinates as Sage symbolic variables
        bounds_list = [] # will contain the coordinate bounds
        for coord_field in coord_list:
            coord_properties = coord_field.split(':')
            coord_symb = coord_properties[0].strip() # the coordinate symbol
            # default values, possibly redefined below:
            coord_latex = None
            xmin = -Infinity ; xmin_included = False
            xmax = +Infinity ; xmax_included = False
            # scan of the properties other than the symbol:
            for prop in coord_properties[1:]:
                prop1 = prop.strip()
                delim_min = prop1[0]
                if delim_min in ['[', ']', '(']:
                    # prop1 is the coordinate's range
                    xmin_str, xmax_str = prop1[1:len(prop1)-1].split(',')
                    if xmin_str not in ['-inf', '-Inf', '-infinity',
                                        '-Infinity', '-oo']:
                        xmin = SR(xmin_str)
                        xmin_included = ( delim_min == '[' )
                    if xmax_str not in ['inf', '+inf', 'Inf', '+Inf',
                                        'infinity', '+infinity', 'Infinity',
                                        '+Infinity', 'oo', '+oo']:
                        xmax = SR(xmax_str)
                        xmax_included = ( prop1[-1] == ']' )
                else:
                    # prop1 is the coordinate's LaTeX symbol
                    coord_latex = prop1
            # Construction of the coordinate as some Sage's symbolic variable:
            coord_var = SR.var(coord_symb, domain='real',
                               latex_name=coord_latex)
            assume(coord_var, 'real')
            if xmin != -Infinity:
                if xmin_included:
                    assume(coord_var >= xmin)
                else:
                    assume(coord_var > xmin)
            if xmax != Infinity:
                if xmax_included:
                    assume(coord_var <= xmax)
                else:
                    assume(coord_var < xmax)
            xx_list.append(coord_var)
            bounds_list.append(((xmin, xmin_included), (xmax, xmax_included)))
        self._xx = tuple(xx_list)
        self._bounds = tuple(bounds_list)
        # End of the treatment of the coordinates

        # Additional restrictions on the coordinates
        self._restrictions = []  # to be set with method add_restrictions()

        # The chart is added to the domain's atlas, as well as to all the
        # atlases of the domain's supersets; moreover the fist defined chart
        # is considered as the default chart
        for sd in self._domain._supersets:
            # the chart is added in the top charts only if its coordinates have
            # not been used:
            for chart in sd._atlas:
                if self._xx == chart._xx:
                    break
            else:
                sd._top_charts.append(self)
            sd._atlas.append(self)
            if sd._def_chart is None:
                sd._def_chart = self
        # The chart is added to the list of the domain's covering charts:
        self._domain._covering_charts.append(self)
        # Initialization of the set of charts that are restrictions of the
        # current chart to subsets of the chart domain:
        self._subcharts = set([self])
        # Initialization of the set of charts which the current chart is a
        # restriction of:
        self._supercharts = set([self])
        #
        self._dom_restrict = {} # dict. of the restrictions of self to
                                # subsets of self._domain, with the
                                # subsets as keys
        # The null function of the coordinates:
        self._zero_function = ZeroFunctionChart(self)
        # Construction of the coordinate frame associated to the chart:
        self._frame = CoordFrame(self)
        self._coframe = self._frame._coframe
        # Expression in self of the zero scalar fields of open sets containing
        # the domain of self:
        for dom in self._domain._supersets:
            if hasattr(dom, '_zero_scalar_field'):
                # dom is an open set
                dom._zero_scalar_field._express[self] = self._zero_function

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = 'chart ' + \
                      '(' + self._domain._name + ', ' + str(self._xx) + ')'
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        description = r'\left(' + latex(self._domain).strip() + ',('
        n = len(self._xx)
        for i in range(n-1):
            description += latex(self._xx[i]).strip() + ', '
        description += latex(self._xx[n-1]).strip() + r')\right)'
        return description

    def _first_ngens(self, n):
        r"""
        Return the list of coordinates.

        This is useful only for the use of Sage preparser::

            sage: preparse("c_cart.<x,y,z> = M.chart()")
            "c_cart = M.chart(names=('x', 'y', 'z',)); (x, y, z,) = c_cart._first_ngens(3)"

        """
        return self[:]


    def __getitem__(self, i):
        r"""
        Access to the coordinates.

        INPUT:

        - ``i`` -- index of the coordinate; if [:] all the coordinates
            are returned

        OUTPUT:

        - the coordinate of index ``i`` or all the coordinates (as a tuple) if
          ``i`` is [:]
        """
        if isinstance(i, slice):
            return self._xx
        else:
            return self._xx[i-self._manifold._sindex]

    def __call__(self, point):
        r"""
        Return the coordinates of a given point.
        """
        return point.coord(self)

    def domain(self):
        r"""
        Return the open subset on which ``self`` is defined.
        """
        return self._domain

    def frame(self):
        r"""
        Return the vector frame (coordinate frame) associated with the chart.

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.vectorframe.CoordFrame`
          representing the coordinate frame.

        EXAMPLE:

        Coordinate frame associated with some chart on a 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_xy.frame()
            coordinate frame (M, (d/dx,d/dy))
            sage: type(c_xy.frame())
            <class 'sage.geometry.manifolds.vectorframe.CoordFrame'>

        Check that c_xy.frame() is indeed the coordinate frame associated with
        (x,y)::

            sage: ex = c_xy.frame()[0] ; ex
            vector field 'd/dx' on the 2-dimensional manifold 'M'
            sage: ey = c_xy.frame()[1] ; ey
            vector field 'd/dy' on the 2-dimensional manifold 'M'
            sage: ex(M.scalar_field(x)).display()
            M --> R
            (x, y) |--> 1
            sage: ex(M.scalar_field(y)).display()
            M --> R
            (x, y) |--> 0
            sage: ey(M.scalar_field(x)).display()
            M --> R
            (x, y) |--> 0
            sage: ey(M.scalar_field(y)).display()
            M --> R
            (x, y) |--> 1

        """
        return self._frame

    def coframe(self):
        r"""
        Return the coframe (basis of coordinate differentials) associated
        with the chart.

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.vectorframe.CoordCoFrame`
          representing the coframe.

        EXAMPLE:

        Coordinate coframe associated with some chart on a 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_xy.coframe()
            coordinate coframe (M, (dx,dy))
            sage: type(c_xy.coframe())
            <class 'sage.geometry.manifolds.vectorframe.CoordCoFrame'>

        Check that c_xy.coframe() is indeed the coordinate coframe associated
        with (x,y)::

            sage: dx = c_xy.coframe()[0] ; dx
            1-form 'dx' on the 2-dimensional manifold 'M'
            sage: dy = c_xy.coframe()[1] ; dy
            1-form 'dy' on the 2-dimensional manifold 'M'
            sage: ex = c_xy.frame()[0] ; ex
            vector field 'd/dx' on the 2-dimensional manifold 'M'
            sage: ey = c_xy.frame()[1] ; ey
            vector field 'd/dy' on the 2-dimensional manifold 'M'
            sage: dx(ex).display()
             dx(d/dx): M --> R
               (x, y) |--> 1
            sage: dx(ey).display()
            dx(d/dy): M --> R
               (x, y) |--> 0
            sage: dy(ex).display()
            dy(d/dx): M --> R
               (x, y) |--> 0
            sage: dy(ey).display()
            dy(d/dy): M --> R
               (x, y) |--> 1

        """
        return self._coframe


    def coord_bounds(self, i=None):
        r"""
        Return the lower and upper bounds of the range of a coordinate.

        For a nicely formatted output, use :meth:`coord_range` instead.

        INPUT:

        - ``i`` -- (default: ``None``)  index of the coordinate; if ``None``,
          the bounds of all the coordinates are returned

        OUTPUT:

        - the coordinate bounds as the tuple
          ((xmin, min_included), (xmax, max_included))
          where

          - xmin is the coordinate lower bound
          - min_included is a Boolean, indicating whether the coordinate can
            take the value xmin, i.e. xmin is a strict lower bound iff
            min_included is False.
          - xmin is the coordinate upper bound
          - max_included is a Boolean, indicating whether the coordinate can
            take the value xmax, i.e. xmax is a strict upper bound iff
            max_included is False.

        EXAMPLES:

        Some coordinate bounds on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart('x y:[0,1)')
            sage: c_xy.coord_bounds(0)  # x in (-oo,+oo) (the default)
            ((-Infinity, False), (+Infinity, False))
            sage: c_xy.coord_bounds(1)  # y in [0,1)
            ((0, True), (1, False))
            sage: c_xy.coord_bounds()
            (((-Infinity, False), (+Infinity, False)), ((0, True), (1, False)))
            sage: c_xy.coord_bounds() == (c_xy.coord_bounds(0), c_xy.coord_bounds(1))
            True

        The coordinate bounds can also be recovered via the method
        :meth:`coord_range`::

            sage: c_xy.coord_range()
            x: (-oo, +oo); y: [0, 1)
            sage: c_xy.coord_range(y)
            y: [0, 1)

        or via Sage's function
        :func:`sage.symbolic.assumptions.assumptions`::

            sage: assumptions(x)
            [x is real]
            sage: assumptions(y)
            [y is real, y >= 0, y < 1]

        """
        if i is None:
            return self._bounds
        else:
            return self._bounds[i-self._manifold._sindex]

    def coord_range(self, xx=None):
        r"""
        Display the range of a coordinate (or all coordinates), as an
        interval.

        INPUT:

        - ``xx`` -- (default: ``None``) symbolic expression corresponding to a
          coordinate of ``self``; if ``None``, the ranges of all coordinates
          are displayed.

        EXAMPLES:

        Ranges of coordinates on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: X.coord_range()
            x: (-oo, +oo); y: (-oo, +oo)
            sage: X.coord_range(x)
            x: (-oo, +oo)
            sage: U = M.open_subset('U', coord_def={X: [x>1, y<pi]})
            sage: XU = X.restrict(U)  # restriction of chart X to U
            sage: XU.coord_range()
            x: (1, +oo); y: (-oo, pi)
            sage: XU.coord_range(x)
            x: (1, +oo)
            sage: XU.coord_range(y)
            y: (-oo, pi)

        The output is LaTeX-formatted for the notebook::

            sage: latex(XU.coord_range(y))
            y :\ \left( -\infty, \pi \right)

        """
        from sage.tensor.modules.format_utilities import FormattedExpansion
        resu_txt = ""
        resu_latex = ""
        if xx is None:
            for x in self._xx:
                if resu_txt != "":
                    resu_txt += "; "
                    resu_latex += r";\quad "
                resu_txt, resu_latex = self._display_coord_range(x, resu_txt,
                                                                 resu_latex)
        else:
            resu_txt, resu_latex = self._display_coord_range(xx, resu_txt,
                                                             resu_latex)
        return FormattedExpansion(resu_txt, resu_latex)

    def _display_coord_range(self, xx, rtxt, rlatex):
        r"""
        Helper function for the display of a coordinate range

        INPUT:

        - ``xx`` -- symbolic expression corresponding to a
          coordinate of ``self``
        - ``rtxt`` -- string; initial state for the text version of the range
        - ``rlatex`` -- string; initial state for the LaTeX version of the
          range

        OUTPUT:

        - tuple of 2 strings containing the final state for the text and
          LaTeX versions of the range

        """
        ind = self._xx.index(xx)
        bounds = self._bounds[ind]
        rtxt += "{}: ".format(xx)
        rlatex += latex(xx) + r":\ "
        if bounds[0][1]:
            rtxt += "["
            rlatex += r"\left["
        else:
            rtxt += "("
            rlatex += r"\left("
        xmin = bounds[0][0]
        if xmin == -Infinity:
            rtxt += "-oo, "
            rlatex += r"-\infty,"
        else:
            rtxt += "{}, ".format(xmin)
            rlatex += latex(xmin) + ","
        xmax = bounds[1][0]
        if xmax == Infinity:
            rtxt += "+oo"
            rlatex += r"+\infty"
        else:
            rtxt += "{}".format(xmax)
            rlatex += latex(xmax)
        if bounds[1][1]:
            rtxt += "]"
            rlatex += r"\right]"
        else:
            rtxt += ")"
            rlatex += r"\right)"
        return rtxt, rlatex

    def add_restrictions(self, restrictions):
        r"""
        Add some restrictions on the coordinates.

        INPUT:

        - ``restrictions`` -- list of restrictions on the
          coordinates, in addition to the ranges declared by the intervals
          specified in the chart constructor.
          A restriction can be any symbolic equality or inequality involving the
          coordinates, such as x>y or x^2+y^2 != 0. The items of the list
          ``restrictions`` are combined with the ``and`` operator; if some
          restrictions are to be combined with the ``or`` operator instead, they
          have to be passed as a tuple in some single item of the list
          ``restrictions``. For example, ``restrictions`` = [x>y, (x!=0, y!=0),
          z^2<x] means (x>y) and ((x!=0) or (y!=0)) and (z^2<x). If the list
          ``restrictions`` contains only one item, this item can be passed as
          such, i.e. writing x>y instead of the single element list [x>y].

        EXAMPLES:

        Cartesian coordinates on the open unit disc in $\RR^2$::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M') # the open unit disc
            sage: X.<x,y> = M.chart()
            sage: X.add_restrictions(x^2+y^2<1)
            sage: X.valid_coordinates(0,2)
            False
            sage: X.valid_coordinates(0,1/3)
            True

        The restrictions are transmitted to subcharts::

            sage: A = M.open_subset('A') # annulus 1/2 < r < 1
            sage: X_A = X.restrict(A, x^2+y^2 > 1/4)
            sage: X_A._restrictions
            [x^2 + y^2 < 1, x^2 + y^2 > (1/4)]
            sage: X_A.valid_coordinates(0,1/3)
            False
            sage: X_A.valid_coordinates(2/3,1/3)
            True

        If appropriate, the restrictions are transformed into bounds on
        the coordinate ranges::

            sage: U = M.open_subset('U')
            sage: X_U = X.restrict(U)
            sage: X_U.coord_range()
            x: (-oo, +oo); y: (-oo, +oo)
            sage: X_U.add_restrictions([x<0, y>1/2])
            sage: X_U.coord_range()
            x: (-oo, 0); y: (1/2, +oo)

        """
        import operator
        if not isinstance(restrictions, list):
            # case of a single condition or conditions to be combined by "or"
            restrictions = [restrictions]
        self._restrictions.extend(restrictions)
        # Update of the coordinate bounds from the restrictions:
        bounds = list(self._bounds) # convert to a list for modifications
        new_restrictions = []
        for restrict in self._restrictions:
            restrict_used = False # determines whether restrict is used
                                  # to set some coordinate bound
            if not isinstance(restrict, tuple): # case of 'or' conditions
                                                # excluded
                operands = restrict.operands()
                left = operands[0]
                right = operands[1]
                right_var = right.variables()
                if left in self._xx:
                    # the l.h.s. of the restriction is a single
                    # coordinate
                    right_coord = [coord for coord in self._xx
                                   if coord in right_var]
                    if not right_coord:
                        # there is no other coordinate in the r.h.s.
                        ind = self._xx.index(left)
                        left_bounds = list(bounds[ind])
                        oper = restrict.operator()
                        oinf = left_bounds[0][0] # old coord inf
                        osup = left_bounds[1][0] # old coord sup
                        if oper == operator.lt:
                            if osup == Infinity or right <= osup:
                                left_bounds[1] = (right, False)
                                restrict_used = True
                        elif oper == operator.le:
                            if osup == Infinity or right < osup:
                                left_bounds[1] = (right, True)
                                restrict_used = True
                        elif oper == operator.gt:
                            if oinf == -Infinity or right >= oinf:
                                left_bounds[0] = (right, False)
                                restrict_used = True
                        elif oper == operator.ge:
                            if oinf == -Infinity or right > oinf:
                                left_bounds[0] = (right, True)
                                restrict_used = True
                        bounds[ind] = tuple(left_bounds)
            if not restrict_used:
                # if restrict has not been used to set a coordinate bound
                # it is maintained in the list of restrictions:
                new_restrictions.append(restrict)
        self._bounds = tuple(bounds)
        self._restrictions = new_restrictions


    def restrict(self, subset, restrictions=None):
        r"""
        Return the restriction of ``self`` to some subset.

        If ``self`` is the chart `(U,\varphi)`, a restriction (or subchart)
        is a chart `(V,\psi)` such that `V\subset U` and `\psi = \varphi |_V`.

        If such subchart has not been defined yet, it is constructed here.

        The coordinates of the subchart bare the same names as the coordinates
        of the mother chart.

        INPUT:

        - ``subset`` -- open subset `V` of the chart domain `U` (must
          be an instance of
          :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`)
        - ``restrictions`` -- (default: ``None``) list of coordinate restrictions
          defining the subset `V`.
          A restriction can be any symbolic equality or
          inequality involving the coordinates, such as x>y or x^2+y^2 != 0.
          The items of the list ``restrictions`` are combined with the ``and``
          operator; if some restrictions are to be combined with the ``or``
          operator instead, they have to be passed as a tuple in some single
          item of the list ``restrictions``. For example, ``restrictions``
          being [x>y, (x!=0, y!=0), z^2<x] means (x>y) and ((x!=0) or (y!=0))
          and (z^2<x). If the list ``restrictions`` contains only one item,
          this item can be passed as such, i.e. writing x>y instead of the
          single element list [x>y]. Note that the argument ``restrictions``
          can be omitted if the subchart has been already initialized by a
          previous call.

        OUTPUT:

        - chart `(V,\psi)`, as an instance of :class:`Chart`.

        EXAMPLES:

        Cartesian coordinates on the unit open disc in `\RR^2` as a subchart
        of the global Cartesian coordinates::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: D = M.open_subset('D') # the unit open disc
            sage: c_cart_D = c_cart.restrict(D, x^2+y^2<1)
            sage: p = M.point((1/2, 0))
            sage: p in D
            True
            sage: q = M.point((1, 2))
            sage: q in D
            False

        Cartesian coordinates on the annulus `1<\sqrt{x^2+y^2}<2`::

            sage: A = M.open_subset('A')
            sage: c_cart_A = c_cart.restrict(A, [x^2+y^2>1, x^2+y^2<4])
            sage: p in A, q in A
            (False, False)
            sage: a = M.point((3/2,0))
            sage: a in A
            True

        """
        if subset == self._domain:
            return self
        if subset not in self._dom_restrict:
            if not subset.is_subset(self._domain):
                raise ValueError("The specified subset is not a subset " +
                                 "of the domain of definition of the chart.")
            coordinates = ""
            for coord in self._xx:
                coordinates += repr(coord) + ' '
            res = Chart(subset, coordinates)
            res._bounds = self._bounds
            res._restrictions.extend(self._restrictions)
            # The coordinate restrictions are added to the result chart and
            # possibly transformed into coordinate bounds:
            if restrictions is not None:
                res.add_restrictions(restrictions)
            # Update of supercharts and subcharts:
            res._supercharts.update(self._supercharts)
            for schart in self._supercharts:
                schart._subcharts.add(res)
                schart._dom_restrict[subset] = res
            # Update of superframes and subframes:
            res._frame._superframes.update(self._frame._superframes)
            for sframe in self._frame._superframes:
                sframe._subframes.add(res._frame)
                sframe._restrictions[subset] = res._frame
            # The subchart frame is not a "top frame" in the supersets
            # (including self._domain):
            for dom in self._domain._supersets:
                if res._frame in dom._top_frames:
                    # it was added by the Chart constructor above
                    dom._top_frames.remove(res._frame)
            # Update of domain restrictions:
            self._dom_restrict[subset] = res
        return self._dom_restrict[subset]

    def valid_coordinates(self, *coordinates, **kwds):
        r"""
        Check whether a tuple of coordinates can be the coordinates of a
        point in the chart domain.

        INPUT:

        - ``*coordinates`` -- coordinate values
        - ``**kwds`` -- options:

          - ``tolerance=0``, to set the absolute tolerance in the test of
            coordinate ranges
          - ``parameters=None``, to set some numerical values to parameters


        OUTPUT:

        - True if the coordinate values are admissible in the chart range.

        """
        n = len(coordinates)
        if n != self._manifold._dim:
            return False
        if 'tolerance' in kwds:
            tolerance = kwds['tolerance']
        else:
            tolerance = 0
        if 'parameters' in kwds:
            parameters = kwds['parameters']
        else:
            parameters = None
        # Check of the coordinate ranges:
        for x, bounds in zip(coordinates, self._bounds):
            xmin = bounds[0][0] - tolerance
            min_included = bounds[0][1]
            xmax = bounds[1][0] + tolerance
            max_included = bounds[1][1]
            if parameters:
                xmin = xmin.subs(parameters)
                xmax = xmax.subs(parameters)
            if min_included:
                if x < xmin:
                    return False
            else:
                if x <= xmin:
                    return False
            if max_included:
                if x > xmax:
                    return False
            else:
                if x >= xmax:
                    return False
        # Check of additional restrictions:
        if self._restrictions != []:
            substitutions = dict([(self._xx[j], coordinates[j]) for j in
                                                                    range(n)])
            if parameters:
                substitutions.update(parameters)
            for restrict in self._restrictions:
                if isinstance(restrict, tuple): # case of or conditions
                    combine = False
                    for expr in restrict:
                        combine = combine or bool(expr.subs(substitutions))
                    if not combine:
                        return False
                else:
                    if not bool(restrict.subs(substitutions)):
                        return False
        # All tests have been passed:
        return True

    def transition_map(self, other, transformations, intersection_name=None,
                       restrictions1=None, restrictions2=None):
        r"""
        Construct the transition map between the current chart,
        `(U,\varphi)` say, and another one, `(V,\psi)` say.

        If `n` is the manifold's dimension, the *transition map* is the
        map

        .. MATH::

            \psi\circ\varphi^{-1}: \varphi(U\cap V) \subset \RR^n
            \rightarrow \psi(U\cap V) \subset \RR^n

        In other words, the
        transition map expresses the coordinates `(y^1,\ldots,y^n)` of
        `(V,\psi)` in terms of the coordinates `(x^1,\ldots,x^n)` of
        `(U,\varphi)` on the open subset where the two charts intersect, i.e.
        on `U\cap V`.

        INPUT:

        - ``other`` -- the chart `(V,\psi)`
        - ``transformations`` -- tuple (Y_1,...,Y_2), where Y_i is a symbolic
          expression expressing the coordinate `y^i` in terms of the
          coordinates `(x^1,\ldots,x^n)`
        - ``intersection_name`` -- (default: ``None``) name to be given to the
          subset `U\cap V` if the latter differs from `U` or `V`
        - ``restrictions1`` -- (default: ``None``) list of conditions on the
          coordinates of the current chart that define `U\cap V` if the
          latter differs from `U`. ``restrictions1`` must be a list of
          of symbolic equalities or inequalities involving the
          coordinates, such as x>y or x^2+y^2 != 0. The items of the list
          ``restrictions1`` are combined with the ``and`` operator; if some
          restrictions are to be combined with the ``or`` operator instead,
          they have to be passed as a tuple in some single item of the list
          ``restrictions1``. For example, ``restrictions1`` = [x>y,
          (x!=0, y!=0), z^2<x] means (x>y) and ((x!=0) or (y!=0)) and (z^2<x).
          If the list ``restrictions1`` contains only one item, this item can
          be passed as such, i.e. writing x>y instead of the single element
          list [x>y].
        - ``restrictions2`` -- (default: ``None``) list of conditions on the
          coordinates of the other chart that define `U\cap V` if the latter
          differs from `V` (see ``restrictions1`` for the syntax)

        OUTPUT:

        - The transition map `\psi\circ\varphi^{-1}` defined on `U\cap V`, as an
          instance of :class:`CoordChange`.

        EXAMPLES:

        Transition map between two stereographic charts on the circle `S^1`::

            sage: M = Manifold(1, 'S^1')
            sage: U = M.open_subset('U') # Complement of the North pole
            sage: cU.<x> = U.chart() # Stereographic chart from the North pole
            sage: V = M.open_subset('V') # Complement of the South pole
            sage: cV.<y> = V.chart() # Stereographic chart from the South pole
            sage: M.declare_union(U,V)   # S^1 is the union of U and V
            sage: trans = cU.transition_map(cV, 1/x, 'W', x!=0, y!=0)
            sage: trans
            coordinate change from chart (W, (x,)) to chart (W, (y,))
            sage: M.list_of_subsets() # the subset W, intersection of U and V, has been created by transition_map()
            [1-dimensional manifold 'S^1',
             open subset 'U' of the 1-dimensional manifold 'S^1',
             open subset 'V' of the 1-dimensional manifold 'S^1',
             open subset 'W' of the 1-dimensional manifold 'S^1']
            sage: W = M.list_of_subsets()[3]
            sage: W is U.intersection(V)
            True
            sage: M.atlas()
            [chart (U, (x,)), chart (V, (y,)), chart (W, (x,)), chart (W, (y,))]

        Transition map between spherical chart and Cartesian chart on `\RR^2`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart()
            sage: U = M.open_subset('U') # the complement of the half line {y=0, x >= 0}
            sage: c_spher.<r,phi> = U.chart(r'r:(0,+oo) phi:(0,2*pi):\phi')
            sage: trans = c_spher.transition_map(c_cart, (r*cos(phi), r*sin(phi)), \
                                                 restrictions2=(y!=0, x<0))
            sage: trans
            coordinate change from chart (U, (r, phi)) to chart (U, (x, y))
            sage: M.list_of_subsets() # in this case, no new subset has been created since U inter M = U
            [2-dimensional manifold 'R^2',
            open subset 'U' of the 2-dimensional manifold 'R^2']
            sage: M.atlas() # ...but a new chart has been created: (U, (x, y))
            [chart (R^2, (x, y)), chart (U, (r, phi)), chart (U, (x, y))]

        """
        dom1 = self._domain
        dom2 = other._domain
        dom = dom1.intersection(dom2, name=intersection_name)
        if dom is dom1:
            chart1 = self
        else:
            chart1 = self.restrict(dom, restrictions1)
        if dom is dom2:
            chart2 = other
        else:
            chart2 = other.restrict(dom, restrictions2)
        if not isinstance(transformations, (tuple, list)):
                transformations = [transformations]
        return CoordChange(chart1, chart2, *transformations)


    def coord_change(self, other, *transformations):
        r"""
        Relate the coordinates of ``self`` to those of another chart.

        The two charts may belong to different manifolds.

        .. NOTE::

            For defining the transition map between two charts on overlapping
            domains, use the method :meth:`transition_map` instead.

        See class :class:`CoordChange` for a complete documentation.

        INPUT:

        - ``other`` -- another chart
        - ``transformations`` -- the coordinate transformations expressed as
          a list of the expressions of the coordinates of ``other`` in terms
          of the coordinates of ``self``

        OUTPUT:

        - instance of class :class:`CoordChange` representing the relation
          between the two sets of coordinates.

        EXAMPLE:

        Connecting spherical (polar) coordinates to Cartesian ones in the
        plane::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x, y> = M.chart() # Cartesian coordinates on the plane
            sage: U = M.open_subset('U') # the complement of the half line {y=0, x>= 0}
            sage: c_spher.<r, ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: spher_to_cart = c_spher.coord_change(c_cart, r*cos(ph), r*sin(ph))
            sage: spher_to_cart
            coordinate change from chart (U, (r, ph)) to chart (R^2, (x, y))
            sage: type(spher_to_cart)
            <class 'sage.geometry.manifolds.chart.CoordChange'>
            sage: spher_to_cart(1, pi/2)
            (0, 1)

        """
        return CoordChange(self, other, *transformations)

    def function(self, expression):
        r"""
        Return a real-valued function of the coordinates.

        If ``self`` is a chart on a `n`-dimensional manifold, a real-valued
        function of the coordinates is a is a function

        .. MATH::

            \begin{array}{llcl}
            f:& U \subset\RR^n & \longrightarrow & \RR \\
              & (x^1,\ldots,x^n) & \longmapsto & f(x^1,\ldots,x^n)
            \end{array}

        where `U` is the subset of `\RR^n` covered by the chart ``self``.

        See class :class:`FunctionChart` for a complete documentation.

        INPUT:

        - ``expression`` -- the coordinate expression of the function

        OUTPUT:

        - instance of class :class:`FunctionChart` representing the function
          `f`

        EXAMPLE:

        Function of two coordinates::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.function(sin(x*y))
            sage: type(f)
            <class 'sage.geometry.manifolds.chart.FunctionChart'>
            sage: f.display()
            (x, y) |--> sin(x*y)
            sage: f(2,3)
            sin(6)

        """
        return FunctionChart(self, expression)

    def zero_function(self):
        r"""
        Return the zero function defined on the chart.

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.chart.ZeroFunctionChart`
          representing the zero function of the chart coordinates.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: X.zero_function()
            0
            sage: X.zero_function().display()
            (x, y) |--> 0
            sage: type(X.zero_function())
            <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>

        The result is cached::

            sage: X.zero_function() is X.zero_function()
            True

        """
        return self._zero_function

    def multifunction(self, *expressions):
        r"""
        Return a `\RR^m`-valued function of the coordinates.

        If ``self`` is a chart on a `n`-dimensional manifold and `m` is
        an integer strictly greater than 0, the returned function is of the
        type

        .. MATH::

            \begin{array}{llcl}
            f:& U \subset\RR^n & \longrightarrow & \RR^m \\
              & (x^1,\ldots,x^n) & \longmapsto & (f_1(x^1,\ldots,x^n),\ldots,
                f_m(x^1,\ldots,x^n))
            \end{array}

        where `U` is the subset of `\RR^n` covered by the chart ``self``.

        See class :class:`MultiFunctionChart` for a complete documentation.

        INPUT:

        - ``*expressions`` -- the list of the coordinate expressions of the `m`
          functions (`m\geq 1`)

        OUTPUT:

        - instance of class :class:`MultiFunctionChart` representing the
          function `f`

        EXAMPLE:

        Function of two coordinates with values in `\RR^3`::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.multifunction(x+cos(y), sin(y)/x^2, x*y) ; f
            functions (x + cos(y), sin(y)/x^2, x*y) on the chart (M, (x, y))
            sage: type(f)
            <class 'sage.geometry.manifolds.chart.MultiFunctionChart'>
            sage: f(2,3,4)
            (cos(3) + 2, 1/4*sin(3), 6)

        """
        return MultiFunctionChart(self, *expressions)

    def plot(self, chart=None, ambient_coords=None, mapping=None,
             fixed_coords=None, ranges=None, max_value=8, nb_values=None,
             steps=None, parameters=None, color='red',  style='-', thickness=1,
             plot_points=75, label_axes=True):
        r"""
        Plot the current chart (``self``) as a "grid" in a Cartesian graph
        based on the coordinates of some ambient chart.

        The "grid" is formed by curves along which a coordinate of ``self``
        varies, the other coordinates being kept fixed; it is drawn in terms of
        two (2D graphics) or three (3D graphics) coordinates of another chart,
        called hereafter the *ambient chart*.

        The ambient chart is related to the current chart (``self``) either by
        a transition map if both charts are defined on the same manifold, or by
        the coordinate expression of some differentiable mapping (typically an
        immersion). In the latter case, the two charts may be defined on two
        different manifolds.

        INPUT:

        - ``chart`` -- (default: ``None``) the ambient chart (see above); if
          ``None``, the ambient chart is set to ``self``
        - ``ambient_coords`` -- (default: ``None``) tuple containing the 2 or 3
          coordinates of the ambient chart in terms of which the plot is
          performed; if ``None``, all the coordinates of the ambient chart are
          considered
        - ``mapping`` -- (default: ``None``) differentiable mapping (instance
          of :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`)
          providing the link between ``self`` and the ambient chart (cf.
          above); if ``None``, both charts are supposed to be defined on the same
          manifold and related by some transition map (see
          :meth:`transition_map`)
        - ``fixed_coords`` -- (default: ``None``) dictionary with keys the
          coordinates of ``self`` that are not drawn and with values the fixed
          value of these coordinates; if ``None``, all the coordinates of ``self``
          are drawn
        - ``ranges`` -- (default: ``None``) dictionary with keys the coordinates
          to be drawn and values tuples ``(x_min,x_max)`` specifying the
          coordinate range for the plot; if ``None``, the entire coordinate range
          declared during the chart construction is considered (with -Infinity
          replaced by ``-max_value`` and +Infinity by ``max_value``)
        - ``max_value`` -- (default: 8) numerical value substituted to
          +Infinity if the latter is the upper bound of the range of a
          coordinate for which the plot is performed over the entire coordinate
          range (i.e. for which no specific plot range has been set in
          ``ranges``); similarly ``-max_value`` is the numerical valued
          substituted for -Infinity
        - ``nb_values`` -- (default: ``None``) either an integer or a dictionary
          with keys the coordinates to be drawn and values the number of
          constant values of the coordinate to be considered; if ``nb_values``
          is a single integer, it represents the number of constant values for all
          coordinates; if ``nb_values`` is ``None``, it is set to 9 for a 2D
          plot and to 5 for a 3D plot
        - ``steps`` -- (default: ``None``) dictionary with keys the coordinates
          to be drawn and values the step between each constant value of
          the coordinate; if ``None``, the step is computed from the coordinate
          range (specified in ``ranges``) and ``nb_values``. On the contrary
          if the step is provided for some coordinate, the corresponding
          number of constant values is deduced from it and the coordinate range.
        - ``parameters`` -- (default: ``None``) dictionary giving the numerical
          values of the parameters that may appear in the relation between
          the two coordinate systems
        - ``color`` -- (default: 'red') either a single color or a dictionary
          of colors, with keys the coordinates to be drawn, representing the
          colors of the lines along which the coordinate varies, the other
          being kept constant; if ``color`` is a single color, it is used for
          all coordinate lines
        - ``style`` -- (default: '-') either a single line style or a dictionary
          of line styles, with keys the coordinates to be drawn, representing
          the style of the lines along which the coordinate varies, the other
          being kept constant; if ``style`` is a single style, it is used for
          all coordinate lines; NB: ``style`` is effective only for 2D plots
        - ``thickness`` -- (default: 1) either a single line thickness or a
          dictionary of line thicknesses, with keys the coordinates to be drawn,
          representing the thickness of the lines along which the coordinate
          varies, the other being kept constant; if ``thickness`` is a single
          value, it is used for all coordinate lines
        - ``plot_points`` -- (default: 75) either a single number of points or
          a dictionary of integers, with keys the coordinates to be drawn,
          representing the number of points to plot the lines along which the
          coordinate varies, the other being kept constant; if ``plot_points``
          is a single integer, it is used for all coordinate lines
        - ``label_axes`` -- (default: ``True``) boolean determining whether the
          labels of the ambient coordinate axes shall be added to the graph;
          can be set to False if the graph is 3D and must be superposed with
          another graph.

        OUTPUT:

        - a graphic object, either an instance of
          :class:`~sage.plot.graphics.Graphics` for a 2D plot (i.e. based on
          2 coordinates of the ambient chart) or an instance of
          :class:`~sage.plot.plot3d.base.Graphics3d` for a 3D plot (i.e.
          based on 3 coordinates of the ambient chart)

        EXAMPLES:

        Grid of polar coordinates in terms of Cartesian coordinates in the
        Euclidean plane::

            sage: R2 = Manifold(2, 'R^2') # the Euclidean plane
            sage: c_cart.<x,y> = R2.chart() # Cartesian coordinates
            sage: U = R2.open_subset('U', coord_def={c_cart: (y!=0, x<0)}) # the complement of the segment y=0 and x>0
            sage: c_pol.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # polar coordinates on U
            sage: pol_to_cart = c_pol.transition_map(c_cart, [r*cos(ph), r*sin(ph)])
            sage: g = c_pol.plot(c_cart)
            sage: type(g)
            <class 'sage.plot.graphics.Graphics'>
            sage: show(g) # graphical display

        Call with non-default values::

            sage: g = c_pol.plot(c_cart, ranges={ph:(pi/4,pi)}, nb_values={r:7, ph:17}, \
            ....:                color={r:'red', ph:'green'}, style={r:'-', ph:'--'})

        A single coordinate line can be drawn::

            sage: g = c_pol.plot(c_cart, fixed_coords={r: 2}) # draw a circle of radius r=2
            sage: g = c_pol.plot(c_cart, fixed_coords={ph: pi/4}) # draw a segment at phi=pi/4

        A chart can be plot in terms of itself, resulting in a rectangular grid::

            sage: g = c_cart.plot()  # equivalent to c_cart.plot(c_cart)
            sage: show(g) # a rectangular grid

        An example with the ambient chart given by the coordinate expression of
        some differentiable mapping: 3D plot of the stereographic charts on the
        2-sphere::

            sage: S2 = Manifold(2, 'S^2') # the 2-sphere
            sage: U = S2.open_subset('U') ; V = S2.open_subset('V') # complement of the North and South pole, respectively
            sage: S2.declare_union(U,V)
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                             intersection_name='W', restrictions1= x^2+y^2!=0, \
                                             restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: R3 = Manifold(3, 'R^3') # the Euclidean space R^3
            sage: c_cart.<X,Y,Z> = R3.chart()  # Cartesian coordinates on R^3
            sage: Phi = S2.diff_mapping(R3, {(c_xy, c_cart): [2*x/(1+x^2+y^2), \
                                              2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)], \
                                             (c_uv, c_cart): [2*u/(1+u^2+v^2), \
                                              2*v/(1+u^2+v^2), (1-u^2-v^2)/(1+u^2+v^2)]}, \
                                        name='Phi', latex_name=r'\Phi') # Embedding of S^2 in R^3
            sage: g = c_xy.plot(c_cart, mapping=Phi)
            sage: show(g) # 3D graphic display
            sage: type(g)
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>

        The same plot without the (X,Y,Z) axes labels::

            sage: g = c_xy.plot(c_cart, mapping=Phi, label_axes=False)

        The North and South stereographic charts on the same plot::

            sage: g2 = c_uv.plot(c_cart, mapping=Phi, color='green')
            sage: show(g+g2)

        South stereographic chart drawned in terms of the North one (we split
        the plot in four parts to avoid the singularity at (u,v)=(0,0))::

            sage: W = U.intersection(V) # the subset common to both charts
            sage: c_uvW = c_uv.restrict(W) # chart (W,(u,v))
            sage: gSN1 = c_uvW.plot(c_xy, ranges={u:[-6.,-0.02], v:[-6.,-0.02]}, nb_values=20, plot_points=100)
            sage: gSN2 = c_uvW.plot(c_xy, ranges={u:[-6.,-0.02], v:[0.02,6.]}, nb_values=20, plot_points=100)
            sage: gSN3 = c_uvW.plot(c_xy, ranges={u:[0.02,6.], v:[-6.,-0.02]}, nb_values=20, plot_points=100)
            sage: gSN4 = c_uvW.plot(c_xy, ranges={u:[0.02,6.], v:[0.02,6.]}, nb_values=20, plot_points=100)
            sage: show(gSN1+gSN2+gSN3+gSN4, xmin=-3, xmax=3, ymin=-3, ymax=3)

        The coordinate line u=1 (red) and the coordinate line v=1 (green) on
        the same plot::

            sage: gu1 = c_uvW.plot(c_xy, fixed_coords={u:1}, max_value=20, plot_points=200)
            sage: gv1 = c_uvW.plot(c_xy, fixed_coords={v:1}, max_value=20, plot_points=200, color='green')
            sage: show(gu1+gv1)

        Note that we have set ``max_value=20`` to have a wider range for the
        coordinates u and v, i.e. to have [-20,20] instead of the default
        [-8,8].

        A 3-dimensional chart plotted in terms of itself results in a 3D
        rectangular grid::

            sage: g = c_cart.plot() # equivalent to c_cart.plot(c_cart)
            sage: show(g)  # a 3D mesh cube

        A 4-dimensional chart plotted in terms of itself (the plot is
        performed for at most 3 coordinates, which must be specified via
        the argument ``ambient_coords``)::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: g = X.plot(ambient_coords=(t,x,y))  # the coordinate z is not depicted
            sage: show(g)  # a 3D mesh cube
            sage: g = X.plot(ambient_coords=(t,y)) # the coordinates x and z are not depicted
            sage: show(g)  # a 2D mesh square

        """
        from sage.misc.functional import numerical_approx
        from sage.plot.graphics import Graphics
        from sage.plot.line import line
        from diffmapping import DiffMapping
        from utilities import set_axes_labels
        if chart is None:
            chart = self
        elif not isinstance(chart, Chart):
            raise TypeError("the argument 'chart' must be a coordinate chart")
        #
        # 1/ Determination of the relation between self and chart
        #    ------------------------------------------------------------
        nc = self._manifold._dim
        if chart is self:
            transf = self.multifunction(*(self._xx))
            if nc > 3:
                if ambient_coords is None:
                    raise TypeError("the argument 'ambient_coords' must be " +
                                    "provided")
                if len(ambient_coords) > 3:
                    raise ValueError("too many ambient coordinates")
                fixed_coords = {}
                for coord in self._xx:
                    if coord not in ambient_coords:
                        fixed_coords[coord] = 0
        else:
            transf = None # to be the MultiFunctionChart relating self to
                          # the ambient chart
            if mapping is None:
                if not self._domain.is_subset(chart._domain):
                    raise TypeError("the domain of {} is not ".format(self) +
                                    "included in that of {}".format(chart))
                coord_changes = chart._domain._coord_changes
                for chart_pair in coord_changes:
                    if chart_pair == (self, chart):
                        transf = coord_changes[chart_pair]._transf
                        break
                else:
                    # Search for a subchart
                    for chart_pair in coord_changes:
                        for schart in chart._subcharts:
                            if chart_pair == (self, schart):
                                transf = coord_changes[chart_pair]._transf
            else:
                if not isinstance(mapping, DiffMapping):
                    raise TypeError("the argument 'mapping' must be a " +
                                    "differentiable mapping")
                if not self._domain.is_subset(mapping._domain):
                    raise TypeError("the domain of {} is not ".format(self) +
                                    "included in that of {}".format(mapping))
                if not chart._domain.is_subset(mapping._codomain):
                    raise TypeError("the domain of {} is not ".format(chart) +
                                    "included in the codomain of {}".format(
                                                                      mapping))
                try:
                    transf = mapping.multi_function_chart(chart1=self,
                                                          chart2=chart)
                except ValueError:
                    pass
            if transf is None:
                raise ValueError("no relation has been found between " +
                                 "{} and {}".format(self, chart))
        #
        # 2/ Treatment of input parameters
        #    -----------------------------
        if fixed_coords is None:
            coords = self._xx
        else:
            fixed_coord_list = fixed_coords.keys()
            coords = []
            for coord in self._xx:
                if coord not in fixed_coord_list:
                    coords.append(coord)
            coords = tuple(coords)
        if ambient_coords is None:
            ambient_coords = chart._xx
        elif not isinstance(ambient_coords, tuple):
            ambient_coords = tuple(ambient_coords)
        nca = len(ambient_coords)
        if nca != 2 and nca !=3:
            raise TypeError("bad number of ambient coordinates: {}".format(nca))
        if ranges is None:
            ranges = {}
        ranges0 = {}
        for coord in coords:
            if coord in ranges:
                ranges0[coord] = (numerical_approx(ranges[coord][0]),
                                  numerical_approx(ranges[coord][1]))
            else:
                bounds = self._bounds[self._xx.index(coord)]
                if bounds[0][0] == -Infinity:
                    xmin = numerical_approx(-max_value)
                elif bounds[0][1]:
                    xmin = numerical_approx(bounds[0][0])
                else:
                    xmin = numerical_approx(bounds[0][0] + 1.e-3)
                if bounds[1][0] == Infinity:
                    xmax = numerical_approx(max_value)
                elif bounds[1][1]:
                    xmax = numerical_approx(bounds[1][0])
                else:
                    xmax = numerical_approx(bounds[1][0] - 1.e-3)
                ranges0[coord] = (xmin, xmax)
        ranges = ranges0
        if nb_values is None:
            if nca == 2: # 2D plot
                nb_values = 9
            else:   # 3D plot
                nb_values = 5
        if not isinstance(nb_values, dict):
            nb_values0 = {}
            for coord in coords:
                nb_values0[coord] = nb_values
            nb_values = nb_values0
        if steps is None:
            steps = {}
        for coord in coords:
            if coord not in steps:
                steps[coord] = (ranges[coord][1] - ranges[coord][0])/ \
                               (nb_values[coord]-1)
            else:
                nb_values[coord] = 1 + int(
                           (ranges[coord][1] - ranges[coord][0])/ steps[coord])
        if not isinstance(color, dict):
            color0 = {}
            for coord in coords:
                color0[coord] = color
            color = color0
        if not isinstance(style, dict):
            style0 = {}
            for coord in coords:
                style0[coord] = style
            style = style0
        if not isinstance(thickness, dict):
            thickness0 = {}
            for coord in coords:
                thickness0[coord] = thickness
            thickness = thickness0
        if not isinstance(plot_points, dict):
            plot_points0 = {}
            for coord in coords:
                plot_points0[coord] = plot_points
            plot_points = plot_points0
        #
        # 3/ Plots
        #    -----
        xx0 = [0] * nc
        if fixed_coords is not None:
            if len(fixed_coords) != nc - len(coords):
                raise ValueError("Bad number of fixed coordinates.")
            for fc, val in fixed_coords.iteritems():
                xx0[self._xx.index(fc)] = val
        ind_a = [chart._xx.index(ac) for ac in ambient_coords]
        resu = Graphics()
        for coord in coords:
            color_c, style_c = color[coord], style[coord]
            thickness_c = thickness[coord]
            rem_coords = list(coords)
            rem_coords.remove(coord)
            xx_list = [xx0]
            if len(rem_coords) >= 1:
                xx_list = self._plot_xx_list(xx_list, rem_coords, ranges,
                                             steps, nb_values)
            xmin, xmax = ranges[coord]
            nbp = plot_points[coord]
            dx = (xmax - xmin) / (nbp-1)
            ind_coord = self._xx.index(coord)
            for xx in xx_list:
                curve = []
                first_invalid = False # initialization
                xc = xmin
                xp = list(xx)
                if parameters is None:
                    for i in range(nbp):
                        xp[ind_coord] = xc
                        if self.valid_coordinates(*xp, tolerance=1e-13):
                            yp = transf(*xp, simplify=False)
                            curve.append( [numerical_approx(yp[j])
                                           for j in ind_a] )
                            first_invalid = True # next invalid point will be
                                                 # the first one
                        else:
                            if first_invalid:
                                # the curve is stopped at previous point and
                                # added to the graph:
                                resu += line(curve, color=color_c,
                                             linestyle=style_c,
                                             thickness=thickness_c)
                                curve = [] # a new curve will start at the
                                           # next valid point
                            first_invalid = False # next invalid point will not
                                                  # be the first one
                        xc += dx
                else:
                    for i in range(nbp):
                        xp[ind_coord] = xc
                        if self.valid_coordinates(*xp, tolerance=1e-13,
                                                  parameters=parameters):
                            yp = transf(*xp, simplify=False)
                            curve.append(
                              [numerical_approx( yp[j].substitute(parameters) )
                               for j in ind_a] )
                            first_invalid = True # next invalid point will be
                                                 # the first one
                        else:
                            if first_invalid:
                                # the curve is stopped at previous point and
                                # added to the graph:
                                resu += line(curve, color=color_c,
                                             linestyle=style_c,
                                             thickness=thickness_c)
                                curve = [] # a new curve will start at the
                                           # next valid point
                            first_invalid = False # next invalid point will not
                                                  # be the first one
                        xc += dx
                if curve != []:
                    resu += line(curve, color=color_c,
                                 linestyle=style_c,
                                 thickness=thickness_c)
        if nca==2:  # 2D graphic
            resu.set_aspect_ratio(1)
            if label_axes:
                # We update the dictionary _extra_kwds (options to be passed
                # to show()), instead of using the method
                # Graphics.axes_labels() since the latter is not robust w.r.t.
                # graph addition
                resu._extra_kwds['axes_labels'] = [r'$'+latex(ac)+r'$'
                                                   for ac in ambient_coords]
        else: # 3D graphic
            resu.aspect_ratio(1)
            if label_axes:
                labels = [str(ac) for ac in ambient_coords]
                resu = set_axes_labels(resu, *labels)
        return resu

    def _plot_xx_list(self, xx_list, rem_coords, ranges, steps, nb_values):
        coord = rem_coords[0]
        xmin = ranges[coord][0]
        sx = steps[coord]
        resu = []
        for xx in xx_list:
            xc = xmin
            for i in range(nb_values[coord]):
                nxx = list(xx)
                nxx[self._xx.index(coord)] = xc
                resu.append(nxx)
                xc += sx
        if len(rem_coords) == 1:
            return resu
        else:
            rem_coords.remove(coord)
            return self._plot_xx_list(resu, rem_coords, ranges, steps,
                                      nb_values)

#*****************************************************************************

class FunctionChart(SageObject):
    r"""
    Real-valued function on a chart codomain.

    Given a chart `(U, \varphi)` on a differentiable manifold `M` of
    dimension `n`, i.e. an open subset `U` of `M` along with an
    homeomorphism `\varphi: U \rightarrow V \subset \RR^n` (`V` being
    an open subset of `\RR^n`), an instance of class :class:`FunctionChart` is
    a function

    .. MATH::

        \begin{array}{llcl}
        f:& V \subset\RR^n & \longrightarrow & \RR \\
          & (x^1,\ldots,x^n) & \longmapsto & f(x^1,\ldots,x^n)
        \end{array}

    Instances of :class:`FunctionChart` differ from Sage symbolic functions,
    implemented as callable symbolic expressions, by the automatic
    simplification of the result of any operation, via the function
    :func:`~sage.geometry.manifolds.utilities.simplify_chain`.

    INPUT:

    - ``chart`` -- the chart `(U, \varphi)`, as an instance of
      class :class:`~sage.geometry.manifolds.chart.Chart`
    - ``expression`` -- a symbolic expression representing `f(x^1,\ldots,x^n)`,
      where `(x^1,\ldots,x^n)` are the coordinates of the chart `(U, \varphi)`

    EXAMPLES:

    Function defined on a 2-dimensional chart::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: f = c_xy.function(x^2+3*y+1)
        sage: type(f)
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        sage: f.display()
        (x, y) |--> x^2 + 3*y + 1
        sage: f(x,y)
        x^2 + 3*y + 1

    The symbolic expression is returned when asking the direct display of
    the function::

        sage: f
        x^2 + 3*y + 1
        sage: latex(f)
        x^{2} + 3 \, y + 1

    or via the method :meth:`expr`::

        sage: f.expr()
        x^2 + 3*y + 1

    The value of the function at specified coordinates is obtained by means
    of the standard parentheses notation::

        sage: f(2,-1)
        2
        sage: var('a b')
        (a, b)
        sage: f(a,b)
        a^2 + 3*b + 1

    An unspecified function on a chart::

        sage: g = c_xy.function(function('G', x, y))
        sage: g
        G(x, y)
        sage: g.display()
        (x, y) |--> G(x, y)
        sage: g.expr()
        G(x, y)
        sage: g(2,3)
        G(2, 3)

    Chart functions can be compared to other values::

        sage: f = c_xy.function(x^2+3*y+1)
        sage: f == 2
        False
        sage: f == x^2 + 3*y + 1
        True
        sage: g = c_xy.function(x*y)
        sage: f == g
        False
        sage: h = c_xy.function(x^2+3*y+1)
        sage: f == h
        True

    .. RUBRIC:: Differences between ``FunctionChart`` and callable symbolic
      expressions

    Callable symbolic expressions are defined directly from symbolic
    expressions of the coordinates::

        sage: f0(x,y) = x^2 + 3*y + 1
        sage: type(f0)
        <type 'sage.symbolic.expression.Expression'>
        sage: f0
        (x, y) |--> x^2 + 3*y + 1
        sage: f0(x,y)
        x^2 + 3*y + 1

    To get an output similar to that of ``f0`` for the chart function ``f``,
    we must use the method :meth:`display`::

        sage: f
        x^2 + 3*y + 1
        sage: f.display()
        (x, y) |--> x^2 + 3*y + 1
        sage: f(x,y)
        x^2 + 3*y + 1

    More importantly, instances of :class:`FunctionChart` differ from
    callable symbolic expression by the automatic simplifications in all
    operations. For instance, adding the two callable symbolic expressions::

        sage: f0(x,y,z) = cos(x)^2 ; g0(x,y,z) = sin(x)^2

    results in::

        sage: f0 + g0
        (x, y, z) |--> cos(x)^2 + sin(x)^2

    To get 1,  one has to call
    :meth:`~sage.symbolic.expression.Expression.simplify_trig`::

        sage: (f0 + g0).simplify_trig()
        (x, y, z) |--> 1

    On the contrary, the sum of the corresponding :class:`FunctionChart`
    instances is automatically simplified (see
    :func:`~sage.geometry.manifolds.utilities.simplify_chain` for details)::

        sage: f = c_xy.function(cos(x)^2) ; g = c_xy.function(sin(x)^2)
        sage: f + g
        1

    Another difference regards the display of partial derivatives: for callable
    symbolic functions, it relies on Pynac notation ``D[0]``, ``D[1]``, etc.::

        sage: g = function('g', x, y)
        sage: f0(x,y) = diff(g, x) + diff(g, y)
        sage: f0
        (x, y) |--> D[0](g)(x, y) + D[1](g)(x, y)

    while for chart functions, the display is more "textbook" like::

        sage: f = c_xy.function(diff(g, x) + diff(g, y))
        sage: f
        d(g)/dx + d(g)/dy

    The difference is even more dramatic on LaTeX outputs::

        sage: latex(f0)
        \left( x, y \right) \ {\mapsto} \ D[0]\left(g\right)\left(x, y\right) + D[1]\left(g\right)\left(x, y\right)
        sage: latex(f)
        \frac{\partial\,g}{\partial x} + \frac{\partial\,g}{\partial y}

    One can switch to Pynac notation via the command
    :func:`~sage.geometry.manifolds.utilities.nice_derivatives`::

        sage: nice_derivatives(False)
        sage: latex(f)
        D[0]\left(g\right)\left(x, y\right) + D[1]\left(g\right)\left(x, y\right)
        sage: nice_derivatives(True)
        sage: latex(f)
        \frac{\partial\,g}{\partial x} + \frac{\partial\,g}{\partial y}

    """

    nice_output = True # static flag for textbook-like output instead of the
                       # Pynac output for derivatives

    omit_fargs  = False # static flag to govern whether or not
                        # the arguments of symbolic functions are printed



    def __init__(self, chart, expression):
        from sage.symbolic.ring import SR
        self._chart = chart
        self._express = SR(expression)
        self._nc = len(self._chart._xx)    # number of coordinates
        # Derived quantities:
        self._der = None  # partial derivatives

    def _repr_(self):
        r"""
        String representation of the object.
        """
        if FunctionChart.nice_output:
            return str(ExpressionNice(self._express))
        else:
            return str(self._express)

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        if FunctionChart.nice_output:
            return latex(ExpressionNice(self._express))
        else:
            return latex(self._express)

    def expr(self):
        r"""
        Return the expression of the image of the function.

        OUTPUT:

        - symbolic expression, involving the chart coordinates.

        EXAMPLES:

        Function on some chart of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.function(x^2+3*y+1)
            sage: f
            x^2 + 3*y + 1
            sage: f.expr()
            x^2 + 3*y + 1
            sage: type(f.expr())
            <type 'sage.symbolic.expression.Expression'>
            sage: f.expr() is f._express
            True

        The method :meth:`expr` is useful for accessing to all the
        symbolic expression functionalities in Sage; for instance::

            sage: var('a')
            a
            sage: f = c_xy.function(a*x*y)
            sage: f.expr()
            a*x*y
            sage: f.expr().subs(a=2)
            2*x*y

        Note that for substituting the value of a coordinate, the function call
        can be used as well::

            sage: f(x,3)
            3*a*x
            sage: bool( f(x,3) == f.expr().subs(y=3) )
            True

        """
        return self._express

    def display(self):
        r"""
        Display the function in arrow notation.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        EXAMPLE:

        Function on a 2-dimensional chart::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.function(x^2+3*y+1)
            sage: f.display()
            (x, y) |--> x^2 + 3*y + 1
            sage: latex(f.display())
             \left(x, y\right) \mapsto x^{2} + 3 \, y + 1

        A shortcut is ``disp()``::

            sage: f.disp()
            (x, y) |--> x^2 + 3*y + 1

        """
        from sage.tensor.modules.format_utilities import FormattedExpansion
        resu_txt = repr((self._chart)[:]) + ' |--> ' + repr(self._express)
        resu_latex = latex(self._chart[:]) + r' \mapsto' + \
                     latex(self._express)
        return FormattedExpansion(resu_txt, resu_latex)

    disp = display

    def view(self):
        r"""
        Deprecated method.

        Use method :meth:`display` instead.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.function(x^2+3*y+1)
            sage: f.view()
            doctest:...: DeprecationWarning: Use function display() instead.
            See http://trac.sagemath.org/15916 for details.
            (x, y) |--> x^2 + 3*y + 1
            sage: f.display()
            (x, y) |--> x^2 + 3*y + 1

        """
        from sage.misc.superseded import deprecation
        deprecation(15916, 'Use function display() instead.')
        return self.display()


    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        self._der = None

    def copy(self):
        r"""
        Returns an exact copy of ``self``.

        The derived quantities are not copied, because they can be
        reconstructed if necessary.

        EXAMPLES:

        Copy on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.function(x^2+3*y+1)
            sage: g = f.copy()
            sage: type(g)
            <class 'sage.geometry.manifolds.chart.FunctionChart'>
            sage: g
            x^2 + 3*y + 1
            sage: g == f    # g is mathematically equal to f:
            True
            sage: g is f    # but differs in computer memory:
            False

        """
        return FunctionChart(self._chart, self._express)

    def __call__(self, *coords, **options):
        r"""
        Computes the value of the function at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates `(x^1,...,x^n)` where the
          function `f` is to be evaluated
        - ``**options`` -- allows to pass ``simplify=False`` to disable the
          call of the simplification chain on the result

        OUTPUT:

        - the value `f(x^1,...,x^n)`

        """
        #!# This should be the Python 2.7 form:
        # substitutions = {self._chart._xx[j]: coords[j] for j in
        #                                                      range(self._nc)}
        #
        # Here we use a form compatible with Python 2.6:
        substitutions = dict([(self._chart._xx[j], coords[j]) for j in
                                                              range(self._nc)])
        resu = self._express.subs(substitutions)
        if 'simplify' in options:
            if options['simplify']:
                return simplify_chain(resu)
            else:
                return resu
        else:
            return simplify_chain(resu)


    def diff(self, coord):
        r"""
        Partial derivative with respect to a coordinate.

        INPUT:

        - ``coord`` -- either the coordinate `x^i` with respect
          to which the derivative of the function `f` is to be taken, or the
          index `i` labelling this coordinate

        OUTPUT:

        - the partial derivative `\frac{\partial f}{\partial x^i}`, as an
          instance of :class:`FunctionChart`

        EXAMPLES:

        Partial derivatives of a function defined on a 2-dimensional chart::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.function(x^2+3*y+1) ; f
            x^2 + 3*y + 1
            sage: f.diff(x)
            2*x
            sage: f.diff(y)
            3

        The partial derivatives are instances of the class
        :class:`FunctionChart`::

            sage: type(f.diff(x))
            <class 'sage.geometry.manifolds.chart.FunctionChart'>

        An index can be used instead of the coordinate symbol::

            sage: f.diff(0)
            2*x
            sage: f.diff(0) is f.diff(x)
            True

        The index range depends on the convention used on the manifold::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.function(x^2+3*y+1)
            sage: f.diff(1)
            2*x
            sage: f.diff(1) is f.diff(x)
            True

        """
        from sage.calculus.functional import diff
        if self._der is None:
            # the partial derivatives have to be updated
            self._der = [FunctionChart(self._chart,
                         simplify_chain(diff(self._express, self._chart._xx[j])))
                                                    for j in range(self._nc) ]
        if isinstance(coord, (int, Integer)):
            return self._der[coord - self._chart._manifold._sindex]
        else:
            return self._der[self._chart._xx.index(coord)]

    def is_zero(self):
        r"""
        Return ``True`` if the function is zero and ``False`` otherwise.

        EXAMPLES:

        Functions on a 2-dimensional chart::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.function(x^2+3*y+1)
            sage: f.is_zero()
            False
            sage: g = c_xy.function(0)
            sage: g.is_zero()
            True

        """
        return self._express.is_zero()

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - True if ``self`` is equal to ``other``,  or False otherwise

        """
        if isinstance(other, FunctionChart):
            if other._chart != self._chart:
                return False
            else:
                return bool(other._express == self._express)
        else:
            return bool(self._express == other)

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - True if ``self`` is different from ``other``,  or False otherwise

        """
        return not self.__eq__(other)

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of ``self``

        """
        return FunctionChart(self._chart, self._express)

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the opposite of the function ``self``

        """
        return FunctionChart(self._chart, simplify_chain(-self._express))

    def __add__(self, other):
        r"""
        Addition operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - chart function resulting from the addition of ``self`` and ``other``

        """
        if isinstance(other, FunctionChart):
            if other._chart != self._chart:
                raise TypeError("Two functions not defined on the same " +
                                "chart cannot be added.")
            if isinstance(other, ZeroFunctionChart):
                return self.copy()
            res = simplify_chain(self._express + other._express)
        elif isinstance(other, (int, RingElement)):  #!# check
            res = simplify_chain(self._express + other)
        else:
            return other.__radd__(self)
        if res == 0:
            return self._chart._zero_function
        else:
            return FunctionChart(self._chart, res)

    def __radd__(self, other):
        r"""
        Addition on the left with ``other``.

        """
        return self.__add__(other)

    def __iadd__(self, other):
        r"""
        In-place addition operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        """
        return self.__add__(other)

    def __sub__(self, other):
        r"""
        Subtraction operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - chart function resulting from the subtraction of ``other`` from
          ``self``

        """
        if isinstance(other, FunctionChart):
            if other._chart != self._chart:
                raise TypeError("Two functions not defined on the same " +
                                "chart cannot be subtracted.")
            if isinstance(other, ZeroFunctionChart):
                return self.copy()
            res = simplify_chain(self._express - other._express)
        elif isinstance(other, (int, RingElement)):  #!# check
            res = simplify_chain(self._express - other)
        else:
            return other.__rsub__(self)
        if res == 0:
            return self._chart._zero_function
        else:
            return FunctionChart(self._chart, res)

    def __rsub__(self, other):
        r"""
        Subtraction from ``other``.

        """
        return (-self).__add__(other)

    def __isub__(self, other):
        r"""
        In-place subtraction operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        """
        return self.__sub__(other)


    def __mul__(self, other):
        r"""
        Multiplication operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - chart function resulting from the multiplication of ``self`` and
          ``other`

        """
        if isinstance(other, FunctionChart):
            if other._chart != self._chart:
                raise TypeError("Two functions not defined on the same " +
                                "chart cannot be multiplied.")
            if isinstance(other, ZeroFunctionChart):
                return self._chart._zero_function
            res = simplify_chain(self._express * other._express)
        elif isinstance(other, (int, RingElement)):  #!# check
            res = simplify_chain(self._express * other)
        else:
            return other.__rmul__(self)
        if res == 0:
            return self._chart._zero_function
        else:
            return FunctionChart(self._chart, res)

    def __rmul__(self, other):
        r"""
        Multiplication on the left by ``other``.

        """
        return self.__mul__(other)

    def __imul__(self, other):
        r"""
        In-place multiplication operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        """
        return self.__mul__(other)


    def __div__(self, other):
        r"""
        Division operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - chart function resulting from the division of ``self`` by
          ``other`

        """
        if isinstance(other, FunctionChart):
            if other._chart != self._chart:
                raise TypeError("Two functions not defined on the same " +
                                "chart cannot be divided.")
            if isinstance(other, ZeroFunctionChart):
                raise ZeroDivisionError("Division of a FunctionChart by zero.")
            res = simplify_chain(self._express / other._express)
        elif isinstance(other, (int, RingElement)):  #!# check
            res = simplify_chain(self._express / other)
        else:
            if other == 0:
                raise ZeroDivisionError("Division of a FunctionChart by zero.")
            return other.__rdiv__(self)
        if res == 0:
            return self._chart._zero_function
        else:
            return FunctionChart(self._chart, res)

    def __rdiv__(self, other):
        r"""
        Division of ``other`` by ``self``.

        """
        #!# to be improved
        res = simplify_chain(other / self._express)
        return FunctionChart(self._chart, res)


    def __idiv__(self, other):
        r"""
        In-place division operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        """
        return self.__div__(other)

    def factor(self):
        r"""
        Factorize the coordinate expression.

        OUTPUT:

        - ``self``, with ``self._express`` factorized

        EXAMPLES:

        Factorization on a 2-dimensional chart::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2 + 2*x*y + y^2)
            sage: f
            x^2 + 2*x*y + y^2
            sage: f.factor()
            (x + y)^2

        The method factor() has changed f::

            sage: f
            (x + y)^2

        """
        self._express = self._express.factor()
        self._del_derived()
        return self

    def simplify(self):
        r"""
        Simplifies the coordinate expression.

        OUTPUT:

        - ``self``, with ``self._express`` simplifyed

        EXAMPLES:

        Simplification on a 2-dimensional chart::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x)^2+sin(x)^2 + sqrt(x^2))
            sage: f
            cos(x)^2 + sin(x)^2 + sqrt(x^2)
            sage: f.simplify()
            abs(x) + 1

        The method simplify() has changed f::

            sage: f
            abs(x) + 1


        Another example::

            sage: f = X.function((x^2-1)/(x+1))
            sage: f
            (x^2 - 1)/(x + 1)
            sage: f.simplify()
            x - 1

        Examples taking into account the declared range of a coordinate::

            sage: M =  Manifold(2, 'M_1')
            sage: X.<x,y> = M.chart('x:(0,+oo) y')
            sage: f = X.function(sqrt(x^2))
            sage: f
            sqrt(x^2)
            sage: f.simplify()
            x

        ::

            sage: forget()  # to clear the previous assumption on x
            sage: M =  Manifold(2, 'M_2')
            sage: X.<x,y> = M.chart('x:(-oo,0) y')
            sage: f = X.function(sqrt(x^2))
            sage: f
            sqrt(x^2)
            sage: f.simplify()
            -x

        """
        self._express = simplify_chain(self._express)
        self._del_derived()
        return self

    def expand(self):
        r"""
        Expands the coordinate expression.

        OUTPUT:

        - ``self``, with ``self._express`` expanded

        EXAMPLES:

        Method `expand()` on a 2-dimensional chart function::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function((x - y)^2)
            sage: f
            (x - y)^2
            sage: f.expand()
            x^2 - 2*x*y + y^2

        Function f is changed after `expand()`::

            sage: f
            x^2 - 2*x*y + y^2
        """
        self._express = self._express.expand()
        self._del_derived()
        return self

    def collect(self, s):
        r"""
        Collect the coefficients of `s` into a group.

        OUTPUT:

        - ``self``, with ``self._express`` rearanged and coefficients of `s`
          collected.

        EXAMPLES:

        Method `collect()` on a 2-dimensional chart function::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2*y + x*y + (x*y)^2)
            sage: f
            x^2*y^2 + x^2*y + x*y
            sage: f.collect(y)
            x^2*y^2 + (x^2 + x)*y

        Function f is changed after `collect()`::

            sage: f
            x^2*y^2 + (x^2 + x)*y

        """
        self._express = self._express.collect(s)
        self._del_derived()
        return self

    def collect_common_factors(self):
        r"""
        Collect common factors (does not perform a full factorization,
        but only looks for factors which are already explicitly present).

        OUTPUT:

        - ``self``, with ``self._express`` rearanged and common factors
          collected.

        EXAMPLES:

        Method `collect_common_factors()` on a 2-dimensional chart function::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x/(x^2*y + x*y))
            sage: f
            x/(x^2*y + x*y)
            sage: f.collect_common_factors()
            1/((x + 1)*y)

        Function f is changed after `collect_common_factors()`::

            sage: f
            1/((x + 1)*y)

        """
        self._express = self._express.collect_common_factors()
        self._del_derived()
        return self


    def scalar_field(self, name=None, latex_name=None):
        r"""
        Construct the scalar field that has ``self`` as coordinate expression.

        The domain of the scalar field is the open subset covered by the chart
        on which ``self`` is defined.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the scalar field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the scalar
          field; if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of class
          :class:`~sage.geometry.manifolds.scalarfield.ScalarField`

        EXAMPLES:

        Construction of a scalar field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: fc = c_xy.function(x+2*y^3)
            sage: f = fc.scalar_field() ; f
            scalar field on the 2-dimensional manifold 'M'
            sage: f.display()
            M --> R
            (x, y) |--> 2*y^3 + x
            sage: f.function_chart(c_xy) is fc
            True

        """
        return self._chart._domain.scalar_field_algebra().element_class(
                     self._chart._domain, coord_expression={self._chart: self},
                     name=name, latex_name=latex_name)

    def exp(self):
        r"""
        Exponential of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.exp()
            e^(x + y)
            sage: exp(f) # equivalent to f.exp()
            e^(x + y)
            sage: exp(f).display()
            (x, y) |--> e^(x + y)
            sage: exp(X.zero_function())
            1

        """
        return FunctionChart(self._chart, simplify_chain(self._express.exp()))

    def log(self, base=None):
        r"""
        Logarithm of the function.

        INPUT:

        - ``base`` -- (default: ``None``) base of the logarithm; if None, the
          natural logarithm (i.e. logarithm to base e) is returned

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.log()
            log(x + y)
            sage: log(f) # equivalent to f.log()
            log(x + y)
            sage: log(f).display()
            (x, y) |--> log(x + y)
            sage: f.log(2)
            log(x + y)/log(2)
            sage: log(f, 2)
            log(x + y)/log(2)

        """
        return FunctionChart(self._chart,
                             simplify_chain(self._express.log(base)))

    def __pow__(self, exponent):
        r"""
        Power of the function.

        INPUT:

        - ``exponent`` -- the exponent

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.__pow__(3)
            x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: f^3  # equivalent to f.__pow__(3)
            x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: f.__pow__(3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: pow(f,3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: (f^3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: pow(X.zero_function(), 3).display()
            (x, y) |--> 0

        """
        return FunctionChart(self._chart,
                             simplify_chain(pow(self._express, exponent)))

    def sqrt(self):
        r"""
        Square root of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.sqrt()
            sqrt(x + y)
            sage: sqrt(f)  # equivalent to f.sqrt()
            sqrt(x + y)
            sage: sqrt(f).display()
            (x, y) |--> sqrt(x + y)
            sage: sqrt(X.zero_function()).display()
            (x, y) |--> 0

        """
        return FunctionChart(self._chart, simplify_chain(self._express.sqrt()))

    def cos(self):
        r"""
        Cosine of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.cos()
            cos(x*y)
            sage: cos(f)  # equivalent to f.cos()
            cos(x*y)
            sage: cos(f).display()
            (x, y) |--> cos(x*y)
            sage: cos(X.zero_function()).display()
            (x, y) |--> 1

        """
        return FunctionChart(self._chart, simplify_chain(self._express.cos()))

    def sin(self):
        r"""
        Sine of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.sin()
            sin(x*y)
            sage: sin(f)  # equivalent to f.sin()
            sin(x*y)
            sage: sin(f).display()
            (x, y) |--> sin(x*y)
            sage: sin(X.zero_function()) == X.zero_function()
            True

        """
        return FunctionChart(self._chart, simplify_chain(self._express.sin()))

    def tan(self):
        r"""
        Tangent of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.tan()
            sin(x*y)/cos(x*y)
            sage: tan(f)  # equivalent to f.tan()
            sin(x*y)/cos(x*y)
            sage: tan(f).display()
            (x, y) |--> sin(x*y)/cos(x*y)
            sage: tan(X.zero_function()) == X.zero_function()
            True

        """
        return FunctionChart(self._chart, simplify_chain(self._express.tan()))

    def arccos(self):
        r"""
        Arc cosine of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arccos()
            arccos(x*y)
            sage: arccos(f)  # equivalent to f.arccos()
            arccos(x*y)
            sage: acos(f)  # equivalent to f.arccos()
            arccos(x*y)
            sage: arccos(f).display()
            (x, y) |--> arccos(x*y)
            sage: arccos(X.zero_function()).display()
            (x, y) |--> 1/2*pi

        """
        return FunctionChart(self._chart,
                             simplify_chain(self._express.arccos()))

    def arcsin(self):
        r"""
        Arc sine of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arcsin()
            arcsin(x*y)
            sage: arcsin(f)  # equivalent to f.arcsin()
            arcsin(x*y)
            sage: asin(f)  # equivalent to f.arcsin()
            arcsin(x*y)
            sage: arcsin(f).display()
            (x, y) |--> arcsin(x*y)
            sage: arcsin(X.zero_function()) == X.zero_function()
            True

        """
        return FunctionChart(self._chart,
                             simplify_chain(self._express.arcsin()))

    def arctan(self):
        r"""
        Arc tangent of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arctan()
            arctan(x*y)
            sage: arctan(f)  # equivalent to f.arctan()
            arctan(x*y)
            sage: atan(f)  # equivalent to f.arctan()
            arctan(x*y)
            sage: arctan(f).display()
            (x, y) |--> arctan(x*y)
            sage: arctan(X.zero_function()) == X.zero_function()
            True

        """
        return FunctionChart(self._chart,
                             simplify_chain(self._express.arctan()))

    def cosh(self):
        r"""
        Hyperbolic cosine of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.cosh()
            cosh(x*y)
            sage: cosh(f)  # equivalent to f.cosh()
            cosh(x*y)
            sage: cosh(f).display()
            (x, y) |--> cosh(x*y)
            sage: cosh(X.zero_function()).display()
            (x, y) |--> 1

        """
        return FunctionChart(self._chart, simplify_chain(self._express.cosh()))

    def sinh(self):
        r"""
        Hyperbolic sine of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.sinh()
            sinh(x*y)
            sage: sinh(f)  # equivalent to f.sinh()
            sinh(x*y)
            sage: sinh(f).display()
            (x, y) |--> sinh(x*y)
            sage: sinh(X.zero_function()) == X.zero_function()
            True

        """
        return FunctionChart(self._chart, simplify_chain(self._express.sinh()))

    def tanh(self):
        r"""
        Hyperbolic tangent of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.tanh()
            sinh(x*y)/cosh(x*y)
            sage: tanh(f)  # equivalent to f.tanh()
            sinh(x*y)/cosh(x*y)
            sage: tanh(f).display()
            (x, y) |--> sinh(x*y)/cosh(x*y)
            sage: tanh(X.zero_function()) == X.zero_function()
            True

        """
        return FunctionChart(self._chart, simplify_chain(self._express.tanh()))

    def arccosh(self):
        r"""
        Inverse hyperbolic cosine of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arccosh()
            arccosh(x*y)
            sage: arccosh(f)  # equivalent to f.arccosh()
            arccosh(x*y)
            sage: acosh(f)  # equivalent to f.arccosh()
            arccosh(x*y)
            sage: arccosh(f).display()
            (x, y) |--> arccosh(x*y)
            sage: arccosh(X.function(1)) == X.zero_function()
            True

        """
        return FunctionChart(self._chart,
                             simplify_chain(self._express.arccosh()))

    def arcsinh(self):
        r"""
        Inverse hyperbolic sine of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arcsinh()
            arcsinh(x*y)
            sage: arcsinh(f)  # equivalent to f.arcsinh()
            arcsinh(x*y)
            sage: asinh(f)  # equivalent to f.arcsinh()
            arcsinh(x*y)
            sage: arcsinh(f).display()
            (x, y) |--> arcsinh(x*y)
            sage: arcsinh(X.zero_function()) == X.zero_function()
            True

        """
        return FunctionChart(self._chart,
                             simplify_chain(self._express.arcsinh()))

    def arctanh(self):
        r"""
        Inverse hyperbolic tangent of the function.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arctanh()
            arctanh(x*y)
            sage: arctanh(f)  # equivalent to f.arctanh()
            arctanh(x*y)
            sage: atanh(f)  # equivalent to f.arctanh()
            arctanh(x*y)
            sage: arctanh(f).display()
            (x, y) |--> arctanh(x*y)
            sage: arctanh(X.zero_function()) == X.zero_function()
            True

        """
        return FunctionChart(self._chart,
                             simplify_chain(self._express.arctanh()))


#*****************************************************************************

class ZeroFunctionChart(FunctionChart):
    r"""
    Null function defined on the codomain of a chart.

    Given a chart `(U, \varphi)` on a differentiable manifold `M` of
    dimension `n`, i.e. an open subset `U` of `M` along with an
    homeomorphism `\varphi: U \rightarrow V \subset \RR^n` (`V` being
    an open subset of `\RR^n`), the function

    .. MATH::

        \begin{array}{llcl}
        f:& V \subset\RR^n & \longrightarrow & \RR \\
          & (x^1,\ldots,x^n) & \longmapsto & 0
        \end{array}

    is implemented as an instance of class :class:`ZeroFunctionChart`.

    INPUT:

    - ``chart`` -- the chart `(U, \varphi)`

    EXAMPLES:

    Null function defined on a 2-dimensional chart::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: from sage.geometry.manifolds.chart import ZeroFunctionChart
        sage: f = ZeroFunctionChart(c_xy) ; f
        0
        sage: f.display()
        (x, y) |--> 0
        sage: f.expr()
        0
        sage: f.is_zero()
        True
        sage: f(1,2)
        0

    Each chart is equipped with its zero function::

        sage: c_xy.zero_function()
        0
        sage: type(c_xy.zero_function())
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        sage: f == c_xy.zero_function()
        True

    Arithmetics between instances of :class:`ZeroFunctionChart`::

        sage: g = ZeroFunctionChart(c_xy)
        sage: s = f+g ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f-g ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f*g ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f/g ; type(s) ; s
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a ZeroFunctionChart by zero.

    Arithmetics with a nonzero instance of :class:`FunctionChart`::

        sage: g = c_xy.function(x+y)
        sage: s = f+g ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x + y
        sage: s = g+f ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x + y
        sage: s = f-g ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        -x - y
        sage: s = g-f ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x + y
        sage: s = f*g ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = g*f ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f/g ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = g/f ; type(s) ; s
        Traceback (most recent call last):
        ...
        ZeroDivisionError: Division of a FunctionChart by zero.

    Arithmetics with a symbolic expression::

        sage: s = f+x ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x
        sage: s = x+f ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x
        sage: s = f-x ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        -x
        sage: s = x-f ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        x
        sage: s = f*x ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = x*f ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0
        sage: s = f/x ; type(s) ; s
        <class 'sage.geometry.manifolds.chart.ZeroFunctionChart'>
        0

    """
    def __init__(self, chart):
        FunctionChart.__init__(self, chart, 0)

    def copy(self):
        r"""
        Returns an exact copy of ``self``.

        The derived quantities are not copied, because they can be reconstructed
        if necessary.

        """
        return ZeroFunctionChart(self._chart)

    def __call__(self, *coords):
        r"""
        Computes the value of the function at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates `(x^1,...,x^n)` where the
          function `f` is to be evaluated

        OUTPUT:

        - the value `f(x^1,...,x^n)`

        """
        return 0    #!# SR(0) instead ?

    def diff(self, coord):
        r"""
        Partial derivative with respect to a coordinate.

        INPUT:

        - ``coord`` -- the coordinate `x^i` with respect
          to which the derivative of the function `f` is to be taken, or the
          index `i` labelling this coordinate

        OUTPUT:

        - the partial derivative `\frac{\partial f}{\partial x^i}`, as an
          instance of :class:`ZeroFunctionChart`

        """
        if self._der is None:
            self._der = [self._chart._zero_function for j in range(self._nc)]
        return self._der[0]

    def is_zero(self):
        r"""
        Always return ``True``.
        """
        return True

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - True if ``self`` is equal to ``other``,  or False otherwise

        """
        if isinstance(other, FunctionChart):
            if other._chart != self._chart:
                return False
            else:
                return other.is_zero()
        else:
            return bool(isinstance(other, (int, Integer)) and other==0)

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - ``self``

        """
        return self

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - `self`` (since ``self`` is zero)

        """
        return self

    def __add__(self, other):
        r"""
        Addition operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - chart function resulting from the addition of ``self`` and ``other``

        """
        if isinstance(other, FunctionChart):
            if other._chart != self._chart:
                raise TypeError("Two functions not defined on the same chart " +
                                "cannot be added.")
            return other.copy()
        elif isinstance(other, (int, RingElement)):  #!# check
            if other == 0:
                return self
            else:
                return FunctionChart(self._chart, other)
        else:
            return other.__radd__(self)

    def __sub__(self, other):
        r"""
        Subtraction operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - chart function resulting from the subtraction of ``other`` from
          ``self``

        """
        if isinstance(other, FunctionChart):
            if other._chart != self._chart:
                raise TypeError("Two functions not defined on the same chart " +
                                "cannot be subtracted.")
            return -other
        elif isinstance(other, (int, RingElement)):  #!# check
            if other == 0:
                return self
            else:
                return FunctionChart(self._chart, -other)
        else:
            return other.__rsub__(self)

    def __mul__(self, other):
        r"""
        Multiplication operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - chart function resulting from the multiplication of ``self`` and
          ``other`

        """
        if isinstance(other, (int, RingElement, FunctionChart)):  #!# check
            return self
        else:
            return other.__rmul__(self)

    def __div__(self, other):
        r"""
        Division operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - chart function resulting from the division of ``self`` by
          ``other`

        """
        if isinstance(other, (int, RingElement, FunctionChart)):  #!# check
            if other == 0:
                raise ZeroDivisionError("Division of a ZeroFunctionChart by " +
                                        "zero.")
            else:
                return self
        else:
            return other.__rdiv__(self)


    def scalar_field(self, name=None, latex_name=None):
        r"""
        Return the zero scalar field on the open subset covered by the chart on
        which ``self`` is defined.

        INPUT:

        - ``name`` -- (default: ``None``) unused
        - ``latex_name`` -- (default: ``None``) unused

        OUTPUT:

        - zero instance of class
          :class:`~sage.geometry.manifolds.scalarfield.ScalarField`

        EXAMPLES:

        Construction of a zero scalar field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: fc = c_xy.zero_function()
            sage: f = fc.scalar_field() ; f
            scalar field 'zero' on the 2-dimensional manifold 'M'
            sage: f.display()
            zero: M --> R
               (x, y) |--> 0

        """
        return self._chart._domain._zero_scalar_field

#*****************************************************************************

class MultiFunctionChart(SageObject):
    r"""
    Function from a chart codomain to `\RR^m`.

    Given a chart `(U, \varphi)` on a differentiable manifold `M` of
    dimension `n`, i.e. an open subset `U` of `M` along with an
    homeomorphism `\varphi: U \rightarrow V` (`V` = open subset of `\RR^n`),
    an instance of class :class:`MultiFunctionChart` is a function

    .. MATH::

        \begin{array}{llcl}
        f:& V \subset\RR^n & \longrightarrow & \RR^m \\
          & (x^1,\ldots,x^n) & \longmapsto & (f_1(x^1,\ldots,x^n),\ldots,
            f_m(x^1,\ldots,x^n))
        \end{array}

    Each function `f_i` is stored as an instance of
    :class:`~sage.geometry.manifolds.chart.FunctionChart`.

    INPUT:

    - ``chart`` -- the chart `(U, \varphi)`
    - ``*expressions`` -- the list of the coordinate expressions of the `m`
      functions `f_i` (`1\leq i \leq m`)

    EXAMPLES:

    A function `V\subset \RR^2 \longrightarrow \RR^3`::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y>  = M.chart()
        sage: f = c_xy.multifunction(x-y, x*y, cos(x)*exp(y)) ; f
        functions (x - y, x*y, cos(x)*e^y) on the chart (M, (x, y))
        sage: type(f)
        <class 'sage.geometry.manifolds.chart.MultiFunctionChart'>
        sage: f._functions
        (x - y, x*y, cos(x)*e^y)
        sage: f(x,y)
        (x - y, x*y, cos(x)*e^y)
        sage: latex(f)
        \left(x - y, x y, \cos\left(x\right) e^{y}\right)

    Each real-valued function `f_i` (`1\leq i \leq m`) composing `f` can be
    accessed via the square-bracket operator, by providing `i-1` as an
    argument::

        sage: f[0]
        x - y
        sage: f[1]
        x*y
        sage: f[2]
        cos(x)*e^y

    Each f[i-1] is an instance of :class:`FunctionChart`::

        sage: type(f[0])
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        sage: f[0].display()
        (x, y) |--> x - y

    An instance of class :class:`MultiFunctionChart` can represent a
    real-valued function (case `m=1`), although one should
    rather employ the class :class:`FunctionChart` for this purpose::

        sage: g = c_xy.multifunction(x*y^2)
        sage: g._functions
        (x*y^2,)

    Evaluating the functions at specified coordinates::

        sage: f(1,2)
        (-1, 2, cos(1)*e^2)
        sage: (a, b) = var('a b')
        sage: f(a,b)
        (a - b, a*b, cos(a)*e^b)
        sage: g(1,2)
        (4,)

    The Jacobian matrix::

        sage: f.jacobian()
        [[1, -1], [y, x], [-e^y*sin(x), cos(x)*e^y]]
        sage: g.jacobian()
        [[y^2, 2*x*y]]

    If `m=n`, the Jacobian determinant can be evaluated::

        sage: h = c_xy.multifunction(x-y, x*y)
        sage: h.jacobian_det()
        x + y

    """
    def __init__(self, chart, *expressions):
        if not isinstance(chart, Chart):
            raise TypeError("The first argument must be a chart.")
        self._chart = chart
        self._nc = len(self._chart._xx)    # number of coordinates
        self._nf = len(expressions)      # number of functions
        self._functions = tuple(FunctionChart(chart, expressions[i])
                                for i in range(self._nf))
        self._jacob = None
        self._jacob_matrix = None
        self._jacob_det = None

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "functions " + str(self._functions) + " on the " + \
                      str(self._chart)
        return description

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        return latex(self._functions)

    def expr(self):
        r"""
        Return the symbolic expression of the image of the `m` functions, as

        .. MATH::

            (f_1(x^1,\ldots,x^n),\ldots, f_m(x^1,\ldots,x^n))

        OUTPUT:

        - tuple of symbolic expressions corresponding to the above formula

        EXAMPLES:

        A set of 3 functions of 2 coordinates::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f.expr()
            (x - y, x*y, cos(x)*e^y)
            sage: type(f.expr()[0])
            <type 'sage.symbolic.expression.Expression'>
            sage: f.expr() == f(x,y)
            True

        """
        return tuple( self._functions[i]._express for i in range(self._nf) )

    def copy(self):
        r"""
        Return an exact copy of ``self``.

        The derived quantities (Jacobian matrix) are not copied, because they
        can be reconstructed if necessary.

        EXAMPLE:

        Copy of a set of 3 functions of 2 coordinates::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: g = f.copy() ; g
            functions (x - y, x*y, cos(x)*e^y) on the chart (M, (x, y))

        """
        return MultiFunctionChart(self._chart, *(self.expr()))

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- another instance of :class:`MultiFunctionChart`

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        EXAMPLES:

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f.__eq__(c_xy)  # other has a different type
            False
            sage: h = c_xy.multifunction(x-y, x*y)
            sage: f.__eq__(h)  # other has a different length
            False
            sage: h = c_xy.multifunction(x-y, x*y, 3+x)
            sage: f.__eq__(h)  # one of the functions of other is different
            False
            sage: h = c_xy.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f.__eq__(h)
            True
            sage: h = c_xy.multifunction(-y+x, y*x, cos(x)/exp(-y))
            sage: f.__eq__(h)  # equivalent expressions
            True

        """
        if not isinstance(other, MultiFunctionChart):
            return False
        if other._chart != self._chart:
            return False
        if other._nf != self._nf:
            return False
        for i in range(self._nf):
            if other._functions[i] != self._functions[i]:
                return False
        return True

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- another instance of :class:`FunctionChart` or a value

        OUTPUT:

        - True if ``self`` is different from ``other``,  or False otherwise

        EXAMPLES:

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f.__ne__(c_xy)  # other has a different type
            True
            sage: h = c_xy.multifunction(x-y, x*y)
            sage: f.__ne__(h)  # other has a different length
            True
            sage: h = c_xy.multifunction(x-y, x*y, 3+x)
            sage: f.__ne__(h)  # one of the functions of other is different
            True
            sage: h = c_xy.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f.__ne__(h)
            False
            sage: h = c_xy.multifunction(-y+x, y*x, cos(x)/exp(-y))
            sage: f.__ne__(h)  # equivalent expressions
            False

        """
        return not self.__eq__(other)

    def __getitem__(self, index):
        r"""
        Return a specified function of the set represented by ``self``.

        INPUT:

        -- ``index`` -- index `i` of the function (`0\leq i \leq m-1`)

        OUTPUT

        -- instance of :class:`FunctionChart` representing the function

        """
        return self._functions[index]

    def __call__(self, *coords, **options):
        r"""
        Compute the values of the functions at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates where the functions are to be
          evaluated
        - ``**options`` -- allows to pass ``simplify=False`` to disable the
          call of the simplification chain on the result

        OUTPUT:

        - the values of the `m` functions.

        """
        return tuple( self._functions[i](*coords, **options) for i in
                                                              range(self._nf) )

    def jacobian(self):
        r"""
        Return the Jacobian matrix of the system of functions.

        ``jacobian()`` is a 2-dimensional array of size `m\times n`
        where `m` is the number of functions and `n` the number of coordinates,
        the generic element being `J_{ij} = \frac{\partial f_i}{\partial x^j}`
        with `1\leq i \leq m` (row index) and `1\leq j \leq n` (column index).

        Each `J_{ij}` is an instance of :class:`FunctionChart`.

        OUTPUT:

        - Jacobian matrix as a 2-dimensional array J of FunctionChart's,
          J[i-1][j-1] being `J_{ij} = \frac{\partial f_i}{\partial x^j}`
          for `1\leq i \leq m` and `1\leq j \leq n`

        EXAMPLES:

        Jacobian of a set of 3 functions of 2 coordinates::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.multifunction(x-y, x*y, cos(x)*exp(y))
            sage: f.jacobian()
            [[1, -1], [y, x], [-e^y*sin(x), cos(x)*e^y]]
            sage: f.jacobian()[0][1]
            -1
            sage: type(f.jacobian()[0][1])
            <class 'sage.geometry.manifolds.chart.FunctionChart'>
            sage: f.jacobian()[0][1].display()
            (x, y) |--> -1

        """
        from sage.matrix.constructor import matrix
        from sage.calculus.functional import diff
        if self._jacob is None:
            self._jacob = [[ FunctionChart(self._chart,
                            simplify_chain(diff(self._functions[i]._express,
                                                self._chart._xx[j])) )
                    for j in range(self._nc) ] for i in range(self._nf) ]
            self._jacob_matrix = matrix( [[ self._jacob[i][j]._express
                    for j in range(self._nc) ] for i in range(self._nf) ] )
        return self._jacob

    def jacobian_det(self):
        r"""
        Return the Jacobian determinant of the system of functions.

        The number `m` of functions must equal the number `n` of
        coordinates.

        OUTPUT:

        - instance of :class:`FunctionChart` representing the determinant

        EXAMPLE:

        Jacobian determinant of a set of 2 functions of 2 coordinates::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = c_xy.multifunction(x-y, x*y)
            sage: f.jacobian_det()
            x + y

        The output of :meth:`jacobian_det` is an instance of
        :class:`FunctionChart` and can therefore be called on specific values
        of the coordinates, e.g. (x,y)=(1,2)::

            sage: type(f.jacobian_det())
            <class 'sage.geometry.manifolds.chart.FunctionChart'>
            sage: f.jacobian_det()(1,2)
            3

        """
        from utilities import simple_determinant
        if self._jacob_det is None:
            if (self._nf != self._nc):
                raise ValueError("The Jacobian matrix is not square.")
            self.jacobian() # to force the computation of self._jacob_matrix
            #!# the following is a workaround for a bug in Sage (cf. trac ticket #14403)
            self._jacob_det = FunctionChart(self._chart,
                       simplify_chain(simple_determinant(self._jacob_matrix)) )
            # the proper writing should be this:
            # self._jacob_det = FunctionChart(self._chart, simplify_chain(self._jacob_matrix.det()) )
        return self._jacob_det


#*****************************************************************************

class CoordChange(SageObject):
    r"""
    Change of coordinates (transition map between two charts).

    The two charts may belong to different manifolds.

    INPUT:

    - ``chart1`` -- initial chart
    - ``chart2`` -- final chart
    - ``transformations`` -- the coordinate transformations expressed as a list
      of the expressions of the "new" coordinates in terms of the "old" ones

    EXAMPLES:

    Change from spherical to Cartesian coordinates on `\RR^3`::

        sage: M = Manifold(3, 'R3', r'\mathbb{R}^3')
        sage: c_spher.<r,th,ph> = M.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
        sage: c_cart.<x,y,z> = M.chart()
        sage: ch = c_spher.coord_change(c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        sage: ch
        coordinate change from chart (R3, (r, th, ph)) to chart (R3, (x, y, z))
        sage: latex(ch)
        \left(\mathbb{R}^3,(r, {\theta}, {\phi})\right) \rightarrow
         \left(\mathbb{R}^3,(x, y, z)\right)
        sage: ch.display()
        x = r*cos(ph)*sin(th)
        y = r*sin(ph)*sin(th)
        z = r*cos(th)
        sage: type(ch)
        <class 'sage.geometry.manifolds.chart.CoordChange'>

    Each created coordinate change is automatically added to the manifold's
    dictionary :attr:`~sage.geometry.manifolds.domain.ManifoldSubset._coord_changes`;
    this dictionary is accessed via the method
    :meth:`~sage.geometry.manifolds.domain.ManifoldSubset.coord_change`::

        sage: M.coord_change(c_spher, c_cart)
        coordinate change from chart (R3, (r, th, ph)) to chart (R3, (x, y, z))

    It also generates a new entry in the manifold's dictionary
    :attr:`~sage.geometry.manifolds.domain.ManifoldSubset._frame_changes`,
    containing the relevant change-of-basis matrix;
    this dictionary is accessed via the method
    :meth:`~sage.geometry.manifolds.domain.ManifoldSubset.frame_change`::

        sage: M.frame_change(c_cart.frame(), c_spher.frame())
        field of tangent-space automorphisms on the 3-dimensional manifold 'R3'
        sage: M.frame_change(c_cart.frame(), c_spher.frame())[:]
        [   cos(ph)*sin(th)  r*cos(ph)*cos(th) -r*sin(ph)*sin(th)]
        [   sin(ph)*sin(th)  r*cos(th)*sin(ph)  r*cos(ph)*sin(th)]
        [           cos(th)         -r*sin(th)                  0]

    The coordinate change can be called directly on a set of "old" coordinates
    to get the "new" ones::

        sage: ch(1,pi/2,0)
        (1, 0, 0)

    The Jacobian matrix of the coordinate change::

        sage: ch._jacobian
        [[cos(ph)*sin(th), r*cos(ph)*cos(th), -r*sin(ph)*sin(th)],
         [sin(ph)*sin(th), r*cos(th)*sin(ph), r*cos(ph)*sin(th)],
         [cos(th), -r*sin(th), 0]]
        sage: ch._jacobian_det  # Jacobian determinant
        r^2*sin(th)

    Two successive change of coordinates can be composed by means of the
    operator \*, which in the present context stands for `\circ`::

        sage: c_cart2.<u,v,w> = M.chart()
        sage: ch2 = c_cart.coord_change(c_cart2, x+y, x-y, z-x-y)
        sage: ch3 = ch2 * ch ; ch3
        coordinate change from chart (R3, (r, th, ph)) to chart (R3, (u, v, w))
        sage: ch3(r,th,ph)
        (r*(cos(ph) + sin(ph))*sin(th),
         r*(cos(ph) - sin(ph))*sin(th),
         -r*(cos(ph) + sin(ph))*sin(th) + r*cos(th))
        sage: ch3 is M.coord_change(c_spher, c_cart2)
        True

    """
    def __init__(self, chart1, chart2, *transformations):
        n1 = len(chart1._xx)
        n2 = len(chart2._xx)
        if len(transformations) != n2:
            raise ValueError(str(n2) +
                             " coordinate transformations must be provided.")
        self._chart1 = chart1
        self._chart2 = chart2
        self._transf = MultiFunctionChart(chart1, *transformations)
        self._inverse = None
        # Jacobian matrix:
        self._jacobian  = self._transf.jacobian()
        # Jacobian determinant:
        if n1 == n2:
            self._jacobian_det = self._transf.jacobian_det()
        # If the two charts are on the same open subset, the coordinate change
        # is added to the subset (and supersets) dictionary and the
        # Jacobian matrix is added to the dictionary of changes of frame:
        if chart1._domain == chart2._domain:
            domain = chart1._domain
            for sdom in domain._supersets:
                sdom._coord_changes[(chart1, chart2)] = self
            frame1 = chart1._frame
            frame2 = chart2._frame
            vf_module = domain.vector_field_module()
            ch_basis = vf_module.automorphism()
            ch_basis.add_comp(frame1)[:, chart1] = self._jacobian
            ch_basis.add_comp(frame2)[:, chart1] = self._jacobian
            vf_module._basis_changes[(frame2, frame1)] = ch_basis
            for sdom in domain._supersets:
                sdom._frame_changes[(frame2, frame1)] = ch_basis
            # The inverse is computed only if it does not exist already
            # (because if it exists it may a simpler expression than that
            #  obtained from the matrix inverse)
            if (frame1, frame2) not in vf_module._basis_changes:
                ch_basis_inv = ch_basis.inverse()
                vf_module._basis_changes[(frame1, frame2)] = ch_basis_inv
                for sdom in domain._supersets:
                    sdom._frame_changes[(frame1, frame2)] = ch_basis_inv

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "coordinate change from " + str(self._chart1) + " to " + \
                      str(self._chart2)
        return description

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        return latex(self._chart1) + r' \rightarrow ' + latex(self._chart2)

    def __call__(self, *old_coords):
        r"""
        Computes the new coordinates from old ones.
        """
        return self._transf(*old_coords)

    def inverse(self):
        r"""
        Computes the inverse coordinate transformation, when the latter is
        invertible.

        OUTPUT:

        - an instance of :class:`CoordChange` representing the inverse of
          ``self``.

        EXAMPLES:

        Inverse of a coordinate transformation corresponding to a pi/3-rotation
        in the plane::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: ch_to_uv = c_xy.coord_change(c_uv, (x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2)
            sage: M._coord_changes
            {(chart (M, (x, y)), chart (M, (u, v))): coordinate change from chart (M, (x, y)) to chart (M, (u, v))}
            sage: ch_to_xy = ch_to_uv.inverse() ; ch_to_xy
            coordinate change from chart (M, (u, v)) to chart (M, (x, y))
            sage: ch_to_xy._transf
            functions (1/2*sqrt(3)*v + 1/2*u, -1/2*sqrt(3)*u + 1/2*v) on the chart (M, (u, v))
            sage: M._coord_changes # optional - dictionary_output
            {(chart (M, (u, v)), chart (M, (x, y))): coordinate change from chart (M, (u, v)) to chart (M, (x, y)),
            (chart (M, (x, y)), chart (M, (u, v))): coordinate change from chart (M, (x, y)) to chart (M, (u, v))}

        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        if self._inverse is not None:
            return self._inverse
        # The computation is necessary:
        x1 = self._chart1._xx  # list of coordinates in chart1
        x2 = self._chart2._xx  # list of coordinates in chart2
        n1 = len(x1)
        n2 = len(x2)
        if n1 != n2:
            raise TypeError("The change of coordinates is not invertible " +
                            "(different number of coordinates in the two " +
                            "charts).")
        # New symbolic variables (different from x2 to allow for a
        #  correct solution even when chart2 = chart1):
        coord_domain = ['real' for i in range(n2)]
        for i in range(n2):
            if x2[i].is_positive():
                coord_domain[i] = 'positive'
        xp2 = [ SR.var('xxxx' + str(i), domain=coord_domain[i])
                                                           for i in range(n2) ]
        equations = [ xp2[i] == self._transf._functions[i]._express
                                                           for i in range(n2) ]
        try:
            solutions = solve(equations, *x1, solution_dict=True)
        except RuntimeError:
            raise RuntimeError("The system could not be solved; use " +
                               "set_inverse() to set the inverse manually.")
        #!# This should be the Python 2.7 form:
        #           substitutions = {xp2[i]: x2[i] for i in range(n2)}
        # Here we use a form compatible with Python 2.6:
        substitutions = dict([(xp2[i], x2[i]) for i in range(n2)])
        if len(solutions) == 1:
            x2_to_x1 = [solutions[0][x1[i]].subs(substitutions)
                                                            for i in range(n1)]
            for transf in x2_to_x1:
                try:
                    transf = simplify_chain(transf)
                except AttributeError:
                    pass
        else:
            list_x2_to_x1 = []
            for sol in solutions:
                if x2[0] in sol:
                    raise ValueError("The system could not be solved; use " +
                                     "set_inverse() to set the inverse " +
                                     "manually.")
                x2_to_x1 = [sol[x1[i]].subs(substitutions) for i in range(n1)]
                for transf in x2_to_x1:
                    try:
                        transf = simplify_chain(transf)
                    except AttributeError:
                        pass
                if self._chart1.valid_coordinates(*x2_to_x1):
                    list_x2_to_x1.append(x2_to_x1)
            if len(list_x2_to_x1) == 0:
                raise ValueError("No solution found; use set_inverse() to " +
                                 "set the inverse manually.")
            if len(list_x2_to_x1) > 1:
                print "Multiple solutions found: "
                print list_x2_to_x1
                raise ValueError(
                   "Non-unique solution to the inverse coordinate " +
                   "transformation; use set_inverse() to set the inverse " +
                   "manually.")
            x2_to_x1 = list_x2_to_x1[0]
        self._inverse = CoordChange(self._chart2, self._chart1, *x2_to_x1)
        #
        # Update of chart expressions of the changes of frames:
        if self._chart1._domain == self._chart2._domain:
            vf_module = self._chart1._domain.vector_field_module()
            frame1 = self._chart1._frame
            frame2 = self._chart2._frame
            fr_change12 = vf_module._basis_changes[(frame1,frame2)]
            fr_change21 = vf_module._basis_changes[(frame2,frame1)]
            for comp in fr_change12._components[frame1]._comp.itervalues():
                comp.function_chart(self._chart1, from_chart=self._chart2)
            for comp in fr_change12._components[frame2]._comp.itervalues():
                comp.function_chart(self._chart1, from_chart=self._chart2)
            for comp in fr_change21._components[frame1]._comp.itervalues():
                comp.function_chart(self._chart2, from_chart=self._chart1)
            for comp in fr_change21._components[frame2]._comp.itervalues():
                comp.function_chart(self._chart2, from_chart=self._chart1)
        return self._inverse


    def set_inverse(self, *transformations, **kwds):
        r"""
        Sets the inverse of the coordinate transformation.

        This is usefull when the automatic computation via :meth:`inverse()`
        fails.

        INPUT:

        - ``transformations`` -- the inverse transformations expressed as a
          list of the expressions of the "old" coordinates in terms of the
          "new" ones
        - ``kwds`` -- keyword arguments: only ``check=True`` (default) or
          ``check=False`` are meaningfull; it determines whether the provided
          transformations are checked to be indeed the inverse coordinate
          transformations.

        EXAMPLES:

        From spherical coordinates to Cartesian ones in the plane::

            sage: M = Manifold(2, 'R^2')
            sage: U = M.open_subset('U') # the complement of the half line {y=0, x>= 0}
            sage: c_cart.<x,y> = U.chart()
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: spher_to_cart = c_spher.coord_change(c_cart, r*cos(ph), r*sin(ph))
            sage: spher_to_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x))
            Check of the inverse coordinate transformation:
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
               x == x
               y == y
            sage: spher_to_cart.inverse()
            coordinate change from chart (U, (x, y)) to chart (U, (r, ph))
            sage: M._coord_changes # random output order
            {(chart (U, (x, y)),
              chart (U, (r, ph))): coordinate change from chart (U, (x, y)) to chart (U, (r, ph)),
             (chart (U, (r, ph)),
              chart (U, (x, y))): coordinate change from chart (U, (r, ph)) to chart (U, (x, y))}

        Introducing a wrong inverse transformation is revealed by the check::

            sage: spher_to_cart.set_inverse(sqrt(x^3+y^2), atan2(y,x)) # note the x^3 typo
            Check of the inverse coordinate transformation:
               r == sqrt(r*cos(ph)^3 + sin(ph)^2)*r
               ph == arctan2(r*sin(ph), r*cos(ph))
               x == sqrt(x^3 + y^2)*x/sqrt(x^2 + y^2)
               y == sqrt(x^3 + y^2)*y/sqrt(x^2 + y^2)
            sage: # the check clearly fails

        """
        if 'check' in kwds:
            check = kwds['check']
        else:
            check = True
        self._inverse = CoordChange(self._chart2, self._chart1, *transformations)
        if check:
            print "Check of the inverse coordinate transformation:"
            x1 = self._chart1._xx
            x2 = self._chart2._xx
            n1 = len(x1)
            for i in range(n1):
                print "  ", x1[i], '==' , self._inverse(*(self(*x1)))[i]
            for i in range(n1):
                print "  ", x2[i], '==', self(*(self._inverse(*x2)))[i]
        # Update of chart expressions of the changes of frames:
        if self._chart1._domain == self._chart2._domain:
            vf_module = self._chart1._domain.vector_field_module()
            frame1 = self._chart1._frame
            frame2 = self._chart2._frame
            fr_change12 = vf_module._basis_changes[(frame1,frame2)]
            fr_change21 = vf_module._basis_changes[(frame2,frame1)]
            for comp in fr_change12._components[frame1]._comp.itervalues():
                comp.function_chart(self._chart1, from_chart=self._chart2)
            for comp in fr_change12._components[frame2]._comp.itervalues():
                comp.function_chart(self._chart1, from_chart=self._chart2)
            for comp in fr_change21._components[frame1]._comp.itervalues():
                comp.function_chart(self._chart2, from_chart=self._chart1)
            for comp in fr_change21._components[frame2]._comp.itervalues():
                comp.function_chart(self._chart2, from_chart=self._chart1)

    def __mul__(self, other):
        r"""
        Composition with another change of coordinates

        INPUT:

        - ``other`` -- another change of coordinate, the final chart of
          it is the initial chart of ``self``

        OUTPUT:

        - the change of coordinates X_1 --> X_3, where X_1 is the initial
          chart of ``other`` and X_3 is the final chart of ``self``

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: U.<u,v> = M.chart()
            sage: X_to_U = X.transition_map(U, (x+y, x-y))
            sage: W.<w,z> = M.chart()
            sage: U_to_W = U.transition_map(W, (u+cos(u)/2, v-sin(v)/2))
            sage: X_to_W = U_to_W * X_to_U ; X_to_W
            coordinate change from chart (M, (x, y)) to chart (M, (w, z))
            sage: X_to_W(x,y)
            (1/2*cos(x)*cos(y) - 1/2*sin(x)*sin(y) + x + y,
             -1/2*cos(y)*sin(x) + 1/2*cos(x)*sin(y) + x - y)

        """
        if not isinstance(other, CoordChange):
            raise TypeError(str(other) + " is not a change of coordinate.")
        if other._chart2 != self._chart1:
            raise ValueError("Composition not possible: " +
                             str(other._chart2) + " is different from "
                             + str(self._chart1))
        transf = self(*(other._transf.expr()))
        return CoordChange(other._chart1, self._chart2, *transf)

    def restrict(self, dom1, dom2=None):
        r"""
        Restriction to subsets.
        """
        if dom2 is None:
            dom2 = dom1
        return CoordChange(self._chart1.restrict(dom1),
                           self._chart2.restrict(dom2),
                           *(self._transf.expr()))

    def display(self):
        r"""
        Display of the coordinate transformation.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        EXAMPLE:

        From spherical coordinates to Cartesian ones in the plane::

            sage: M = Manifold(2, 'R^2')
            sage: U = M.open_subset('U') # the complement of the half line {y=0, x>= 0}
            sage: c_cart.<x,y> = U.chart()
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: spher_to_cart = c_spher.coord_change(c_cart, r*cos(ph), r*sin(ph))
            sage: spher_to_cart.display()
            x = r*cos(ph)
            y = r*sin(ph)
            sage: latex(spher_to_cart.display())
            \left\{\begin{array}{lcl} x & = & r \cos\left({\phi}\right) \\
             y & = & r \sin\left({\phi}\right) \end{array}\right.

        """
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion
        coords2 = self._chart2[:]
        n2 = len(coords2)
        expr = self._transf.expr()
        rtxt = ""
        if n2 == 1:
            rlatex = r"\begin{array}{lcl}"
        else:
            rlatex = r"\left\{\begin{array}{lcl}"
        for i in range(n2):
            x2 = coords2[i]
            x2f = expr[i]
            rtxt += repr(x2) + " = " + repr(x2f) + "\n"
            rlatex += latex(x2) + r" & = & " + latex(x2f) + r"\\"
        rtxt = rtxt[:-1]  # remove the last new line
        rlatex = rlatex[:-2] + r"\end{array}"
        if n2 > 1:
            rlatex += r"\right."
        return FormattedExpansion(rtxt, rlatex)



