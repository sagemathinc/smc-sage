r"""
Tangent spaces

The class :class:`TangentSpace` implements vector spaces tangent to a
differentiable manifold and the class :class:`TangentVector` implements the
corresponding tangent vectors.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version

REFERENCES:

- Chap. 3 of J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer
  (New York) (2013)

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

from sage.symbolic.ring import SR
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.tensor.modules.free_module_tensor import FiniteRankFreeModuleElement

class TangentVector(FiniteRankFreeModuleElement):
    r"""
    Tangent vector to a differentiable manifold at a given point.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.geometry.manifolds.tangentspace.TangentSpace`.
    It inherits from
    :class:`~sage.tensor.modules.free_module_tensor.FiniteRankFreeModuleElement`
    since :class:`~sage.geometry.manifolds.tangentspace.TangentSpace`
    inherits from
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`.

    INPUT:

    - ``parent`` -- the tangent space to which the vector belongs (must be an
      instance of
      :class:`~sage.geometry.manifolds.tangentspace.TangentSpace`)
    - ``name`` -- (default: ``None``) string; symbol given to the vector
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      the vector; if none is provided, ``name`` will be used

    EXAMPLES:

    Tangent vector on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: p = M.point((2,3), name='p')
        sage: Tp = p.tangent_space()
        sage: v = Tp((-2,1), name='v') ; v
        tangent vector v at point 'p' on 2-dimensional manifold 'M'
        sage: v.display()
        v = -2 d/dx + d/dy
        sage: v.parent()
        tangent space at point 'p' on 2-dimensional manifold 'M'
        sage: v in Tp
        True

    """
    def __init__(self, parent, name=None, latex_name=None):
        FiniteRankFreeModuleElement.__init__(self, parent, name=name, latex_name=latex_name)
        # Extra data (with respect to FiniteRankFreeModuleElement):
        self._point = parent._point

    def _repr_(self):
        r"""
        String representation of the object.
        """
        desc = "tangent vector"
        if self._name:
            desc += " " + str(self._name)
        desc += " at " + str(self._point)
        return desc

    def plot(self, chart=None, ambient_coords=None, mapping=None, scale=1,
             color='blue', print_label=True, label=None,  label_color=None,
             fontsize=10, label_offset=0.1, parameters=None, **extra_options):
        r"""
        Plot the vector in a Cartesian graph based on the coordinates of some
        ambient chart.

        The vector is drawn in terms of two (2D graphics) or three (3D graphics)
        coordinates of a given chart, called hereafter the *ambient chart*.
        The vector's base point `p` (or its image `\Phi(p)` by some
        differentiable mapping `\Phi`) must lie in the ambient chart's domain.
        If `\Phi` is different from the identity mapping, the vector
        actually depicted is `\mathrm{d}\Phi_p(v)`, where `v` is the current
        vector (``self``) (see the example of a vector tangent to the
        2-sphere below, where `\Phi: S^2 \rightarrow \RR^3`).

        INPUT:

        - ``chart`` -- (default: ``None``) the ambient chart (see above); if
          ``None``, it is set to the default chart of the open set containing
          the point at which the vector (or the vector image via the
          differential `\mathrm{d}\Phi_p` of ``mapping``) is defined
        - ``ambient_coords`` -- (default: ``None``) tuple containing the 2 or 3
          coordinates of the ambient chart in terms of which the plot is
          performed; if ``None``, all the coordinates of the ambient chart are
          considered
        - ``mapping`` -- (default: ``None``) differentiable mapping `\Phi`
          (instance of
          :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`)
          providing the link between the point `p` at which the vector is
          defined and the ambient chart ``chart``: the domain of ``chart`` must
          contain `\Phi(p)`; if ``None``, the identity mapping is assumed
        - ``scale`` -- (default: 1) value by which the length of the arrow
          representing the vector is multiplied
        - ``color`` -- (default: 'blue') color of the arrow representing the
          vector
        - ``print_label`` -- (boolean; default: ``True``) determines whether a
          label is printed next to the arrow representing the vector
        - ``label`` -- (string; default: ``None``) label printed next to the
          arrow representing the vector; if ``None``, the vector's symbol is
          used, if any
        - ``label_color`` -- (default: ``None``) color to print the label;
          if ``None``, the value of ``color`` is used
        - ``fontsize`` -- (default: 10) size of the font used to print the
          label
        - ``label_offset`` -- (default: 0.1) determines the separation between
          the vector arrow and the label
        - ``parameters`` -- (default: ``None``) dictionary giving the numerical
          values of the parameters that may appear in the coordinate expression
          of ``self`` (see example below)
        - ``**extra_options`` -- extra options for the arrow plot, like
          ``linestyle``, ``width`` or ``arrowsize`` (see
          :func:`~sage.plot.arrow.arrow2d` and
          :func:`~sage.plot.plot3d.shapes.arrow3d` for details)

        OUTPUT:

        - a graphic object, either an instance of
          :class:`~sage.plot.graphics.Graphics` for a 2D plot (i.e. based on
          2 coordinates of ``chart``) or an instance of
          :class:`~sage.plot.plot3d.base.Graphics3d` for a 3D plot (i.e.
          based on 3 coordinates of ``chart``)

        EXAMPLES:

        Vector tangent to a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,2), name='p')
            sage: Tp = p.tangent_space()
            sage: v = Tp((2, 1), name='v') ; v
            tangent vector v at point 'p' on 2-dimensional manifold 'M'

        Plot of the vector alone (arrow + label)::

            sage: v.plot()
            Graphics object consisting of 2 graphics primitives

        Plot atop of the chart grid::

            sage: show(X.plot() + v.plot())

        Plots with various options::

            sage: show(X.plot() + v.plot(color='green', scale=2, label='V'))
            sage: show(X.plot() + v.plot(print_label=False))
            sage: show(X.plot() + v.plot(color='green', label_color='black',
            ....:                         fontsize=20, label_offset=0.2))

        Plot with extra options::

            sage: show(X.plot() + v.plot(linestyle=':', width=4, arrowsize=8))

        Plot with specific values of some free parameters::

            sage: var('a b')
            (a, b)
            sage: v = Tp((1+a, -b^2), name='v') ; v.display()
            v = (a + 1) d/dx - b^2 d/dy
            sage: show(X.plot() + v.plot(parameters={a: -2, b: 3}))

        Special case of the zero vector::

            sage: v = Tp.zero() ; v
            tangent vector zero at point 'p' on 2-dimensional manifold 'M'
            sage: show(X.plot() + v.plot())

        Vector tangent to a 4-dimensional manifold::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: p = M((0,1,2,3), name='p')
            sage: Tp = p.tangent_space()
            sage: v = Tp((5,4,3,2), name='v') ; v
            tangent vector v at point 'p' on 4-dimensional manifold 'M'

        We cannot make a 4D plot directly::

            sage: v.plot()
            Traceback (most recent call last):
            ...
            ValueError: The number of coordinates involved in the plot must be either 2 or 3, not 4

        Rather, we have to select some chart coordinates for the plot, via
        the argument ``ambient_coords``. For instance, for a 2-dimensional plot
        in terms of the coordinates `(x,y)`::

            sage: v.plot(ambient_coords=(x,y))
            Graphics object consisting of 2 graphics primitives

        This plot involves only the components `v^x` and `v^y` of `v`.
        Similarly, for a 3-dimensional plot in terms of the coordinates
        `(t,x,y)`::

            sage: v.plot(ambient_coords=(t,x,z))
            Graphics3d Object

        This plot involves only the components `v^t`,  `v^x` and `v^z` of `v`.
        A nice 3D view atop the coordinate grid is obtained via::

            sage: show(X.plot(ambient_coords=(t,x,z)) +
            ....:      v.plot(ambient_coords=(t,x,z), label_offset=0.5, width=6))

        An example of plot via a differential mapping: plot of a vector tangent
        to a 2-sphere viewed in `\RR^3`::

            sage: S2 = Manifold(2, 'S^2')
            sage: U = S2.open_subset('U') # the open set covered by spherical coord.
            sage: XS.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: R3 = Manifold(3, 'R^3')
            sage: X3.<x,y,z> = R3.chart()
            sage: F = S2.diff_mapping(R3, {(XS, X3): [sin(th)*cos(ph),
            ....:                     sin(th)*sin(ph), cos(th)]}, name='F')
            sage: F.display() # the standard embedding of S^2 into R^3
            F: S^2 --> R^3
            on U: (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: p = U.point((pi/4, pi/4), name='p')
            sage: v = XS.frame()[1].at(p) ; v
            tangent vector d/dph at point 'p' on 2-dimensional manifold 'S^2'
            sage: graph_v = v.plot(mapping=F)
            sage: graph_S2 = XS.plot(chart=X3, mapping=F, nb_values=9)
            sage: show(graph_v + graph_S2)

        """
        from sage.plot.arrow import arrow2d
        from sage.plot.text import text
        from sage.plot.graphics import Graphics
        from sage.plot.plot3d.shapes import arrow3d
        from sage.plot.plot3d.shapes2 import text3d
        from sage.misc.functional import numerical_approx
        from sage.geometry.manifolds.chart import Chart
        #
        # The "effective" vector to be plotted
        #
        if mapping is None:
            eff_vector = self
            base_point = self._point
        else:
            # For efficiency, the method FiniteRankFreeModuleMorphism._call_()
            # is called instead of FiniteRankFreeModuleMorphism.__call__()
            eff_vector = mapping.differential(self._point)._call_(self)
            base_point = mapping(self._point)
        #
        # The chart w.r.t. which the vector is plotted
        #
        if chart is None:
            chart = base_point.containing_set().default_chart()
        elif not isinstance(chart, Chart):
            raise TypeError("{} is not a chart".format(chart))
        #
        # Coordinates of the above chart w.r.t. which the vector is plotted
        #
        if ambient_coords is None:
            ambient_coords = chart[:]  # all chart coordinates are used
        n_pc = len(ambient_coords)
        if n_pc != 2 and n_pc !=3:
            raise ValueError("The number of coordinates involved in the " +
                             "plot must be either 2 or 3, not {}".format(n_pc))
        # indices coordinates involved in the plot:
        ind_pc = [chart[:].index(pc) for pc in ambient_coords]
        #
        # Components of the vector w.r.t. the chart frame
        #
        basis = chart.frame().at(base_point)
        vcomp = eff_vector.comp(basis=basis)[:]
        xp = base_point.coord(chart=chart)
        #
        # The arrow
        #
        resu = Graphics()
        if parameters is None:
            coord_tail = [numerical_approx(xp[i]) for i in ind_pc]
            coord_head = [numerical_approx(xp[i] + scale*vcomp[i])
                          for i in ind_pc]
        else:
            coord_tail = [numerical_approx(
                           xp[i].substitute(parameters))
                          for i in ind_pc]
            coord_head = [numerical_approx(
                           (xp[i] + scale*vcomp[i]).substitute(parameters))
                          for i in ind_pc]
        if coord_head != coord_tail:
            if n_pc == 2:
                resu += arrow2d(tailpoint=coord_tail, headpoint=coord_head,
                                color=color, **extra_options)
            else:
                resu += arrow3d(coord_tail, coord_head, color=color,
                                **extra_options)
        #
        # The label
        #
        if print_label:
            if label is None:
                if n_pc == 2 and self._latex_name is not None:
                    label = r'$' + self._latex_name + r'$'
                if n_pc == 3 and self._name is not None:
                    label = self._name
            if label is not None:
                xlab = [xh + label_offset for xh in coord_head]
                if label_color is None:
                    label_color = color
                if n_pc == 2:
                    resu += text(label, xlab, fontsize=fontsize,
                                 color=label_color)
                else:
                    resu += text3d(label, xlab, fontsize=fontsize,
                                   color=label_color)
        return resu

#******************************************************************************

class TangentSpace(FiniteRankFreeModule):
    r"""
    Tangent space to a differentiable manifold at a given point.

    Since the tangent space at a given point `p` of a differentiable
    manifold `M` of dimension `n` is a `n`-dimensional vector space over `\RR`
    without any distinguished basis, the class :class:`TangentSpace` inherits
    from
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`,
    which implements free modules of finite rank (hence vector spaces of finite
    dimension) without any distinguished basis.

    This is a Sage *parent* class, the corresponding *element* class being
    :class:`~sage.geometry.manifolds.tangentspace.TangentVector`.

    INPUT:

    - ``point`` -- (instance of
      :class:`~sage.geometry.manifolds.point.ManifoldPoint`) point `p` at which the
      tangent space is defined.

    EXAMPLES:

    Tangent space on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: p = M.point((-1,2), name='p')
        sage: Tp = p.tangent_space() ; Tp
        tangent space at point 'p' on 2-dimensional manifold 'M'

    Tangent spaces actually belong to a dynamically generated subclass of
    :class:`TangentSpace`::

        sage: type(Tp)
        <class 'sage.geometry.manifolds.tangentspace.TangentSpace_with_category'>

    They are free modules of finite rank over Sage's Symbolic Ring (actually
    vector spaces of finite dimension over `\RR`)::

        sage: isinstance(Tp, FiniteRankFreeModule)
        True
        sage: Tp.base_ring()
        Symbolic Ring
        sage: Tp.category()
        Category of vector spaces over Symbolic Ring
        sage: Tp.rank()
        2
        sage: dim(Tp)
        2

    The tangent space is automatically endowed with bases deduced from the
    vector frames around the point::

        sage: Tp.bases()
        [Basis (d/dx,d/dy) on the tangent space at point 'p' on 2-dimensional manifold 'M']
        sage: M.frames()
        [coordinate frame (M, (d/dx,d/dy))]

    At this stage, only one basis has been defined in the tangent space, but
    new bases can be added from vector frames on the manifolds by means of the
    method :meth:`~sage.geometry.manifolds.vectorframe.VectorFrame.at`, for
    instance, from the frame associated with new coordinates on the manifold::

        sage: c_uv.<u,v> = M.chart()
        sage: c_uv.frame().at(p)
        Basis (d/du,d/dv) on the tangent space at point 'p' on 2-dimensional manifold 'M'
        sage: Tp.bases()
        [Basis (d/dx,d/dy) on the tangent space at point 'p' on 2-dimensional manifold 'M',
         Basis (d/du,d/dv) on the tangent space at point 'p' on 2-dimensional manifold 'M']

    All the bases defined on ``Tp`` are on the same footing. Accordingly the
    tangent space is not in the category of modules with a distinguished
    basis::

        sage: Tp in ModulesWithBasis(SR)
        False

    It is simply in the category of modules, as any instance of
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`::

        sage: Tp in Modules(SR)
        True

    Since the base ring is a field, it is actually in the category of
    vector spaces::

        sage: Tp in VectorSpaces(SR)
        True

    A typical element::

        sage: v = Tp.an_element() ; v
        tangent vector at point 'p' on 2-dimensional manifold 'M'
        sage: v.display()
        d/dx + 2 d/dy
        sage: v.parent()
        tangent space at point 'p' on 2-dimensional manifold 'M'

    The zero vector::

        sage: Tp.zero()
        tangent vector zero at point 'p' on 2-dimensional manifold 'M'
        sage: Tp.zero().display()
        zero = 0
        sage: Tp.zero().parent()
        tangent space at point 'p' on 2-dimensional manifold 'M'

    """

    Element = TangentVector

    def __init__(self, point):
        manif = point._manifold
        name = "T_" + str(point._name) + " " + str(manif._name)
        latex_name = r"T_{" + str(point._latex_name) + "}\," + \
                     str(manif._latex_name)
        self._point = point
        self._manif = manif
        FiniteRankFreeModule.__init__(self, SR, manif._dim, name=name,
                                      latex_name=latex_name,
                                      start_index=manif._sindex)
        # Initialization of bases of the tangent space from existing vector
        # frames around the point:
        for frame in point._subset._top_frames:
            if point in frame._domain:
                basis = self.basis(symbol=frame._symbol,
                                   latex_symbol=frame._latex_symbol)
                # Names of basis vectors set to those of the frame vector
                # fields:
                n = manif._dim
                for i in range(n):
                    basis._vec[i]._name = frame._vec[i]._name
                    basis._vec[i]._latex_name = frame._vec[i]._latex_name
                basis._name = "(" + \
                        ",".join([basis._vec[i]._name for i in range(n)]) + ")"
                basis._latex_name = r"\left(" + \
                     ",".join([basis._vec[i]._latex_name for i in range(n)])+ \
                     r"\right)"
                basis._symbol = basis._name
                basis._latex_symbol = basis._latex_name
                # Names of cobasis linear forms set to those of the coframe
                # 1-forms:
                coframe = frame.coframe()
                cobasis = basis.dual_basis()
                for i in range(n):
                    cobasis._form[i]._name = coframe._form[i]._name
                    cobasis._form[i]._latex_name = coframe._form[i]._latex_name
                cobasis._name = "(" + \
                     ",".join([cobasis._form[i]._name for i in range(n)]) + ")"
                cobasis._latex_name = r"\left(" + \
                  ",".join([cobasis._form[i]._latex_name for i in range(n)])+ \
                  r"\right)"
                point._frame_bases[frame] = basis
        # The basis induced by the default frame of the manifold subset
        # in which the point has been created is declared the default
        # basis of self:
        def_frame = point._subset._def_frame
        if def_frame in point._frame_bases:
            self._def_basis = point._frame_bases[def_frame]
        # Initialization of the changes of bases from the existing changes of
        # frames around the point:
        for frame_pair, automorph in point._subset._frame_changes.iteritems():
            frame1 = frame_pair[0] ; frame2 = frame_pair[1]
            fr1, fr2 = None, None
            for frame in point._frame_bases:
                if frame1 in frame._subframes:
                    fr1 = frame
                    break
            for frame in point._frame_bases:
                if frame2 in frame._subframes:
                    fr2 = frame
                    break
            if fr1 is not None and fr2 is not None:
                basis1 = point._frame_bases[fr1]
                basis2 = point._frame_bases[fr2]
                auto = self.automorphism()
                for frame, comp in automorph._components.iteritems():
                    basis = None
                    if frame is frame1:
                        basis = basis1
                    if frame is frame2:
                        basis = basis2
                    if basis is not None:
                        cauto = auto.add_comp(basis)
                        for ind, val in comp._comp.iteritems():
                            cauto._comp[ind] = val(point)
                self._basis_changes[(basis1, basis2)] = auto

    def _repr_(self):
        r"""
        String representation of the object.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((3,-2), name='p')
            sage: Tp = p.tangent_space()
        sage: Tp._repr_()
        "tangent space at point 'p' on 2-dimensional manifold 'M'"

        """
        description = "tangent space at {}".format(self._point)
        return description

    def _an_element_(self):
        r"""
        Construct some (unamed) vector in the tangent space

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((3,-2), name='p')
            sage: Tp = p.tangent_space()
            sage: Tp._an_element_()
            tangent vector at point 'p' on 2-dimensional manifold 'M'
            sage: Tp._an_element_().display()
            d/dx + 2 d/dy

        """
        resu = self.element_class(self)
        if self._def_basis is not None:
            resu.set_comp()[:] = range(1, self._rank+1)
        return resu

    def dimension(self):
        r"""
        Return the vector space dimension of ``self``.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = p.tangent_space()
            sage: Tp.dimension()
            2

        A shortcut is ``dim()``::

            sage: Tp.dim()
            2

        One can also use the global function ``dim``::

            sage: dim(Tp)
            2

        """
        return self._rank

    dim = dimension

    def base_point(self):
        r"""
        Return the manifold point at which ``self`` is defined.

        EXAMPLE::

            sage: from sage.geometry.manifolds.tangentspace import TangentSpace # for doctests only
            sage: TangentSpace._clear_cache_() ; Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = p.tangent_space()
            sage: Tp.base_point()
            point 'p' on 2-dimensional manifold 'M'
            sage: Tp.base_point() is p
            True

        """
        return self._point

    dim = dimension
