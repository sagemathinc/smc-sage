r"""
Manifold subsets

The class :class:`ManifoldSubset` implements subsets of a differentiable
manifold over `\RR`. The subclass :class:`ManifoldOpenSubset` is devoted
to open subsets, with respect to the manifold topology.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version

REFERENCES:

- J.M. Lee : *Introduction to Topological Manifolds*, 2nd ed., Springer (New
  York) (2011)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)

EXAMPLES:

Two subsets on a manifold::

    sage: M = Manifold(2, 'M')
    sage: a = M.subset('A') ; a
    subset 'A' of the 2-dimensional manifold 'M'
    sage: b = M.subset('B') ; b
    subset 'B' of the 2-dimensional manifold 'M'
    sage: M.subsets()  # random (set output)
    {subset 'A' of the 2-dimensional manifold 'M',
     subset 'B' of the 2-dimensional manifold 'M',
     2-dimensional manifold 'M'}

The intersection of the two subsets::

    sage: c = a.intersection(b) ; c
    subset 'A_inter_B' of the 2-dimensional manifold 'M'

Their union::

    sage: d = a.union(b) ; d
    subset 'A_union_B' of the 2-dimensional manifold 'M'

State of various data members after the above operations::

    sage: M.subsets()  # random (set output)
    {subset 'A' of the 2-dimensional manifold 'M',
     subset 'B' of the 2-dimensional manifold 'M',
     subset 'A_inter_B' of the 2-dimensional manifold 'M',
     2-dimensional manifold 'M'}
    sage: a.subsets()  # random (set output)
    set([subset 'A' of the 2-dimensional manifold 'M',
         subset 'A_inter_B' of the 2-dimensional manifold 'M'])
    sage: a._supersets  # random (set output)
    set([subset 'A_union_B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M',
         subset 'A' of the 2-dimensional manifold 'M'])
    sage: c._supersets  # random (set output)
    set([subset 'B' of the 2-dimensional manifold 'M',
         subset 'A_union_B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M',
         subset 'A' of the 2-dimensional manifold 'M',
         subset 'A_inter_B' of the 2-dimensional manifold 'M'])
    sage: c.subsets()  # random (set output)
    set([subset 'A_inter_B' of the 2-dimensional manifold 'M'])
    sage: d.subsets()  # random (set output)
    set([subset 'B' of the 2-dimensional manifold 'M',
         subset 'A_union_B' of the 2-dimensional manifold 'M',
         subset 'A_inter_B' of the 2-dimensional manifold 'M',
         subset 'A' of the 2-dimensional manifold 'M'])
    sage: d._supersets  # random (set output)
    set([subset 'A_union_B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M'])

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

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
from sage.categories.homset import Hom
from sage.rings.infinity import Infinity
from sage.geometry.manifolds.point import ManifoldPoint

class ManifoldSubset(UniqueRepresentation, Parent):
    r"""
    Subset of a differentiable manifold over `\RR`.

    The class :class:`ManifoldSubset` inherits from the generic Sage class
    :class:`~sage.structure.parent.Parent` and is declared to belong to
    the category of facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`).
    The corresponding element class is
    :class:`~sage.geometry.manifolds.point.ManifoldPoint`. A subset acts
    as a facade for the true parent of its points, which is the whole manifold
    (see example below).

    For an open subset, use the class :class:`ManifoldOpenSubset` instead.

    INPUT:

    - ``manifold`` -- manifold on which the subset is defined
    - ``name`` -- string; name (symbol) given to the subset
    - ``latex_name`` --  (default: ``None``) string: LaTeX symbol to denote the
      subset; if none is provided, it is set to ``name``

    EXAMPLES:

    A subset of a manifold::

        sage: Manifold._clear_cache_() # for doctests only
        sage: M = Manifold(2, 'M')
        sage: from sage.geometry.manifolds.domain import ManifoldSubset
        sage: A = ManifoldSubset(M, 'A', latex_name=r'\mathcal{A}') ; A
        subset 'A' of the 2-dimensional manifold 'M'
        sage: latex(A)
        \mathcal{A}
        sage: A.is_subset(M)
        True

    Instead of importing :class:`ManifoldSubset` in the global namespace,
    it is recommended to use the method :meth:`subset` to create a new subset::

        sage: B = M.subset('B', latex_name=r'\mathcal{B}') ; B
        subset 'B' of the 2-dimensional manifold 'M'
        sage: M.subsets()  # random (set output)
        {subset 'A' of the 2-dimensional manifold 'M',
         subset 'B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M'}

    The manifold is itself a subset::

        sage: isinstance(M, ManifoldSubset)
        True

    Actually, it is an instance of the subclass
    :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`, for it is
    (by definition) open::

        sage: isinstance(M, sage.geometry.manifolds.domain.ManifoldOpenSubset)
        True

    Instances of :class:`ManifoldSubset` are Sage's facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`):
    their elements are manifold points
    (class :class:`~sage.geometry.manifolds.point.ManifoldPoint`),
    which have the manifold (and not the subset) as parent::

        sage: isinstance(A, Parent)
        True
        sage: A.category()
        Category of facade sets
        sage: A.facade_for()
        (2-dimensional manifold 'M',)
        sage: p = A.an_element() ; p
        point on 2-dimensional manifold 'M'
        sage: p.parent()
        2-dimensional manifold 'M'
        sage: p in A
        True
        sage: p in M
        True

    """

    Element = ManifoldPoint

    def __init__(self, manifold, name, latex_name=None):
        r"""
        Construct a manifold subset.

        TESTS::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: A = M.subset('A') ; A
            subset 'A' of the 2-dimensional manifold 'M'

        """
        # Except for the manifold itself, the subsets are facade sets:
        if self is manifold:
            Parent.__init__(self, category=Sets())
        else:
            Parent.__init__(self, category=Sets(), facade=manifold)
            for dom in manifold._subsets:
                if name == dom._name:
                    raise ValueError("The name '" + name +
                                     "' is already used for another " +
                                     "subset of the {}".format(manifold))
            manifold._subsets.add(self)
        self._manifold = manifold
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._supersets = set([manifold, self]) # subsets containing self
        self._subsets = set([self]) # subsets of self
        self._top_subsets = set([self]) # subsets contained in self but not
                                        # in another strict subset of self
        self._intersections = {} # dict. of intersections with other subsets
                                 # (key: subset name)
        self._unions = {} # dict. of unions with other subsets (key: subset
                          # name)
        self._atlas = []  # list of charts defined on subsets of self
        self._top_charts = []  # list of charts defined on subsets of self
                        # that are not subcharts of charts on larger subsets
        self._def_chart = None  # default chart
        self._coord_changes = {} # dictionary of transition maps
        self._frames = []  # list of vector frames defined on subsets of self
        self._top_frames = []  # list of vector frames defined on subsets
               # of self that are not subframes of frames on larger subsets
        self._def_frame = None  # default frame
        self._frame_changes = {} # dictionary of changes of frames
        self._coframes = []  # list of coframes defined on subsets of self
        self._parallelizable_parts = set() # parallelizable subsets contained in self

    #### Methods required for any Parent in the category of sets:
    def _element_constructor_(self, coords=None, chart=None, name=None,
                              latex_name=None, check_coords=True):
        r"""
        Construct a point on ``self`` from its coordinates in some chart.
        """
        if isinstance(coords, ManifoldPoint):
            point = coords # for readability
            if point._subset is self:
                return point
            if point in self:
                resu = self.element_class(self, name=point._name,
                                          latex_name=point._latex_name)
                for chart, coords in point._coordinates.iteritems():
                    resu._coordinates[chart] = coords
                return resu
            else:
                raise ValueError("the {}".format(point) +
                                 " is not in {}".format(self))
        return self.element_class(self, coords=coords, chart=chart, name=name,
                              latex_name=latex_name, check_coords=check_coords)

    def _an_element_(self):
        r"""
        Construct some (unamed) point on ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M._an_element_() ; p
            point on 2-dimensional manifold 'M'
            sage: p.coord()
            (0, 0)
            sage: U = M.open_subset('U', coord_def={X: y>1}) ; U
            open subset 'U' of the 2-dimensional manifold 'M'
            sage: p = U._an_element_() ; p
            point on 2-dimensional manifold 'M'
            sage: p in U
            True
            sage: p.coord()
            (0, 2)
            sage: V = U.open_subset('V', coord_def={X.restrict(U): x<-pi})
            sage: p = V._an_element_() ; p
            point on 2-dimensional manifold 'M'
            sage: p in V
            True
            sage: p.coord()
            (-pi - 1, 2)

        """
        if self._def_chart is None:
            return self.element_class(self)
        # Attempt to construct a point in the domain of the default chart
        chart = self._def_chart
        coords = []
        for coord_range in chart._bounds:
            xmin = coord_range[0][0]
            xmax = coord_range[1][0]
            if xmin == -Infinity:
                if xmax == Infinity:
                    x = 0
                else:
                    x = xmax - 1
            else:
                if xmax == Infinity:
                    x = xmin + 1
                else:
                    x = (xmin + xmax)/2
            coords.append(x)
        if not chart.valid_coordinates(*coords):
            # Attempt to construct a point in the domain of other charts
            for ch in self._atlas:
                if ch is self._def_chart:
                    continue # since this case has already been attempted
                coords = []
                for coord_range in ch._bounds:
                    xmin = coord_range[0][0]
                    xmax = coord_range[1][0]
                    if xmin == -Infinity:
                        if xmax == Infinity:
                            x = 0
                        else:
                            x = xmax - 1
                    else:
                        if xmax == Infinity:
                            x = xmin + 1
                        else:
                            x = (xmin + xmax)/2
                    coords.append(x)
                if ch.valid_coordinates(*coords):
                    chart = ch
                    break
            else:
                # A generic element with specific coordinates could not be
                # automatically generated, due to too complex cooordinate
                # conditions. An element without any coordinate set is
                # returned instead:
                return self.element_class(self)
        # The point is constructed with check_coords=False since the check
        # has just been performed above:
        return self.element_class(self, coords=coords, chart=chart,
                                  check_coords=False)

    #### End of methods required for any Parent in the category of sets

    def _repr_(self):
        r"""
        String representation of ``self``.
        """
        return "subset '" + self._name + "' of the " + str(self._manifold)

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.
        """
        return self._latex_name

    def manifold(self):
        r"""
        Return the manifold of which ``self`` is a subset.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: A = M.subset('A')
            sage: A.manifold()
            2-dimensional manifold 'M'
            sage: B = A.subset('B')
            sage: B.manifold()
            2-dimensional manifold 'M'
            sage: M.manifold() is M
            True

        """
        return self._manifold

    def subsets(self):
        r"""
        Return the set of subsets that have been defined on ``self``.

        OUTPUT:

        - A Python set containing all the subsets that have been defined on
          ``self``.

        .. NOTE::

            To get the subsets as a list, used the method
            :meth:`list_of_subsets` instead.

        EXAMPLE:

        Subsets of a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: M.subsets()  # random (set output)
             {open subset 'U' of the 2-dimensional manifold 'M',
              subset 'V' of the 2-dimensional manifold 'M',
              2-dimensional manifold 'M'}
            sage: U in M.subsets()
            True

        The method :meth:`list_of_subsets` returns a list instead of a set::

            sage: M.list_of_subsets()
            [2-dimensional manifold 'M',
             open subset 'U' of the 2-dimensional manifold 'M',
             subset 'V' of the 2-dimensional manifold 'M']

        """
        return self._subsets

    def list_of_subsets(self):
        r"""
        Return the list of subsets that have been defined on ``self``.

        The list is sorted by the alphabetical names of the subsets.

        OUTPUT:

        - A list containing all the subsets that have been defined on
          ``self``.

        .. NOTE::

            To get the subsets as a Python set, used the method
            :meth:`subsets` instead.

        EXAMPLE:

        Subsets of a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: M.list_of_subsets()
            [2-dimensional manifold 'M',
             open subset 'U' of the 2-dimensional manifold 'M',
             subset 'V' of the 2-dimensional manifold 'M']

        The method :meth:`subsets` returns a set instead of a list::

            sage: M.subsets()  # random (set output)
            {open subset 'U' of the 2-dimensional manifold 'M',
             subset 'V' of the 2-dimensional manifold 'M',
             2-dimensional manifold 'M'}

        """
        return sorted(self._subsets, key = lambda x: x._name)

    def atlas(self):
        r"""
        Return the atlas of ``self``.

        OUTPUT:

        - list of charts defined on open subsets of ``self``.

        EXAMPLES:

        Charts on subsets of `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: M.atlas()
            [chart (R^2, (x, y))]
            sage: U = M.open_subset('U', coord_def={c_cart: (y!=0,x<0)}) # U = R^2 \ half line {y=0,x>=0}
            sage: U.atlas()
            [chart (U, (x, y))]
            sage: M.atlas()
            [chart (R^2, (x, y)), chart (U, (x, y))]
            sage: c_spher.<r, ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical (polar) coordinates on U
            sage: U.atlas()
            [chart (U, (x, y)), chart (U, (r, ph))]
            sage: M.atlas()
            [chart (R^2, (x, y)), chart (U, (x, y)), chart (U, (r, ph))]

        """
        return self._atlas

    def frames(self):
        r"""
        Return the list of vector frames defined on subsets of ``self``.

        OUTPUT:

        - list of vector frames defined on open subsets of ``self``.

        EXAMPLES:

        Vector frames on subsets of `\RR^2`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: M.frames()
            [coordinate frame (R^2, (d/dx,d/dy))]
            sage: e = M.vector_frame('e')
            sage: M.frames()
            [coordinate frame (R^2, (d/dx,d/dy)), vector frame (R^2, (e_0,e_1))]
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}) # Unit disk
            sage: U.frames()
            [coordinate frame (U, (d/dx,d/dy))]
            sage: M.frames()
            [coordinate frame (R^2, (d/dx,d/dy)),
             vector frame (R^2, (e_0,e_1)),
             coordinate frame (U, (d/dx,d/dy))]

        """
        return self._frames

    def coframes(self):
        r"""
        Return the list of coframes defined on subsets of ``self``.

        OUTPUT:

        - list of coframes defined on open subsets of ``self``.

        EXAMPLES:

        Coframes on subsets of `\RR^2`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: M.coframes()
            [coordinate coframe (R^2, (dx,dy))]
            sage: e = M.vector_frame('e')
            sage: M.coframes()
            [coordinate coframe (R^2, (dx,dy)), coframe (R^2, (e^0,e^1))]
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}) # Unit disk
            sage: U.coframes()
            [coordinate coframe (U, (dx,dy))]
            sage: e.restrict(U)
            vector frame (U, (e_0,e_1))
            sage: U.coframes()
            [coordinate coframe (U, (dx,dy)), coframe (U, (e^0,e^1))]
            sage: M.coframes()
            [coordinate coframe (R^2, (dx,dy)),
             coframe (R^2, (e^0,e^1)),
             coordinate coframe (U, (dx,dy)),
             coframe (U, (e^0,e^1))]

        """
        return self._coframes

    def coord_changes(self):
        r"""
        Return the changes of coordinates defined on ``self``.
        """
        return self._coord_changes

    def frame_changes(self):
        r"""
        Return the changes of vector frames defined on ``self``.
        """
        return self._frame_changes

    def subset(self, name, latex_name=None, is_open=False):
        r"""
        Create a subset of ``self``.

        INPUT:

        - ``name`` -- name given to the subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``is_open`` -- (default: False) if True, the created subset is
          assumed to be open with respect to the manifold's topology

        OUTPUT:

        - the subset, as an instance of :class:`ManifoldSubset`, or of
          :class:`ManifoldOpenSubset` if ``is_open`` is True.

        EXAMPLES:

        Creating a subset of a manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: a = M.subset('A') ; a
            subset 'A' of the 2-dimensional manifold 'M'

        Creating a subset of A::

            sage: b = a.subset('B', latex_name=r'\mathcal{B}') ; b
            subset 'B' of the 2-dimensional manifold 'M'
            sage: latex(b)
            \mathcal{B}

        B is then a subset of A and A is a superset of B::

            sage: a._subsets  # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])
            sage: b._supersets  # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])

        Creating an open subset of A::

            sage: c = a.subset('C', is_open=True) ; c
            open subset 'C' of the 2-dimensional manifold 'M'

        """
        if is_open:
            res = ManifoldOpenSubset(self._manifold, name, latex_name)
        else:
            res = ManifoldSubset(self._manifold, name, latex_name)
        res._supersets.update(self._supersets)
        for sd in self._supersets:
            sd._subsets.add(res)
        self._top_subsets.add(res)
        return res

    def superset(self, name, latex_name=None, is_open=False):
        r"""
        Create a superset of ``self``.

        A *superset* is a manifold subset in which ``self`` is included.

        INPUT:

        - ``name`` -- name given to the superset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          superset; if none is provided, it is set to ``name``
        - ``is_open`` -- (default: False) if True, the created subset is
          assumed to be open with respect to the manifold's topology

        OUTPUT:

        - the superset, as an instance of :class:`ManifoldSubset` or of
          :class:`ManifoldOpenSubset` if ``is_open==True``.

        EXAMPLES:

        Creating some superset of a given subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: a = M.subset('A')
            sage: b = a.superset('B') ; b
            subset 'B' of the 2-dimensional manifold 'M'
            sage: b._subsets # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])
            sage: a._supersets # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])

        The superset of the manifold is itself::

            sage: M.superset('SM') is M
            True

        Two supersets of a given subset are a priori different::

            sage: c = a.superset('C')
            sage: c == b
            False

        """
        if self is self._manifold:
            return self
        if is_open:
            res = ManifoldOpenSubset(self._manifold, name, latex_name)
        else:
            res = ManifoldSubset(self._manifold, name, latex_name)
        res._subsets.update(self._subsets)
        for sd in self._subsets:
            sd._supersets.add(res)
        res._atlas = list(self._atlas)
        res._top_charts = list(self._top_charts)
        res._coord_changes = dict(self._coord_changes)
        res._frames = list(self._frames)
        res._top_frames = list(self._top_frames)
        res._frame_changes = dict(self._frame_changes)
        res._coframes = list(self._coframes)
        res._def_chart = self._def_chart
        res._def_frame = self._def_frame
        return res


    def intersection(self, other, name=None, latex_name=None):
        r"""
        Return the intersection of ``self`` with another subset.

        INPUT:

        - ``other`` -- another subset of the same manifold
        - ``name`` -- (default: ``None``) name given to the intersection in the
          case the latter has to be created; the default is
          ``self._name`` inter ``other._name``
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          intersection in the case the latter has to be created; the default
          is built upon the symbol `\cap`

        OUTPUT:

        - instance of :class:`ManifoldSubset` (or :class:`ManifoldOpenSubset`
          if both subsets are open) representing the subset that is the
          intersection of ``self`` with ``other``

        EXAMPLES:

        Intersection of two subsets::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: a = M.subset('A')
            sage: b = M.subset('B')
            sage: c = a.intersection(b) ; c
            subset 'A_inter_B' of the 2-dimensional manifold 'M'
            sage: a._subsets  # random (set output)
            set([subset 'A_inter_B' of the 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])
            sage: b._subsets  # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 subset 'A_inter_B' of the 2-dimensional manifold 'M'])
            sage: c._supersets  # random (set output)
            set([subset 'A' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A_inter_B' of the 2-dimensional manifold 'M',
                 subset 'B' of the 2-dimensional manifold 'M'])

        Some checks::

            sage: (a.intersection(b)).is_subset(a)
            True
            sage: (a.intersection(b)).is_subset(a)
            True
            sage: a.intersection(b) is b.intersection(a)
            True
            sage: a.intersection(a.intersection(b)) is a.intersection(b)
            True
            sage: (a.intersection(b)).intersection(a) is a.intersection(b)
            True
            sage: M.intersection(a) is a
            True
            sage: a.intersection(M) is a
            True

        """
        if other._manifold != self._manifold:
            raise TypeError(
                "The two subsets do not belong to the same manifold.")
        # Particular cases:
        if self is self._manifold:
            return other
        if other is self._manifold:
            return self
        if self in other._subsets:
            return self
        if other in self._subsets:
            return other
        # Generic case:
        if other._name in self._intersections:
            # the intersection has already been created:
            return self._intersections[other._name]
        else:
            # the intersection must be created:
            if latex_name is None:
                if name is None:
                    latex_name = self._latex_name + r'\cap ' + other._latex_name
                else:
                    latex_name = name
            if name is None:
                name = self._name + "_inter_" + other._name
            if isinstance(self, ManifoldOpenSubset) and isinstance(other, ManifoldOpenSubset):
                res = self.open_subset(name, latex_name)
            else:
                res = self.subset(name, latex_name)
            res._supersets.update(other._supersets)
            for sd in other._supersets:
                sd._subsets.add(res)
            other._top_subsets.add(res)
            self._intersections[other._name] = res
            other._intersections[self._name] = res
            return res

    def union(self, other, name=None, latex_name=None):
        r"""
        Return the union of ``self`` with another subset.

        INPUT:

        - ``other`` -- another subset of the same manifold
        - ``name`` -- (default: ``None``) name given to the union in the
          case the latter has to be created; the default is
          ``self._name`` union ``other._name``
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          union in the case the latter has to be created; the default
          is built upon the symbol `\cup`

        OUTPUT:

        - instance of :class:`ManifoldSubset` (or :class:`ManifoldOpenSubset`
          if both subsets are open) representing the subset that is the
          union of ``self`` with ``other``

        EXAMPLES:

        Union of two subsets::

            sage: M = Manifold(2, 'M')
            sage: a = M.subset('A')
            sage: b = M.subset('B')
            sage: c = a.union(b) ; c
            subset 'A_union_B' of the 2-dimensional manifold 'M'
            sage: a._supersets  # random (set output)
            set([subset 'A_union_B' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])
            sage: b._supersets  # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A_union_B' of the 2-dimensional manifold 'M'])
            sage: c._subsets  # random (set output)
            set([subset 'A_union_B' of the 2-dimensional manifold 'M',
                subset 'A' of the 2-dimensional manifold 'M',
                subset 'B' of the 2-dimensional manifold 'M'])

        Some checks::

            sage: a.is_subset(a.union(b))
            True
            sage: b.is_subset(a.union(b))
            True
            sage: a.union(b) is b.union(a)
            True
            sage: a.union(a.union(b)) is a.union(b)
            True
            sage: (a.union(b)).union(a) is a.union(b)
            True
            sage: a.union(M) is M
            True
            sage: M.union(a) is M
            True

        """
        if other._manifold != self._manifold:
            raise TypeError(
                "The two subsets do not belong to the same manifold.")
        # Particular cases:
        if (self is self._manifold) or (other is self._manifold):
            return self._manifold
        if self in other._subsets:
            return other
        if other in self._subsets:
            return self
        # Generic case:
        if other._name in self._unions:
            # the union has already been created:
            return self._unions[other._name]
        else:
            # the union must be created:
            if latex_name is None:
                if name is None:
                    latex_name = self._latex_name + r'\cup ' + other._latex_name
                else:
                    latex_name = name
            if name is None:
                name = self._name + "_union_" + other._name
            res_open = isinstance(self, ManifoldOpenSubset) and \
                       isinstance(other, ManifoldOpenSubset)
            res = self.superset(name, latex_name, is_open=res_open)
            res._subsets.update(other._subsets)
            res._top_subsets.add(self)
            res._top_subsets.add(other)
            for sd in other._subsets:
                sd._supersets.add(res)
            for chart in other._atlas:
                if chart not in res._atlas:
                    res._atlas.append(chart)
            for chart in other._top_charts:
                if chart not in res._top_charts:
                    res._top_charts.append(chart)
            res._coord_changes.update(other._coord_changes)
            for frame in other._frames:
                if frame not in res._frames:
                    res._frames.append(frame)
            for frame in other._top_frames:
                if frame not in res._top_frames:
                    res._top_frames.append(frame)
            res._frame_changes.update(other._frame_changes)
            for coframe in other._coframes:
                if coframe not in res._coframes:
                    res._coframes.append(coframe)
            self._unions[other._name] = res
            other._unions[self._name] = res
            return res

    def declare_union(self, dom1, dom2):
        r"""
        Declare that ``self`` is the union of two subsets, i.e.
        that

        .. MATH::

            U = U_1 \cup U_2

        where `U` is ``self``,  `U_1\subset U` and `U_2\subset U`.

        INPUT:

        - ``dom1`` -- subset `U_1`
        - ``dom2`` -- subset `U_2`

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: M.declare_union(A, B)
            sage: A.union(B)
            2-dimensional manifold 'M'

        """
        if not dom1.is_subset(self):
            raise TypeError("The " + str(dom1) + " is not a subset of " +
                            "the " + str(self) + ".")
        if not dom2.is_subset(self):
            raise TypeError("The " + str(dom2) + " is not a subset of " +
                            "the " + str(self) + ".")
        dom1._unions[dom2._name] = self
        dom2._unions[dom1._name] = self

    def is_subset(self, other):
        r"""
        Return ``True`` iff ``self`` is included in ``other``.

        EXAMPLES:

        Subsubsets on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: a = M.subset('A')
            sage: b = a.subset('B')
            sage: c = M.subset('C')
            sage: a.is_subset(M)
            True
            sage: b.is_subset(a)
            True
            sage: b.is_subset(M)
            True
            sage: a.is_subset(b)
            False
            sage: c.is_subset(a)
            False

        """
        return self in other._subsets

    def __contains__(self, point):
        r"""
        Check whether a point is contained in ``self``.
        """
        # for efficiency, a quick test first:
        if point._subset is self:
            return True
        if point._subset.is_subset(self):
            return True
        for chart in self._atlas:
            if chart in point._coordinates:
                if chart.valid_coordinates( *(point._coordinates[chart]) ):
                    return True
        for chart in point._coordinates:
            for schart in chart._subcharts:
                if schart in self._atlas and schart.valid_coordinates(
                                          *(point._coordinates[chart]) ):
                    return True
        return False

    def point(self, coords=None, chart=None, name=None, latex_name=None):
        r"""
        Define a point in ``self``.

        See :class:`~sage.geometry.manifolds.point.ManifoldPoint` for a
        complete documentation.

        INPUT:

        - ``coords`` -- the point coordinates (as a tuple or a list) in the
          chart specified by ``chart``
        - ``chart`` -- (default: ``None``) chart in which the point coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the default chart of ``self``
        - ``name`` -- (default: ``None``) name given to the point
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the point;
          if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - the declared point, as an instance of
          :class:`~sage.geometry.manifolds.point.ManifoldPoint`.

        EXAMPLES:

        Points on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: p = M.point((1,2), name='p') ; p
            point 'p' on 2-dimensional manifold 'M'
            sage: p in M
            True
            sage: a = M.open_subset('A')
            sage: c_uv.<u,v> = a.chart()
            sage: q = a.point((-1,0), name='q') ; q
            point 'q' on 2-dimensional manifold 'M'
            sage: q in a
            True
            sage: p._coordinates
            {chart (M, (x, y)): (1, 2)}
            sage: q._coordinates
            {chart (A, (u, v)): (-1, 0)}

        """
        return self.element_class(self, coords=coords, chart=chart, name=name,
                                  latex_name=latex_name)

    def default_chart(self):
        r"""
        Return the default chart defined on ``self``.

        Unless changed via :meth:`set_default_chart`, the *default chart*
        is the first one defined on a subset of ``self`` (possibly itself).

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.chart.Chart`
          representing the default chart.

        EXAMPLES:

        Default chart on a 2-dimensional manifold and on some subsets::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: M.chart('x y')
            chart (M, (x, y))
            sage: M.chart('u v')
            chart (M, (u, v))
            sage: M.default_chart()
            chart (M, (x, y))
            sage: a = M.subset('A')
            sage: b = a.subset('B', is_open=True)
            sage: b.chart('t z')
            chart (B, (t, z))
            sage: a.default_chart()
            chart (B, (t, z))
            sage: b.default_chart()
            chart (B, (t, z))

        """
        return self._def_chart

    def set_default_chart(self, chart):
        r"""
        Changing the default chart on ``self``.

        INPUT:

        - ``chart`` -- a chart (must be defined on some subset ``self``)

        EXAMPLES:

        Charts on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: M.default_chart()
            chart (M, (x, y))
            sage: M.set_default_chart(c_uv)
            sage: M.default_chart()
            chart (M, (u, v))

        """
        from chart import Chart
        if not isinstance(chart, Chart):
            raise TypeError(str(chart) + " is not a chart.")
        if chart._domain is not self:
            if self.is_manifestly_coordinate_domain():
                raise TypeError("The chart domain must coincide with the " +
                                str(self) + ".")
            if chart not in self._atlas:
                raise ValueError("The chart must be defined on the " +
                                 str(self))
        self._def_chart = chart

    def coord_change(self, chart1, chart2):
        r"""
        Return a change of coordinates (transition map) defined on some
        subset of ``self``.

        INPUT:

        - ``chart1`` -- chart 1
        - ``chart2`` -- chart 2

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.chart.CoordChange`
          representing the transition map from chart 1 to chart 2

        EXAMPLES:

        Change of coordinates on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: c_xy.transition_map(c_uv, (x+y, x-y)) # defines the coordinate change
            coordinate change from chart (M, (x, y)) to chart (M, (u, v))
            sage: M.coord_change(c_xy, c_uv) # returns the coordinate change defined above
            coordinate change from chart (M, (x, y)) to chart (M, (u, v))

        """
        if (chart1, chart2) not in self._coord_changes:
            raise TypeError("The change of coordinates from " + str(chart1) +
                            " to " + str(chart2) + " has not been " +
                            "defined on the " + str(self))
        return self._coord_changes[(chart1, chart2)]


    def default_frame(self):
        r"""
        Return the default vector frame defined on ``self``.

        By 'vector frame' it is meant a field on ``self`` that provides,
        at each point p, a vector basis of the tangent space at p.

        Unless changed via :meth:`set_default_frame`, the default frame is the
        first one defined on ``self``, usually implicitely as the coordinate
        basis associated with the first chart defined on ``self``.

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.vectorframe.VectorFrame`
          representing the default vector frame.

        EXAMPLES:

        The default vector frame is often the coordinate frame associated
        with the first chart defined on ``self``::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: M.default_frame()
            coordinate frame (M, (d/dx,d/dy))

        """
        return self._def_frame

    def set_default_frame(self, frame):
        r"""
        Changing the default vector frame on ``self``.

        INPUT:

        - ``frame`` -- a vector frame (instance of
          :class:`~sage.geometry.manifolds.vectorframe.VectorFrame`) defined
          on ``self``

        EXAMPLE:

        Changing the default frame on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: e = M.vector_frame('e')
            sage: M.default_frame()
            coordinate frame (M, (d/dx,d/dy))
            sage: M.set_default_frame(e)
            sage: M.default_frame()
            vector frame (M, (e_0,e_1))

        """
        from vectorframe import VectorFrame
        if not isinstance(frame, VectorFrame):
            raise TypeError(str(frame) + " is not a vector frame.")
        if frame._domain is not self:
            if self.is_manifestly_parallelizable():
                raise TypeError("the frame domain must coincide with the " +
                                str(self))
            if not frame._domain.is_subset(self):
                raise TypeError("the frame must be defined on " + str(self))
        self._def_frame = frame
        frame._fmodule.set_default_basis(frame)

    def frame_change(self, frame1, frame2):
        r"""
        Return a change of vector frames defined on ``self``.

        INPUT:

        - ``frame1`` -- vector frame 1
        - ``frame2`` -- vector frame 2

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField`
          representing, at each point, the vector space automorphism `P`
          that relates frame 1, `(e_i)` say, to frame 2, `(n_i)` say,
          according to `n_i = P(e_i)`

        EXAMPLES:

        Change of vector frames induced by a change of coordinates::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: c_xy.transition_map(c_uv, (x+y, x-y))
            coordinate change from chart (M, (x, y)) to chart (M, (u, v))
            sage: M.frame_change(c_xy.frame(), c_uv.frame())
            field of tangent-space automorphisms on the 2-dimensional manifold 'M'
            sage: M.frame_change(c_xy.frame(), c_uv.frame())[:]
            [ 1/2  1/2]
            [ 1/2 -1/2]
            sage: M.frame_change(c_uv.frame(), c_xy.frame())
            field of tangent-space automorphisms on the 2-dimensional manifold 'M'
            sage: M.frame_change(c_uv.frame(), c_xy.frame())[:]
            [ 1  1]
            [ 1 -1]
            sage: M.frame_change(c_uv.frame(), c_xy.frame()) ==  M.frame_change(c_xy.frame(), c_uv.frame()).inverse()
            True

        In the present example, the manifold M is parallelizable, so that the
        module X(M) of vector fields on M is free. A change of frame on M is
        then identical to a change of basis in X(M)::

            sage: XM = M.vector_field_module() ; XM
            free module X(M) of vector fields on the 2-dimensional manifold 'M'
            sage: XM.print_bases()
            Bases defined on the free module X(M) of vector fields on the 2-dimensional manifold 'M':
             - (M, (d/dx,d/dy)) (default basis)
             - (M, (d/du,d/dv))
            sage: XM.change_of_basis(c_xy.frame(), c_uv.frame())
            field of tangent-space automorphisms on the 2-dimensional manifold 'M'
            sage: M.frame_change(c_xy.frame(), c_uv.frame()) is XM.change_of_basis(c_xy.frame(), c_uv.frame())
            True

        """
        if (frame1, frame2) not in self._frame_changes:
            raise TypeError("The change of frame from '" + repr(frame1) +
                            "' to '" + repr(frame2) + "' has not been " +
                            "defined on the " + repr(self))
        return self._frame_changes[(frame1, frame2)]

    def is_manifestly_coordinate_domain(self):
        r"""
        Returns True if ``self`` is known to be the domain of some coordinate
        chart and False otherwise.

        If False is returned, either ``self`` cannot be the domain of some
        coordinate chart or no such chart has been declared yet.
        """
        if not isinstance(self, ManifoldOpenSubset):
            return False
        return not self._covering_charts == []

    def is_manifestly_parallelizable(self):
        r"""
        Returns True if self is known to be a parallelizable subset and False
        otherwise.

        If False is returned, either self is not parallelizable or no vector
        frame has been defined on self yet.
        """
        if not isinstance(self, ManifoldOpenSubset):
            return False
        return not self._covering_frames == []

#******************************************************************************

class ManifoldOpenSubset(ManifoldSubset):
    r"""
    Open subset of a differentiable manifold over `\RR`.

    The class :class:`ManifoldOpenSubset` inherits naturally from the class
    :class:`ManifoldSubset`. Via the latter, it inherits from the generic
    Sage class :class:`~sage.structure.parent.Parent` and is declared to belong
    to the category of facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`).
    The corresponding element class is
    :class:`~sage.geometry.manifolds.point.ManifoldPoint`. An open subset acts
    as a facade for the true parent of its points, which is the whole manifold
    (see example below).

    INPUT:

    - ``manifold`` -- manifold on which the open subset is defined
    - ``name`` -- string; name (symbol) given to the open subset
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the open
      subset; if none is provided, it is set to ``name``

    EXAMPLES:

    An open subset of a 2-dimensional manifold::

        sage: Manifold._clear_cache_() # for doctests only
        sage: M = Manifold(2, 'M')
        sage: from sage.geometry.manifolds.domain import ManifoldOpenSubset
        sage: A = ManifoldOpenSubset(M, 'A', latex_name=r'\mathcal{A}') ; A
        open subset 'A' of the 2-dimensional manifold 'M'
        sage: latex(A)
        \mathcal{A}

    Instead of importing :class:`ManifoldOpenSubset` in the global namespace,
    it is recommended to use the method :meth:`open_subset` to create a new
    open subset::

        sage: B = M.open_subset('B', latex_name=r'\mathcal{B}') ; B
        open subset 'B' of the 2-dimensional manifold 'M'
        sage: M.subsets()  # random (set output)
        {open subset 'A' of the 2-dimensional manifold 'M',
         open subset 'B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M'}

    The manifold is itself an open subset (by definition!)::

        sage: isinstance(M, ManifoldOpenSubset)
        True

    Instances of :class:`ManifoldOpenSubset` are Sage's facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`):
    their elements are manifold points
    (class :class:`~sage.geometry.manifolds.point.ManifoldPoint`),
    which have the manifold (and not the open subset) as parent::

        sage: A.category()
        Category of facade sets
        sage: A.facade_for()
        (2-dimensional manifold 'M',)
        sage: p = A.an_element() ; p
        point on 2-dimensional manifold 'M'
        sage: p.parent()
        2-dimensional manifold 'M'

    Points can be created by providing their coordinates in some
    chart via the operator () applied to the subset::

        sage: X.<x,y> = A.chart()
        sage: p = A((-2,3)) ; p
        point on 2-dimensional manifold 'M'
        sage: p.coord()
        (-2, 3)

    Other arguments can be specified::

        sage: p = A((-2,3), chart=X, name='p') ; p
        point 'p' on 2-dimensional manifold 'M'

    It is equivalent to use the method
    :meth:`~sage.geometry.manifolds.domain.ManifoldSubset.point`::

        sage: A((-2,3)) == A.point((-2,3))
        True

    An open subset can be defined by some coordinate conditions::

        sage: U = A.open_subset('U', coord_def={X: x<0}) ; U
        open subset 'U' of the 2-dimensional manifold 'M'
        sage: U.is_subset(A)
        True
        sage: p in U  # since p's coordinates in chart X are (x,y)=(-2,3)
        True

    The coordinate conditions are taken into account when asking for a
    generic element::

        sage: q = U.an_element() ; q
        point on 2-dimensional manifold 'M'
        sage: q in U
        True
        sage: q.coord() # coordinates in U's default chart (x,y)
        (-1, 0)

    An open subset can be called to construct a point from another one,
    provided that the latter belongs to the open subset::

        sage: p1 = U(p) ; p1
        point 'p' on 2-dimensional manifold 'M'

    The two points simply differ by the returned value of
    :meth:`~sage.geometry.manifolds.point.ManifoldPoint.containing_set`::

        sage: p1 == p
        True
        sage: p.containing_set()
        open subset 'A' of the 2-dimensional manifold 'M'
        sage: p1.containing_set()
        open subset 'U' of the 2-dimensional manifold 'M'

    """
    def __init__(self, manifold, name, latex_name=None):
        r"""
        Construct an open subset of a manifold.

        TESTS::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x^2+y^2<1}) ; U
            open subset 'U' of the 2-dimensional manifold 'M'
            sage: U.category()
            Category of facade sets
            sage: TestSuite(U).run()

        """
        from sage.geometry.manifolds.scalarfield_algebra import \
                                                             ScalarFieldAlgebra
        ManifoldSubset.__init__(self, manifold, name, latex_name)
        # list of charts that individually cover self, i.e. whose
        # domains are self (if non-empty, self is a coordinate domain):
        self._covering_charts = []
        # list of vector frames that individually cover self, i.e. whose
        # domains are self (if non-empty, self is parallelizable):
        self._covering_frames = []
        # algebra of scalar fields defined on self (not contructed yet)
        self._scalar_field_algebra = ScalarFieldAlgebra(self)
        # The zero scalar field is constructed:
        self._zero_scalar_field = self._scalar_field_algebra.zero()
        # dict. of vector field modules along self
        # (keys = diff. mapping from self to an open set (possibly the identity
        #  map))
        self._vector_field_modules = {}
        # the identity map on self
        self._identity_map = Hom(self, self).one()

    def _repr_(self):
        r"""
        String representation of ``self``.
        """
        return "open subset '" + self._name + "' of the " + str(self._manifold)

    def open_subset(self, name, latex_name=None, coord_def={}):
        r"""
        Create an open subset of ``self``.

        An open subset is a set that is (i) included in ``self`` and (ii)
        open with respect to the manifold's topology.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts on ``self`` and values the symbolic expressions formed by the
          coordinates to define the subset.

        OUTPUT:

        - the open subset, as an instance of :class:`ManifoldOpenSubset`.

        EXAMPLES:

        Creating an open subset of a manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: a = M.open_subset('A') ; a
            open subset 'A' of the 2-dimensional manifold 'M'

        Creating an open subset of A::

            sage: b = a.open_subset('B') ; b
            open subset 'B' of the 2-dimensional manifold 'M'

        B is then a subset of A and A is a superset of B::

            sage: a._subsets # random (set output)
            set([open subset 'A' of the 2-dimensional manifold 'M',
                 open subset 'B' of the 2-dimensional manifold 'M'])
            sage: b._supersets # random (set output)
            set([open subset 'A' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 open subset 'B' of the 2-dimensional manifold 'M'])

        Defining an open subset by some coordinate restrictions: the open
        unit disk in `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}) ; U
            open subset 'U' of the 2-dimensional manifold 'R^2'

        Since the argument ``coord_def`` has been set, U is automatically
        provided with a chart, which is the restriction of the Cartesian one
        to U::

            sage: U.atlas()
            [chart (U, (x, y))]

        Therefore, one can immediately check whether a point belongs to U::

            sage: M.point((0,0)) in U
            True
            sage: M.point((1/2,1/3)) in U
            True
            sage: M.point((1,2)) in U
            False

        """
        resu = self.subset(name, latex_name=latex_name, is_open=True)
        for chart, restrictions in coord_def.iteritems():
            if chart not in self._atlas:
                raise ValueError("The " + str(chart) + "does not belong to " +
                    "the atlas of " + str(self))
            chart.restrict(resu, restrictions)
        return resu

    def chart(self, coordinates='', names=None):
        r"""
        Define a chart on ``self``.

        A *chart* is a pair `(U,\varphi)`, where `U` is the open subset
        represented by ``self`` and `\varphi: U \rightarrow V \subset \RR^n`
        is a homeomorphism from `U` to an open subset `V` of `\RR^n`.

        The components `(x^1,\ldots,x^n)` of `\varphi`, defined by
        `\varphi(p) = (x^1(p),\ldots,x^n(p))`, are called the *coordinates*
        of the chart `(U,\varphi)`.

        See :class:`~sage.geometry.manifolds.chart.Chart` for a complete
        documentation.

        INPUT:

        - ``coordinates`` -- single string defining the coordinate symbols and
          ranges: the coordinates are separated by ' ' (space) and each
          coordinate has at most three fields, separated by ':':

          1. The coordinate symbol (a letter or a few letters)
          2. (optional) The interval `I` defining the coordinate range: if not
             provided, the coordinate is assumed to span all `\RR`; otherwise
             `I` must be provided in the form (a,b) (or equivalently ]a,b[)
             The bounds a and b can be +/-Infinity, Inf, infinity, inf or oo.
             For *singular* coordinates, non-open intervals such as [a,b] and
             (a,b] (or equivalently ]a,b]) are allowed.
             Note that the interval declaration must not contain any space
             character.
          3. (optional) The LaTeX spelling of the coordinate; if not provided the
             coordinate symbol given in the first field will be used.

          The order of the fields 2 and 3 does not matter and each of them can
          be omitted.
          If it contains any LaTeX expression, the string ``coordinates`` must
          be declared with the prefix 'r' (for "raw") to allow for a proper
          treatment of the backslash character (see examples below).
          If no interval range and no LaTeX spelling is to be provided for any
          coordinate, the argument ``coordinates`` can be omitted when the
          shortcut operator <,> is used via Sage preparser (see examples below)
        - ``names`` -- (default: ``None``) unused argument, except if
          ``coordinates`` is not provided; it must then be a tuple containing
          the coordinate symbols (this is guaranted if the shortcut operator <,>
          is used).

        OUTPUT:

        - the created chart, as an instance of
          :class:`~sage.geometry.manifolds.chart.Chart`.

        EXAMPLES:

        Chart on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: X = U.chart('x y') ; X
            chart (U, (x, y))
            sage: X[0]
            x
            sage: X[1]
            y
            sage: X[:]
            (x, y)

        The declared coordinates are not known at the global level::

            sage: y
            Traceback (most recent call last):
            ...
            NameError: name 'y' is not defined

        They can be recovered by the Chart method [:]::

            sage: (x, y) = X[:]
            sage: y
            y
            sage: type(y)
            <type 'sage.symbolic.expression.Expression'>

        But a shorter way to proceed is to use the operator <,> in the
        left-hand side of the chart declaration (there is then no need to
        pass the string 'x y' to chart())::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: X.<x,y> = U.chart() ; X
            chart (U, (x, y))

        Indeed, the declared coordinates are then known at the global level::

            sage: y
            y
            sage: (x,y) == X[:]
            True

        Actually the instruction ``X.<x,y> = U.chart()`` is
        equivalent to the two instructions ``X = U.chart('x y')``
        and ``(x,y) = X[:]``.

        See the documentation of class
        :class:`~sage.geometry.manifolds.chart.Chart` for more examples,
        especially regarding the coordinates ranges and restrictions.

        """
        from chart import Chart
        return Chart(self, coordinates, names)

    def vector_frame(self, symbol=None, latex_symbol=None, dest_map=None,
                     from_frame=None):
        r"""
        Define a vector frame on ``self``.

        A *vector frame* is a field on ``self`` that provides, at each point
        p of ``self``, a vector basis of the tangent space at p.

        See :class:`~sage.geometry.manifolds.vectorframe.VectorFrame` for a
        complete documentation.

        INPUT:

        - ``symbol`` -- (default: ``None``) a letter (of a few letters) to denote a
          generic vector of the frame; can be set to None if the parameter
          ``from_frame`` is filled.
        - ``latex_symbol`` -- (default: ``None``) symbol to denote a generic vector
          of the frame; if None, the value of ``symbol`` is used.
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`
          (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`);
          if none is provided, the identity is assumed (case of a vector frame
          *on* `U`)
        - ``from_frame`` -- (default: ``None``) vector frame `\tilde e` on the
          codomain `V` of the destination map `\Phi`; the returned frame `e` is
          then such that `\forall p \in U, e(p) = \tilde e(\Phi(p))`

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.vectorframe.VectorFrame`
          representing the defined vector frame.

        EXAMPLES:

        Setting a vector frame on a 3-dimensional open subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: A = M.open_subset('A', latex_name=r'\mathcal{A}'); A
            open subset 'A' of the 3-dimensional manifold 'M'
            sage: c_xyz.<x,y,z> = A.chart()
            sage: e = A.vector_frame('e'); e
            vector frame (A, (e_0,e_1,e_2))
            sage: e[0]
            vector field 'e_0' on the open subset 'A' of the 3-dimensional manifold 'M'

        See the documentation of class
        :class:`~sage.geometry.manifolds.vectorframe.VectorFrame` for more
        examples.

        """
        from vectorframe import VectorFrame
        return VectorFrame(self.vector_field_module(dest_map=dest_map,
                                                    force_free=True),
                           symbol=symbol, latex_symbol=latex_symbol,
                           from_frame=from_frame)

    def _set_covering_frame(self, frame):
        r"""
        Declare a frame covering ``self``.
        """
        self._covering_frames.append(frame)
        self._parallelizable_parts = set([self])
        # if self cotained smaller parallelizable parts, they are forgotten
        for sd in self._supersets:
            if not sd.is_manifestly_parallelizable():
                sd._parallelizable_parts.add(self)

    def scalar_field_algebra(self):
        r"""
        Returns the algebra of scalar fields defined on ``self``.

        See :class:`~sage.geometry.manifolds.scalarfield_algebra.ScalarFieldAlgebra`
        for a complete documentation.

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.scalarfield_algebra.ScalarFieldAlgebra`
          representing the algebra `C^\infty(U)` of all scalar fields defined
          on `U` = ``self``.

        EXAMPLE:

        Scalar algebra of a 3-dimensional open subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: CU = U.scalar_field_algebra() ; CU
            algebra of scalar fields on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: CU.category()
            Category of commutative algebras over Symbolic Ring
            sage: CU.zero()
            scalar field 'zero' on the open subset 'U' of the 3-dimensional manifold 'M'

        """
        return self._scalar_field_algebra

    def vector_field_module(self, dest_map=None, force_free=False):
        r"""
        Returns the set of vector fields defined on ``self``, possibly
        within some ambient manifold, as a module over the algebra of scalar
        fields defined on ``self``.

        See :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldModule`
        for a complete documentation.

        INPUT:

        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `U` is ``self``
          (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`);
          if none is provided, the identity map is assumed (case of vector
          fields *on* `U`)
        - ``force_free`` -- (default: False) if set to True, force the
          construction of a *free* module (this implies that `V` is
          parallelizable)

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldModule`
          (or of
          :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldFreeModule`
          if `V` is parallelizable)
          representing the module `\mathcal{X}(U,\Phi)` of vector fields on the
          open subset `U`=``self`` taking values on
          `\Phi(U)\subset V\subset M`.

        EXAMPLES:

        Vector field module `\mathcal{X}(U)` of the complement `U` of the two
        poles on the sphere `\mathbb{S}^2`::

            sage: S2 = Manifold(2, 'S^2')
            sage: U = S2.open_subset('U')  # the complement of the two poles
            sage: spher_coord.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coordinates
            sage: XU = U.vector_field_module() ; XU
            free module X(U) of vector fields on the open subset 'U' of the 2-dimensional manifold 'S^2'
            sage: XU.category()
            Category of modules over algebra of scalar fields on the open subset 'U' of the 2-dimensional manifold 'S^2'
            sage: XU.base_ring()
            algebra of scalar fields on the open subset 'U' of the 2-dimensional manifold 'S^2'
            sage: XU.base_ring() is U.scalar_field_algebra()
            True

        `\mathcal{X}(U)` is a free module because `U` is parallelizable (being
        a chart domain)::

            sage: U.is_manifestly_parallelizable()
            True

        Its rank is the manifold's dimension::

            sage: XU.rank()
            2

        The elements of `\mathcal{X}(U)` are vector fields on `U`::

            sage: XU.an_element()
            vector field on the open subset 'U' of the 2-dimensional manifold 'S^2'
            sage: XU.an_element().display()
            2 d/dth + 2 d/dph

        Vector field module `\mathcal{X}(U,\Phi)` of the
        `\mathbb{R}^3`-valued vector fields along `U`, associated with the
        embedding `\Phi` of `\mathbb{S}^2` into `\mathbb{R}^3`::

            sage: R3 = Manifold(3, 'R^3')
            sage: cart_coord.<x, y, z> = R3.chart()
            sage: Phi = U.diff_mapping(R3, [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)], name='Phi')
            sage: XU_R3 = U.vector_field_module(dest_map=Phi) ; XU_R3
            free module X(U,Phi) of vector fields along the open subset 'U' of the 2-dimensional manifold 'S^2' mapped into the 3-dimensional manifold 'R^3'
            sage: XU_R3.base_ring()
            algebra of scalar fields on the open subset 'U' of the 2-dimensional manifold 'S^2'

        `\mathcal{X}(U,\mathbb{R}^3)` is a free module because `\mathbb{R}^3`
        is parallelizable and its rank is 3::

            sage: XU_R3.rank()
            3

        """
        from vectorfield_module import VectorFieldModule, VectorFieldFreeModule
        if dest_map is None:
            dest_map = self._identity_map
        codomain = dest_map._codomain
        if dest_map not in self._vector_field_modules:
            if codomain.is_manifestly_parallelizable() or force_free:
                self._vector_field_modules[dest_map] = \
                                 VectorFieldFreeModule(self, dest_map=dest_map)
            else:
                self._vector_field_modules[dest_map] = \
                                     VectorFieldModule(self, dest_map=dest_map)
        return self._vector_field_modules[dest_map]

    def tensor_field_module(self, tensor_type, dest_map=None):
        r"""
        Return the set of tensor fields of a given type defined on ``self``,
        possibly within some ambient manifold, as a module over the algebra of
        scalar fields defined on ``self``.

        See :class:`~sage.geometry.manifolds.tensorfield_module.TensorFieldModule`
        for a complete documentation.

        INPUT:

        - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant
          rank and `l` the covariant rank
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `U` is ``self``
          (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`);
          if none is provided, the identity map is assumed (case of tensor
          fields *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.tensorfield_module.TensorFieldModule`
          representing the module `\mathcal{T}^{(k,l)}(U,\Phi)` of type-`(k,l)`
          tensor fields on the open subset `U` = ``self`` taking values on
          `\Phi(U)\subset V\subset M`.

        EXAMPLE:

        Module of type-(2,1) tensor fields on a 3-dimensional open subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: c_xyz.<x,y,z> = U.chart()
            sage: TU = U.tensor_field_module((2,1)) ; TU
            free module T^(2,1)(U) of type-(2,1) tensors fields on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: TU.category()
            Category of modules over algebra of scalar fields on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: TU.base_ring()
            algebra of scalar fields on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: TU.base_ring() is U.scalar_field_algebra()
            True
            sage: TU.an_element()
            tensor field of type (2,1) on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: TU.an_element().display()
            2 d/dx*d/dx*dx

        """
        return self.vector_field_module(dest_map=dest_map).tensor_module(
                                                                  *tensor_type)

    def diff_form_module(self, degree, dest_map=None):
        r"""
        Return the set of differential forms of a given degree defined on
        ``self``, possibly within some ambient manifold, as a module over the
        algebra of scalar fields defined on ``self``.

        See :class:`~sage.geometry.manifolds.diffform_module.DiffFormModule`
        for a complete documentation.

        INPUT:

        - ``degree`` -- positive integer; the degree `p` of the differential forms
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `U` is ``self``
          (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`);
          if none is provided, the identity map is assumed (case of
          differential forms *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.diffform_module.DiffFormModule`
          representing the module `\Lambda^p(U,\Phi)` of `p`-forms on the open
          subset `U` = ``self`` taking values on `\Phi(U)\subset V\subset M`.

        EXAMPLE:

        """
        return self.vector_field_module(dest_map=dest_map).dual_exterior_power(
                                                                        degree)

    def automorphism_field_group(self, dest_map=None):
        r"""
        Return the group of tangent-space automorphism fields defined on
        ``self``, possibly within some ambient manifold.

        If `U` stands for the open subset ``self`` and `\Phi` is a
        differentiable mapping `\Phi: U \rightarrow V = \Phi(U) \subset M`,
        where `M` is a manifold, this method called with ``dest_map``
        being `\Phi` returns the general linear group
        `\mathrm{GL}(\mathcal{X}(U,\Phi))` of the module
        `\mathcal{X}(U,\Phi)` of vector fields along `U` with values in
        `V=\Phi(U)`.

        INPUT:

        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ U \rightarrow V`, where `U` is ``self``
          (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`);
          if none is provided, the identity map is assumed

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.automorphismfield_group.AutomorphismFieldParalGroup`
          (if ``self`` is parallelizable) or of
          :class:`~sage.geometry.manifolds.automorphismfield_group.AutomorphismFieldGroup`
          (if ``self`` is not parallelizable) representing
          `\mathrm{GL}(\mathcal{X}(U,\Phi))`

        EXAMPLE:

        """
        return self.vector_field_module(dest_map=dest_map).general_linear_group()

    def _Hom_(self, other, category=None):
        r"""
        Construct the set of morphisms (i.e. differentiable maps)
        ``self`` --> ``other``.

        INPUT:

        - ``other`` -- an open subset of some manifold
        - ``category`` -- (default: ``None``) not used here (to ensure
          compatibility with generic hook ``_Hom_``)

        OUTPUT:

        - the homset Hom(U,V), where U is ``self`` and V is ``other``

        See class
        :class:`~sage.geometry.manifolds.manifold_homset.ManifoldHomset`
        for more documentation.

        """
        from sage.geometry.manifolds.manifold_homset import ManifoldHomset
        return ManifoldHomset(self, other)

    def scalar_field(self, coord_expression=None, chart=None, name=None,
                     latex_name=None):
        r"""
        Define a scalar field on the open set.

        See :class:`~sage.geometry.manifolds.scalarfield.ScalarField` for a
        complete documentation.

        INPUT:

        - ``coord_expression`` -- (default: ``None``) coordinate expression(s)
          of the scalar field; this can be either

          - a single coordinate expression; if the argument ``chart`` is
            ``'all'``, this expression is set to all the charts defined
            on the open set; otherwise, the expression is set in the
            specific chart provided by the argument ``chart``
          - a dictionary of coordinate expressions, with the charts as keys.

          If ``coord_expression`` is ``None`` or does not fully specified the
          scalar field, other coordinate expressions can be added subsequently
          by means of the methods
          :meth:`~sage.geometry.manifolds.scalarfield.ScalarField.add_expr`,
          :meth:`~sage.geometry.manifolds.scalarfield.ScalarField.add_expr_by_continuation`,
          or :meth:`~sage.geometry.manifolds.scalarfield.ScalarField.set_expr`
        - ``chart`` -- (default: ``None``) chart defining the coordinates used
          in ``coord_expression`` when the latter is a single coordinate
          expression; if none is provided (default), the default chart of the
          open set is assumed. If ``chart=='all'``, ``coord_expression`` is
          assumed to be independent of the chart (constant scalar field).
        - ``name`` -- (default: ``None``) name given to the scalar field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the scalar
          field; if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.scalarfield.ScalarField`
          representing the defined scalar field.

        EXAMPLES:

        A scalar field defined by its coordinate expression in the open
        set's default chart::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: c_xyz.<x,y,z> = U.chart()
            sage: f = U.scalar_field(sin(x)*cos(y) + z, name='F'); f
            scalar field 'F' on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: f.display()
            F: U --> R
               (x, y, z) |--> cos(y)*sin(x) + z
            sage: f.parent()
            algebra of scalar fields on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: f in U.scalar_field_algebra()
            True

        Equivalent definition with the chart specified::

            sage: f = U.scalar_field(sin(x)*cos(y) + z, chart=c_xyz, name='F')
            sage: f.display()
            F: U --> R
               (x, y, z) |--> cos(y)*sin(x) + z

        Equivalent definition with a dictionary of coordinate expression(s)::

            sage: f = U.scalar_field({c_xyz: sin(x)*cos(y) + z}, name='F')
            sage: f.display()
            F: U --> R
               (x, y, z) |--> cos(y)*sin(x) + z

        See the documentation of class
        :class:`~sage.geometry.manifolds.scalarfield.ScalarField` for more
        examples.

        """
        if isinstance(coord_expression, dict):
            # check validity of entry
            for chart in coord_expression:
                if not chart._domain.is_subset(self):
                    raise ValueError("The " + str(chart) + " is not defined " +
                                     "on some subset of the " + str(self))
        elif coord_expression is not None and chart != 'all':
            # coord_expression is valid only in a specific chart
            if chart is None:
                chart = self._def_chart
            coord_expression = {chart: coord_expression}
        return self.scalar_field_algebra()._element_constructor_(
                                            coord_expression=coord_expression,
                                            name=name, latex_name=latex_name)

    def zero_scalar_field(self):
        r"""
        Return the zero scalar field defined on the open subset.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.zero_scalar_field() ; f
            scalar field 'zero' on the 2-dimensional manifold 'M'
            sage: f.display()
            zero: M --> R
               (x, y) |--> 0
            sage: f.parent()
            algebra of scalar fields on the 2-dimensional manifold 'M'
            sage: f is M.scalar_field_algebra().zero()
            True

        """
        return self._zero_scalar_field

    def vector_field(self, name=None, latex_name=None, dest_map=None):
        r"""
        Define a vector field on ``self``.

        See :class:`~sage.geometry.manifolds.vectorfield.VectorField` for a
        complete documentation.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the vector field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the vector
          field; if none is provided, the LaTeX symbol is set to ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of a vector field *on* ``self``)

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.vectorfield.VectorField`
          representing the defined vector field.

        EXAMPLES:

        A vector field on a 3-dimensional open subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: c_xyz.<x,y,z> = U.chart()
            sage: v = U.vector_field('v'); v
            vector field 'v' on the open subset 'U' of the 3-dimensional manifold 'M'

        Vector fields on `U` form the set `\mathcal{X}(U)`, which is a module
        over the algebra `C^\infty(U)` of smooth scalar fields on `U`::

            sage: v.parent()
            free module X(U) of vector fields on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: v in U.vector_field_module()
            True

        See the documentation of class
        :class:`~sage.geometry.manifolds.vectorfield.VectorField` for more
        examples.

        """
        vmodule = self.vector_field_module(dest_map)  # the parent
        return vmodule.element_class(vmodule, name=name, latex_name=latex_name)

    def tensor_field(self, k, l, name=None, latex_name=None, sym=None,
        antisym=None, dest_map=None):
        r"""
        Define a tensor field on ``self``.

        See :class:`~sage.geometry.manifolds.tensorfield.TensorField` for a
        complete documentation.

        INPUT:

        - ``k`` -- the contravariant rank, the tensor type being (k,l)
        - ``l`` -- the covariant rank, the tensor type being (k,l)
        - ``name`` -- (default: ``None``) name given to the tensor field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the tensor
          field; if none is provided, the LaTeX symbol is set to ``name``
        - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries among
          the tensor arguments: each symmetry is described by a tuple containing
          the positions of the involved arguments, with the convention position=0
          for the first argument. For instance:

          * sym=(0,1) for a symmetry between the 1st and 2nd arguments
          * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments.

        - ``antisym`` -- (default: ``None``) antisymmetry or list of antisymmetries
          among the arguments, with the same convention as for ``sym``.
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of a tensor field *on* ``self``)

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.tensorfield.TensorField`
          representing the defined tensor field.

        EXAMPLES:

        A tensor field of type (2,0) on a 3-dimensional open subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: c_xyz.<x,y,z> = U.chart()
            sage: t = U.tensor_field(2, 0, 'T'); t
            tensor field 'T' of type (2,0) on the open subset 'U' of the 3-dimensional manifold 'M'

        Type-(2,0) tensor fields on `U` form the set `\mathcal{T}^{(2,0)}(U)`,
        which is a module over the algebra `C^\infty(U)` of smooth scalar
        fields on `U`::

            sage: t.parent()
            free module T^(2,0)(U) of type-(2,0) tensors fields on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: t in U.tensor_field_module((2,0))
            True

        See the documentation of class
        :class:`~sage.geometry.manifolds.tensorfield.TensorField` for more
        examples.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.tensor((k,l), name=name, latex_name=latex_name, sym=sym,
                              antisym=antisym)

    def sym_bilin_form_field(self, name=None, latex_name=None, dest_map=None):
        r"""
        Define a field of symmetric bilinear forms on ``self``.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the field;
          if none is provided, the LaTeX symbol is set to ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of a field of symmetric bilinear forms *on* ``self``)

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.tensorfield.TensorField`
          of tensor type (0,2) and symmetric representing the defined
          symmetric bilinear form field.

        EXAMPLE:

        A field of symmetric bilinear forms on a 3-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: t = M.sym_bilin_form_field('T'); t
            field of symmetric bilinear forms 'T' on the 3-dimensional manifold 'M'

        Such a object is a tensor field of rank 2 and type (0,2)::

            sage: t.parent()
            free module T^(0,2)(M) of type-(0,2) tensors fields on the 3-dimensional manifold 'M'
            sage: t._tensor_rank
            2
            sage: t._tensor_type
            (0, 2)

        The LaTeX symbol is deduced from the name or can be specified when
        creating the object::

            sage: latex(t)
            T
            sage: om = M.sym_bilin_form_field('Omega', r'\Omega')
            sage: latex(om)
            \Omega

        Components with respect to some vector frame::

            sage: e = M.vector_frame('e') ; M.set_default_frame(e)
            sage: t.set_comp()
            Fully symmetric 2-indices components w.r.t. vector frame (M, (e_0,e_1,e_2))
            sage: type(t.comp())
            <class 'sage.tensor.modules.comp.CompFullySym'>

        For the default frame, the components are accessed with the
        square brackets::

            sage: t[0,0], t[0,1], t[0,2] = (1, 2, 3)
            sage: t[1,1], t[1,2] = (4, 5)
            sage: t[2,2] = 6

        The other components are deduced by symmetry::

            sage: t[1,0], t[2,0], t[2,1]
            (2, 3, 5)
            sage: t[:]
            [1 2 3]
            [2 4 5]
            [3 5 6]

        A symmetric bilinear form acts on vector pairs::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: t = M.sym_bilin_form_field('T')
            sage: t[0,0], t[0,1], t[1,1] = (-1, x, y*x)
            sage: v1 = M.vector_field('V_1')
            sage: v1[:] = (y,x)
            sage: v2 = M.vector_field('V_2')
            sage: v2[:] = (x+y,2)
            sage: s = t(v1,v2) ; s
            scalar field 'T(V_1,V_2)' on the 2-dimensional manifold 'M'
            sage: s.expr()
            x^3 + (3*x^2 + x)*y - y^2
            sage: s.expr() - t[0,0]*v1[0]*v2[0] - t[0,1]*(v1[0]*v2[1]+v1[1]*v2[0]) - t[1,1]*v1[1]*v2[1]
            0
            sage: latex(s)
            T\left(V_1,V_2\right)

        Adding two symmetric bilinear forms results in another symmetric
        bilinear form::

            sage: a = M.sym_bilin_form_field()
            sage: a[0,0], a[0,1], a[1,1] = (1,2,3)
            sage: b = M.sym_bilin_form_field()
            sage: b[0,0], b[0,1], b[1,1] = (-1,4,5)
            sage: s = a + b ; s
            field of symmetric bilinear forms on the 2-dimensional manifold 'M'
            sage: s[:]
            [0 6]
            [6 8]

        But adding a symmetric bilinear from with a non-symmetric bilinear form
        results in a generic type (0,2) tensor::

            sage: c = M.tensor_field(0,2)
            sage: c[:] = [[-2, -3], [1,7]]
            sage: s1 = a + c ; s1
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: s1[:]
            [-1 -1]
            [ 3 10]
            sage: s2 = c + a ; s2
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: s2[:]
            [-1 -1]
            [ 3 10]

        """
        return self.tensor_field(0, 2, name=name, latex_name=latex_name,
                                 sym=(0,1))

    def automorphism_field(self, name=None, latex_name=None,
                           dest_map=None):
        r"""
        Define a field of automorphisms (invertible endomorphisms in each
        tangent space) on ``self``.

        See :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField` for
        a complete documentation.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the field;
          if none is provided, the LaTeX symbol is set to ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of an automorphism field *on* ``self``)

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField`
          (or of
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismFieldParal`
          if the codomain is parallelizable) representing the defined field of
          automorphisms.

        EXAMPLE:

        A field of automorphisms on a 3-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3,'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: a = M.automorphism_field('A') ; a
            field of tangent-space automorphisms 'A' on the 3-dimensional
             manifold 'M'
            sage: a.parent()
            General linear group of the free module X(M) of vector fields on
             the 3-dimensional manifold 'M'

        See the documentation of class
        :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField` for more
        examples.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.automorphism(name=name, latex_name=latex_name)

    def tangent_identity_field(self, name='Id', latex_name=None,
                               dest_map=None):
        r"""
        Return the field of identity maps in the tangent spaces on ``self``.

        INPUT:

        - ``name`` -- (string; default: 'Id') name given to the field of
          identity maps
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the field of identity map; if none is provided, the LaTeX symbol is
          set to '\mathrm{Id}' if ``name`` is 'Id' and to ``name`` otherwise
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of a tangent-space identity map field *on* ``self``)

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismField`
          (or of
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismFieldParal`
          if the codomain is parallelizable) representing the field of identity
          maps.

        EXAMPLE:

        Field of tangent-space identity maps on a 3-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: a = M.tangent_identity_field(); a
            field of tangent-space identity maps on the 3-dimensional
             manifold 'M'
            sage: a.comp()
            Kronecker delta of size 3x3

        See the documentation of class
        :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismFieldField`
        for more examples.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.identity_map(name=name, latex_name=latex_name)

    def curve(self, coord_expression, param, chart=None, name=None,
              latex_name=None):
        r"""
        Define a differentiable curve in ``self``.

        See :class:`~sage.geometry.manifolds.curve.ManifoldCurve` for details.

        INPUT:

        - ``coord_expression`` -- either

          - (i) a dictionary whose keys are charts on ``self`` and values
            the coordinate expressions (as lists or tuples) of the curve in
            the given chart
          - (ii) a single coordinate expression in a given chart on ``self``,
            the latter being provided by the argument ``chart``

          In both cases, if the dimension of the arrival manifold is 1,
          a single coordinate expression can be passed instead of a tuple with
          a single element
        - ``param`` -- a tuple of the type ``(t, t_min, t_max)``, where ``t``
          is the curve parameter used in ``coord_expression``, ``t_min`` is its
          minimal value and ``t_max`` its maximal value; if ``t_min=-Infinity``
          and ``t_max=+Infinity``, they can be omitted and ``t`` can be passed
          for ``param``, instead of the tuple ``(t, t_min, t_max)``
        - ``chart`` -- (default: ``None``) chart on ``self`` used for case (ii)
          above; if ``None`` the default chart of ``self`` is assumed.
        - ``name`` -- (default: ``None``) string; symbol given to the curve
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
          the curve; if none is provided, ``name`` will be used

        OUTPUT:

        - instance of
          :class:`~sage.geometry.manifolds.curve.ManifoldCurve`

        EXAMPLES:

        The lemniscate of Gerono in the 2-dimensional Euclidean plane::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c') ; c
            Curve 'c' in the 2-dimensional manifold 'M'

        The same definition with the coordinate expression passed as a
        dictionary::

            sage: c = M.curve({X: [sin(t), sin(2*t)/2]}, (t, 0, 2*pi), name='c') ; c
            Curve 'c' in the 2-dimensional manifold 'M'

        An example of definition with ``t_min`` and ``t_max`` omitted: a helix
        in `\RR^3`::

            sage: R3 = Manifold(3, 'R^3')
            sage: X.<x,y,z> = R3.chart()
            sage: c = R3.curve([cos(t), sin(t), t], t, name='c') ; c
            Curve 'c' in the 3-dimensional manifold 'R^3'
            sage: c.domain() # check that t is unbounded
            field R of real numbers

        See the documentation of
        :class:`~sage.geometry.manifolds.curve.ManifoldCurve` for more
        examples.

        """
        from sage.geometry.manifolds.manifold import RealLine
        if not isinstance(param, (tuple, list)):
            param = (param, -Infinity, Infinity)
        elif len(param) != 3:
            raise TypeError("the argument 'param' must be of the type " +
                            "(t, t_min, t_max)")
        t = param[0]
        t_min = param[1]
        t_max = param[2]
        real_field = RealLine(names=(repr(t),))
        interval = real_field.open_interval(t_min, t_max)
        curve_set = Hom(interval, self)
        if not isinstance(coord_expression, dict):
            # Turn coord_expression into a dictionary:
            if chart is None:
                chart = self._def_chart
            elif chart not in self._atlas:
                raise ValueError("the {} has not been".format(chart) +
                                     " defined on the {}".format(self))
            if isinstance(coord_expression, (tuple, list)):
                coord_expression = {chart: coord_expression}
            else:
                # case self.dim()=1
                coord_expression = {chart: (coord_expression,)}
        return curve_set(coord_expression, name=name, latex_name=latex_name)

    def metric(self, name, signature=None, latex_name=None, dest_map=None):
        r"""
        Define a pseudo-Riemannian metric on ``self``.

        See :class:`~sage.geometry.manifolds.metric.Metric` for a complete
        documentation.

        INPUT:

        - ``name`` -- name given to the metric
        - ``signature`` -- (default: ``None``) signature `S` of the metric as a single
          integer: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the number of
          positive terms (resp. number of negative terms) in any diagonal writing
          of the metric components; if ``signature`` is not provided, `S` is set to
          the manifold's dimension (Riemannian signature)
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the metric; if
          none, it is formed from ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of a metric field *on* ``self``)

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.metric.Metric`
          representing the defined pseudo-Riemannian metric.

        EXAMPLE:

        Metric on a 3-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: g = M.metric('g'); g
            pseudo-Riemannian metric 'g' on the 3-dimensional manifold 'M'

        See the documentation of class
        :class:`~sage.geometry.manifolds.metric.Metric` for more examples.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.metric(name, signature=signature, latex_name=latex_name)

    def riemann_metric(self, name, latex_name=None, dest_map=None):
        r"""
        Define a Riemannian metric on ``self``.

        A Riemannian metric is a field of positive-definite symmetric bilinear
        forms on ``self``.

        See :class:`~sage.geometry.manifolds.metric.RiemannMetric` for a
        complete documentation.

        INPUT:

        - ``name`` -- name given to the metric
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the metric; if
          none, it is formed from ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of a metric field *on* ``self``)

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.metric.RiemannMetric`
          representing the defined Riemannian metric.

        EXAMPLE:

        Standard metric on the 2-sphere `S^2`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: c_spher.<th,ph> = M.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: g = M.riemann_metric('g'); g
            Riemannian metric 'g' on the 2-dimensional manifold 'S^2'
            sage: g[1,1], g[2,2] = 1, sin(th)^2
            sage: g.display()
            g = dth*dth + sin(th)^2 dph*dph
            sage: g.signature()
            2

        See the documentation of class
        :class:`~sage.geometry.manifolds.metric.RiemannMetric` for more examples.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.riemann_metric(name, latex_name=latex_name)

    def lorentz_metric(self, name, signature='positive', latex_name=None,
                       dest_map=None):
        r"""
        Define a Lorentzian metric on ``self``.

        A Lorentzian metric is a field of nondegenerate symmetric bilinear
        forms with signature `(-,+,\cdots,+)` or `(+,-,\cdots,-)`.

        See :class:`~sage.geometry.manifolds.metric.LorentzMetric` for a
        complete documentation.

        INPUT:

        - ``name`` -- name given to the metric
        - ``signature`` -- (default: 'positive') sign of the metric
          signature:

          * if set to 'positive', the signature is n-2, where n is the manifold's
            dimension, i.e. `(-,+,\cdots,+)`
          * if set to 'negative', the signature is -n+2, i.e. `(+,-,\cdots,-)`

        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the metric; if
          none, it is formed from ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of a metric field *on* ``self``)

        OUTPUT:

        - instance of :class:`~sage.geometry.manifolds.metric.LorentzMetric`
          representing the defined Lorentzian metric.

        EXAMPLE:

        Metric of Minkowski spacetime::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(4, 'M')
            sage: c_cart.<t,x,y,z> = M.chart()
            sage: g = M.lorentz_metric('g'); g
            Lorentzian metric 'g' on the 4-dimensional manifold 'M'
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
            sage: g.display()
            g = -dt*dt + dx*dx + dy*dy + dz*dz
            sage: g.signature()
            2

        See the documentation of class
        :class:`~sage.geometry.manifolds.metric.LorentzMetric` for more
        examples.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.lorentz_metric(name, signature=signature,
                                      latex_name=latex_name)

    def diff_form(self, degree, name=None, latex_name=None,
                  dest_map=None):
        r"""
        Define a differential form on ``self``.

        See :class:`~sage.geometry.manifolds.diffform.DiffForm` for a complete
        documentation.

        INPUT:

        - ``degree`` -- the degree `p` of the differential form (i.e. its
          tensor rank)
        - ``name`` -- (default: ``None``) name given to the differential form
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          differential form; if none is provided, the LaTeX symbol is set to ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of a differential form *on* ``self``)

        OUTPUT:

        - the `p`-form, as an instance of
          :class:`~sage.geometry.manifolds.diffform.DiffForm`

        EXAMPLE:

        A 2-form on a 4-dimensional open subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(4, 'M')
            sage: A = M.open_subset('A', latex_name=r'\mathcal{A}'); A
            open subset 'A' of the 4-dimensional manifold 'M'
            sage: c_xyzt.<x,y,z,t> = A.chart()
            sage: f = A.diff_form(2, 'F'); f
            2-form 'F' on the open subset 'A' of the 4-dimensional manifold 'M'

        See the documentation of class
        :class:`~sage.geometry.manifolds.diffform.DiffForm` for more examples.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.alternating_form(degree, name=name,
                                        latex_name=latex_name)

    def one_form(self, name=None, latex_name=None, dest_map=None):
        r"""
        Define a 1-form on ``self``.

        See :class:`~sage.geometry.manifolds.diffform.DiffForm` for a complete
        documentation.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the 1-form
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the 1-form;
          if none is provided, the LaTeX symbol is set to ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`
          representing the destination map `\Phi:\ U \rightarrow V`, where `U`
          is ``self``; if none is provided, the identity map is assumed (case
          of a 1-form *on* ``self``)

        OUTPUT:

        - the 1-form, as an instance of
          :class:`~sage.geometry.manifolds.diffform.DiffForm`

        EXAMPLE:

        A 1-form on a 3-dimensional open subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: A = M.open_subset('A', latex_name=r'\mathcal{A}')
            sage: X.<x,y,z> = A.chart()
            sage: om = A.one_form('omega', r'\omega') ; om
            1-form 'omega' on the open subset 'A' of the 3-dimensional manifold 'M'
            sage: om.parent()
            Free module /\^1(A) of 1-forms on the open subset 'A' of the
             3-dimensional manifold 'M'

        See the documentation of class
        :class:`~sage.geometry.manifolds.diffform.DiffForm` for more examples.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.linear_form(name=name, latex_name=latex_name)

    def diff_mapping(self, codomain, coord_functions=None, chart1=None,
                     chart2=None, name=None, latex_name=None):
        r"""
        Define a differentiable mapping between ``self`` and another
        subset (possibly on another manifold).

        See :class:`~sage.geometry.manifolds.diffmapping.DiffMapping` for a
        complete documentation.

        INPUT:

        - ``codomain`` -- mapping's codomain (the arrival manifold or some
          subset of it)
        - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
          either

          - (i) a dictionary of
            the coordinate expressions (as lists (or tuples) of the
            coordinates of the image expressed in terms of the coordinates of
            the considered point) with the pairs of charts (chart1, chart2)
            as keys (chart1 being a chart on ``self`` and chart2 a chart on
            ``codomain``)
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``

          In both cases, if the dimension of the arrival manifold is 1,
          a single coordinate expression can be passed instead of a tuple with
          a single element
        - ``chart1`` -- (default: ``None``; used only in case (ii) above) chart
          on ``self`` defining the start coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``self``
        - ``chart2`` -- (default: ``None``; used only in case (ii) above) chart
          on ``codomain`` defining the arrival coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``codomain``
        - ``name`` -- (default: ``None``) name given to the differentiable
          mapping
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          differentiable mapping; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - the differentiable mapping, as an instance of
          :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`

        EXAMPLES:

        A mapping between an open subset of `S^2` covered by regular spherical
        coordinates and `\RR^3`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'S^2')
            sage: U = M.open_subset('U')
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: N = Manifold(3, 'R^3', r'\RR^3')
            sage: c_cart.<x,y,z> = N.chart()  # Cartesian coord. on R^3
            sage: Phi = U.diff_mapping(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)),
            ....:                      name='Phi', latex_name=r'\Phi') ; Phi
            differentiable mapping 'Phi' from the open subset 'U' of the
             2-dimensional manifold 'S^2' to the 3-dimensional manifold 'R^3'

        The same definition, but with a dictionary with pairs of charts as
        keys (case (i) above)::

            sage: Phi1 = U.diff_mapping(N,
            ....:        {(c_spher, c_cart): (sin(th)*cos(ph), sin(th)*sin(ph), cos(th))},
            ....:        name='Phi', latex_name=r'\Phi')
            sage: Phi1 == Phi
            True

        The differentiable mapping acting on a point::

            sage: p = U.point((pi/2, pi)) ; p
            point on 2-dimensional manifold 'S^2'
            sage: Phi(p)
            point on 3-dimensional manifold 'R^3'
            sage: Phi(p).coord(c_cart)
            (-1, 0, 0)
            sage: Phi1(p) == Phi(p)
            True

        See the documentation of class
        :class:`~sage.geometry.manifolds.diffmapping.DiffMapping` for more
        examples.

        """
        homset = Hom(self, codomain)
        if coord_functions is None:
            coord_functions = {}
        if not isinstance(coord_functions, dict):
            # Turn coord_functions into a dictionary:
            if chart1 is None:
                chart1 = self._def_chart
            elif chart1 not in self._atlas:
                raise ValueError("{} is not a chart ".format(chart1) +
                                 "defined on the {}".format(self))
            if chart2 is None:
                chart2 = codomain._def_chart
            elif chart2 not in codomain._atlas:
                raise ValueError("{} is not a chart ".format(chart2) +
                                 " defined on the {}".format(codomain))
            coord_functions = {(chart1, chart2): coord_functions}
        return homset(coord_functions, name=name, latex_name=latex_name)

    def diffeomorphism(self, codomain, coord_functions=None, chart1=None,
                       chart2=None, name=None, latex_name=None):
        r"""
        Define a diffeomorphism between ``self`` and another open subset
        (possibly on another manifold).

        See :class:`~sage.geometry.manifolds.diffmapping.DiffMapping` for a
        complete documentation.

        INPUT:

        - ``codomain`` -- mapping's codomain (the arrival manifold or some
          subset of it)
        - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
          either

          - (i) a dictionary of
            the coordinate expressions (as lists (or tuples) of the
            coordinates of the image expressed in terms of the coordinates of
            the considered point) with the pairs of charts (chart1, chart2)
            as keys (chart1 being a chart on ``self`` and chart2 a chart on
            ``codomain``)
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``

          In both cases, if the dimension of the arrival manifold is 1,
          a single coordinate expression can be passed instead of a tuple with
          a single element
        - ``chart1`` -- (default: ``None``; used only in case (ii) above) chart
          on ``self`` defining the start coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``self``
        - ``chart2`` -- (default: ``None``; used only in case (ii) above) chart
          on ``codomain`` defining the arrival coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``codomain``
        - ``name`` -- (default: ``None``) name given to the diffeomorphism
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          diffeomorphism; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - the diffeomorphism, as an instance of
          :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`

        EXAMPLE:

        A diffeomorphism between two 2-dimensional subsets::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M', r'{\cal M}')
            sage: U = M.open_subset('U')
            sage: c_xv.<x,y> = U.chart(r'x:(-pi/2,+oo) y:(-pi/2,+oo)')
            sage: N = Manifold(2, 'N', r'{\cal N}')
            sage: V = N.open_subset('V')
            sage: c_zt.<z,t> = V.chart(r'z t')
            sage: Phi = U.diffeomorphism(V, (arctan(x), arctan(y)), name='Phi',
            ....:                        latex_name=r'\Phi')

        See the documentation of class
        :class:`~sage.geometry.manifolds.diffmapping.DiffMapping` for more
        examples.

        """
        homset = Hom(self, codomain)
        if coord_functions is None:
            coord_functions = {}
        if not isinstance(coord_functions, dict):
            # Turn coord_functions into a dictionary:
            if chart1 is None:
                chart1 = self._def_chart
            elif chart1 not in self._atlas:
                raise ValueError("{} is not a chart ".format(chart1) +
                                 "defined on the {}".format(self))
            if chart2 is None:
                chart2 = codomain._def_chart
            elif chart2 not in codomain._atlas:
                raise ValueError("{} is not a chart ".format(chart2) +
                                 " defined on the {}".format(codomain))
            coord_functions = {(chart1, chart2): coord_functions}
        return homset(coord_functions, name=name, latex_name=latex_name,
                      is_diffeomorphism=True)

    def identity_map(self):
        r"""
        Identity map on ``self``.

        See :class:`~sage.geometry.manifolds.diffmapping.DiffMapping` for a
        complete documentation.

        OUTPUT:

        - the identity map, as an instance of
          :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`

        """
        return self._identity_map

    def aff_connection(self, name, latex_name=None):
        r"""
        Define an affine connection on ``self``.

        See :class:`~sage.geometry.manifolds.connection.AffConnection` for a
        complete documentation.

        INPUT:

        - ``name`` -- name given to the affine connection
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the affine
          connection

        OUTPUT:

        - the affine connection, as an instance of
          :class:`~sage.geometry.manifolds.connection.AffConnection`

        EXAMPLE:

        Affine connection on a 3-dimensional subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M', start_index=1)
            sage: A = M.open_subset('A', latex_name=r'\mathcal{A}')
            sage: #!# nab = A.aff_connection('nabla', r'\nabla') ; nab
            affine connection 'nabla' on the subset 'A' of the 3-dimensional manifold 'M'

        See the documentation of class
        :class:`~sage.geometry.manifolds.connection.AffConnection` for more
        examples.

        """
        from connection import AffConnection
        return AffConnection(self, name, latex_name)

    def set_frame_change(self, frame1, frame2, change_of_frame,
                         compute_inverse=True):
        r"""
        Relates two vector frames by an automorphism.

        This updates the internal dictionary ``self._frame_changes``.

        INPUT:

        - ``frame1`` -- frame 1, denoted `(e_i)`  below
        - ``frame2`` -- frame 2, denoted `(f_i)`  below
        - ``change_of_frame`` -- instance of class
          :class:`~sage.geometry.manifolds.automorphismfield.AutomorphismFieldParal`
          describing the automorphism `P` that relates the basis `(e_i)` to
          the basis `(f_i)` according to `f_i = P(e_i)`
        - ``compute_inverse`` (default: True) -- if set to True, the inverse
          automorphism is computed and the change from basis `(f_i)` to `(e_i)`
          is set to it in the internal dictionary ``self._frame_changes``

        EXAMPLE:

        Connecting two vector frames on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: e = M.vector_frame('e')
            sage: f = M.vector_frame('f')
            sage: a = M.automorphism_field()
            sage: a[e,:] = [[1,2],[0,3]]
            sage: M.set_frame_change(e, f, a)
            sage: f[0].display(e)
            f_0 = e_0
            sage: f[1].display(e)
            f_1 = 2 e_0 + 3 e_1
            sage: e[0].display(f)
            e_0 = f_0
            sage: e[1].display(f)
            e_1 = -2/3 f_0 + 1/3 f_1
            sage: M.frame_change(e,f)[e,:]
            [1 2]
            [0 3]

        """
        from automorphismfield import AutomorphismFieldParal
        fmodule = frame1._fmodule
        if frame2._fmodule != fmodule:
            raise ValueError("The two frames are not defined on the same " +
                             "vector field module.")
        if not isinstance(change_of_frame, AutomorphismFieldParal):
            raise TypeError("The argument change_of_frame must be some " +
                            "instance of AutomorphismFieldParal.")
        fmodule.set_change_of_basis(frame1, frame2, change_of_frame,
                                    compute_inverse=compute_inverse)
        for sdom in self._supersets:
            sdom._frame_changes[(frame1, frame2)] = change_of_frame
        if compute_inverse:
            if (frame2, frame1) not in self._frame_changes:
                for sdom in self._supersets:
                    sdom._frame_changes[(frame2, frame1)] = \
                                                      change_of_frame.inverse()
