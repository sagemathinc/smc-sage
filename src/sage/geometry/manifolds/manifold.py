r"""
Differentiable manifolds

The class :class:`Manifold` implements differentiable manifolds over `\RR`. 

Ideally this class should inherit from a class describing topological 
manifolds or at least topological spaces. Since such classes do not
exist in Sage yet, the class :class:`Manifold` inherits from the 
class :class:`~sage.geometry.manifolds.domain.OpenDomain`. 
Via the latter, the class :class:`Manifold` inherits 
from the generic Sage class :class:`~sage.structure.parent.Parent` 
and is declared to belong to the category of sets (Sage category 
:class:`~sage.categories.sets_cat.Sets`).
The corresponding Sage :class:`~sage.structure.element.Element`'s are 
implemented via the class :class:`~sage.geometry.manifolds.point.Point`. 

The derived class :class:`RealLine` implements the field of real numbers
`\RR` as a manifold of dimension one. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014): initial version

EXAMPLES:
    
    The sphere `S^2` as a 2-dimensional manifold::
    
        sage: M = Manifold(2, 'S^2')
        sage: M
        2-dimensional manifold 'S^2'
        sage: M.dim()
        2

    Let us consider the complement of the North pole; it is an open domain
    of `S^2`, which we call U::
        
        sage: U = M.open_domain('U') ; U
        open domain 'U' on the 2-dimensional manifold 'S^2'
        
    A standard chart on U is provided by the stereographic projection from the
    North pole to the equatorial plane::
    
        sage: stereoN.<x,y> = U.chart() ; stereoN
        chart (U, (x, y))
        
    Thanks to the operator <x,y> on the left-hand side, the coordinates 
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

    Let us call V the domain that is the complement of the South pole and let
    us introduce on it the chart induced by the stereographic projection from
    the South pole to the equatorial plane::
    
        sage: V = M.open_domain('V') ; V
        open domain 'V' on the 2-dimensional manifold 'S^2'
        sage: stereoS.<u,v> = V.chart() ; stereoS
        chart (V, (u, v))

    The North pole is the point of coordinates `(u,v)=(0,0)` in this chart::
    
        sage: N = V.point((0,0), name='N') ; N
        point 'N' on 2-dimensional manifold 'S^2'

    At this stage, the manifold's atlas contains two charts::
    
        sage: M.atlas()
        [chart (U, (x, y)), chart (V, (u, v))]

    To finalize things, we must declare the transition map between these two
    charts: calling W the intersection of U and V, (W is the subdomain of U 
    defined by `x^2+y^2\not=0`, as well as the subdomain of V defined by 
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
           
    At this stage, we have four open domains on `S^2`::
        
        sage: M.domains()
        [2-dimensional manifold 'S^2',
         open domain 'U' on the 2-dimensional manifold 'S^2',
         open domain 'V' on the 2-dimensional manifold 'S^2',
         open domain 'W' on the 2-dimensional manifold 'S^2']

    W is the open domain that is the complement of the two poles::
    
        sage: N in W
        False
        sage: S in W
        False

    The North pole lies in `V\setminus U` and the South pole in `U\setminus V`::
    
        sage: N in V, N in U
        (True, False)
        sage: S in U, S in V
        (True, False)

    Four charts have been defined on the manifold::
    
        sage: M.atlas()
        [chart (U, (x, y)), chart (V, (u, v)), chart (W, (x, y)), chart (W, (u, v))]
         
    The first defined chart is considered as the default chart on the 
    manifold (unless it is changed by the method 
    :meth:`~sage.geometry.manifolds.domain.Domain.set_default_chart`)::
    
        sage: M.default_chart()
        chart (U, (x, y))

    Being the *default chart* means that its mention can be omitted when 
    specifying some point coordinates::
    
        sage: p = M.point((1,2))  # a point is created with coordinates (1,2)
        sage: p._coordinates # random (dictionary output):
        {chart (W, (x, y)): (1, 2), chart (U, (x, y)): (1, 2)}
        sage: p.coord() # if the chart is not specified, the default chart coordinates are returned:
        (1, 2)
        sage: p.coord(stereoS_W) # the coordinates in the chart stereoS_W are computed by means of the transition map:
        (1/5, 2/5)
        
    Manifolds are 'Parent' Sage objects, whose elements are the points::
    
        sage: p.parent()
        2-dimensional manifold 'S^2'
        sage: p in M
        True
        sage: p == M((1,2))
        True
        
    A manifold has a default vector frame, which, unless otherwise specified, 
    is the coordinate frame associated with the first defined chart::
    
        sage: M.default_frame()
        coordinate frame (U, (d/dx,d/dy))
        sage: latex(M.default_frame())
        \left(U ,\left(\frac{\partial}{\partial x },\frac{\partial}{\partial y }\right)\right)
        
    A manifold has a predefined zero scalar field, mapping all the points to 0; 
    it is an instance of 
    :class:`~sage.geometry.manifolds.scalarfield.ZeroScalarField`::
    
        sage: M._zero_scalar_field
        zero scalar field on the 2-dimensional manifold 'S^2'
        sage: M._zero_scalar_field(p)
        0
        sage: M._zero_scalar_field(N)
        0
        sage: M._zero_scalar_field(S)
        0
        
"""

#*****************************************************************************
#       Copyright (C) 2013, 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013, 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from domain import OpenDomain

class Manifold(OpenDomain):
    r"""  
    Base class for differentiable manifolds.
    
    This class implements differentiable manifolds over `\RR`. Ideally it 
    should inherit from a class describing topological manifolds, or at 
    least, topological spaces (not existing yet in Sage!). 
    
    INPUT:
    
    - ``n`` -- dimension of the manifold
    - ``name`` -- name given to the manifold 
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the manifold; if 
      none is provided, it is set to ``name``
    - ``start_index`` -- (default: 0) lower bound of the range of indices on the
      manifold
    
    EXAMPLES:

    A 4-dimensional manifold::
    
        sage: M = Manifold(4, 'M', latex_name=r'\mathcal{M}')
        sage: M
        4-dimensional manifold 'M'
        sage: latex(M)
        \mathcal{M}
                
    The input parameter ``start_index`` becomes the attribute :attr:`sindex`
    of the manifold::
    
        sage: M = Manifold(4, 'M')  # default value of start_index is 0
        sage: M._sindex
        0
        sage: M = Manifold(4, 'M', start_index=1)
        sage: M._sindex
        1
        
    It defines the range of indices on the manifold::
    
        sage: M = Manifold(4, 'M')
        sage: list(M.irange())
        [0, 1, 2, 3]
        sage: M = Manifold(4, 'M', start_index=2)
        sage: list(M.irange())
        [2, 3, 4, 5]

    A manifold is a Sage *Parent* object, in the category of sets::
    
        sage: M.category()
        Category of sets
        sage: M in Sets()
        True

    The corresponding Sage *Elements* are points::
    
        sage: M.chart('t x y z')
        chart (M, (t, x, y, z))
        sage: p = M.an_element(); p
        point on 4-dimensional manifold 'M'
        sage: p.parent()
        4-dimensional manifold 'M'
        sage: p in M
        True

    The manifold's points are instances of class 
    :class:`~sage.geometry.manifolds.point.Point`::
    
        sage: isinstance(p, sage.geometry.manifolds.point.Point)
        True

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
        OpenDomain.__init__(self, self, name, latex_name)
        self._sindex = start_index
        self._domains = [self]
        
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

    def dim(self):
        r"""
        Return the dimension of the manifold.
        
        EXAMPLE::
        
            sage: M = Manifold(2, 'M')
            sage: M.dim()
            2

        """
        return self._dim

    def domains(self):
        r"""
        Return the list of subdomains that have been defined on the manifold.
        
        EXAMPLE:
        
        Domains on a 2-dimensional manifold::
        
            sage: M = Manifold(2, 'M')
            sage: U = M.open_domain('U')
            sage: V = M.domain('V')
            sage: M.domains()
            [2-dimensional manifold 'M',
             open domain 'U' on the 2-dimensional manifold 'M',
             domain 'V' on the 2-dimensional manifold 'M']
            sage: U is M.domains()[1]
            True

        """
        return self._domains


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
    canonical coordinate. 
        
    INPUT: 

    - ``coordinate`` -- (default: None) string defining the symbol of the 
      canonical coordinate set on the real line; if none is provided and 
      ``names`` is none, the symbol 't' is used
    - ``name`` -- (default: 'R') name given to the real line
    - ``latex_name`` -- (default: r'\\RR') LaTeX symbol to denote the real line
    - ``start_index`` -- (default: 0) unique value of the index for vectors and
      forms on the real line. 
    - ``names`` -- (default: None) used only when ``coordinate`` is None: it 
      must be a single-element tuple containing the canonical coordinate
      symbol (this is guaranted if the shortcut operator <,> is used, see
      examples below). 
    
    EXAMPLES:
                
    Constructing the real line without any argument::
    
        sage: R = RealLine() ; R
        field R of real numbers
        sage: latex(R)
        \RR

    R is a 1-dimensional manifold::
    
        sage: isinstance(R, Manifold)
        True
        sage: R.dim()
        1
    
    It is endowed with a default chart (canonical coordinate)::
    
        sage: R.default_chart()
        chart (R, (t,))
        sage: R.atlas()
        [chart (R, (t,))]

    The instance is unique (as long as the constructor arguments are the same)::
    
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
    by means of the shortcut operator <>::
    
        sage: R.<t> = RealLine()
        sage: t
        t
        sage: type(t)
        <type 'sage.symbolic.expression.Expression'>

    The trick is performed by Sage preparser::
    
        sage: preparse("R.<t> = RealLine()")
        "R = RealLine(names=('t',)); (t,) = R._first_ngens(1)"

    In particular the <> operator is to be used to set a canonical 
    coordinate symbol different from 't'::
    
        sage: R.<u> = RealLine()
        sage: R.atlas()
        [chart (R, (u,))]
        sage: R.canonical_coordinate()
        u        
    
    The LaTeX symbol of the canonical coordinate can be adjusted via the same
    syntax as a chart declaration (see 
    :class:`~sage.geometry.manifolds.chart.Chart`)::
    
        sage: R.<x> = RealLine(coordinate=r'x:\xi')
        sage: latex(x)
        \xi
        sage: latex(R.default_chart())
        (\RR,(\xi))

    The LaTeX symbol of the real line itself can also be customized::
    
        sage: R.<x> = RealLine(latex_name=r'\mathbb{R}')
        sage: latex(R)
        \mathbb{R}
        
    """
    def __init__(self, coordinate=None, name='R', latex_name=r'\RR', 
                 start_index=0, names=None):
        from chart import Chart
        Manifold.__init__(self, 1, name, latex_name=latex_name, 
                          start_index=start_index)
        if coordinate is None:
            if names is None:
                coordinate = 't'
            else:
                coordinate = names[0]
        Chart(self, coordinates=coordinate)

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "field " + self._name + " of real numbers"

    def _first_ngens(self, n):
        r"""
        Return the coordinate of the default chart
        
        This is useful only for the use of Sage preparser:

        """
        return self._def_chart[:]

    def canonical_coordinate(self):
        r""" 
        Return the canonical coordinate defined on ``self``. 
        
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
        return self._def_chart._xx[0]
