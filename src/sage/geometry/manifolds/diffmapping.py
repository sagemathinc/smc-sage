r"""
Differentiable mappings between manifolds

The class :class:`DiffMapping` implements differentiable mappings from an open
domain `U` of a differentiable manifold `\mathcal{M}` to a differentiable
manifold `\mathcal{N}`: 

.. MATH::

    \Phi: U\subset \mathcal{M} \longrightarrow \mathcal{N}
    
In what follows, `\mathcal{M}` is called the *start manifold* and `\mathcal{N}`
the *arrival manifold*. The case `\mathcal{N}=\mathcal{M}` is allowed. 

The special case of *diffeomorphisms*, i.e. of invertible mappings such that
both `\Phi` and `\Phi^{-1}` are differentiable, is implemented through the 
class :class:`Diffeomorphism`, which inherits from :class:`DiffMapping`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014): initial version

EXAMPLES: 

    A mapping between the sphere `S^2` and `\RR^3`::

        sage: M = Manifold(2, 'S^2')
        sage: U = M.open_domain('U') # the subdomain of S^2 covered by regular spherical coordinates
        sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        sage: N = Manifold(3, 'R^3', r'\RR^3')
        sage: c_cart.<x,y,z> = N.chart('x y z')  # Cartesian coord. on R^3
        sage: Phi = U.diff_mapping(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
        sage: Phi.view()
        Phi: U --> R^3, (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
        
    The mapping acting on a point::
       
        sage: p = U.point((pi/2,pi/2), name='P')
        sage: q = Phi(p) ; q
        point 'Phi(P)' on 3-dimensional manifold 'R^3'
        sage: q.coord()
        (0, 1, 0)
        sage: (u, v) = var('u v')
        sage: a = U.point((u,v)) # a (unnamed) point defined by a symbolic expression
        sage: Phi(a)
        point on 3-dimensional manifold 'R^3'
        sage: Phi(a).coord()
        (cos(v)*sin(u), sin(u)*sin(v), cos(u))
        
    An example of diffeomorphism: a rotation in the Euclidean plane::

        sage: M = Manifold(2, 'R^2', r'\RR^2')
        sage: c_cart.<x,y> = M.chart('x y') # Cartesian coordinates
        sage: # A pi/3 rotation around the origin:
        sage: rot = M.diffeomorphism(M, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
        sage: p = M.point((1,2), name='p')
        sage: q = rot(p) ; q
        point 'R(p)' on 2-dimensional manifold 'R^2'
        sage: q.coord()
        (-sqrt(3) + 1/2, 1/2*sqrt(3) + 1)
 
    The inverse diffeormorphism::

        sage: rot.inverse() 
        diffeomorphism 'R^(-1)' on the 2-dimensional manifold 'R^2'
        sage: rot.inverse().view()
        R^(-1): R^2 --> R^2, (x, y) |--> (1/2*sqrt(3)*y + 1/2*x, -1/2*sqrt(3)*x + 1/2*y)
        sage: p1 = rot.inverse()(q)
        sage: p1
        point 'R^(-1)(R(p))' on 2-dimensional manifold 'R^2'
        sage: p1 == p 
        True
        
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

from sage.structure.sage_object import SageObject
from domain import Domain
from chart import Chart, FunctionChart, MultiFunctionChart, CoordChange
from point import Point
     
class DiffMapping(SageObject):
    r"""
    Class for differentiable mappings between manifolds.

    This class implements differentiable mappings of the type
    
    .. MATH::
    
        \Phi: U\subset \mathcal{M} \longrightarrow \mathcal{N}
    
    where  `U` is a open subset of some differentiable manifold `\mathcal{M}` 
    and `\mathcal{N}` is a differentiable manifold.
    In what follows, `\mathcal{M}` is called the *start manifold* and 
    `\mathcal{N}` the *arrival manifold*. 

    INPUT:
    
    - ``domain`` -- mapping's domain `U` (open subset of the start 
      manifold)
    - ``codomain`` -- mapping's codomain (the arrival manifold or some subset
      of it)
    - ``coord_functions`` -- (default: None) the coordinate symbolic expression 
      of the mapping in some pair of charts: list (or tuple) of the 
      coordinates of the image expressed in terms of the coordinates of 
      the considered point; if the dimension of the arrival manifold is 1, 
      a single expression is expected (not a list with a single element)
    - ``chart1`` -- (default: None) chart on domain `U` in which the 
      coordinates are given for ``coord_functions``; if none is provided, the 
      coordinates are assumed to refer to domain's default chart
    - ``chart2`` -- (default: None) chart on the codomain for the coordinates
      on the arrival manifold for ``coord_functions``; if none is provided, the 
      coordinates are assumed to refer to the codomain's default chart
    - ``name`` -- (default: None) name given to the differentiable mapping
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
      differentiable mapping; if none is provided, the LaTeX symbol is set to 
      ``name``
    
    .. NOTE:
    
        If ``chart1`` does not cover the entire domain `U`, the argument
        ``coord_functions` is not not sufficient to fully specify the 
        differential mapping; further coordinate expressions, in other charts,
        can be subsequently added by means of the method :meth:`add_expr`
    
    EXAMPLES:
    
    A mapping between the sphere `S^2` and `\RR^3`::

        sage: M = Manifold(2, 'S^2')
        sage: U = M.open_domain('U') # the subdomain of S^2 covered by regular spherical coordinates
        sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
        sage: N = Manifold(3, 'R^3', r'\RR^3')
        sage: c_cart.<x,y,z> = N.chart('x y z')  # Cartesian coord. on R^3
        sage: Phi = U.diff_mapping(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
        sage: Phi
        differentiable mapping 'Phi' from open domain 'U' on the 2-dimensional manifold 'S^2' to 3-dimensional manifold 'R^3'
        sage: type(Phi)
        <class 'sage.geometry.manifolds.diffmapping.DiffMapping'>
        sage: Phi.view()
        Phi: U --> R^3, (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
        
    The mapping acting on a point::
       
        sage: p = U.point((pi/2,pi/2), name='P')
        sage: q = Phi(p) ; q
        point 'Phi(P)' on 3-dimensional manifold 'R^3'
        sage: q.coord()
        (0, 1, 0)
        sage: (u, v) = var('u v')
        sage: a = U.point((u,v)) # a point defined by a symbolic expression
        sage: Phi(a).coord()
        (cos(v)*sin(u), sin(u)*sin(v), cos(u))

    If the arrival manifold is 1-dimensional, the mapping is defined by a
    single symbolic expression and not a list with a single element::

        sage: N = Manifold(1, 'N')
        sage: chart_n = N.chart('x')
        sage: Phi = M.diff_mapping(N, sin(th)*cos(ph)) # and not ...,(sin(th)*cos(ph),))
        
    If the arrival manifold is the field of real numbers `\RR` (the Sage object
    :data:`RealLine`), the action on a point returns a real number, i.e. the 
    canonical coordinate of the image point, and not the image point itself::

        sage: Phi = M.diff_mapping(RealLine, sin(th)*cos(ph))
        sage: p = U.point((pi/2,pi))
        sage: Phi(p)      
        -1

    """
    def __init__(self, domain, codomain, coord_functions=None, chart1=None, 
                 chart2=None, name=None, latex_name=None): 
        if not isinstance(domain, Domain):
            raise TypeError("The argument domain must be a domain.")
        if not isinstance(codomain, Domain):
            raise TypeError("The argument codomain must be a domain.")
        self.domain = domain
        self.codomain = codomain
        if coord_functions is not None:
            if chart1 is None: chart1 = domain.def_chart
            if chart2 is None: chart2 = codomain.def_chart
            if chart1 not in self.domain.atlas:
                raise ValueError("The " + str(chart1) +
                                    " has not been defined on the " + 
                                    str(self.domain))
            if chart2 not in self.codomain.atlas:
                raise ValueError("The " + str(chart2) +
                                    " has not been defined on the " + 
                                    str(self.codomain))
            n2 = self.codomain.manifold.dim
            if n2 > 1:
                if len(coord_functions) != n2:
                    raise ValueError(str(n2) + 
                                     " coordinate function must be provided.")
                self.coord_expression = {(chart1, chart2): 
                                        MultiFunctionChart(chart1, *coord_functions)}
            else:
                self.coord_expression = {(chart1, chart2): 
                                        MultiFunctionChart(chart1, coord_functions)}
        else: # case coord_functions is None:
            self.coord_expression = {}
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        # Initialization of derived quantities:
        DiffMapping._init_derived(self)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "differentiable mapping"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " from " + str(self.domain) + " to " + \
                       str(self.codomain)
        return description
        
    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        if self.latex_name is None:
            return r'\mbox{no symbol}'
        else:
           return self.latex_name

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        pass # no derived quantity yet

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        pass # no derived quantity yet


    def view(self, chart1=None, chart2=None):
        r""" 
        Display the expression of the differentiable mapping in a given 
        pair of charts. 
        
        If the expression is not known already, it is computed from some
        expression in other charts by means of change-of-chart formulas.
        
        INPUT:
        
        - ``chart1`` -- (default: None) chart on the mapping's domain; if None, 
          the domain's default chart will be used
        - ``chart2`` -- (default: None) chart on the mapping's codomain; if 
          None, the codomain's default chart will be used
          
        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode). 
        
        EXAMPLES:
        
        Standard embedding of the sphere `S^2` in `\RR^3`::
    
            sage: M = Manifold(2, 'S^2')
            sage: c_spher.<th,ph> = M.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: N = Manifold(3, 'R^3', r'\RR^3')
            sage: c_cart.<x,y,z> = N.chart('x y z')
            sage: Phi = M.diff_mapping(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
            sage: Phi.view()
            Phi: S^2 --> R^3, (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: latex(Phi.view())
            \begin{array}{llcl} \Phi:& S^2 & \longrightarrow & \RR^3 \\ & \left(\theta, \phi\right) & \longmapsto & \left(x, y, z\right) = \left(\cos\left(\phi\right) \sin\left(\theta\right), \sin\left(\phi\right) \sin\left(\theta\right), \cos\left(\theta\right)\right) \end{array}

        """
        from sage.misc.latex import latex
        from utilities import FormattedExpansion
        result = FormattedExpansion(self)
        if chart1 is None:
            chart1 = self.domain.def_chart
        if chart2 is None:
            chart2 = self.codomain.def_chart
        expression = self.expr(chart1, chart2)
        if self.name is None:
            symbol = ""
        else:
            symbol = self.name + ": "
        result.txt = symbol + self.domain.name + " --> " + \
                     self.codomain.name + ", " + repr(chart1[:]) + " |--> " 
        if chart2 == chart1:
            result.txt += repr(expression)
        else:
            result.txt += repr(chart2[:]) + " = " + repr(expression)
        if self.latex_name is None:
            symbol = ""
        else:
            symbol = self.latex_name + ":"
        result.latex = r"\begin{array}{llcl} " + symbol + r"&" + \
                       latex(self.domain) + r"& \longrightarrow & " + \
                       latex(self.codomain) + r"\\ &" + latex(chart1[:]) + \
                       r"& \longmapsto & " 
        if chart2 == chart1:
            result.latex += latex(expression) + r"\end{array}"
        else:
            result.latex += latex(chart2[:]) + " = " + latex(expression) + \
                            r"\end{array}"
        return result


    def multi_function_chart(self, chart1=None, chart2=None):
        r""" 
        Return the functions of the coordinates representing the differentiable
        mapping in a given pair of charts.
        
        If these functions are not already known, they are computed from known 
        ones by means of change-of-chart formulas. 
        
        INPUT:
        
        - ``chart1`` -- (default: None) chart on the mapping's domain; if None, 
          the domain's default chart is assumed
        - ``chart2`` -- (default: None) chart on the mapping's codomain; if 
          None,  the codomain's default chart is assumed

        OUTPUT:
        
        - instance of class 
          :class:`~sage.geometry.manifolds.chart.MultiFunctionChart` 
          representing the differentiable mapping in the above two charts

        EXAMPLES:

        Differential mapping from a 2-dimensional manifold to a 3-dimensional 
        one::
        
            sage: M = Manifold(2, 'M')
            sage: N = Manifold(3, 'N')
            sage: c_uv.<u,v> = M.chart('u v')
            sage: c_xyz.<x,y,z> = N.chart('x y z')
            sage: Phi = M.diff_mapping(N, (u*v, u/v, u+v), name='Phi', latex_name=r'\Phi')
            sage: Phi.view()
            Phi: M --> N, (u, v) |--> (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.multi_function_chart(c_uv, c_xyz)
            functions (u*v, u/v, u + v) on the chart (M, (u, v))
            sage: Phi.multi_function_chart() # equivalent to above since 'uv' and 'xyz' are default charts
            functions (u*v, u/v, u + v) on the chart (M, (u, v))
            sage: type(Phi.multi_function_chart())
            <class 'sage.geometry.manifolds.chart.MultiFunctionChart'>

        Representation in other charts::
        
            sage: c_UV.<U,V> = M.chart('U V')  # new chart on M
            sage: ch_uv_UV = c_uv.coord_change(c_UV, u-v, u+v)
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ.<X,Y,Z> = N.chart('X Y Z') # new chart on N
            sage: ch_xyz_XYZ = c_xyz.coord_change(c_XYZ, 2*x-3*y+z, y+z-x, -x+2*y-z)
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.multi_function_chart(c_UV, c_xyz)
            functions (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V) on the chart (M, (U, V))
            sage: Phi.multi_function_chart(c_uv, c_XYZ)
            functions (((2*u + 1)*v^2 + u*v - 3*u)/v, -((u - 1)*v^2 - u*v - u)/v, -((u + 1)*v^2 + u*v - 2*u)/v) on the chart (M, (u, v))
            sage: Phi.multi_function_chart(c_UV, c_XYZ)
            functions (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V), 1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V), 1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V)) on the chart (M, (U, V))

        """
        dom1 = self.domain; dom2 = self.codomain
        def_chart1 = dom1.def_chart; def_chart2 = dom2.def_chart
        if chart1 is None:
            chart1 = def_chart1
        if chart2 is None:
            chart2 = def_chart2
        if (chart1, chart2) not in self.coord_expression:
            # some change of coordinates must be performed
            change_start = [] ; change_arrival = []
            for (ochart1, ochart2) in self.coord_expression:
                if chart1 == ochart1:
                    change_arrival.append(ochart2)
                if chart2 == ochart2:
                    change_start.append(ochart1)
            # 1/ Trying to make a change of chart only on the arrival domain:
            # the arrival default chart is privileged:
            sel_chart2 = None # selected chart2
            if def_chart2 in change_arrival \
                    and (def_chart2, chart2) in dom2.coord_changes:
                sel_chart2 = def_chart2
            else:
                for ochart2 in change_arrival:
                    if (ochart2, chart2) in dom2.coord_changes:
                        sel_chart2 = ochart2
                        break 
            if sel_chart2 is not None:
                oexpr = self.coord_expression[(chart1, sel_chart2)]
                chg2 = dom2.coord_changes[(sel_chart2, chart2)]
                self.coord_expression[(chart1, chart2)] = \
                    MultiFunctionChart(chart1, *(chg2(*(oexpr.expr()))) )
                return self.coord_expression[(chart1, chart2)]

            # 2/ Trying to make a change of chart only on the start domain:
            # the start default chart is privileged:
            sel_chart1 = None # selected chart1
            if def_chart1 in change_start \
                    and (chart1, def_chart1) in dom1.coord_changes:
                sel_chart1 = def_chart1
            else:
                for ochart1 in change_start:
                    if (chart1, ochart1) in dom1.coord_changes:
                        sel_chart1 = ochart1
                        break
            if sel_chart1 is not None:
                oexpr = self.coord_expression[(sel_chart1, chart2)]
                chg1 = dom1.coord_changes[(chart1, sel_chart1)]
                self.coord_expression[(chart1, chart2)] = \
                    MultiFunctionChart(chart1, 
                                       *(oexpr( *(chg1.transf.expr()) )) )
                return self.coord_expression[(chart1, chart2)]
                    
            # 3/ If this point is reached, it is necessary to perform some 
            # coordinate change both on the start domain and the arrival one
            # the default charts are privileged:
            if (def_chart1, def_chart2) in self.coord_expression \
                    and (chart1, def_chart1) in dom1.coord_changes \
                    and (def_chart2, chart2) in dom2.coord_changes:
                sel_chart1 = def_chart1
                sel_chart2 = def_chart2
            else:
                for (ochart1, ochart2) in self.coord_expression:
                    if (chart1, ochart1) in dom1.coord_changes \
                        and (ochart2, chart2) in dom2.coord_changes:
                        sel_chart1 = ochart1
                        sel_chart2 = ochart2
                        break
            if (sel_chart1 is not None) and (sel_chart2 is not None):
                oexpr = self.coord_expression[(sel_chart1, sel_chart2)]
                chg1 = dom1.coord_changes[(chart1, sel_chart1)]
                chg2 = dom2.coord_changes[(sel_chart2, chart2)]
                self.coord_expression[(chart1, chart2)] = \
                     MultiFunctionChart(chart1, 
                                *(chg2( *(oexpr(*(chg1.transf.expr()))) )) )
                return self.coord_expression[(chart1, chart2)]
                
            # 4/ If this point is reached, the demanded value cannot be
            # computed 
            raise ValueError("The expression of the mapping in the pair of " +
                "charts (" + str(chart1) + ", " + str(chart2) + ") cannot " + 
                "be computed by means of known changes of charts.")
                
        return self.coord_expression[(chart1, chart2)]
            

    def expr(self, chart1=None, chart2=None):
        r""" 
        Return the expression of the differentiable mapping in terms of
        specified coordinates.
        
        If the expression is not already known, it is computed from some known 
        expression by means of change-of-chart formulas. 
        
        INPUT:
        
        - ``chart1`` -- (default: None) chart on the mapping's domain; if None, 
          the domain's default chart is assumed
        - ``chart2`` -- (default: None) chart on the mapping's codomain; if 
          None, the codomain's default chart is assumed

        OUTPUT:
        
        - symbolic expression representing the differentiable mapping in the 
          above two charts

        EXAMPLES:
        
        Differential mapping from a 2-dimensional manifold to a 3-dimensional 
        one::
                
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: N = Manifold(3, 'N')
            sage: c_uv.<u,v> = M.chart('u v')
            sage: c_xyz.<x,y,z> = N.chart('x y z')
            sage: Phi = M.diff_mapping(N, (u*v, u/v, u+v), name='Phi', latex_name=r'\Phi')
            sage: Phi.view()
            Phi: M --> N, (u, v) |--> (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.expr(c_uv, c_xyz)
            (u*v, u/v, u + v)
            sage: Phi.expr()  # equivalent to above since 'uv' and 'xyz' are default charts
            (u*v, u/v, u + v)
            sage: type(Phi.expr()[0])
            <type 'sage.symbolic.expression.Expression'>

        Expressions in other charts::
        
            sage: c_UV.<U,V> = M.chart('U V')  # new chart on M
            sage: ch_uv_UV = c_uv.coord_change(c_UV, u-v, u+v)
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ.<X,Y,Z> = N.chart('X Y Z') # new chart on N
            sage: ch_xyz_XYZ = c_xyz.coord_change(c_XYZ, 2*x-3*y+z, y+z-x, -x+2*y-z)
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.expr(c_UV, c_xyz)
            (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V)
            sage: Phi.expr(c_uv, c_XYZ)
            (((2*u + 1)*v^2 + u*v - 3*u)/v,
             -((u - 1)*v^2 - u*v - u)/v,
             -((u + 1)*v^2 + u*v - 2*u)/v)
            sage: Phi.expr(c_UV, c_XYZ)
             (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V), 1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V), 1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V))

        A rotation in some Euclidean plane::
        
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M') # the plane
            sage: c_spher.<r,ph> = M.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates on the plane
            sage: rot = M.diff_mapping(M, (r, ph+pi/3), name='R') # pi/3 rotation around r=0
            sage: rot.expr()
            (r, 1/3*pi + ph)

        Expression of the rotation in terms of Cartesian coordinates::
        
            sage: c_cart.<x,y> = M.chart('x y') # Declaration of Cartesian coordinates
            sage: ch_spher_cart = c_spher.coord_change(c_cart, r*cos(ph), r*sin(ph)) # relation to spherical coordinates
            sage: ch_spher_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x))              
            Check of the inverse coordinate transformation:
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
               x == x
               y == y
            sage: rot.expr(c_cart, c_cart)                            
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        """
        return self.multi_function_chart(chart1, chart2).expr()
            
    def set_expr(self, chart1, chart2, coord_functions): 
        r"""
        Set a new coordinate representation of the mapping.

        The expressions with respect to other charts are deleted, in order to 
        avoid any inconsistency. To keep them, use :meth:`add_expr` instead.

        INPUT:
    
        - ``chart1`` -- chart for the coordinates on the mapping's domain
        - ``chart2`` -- chart for the coordinates on the mapping's codomain
        - ``coord_functions`` -- the coordinate symbolic expression of the 
          mapping in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single 
          expression is expected (not a list with a single element)
        
        EXAMPLES:
        
        Polar representation of a planar rotation initally defined in 
        Cartesian coordinates::
            
            sage: M = Manifold(2, 'R^2', r'\RR^2')   # Euclidean plane
            sage: c_cart.<x,y> = M.chart('x y')  # Cartesian coordinates
            sage: c_spher.<r,ph> = M.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates
            sage: # Links between spherical coordinates and Cartesian ones:
            sage: ch_cart_spher = c_cart.coord_change(c_spher, sqrt(x*x+y*y), atan2(y,x))
            sage: ch_cart_spher.set_inverse(r*cos(ph), r*sin(ph))
            Check of the inverse coordinate transformation:
               x == x
               y == y
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
            sage: # A pi/3 rotation around the origin defined in terms of Cartesian coordinates:
            sage: rot = M.diff_mapping(M, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
            sage: rot.expr()
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
        
        Let us use the method :meth:`set_expr` to set the 
        spherical-coordinate expression by hand::

            sage: rot.set_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.expr(c_spher, c_spher) 
            (r, 1/3*pi + ph)
        
        The expression in Cartesian coordinates has been lost in the dictionary :attr:`coord_expression` that stores the various representations of 
        the differentiable mapping::
             
            sage: rot.coord_expression
            {(chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph))}
            
        It is recovered by a call to :meth:`expr`::
        
            sage: rot.expr(c_cart, c_cart)
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
            sage: rot.coord_expression
            {(chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph)), (chart (R^2, (x, y)), chart (R^2, (x, y))): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart (R^2, (x, y))}
            
        The rotation can be applied to a point by means of either coordinate 
        system::
            
            sage: p = M.point((1,2))  #  p defined by its Cartesian coord.
            sage: q = rot(p)  # q is computed by means of Cartesian coord.
            sage: p.coord(c_spher) # the spherical coord. of p are evaluated
            (sqrt(5), arctan(2))
            sage: q1 = rot(p, c_spher, c_spher) # q1 is computed by means of spherical coord.
            sage: q.coord(c_spher) ; # the spherical coord. of q are evaluated
            (sqrt(5), pi - arctan((sqrt(3) + 2)/(2*sqrt(3) - 1)))
            sage: q1 == q
            True
                
        """
        if chart1 not in self.domain.atlas:
            raise ValueError("The " + str(chart1) +
               " has not been defined on the " + str(self.domain))
        if chart2 not in self.codomain.atlas:
            raise ValueError("The " + str(chart2) +
              " has not been defined on the " + str(self.codomain))
        self.coord_expression.clear()
        self._del_derived()
        n2 = self.codomain.manifold.dim
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError(str(n2) + 
                                 " coordinate function must be provided.")
            self.coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
        else:
            self.coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, coord_functions)

    def add_expr(self, chart1, chart2, coord_functions): 
        r"""
        Set a new coordinate representation of the mapping.

        The previous expressions with respect to other charts are kept. To
        clear them, use :meth:`set_expr` instead. 

        INPUT:
    
        - ``chart1`` -- chart for the coordinates on the mapping's domain
        - ``chart2`` -- chart for the coordinates on the mapping's codomain
        - ``coord_functions`` -- the coordinate symbolic expression of the 
          mapping in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single 
          expression is expected (not a list with a single element)
          
        .. WARNING::
        
            If the mapping has already expressions in other charts, it 
            is the user's responsability to make sure that the expression
            to be added is consistent with them.         
        
        EXAMPLES:
        
        Polar representation of a planar rotation initally defined in 
        Cartesian coordinates::
            
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'R^2', r'\RR^2')   # Euclidean plane
            sage: c_cart.<x,y> = M.chart('x y')  # Cartesian coordinates
            sage: c_spher.<r,ph> = M.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates
            sage: # Links between spherical coordinates and Cartesian ones:
            sage: ch_cart_spher = c_cart.coord_change(c_spher, sqrt(x*x+y*y), atan2(y,x))
            sage: ch_cart_spher.set_inverse(r*cos(ph), r*sin(ph))
            Check of the inverse coordinate transformation:
               x == x
               y == y
               r == r
               ph == arctan2(r*sin(ph), r*cos(ph))
            sage: # A pi/3 rotation around the origin defined in terms of Cartesian coordinates:
            sage: rot = M.diff_mapping(M, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
            sage: rot.expr()
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        If we try to make Sage calculate the expression in terms of spherical
        coordinates, via the method :meth:`expr`, we notice some difficulties
        in arctan2 simplifications::
        
            sage: rot.expr(c_spher, c_spher) # correct output but could be simplified !
            (r,
             arctan2(1/2*(sqrt(3)*cos(ph) + sin(ph))*r, -1/2*(sqrt(3)*sin(ph) - cos(ph))*r))
        
        Therefore, we use the method :meth:`add_expr` to set the 
        spherical-coordinate expression by hand::

            sage: rot.add_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.expr(c_spher, c_spher)  # the output is now satisfactory
            (r, 1/3*pi + ph)
        
        The expression in Cartesian coordinates has been kept in the 
        dictionary :attr:`coord_expression` 
        that stores the various representations of the differentiable mapping::
             
            sage: rot.coord_expression
            {(chart (R^2, (x, y)), chart (R^2, (x, y))): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart (R^2, (x, y)), (chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph))}
             
        If, on the contrary, we use :meth:`set_expr`, the expression in 
        Cartesian coordinates is lost::
        
            sage: rot.set_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.coord_expression
            {(chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph))}
 
        It is recovered by a call to :meth:`expr`::
        
            sage: rot.expr(c_cart,c_cart)
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
            sage: rot.coord_expression
            {(chart (R^2, (r, ph)), chart (R^2, (r, ph))): functions (r, 1/3*pi + ph) on the chart (R^2, (r, ph)), (chart (R^2, (x, y)), chart (R^2, (x, y))): functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y) on the chart (R^2, (x, y))}

        The rotation can be applied to a point by means of either coordinate 
        system::
            
            sage: p = M.point((1,2))  #  p defined by its Cartesian coord.
            sage: q = rot(p)  # q is computed by means of Cartesian coord.
            sage: p.coord(c_spher) # the spherical coord. of p are evaluated
            (sqrt(5), arctan(2))
            sage: q1 = rot(p, c_spher, c_spher) # q1 is computed by means of spherical coord.
            sage: q.coord(c_spher) ; # the spherical coord. of q are evaluated
            (sqrt(5), pi - arctan((sqrt(3) + 2)/(2*sqrt(3) - 1)))
            sage: q1 == q
            True
                
        """
        if chart1 not in self.domain.atlas:
            raise ValueError("The " + str(chart1) +
               " has not been defined on the " + str(self.domain))
        if chart2 not in self.codomain.atlas:
            raise ValueError("The " + str(chart2) +
              " has not been defined on the " + str(self.codomain))
        n2 = self.codomain.manifold.dim
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError(str(n2) + 
                                 " coordinate function must be provided.")
            self.coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
        else:
            self.coord_expression[(chart1, chart2)] = \
                                   MultiFunctionChart(chart1, coord_functions)


    def __call__(self, p, chart1=None, chart2=None):
        r"""
        Compute the image of a point.

        INPUT:
    
        - ``p`` -- point on the mapping's domain (type: 
          :class:`~sage.geometry.manifolds.point.Point`)
        - ``chart1`` -- (default: None) chart in which the coordinates of p 
          are to be considered; if none is provided, a chart in which both p's 
          coordinates and the expression of ``self`` are known is searched, 
          starting from the default chart of self.domain will be used
        - ``chart2`` -- (default: None) chart in which the coordinates of the 
          image of p will be computed; if none is provided, the default chart 
          of self.codomain is assumed.
        
        OUTPUT:

        - image of the point by the mapping (type: 
          :class:`~sage.geometry.manifolds.point.Point`)

        EXAMPLES:
        
        Planar rotation acting on a point::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'R^2', r'\RR^2') # Euclidean plane
            sage: c_cart.<x,y> = M.chart('x y') # Cartesian coordinates
            sage: # A pi/3 rotation around the origin defined in Cartesian coordinates:
            sage: rot = M.diff_mapping(M, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
            sage: p = M.point((1,2), name='p')
            sage: q = rot(p) ; q
            point 'R(p)' on 2-dimensional manifold 'R^2'
            sage: q.coord()
            (-sqrt(3) + 1/2, 1/2*sqrt(3) + 1)
            
        Image computed by means of coordinates different from the default 
        ones::
        
            sage: # Spherical coord. on the plane:
            sage: c_spher.<r,ph> = M.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: ch = c_cart.coord_change(c_spher, sqrt(x*x+y*y), atan2(y,x))
            sage: rot.add_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: p.coord(c_spher) # the spherical coord. of p are evaluated
            (sqrt(5), arctan(2))
            sage: q1 = rot(p, c_spher, c_spher) # q1 is computed by means of spherical coord.
            sage: q.coord(c_spher) ; # the spherical coord. of q are evaluated
            (sqrt(5), pi - arctan((sqrt(3) + 2)/(2*sqrt(3) - 1)))
            sage: q1 == q
            True
    
        """
        from manifold import RealLine
        if p not in self.domain.manifold: 
            raise ValueError("The point " + str(p) +
                  " does not belong to the " + str(self.domain.manifold))
        if chart2 is None: 
            chart2 = self.codomain.def_chart
        if chart1 is None: 
            def_chart1 = self.domain.def_chart
            if def_chart1 in p.coordinates and \
                        (def_chart1, chart2) in self.coord_expression:
                chart1 = def_chart1
            else:
                for chart in p.coordinates:
                    if (chart, chart2) in self.coord_expression:
                        chart1 = chart
                        break
        if chart1 is None:
            raise ValueError("No common chart has been found to evaluate " \
                "the action of " + str(self) + " on the " + str(p) + ".")

        coord_map = self.coord_expression[(chart1, chart2)]
        y = coord_map(*(p.coordinates[chart1])) 
        
        if self.codomain.manifold is RealLine:   # special case of a mapping to R
            return y[0]
        else:
            if p.name is None or self.name is None:
                res_name = None
            else:
                res_name = self.name + '(' + p.name + ')'
            if p.latex_name is None or self.latex_name is None:
                res_latex_name = None
            else:
                res_latex_name = self.latex_name + r'\left(' + p.latex_name + \
                                 r'\right)'
            
            return Point(self.codomain.manifold, y, chart2, name=res_name, 
                         latex_name=res_latex_name)  #!# check
        
    def pullback(self, tensor):
        r""" 
        Pullback operator associated with the differentiable mapping. 
        
        INPUT:
        
        - ``tensor`` -- instance of class 
          :class:`~sage.geometry.manifolds.tensorfield.TensorField` 
          representing a fully covariant tensor field `T` on the mapping's
          codomain, i.e. a tensor field of type (0,p), with p a positive or 
          zero integer. The case p=0 corresponds to a scalar field.
          
        OUTPUT:
        
        - instance of class
          :class:`~sage.geometry.manifolds.tensorfield.TensorField` 
          representing a fully covariant tensor field on the mapping's domain 
          that is the pullback of `T` given by ``self``. 
          
        EXAMPLES:
        
        Pullback on `S^2` of a scalar field defined on `R^3`::
        
            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: c_spher.<th,ph> = M.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coord. on S^2
            sage: N = Manifold(3, 'R^3', r'\RR^3', start_index=1)
            sage: c_cart.<x,y,z> = N.chart('x y z') # Cartesian coord. on R^3
            sage: Phi = M.diff_mapping(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)), name='Phi', latex_name=r'\Phi')
            sage: f = N.scalar_field(x*y*z, name='f') ; f
            scalar field 'f' on the 3-dimensional manifold 'R^3'
            sage: f.view()
            f: (x, y, z) |--> x*y*z
            sage: pf = Phi.pullback(f) ; pf
            scalar field 'Phi_*(f)' on the 2-dimensional manifold 'S^2'
            sage: pf.view()
            Phi_*(f): (th, ph) |--> cos(ph)*cos(th)*sin(ph)*sin(th)^2
            
        Pullback on `S^2` of the standard Euclidean metric on `R^3`::
                
            sage: g = N.sym_bilin_form_field('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: g.view()
            g = dx*dx + dy*dy + dz*dz
            sage: pg = Phi.pullback(g) ; pg
            field of symmetric bilinear forms 'Phi_*(g)' on the 2-dimensional manifold 'S^2'
            sage: pg.view()
            Phi_*(g) = dth*dth + sin(th)^2 dph*dph

        Pullback on `S^2` of a 3-form on `R^3`::
                
            sage: a = N.diff_form(3, 'A')
            sage: a[1,2,3] = f 
            sage: a.view()
            A = x*y*z dx/\dy/\dz
            sage: pa = Phi.pullback(a) ; pa
            3-form 'Phi_*(A)' on the 2-dimensional manifold 'S^2'
            sage: pa.view() # should be zero (as any 3-form on a 2-dimensional manifold)
            Phi_*(A) = 0

        """
        from scalarfield import ScalarField
        from vectorframe import CoordFrame
        from diffform import OneFormParal, DiffFormParal
        from rank2field import SymBilinFormFieldParal
        from tensorfield import TensorFieldParal
        from sage.tensor.modules.comp import Components, CompWithSym, \
                                                 CompFullySym, CompFullyAntiSym

        #!# if not isinstance(tensor, TensorField):
        #    raise TypeError("The argument 'tensor' must be a tensor field.")
        dom1 = self.domain
        dom2 = self.codomain
        if not tensor.domain.is_subdomain(dom2):
            raise TypeError("The tensor field is not defined on the mapping " +
                            "arrival domain.")
        (ncon, ncov) = tensor.tensor_type
        if ncon != 0:
            raise TypeError("The pullback cannot be taken on a tensor " + 
                            "with some contravariant part.")
        resu_name = None ; resu_latex_name = None
        if self.name is not None and tensor.name is not None:
            resu_name = self.name + '_*(' + tensor.name + ')'
        if self.latex_name is not None and tensor.latex_name is not None:
            resu_latex_name = self.latex_name + '_*' + tensor.latex_name                
        if ncov == 0:
            # Case of a scalar field
            # ----------------------
            resu = ScalarField(dom1, name=resu_name, 
                               latex_name=resu_latex_name)
            for chart2 in tensor.express:
                for chart1 in dom1.atlas:
                    if (chart1, chart2) in self.coord_expression:
                        phi = self.coord_expression[(chart1, chart2)]
                        coord1 = chart1.xx
                        ff = tensor.express[chart2]
                        resu.add_expr( ff(*(phi(*coord1))), chart1)
            return resu
        else:
            # Case of tensor field of rank >= 1
            # ---------------------------------
            #!# What follows is valid only if dom2 and dom1 are parallelizable:
            if not dom1.is_manifestly_parallelizable():
                raise NotImplementedError("Pullback to non-parallelizable " + 
                                          "domains not implemented yet")
            fmodule1 = dom1.vector_field_module()
            ring1 = fmodule1.ring
            si1 = fmodule1.sindex
            of1 = fmodule1.output_formatter
            si2 = dom2.manifold.sindex
            if isinstance(tensor, OneFormParal):
                resu = OneFormParal(fmodule1, name=resu_name, 
                                                    latex_name=resu_latex_name)
            elif isinstance(tensor, DiffFormParal):
                resu = DiffFormParal(fmodule1, ncov, name=resu_name, 
                                                    latex_name=resu_latex_name)
            elif isinstance(tensor, SymBilinFormFieldParal):
                resu = SymBilinFormFieldParal(fmodule1, name=resu_name, 
                                                    latex_name=resu_latex_name)                
            else:
                resu = TensorFieldParal(fmodule1, (0,ncov), name=resu_name, 
                                    latex_name=resu_latex_name, sym=tensor.sym,
                                    antisym=tensor.antisym)
            for frame2 in tensor.components:
                if isinstance(frame2, CoordFrame):
                    chart2 = frame2.chart
                    for chart1 in dom1.atlas:
                        if (chart1, chart2) in self.coord_expression:
                            # Computation at the component level:
                            frame1 = chart1.frame
                            tcomp = tensor.components[frame2]
                            if isinstance(tcomp, CompFullySym):
                                ptcomp = CompFullySym(ring1, frame1, ncov, 
                                                      start_index=si1, 
                                                      output_formatter=of1)
                            elif isinstance(tcomp, CompFullyAntiSym):
                                ptcomp = CompFullyAntiSym(ring1, frame1, ncov,
                                                          start_index=si1, 
                                                          output_formatter=of1)
                            elif isinstance(tcomp, CompWithSym):
                                ptcomp = CompWithSym(ring1, frame1, ncov, 
                                                     start_index=si1, 
                                                     output_formatter=of1,
                                                     sym=tcomp.sym, 
                                                     antisym=tcomp.antisym)
                            else:
                                ptcomp = Components(ring1, frame1, ncov,
                                                    start_index=si1, 
                                                    output_formatter=of1)
                            phi = self.coord_expression[(chart1, chart2)]
                            jacob = phi.jacobian()
                            # X2 coordinates expressed in terms of X1 ones via the mapping:
                            coord2_1 = phi(*(chart1.xx)) 
                            for ind_new in ptcomp.non_redundant_index_generator(): 
                                res = 0 
                                for ind_old in dom2.manifold.index_generator(ncov): 
                                    ff = tcomp[[ind_old]].function_chart(chart2)
                                    t = FunctionChart(chart1, ff(*coord2_1))
                                    for i in range(ncov):
                                        t *= jacob[ind_old[i]-si2][ind_new[i]-si1]
                                    res += t
                                ptcomp[ind_new] = res
                            resu.components[frame1] = ptcomp
            return resu

        
#*****************************************************************************

class Diffeomorphism(DiffMapping):
    r"""
    Class for manifold diffeomorphisms.

    A *diffeomorphism* is a differential mapping 
    
    .. MATH::
    
        \Phi: U\subset \mathcal{M} \longrightarrow \Phi(U)\subset\mathcal{N}
        
    where  `U` is a open subset of some differentiable manifold `\mathcal{M}` 
    and `\mathcal{N}` is a differentiable manifold, such that $\Phi(U) is 
    an open subset of `\mathcal{N}`, $\Phi$ is invertible on its image and
    both `\Phi` and `\Phi^{-1}` are differentiable.

    INPUT:
    
    - ``domain`` -- domain `U` of the diffeomorphism (open subset of the start 
      manifold)
    - ``codomain`` -- codomain of the diffeomorphism (the arrival manifold or 
      some subset of it)
    - ``coord_functions`` -- (default: None) the coordinate symbolic expression 
      of the mapping in some pair of charts: list (or tuple) of the 
      coordinates of the image expressed in terms of the coordinates of 
      the considered point; if the dimension of the arrival manifold is 1, 
      a single expression is expected (not a list with a single element)
    - ``chart1`` -- (default: None) chart on domain `U` in which the 
      coordinates are given for ``coord_functions``; if none is provided, the 
      coordinates are assumed to refer to domain's default chart
    - ``chart2`` -- (default: None) chart on the codomain for the coordinates
      on the arrival manifold for ``coord_functions``; if none is provided, the 
      coordinates are assumed to refer to the codomain's default chart
    - ``name`` -- (default: None) name given to the diffeomorphism
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the 
      diffeomorphism; if none is provided, the LaTeX symbol is set to 
      ``name``
    
    .. NOTE:
    
        If ``chart1`` does not cover the entire domain `U`, the argument
        ``coord_functions` is not not sufficient to fully specify the 
        diffeomorphism; further coordinate expressions, in other charts,
        can be subsequently added by means of the method :meth:`add_expr`
    
    """
    def __init__(self, domain, codomain, coord_functions=None, 
                 chart1=None, chart2=None, name=None, latex_name=None): 
        DiffMapping.__init__(self, domain, codomain, coord_functions, chart1, 
                             chart2, name, latex_name)
        if self.domain.manifold.dim != self.codomain.manifold.dim:
            raise ValueError("The manifolds " + str(self.domain.manifold) + 
                             " and " + str(self.codomain.manifold) + 
                             " do not have the same dimension.")
        # Initialization of derived quantities:
        Diffeomorphism._init_derived(self)
    
    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "diffeomorphism"
        if self.name is not None:
            description += " '%s'" % self.name
        if self.domain == self.codomain:
            description += " on the " + str(self.domain)
        else:
            description += " between the " + str(self.domain) + \
                           " and the " + str(self.codomain)
        return description

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        DiffMapping._init_derived(self) # derived quantities of the mother class
        self._inverse = None

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        DiffMapping._del_derived(self) # derived quantities of the mother class
        self._inverse = None


    def inverse(self, chart1=None, chart2=None): 
        r"""
        Return the inverse diffeomorphism. 
        
        INPUT:
    
        - ``chart1`` -- (default: None) chart in which the computation of the 
          inverse is performed if necessary; if none is provided, the default 
          chart of the start domain will be used
        - ``chart2`` -- (default: None) chart in which the computation of the 
          inverse is performed if necessary; if none is provided, the default 
          chart of the arrival domain will be used
        
        OUTPUT:
        
        - the inverse diffeomorphism
        
        EXAMPLES:
        
        The inverse of a rotation in the Euclidean plane::
        
            sage: M = Manifold(2, 'R^2', r'\RR^2')
            sage: c_cart.<x,y> = M.chart('x y')
            sage: # A pi/3 rotation around the origin:
            sage: rot = M.diffeomorphism(M, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2), name='R')
            sage: rot.inverse() 
            diffeomorphism 'R^(-1)' on the 2-dimensional manifold 'R^2'
            sage: rot.inverse().view()
            R^(-1): R^2 --> R^2, (x, y) |--> (1/2*sqrt(3)*y + 1/2*x, -1/2*sqrt(3)*x + 1/2*y)

        Checking that applying successively the diffeomorphism and its 
        inverse results in the identity::
        
            sage: (a, b) = var('a b')
            sage: p = M.point((a,b)) # a generic point on M
            sage: q = rot(p)
            sage: p1 = rot.inverse()(q)
            sage: p1 == p 
            True

        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        from utilities import simplify_chain
        if self._inverse is not None:
            return self._inverse
            
        if chart1 is None: chart1 = self.domain.def_chart
        if chart2 is None: chart2 = self.codomain.def_chart
        coord_map = self.coord_expression[(chart1, chart2)]
        n1 = len(chart1.xx)
        n2 = len(chart2.xx)
        
        # New symbolic variables (different from chart2.xx to allow for a 
        #  correct solution even when chart2 = chart1):
        x2 = [ SR.var('xxxx' + str(i)) for i in range(n2) ]
        equations = [ x2[i] == coord_map.functions[i].express 
                      for i in range(n2) ]
        solutions = solve(equations, chart1.xx, solution_dict=True)
        if len(solutions) == 0: 
            raise ValueError("No solution found")
        if len(solutions) > 1: 
            raise ValueError("Non-unique solution found")
            
        #!# This should be the Python 2.7 form: 
        # substitutions = {x2[i]: chart2.xx[i] for i in range(n2)}
        #
        # Here we use a form compatible with Python 2.6:
        substitutions = dict([(x2[i], chart2.xx[i]) for i in range(n2)])
       
        inv_functions = [solutions[0][chart1.xx[i]].subs(substitutions) 
                           for i in range(n1)]
        for i in range(n1):
            x = inv_functions[i]
            try:
                inv_functions[i] = simplify_chain(x)
            except AttributeError:
                pass
        if self.name is None:
            name = None
        else:
            name = self.name + '^(-1)'
        
        if self.latex_name is None:
            latex_name = None
        else:
            latex_name = self.latex_name + r'^{-1}'
        self._inverse = Diffeomorphism(self.codomain, self.domain, 
                                       inv_functions, chart2, chart1,
                                       name=name, latex_name=latex_name)
        return self._inverse

        
#******************************************************************************

class IdentityMapping(Diffeomorphism):
    r"""
    Class for identity mapping on an open subset of some differentiable 
    manifold.

    INPUT:
    
    - ``domain`` -- open subset of some differentiable manifold
    - ``name`` -- (default: None) name given to the identity mapping; if None,
      it is set to 'Id_U', where 'U' is the domain's name.
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the identity 
      mapping; if None, it is set to `\mathrm{Id}_U`, where `U` is the symbol
      denoting the domain. 
      
    EXAMPLES:
    
    Identity mapping on a open subset of a 2-dimensional manifold::
    
        sage: M = Manifold(2, 'M')
        sage: U = M.open_domain('U')
        sage: c_xy.<x, y> = U.chart('x y')
        sage: i = U.identity_mapping() ; i
        identity mapping 'Id_U' on the open domain 'U' on the 2-dimensional manifold 'M'
        sage: latex(i)
        \mathrm{Id}_{U}

    The identity mapping acting on a point::
    
        sage: p = U.point((1,-2), name='p')
        sage: i(p)
        point 'p' on 2-dimensional manifold 'M'
        sage: i(p) == p
        True
        sage: i(p) is p
        True
    
    The coordinate expression of the identity mapping::
    
        sage: i.view()
        Id_U: U --> U, (x, y) |--> (x, y)

    """
    def __init__(self, domain, name=None, latex_name=None):
        if name is None:
            name = 'Id_' + domain.name
        if latex_name is None:
            latex_name = r'\mathrm{Id}_{' + domain.latex_name + r'}'
        Diffeomorphism.__init__(self, domain, domain, name=name, 
                                latex_name=latex_name)
        for chart in domain.atlas:
            coord_functions = chart[:]
            self.coord_expression[(chart, chart)] = \
                                    MultiFunctionChart(chart, *coord_functions)
        self._inverse = self 
    
    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "identity mapping '%s'" % self.name
        description += " on the " + str(self.domain)
        return description

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        DiffMapping._del_derived(self) # derived quantities of the mother class
        self._inverse = None
        
    def set_expr(self, chart1, chart2, coord_functions): 
        r"""
        Redefinition of :meth:`DiffMapping.set_expr`: should not be used
        """
        raise NotImplementedError("IdentityMapping.set_expr must not be used.")

    def multi_function_chart(self, chart1=None, chart2=None):
        r""" 
        Return the functions of the coordinates representing the differentiable
        mapping in a given pair of charts.
        
        If these functions are not already known, they are computed from known 
        ones by means of change-of-chart formulas. 
        
        INPUT:
        
        - ``chart1`` -- (default: None) chart on the mapping's domain; if None, 
          the domain's default chart is assumed
        - ``chart2`` -- (default: None) chart on the mapping's codomain; if 
          None,  the codomain's default chart is assumed

        OUTPUT:
        
        - instance of class 
          :class:`~sage.geometry.manifolds.chart.MultiFunctionChart` 
          representing the identity mapping in the above two charts

        EXAMPLES:

        """
        def_chart = self.domain.def_chart
        if chart1 is None:
            chart1 = def_chart
        if chart2 is None:
            chart2 = def_chart
        if (chart1, chart2) not in self.coord_expression:
            if chart1 == chart2:
                coord_functions = chart1[:]
                self.coord_expression[(chart1, chart1)] = \
                                   MultiFunctionChart(chart1, *coord_functions)
            else:
                return DiffMapping.multi_function_chart(self, chart1, chart2)
        return self.coord_expression[(chart1, chart2)]


    def inverse(self, chart1=None, chart2=None): 
        r"""
        Return the inverse diffeomorphism, i.e. itself !

        This is a redefinition of :meth:`Diffeomorphism.inverse`
        
        INPUT:
    
        - ``chart1`` -- (default: None) unsued
        - ``chart2`` -- (default: None) unsued
        
        OUTPUT:
        
        - the identity mapping
        
        """
        return self
            
    def __call__(self, p, chart1=None, chart2=None):
        r"""
        Image of a point.

        This is a redefinition of :meth:`DiffMapping.__call__`

        INPUT:
    
        - ``p`` -- point on the mapping's domain (type: 
          :class:`~sage.geometry.manifolds.point.Point`)
        - ``chart1`` -- (default: None) unused 
        - ``chart2`` -- (default: None) unused
        
        OUTPUT:

        - point ``p`` (since ``self`` is the identity mapping
        
        """
        # no test of p being a point in the domain (for efficiency)
        return p

    def pullback(self, tensor):
        r""" 
        Pullback operator associated with the identity mapping.
        
        This is a redefinition of :meth:`DiffMapping.pullback`
        """
        # no test for efficiency
        return tensor

        
