r"""
Tangent space module

Defines the tangent space over a point `p` at some open subset `U`, 
possibly belonging to a manifold `M`. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version
"""

#******************************************************************************
#       Copyright (C) 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule, FiniteRankFreeModuleElement

class TangentVector(FiniteRankFreeModuleElement):
    r""" 
    Tangent vector `v` beloning to the tangent space `T` at `p`

    EXAMPLES:

    Tangent space on `\RR^2` and some tangent vector::

        sage: M=Manifold(2, "R^2")
        sage: c.<x,y> = M.chart()
        sage: p = M.point((1, 2), name='P'); p
        point 'P' on 2-dimensional manifold 'R^2'
        sage: T = p.tangent_space(); T
        tangent space at a point 'P' on 2-dimensional manifold 'R^2'
        sage: v = T.tangent_vector('V'); v
        tangent vector V at point 'P' on 2-dimensional manifold 'R^2'
        sage: w = T._an_element_()
        sage: w
        tangent vector at point 'P' on 2-dimensional manifold 'R^2' 
    """
    
    def __init__(self, parent, name=None, latex_name=None):
        from sage.symbolic.ring import SR

        self._point = parent._point
        self._dim = parent._dim    
        self.name = name 
        self.latex_name = latex_name
        vmodule = FiniteRankFreeModule(SR, self._dim) 
        FiniteRankFreeModuleElement.__init__(self, vmodule, name=name, latex_name=latex_name)

    def _repr_(self): 
        r"""
        String representation of the object.
        """

        desc = "tangent vector"
        if self.name:
            desc += " " + str(self.name) 
        desc += " at " + str(self._point)
        return desc

class TangentSpace(FiniteRankFreeModule):
    r"""
    Tangent space at a point `p` beloning to 
    an open domain `U` of some manifold `M`. 

    INPUT:
    
    - ``point`` -- point `p` beloning to `U` at which the tangent space is defined. 

    EXAMPLES:

    Tangent space on `\RR^2`::

        sage: M=Manifold(2, "R^2")
        sage: c.<x,y> = M.chart()
        sage: p = M.point((1, 2), name='P'); p
        point 'P' on 2-dimensional manifold 'R^2'
        sage: T = p.tangent_space(); T
        tangent space at a point 'P' on 2-dimensional manifold 'R^2'
    """

    Element = TangentVector

    def __init__(self, point):
        self._point = point 
        self._dim = (point._manifold)._dim

    def _repr_(self):
        r"""
        String representation of the object.
        """

        description = "tangent space at a " + str(self._point)  
        return description
       
    def _element_constructor_(self, point, name=None, latex_name=None):
        r"""
        Construct an element of the module
        """
        return self.element_class(self, point, name=name, latex_name=latex_name)

    def _an_element_(self):
        r"""
        Construct some (unamed) element
        """
        return self.element_class(self)

    def tangent_vector(self, name):
        r"""
        Construct an element with name
        """
        return self.element_class(self, name=name)

