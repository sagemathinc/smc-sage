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

from sage.symbolic.ring import SR
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule, \
                                                    FiniteRankFreeModuleElement

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

    Check if vector `v` belongs to `T`::

        sage: v in T
        True

    Checking the object's type::  

        sage: type(v)
        <class 'sage.geometry.manifolds.tangentspace.TangentSpace_with_category.element_class'>

    Unnamed element of the tangent space:: 

        sage: w = T._an_element_()
        sage: w
        tangent vector at point 'P' on 2-dimensional manifold 'R^2' 
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

#******************************************************************************

class TangentSpace(FiniteRankFreeModule):
    r"""
    Tangent space at a given point on a differentiable manifold.

    INPUT:
    
    - ``point`` -- point at which the tangent space is defined. 

    EXAMPLES:

    Tangent space on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: p = M.point((-1,2), name='p')
        sage: T = p.tangent_space() ; T
        tangent space at point 'p' on 2-dimensional manifold 'M'

    Tangent spaces belong to a subclass of :class:`TangentSpace`::
    
        sage: type(T)
        <class 'sage.geometry.manifolds.tangentspace.TangentSpace_with_category'>

    They are free modules of finite rank over Sage Symbolic Rank (actually
    vector space of finite dimension over `\RR`)::
    
        sage: isinstance(T, FiniteRankFreeModule)
        True
        sage: T.base_ring()
        Symbolic Ring
        sage: T.rank()
        2
        sage: T.dim()
        2

    """

    Element = TangentVector

    def __init__(self, point):
        manif = point._manifold
        name = "T_" + str(point._name) + " " + str(manif._name)
        latex_name = r"T_{" + str(point._name) + "}\," + str(manif._latex_name)
        self._point = point
        self._manif = manif
        FiniteRankFreeModule.__init__(self, SR, manif._dim, name=name, 
                                      latex_name=latex_name, 
                                      start_index=manif._sindex)
        # Initialization of bases of the tangent space from existing vector
        # frames around the point:
        for frame in manif._frames:
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
                self._point._frame_bases[frame] = basis

    def _repr_(self):
        r"""
        String representation of the object.
        """

        description = "tangent space at " + str(self._point)  
        return description
       
    def dim(self):
        r"""
        Return the vector space dimension of ``self``. 
        """
        return self._rank
    
