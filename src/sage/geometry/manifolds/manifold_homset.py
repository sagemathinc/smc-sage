r"""
Sets of morphisms between manifolds

The class :class:`ManifoldHomset` implements sets of morphisms between
two differentiable manifolds over `\RR`, a morphism being a *differentiable
mapping* in the present context. 

AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

- S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*, vol. 1,
  Interscience Publishers (New York) (1963)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.homset import Homset
from sage.categories.sets_cat import Sets
from sage.geometry.manifolds.diffmapping import DiffMapping, Diffeomorphism

class ManifoldHomset(Homset):
    r"""
    Set of differentiable mappings between two differentiable manifolds. 

    Given two differentiable manifolds `M` and `N` over `\RR`, the class
    :class:`ManifoldHomset` implements the set `\mathrm{Hom}(U,V)` of morphisms
    (i.e. differentiable mappings) `U\rightarrow V`, where `U` is an open
    subset of `M` and `V` an open subset of `N`. Note that, as open subsets of
    manifolds, `U` and `V` are manifolds by themselves. 

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`.

    INPUT:

    - ``domain`` -- open subset `U\subset M` (domain of the morphisms),
      as an instance of
      :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`
    - ``codomain`` -- open subset `V\subset N` (codomain of the morphisms),
      as an instance of
      :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`
    - ``name`` -- (default: ``None``) string; name given to the homset; if
      none is provided, Hom(U,V) will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      homset; if none is provided, `\mathrm{Hom}(U,V)` will be used

    EXAMPLES:

    """

    Element = DiffMapping

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        TESTS

        """
        from sage.geometry.manifolds.domain import ManifoldOpenSubset
        if not isinstance(domain, ManifoldOpenSubset):
            raise TypeError("domain = {} is not an ".format(domain) +
                            "instance of ManifoldOpenSubset")
        if not isinstance(codomain, ManifoldOpenSubset):
            raise TypeError("codomain = {} is not an ".format(codomain) +
                            "instance of ManifoldOpenSubset")
        Homset.__init__(self, domain, codomain)
        if name is None:
            self._name = "Hom(" + domain._name + "," + codomain._name + ")"
        else:
            self._name = name
        if latex_name is None:
            self._latex_name = \
                    r"\mathrm{Hom}\left(" + domain._latex_name + "," + \
                    codomain._latex_name + r"\right)"
        else:
            self._latex_name = latex_name

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def __call__(self, *args, **kwds):
        r"""
        To bypass Homset.__call__, enforcing Parent.__call__ instead.

        EXAMPLES:
        """
        return Parent.__call__(self, *args, **kwds)


    #### Methods required for any Parent

    def _element_constructor_(self, coord_functions, chart1=None,
                              chart2=None, name=None, latex_name=None):
        r"""
        Construct an element of ``self``, i.e. a differentiable mapping
        U --> V, where U is the domain of ``self`` and V its codomain.

        INPUT:

        - ``coord_functions`` -- must be either
    
          - (i) a dictionary of
            the coordinate expressions (as lists (or tuples) of the
            coordinates of the image expressed in terms of the coordinates of
            the considered point) with the pairs of charts (chart1, chart2)
            as keys (chart1 being a chart on `U` and chart2 a chart on `N`)
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``
    
          In both cases, if the dimension of the arrival manifold is 1,
          a single coordinate expression is expected (not a list or tuple with a
          single element)
        - ``chart1`` -- (default: None; used only in case (ii) above) chart on
          domain `U` defining the start coordinates involved in ``coord_functions``
          for case (ii); if none is provided, the coordinates are assumed to
          refer to `U`'s default chart
        - ``chart2`` -- (default: None; used only in case (ii) above) chart on the
          codomain defining the arrival coordinates involved in ``coord_functions``
          for case (ii); if none is provided, the coordinates are assumed to
          refer to the codomain's default chart
        - ``name`` -- (default: None) name given to the differentiable mapping
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the
          differentiable mapping; if none is provided, the LaTeX symbol is set to
          ``name``
    
        .. NOTE::
    
            If the information passed by means of the argument
            ``coord_functions`` is not sufficient to fully specify the
            differential mapping (for instance case (ii) with ``chart1`` not
            covering the entire domain `U`), further coordinate expressions,
            in other charts, can be subsequently added by means of the method
            :meth:`~sage.geometry.manifolds.diffmapping.DiffMapping.add_expr`

        """
        # Standard construction
        return self.element_class(self, coord_functions=coord_functions,
                                  chart1=chart1, chart2=chart2, name=name,
                                  latex_name=latex_name)

    def _an_element_(self):
        r"""
        Construct some element.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M,N) 
            sage: f = H._an_element_() ; f
            differentiable mapping from 2-dimensional manifold 'M' to
             3-dimensional manifold 'N'
            sage: f.display()
            M --> N
               (x, y) |--> (u, v, w) = (0, 0, 0)
            sage: p = M((-2,3)) ; p
            point on 2-dimensional manifold 'M'
            sage: f(p)
            point on 3-dimensional manifold 'N'
            sage: f(p).coord(Y)
            (0, 0, 0)
            sage: TestSuite(f).run()

        """
        dom = self.domain()
        codom = self.codomain()
        # A constant diff. mapping is constructed:
        chart2 = codom.default_chart()
        target_point = chart2.domain().an_element()
        target_coord = target_point.coord(chart2)
        coord_functions = {}
        for chart in dom.atlas():
            coord_functions[(chart, chart2)] = target_coord
        return self.element_class(self, coord_functions=coord_functions)

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent.

        EXAMPLES:


        """
        #!# for the time being:
        return False


    #### End of methods required for any Parent


#******************************************************************************

class DiffeoSet(UniqueRepresentation, Parent):
    r"""
    Set of diffeomorphisms between two differentiable manifolds. 

    Given two differentiable manifolds `M` and `N` over `\RR`, the class
    :class:`DiffeoSet` implements the set `\mathrm{Diff}(U,V)` of
    diffeomorphisms `U\rightarrow V`, where `U` is an open
    subset of `M` and `V` an open subset of `N`,  a *diffeomorphism* being a
    differentiable mapping `U\rightarrow V` that is bijective and whose inverse
    is differentiable. 

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.geometry.manifolds.diffmapping.Diffeomorphism`.

    INPUT:

    - ``domain`` -- open subset `U\subset M` (domain of the diffeomorphisms),
      as an instance of
      :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`
    - ``codomain`` -- open subset `V\subset N` (codomain of the
      diffeommorphisms),
      as an instance of
      :class:`~sage.geometry.manifolds.domain.ManifoldOpenSubset`
    - ``name`` -- (default: ``None``) string; name given to the set of
      diffeomorphism; if none is provided, Diff(U,V) will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      set of diffeomorphism; if none is provided, `\mathrm{Diff}(U,V)` will be
      used

    EXAMPLES:

    """

    Element = Diffeomorphism

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        TESTS

        """
        from sage.geometry.manifolds.domain import ManifoldOpenSubset
        if not isinstance(domain, ManifoldOpenSubset):
            raise TypeError("domain = {} is not an ".format(domain) +
                            "instance of ManifoldOpenSubset")
        if not isinstance(codomain, ManifoldOpenSubset):
            raise TypeError("codomain = {} is not an ".format(codomain) +
                            "instance of ManifoldOpenSubset")
        Parent.__init__(self, category=Sets())
        self._domain = domain
        self._codomain = codomain
        if name is None:
            self._name = "Diff(" + domain._name + "," + codomain._name + ")"
        else:
            self._name = name
        if latex_name is None:
            self._latex_name = \
                    r"\mathrm{Diff}\left(" + domain._latex_name + "," + \
                    codomain._latex_name + r"\right)"
        else:
            self._latex_name = latex_name

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name


    #### Methods required for any Parent

    def _element_constructor_(self, coord_functions, chart1=None,
                              chart2=None, name=None, latex_name=None,
                              is_identity=False):
        r"""
        Construct an element of ``self``, i.e. a diffeomorphism
        U --> V, where U is the domain of ``self`` and V its codomain.

        INPUT:
    
          - ``coord_functions`` -- (default: None) if not None, must be either
    
          - (i) a dictionary of
            the coordinate expressions (as lists (or tuples) of the
            coordinates of the image expressed in terms of the coordinates of
            the considered point) with the pairs of charts (chart1, chart2)
            as keys (chart1 being a chart on `U` and chart2 a chart on `N`)
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``
    
          In both cases, if the dimension of the manifolds is 1,
          a single coordinate expression is expected (not a list or tuple with
          a single element)
        - ``chart1`` -- (default: None; used only in case (ii) above) chart on
          domain `U` defining the start coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to `U`'s default chart
        - ``chart2`` -- (default: None; used only in case (ii) above) chart on
          the codomain defining the arrival coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the codomain's default chart
        - ``name`` -- (default: None) name given to the diffeomorphism
        - ``latex_name`` -- (default: None) LaTeX symbol to denote the
          diffeomorphism; if none is provided, the LaTeX symbol is set to
          ``name``
        - ``is_identity`` -- (default: ``False``) determines whether the
          constructed diffeomorphism is the identity map; if set to ``True``,
          then V must be U and the entries ``coord_functions``, ``chart1``
          and ``chart2`` are not used.
    
        .. NOTE::
    
            If the information passed by means of the argument
            ``coord_functions`` is not sufficient to fully specify the
            diffeomorphism (for instance case (ii) with ``chart1`` not
            covering the entire domain `U`), further coordinate expressions,
            in other charts, can be subsequently added by means of the method
            :meth:`~sage.geometry.manifolds.diffmapping.DiffMapping.add_expr`

        """
        # Standard construction
        return self.element_class(self, coord_functions=coord_functions,
                                  chart1=chart1, chart2=chart2, name=name,
                                  latex_name=latex_name,
                                  is_identity=is_identity)

    def _an_element_(self):
        r"""
        Construct some element.

        EXAMPLE:

        """
        if self.domain() == self.codomain():
            # return the identity map
            return self.element_class(self, is_identity=True)
        else:
            raise NotImplementedError("generation of a generic element is " +
                "not implemented when the domain and codomain do not coincide")

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent.

        EXAMPLES:


        """
        #!# for the time being:
        return False


    #### End of methods required for any Parent

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain








