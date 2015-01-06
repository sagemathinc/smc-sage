r"""
Exterior powers of dual free modules

Given a free module `M` of finite rank over a commutative ring `R`
and a positive integer `p`, the *p-th exterior power* of the dual of `M` is the
set `\Lambda^p(M^*)` of all alternating forms of degree `p` on `M`, i.e. of
all multilinear maps

.. MATH::

    \underbrace{M\times\cdots\times M}_{p\ \; \mbox{times}}
    \longrightarrow R

that vanish whenever any of two of their arguments are equals.
Note that `\Lambda^1(M^*) = M^*` (the dual of `M`). 

`\Lambda^p(M^*)` is a free module of rank `\left({n\atop p}\right)` over `R`,
where `n` is the rank of `M`. 
Accordingly, exterior powers of free modules are implemented by a class,
:class:`ExtPowerFreeModule`, which inherits from the class
:class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.tensor.modules.free_module_alt_form import FreeModuleAltForm

class ExtPowerFreeModule(FiniteRankFreeModule):
    r"""
    Class for the exterior powers of the dual of a free module of finite rank
    over a commutative ring.

    Given a free module `M` of finite rank over a commutative ring `R`
    and a positive integer `p`, the *p-th exterior power* of the dual of `M` is
    the set `\Lambda^p(M^*)` of all alternating forms of degree `p` on `M`,
    i.e. of all multilinear maps

    .. MATH::

        \underbrace{M\times\cdots\times M}_{p\ \; \mbox{times}}
        \longrightarrow R

    that vanish whenever any of two of their arguments are equals.
    Note that `\Lambda^1(M^*) = M^*` (the dual of `M`). 

    `\Lambda^p(M^*)` is a free module of rank `\left({n\atop p}\right)` over
    `R`, where `n` is the rank of `M`. 
    Accordingly, the class :class:`ExtPowerFreeModule` inherits from the class
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`.

    The class :class:`ExtPowerFreeModule` is a Sage *parent* class, whose
    *element* class is
    :class:`~sage.tensor.modules.free_module_alt_form.FreeModuleAltForm`.

    INPUT:

    - ``fmodule`` -- free module `M` of finite rank (must be an instance of
      :class:`FiniteRankFreeModule`)
    - ``degree`` -- degree `p` of the alternating forms
    - ``name`` -- (default: ``None``) string; name given to `\Lambda^p(M^*)`
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      `\Lambda^p(M^*)`; if none is provided, it is set to ``name``

    """

    Element = FreeModuleAltForm

    def __init__(self, fmodule, degree, name=None, latex_name=None):
        r"""
        TEST::

            sage: from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule
            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')

        """
        from sage.functions.other import binomial
        self._fmodule = fmodule
        self._degree = degree
        rank = binomial(fmodule._rank, degree)
        self._zero_element = 0 # provisory (to avoid infinite recursion in what
                               # follows)
        if degree == 1:  # case of the dual
            if name is None and fmodule._name is not None:
                name = fmodule._name + '*'
            if latex_name is None and fmodule._latex_name is not None:
                latex_name = fmodule._latex_name + r'^*'
        FiniteRankFreeModule.__init__(self, fmodule._ring, rank, name=name,
                                      latex_name=latex_name,
                                      start_index=fmodule._sindex,
                                    output_formatter=fmodule._output_formatter)
        # Unique representation:
        if self._degree in self._fmodule._dual_exterior_powers:
            raise ValueError("the " + str(degree) + "th exterior power of " +
                             "the dual of " + str(self._fmodule) +
                             " has already been created")
        else:
            self._fmodule._dual_exterior_powers[self._degree] = self
        # Zero element
        self._zero_element = self._element_constructor_(name='zero',
                                                        latex_name='0')
        for basis in self._fmodule._known_bases:
            self._zero_element._components[basis] = \
                                            self._zero_element._new_comp(basis)
            # (since new components are initialized to zero)

    #### Methods required for any Parent
    
    def _element_constructor_(self, comp=[], basis=None, name=None,
                              latex_name=None):
        r"""
        Construct an alternating form.

        EXAMPLES::


        """
        if comp == 0:
            return self._zero_element
        resu = self.element_class(self._fmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp:
            resu.set_comp(basis)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) alternating form.

        EXAMPLES::

        """
        resu = self.element_class(self._fmodule, self._degree)
        if self._fmodule._def_basis is not None:
            sindex = self._fmodule._sindex
            ind = [sindex + i for i in range(resu._tensor_rank)]
            resu.set_comp()[ind] = self._fmodule._ring.an_element()
        return resu


    #### End of methods required for any Parent

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        """
        if self._degree == 1:
            return "Dual of the {}".format(self._fmodule)
        description = "{}".format(self._degree)
        if self._degree == 2:
            description += "nd"
        elif self._degree == 3:
            description += "rd"
        else:
            description += "th"
        description += " exterior power of the dual of the {}".format(
                                                                 self._fmodule)
        return description

    def base_module(self):
        r"""
        Return the free module on which ``self`` is constructed.

        OUTPUT:

        - instance of :class:`FiniteRankFreeModule` representing the free
          module on which the exterior power is defined.

        EXAMPLE:

        """
        return self._fmodule

    def degree(self):
        r"""
        Return the degree of ``self``.

        OUTPUT:

        - integer `p` such that ``self`` is the exterior power `\Lambda^p(M^*)`

        EXAMPLE:

        """
        return self._degree
