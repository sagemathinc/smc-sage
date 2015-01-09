r"""
Differential form modules

The set `\Lambda^p(U,\Phi)` of `p`-forms along an open subset `U` of
some manifold `S` with values in a open subset `V=\Phi(U)` of a manifold `M`
(`\Phi` being a differentiable mapping `U\rightarrow M`; possibly
`\Phi=\mathrm{Id}`, `S=M` and `U=V`) is a module over the algebra `C^\infty(U)`
of differentiable scalar fields on `U`. It is a free module iff `V` is
parallelizable. Accordingly, two classes are implementing `\Lambda^p(U,\Phi)`:

- :class:`DiffFormModule` for differential forms with values in a generic (in
  practice, not parallelizable) open set `V`
- :class:`DiffFormFreeModule` for differential forms with values in a
  parallelizable open set `V`

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

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.modules import Modules
from sage.tensor.modules.ext_pow_free_module import ExtPowerFreeModule
from diffform import DiffForm, DiffFormParal
from tensorfield import TensorField, TensorFieldParal

class DiffFormModule(UniqueRepresentation, Parent):
    r"""
    Module of differential forms of a given degree `p` (`p`-forms) along an
    open subset `U` of some manifold `S` with values in an open subset `V` of
    a manifold `M`.

    Given an open subset `U` of a manifold `S` and a differentiable mapping
    `\Phi: U \rightarrow V = \Phi(U) \subset M`, where `M` is a differentiable
    manifold, the set `\Lambda^p(U,\Phi)` of `p`-forms along `U` with values
    in `V` is a module over `C^\infty(U)`, the commutative algebra of
    differentiable scalar fields on `U`.
    The standard case of `p`-forms *on* a manifold corresponds to
    `\Phi=\mathrm{Id}`, `U=V` and `S=M`. Another common case is `\Phi` being an
    immersion.

    This class implements `\Lambda^p(U,\Phi)` in the case where `V=\Phi(U)` is
    not assumed to be parallelizable; the module `\Lambda^p(U,\Phi)` is then
    not necessarily free. If `V` is parallelizable, the class
    :class:`DiffFormFreeModule` must be used instead. 
    
    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.geometry.manifolds.diffform.DiffForm`.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the mapping `\Phi:\; U \rightarrow V`.
    - ``degree`` -- positive integer; the degree `p` of the differential forms

    """

    Element = DiffForm

    def __init__(self, vector_field_module, degree):
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        name = "/\^{}(".format(degree) + domain._name
        latex_name = r"\Lambda^{" + str(degree) + r"}\left(" + \
                     domain._latex_name
        if dest_map is domain._identity_map:
            name += ")"
            latex_name += r"\right)"
        else:
            name += "," + dest_map._name + ")"
            latex_name += "," + dest_map._latex_name + r"\right)"
        self._vmodule = vector_field_module
        self._degree = degree
        self._name = name
        self._latex_name = latex_name
        # the member self._ring is created for efficiency (to avoid calls to
        # self.base_ring()):
        self._ring = domain.scalar_field_algebra()
        Parent.__init__(self, base=self._ring, category=Modules(self._ring))
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain
        # NB: self._zero_element is not constructed here, since no element
        # can be constructed here, to avoid some infinite recursion.

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None):
        r"""
        Construct a differential form.

        """
        if comp == 0:
            if not hasattr(self, '_zero_element'):
                self._zero_element = self._element_constructor_(name='zero',
                                                                latex_name='0')
                for frame in self._domain._frames:
                    if self._dest_map.restrict(frame._domain) == \
                                                               frame._dest_map:
                        self._zero_element.add_comp(frame)
                        # (since new components are initialized to zero)
            return self._zero_element
        if isinstance(comp, DiffForm):
            # coercion by domain restriction
            if self._degree == comp._tensor_type[1] and \
               self._domain.is_subset(comp._domain) and \
               self._ambient_domain.is_subset(comp._ambient_domain):
                return comp.restrict(self._domain)
            else:
                raise TypeError("cannot coerce the {}".format(comp) +
                                " to an element of {}".format(self))
        if isinstance(comp, TensorField):
            # coercion of a tensor of type (0,1) to a linear form
            tensor = comp # for readability
            if tensor.tensor_type() == (0,1) and self._degree == 1 and \
                                         tensor._vmodule is self._vmodule:
                resu = self.element_class(self._vmodule, 1, name=tensor._name,
                                          latex_name=tensor._latex_name)
                for dom, rst in tensor._restrictions.iteritems():
                    resu._restrictions[dom] = dom.diff_form_module(1)(rst)
                return resu
            else:
                raise TypeError("cannot coerce the {}".format(tensor) +
                                " to an element of {}".format(self))
        # standard construction
        resu = self.element_class(self._vmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _an_element_(self):
        r"""
        Construct some (unamed) differential form.

        """
        resu = self.element_class(self._vmodule, self._degree)
        #!# a zero element is returned
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.
        
        """
        from tensorfield_module import TensorFieldModule
        if isinstance(other, (DiffFormModule, DiffFormFreeModule)):
            # coercion by domain restriction
            return self._degree == other._degree and \
                   self._domain.is_subset(other._domain) and \
                   self._ambient_domain.is_subset(other._ambient_domain)
        if isinstance(other, TensorFieldModule):
            # coercion of a type-(0,1) tensor to a linear form
            return self._vmodule is other._vmodule and self._degree == 1 and \
               other.tensor_type() == (0,1)
        return False

    #### End of Parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        """
        description = "Module "
        if self._name is not None:
            description += self._name + " "
        description += "of {}-forms ".format(self._degree)
        if self._dest_map is self._domain._identity_map:
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} mapped into the {}".format(
                                            self._domain, self._ambient_domain)
        return description

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.
        
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def base_module(self):
        r"""
        Return the vector field module on which ``self`` is constructed.

        OUTPUT:

        - instance of class
          :class:`~sage.geometry.manifolds.vectorfield_module.VectorFieldModule`
          representing the module on which the tensor module is defined.

        """
        return self._vmodule


#******************************************************************************

class DiffFormFreeModule(ExtPowerFreeModule):
    r"""
    Module of differential forms of a given degree `p` (`p`-forms) along an
    open subset `U` of some manifold `S` with values in a parallelizable open
    subset `V` of a manifold `M`.

    Given an open subset `U` of a manifold `S` and a differentiable mapping
    `\Phi: U \rightarrow V = \Phi(U) \subset M`, where `M` is a differentiable
    manifold, the set `\Lambda^p(U,\Phi)` of `p`-forms along `U` with values
    in `V` is a module over `C^\infty(U)`, the commutative algebra of
    differentiable scalar fields on `U`.
    The standard case of `p`-forms *on* a manifold corresponds to
    `\Phi=\mathrm{Id}`, `U=V` and `S=M`. Another common case is `\Phi` being an
    immersion.

    This class implements `\Lambda^p(U,\Phi)` in the case where `V=\Phi(U)` is
    parallelizable; `\Lambda^p(U,\Phi)` is then a *free* module. If `V` is not
    parallelizable, the class :class:`DiffFormModule` must be used instead. 
    
    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.geometry.manifolds.diffform.DiffFormParal`.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` associated with the mapping `\Phi:\; U \rightarrow V`.
    - ``degree`` -- positive integer; the degree `p` of the differential forms

    """

    Element = DiffFormParal

    def __init__(self, vector_field_module, degree):
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        name = "/\^{}(".format(degree) + domain._name
        latex_name = r"\Lambda^{" + str(degree) + r"}\left(" + \
                     domain._latex_name
        if dest_map is domain._identity_map:
            name += ")"
            latex_name += r"\right)"
        else:
            name += "," + dest_map._name + ")"
            latex_name += "," + dest_map._latex_name + r"\right)"
        ExtPowerFreeModule.__init__(self, vector_field_module, degree,
                                    name=name, latex_name=latex_name)
        self._domain = domain
        self._dest_map = dest_map
        self._ambient_domain = vector_field_module._ambient_domain

    #### Parent methods

    def _element_constructor_(self, comp=[], frame=None, name=None,
                              latex_name=None, sym=None, antisym=None):
        r"""
        Construct a differential form.

        """
        if comp == 0:
            return self._zero_element
        if isinstance(comp, DiffForm):
            # coercion by domain restriction
            if self._degree == comp._tensor_type[1] and \
               self._domain.is_subset(comp._domain) and \
               self._ambient_domain.is_subset(comp._ambient_domain):
                return comp.restrict(self._domain)
            else:
                raise TypeError("Cannot coerce the {} ".format(comp) +
                                "to a differential form in {}".format(self))
        if isinstance(comp, TensorFieldParal):
            # coercion of a tensor of type (0,1) to a linear form
            tensor = comp # for readability
            if tensor.tensor_type() == (0,1) and self._degree == 1 and \
                                         tensor._fmodule is self._fmodule:
                resu = self.element_class(self._fmodule, 1, name=tensor._name,
                                          latex_name=tensor._latex_name)
                for frame, comp in tensor._components.iteritems():
                    resu._components[frame] = comp.copy()
                return resu
            else:
                raise TypeError("cannot coerce the {}".format(tensor) +
                                " to an element of {}".format(self))
        resu = self.element_class(self._fmodule, self._degree, name=name,
                                  latex_name=latex_name)
        if comp != []:
            resu.set_comp(frame)[:] = comp
        return resu

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.
        
        """
        from tensorfield_module import TensorFieldFreeModule
        if isinstance(other, (DiffFormModule, DiffFormFreeModule)):
            # coercion by domain restriction
            return self._degree == other._degree and \
                   self._domain.is_subset(other._domain) and \
                   self._ambient_domain.is_subset(other._ambient_domain)
        if isinstance(other, TensorFieldFreeModule):
            # coercion of a type-(0,1) tensor to a linear form
            return self._fmodule is other._fmodule and self._degree == 1 and \
               other.tensor_type() == (0,1)
        return False

    #### End of Parent methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        """
        description = "Free module "
        if self._name is not None:
            description += self._name + " "
        description += "of {}-forms ".format(self._degree)
        if self._dest_map is self._domain._identity_map:
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} mapped into the {}".format(
                                            self._domain, self._ambient_domain)
        return description

