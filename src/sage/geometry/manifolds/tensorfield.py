r"""
Tensor fields

The class :class:`TensorField` implements tensor fields on differentiable 
manifolds over `\RR`. 

A tensor field of type `(k,\ell)` is a field of multilinear maps:

.. MATH::

    \underbrace{T_p^*M\times\cdots\times T_p^*M}_{k\ \; \mbox{times}}
    \times \underbrace{T_p M\times\cdots\times T_p M}_{\ell\ \; \mbox{times}}
    \longrightarrow \RR
    
where `T_p M` stands for the tangent space at the point `p` on the
manifold `M` and `T_p^*M` for its dual vector space. The integer `k+\ell`
is called the tensor rank. 

Various derived classes of :class:`TensorField` are devoted to specific tensor
fields:

* :class:`~sage.geometry.manifolds.vectorfield.VectorField` for vector fields 
  (rank-1 contravariant tensor fields)
* :class:`~sage.geometry.manifolds.diffform.OneForm` for 1-forms (rank-1 
  covariant tensor fields)
* :class:`~sage.geometry.manifolds.rank2field.EndomorphismField` for fields of 
  endomorphisms (type (1,1) tensor fields)
* :class:`~sage.geometry.manifolds.rank2field.SymBilinFormField` for fields of 
  symmetric bilinear forms (rank-2 symmetric covariant tensor fields)
* :class:`~sage.geometry.manifolds.diffform.DiffForm` for differential forms 
  (fully antisymmetric covariant tensor fields)


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014) : initial version

EXAMPLES:

    A tensor field of type (1,1) on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M', start_index=1)
        sage: c_xy.<x,y> = M.chart('x y')
        sage: t = M.tensor_field(1, 1, 'T') ; t
        tensor field 'T' of type (1,1) on the 2-dimensional manifold 'M'
        sage: t.tensor_rank
        2

    A just-created tensor field has no components::

        sage: t.components
        {}

    Components w.r.t. the manifold's default frame are created by providing the
    relevant indices inside square brackets::

        sage: t[1,1] = x^2

    Unset components are initialized to zero::

        sage: t[:]  # list of components w.r.t. the manifold's default vector frame
        [x^2   0]
        [  0   0]

    The full set of components w.r.t. a given vector frame is returned by the 
    method :meth:`comp`; it is an instance of the class 
    :class:`~sage.tensor.modules.comp.Components`::
    
        sage: t.comp(c_xy.frame)
        2-indices components w.r.t. coordinate frame (M, (d/dx,d/dy)) 
        sage: print type(t.comp(c_xy.frame))
        <class 'sage.tensor.modules.comp.Components'>

    The vector frame can be skipped, it is then assumed to be the
    manifold's default frame::
    
        sage: M.default_frame()
        coordinate frame (M, (d/dx,d/dy))
        sage: t.comp() is t.comp(c_xy.frame)
        True

    Individual components w.r.t. the manifold's default frame are accessed by 
    listing their indices inside double square brackets; they are scalar
    fields on the manifold, and therefore instances of the class 
    :class:`~sage.geometry.manifolds.scalarfield.ScalarField`::

        sage: t[[1,1]]
        scalar field on the 2-dimensional manifold 'M'
        sage: t[[1,1]].expr()
        x^2
        sage: t[[1,2]]
        zero scalar field on the 2-dimensional manifold 'M'
        sage: t[[1,2]].expr()
        0
        
    A direct access to the coordinate expression of some component is obtained
    via the single square brackets::
    
        sage: t[1,1] 
        x^2
        sage: t[1,1] is t[[1,1]].function_chart() # the coordinate function
        True
        sage: t[1,1] is t[[1,1]].function_chart(c_xy)
        True
        sage: t[1,1] == t[[1,1]].expr() # check the value of the coordinate function
        True
        sage: t[1,1].expr() is t[[1,1]].expr() # the symbolic expression
        True

    In other words, the single square brackets return an instance of 
    :class:`~sage.geometry.manifolds.chart.FunctionChart` that is the 
    coordinate function representing the component in some chart (by default, 
    the manifold's default chart)::
    
        sage: print type(t[1,1])    # single bracket --> FunctionChart
        <class 'sage.geometry.manifolds.chart.FunctionChart'>
        sage: print type(t[[1,1]])  # double bracket --> ScalarField
        <class 'sage.geometry.manifolds.scalarfield.ScalarFieldRing_with_category.element_class'>
    
    Expressions in a chart different from the manifold's default one are 
    obtained by specifying the chart as the last argument inside the
    single square brackets::
    
        sage: c_uv.<u,v> = M.chart('u v')
        sage: xy_to_uv = c_xy.coord_change(c_uv, x+y, x-y)  
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: t[1,1, c_uv] 
        1/4*u^2 + 1/2*u*v + 1/4*v^2

    Note that ``t[1,1, c_uv]`` is the component of the tensor t w.r.t. to 
    the coordinate frame associated to the chart (x,y) expressed in terms of 
    the coordinates (u,v). Indeed, ``t[1,1, c_uv]`` is a shortcut for 
    ``t.comp(c_xy.frame)[[1,1]].function_chart(c_uv)``::
        
        sage: t[1,1, c_uv] is t.comp(c_xy.frame)[[1,1]].function_chart(c_uv)
        True
    
    Similarly, ``t[1,1]`` is a shortcut for 
    ``t.comp(c_xy.frame)[[1,1]].function_chart(c_xy)``::
    
        sage: t[1,1] is t.comp(c_xy.frame)[[1,1]].function_chart(c_xy)            
        True
        sage: t[1,1] is t.comp()[[1,1]].function_chart()  # since c_xy.frame and c_xy are the manifold's default values                
        True

    Internally, the components are stored as a dictionary (attribute 
    :attr:`_comp` of the class 
    :class:`~sage.tensor.modules.comp.Components`) whose
    keys are the indices. Only the non-zero components and non-redundant
    components (in case of symmetries) are stored::

        sage: t.comp()._comp
        {(1, 1): scalar field on the 2-dimensional manifold 'M'}

    All the components can be set at once via [:]::
    
        sage: t[:] = [[1, -x], [x*y, 2]]
        sage: t[:]
        [  1  -x]
        [x*y   2]
        
    The different sets of components, corresponding to representations of the
    tensor in different vector frames, are stored in the dictionary 
    :attr:`components`, each item being an instance of the class 
    :class:`~sage.tensor.modules.comp.Components`::
    
        sage: t.components
        {coordinate frame (M, (d/dx,d/dy)): 2-indices components w.r.t. coordinate frame (M, (d/dx,d/dy))}
        sage: print type(t.components[c_xy.frame])
        <class 'sage.tensor.modules.comp.Components'>
        sage: print type(t.comp())
        <class 'sage.tensor.modules.comp.Components'>
        sage: t.comp() is t.components[c_xy.frame]
        True

    To set the components in a vector frame different from the manifold's 
    default one, the method :meth:`set_comp` must be employed::

        sage: e = M.vector_frame('e')
        sage: t.set_comp(e)[1,1], t.set_comp(e)[1,2] = (x+y, 0)
        sage: t.set_comp(e)[2,1], t.set_comp(e)[2,2] = (y, -3*x)
        sage: t.comp(e)
        2-indices components w.r.t. vector frame (M, (e_1,e_2))
        sage: t.comp(e)[:]
        [x + y     0]
        [    y  -3*x]

    All the components in some frame can be set at once, via the operator
    [:]::

        sage: t.set_comp(e)[:] = [[x+y, 0], [y, -3*x]]
        sage: t.comp(e)[:]  # same as above:
        [x + y     0]
        [    y  -3*x]
    
    To avoid any insconstency between the various components, the method 
    :meth:`set_comp` clears the components in other frames. 
    Accordingly, the components in the frame c_xy.frame have been deleted::
    
        sage: t.components
        {vector frame (M, (e_1,e_2)): 2-indices components w.r.t. vector frame (M, (e_1,e_2))}

    To keep the other components, one must use the method :meth:`add_comp`::
    
        sage: t = M.tensor_field(1, 1, 'T')  # Let us restart 
        sage: t[:] = [[1, -x], [x*y, 2]]  # by first setting the components in the frame c_xy.frame
        sage: # We now set the components in the frame e with add_comp:
        sage: t.add_comp(e)[:] = [[x+y, 0], [y, -3*x]]
        sage: t.components  # Both set of components are present:
        {coordinate frame (M, (d/dx,d/dy)): 2-indices components w.r.t. coordinate frame (M, (d/dx,d/dy)), vector frame (M, (e_1,e_2)): 2-indices components w.r.t. vector frame (M, (e_1,e_2))}

    The expansion of the tensor field in a given frame is displayed via the 
    method :meth:`view` (the symbol * stands for tensor product)::
    
        sage: t.view()  # expansion in the manifold's default frame
        T = d/dx*dx - x d/dx*dy + x*y d/dy*dx + 2 d/dy*dy
        sage: t.view(e)
        T = (x + y) e_1*e^1 + y e_2*e^1 - 3*x e_2*e^2

    A tensor field acts as a multilinear map on 1-forms and vector fields; 
    in the present case, T being of type (1,1), it acts on pairs 
    (1-form, vector)::
    
        sage: a = M.one_form('a')
        sage: a[:] = (1, x)
        sage: v = M.vector_field('V')
        sage: v[:] = (y, 2)
        sage: t(a,v)
        scalar field 'T(a,V)' on the 2-dimensional manifold 'M'
        sage: t(a,v).expr()
        x^2*y^2 + 2*x + y
        sage: latex(t(a,v))
        T\left(a,V\right)
    
    Check by means of the component expression of t(a,v)::
    
        sage: t[1,1]*a[1]*v[1] + t[1,2]*a[1]*v[2] + t[2,1]*a[2]*v[1] + t[2,2]*a[2]*v[2] - t(a,v).expr()
        0

    A scalar field (rank-0 tensor field)::
    
        sage: f = M.scalar_field(x*y + 2, name='f') ; f 
        scalar field 'f' on the 2-dimensional manifold 'M'
        sage: f.tensor_type
        (0, 0)
        
    A scalar field acts on points on the manifold::
    
        sage: p = M.point((1,2))
        sage: f(p)
        4
        
    A vector field (rank-1 contravariant tensor field)::
    
        sage: v = M.vector_field('v') ; v
        vector field 'v' on the 2-dimensional manifold 'M'
        sage: v.tensor_type
        (1, 0)
        sage: v[1], v[2] = -x, y
        sage: v.view()
        v = -x d/dx + y d/dy        

    A field of symmetric bilinear forms::
    
        sage: q = M.sym_bilin_form_field('Q') ; q
        field of symmetric bilinear forms 'Q' on the 2-dimensional manifold 'M'
        sage: q.tensor_type
        (0, 2)

    The components of a symmetric bilinear form are dealt by the subclass 
    :class:`~sage.tensor.modules.comp.CompFullySym` of the class 
    :class:`~sage.tensor.modules.comp.Components`, which takes into 
    account the symmetry between the two indices::
    
        sage: q[1,1], q[1,2], q[2,2] = (0, -x, y) # no need to set the component (2,1)
        sage: print type(q.comp())
        <class 'sage.tensor.modules.comp.CompFullySym'>
        sage: q[:] # note that the component (2,1) is equal to the component (1,2)
        [ 0 -x]
        [-x  y]
        sage: q.view()
        Q = -x dx*dy - x dy*dx + y dy*dy
    
    Internally (dictionary :attr:`_comp` of the class 
    :class:`~sage.tensor.modules.comp.Components`), only
    the non-zero and non-redundant components are stored::
    
        sage: q.comp()._comp
        {(1, 2): scalar field on the 2-dimensional manifold 'M',
        (2, 2): scalar field on the 2-dimensional manifold 'M'}
        sage: q.comp()._comp[(1,2)].expr()
        -x
        sage: q.comp()._comp[(2,2)].expr()
        y

    More generally, tensor symmetries or antisymmetries can be specified via
    the keywords ``sym`` and ``antisym``. For instance a rank-4 covariant 
    tensor symmetric with respect to its first two arguments and 
    antisymmetric with respect to its last two ones is declared as follows::
    
        sage: t = M.tensor_field(0, 4, 'T', sym=(0,1), antisym=(2,3))
        sage: t[1,2,1,2] = 3
        sage: t[2,1,1,2] # check of the symmetry with respect to the first 2 indices
        3
        sage: t[1,2,2,1] # check of the antisymmetry with respect to the last 2 indices
        -3

"""

#******************************************************************************
#       Copyright (C) 2013, 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013, 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.rings.integer import Integer
from sage.structure.element import ModuleElement  
from sage.tensor.modules.free_module_tensor import FreeModuleTensor

class TensorField(ModuleElement):
    r"""
    Base class for tensor fields on an open set of a differentiable manifold, 
    with values on an open subset of a differentiable manifold. 
    
    An instance of this class is a tensor field along an open subset `U` 
    of some manifold `S` with values in an open subset `V` 
    of a manifold `M`, via a differentiable mapping `\Phi: U \rightarrow V`. 
    The standard case of a tensor field *on* a manifold corresponds to `S=M`, 
    `U=V` and `\Phi = \mathrm{Id}`. Another common case is `\Phi` being an
    immersion.

    If `V` is parallelizable, the class :class:`TensorFieldParal` should be
    used instead. 
    
    A tensor field of type `(k,\ell)` is a field `t` on `U`, so that at each 
    point `p\in U`, `t(p)` is a multilinear map of the type:

    .. MATH::

        t(p):\ \underbrace{T_p^*M\times\cdots\times T_p^*M}_{k\ \; \mbox{times}}
        \times \underbrace{T_p M\times\cdots\times T_p M}_{\ell\ \; \mbox{times}}
        \longrightarrow \RR
    
    where `T_p M` stands for the tangent space at the point `p` on the
    manifold `M` and `T_p^* M` for its dual vector space. The integer `k+\ell`
    is called the tensor rank. 
    
    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector 
      fields along `U` with values on `\Phi(U)\subset V \subset M`
    - ``tensor_type`` -- pair (k,l) with k being the contravariant rank and l 
      the covariant rank
    - ``name`` -- (default: None) name given to the tensor field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor field; 
      if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: None) a symmetry or a list of symmetries among the 
      tensor arguments: each symmetry is described by a tuple containing 
      the positions of the involved arguments, with the convention position=0
      for the first argument. For instance:
        * sym=(0,1) for a symmetry between the 1st and 2nd arguments 
        * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
          arguments and a symmetry between the 2nd, 4th and 5th arguments.
    - ``antisym`` -- (default: None) antisymmetry or list of antisymmetries 
      among the arguments, with the same convention as for ``sym``. 


    EXAMPLES:

    """
    def __init__(self, vector_field_module, tensor_type, name=None, 
                 latex_name=None, sym=None, antisym=None):
        ModuleElement.__init__(self, 
                               vector_field_module.tensor_module(*tensor_type))
        self.vmodule = vector_field_module
        self.tensor_type = tuple(tensor_type)
        self.tensor_rank = self.tensor_type[0] + self.tensor_type[1]
        self.name = name
        if latex_name is None:
            self.latex_name = self.name
        else:
            self.latex_name = latex_name
        self.domain = vector_field_module.domain
        self.ambient_domain = vector_field_module.ambient_domain
        self.restrictions = {} # dict. of restrictions on subdomains of self.domain
        # Treatment of symmetry declarations:
        self.sym = []
        if sym is not None and sym != []:
            if isinstance(sym[0], (int, Integer)):  
                # a single symmetry is provided as a tuple -> 1-item list:
                sym = [tuple(sym)]
            for isym in sym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self.tensor_rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                 " not in [0," + str(self.tensor_rank-1) + "]")
                    self.sym.append(tuple(isym))       
        self.antisym = []
        if antisym is not None and antisym != []:
            if isinstance(antisym[0], (int, Integer)):  
                # a single antisymmetry is provided as a tuple -> 1-item list:
                antisym = [tuple(antisym)]
            for isym in antisym:
                if len(isym) > 1:
                    for i in isym:
                        if i<0 or i>self.tensor_rank-1:
                            raise IndexError("Invalid position: " + str(i) +
                                " not in [0," + str(self.tensor_rank-1) + "]")
                    self.antisym.append(tuple(isym))
        # Final consistency check:
        index_list = []
        for isym in self.sym:
            index_list += isym
        for isym in self.antisym:
            index_list += isym
        if len(index_list) != len(set(index_list)):
            # There is a repeated index position:
            raise IndexError("Incompatible lists of symmetries: the same " + 
                             "position appears more than once.")
        # Initialization of derived quantities:
        self._init_derived() 

    ####### Required methods for ModuleElement (beside arithmetic) #######
    
    def __nonzero__(self):
        r"""
        Return True if ``self`` is nonzero and False otherwise. 
        
        This method is called by self.is_zero(). 
        
        EXAMPLE:
        
        Tensor field defined by parts on a 2-dimensional manifold::
        
            sage: M = Manifold(2, 'M')
            sage: U = M.open_domain('U')
            sage: c_xy.<x, y> = U.chart('x y')
            sage: V = M.open_domain('V')
            sage: c_uv.<u, v> = V.chart('u v')
            sage: t = M.tensor_field(1, 2, name='t')
            sage: tu = U.tensor_field(1, 2, name='t')
            sage: t.set_restriction(tu)
            sage: tv = V.tensor_field(1, 2, name='t')
            sage: t.set_restriction(tv)
            sage: tu[0,0,0] = 0
            sage: tv[0,0,0] = 0
            sage: t.is_zero()
            True
            sage: tv[0,0,0] = 1
            sage: t.is_zero()
            False

        """
        resu = False
        for rst in self.restrictions.values():
            resu = resu or rst.__nonzero__()
        return resu
                
    ####### End of required methods for ModuleElement (beside arithmetic) #######

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "tensor field"
        if self.name is not None:
            description += " '%s'" % self.name
        description += " of type (%s,%s) " % (str(self.tensor_type[0]), 
                                             str(self.tensor_type[1]))
        return self._final_repr(description)

    def _latex_(self):
        r"""
        LaTeX representation of the object.
        """
        if self.latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self.latex_name

    def _new_instance(self):
        r"""
        Create a :class:`TensorField` instance of the same tensor type, 
        with the same symmetries and on the same domain.

        This method must be redefined by derived classes of 
        :class:`TensorField`.
        
        """
        return TensorField(self.vmodule, self.tensor_type, sym=self.sym, 
                           antisym=self.antisym)

    def _final_repr(self, description):
        r"""
        Part of string representation common to all tensor fields.
        """
        if self.domain == self.ambient_domain:
            description += "on the " + str(self.domain)
        else:
            description += "along the " + str(self.domain) + \
                           " with values on the " + str(self.ambient_domain)
        return description

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        self._lie_derivatives = {} # collection of Lie derivatives of self

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        # First deletes any reference to self in the vectors' dictionary:
        for vid, val in self._lie_derivatives.items():
            del val[0]._lie_der_along_self[id(self)]
        self._lie_derivatives.clear()

    def set_restriction(self, rst):
        r"""
        Define a restriction of ``self`` to some subdomain.

        INPUT:
        
        - ``rst`` -- tensor field of the same type and symmetries as ``self``, 
          defined on a subdomain of ``self.domain`` (must be an instance of
          :class:`TensorField`)
          
        """
        if not isinstance(rst, TensorField):
            raise TypeError("The argument must be a tensor field.")
        if not rst.domain.is_subdomain(self.domain):
            raise ValueError("The domain of the declared restriction is not " +
                             "a subdomain of the current field's domain.")
        if not rst.ambient_domain.is_subdomain(self.ambient_domain):
            raise ValueError("The ambient domain of the declared " + 
                             "restriction is not a subdomain of the current " + 
                             "field's ambient domain.")
        if rst.tensor_type != self.tensor_type:
            raise ValueError("The declared restriction has not the same " + 
                             "tensor type as the current tensor field.")
        if rst.tensor_type != self.tensor_type:
            raise ValueError("The declared restriction has not the same " + 
                             "tensor type as the current tensor field.")
        if rst.sym != self.sym:
            raise ValueError("The declared restriction has not the same " +
                             "symmetries as the current tensor field.")
        if rst.antisym != self.antisym:
            raise ValueError("The declared restriction has not the same " +
                             "antisymmetries as the current tensor field.")
        self.restrictions[rst.domain] = rst

    def restriction(self, subdomain, dest_map=None):
        r"""
        Return the restriction of ``self`` to some subdomain.
        
        If such restriction has not been defined yet, it is constructed here.

        INPUT:
        
        - ``subdomain`` -- open subset `U` of ``self.domain`` (must be an instance 
          of :class:`~sage.geometry.manifolds.domain.OpenDomain`)
        - ``dest_map`` -- (default: None) destination map 
          `\Phi:\ U \rightarrow V`, where `V` is a subdomain of 
          ``self.codomain``
          (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`)
          If None, the restriction of ``self.vmodule.dest_map`` to `U` is 
          used.
          
        OUTPUT:
        
        - instance of :class:`TensorField` representing the restriction.
          
        """
        if subdomain == self.domain:
            return self
        if subdomain not in self.restrictions:
            if not subdomain.is_subdomain(self.domain):
                raise ValueError("The provided domain is not a subdomain of " + 
                                 "the current field's domain.")
            if dest_map is None:
                if self.vmodule.dest_map is not None:
                    dest_map = self.vmodule.dest_map.restriction(subdomain)
            elif not dest_map.codomain.is_subdomain(self.ambient_domain):
                raise ValueError("Argument dest_map not compatible with " + 
                                 "self.ambient_domain")
            smodule = subdomain.vector_field_module(dest_map=dest_map)
            self.restrictions[subdomain] = smodule.tensor(self.tensor_type, 
                                                    name=self.name, 
                                                    latex_name=self.latex_name, 
                                                    sym=self.sym, 
                                                    antisym=self.antisym)
        return self.restrictions[subdomain]

    def set_comp(self, basis=None):
        r"""
        Return the components of the tensor field in a given vector frame 
        for assignment.
        
        The components with respect to other frames on the same domain are 
        deleted, in order to avoid any inconsistency. To keep them, use the 
        method :meth:`add_comp` instead.
        
        INPUT:
        
        - ``basis`` -- (default: None) vector frame in which the components are
          defined; if none is provided, the components are assumed to refer to 
          the tensor field domain's default frame.
         
        OUTPUT: 
        
        - components in the given frame, as an instance of the 
          class :class:`~sage.tensor.modules.comp.Components`; if such 
          components did not exist previously, they are created.  
        
        EXAMPLES:
        
          

        """
        if self is self.parent().zero(): #!# this is maybe not very efficient
            raise ValueError("The zero tensor field cannot be changed.")
        if basis is None: 
            basis = self.domain.def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restriction(basis.domain, dest_map=basis.dest_map)
        return rst.set_comp(basis)

    def add_comp(self, basis=None):
        r"""
        Return the components of the tensor field in a given vector frame 
        for assignment.
        
        The components with respect to other frames on the same domain are 
        kept. To delete them them, use the method :meth:`set_comp` instead.
        
        INPUT:
        
        - ``basis`` -- (default: None) vector frame in which the components are
          defined; if none is provided, the components are assumed to refer to 
          the tensor field domain's default frame.
         
        OUTPUT: 
        
        - components in the given frame, as an instance of the 
          class :class:`~sage.tensor.modules.comp.Components`; if such 
          components did not exist previously, they are created.  
        
        EXAMPLES:
        
          

        """
        if self is self.parent().zero(): #!# this is maybe not very efficient
            raise ValueError("The zero tensor field cannot be changed.")
        if basis is None: 
            basis = self.domain.def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restriction(basis.domain, dest_map=basis.dest_map)
        return rst.add_comp(basis)


    def comp(self, basis=None, from_basis=None):
        r"""
        Return the components in a given vector frame.
        
        If the components are not known already, they are computed by the tensor
        change-of-basis formula from components in another vector frame. 
        
        INPUT:
        
        - ``basis`` -- (default: None) vector frame in which the components are 
          required; if none is provided, the components are assumed to refer to
          the tensor field domain's default frame
        - ``from_basis`` -- (default: None) vector frame from which the
          required components are computed, via the tensor change-of-basis 
          formula, if they are not known already in the basis ``basis``
          
        OUTPUT: 
        
        - components in the vector frame ``basis``, as an instance of the 
          class :class:`~sage.tensor.modules.comp.Components` 

        """
        if basis is None: 
            basis = self.domain.def_frame
        self._del_derived() # deletes the derived quantities
        rst = self.restriction(basis.domain, dest_map=basis.dest_map)
        return rst.comp(basis=basis, from_basis=from_basis)

    def __getitem__(self, args):
        r"""
        Return a component w.r.t. some frame.

        INPUT:
        
        - ``args`` -- list of indices defining the component; if [:] is 
          provided, all the components are returned. The frame can be passed
          as the first item of ``args``; if not, the default frame of the 
          tensor field's domain is assumed. 
    
        """
        if isinstance(args, list):  # case of [[...]] syntax
            if not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self.domain.def_frame
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self.domain.def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self.domain.def_frame
        return self.comp(frame)[args]

    def __setitem__(self, args, value):
        r"""
        Sets a component w.r.t to some vector frame. 

        INPUT:
        
       - ``args`` -- list of indices; if [:] is provided, all the components 
          are set. The frame can be passed as the first item of ``args``; if 
          not, the default frame of the tensor field's domain is assumed. 
        - ``value`` -- the value to be set or a list of values if ``args``
          == ``[:]``
    
        """
        if isinstance(args, list):  # case of [[...]] syntax
            if not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self.domain.def_frame
        else:
            if isinstance(args, (int, Integer, slice)):
                frame = self.domain.def_frame
            elif not isinstance(args[0], (int, Integer, slice)):
                frame = args[0]
                args = args[1:]
            else:
                frame = self.domain.def_frame
        self.set_comp(frame)[args] = value


    def copy(self):
        r"""
        Return an exact copy of ``self``.
        
        The name and the derived quantities are not copied. 
        
        EXAMPLES:
        

        # If the original tensor field is modified, the copy is not::
        

        """
        resu = self._new_instance()
        for dom, rst in self.restrictions.items():
            resu.restrictions[dom] = rst.copy()
        return resu

    def _common_subdomains(self, other):
        r"""
        Returns the list of subdomains of self.domain on which both ``self``
        and ``other`` have known restrictions.
        """
        resu = []
        for dom in self.restrictions:
            if dom in other.restrictions:
                resu.append(dom)
        return resu

    def __eq__(self, other):
        r"""
        Comparison (equality) operator. 
        
        INPUT:
        
        - ``other`` -- a tensor field or 0
        
        OUTPUT:
        
        - True if ``self`` is equal to ``other`` and False otherwise
        
        """
        if isinstance(other, (int, Integer)): # other should be 0
            if other == 0:
                return self.is_zero()
            else:
                return False
        elif not isinstance(other, TensorField):
            return False
        else: # other is another tensor field
            if other.vmodule != self.vmodule:
                return False
            if other.tensor_type != self.tensor_type:
                return False
            resu = True
            for dom, rst in self.restrictions.items():
                if dom in other.restrictions:
                    resu = resu and bool(rst == other.restrictions[dom])
            return resu

    def __ne__(self, other):
        r"""
        Inequality operator. 
        
        INPUT:
        
        - ``other`` -- a tensor field or 0
        
        OUTPUT:
        
        - True if ``self`` is different from ``other`` and False otherwise
        
        """
        return not self.__eq__(other)

    def __pos__(self):
        r"""
        Unary plus operator. 
        
        OUTPUT:
        
        - an exact copy of ``self``
    
        """
        resu = self._new_instance()
        for dom, rst in self.restrictions.items():
            resu.restrictions[dom] = + rst
        if self.name is not None:
            resu.name = '+' + self.name 
        if self.latex_name is not None:
            resu.latex_name = '+' + self.latex_name
        return resu

    def __neg__(self):
        r"""
        Unary minus operator. 
        
        OUTPUT:
        
        - the tensor field `-T`, where `T` is ``self``
    
        """
        resu = self._new_instance()
        for dom, rst in self.restrictions.items():
            resu.restrictions[dom] = - rst
        if self.name is not None:
            resu.name = '-' + self.name 
        if self.latex_name is not None:
            resu.latex_name = '-' + self.latex_name
        return resu

    ######### ModuleElement arithmetic operators ########
    
    def _add_(self, other):
        r"""
        Tensor field addition. 
        
        INPUT:
        
        - ``other`` -- a tensor field, in the same tensor module as ``self`` 
        
        OUPUT:
        
        - the tensor field resulting from the addition of ``self`` and ``other``
        
        """
        if other == 0:
            return +self
        if not isinstance(other, TensorField):
            raise TypeError("For the addition, other must be a tensor field.")
        if other.vmodule != self.vmodule:
            raise ValueError("The two tensor fields do not belong to the " + 
                             "same module.")
        if other.tensor_type != self.tensor_type:
            raise TypeError("The two tensor fields are not of the same type.")
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self.restrictions[dom] + other.restrictions[dom]
        some_rst = resu_rst.values()[0]
        resu_sym = some_rst.sym
        resu_antisym = some_rst.antisym
        resu = self.vmodule.tensor(self.tensor_type, sym=resu_sym, 
                                   antisym=resu_antisym)
        resu.restrictions = resu_rst
        if self.name is not None and other.name is not None:
            resu.name = self.name + '+' + other.name
        if self.latex_name is not None and other.latex_name is not None:
            resu.latex_name = self.latex_name + '+' + other.latex_name
        return resu

    def _sub_(self, other):
        r"""
        Tensor field subtraction. 
        
        INPUT:
        
        - ``other`` -- a tensor field, in the same tensor module as ``self``
        
        OUPUT:
        
        - the tensor field resulting from the subtraction of ``other`` from 
          ``self``
                
        """
        if other == 0:
            return +self
        if not isinstance(other, TensorField):
            raise TypeError("For the subtraction, other must be a tensor " + 
                            "field.")
        if other.vmodule != self.vmodule:
            raise ValueError("The two tensor fields do not belong to the " + 
                             "same module.")
        if other.tensor_type != self.tensor_type:
            raise TypeError("The two tensor fields are not of the same type.")
        resu_rst = {}
        for dom in self._common_subdomains(other):
            resu_rst[dom] = self.restrictions[dom] - other.restrictions[dom]
        some_rst = resu_rst.values()[0]
        resu_sym = some_rst.sym
        resu_antisym = some_rst.antisym
        resu = self.vmodule.tensor(self.tensor_type, sym=resu_sym, 
                                   antisym=resu_antisym)
        resu.restrictions = resu_rst
        if self.name is not None and other.name is not None:
            resu.name = self.name + '-' + other.name
        if self.latex_name is not None and other.latex_name is not None:
            resu.latex_name = self.latex_name + '-' + other.latex_name
        return resu

    def _rmul_(self, other):
        r"""
        Multiplication on the left by a scalar field (``other``)
        
        """
        #!# The following test is probably not necessary:
        if isinstance(other, TensorField):
            raise NotImplementedError("Left tensor product not implemented.")
        # Left multiplication by a scalar field: 
        resu = self._new_instance()
        for dom, rst in self.restrictions.items():
            resu.restrictions[dom] = other.restriction(dom) * rst
        return resu

    ######### End of ModuleElement arithmetic operators ########

    def __radd__(self, other):
        r"""
        Addition on the left with ``other``. 
        
        This allows to write "0 + t", where "t" is a tensor
        
        """
        return self.__add__(other)

    def __rsub__(self, other):
        r"""
        Subtraction from ``other``. 

        This allows to write "0 - t", where "t" is a tensor
        
        """
        return (-self).__add__(other)


    def up(self, metric, pos=None):
        r"""
        Compute a metric dual by raising some index with a given metric.
        
        If ``self`` is a tensor field `T` of type `(k,\ell)` and `p` is the 
        position of a covariant index (i.e. `k\leq p < k+\ell`), 
        the output with ``pos`` `=p` is the tensor field `T^\sharp` of type 
        `(k+1,\ell-1)` whose components are

        .. MATH::

            (T^\sharp)^{a_1\ldots a_{k+1}}_{\qquad\quad b_1 \ldots b_{\ell-1}}
            = g^{a_{k+1} i} \, 
            T^{a_1\ldots a_k}_{\qquad\ \  b_1 \ldots b_{p-k} \, i \, b_{p-k+1} \ldots b_{\ell-1}},
            
        `g^{ab}` being the components of the inverse metric. 

        The reverse operation is :meth:`TensorFieldParal.down`

        INPUT:
        
        - ``metric`` -- metric `g`, as an instance of 
          :class:`~sage.geometry.manifolds.metric.Metric`
        - ``pos`` -- (default: None) position of the index (with the
          convention ``pos=0`` for the first index); if none, the raising is 
          performed over all the covariant indices, starting from the first one
         
        OUTPUT:
        
        - the tensor field `T^\sharp` resulting from the index raising operation

        EXAMPLES:
        
        Raising the index of a 1-form results in a vector field::
        
            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart('x y')
            sage: g = M.metric('g')
            sage: g[1,1], g[1,2], g[2,2] = 1+x, x*y, 1-y
            sage: w = M.one_form()
            sage: w[:] = [-1, 2]
            sage: v = w.up(g) ; v
            vector field on the 2-dimensional manifold 'M'
            sage: v.view()
            ((2*x - 1)*y + 1)/(x^2*y^2 + (x + 1)*y - x - 1) d/dx - (x*y + 2*x + 2)/(x^2*y^2 + (x + 1)*y - x - 1) d/dy
            sage: g.inverse()[:]
            [ (y - 1)/(x^2*y^2 + (x + 1)*y - x - 1)      x*y/(x^2*y^2 + (x + 1)*y - x - 1)]
            [     x*y/(x^2*y^2 + (x + 1)*y - x - 1) -(x + 1)/(x^2*y^2 + (x + 1)*y - x - 1)]
            sage: w1 = v.down(g) ; w1   # the reverse operation
            1-form on the 2-dimensional manifold 'M'
            sage: w1.view()
            -dx + 2 dy
            sage: w1 == w
            True

        Raising the indices of a tensor field of type (0,2)::

            sage: t = M.tensor_field(0, 2)
            sage: t[:] = [[1,2], [3,4]]
            sage: tu0 = t.up(g, 0) ; tu0  # raising the first index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: tu0[:]
            [  ((3*x + 1)*y - 1)/(x^2*y^2 + (x + 1)*y - x - 1) 2*((2*x + 1)*y - 1)/(x^2*y^2 + (x + 1)*y - x - 1)]
            [    (x*y - 3*x - 3)/(x^2*y^2 + (x + 1)*y - x - 1)   2*(x*y - 2*x - 2)/(x^2*y^2 + (x + 1)*y - x - 1)]
            sage: tuu0 = tu0.up(g) ; tuu0 # the two indices have been raised, starting from the first one
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: tu1 = t.up(g, 1) ; tu1 # raising the second index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: tu1[:]
            [((2*x + 1)*y - 1)/(x^2*y^2 + (x + 1)*y - x - 1) ((4*x + 3)*y - 3)/(x^2*y^2 + (x + 1)*y - x - 1)]
            [  (x*y - 2*x - 2)/(x^2*y^2 + (x + 1)*y - x - 1) (3*x*y - 4*x - 4)/(x^2*y^2 + (x + 1)*y - x - 1)]
            sage: tuu1 = tu1.up(g) ; tuu1 # the two indices have been raised, starting from the second one
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: tuu0 == tuu1 # the order of index raising is important
            False
            sage: tuu = t.up(g) ; tuu # both indices are raised, starting from the first one
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: tuu0 == tuu # the same order for index raising has been applied
            True
            sage: tuu1 == tuu # to get tuu1, indices have been raised from the last one, contrary to tuu 
            False
            sage: d0tuu = tuu.down(g, 0) ; d0tuu # the first index is lowered again
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: dd0tuu = d0tuu.down(g) ; dd0tuu  # the second index is then lowered
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: d1tuu = tuu.down(g, 1) ; d1tuu # lowering operation, starting from the last index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: dd1tuu = d1tuu.down(g) ; dd1tuu
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: ddtuu = tuu.down(g) ; ddtuu # both indices are lowered, starting from the last one
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: ddtuu == t # should be true
            True
            sage: dd0tuu == t # not true, because of the order of index lowering to get dd0tuu
            False
            sage: dd1tuu == t # should be true
            True

        """
        from sage.rings.integer import Integer
        n_con = self.tensor_type[0] # number of contravariant indices = k 
        if pos is None:
            result = self
            for p in range(n_con, self.tensor_rank):
                k = result.tensor_type[0]
                result = result.up(metric, k)
            return result
        if not isinstance(pos, (int, Integer)):
            raise TypeError("The argument 'pos' must be an integer.")
        if pos<n_con or pos>self.tensor_rank-1:
            print "pos = ", pos
            raise ValueError("Position out of range.")
        return self.contract(pos, metric.inverse(), 0)

    def down(self, metric, pos=None):
        r"""
        Compute a metric dual by lowering some index with a given metric.
        
        If ``self`` is a tensor field `T` of type `(k,\ell)` and `p` is the 
        position of a contravariant index (i.e. `0\leq p < k`), the output with
        ``pos`` `=p` is the tensor field `T^\flat` of type `(k-1,\ell+1)` whose 
        components are

        .. MATH::

            (T^\flat)^{a_1\ldots a_{k-1}}_{\qquad\quad b_1 \ldots b_{\ell+1}}
            = g_{b_1 i} \, 
            T^{a_1\ldots a_{p} \, i \, a_{p+1}\ldots a_{k-1}}_{\qquad\qquad\qquad\quad b_2 \ldots b_{\ell+1}},
            
        `g_{ab}` being the components of the metric tensor. 

        The reverse operation is :meth:`TensorFieldParal.up`
        
        INPUT:
        
        - ``metric`` -- metric `g`, as an instance of 
          :class:`~sage.geometry.manifolds.metric.Metric`
        - ``pos`` -- (default: None) position of the index (with the 
          convention ``pos=0`` for the first index); if none, the lowering is 
          performed over all the contravariant indices, starting from the last 
          one
         
        OUTPUT:
        
        - the tensor field `T^\flat` resulting from the index lowering operation
        
        EXAMPLES:
        
        Lowering the index of a vector field results in a 1-form::
        
            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart('x y')
            sage: g = M.metric('g')
            sage: g[1,1], g[1,2], g[2,2] = 1+x, x*y, 1-y
            sage: v = M.vector_field()
            sage: v[:] = [-1,2]
            sage: w = v.down(g) ; w
            1-form on the 2-dimensional manifold 'M'
            sage: w.view()
            (2*x*y - x - 1) dx + (-(x + 2)*y + 2) dy
            sage: v1 = w.up(g) ; v1  # the reverse operation
            vector field on the 2-dimensional manifold 'M'
            sage: v1 == v
            True

        Lowering the indices of a tensor field of type (2,0)::
        
            sage: t = M.tensor_field(2, 0)
            sage: t[:] = [[1,2], [3,4]]
            sage: td0 = t.down(g, 0) ; td0  # lowering the first index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: td0[:]
            [  3*x*y + x + 1   (x - 3)*y + 3]
            [4*x*y + 2*x + 2 2*(x - 2)*y + 4]
            sage: tdd0 = td0.down(g) ; tdd0 # the two indices have been lowered, starting from the first one
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: tdd0[:]
            [      4*x^2*y^2 + x^2 + 5*(x^2 + x)*y + 2*x + 1 2*(x^2 - 2*x)*y^2 + (x^2 + 2*x - 3)*y + 3*x + 3]
            [(3*x^2 - 4*x)*y^2 + (x^2 + 3*x - 2)*y + 2*x + 2           (x^2 - 5*x + 4)*y^2 + (5*x - 8)*y + 4]
            sage: td1 = t.down(g, 1) ; td1  # lowering the second index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: td1[:]
            [  2*x*y + x + 1   (x - 2)*y + 2]
            [4*x*y + 3*x + 3 (3*x - 4)*y + 4]
            sage: tdd1 = td1.down(g) ; tdd1 # the two indices have been lowered, starting from the second one
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: tdd1[:]
            [      4*x^2*y^2 + x^2 + 5*(x^2 + x)*y + 2*x + 1 (3*x^2 - 4*x)*y^2 + (x^2 + 3*x - 2)*y + 2*x + 2]
            [2*(x^2 - 2*x)*y^2 + (x^2 + 2*x - 3)*y + 3*x + 3           (x^2 - 5*x + 4)*y^2 + (5*x - 8)*y + 4]
            sage: tdd1 == tdd0   # the order of index lowering is important
            False
            sage: tdd = t.down(g) ; tdd  # both indices are lowered, starting from the last one
            tensor field of type (0,2) on the 2-dimensional manifold 'M'
            sage: tdd[:]
            [      4*x^2*y^2 + x^2 + 5*(x^2 + x)*y + 2*x + 1 (3*x^2 - 4*x)*y^2 + (x^2 + 3*x - 2)*y + 2*x + 2]
            [2*(x^2 - 2*x)*y^2 + (x^2 + 2*x - 3)*y + 3*x + 3           (x^2 - 5*x + 4)*y^2 + (5*x - 8)*y + 4]
            sage: tdd0 == tdd  # to get tdd0, indices have been lowered from the first one, contrary to tdd 
            False
            sage: tdd1 == tdd  # the same order for index lowering has been applied
            True
            sage: u0tdd = tdd.up(g, 0) ; u0tdd # the first index is raised again
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: uu0tdd = u0tdd.up(g) ; uu0tdd # the second index is then raised
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: u1tdd = tdd.up(g, 1) ; u1tdd  # raising operation, starting from the last index
            field of endomorphisms on the 2-dimensional manifold 'M'
            sage: uu1tdd = u1tdd.up(g) ; uu1tdd
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: uutdd = tdd.up(g) ; uutdd  # both indices are raised, starting from the first one
            tensor field of type (2,0) on the 2-dimensional manifold 'M'
            sage: uutdd == t  # should be true
            True
            sage: uu0tdd == t # should be true
            True
            sage: uu1tdd == t # not true, because of the order of index raising to get uu1tdd
            False
 
        """
        from sage.rings.integer import Integer
        n_con = self.tensor_type[0] # number of contravariant indices = k 
        if pos is None:
            result = self
            for p in range(0, n_con):
                k = result.tensor_type[0]
                result = result.down(metric, k-1)
            return result
        if not isinstance(pos, (int, Integer)):
            raise TypeError("The argument 'pos' must be an integer.")
        if pos<0 or pos>=n_con:
            print "pos = ", pos
            raise ValueError("Position out of range.")
        return metric.contract(1, self, pos)


#******************************************************************************

class TensorFieldParal(FreeModuleTensor, TensorField):
    r"""
    Base class for tensor fields on an open set of a differentiable manifold, 
    with values on parallelizable open subset of a differentiable manifold. 
    
    An instance of this class is a tensor field along an open subset `U` 
    of some manifold `S` with values in a parallelizable open subset `V` 
    of a manifold `M`, via a differentiable mapping `\Phi: U \rightarrow V`. 
    The standard case of a tensor field *on* a manifold corresponds to `S=M`, 
    `U=V` and `\Phi = \mathrm{Id}`. Another common case is `\Phi` being an
    immersion.

    A tensor field of type `(k,\ell)` is a field `t` on `U`, so that at each 
    point `p\in U`, `t(p)` is a multilinear map of the type:

    .. MATH::

        t(p):\ \underbrace{T_p^*M\times\cdots\times T_p^*M}_{k\ \; \mbox{times}}
        \times \underbrace{T_p M\times\cdots\times T_p M}_{\ell\ \; \mbox{times}}
        \longrightarrow \RR
    
    where `T_p M` stands for the tangent space at the point `p` on the
    manifold `M` and `T_p^* M` for its dual vector space. The integer `k+\ell`
    is called the tensor rank. 
    
    INPUT:
    
    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector 
      fields along `U` with values on `\Phi(U)\subset V \subset M`
    - ``tensor_type`` -- pair (k,l) with k being the contravariant rank and l 
      the covariant rank
    - ``name`` -- (default: None) name given to the tensor field
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the tensor field; 
      if none is provided, the LaTeX symbol is set to ``name``
    - ``sym`` -- (default: None) a symmetry or a list of symmetries among the 
      tensor arguments: each symmetry is described by a tuple containing 
      the positions of the involved arguments, with the convention position=0
      for the first argument. For instance:
        * sym=(0,1) for a symmetry between the 1st and 2nd arguments 
        * sym=[(0,2),(1,3,4)] for a symmetry between the 1st and 3rd
          arguments and a symmetry between the 2nd, 4th and 5th arguments.
    - ``antisym`` -- (default: None) antisymmetry or list of antisymmetries 
      among the arguments, with the same convention as for ``sym``. 


    EXAMPLES:

    A tensor field of type (2,0) on a 3-dimensional manifold::
    
        sage: M = Manifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart('x y z')
        sage: t = M.tensor_field(2, 0, 'T') ; t
        tensor field 'T' of type (2,0) on the 3-dimensional manifold 'M'

    Tensor fields are considered as elements of a module over the ring 
    `C^\infty(M)` of scalar fields on `M`::
    
        sage: t.parent()
        free module TF^(2,0)(M) of type-(2,0) tensors fields on the 3-dimensional manifold 'M'
        sage: t.parent().base_ring()
        ring of scalar fields on the 3-dimensional manifold 'M'

    The components with respect to the manifold's default frame are set or read
    by means of square brackets::
    
        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: for i in M.irange():
        ...       for j in M.irange():
        ...           t[i,j] = (i+1)**(j+1)
        ...
        sage: [[ t[i,j] for j in M.irange()] for i in M.irange()]
        [[1, 1, 1], [2, 4, 8], [3, 9, 27]]
    
    A shortcut for the above is using [:]::
    
        sage: t[:]
        [ 1  1  1]
        [ 2  4  8]
        [ 3  9 27]

    The components with respect to another frame are set via the method
    :meth:`set_comp` and read via the method :meth:`comp`; both return an 
    instance of :class:`~sage.tensor.modules.comp.Components`::
    
        sage: f = M.vector_frame('f')  # a new frame defined on M, in addition to e
        sage: t.set_comp(f)[0,0] = -3
        sage: t.comp(f)
        2-indices components w.r.t. vector frame (M, (f_0,f_1,f_2))
        sage: t.comp(f)[0,0]
        -3
        sage: t.comp(f)[:]  # the full list of components
        [-3  0  0]
        [ 0  0  0]
        [ 0  0  0]

    To avoid any insconstency between the various components, the method 
    :meth:`set_comp` deletes the components in other frames. 
    Accordingly, the components in the frame e have been deleted::
    
        sage: t.components
        {vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. vector frame (M, (f_0,f_1,f_2))}

    To keep the other components, one must use the method :meth:`add_comp`::
    
        sage: t = M.tensor_field(2, 0, 'T')  # Let us restart
        sage: t[0,0] = 2                     # sets the components in the frame e
        sage: # We now set the components in the frame f with add_comp:
        sage: t.add_comp(f)[0,0] = -3
        sage: # The components w.r.t. frame e have been kept: 
        sage: t.components
        {vector frame (M, (e_0,e_1,e_2)): 2-indices components w.r.t. vector frame (M, (e_0,e_1,e_2)), vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. vector frame (M, (f_0,f_1,f_2))}

    The basic attributes of :class:`TensorField` are :attr:`domain`, 
    :attr:`tensor_type` (the pair (k,l)), :attr:`tensor_rank` (the integer k+l),  
    and :attr:`components` (the dictionary of the components w.r.t. various 
    frames)::

        sage: t.domain
        3-dimensional manifold 'M'
        sage: t.tensor_type
        (2, 0)
        sage: t.tensor_rank
        2
        sage: t.components
        {vector frame (M, (e_0,e_1,e_2)): 2-indices components w.r.t. vector frame (M, (e_0,e_1,e_2)), vector frame (M, (f_0,f_1,f_2)): 2-indices components w.r.t. vector frame (M, (f_0,f_1,f_2))}

    Symmetries and antisymmetries are declared via the keywords ``sym`` and
    ``antisym``. For instance, a rank-6 covariant tensor that is symmetric with
    respect to its 1st and 3rd arguments and antisymmetric with respect to the 
    2nd, 5th and 6th arguments is set up as follows::
    
        sage: a = M.tensor_field(0, 6, 'T', sym=(0,2), antisym=(1,4,5))
        sage: a[0,0,1,0,1,2] = 3
        sage: a[1,0,0,0,1,2] # check of the symmetry
        3
        sage: a[0,1,1,0,0,2], a[0,1,1,0,2,0] # check of the antisymmetry
        (-3, 3)  

    Multiple symmetries or antisymmetries are allowed; they must then be 
    declared as a list. For instance, a rank-4 covariant tensor that is 
    antisymmetric with respect to its 1st and 2nd arguments and with respect to
    its 3rd and 4th argument must be declared as::
    
        sage: r = M.tensor_field(0, 4, 'T', antisym=[(0,1), (2,3)])
        sage: r[0,1,2,0] = 3
        sage: r[1,0,2,0] # first antisymmetry
        -3
        sage: r[0,1,0,2] # second antisymmetry
        -3
        sage: r[1,0,0,2] # both antisymmetries acting
        3
    
    Tensor fields of the same type can be added and subtracted::
    
        sage: a = M.tensor_field(2, 0)
        sage: a[0,0], a[0,1], a[0,2] = (1,2,3)
        sage: b = M.tensor_field(2, 0)
        sage: b[0,0], b[1,1], b[2,2], b[0,2] = (4,5,6,7)
        sage: s = a + 2*b ; s
        tensor field of type (2,0) on the 3-dimensional manifold 'M'
        sage: a[:], (2*b)[:], s[:]
        (
        [1 2 3]  [ 8  0 14]  [ 9  2 17]
        [0 0 0]  [ 0 10  0]  [ 0 10  0]
        [0 0 0], [ 0  0 12], [ 0  0 12]
        )
        sage: s = a - b ; s
        tensor field of type (2,0) on the 3-dimensional manifold 'M'
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [-3  2 -4]
        [0 0 0]  [0 5 0]  [ 0 -5  0]
        [0 0 0], [0 0 6], [ 0  0 -6]
        )

    Symmetries are preserved by the addition whenever it is possible::
    
        sage: a = M.tensor_field(2, 0, sym=(0,1))
        sage: a[0,0], a[0,1], a[0,2] = (1,2,3)
        sage: s = a + b
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [ 5  2 10]
        [2 0 0]  [0 5 0]  [ 2  5  0]
        [3 0 0], [0 0 6], [ 3  0  6]
        )
        sage: a.symmetries()
        symmetry: (0, 1);  no antisymmetry
        sage: b.symmetries()
        no symmetry;  no antisymmetry
        sage: s.symmetries()
        no symmetry;  no antisymmetry
        sage: # let us now make b symmetric:
        sage: b = M.tensor_field(2, 0, sym=(0,1))
        sage: b[0,0], b[1,1], b[2,2], b[0,2] = (4,5,6,7)
        sage: s = a + b
        sage: a[:], b[:], s[:]
        (
        [1 2 3]  [4 0 7]  [ 5  2 10]
        [2 0 0]  [0 5 0]  [ 2  5  0]
        [3 0 0], [7 0 6], [10  0  6]
        )
        sage: s.symmetries()  # s is symmetric because both a and b are
        symmetry: (0, 1);  no antisymmetry

    The tensor product is taken with the operator \*::
    
        sage: c = a*b ; c
        tensor field of type (4,0) on the 3-dimensional manifold 'M'
        sage: c.symmetries()  # since a and b are both symmetric, a*b has two symmetries:
        symmetries: [(0, 1), (2, 3)];  no antisymmetry

    The tensor product of two fully contravariant tensors is not symmetric in 
    general::
    
        sage: a*b == b*a
        False

    The tensor product of a fully contravariant tensor by a fully covariant one
    is symmetric::
    
        sage: d = M.diff_form(2)  # a fully covariant tensor field
        sage: d[0,1], d[0,2], d[1,2] = (3, 2, 1)
        sage: s = a*d ; s 
        tensor field of type (2,2) on the 3-dimensional manifold 'M'
        sage: s.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: s1 = d*a ; s1 
        tensor field of type (2,2) on the 3-dimensional manifold 'M'
        sage: s1.symmetries()
        symmetry: (0, 1);  antisymmetry: (2, 3)
        sage: d*a == a*d
        True

    """
    def __init__(self, vector_field_module, tensor_type, name=None, 
                 latex_name=None, sym=None, antisym=None):
        FreeModuleTensor.__init__(self, vector_field_module, tensor_type,
                                  name=name, latex_name=latex_name,
                                  sym=sym, antisym=antisym)
        self.domain = vector_field_module.domain
        self.ambient_domain = vector_field_module.ambient_domain
        # TensorField attributes:
        self.vmodule = vector_field_module
        self.restrictions = {} # dict. of restrictions on subdomains of self.domain        
        # Initialization of derived quantities:
        self._init_derived() 

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return TensorField._repr_(self)
        
    def _new_instance(self):
        r"""
        Create a :class:`TensorFieldParal` instance of the same tensor type, 
        with the same symmetries and on the same domain.

        This method must be redefined by derived classes of 
        :class:`TensorFieldParal`.
        
        """
        return TensorFieldParal(self.fmodule, self.tensor_type, sym=self.sym, 
                                antisym=self.antisym)

    def _init_derived(self):
        r"""
        Initialize the derived quantities
        """
        FreeModuleTensor._init_derived(self) 
        TensorField._init_derived(self)

    def _del_derived(self):
        r"""
        Delete the derived quantities
        """
        FreeModuleTensor._del_derived(self) 
        TensorField._del_derived(self)

    def common_coord_frame(self, other):
        r"""
        Find a common coordinate frame for the components of ``self`` and 
        ``other``. 
        
        In case of multiple common bases, the domain's default coordinate
        basis is privileged. 
        If the current components of ``self`` and ``other`` are all relative to
        different frames, a common frame is searched by performing a component
        transformation, via the transformations listed in 
        ``self.domain.frame_changes``, still privileging transformations to 
        the domain's default frame.
        
        INPUT:
        
        - ``other`` -- a tensor field (instance of :class:`TensorFieldParal`)
        
        OUPUT:
        
        - common coordinate frame; if no common basis is found, None is 
          returned. 
        
        """
        from vectorframe import CoordFrame
        # Compatibility checks:
        if not isinstance(other, TensorFieldParal):
            raise TypeError("The argument must be of type TensorFieldParal.")
        dom = self.domain
        def_frame = dom.def_frame
        #
        # 1/ Search for a common frame among the existing components, i.e. 
        #    without performing any component transformation. 
        #    -------------------------------------------------------------
        # 1a/ Direct search
        if def_frame in self.components and \
           def_frame in other.components and \
           isinstance(dom.def_frame, CoordFrame):
            return def_frame # the domain's default frame is privileged
        for frame1 in self.components:
            if frame1 in other.components and \
               isinstance(frame1, CoordFrame):
                return frame1
        # 1b/ Search involving subframes
        dom2 = other.domain
        for frame1 in self.components:
            if not isinstance(frame1, CoordFrame):
                continue
            for frame2 in other.components:
                if not isinstance(frame2, CoordFrame):
                    continue
                if frame2 in frame1.subframes:
                    self.comp(frame2)
                    return frame2
                if frame1 in frame2.subframes:
                    other.comp(frame1)
                    return frame1
        #
        # 2/ Search for a common frame via one component transformation
        #    ----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at least 
        # one component transformation to get a common frame
        if isinstance(dom.def_frame, CoordFrame):
            if def_frame in self.components:
                for oframe in other.components:
                    if (oframe, def_frame) in dom.frame_changes:
                        other.comp(def_frame, from_frame=oframe)
                        return def_frame
            if def_frame in other.components:
                for sframe in self.components:
                    if (sframe, def_frame) in dom.frame_changes:
                        self.comp(def_frame, from_frame=sframe)
                        return def_frame
        # If this point is reached, then def_frame cannot be a common frame
        # via a single component transformation
        for sframe in self.components:
            if not instance(sframe, CoordFrame):
                continue
            for oframe in other.components:
                if not instance(oframe, CoordFrame):
                    continue
                if (oframe, sframe) in dom.frame_changes:
                    other.comp(sframe, from_frame=oframe)
                    return sframe
                if (sframe, oframe) in dom.frame_changes:
                    self.comp(oframe, from_frame=sframe)
                    return oframe
        #
        # 3/ Search for a common frame via two component transformations
        #    -----------------------------------------------------------
        # If this point is reached, it is indeed necessary to perform at two
        # component transformation to get a common frame
        for sframe in self.components:
            for oframe in other.components:
                if (sframe, def_frame) in dom.frame_changes and \
                   (oframe, def_frame) in dom.frame_changes and \
                   isinstance(def_frame, CoordFrame):
                    self.comp(def_frame, from_frame=sframe)
                    other.comp(def_frame, from_frame=oframe)
                    return def_frame
                for frame in dom.frames:
                    if (sframe, frame) in dom.frame_changes and \
                       (oframe, frame) in dom.frame_changes and \
                       isinstance(frame, CoordFrame):
                        self.comp(frame, from_frame=sframe)
                        other.comp(frame, from_frame=oframe)
                        return frame
        #
        # If this point is reached, no common frame could be found, even at 
        # the price of component transformations:
        return None
    

    def lie_der(self, vector):
        r"""
        Computes the Lie derivative with respect to a vector field.
        
        The Lie derivative is stored in the dictionary 
        :attr:`_lie_derivatives`, so that there is no need to 
        recompute it at the next call if neither ``self`` nor ``vector``
        have been modified meanwhile. 
        
        INPUT:
        
        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken
          
        OUTPUT:
        
        - the tensor field that is the Lie derivative of ``self`` with respect 
          to ``vector``
        
        EXAMPLES:
        
        Lie derivative of a vector::
        
            sage: M = Manifold(2, 'M', start_index=1)
            sage: c_xy.<x,y> = M.chart('x y')
            sage: v = M.vector_field('v')
            sage: v[:] = (-y, x)
            sage: w = M.vector_field()
            sage: w[:] = (2*x+y, x*y)
            sage: w.lie_der(v)
            vector field on the 2-dimensional manifold 'M'
            sage: w.lie_der(v).view()
            ((x - 2)*y + x) d/dx + (x^2 - y^2 - 2*x - y) d/dy

        The Lie derivative is antisymmetric::
        
            sage: w.lie_der(v) == -v.lie_der(w)
            True
            
        For vectors, it coincides with the commutator::

            sage: f = M.scalar_field(x^3 + x*y^2)
            sage: w.lie_der(v)(f).view()
            (x, y) |--> -(x + 2)*y^3 + 3*x^3 - x*y^2 + 5*(x^3 - 2*x^2)*y
            sage: w.lie_der(v)(f) == v(w(f)) - w(v(f))  # rhs = commutator [v,w] acting on f
            True
            
        Lie derivative of a 1-form::
        
            sage: om = M.one_form()
            sage: om[:] = (y^2*sin(x), x^3*cos(y))
            sage: om.lie_der(v)
            1-form on the 2-dimensional manifold 'M'
            sage: om.lie_der(v).view()
            (-y^3*cos(x) + x^3*cos(y) + 2*x*y*sin(x)) dx + (-x^4*sin(y) - 3*x^2*y*cos(y) - y^2*sin(x)) dy
            
        Check of Cartan identity::
        
            sage: om.lie_der(v) == v.contract(0,om.exterior_der(),0) + (om(v)).exterior_der()
            True
        
        """
        from vectorfield import VectorFieldParal
        if not isinstance(vector, VectorFieldParal):
            raise TypeError("The argument must be of type VectorFieldParal.")
        if id(vector) not in self._lie_derivatives:
            # A new computation must be performed
            #
            # 1/ Search for a common coordinate frame:
            coord_frame = self.common_coord_frame(vector)
            if coord_frame is None:
                raise TypeError("No common coordinate frame found.")
            chart = coord_frame.chart
            #
            # 2/ Component computation:
            tc = self.components[coord_frame]
            vc = vector.components[coord_frame]
            # the result has the same tensor type and same symmetries as self:
            resc = self._new_comp(coord_frame) 
            n_con = self.tensor_type[0]
            vf_module = vector.fmodule
            for ind in resc.non_redundant_index_generator():
                rsum = 0
                for i in vf_module.irange():
                    rsum += vc[[i]].function_chart(chart) * \
                           tc[[ind]].function_chart(chart).diff(i)
                # loop on contravariant indices:
                for k in range(n_con): 
                    for i in vf_module.irange():
                        indk = list(ind)
                        indk[k] = i  
                        rsum -= tc[[indk]].function_chart(chart) * \
                                vc[[ind[k]]].function_chart(chart).diff(i)
                # loop on covariant indices:
                for k in range(n_con, self.tensor_rank): 
                    for i in vf_module.irange():
                        indk = list(ind)
                        indk[k] = i  
                        rsum += tc[[indk]].function_chart(chart) * \
                                vc[[i]].function_chart(chart).diff(ind[k])
                resc[[ind]] = rsum.scalar_field()
            #
            # 3/ Final result (the tensor)
            res = vf_module.tensor_from_comp(self.tensor_type, resc)
            self._lie_derivatives[id(vector)] = (vector, res)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]


